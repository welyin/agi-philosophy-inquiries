"""1221: a Gauss-matched variational comparison, not a universal channel.

Numerics use genuine original loop characters (trivial / weak adjoint) and
the four original sterile-neutrino CAR modes at two vertices. The finite
matrix is the exact compression of the old H up to a common scalar. Its
Gibbs state is NOT the full Gauss Gibbs state or an invariant-H sector.
The full-model finite-energy theorem is proved in proof.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
TARGET = HERE / "results.json"
SOURCES = [
    "archive_554_584/research_note_574.md",
    "archive_554_584/569/joint_quotient_gauge_completion_results.json",
    "archive_585_628/research_note_589.md",
    "archive_585_628/research_note_598.md",
    "archive_585_628/research_note_617.md",
    "archive_585_628/617/joint_region_energy_gluing_results.json",
    "archive_629_652/research_note_637.md",
    "archive_629_652/research_note_642.md",
    "archive_629_652/research_note_644.md",
    "archive_629_652/644/joint_region_replica_source.py",
]
SHAPE = np.array([[.24, .31, -.12], [.31, -.07, .17], [-.12, .17, -.17]])
BETA, HOP, EPS, T0, SIGMA0 = 1.7, .17, .73, .8, .12


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hermitian_function(a, function):
    w, u = np.linalg.eigh((a + a.conj().T) / 2)
    return (u * function(w)) @ u.conj().T


def entropy(a):
    w = np.linalg.eigvalsh((a + a.conj().T) / 2)
    assert w.min() > -2e-13
    w = w[w > 1e-14]
    return float(-w @ np.log(w))


def log_positive(a):
    w = np.linalg.eigvalsh(a)
    assert w.min() > 1e-12, "These finite diagnostic references are strictly positive."
    return hermitian_function(a, np.log)


def relative(a, b):
    return float(np.trace(a @ (log_positive(a) - log_positive(b))).real)


def gibbs(h):
    w, u = np.linalg.eigh(h)
    p = np.exp(-BETA * (w - w.min()))
    p /= p.sum()
    return (u * p) @ u.conj().T


def expectation(a, h):
    return float(np.trace(a @ h).real)


def annihilators(n):
    # Binary order: the first two modes are A, the last two B.
    result = []
    for k in range(n):
        c = np.zeros((1 << n, 1 << n), complex)
        mask = 1 << (n - 1 - k)
        for j in range(1 << n):
            if j & mask:
                occupied_before = (j >> (n - k)).bit_count()
                c[j ^ mask, j] = (-1) ** occupied_before
        result.append(c)
    return result


def margins(a):
    r = a.reshape(4, 4, 4, 4)
    return np.einsum("abcb->ac", r), np.einsum("abad->bd", r)


def compare(a):
    """Nonlinear state-dependent comparison. NOT a CPTP map on inputs."""
    tau = np.zeros_like(a)
    pinched = np.zeros_like(a)
    ps, conditional_mi = [], 0.
    for k in range(2):
        sl = slice(16 * k, 16 * (k + 1))
        block = a[sl, sl]
        p = float(np.trace(block).real)
        assert p > 0
        q = block / p
        qa, qb = margins(q)
        tau[sl, sl] = p * np.kron(qa, qb)
        pinched[sl, sl] = block
        ps.append(p)
        conditional_mi += p * (entropy(qa) + entropy(qb) - entropy(q))
    return tau, pinched, np.array(ps), float(conditional_mi)


def physical_cut(a):
    # Two crossing original links: Schmidt ranks 1 and 3^2, not new qubits.
    # Regional multiplicities contain the original two sterile CAR modes.
    offsets, ranks, d = [0, 4], [1, 9], 40
    v = np.zeros((d * d, 32), complex)
    for sector, rank in enumerate(ranks):
        for ia in range(4):
            for ib in range(4):
                col = sector * 16 + ia * 4 + ib
                for edge in range(rank):
                    ar = offsets[sector] + ia * rank + edge
                    br = offsets[sector] + ib * rank + edge
                    v[ar * d + br, col] = 1 / np.sqrt(rank)
    assert np.max(abs(v.conj().T @ v - np.eye(32))) < 2e-15
    # Avoid materialising a 1600^2 density. Trace via V's row tensor.
    vt = v.reshape(d, d, 32)
    ra = np.einsum("abi,ij,cbj->ac", vt, a, vt.conj(), optimize=True)
    rb = np.einsum("abi,ij,acj->bc", vt, a, vt.conj(), optimize=True)
    return ra, rb


def coefficients(t=T0, sigma=SIGMA0):
    old = json.loads((RESEARCH / SOURCES[1]).read_text("utf-8"))["evidence"]["units_check"]
    b = old["electric_coefficients_colour_weak_integer_circle"][1]
    alpha = old["residual_positive_coefficients"][1]
    shape = hermitian_function(t * SHAPE, np.exp).real
    prefactor = np.exp(-2 * sigma) / EPS
    # Original four-edge square: the two outgoing edges at its origin give
    # P1.P2=-Casimir. No off-diagonal electric term has been discarded.
    kappa = b * prefactor * (2 * shape[0, 0] + 2 * shape[1, 1] - 2 * shape[0, 1])
    # det(shape)=1: the (12,12) inverse-shape two-form coefficient is shape33.
    g = alpha * prefactor * shape[2, 2]
    return float(kappa), float(g), shape, old


def matrices(t=T0, sigma=SIGMA0):
    kappa, g, shape, old = coefficients(t, sigma)
    fusion = np.array([[0., 1.], [1., 1.]])
    electric = np.diag([0., 2 * kappa])
    magnetic = -g * fusion
    c = annihilators(4)
    car = HOP * sum(c[k].conj().T @ c[k + 2] + c[k + 2].conj().T @ c[k]
                    for k in range(2))
    hbulk = np.kron(electric, np.eye(16))
    hboundary = np.kron(magnetic, np.eye(16)) + np.kron(np.eye(2), car)
    return hbulk + hboundary, hbulk, hboundary, electric, magnetic, car


def group_and_original_compression():
    x, w = np.polynomial.legendre.leggauss(160)
    theta, weight = (x + 1) * np.pi / 2, w * np.sin((x + 1) * np.pi / 2) ** 2
    chi = 1 + 2 * np.cos(2 * theta)  # SU(2) spin-one character.
    basis = np.stack([np.ones_like(chi), chi])
    gram = np.einsum("ik,k,jk->ij", basis, weight, basis)
    fusion = np.einsum("ik,k,jk->ij", basis, weight * chi, basis)
    assert np.max(abs(gram - np.eye(2))) < 2e-13
    assert np.max(abs(fusion - [[0, 1], [1, 1]])) < 2e-13
    # The weak adjoint is a legal irrep of the ORIGINAL quotient group.
    assert (2 * (0 + 2 * 0) + 3 * 2 + 0) % 6 == 0
    c = annihilators(4)
    car_error = max(float(np.max(abs(c[i] @ c[j].conj().T + c[j].conj().T @ c[i]
                    - (np.eye(16) if i == j else 0)))) for i in range(4) for j in range(4))
    assert car_error == 0
    kappa, g, shape, old = coefficients()
    original = json.loads((RESEARCH / SOURCES[5]).read_text("utf-8"))["electric"]
    k = old["electric_coefficients_colour_weak_integer_circle"][1] / EPS * np.exp(-2 * SIGMA0) * shape
    coefficient_error = float(np.max(abs(k - np.array(original["coefficient_matrix"]))))
    assert coefficient_error < 1e-14
    return dict(original_labels=[[0, 0, 0, 0], [0, 0, 2, 0]],
        dimensions=[1, 3], cut_schmidt_ranks=[1, 9], original_car_slots=[30, 31, 62, 63],
        haar_gram_error=float(np.max(abs(gram - np.eye(2)))),
        haar_fusion_error=float(np.max(abs(fusion - [[0, 1], [1, 1]]))),
        car_error=car_error, inherited_electric_coefficients_error=coefficient_error,
        kappa=kappa, magnetic_g=g, beta=BETA, original_declared_hopping=HOP,
        exact_compressed_H_up_to_common_scalar=True,
        full_Gauss_Gibbs=False, H_invariant_subspace=False,
        full_H_approximation_error_claimed=False)


def entropy_and_matching():
    h, hbulk, hb, he, hm, hc = matrices()
    rho = gibbs(h)
    tau, pinched, p, mi = compare(rho)
    ra, rb = physical_cut(rho)
    ta, tb = physical_cut(tau)
    marginal_error = max(float(np.max(abs(ra-ta))), float(np.max(abs(rb-tb))))
    c = relative(rho, tau)
    coh = relative(rho, pinched)
    reverse = relative(tau, rho)
    delta_e = expectation(tau-rho, h)
    bulk_gap = expectation(tau-rho, hbulk)
    center = float(-p @ np.log(p))
    matching_one_side = float(p @ np.log([1., 9.]))
    iext = entropy(ra) + entropy(rb) - entropy(rho)
    errors = [marginal_error, abs(c-(entropy(tau)-entropy(rho))), abs(c-coh-mi),
              abs(c+reverse-BETA*delta_e), abs(bulk_gap),
              abs(iext-center-2*matching_one_side-c)]
    assert max(errors) < 8e-12
    assert min(c, coh, mi, reverse, delta_e) > 1e-6
    # An unrestricted product in the extended spaces leaves the matched image.
    # Within sector lambda the two edge marginals have rank r; singlet overlap
    # of their product is 1/r^2. A/B-sector mismatch is also excluded.
    unmatched = 1 - float(np.sum(p**2 / np.array([1., 9.])**2))
    assert unmatched > .1
    return dict(sector_probabilities=p.tolist(),
        correlation_relative_entropy=c, sector_coherence=coh,
        conditional_multiplicity_mutual_information=mi,
        reverse_relative_entropy=reverse, entropy_gain=entropy(tau)-entropy(rho),
        boundary_energy_gap=delta_e, bulk_energy_gap=bulk_gap,
        center_Shannon=center, one_side_forced_matching_entropy=matching_one_side,
        extended_mutual_information=iext,
        naive_extended_product_outside_matching_weight=unmatched,
        maximum_identity_error=max(errors), comparison_is_universal_CPTP=False)


def complete_boundary_source():
    h, bulk, hb, he, hm, hc = matrices()
    rho = gibbs(h)
    tau = compare(rho)[0]
    diff = tau-rho
    kappa, g, shape, old = coefficients()
    ds = SHAPE @ shape
    common = np.exp(-2*SIGMA0) / EPS
    dk = old["electric_coefficients_colour_weak_integer_circle"][1] * common * (
        2*ds[0, 0] + 2*ds[1, 1] - 2*ds[0, 1])
    dg = old["residual_positive_coefficients"][1] * common * ds[2, 2]
    fusion = np.array([[0., 1.], [1., 1.]])
    full_sources = [np.kron(np.diag([0., 2*dk]), np.eye(16)) - np.kron(dg*fusion, np.eye(16)),
                    -2*np.kron(he+hm, np.eye(16))]
    boundary_sources = [-np.kron(dg*fusion, np.eye(16)), -2*np.kron(hm, np.eye(16))]
    errors, source_gaps = [], []
    for i in range(2):
        args = [T0, SIGMA0]
        step = 2e-6
        args[i] += step
        hp = matrices(*args)[0]
        args[i] -= 2*step
        hn = matrices(*args)[0]
        derivative = (hp-hn)/(2*step)
        errors.append(float(np.max(abs(derivative-full_sources[i]))))
        errors.append(abs(expectation(diff, full_sources[i]-boundary_sources[i])))
        source_gaps.append(expectation(diff, full_sources[i]))
    assert max(errors) < 2e-9
    # Moving the actual comparison pair adds its state derivative.
    def gap(sigma):
        hh = matrices(T0, sigma)[0]
        rr = gibbs(hh)
        tt = compare(rr)[0]
        return expectation(tt-rr, hh), tt-rr
    step = 2e-5
    ep, dp = gap(SIGMA0+step)
    em, dm = gap(SIGMA0-step)
    total = (ep-em)/(2*step)
    state = expectation((dp-dm)/(2*step), h)
    chain_error = abs(total-source_gaps[1]-state)
    assert chain_error < 2e-8
    assert abs(state) > .01
    # Independently specified face area, no entropy definition. It is still
    # given background geometry, NOT an actual internal ruler certificate.
    inverse = np.linalg.inv(shape)
    face = EPS**2*np.exp(4*SIGMA0)*np.sqrt(inverse[0, 0])
    return dict(fixed_state_shear_source_gap=source_gaps[0],
        fixed_state_conformal_source_gap=source_gaps[1],
        coefficient_and_bulk_cancellation_max_error=max(errors),
        moving_pair_total_energy_gap_derivative=total,
        indispensable_state_derivative=state, chain_rule_error=chain_error,
        given_background_dual_face_area=float(face),
        independent_internal_area_measurement_proved=False,
        entropy_area_response_proved=False)


def nonlinear_comparison_witness():
    # The same four original sterile CAR modes, even physical vectors. Only
    # the first spin at each vertex is used; no new Hilbert system is added.
    empty = np.zeros(32, complex)
    both = np.zeros(32, complex)
    empty[0] = 1
    both[10] = 1  # occupations 1010: one neutral sterile mode in A and B.
    # Positive weight on the second original flux sector keeps compare's
    # diagnostic inputs faithful at the sector level.
    base = .2*np.eye(32)/32
    r0 = base + .8*np.outer(empty, empty.conj())
    r1 = base + .8*np.outer(both, both.conj())
    midpoint = compare((r0+r1)/2)[0]
    average = (compare(r0)[0]+compare(r1)[0])/2
    defect = float(np.linalg.norm(midpoint-average))
    assert defect > .1
    return dict(non_affinity_Frobenius_defect=defect,
        same_original_Gauss_and_even_CAR_support=True,
        universal_unknown_state_channel_excluded=True)


def run():
    groups = dict(original_physical_compression=group_and_original_compression(),
                  matched_entropy_energy=entropy_and_matching(),
                  original_boundary_sources=complete_boundary_source(),
                  no_universal_channel=nonlinear_comparison_witness())
    return dict(schema="round1221_gauss_matched_comparison_v1", status="PASS", check_groups=4,
                adopted_axioms=0, empirical_calibrations=0,
                source_sha256={s:digest(RESEARCH/s) for s in SOURCES},
                code_sha256=digest(Path(__file__).resolve()), groups=groups,
                full_model_claim="analytic fixed-graph finite-form-energy theorem; proof.md",
                numerical_scope="exact original physical Galerkin compression, not full thermal approximation",
                TG3_complete=False)


def numeric_close(a, b):
    if isinstance(a, dict):
        return a.keys()==b.keys() and all(numeric_close(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a)==len(b) and all(numeric_close(x,y) for x,y in zip(a,b))
    if isinstance(a, float):
        return abs(a-b) <= 1e-11 * max(1., abs(a), abs(b))
    return a == b


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        TARGET.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    if args.verify:
        assert numeric_close(result, json.loads(TARGET.read_text("utf-8"))), "saved-result mismatch"
    print(json.dumps(result, ensure_ascii=False, indent=2))
