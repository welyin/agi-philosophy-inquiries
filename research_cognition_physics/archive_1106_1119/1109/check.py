"""Finite four-momentum relational-clock boundary instrument; no port simulation.

Default recomputes and compares results without writing. The sole numerical
material witness has N=3; general quantifiers are proved in proof.md.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]


def norm(a):
    return float(np.linalg.norm(a, 2))


def spectral(a, f):
    e, v = np.linalg.eigh(a)
    return (v * f(e)) @ v.conj().T


def mean(a, psi):
    return float(np.vdot(psi, a @ psi).real)


def boost(beta):
    beta = np.asarray(beta, float)
    b2 = float(beta @ beta)
    g = 1 / math.sqrt(1-b2)
    out = np.eye(4)
    out[0, 0] = g
    out[0, 1:] = g * beta
    out[1:, 0] = g * beta
    if b2:
        out[1:, 1:] += (g-1) * np.outer(beta, beta) / b2
    return out


def run():
    nref, omega, mass = 3, 1., 2.
    uc, ua = np.array([1.25, .75, 0., 0.]), np.array([1.25, -.75, 0., 0.])
    eta_metric = np.diag([1., -1., -1., -1.])
    delta = omega * (uc-ua)
    pca = np.array([(mass+omega*c)*uc+(mass+omega*a)*ua
                    for c in range(2) for a in range(2)])
    q = np.array([[2.5, (n-1)*1.5, 0., 0.] for n in range(nref)])
    mref = np.sqrt(np.einsum('ni,ij,nj->n', q, eta_metric, q))
    assert np.allclose(mref, [2., 2.5, 2.]) and np.all(q[:, 0] > 0)
    b = np.zeros((4, 4), complex)
    b[1, 2], b[2, 1] = 1j, -1j
    bare = [norm(b @ np.diag(pca[:, mu])-np.diag(pca[:, mu]) @ b)
            for mu in range(4)]
    assert np.allclose(bare, [0., 1.5, 0., 0.])
    g1 = dict(clock_velocities=[uc.tolist(), ua.tolist()], delta_four_momentum=delta.tolist(),
              bare_readout_commutator_norms=bare,
              energy_only_check_misses_momentum_violation=True)

    dim = 4*nref
    bt = np.zeros((dim, dim), complex)
    for n in range(nref-1):
        a, d = nref+n+1, 2*nref+n
        bt[a, d], bt[d, a] = 1j, -1j
    four = [np.diag((pca[:, mu, None]+q[None, :, mu]).reshape(-1)) for mu in range(4)]
    comm = [norm(bt @ p-p @ bt) for p in four]
    ap, am = [spectral(bt, lambda e, s=s: np.sqrt(np.clip((1+s*e)/2, 0, 1)))
              for s in (1, -1)]
    y = np.array([[0, -1j], [1j, 0]])
    x = np.array([[0, 1], [1, 0]])
    n = np.diag([0., 1.])
    ur = np.kron(ap, np.eye(2))-1j*np.kron(am, y)
    um = np.kron(np.eye(2)-n, np.eye(2))-1j*np.kron(n, x)
    w = np.kron(np.eye(dim), um) @ np.kron(ur, np.eye(2))
    wp = [np.kron(p, np.eye(4)) for p in four]
    wi = [norm(w @ p-p @ w) for p in wp]
    assert max(comm+wi+[norm(w.conj().T @ w-np.eye(4*dim))]) < 1e-12
    g2 = dict(reference_four_momenta=q.tolist(), reference_rest_masses=mref.tolist(),
              reference_rest_internal_energies=(mref-2).tolist(),
              lifted_effect_norm=norm(bt), lifted_four_commutators=comm,
              full_instrument_dimension=4*dim, instrument_four_commutators=wi,
              full_unitarity_error=norm(w.conj().T @ w-np.eye(4*dim)))

    eta = np.ones(nref)/math.sqrt(nref)
    ready = np.array([1., 0., 0., 0.])
    done = np.array([0., 0., 0., 1.])
    j = np.kron(np.eye(4), eta[:, None])
    rho_eta = np.outer(eta, eta)
    eff = (np.eye(dim)+bt)/2
    reduced_b_error = norm(j.conj().T @ bt @ j-(1-1/nref)*b)
    rows = []
    for phi in (-math.pi/6, 0., math.pi/6):
        psi = np.array([0, 1, np.exp(-1j*phi), 0])/math.sqrt(2)
        out = (w @ np.kron(j @ psi, ready)).reshape(dim, 2, 2)
        plus = float(np.sum(abs(out[:, 0, 0])**2))
        mismatch = float(np.sum(abs(out[:, 0, 1])**2)+np.sum(abs(out[:, 1, 0])**2))
        rho_mix = np.kron(np.outer(psi, psi.conj()), np.eye(nref)/nref)
        dephased = float(np.trace(eff @ rho_mix).real)
        expected = (1+(1-1/nref)*np.sin(phi))/2
        assert abs(plus-expected) < 1e-12 and abs(dephased-.5) < 1e-12 and mismatch < 1e-12
        rows.append(dict(phase=phi, plus=plus, expected=float(expected),
                         record_disagreement=mismatch, incoherent_reference_plus=dephased))
    assert reduced_b_error < 1e-12
    g3 = dict(reference_overlap_factor=1-1/nref, all_input_effect_identity_error=reduced_b_error,
              cases=rows, finite_phase_window=[-math.pi/3, math.pi/3],
              minimum_probability_slope=(1-1/nref)*omega/4,
              probability_to_relative_time_error=4/((1-1/nref)*omega),
              two_phase_single_record_discrimination_error=1/3)

    a0p, a0m = [spectral(b, lambda e, s=s: np.sqrt(np.clip((1+s*e)/2, 0, 1)))
                for s in (1, -1)]
    vin = np.kron(j, ready[:, None])
    actual = w @ vin
    ideal = np.column_stack([np.kron(np.kron(a0p[:, a], eta), ready)
                              -1j*np.kron(np.kron(a0m[:, a], eta), done) for a in range(4)])
    overlap = 1-(1-1/math.sqrt(2))/nref
    overlap_matrix = ideal.conj().T @ actual
    overlap_error = norm(overlap_matrix-np.diag([1, overlap, overlap, 1]))
    iso_error = norm(actual-ideal)
    complete_error = math.sqrt(1-overlap**2)
    witness = np.array([0., 1., 0., 0.])
    witness_error = math.sqrt(max(0., 1-abs(np.vdot(ideal @ witness, actual @ witness))**2))
    assert overlap_error < 1e-12 and abs(witness_error-complete_error) < 1e-12
    assert abs(iso_error-math.sqrt((2-math.sqrt(2))/nref)) < 1e-12
    g4 = dict(complete_output_overlap_error=overlap_error, minimum_complete_overlap=overlap,
              complete_isometry_difference=iso_error,
              analytic_complete_half_diamond=complete_error, pure_witness_distance=witness_error,
              effect_probability_error_bound=1/(2*nref),
              reference_and_both_records_kept=True)

    psi = np.kron(np.array([0., 0., 1., 0.]), eta)
    rho_out = sum(a @ np.outer(psi, psi.conj()) @ a.conj().T for a in (ap, am))
    recoil = []
    for mu in range(4):
        pr = np.kron(np.eye(4), np.diag(q[:, mu]))
        pc = np.kron(np.diag(pca[:, mu]), np.eye(nref))
        dr = float(np.trace(pr @ rho_out).real-mean(pr, psi))
        dc = float(np.trace(pc @ rho_out).real-mean(pc, psi))
        assert abs(dr-(1-1/nref)*delta[mu]/2) < 1e-12 and abs(dr+dc) < 1e-12
        recoil.append(dict(component=mu, reference_gain=dr, clock_change=dc))
    lam = boost([.6, 0., 0.])
    assert norm(lam.T @ eta_metric @ lam-eta_metric) < 1e-12
    qt = q @ lam.T
    pt = pca @ lam.T
    boosts_comm = [norm(bt @ np.diag((pt[:, mu, None]+qt[None, :, mu]).reshape(-1))
                        -np.diag((pt[:, mu, None]+qt[None, :, mu]).reshape(-1)) @ bt)
                   for mu in range(4)]
    variance_before = float(np.var(q[:, 0]))
    variance_after = float(np.var(qt[:, 0]))
    assert abs(variance_after-27/32) < 1e-12 and max(boosts_comm) < 1e-12
    g5 = dict(unconditional_recoil_for_10=recoil,
              reference_energy_variance_cm=variance_before,
              reference_momentum_variance_cm=float(np.var(q[:, 1])),
              reference_energy_variance_boosted=variance_after,
              boosted_four_commutators=boosts_comm,
              mass_invariance_error=float(np.max(abs(np.einsum('ni,ij,nj->n', qt, eta_metric, qt)-mref*mref))))

    swap = np.zeros((4, 4))
    for a in range(2):
        for b0 in range(2): swap[2*b0+a, 2*a+b0] = 1
    # Physical factors X and Y are different; Y is initially blank.
    jq = np.kron(np.eye(2), np.array([[1.], [0.]]))
    targetq = np.column_stack([np.kron(np.array([1., 0.]), np.eye(2)[:, a]) for a in range(2)])
    factor_error = norm(swap @ jq-targetq)
    both = np.kron(w, swap)
    # Positive baseline four-momenta of the two payload carriers and records.
    baseline = 7*uc+8*ua+np.array([9., 0., 0., 0.])
    joint_four = [np.kron(p, np.eye(4))+baseline[mu]*np.eye(16*dim)
                  for mu, p in enumerate(wp)]
    conservation = [norm(both @ p-p @ both) for p in joint_four]
    assert factor_error == 0. and max(conservation) < 1e-11
    g6 = dict(actual_payload_factor_swap_error=factor_error,
              complete_local_boundary_dimension=16*dim,
              payload_carrier_and_record_baseline_four_momentum=baseline.tolist(),
              full_joint_four_commutators=conservation,
              all_passive_references_protected_by_factor_identity=True,
              incoming_port_response_implemented=False,
              normal_localized_preparation_implemented=False)
    paths = ['archive_370_428/research_note_421.md', 'archive_935_955/research_note_939.md',
             'archive_956_989/research_note_976.md', 'archive_1106_/1108/proof.md',
             'archive_1106_/_shared/clock_instrument_admission/README.md']
    return dict(round=1109, status='PASS_MOVING_CLOCK_BOUNDARY_INSTRUMENT',
                groups=dict(bare_violation=g1, positive_reference_and_instrument=g2,
                            clock_record=g3, complete_error=g4, recoil_and_boost=g5,
                            joint_payload_receive=g6),
                source_sha256={p: hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in paths},
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                new_adopted_axioms=0, scientific_calibration_increment=0,
                goal_complete=False, actual_SR_common_model_complete=False)


def compare(a, b):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys()
        for k in a: compare(a[k], b[k])
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b)
        for x, y in zip(a, b): compare(x, y)
    elif isinstance(a, float):
        assert math.isfinite(a) and math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-11), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    data = run()
    destination = HERE/'results.json'
    if args.write_results:
        with destination.open('x', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
            f.write('\n')
    compare(data, json.loads(destination.read_text(encoding='utf-8')))
    print(json.dumps(dict(round=1109, status=data['status'], groups=len(data['groups']),
                          saved_results_match=True, scientific_calibration_increment=0)))
