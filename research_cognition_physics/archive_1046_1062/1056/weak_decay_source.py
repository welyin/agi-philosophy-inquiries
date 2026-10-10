"""Round 1056: joint electron record and total outgoing angular energy flux.

Exact fractions certify the finite-window bound. Numerical quadratures are
independent checks in two three-body coordinate systems, not an error bound
for full QFT, a material detector, or finite-frequency gravity.
"""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import argparse
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / "results.json"
HISTORY = [
    "archive_301_341/research_note_326.md",
    "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
    "archive_956_989/981/drafts/common_parent_contract_v1.md",
    "archive_956_989/research_note_985.md",
    "archive_956_989/research_note_986.md",
    "archive_1046_/research_note_1047.md",
    "archive_1046_/research_note_1053.md",
    "archive_1046_/_shared/muon_decay_readout_adoption.md",
    "archive_1046_/_shared/muon_decay_readout_adoption_review.md",
    "archive_1046_/_shared/muon_decay_readout_adoption/results.json",
    "archive_1046_/_shared/operational_adoptions_after1055_checks.json",
]
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def gauss(n, low, high):
    x, w = leggauss(n)
    return low + (x + 1) * (high - low) / 2, w * (high - low) / 2


def exact_certificate():
    lo, hi = F(2, 5), F(1, 2)
    alpha = (hi**3 - hi**4/2) - (lo**3 - lo**4/2)
    zeta = (hi**4/4 - hi**3/6) - (lo**4/4 - lo**3/6)
    assert alpha == F(851, 20000) and zeta == -F(113, 120000)
    assert 0 < alpha-abs(zeta) < alpha+abs(zeta) < 1
    beta_lo, beta_hi = lo/(2-lo), hi/(2-hi)
    assert (beta_lo, beta_hi) == (F(1, 4), F(1, 3))
    # For every beta in this interval: N0 <= 0, N2 <= -1/6,
    # N4 <= -8/27, and 0 < 1-beta^2*c^2 <= 1.
    n2_upper = -F(1, 2)+beta_hi
    n4_upper = beta_hi**3+beta_hi/2-F(1, 2)
    assert n2_upper == -F(1, 6) and n4_upper == -F(8, 27)
    denominator_lower = 1-beta_hi**2
    integral_upper = n2_upper*F(2, 3)  # dropping strictly negative N4
    prefactor_lower = lo**2*(3-2*hi)/16
    k_upper = prefactor_lower*integral_upper
    c_upper = (hi-lo)*k_upper
    assert (denominator_lower, integral_upper, prefactor_lower, c_upper) == (
        F(8, 9), -F(1, 9), F(1, 50), -F(1, 4500))
    source_grid_step = F(1, 100000)
    grid_error = source_grid_step/2
    recovery_lower = -c_upper-grid_error
    assert recovery_lower > F(1, 5000)
    # Energy-weighted Michel polynomials. Values here multiply m/(4*pi).
    charged_energy_scalar = F(3, 4)-F(2, 5)
    charged_energy_spin = F(2, 5)-F(1, 4)
    electron_neutrino_scalar = 6*(F(1, 4)-F(1, 5))
    electron_neutrino_spin = -electron_neutrino_scalar
    assert 2*charged_energy_scalar+electron_neutrino_scalar == 1
    assert 2*charged_energy_spin+electron_neutrino_spin == 0
    # The logarithm check is an exact series interval, not rounded ln(2).
    ln2_lo = 2*sum((F(1, 3)**(2*j+1))/F(2*j+1) for j in range(24))
    ln2_tail = 2*F(1, 3)**49/(49*(1-F(1, 9)))
    fixed_direction_lower = 25-36*(ln2_lo+ln2_tail)
    assert fixed_direction_lower > F(46, 1000)
    values = dict(record_probability=alpha, record_z_coefficient=zeta,
                  beta_low=beta_lo, beta_high=beta_hi,
                  N2_upper=n2_upper, N4_upper=n4_upper,
                  denominator_lower=denominator_lower,
                  integral_upper=integral_upper,
                  prefactor_lower=prefactor_lower,
                  window_source_coefficient_upper=c_upper,
                  joint_moment_gap_lower=-2*c_upper,
                  record_only_TV_error_lower=-c_upper,
                  rounded_source_error_bound=grid_error,
                  rounded_record_only_TV_error_lower=recovery_lower,
                  charged_energy_scalar=charged_energy_scalar,
                  charged_energy_spin=charged_energy_spin,
                  electron_neutrino_energy_scalar=electron_neutrino_scalar,
                  electron_neutrino_energy_spin=electron_neutrino_spin,
                  fixed_direction_coefficient_lower=fixed_direction_lower)
    return {k: str(v) for k, v in values.items()}


def pair_coefficient(n):
    xs, wx = gauss(n, .4, .5)
    cs, wc = gauss(n, -1., 1.)
    values = []
    for x in xs:
        beta = x/(2-x)
        numerator = (beta*(3*beta-1)/2
                     +(-.5+beta-beta**2/2-2*beta**3)*cs**2
                     +(beta**3+beta/2-.5)*cs**4)
        integral = np.dot(wc, numerator/(1-beta**2*cs**2))
        values.append(x*x*(3-2*x)*integral/(16*(1-beta/3)))
    return float(np.dot(wx, values))


def dalitz_k(a, n):
    b, weights = gauss(n, 1-a, 1.)
    c = 2-a-b
    u = 1-2*(a+b-1)/(a*b)
    v2 = 1-u*u
    A = (a+b*u*u+(a+b*u)**2/c)/2
    C = b*v2*(2-a)/(2*c)
    Bv = b*v2*((2-a)*u+a)/(2*c)
    D = Bv/2
    L = u*(A-C/2)-2*D
    # Independent electron-sphere integrals: <nz>_H=1/4,
    # <nx^2*nz>_H=1/16, with normalized dOmega/(4*pi).
    return float(np.dot(weights, -12*b*(1-b)*(L/16+D/4)))


def dalitz_coefficient(n):
    a, w = gauss(n, .4, .5)
    return float(np.dot(w, [dalitz_k(x, n) for x in a]))


def direct_momentum_check():
    """Integrate actual three momenta, rather than the reduced N_beta kernel.

    The non-smooth sign partition is only a finite CP/reference diagnostic;
    the unbinned energy moment has its analytic certificate above.
    """
    az = (np.arange(12)+.5)*2*np.pi/12
    costheta, wt = gauss(10, 0., 1.)
    ct, ph, chi = np.meshgrid(costheta, az, az, indexing="ij")
    st = np.sqrt(1-ct*ct)
    nvec = np.stack((st*np.cos(ph), st*np.sin(ph), ct), axis=-1)
    e_theta = np.stack((ct*np.cos(ph), ct*np.sin(ph), -st), axis=-1)
    e_phi = np.stack((-np.sin(ph), np.cos(ph), np.zeros_like(ph)), axis=-1)
    transverse = np.cos(chi)[..., None]*e_theta+np.sin(chi)[..., None]*e_phi
    angle_weight = np.broadcast_to(wt[:, None, None]/(2*12*12), ct.shape)
    effects = np.zeros(4)
    moment = np.zeros(4)
    positive_q = np.zeros(4)
    residual = 0.
    a_values, wa_values = gauss(20, .4, .5)
    for a, wa in zip(a_values, wa_values):
        b_values, wb_values = gauss(20, 1-a, 1.)
        for b, wb in zip(b_values, wb_values):
            c = 2-a-b
            u = 1-2*(a+b-1)/(a*b)
            nq = u*nvec+math.sqrt(max(0., 1-u*u))*transverse
            nr = -(a*nvec+b*nq)/c
            residual = max(residual, float(np.max(abs(np.sum(nr*nr, axis=-1)-1))),
                           float(np.max(abs(a*nvec+b*nq+c*nr))))
            Qxz = (a*nvec[..., 0]*nvec[..., 2]
                   +b*nq[..., 0]*nq[..., 2]+c*nr[..., 0]*nr[..., 2])/2
            assert float(np.max(abs(Qxz))) <= .5+1e-12
            pauli = np.concatenate((np.ones(ct.shape+(1,)), -nq), axis=-1)
            weight = wa*wb*12*b*(1-b)*angle_weight
            effects += np.sum(weight[..., None]*pauli, axis=(0, 1, 2))
            moment += np.sum((weight*Qxz)[..., None]*pauli, axis=(0, 1, 2))
            positive_q += np.sum((weight*(Qxz >= 0))[..., None]*pauli, axis=(0, 1, 2))
    target = np.array([851/20000, 0, 0, -113/120000])
    target_moment = np.array([0, pair_coefficient(64), 0, 0])
    effect_residual = float(np.max(abs(effects-target)))
    moment_residual = float(np.max(abs(moment-target_moment)))
    assert residual < 2e-13 and effect_residual < 2e-13 and moment_residual < 2e-13
    return effects, positive_q, dict(
        max_on_shell_and_momentum_residual=residual,
        direct_effect_residual=effect_residual,
        direct_source_moment_residual=moment_residual,
        direct_pauli_effect=effects.tolist(), direct_source_moment=moment.tolist())


def matrix(coefficients):
    return coefficients[0]*I2+np.einsum("i,ijk->jk", coefficients[1:], SIGMA)


def reference_diagnostic(total_b, positive_q):
    effects = [matrix(positive_q), matrix(total_b-positive_q), I2-matrix(total_b)]
    eigen_min = min(float(np.linalg.eigvalsh(e)[0]) for e in effects)
    assert eigen_min > -1e-12
    assert np.linalg.norm(sum(effects)-I2) < 1e-12
    rng = np.random.default_rng(1056)
    largest = 0.
    for _ in range(16):
        x = rng.normal(size=(6, 6))+1j*rng.normal(size=(6, 6))
        rho = x@x.conj().T
        rho /= np.trace(rho)
        original_r = np.trace(rho.reshape(2, 3, 2, 3), axis1=0, axis2=2)
        post_r = np.zeros((3, 3), complex)
        for e in effects:
            ev, v = np.linalg.eigh(e)
            sq = (v*np.sqrt(np.maximum(ev, 0)))@v.conj().T
            k = np.kron(sq, np.eye(3))
            branch = np.trace((k@rho@k.conj().T).reshape(2, 3, 2, 3), axis1=0, axis2=2)
            linear_branch = np.trace((np.kron(e, np.eye(3))@rho).reshape(2, 3, 2, 3),
                                     axis1=0, axis2=2)
            assert np.linalg.eigvalsh(branch)[0] > -1e-12
            largest = max(largest, float(np.linalg.norm(branch-linear_branch)))
            post_r += branch
        largest = max(largest, float(np.linalg.norm(post_r-original_r)))
    assert largest < 1e-12
    return dict(samples=16, reference_dimension=3, min_diagnostic_effect_eigenvalue=eigen_min,
                max_reference_identity_residual=largest,
                meaning="cq labels and passive R only; sqrt(E) is not a surviving muon state")


def run():
    exact = exact_certificate()
    comparisons = {str(n): dict(pair=pair_coefficient(n), dalitz=dalitz_coefficient(n))
                   for n in (24, 48, 72)}
    differences = [abs(v["pair"]-v["dalitz"]) for v in comparisons.values()]
    assert max(differences) < 1e-14
    assert all(v["pair"] < -1/4500 for v in comparisons.values())
    special_k = (31-45*math.log(2))/16
    assert abs(dalitz_k(.5, 64)-special_k) < 1e-14
    total_b, positive_q, direct = direct_momentum_check()
    return {
        "round": 1056, "all_passed": True,
        "scope": "mu+ massless-daughter tree V-A; normalized decay branch; joint classical electron bin, total angular energy flux and passive R; not full coherent daughters, finite-time QFT or gravity instrument",
        "exact_certificate": exact,
        "quadrature_diagnostics_not_proof": comparisons,
        "special_energy_k_over_m": special_k,
        "fixed_direction_coefficient_magnitude": 25-36*math.log(2),
        "direct_momentum_diagnostics": direct,
        "reference_diagnostics": reference_diagnostic(total_b, positive_q),
        "historical_sha256": {f: digest(ROOT/f) for f in HISTORY},
        "new_cognitive_axioms": 0,
        "whole_roadmap_completed": False,
    }


def compare(current, saved, location=""):
    if isinstance(current, dict):
        assert current.keys() == saved.keys(), location
        for key in current:
            compare(current[key], saved[key], location+"/"+key)
    elif isinstance(current, list):
        assert len(current) == len(saved), location
        for i, (a, b) in enumerate(zip(current, saved)):
            compare(a, b, location+"/"+str(i))
    elif isinstance(current, float):
        assert math.isclose(current, saved, rel_tol=1e-11, abs_tol=2e-13), location
    else:
        assert current == saved, location


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="exclusively create first result")
    args = parser.parse_args()
    result = run()
    if args.write:
        with RESULT.open("x", encoding="utf-8", newline="\n") as output:
            json.dump(result, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        compare(result, json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"round": 1056, "all_passed": True,
                      "C_b_over_m_diagnostic": result["quadrature_diagnostics_not_proof"]["72"]["pair"],
                      "record_only_error_lower_exact": result["exact_certificate"]["record_only_TV_error_lower"],
                      "history_files": len(HISTORY), "mode": "create" if args.write else "read_only"}))


if __name__ == "__main__":
    main()
