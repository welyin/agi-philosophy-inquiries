"""Round 1062: a conditional collisionless mass-identification certificate.

Default: recompute and compare results.json, without writing anything.
--save-exclusive: create results.json once; never overwrite it.
NumPy quadrature checks finite examples; analysis.md supplies the proofs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
RESULT = HERE / "results.json"
HISTORY = [
    "archive_956_989/research_note_973.md",
    "archive_990_1008/research_note_1007.md",
    "archive_1009_1043/research_note_1034.md",
    "archive_1009_1043/research_note_1035.md",
    "archive_1009_1043/research_note_1036.md",
    "archive_1046_/_shared/neutrino_decoupling_adoption.md",
    "archive_1046_/_shared/neutrino_propagation_adoption.md",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quad(fn, lo: float, hi: float, order: int = 96) -> float:
    z, w = np.polynomial.legendre.leggauss(order)
    q = (hi + lo) / 2 + (hi - lo) * z / 2
    return float(np.dot(w, fn(q)) * (hi - lo) / 2)


def occupation(q, amplitude, dilation=1.0):
    return amplitude / dilation**4 * np.maximum(1 - (q / dilation)**2, 0)**2


def moments(a, mass, amplitude, dilation=1.0, order=96):
    f = lambda q: occupation(q, amplitude, dilation)
    rho = quad(lambda q: q*q * np.sqrt(q*q + (a*mass)**2) * f(q),
               0, dilation, order) / (2 * math.pi**2 * a**4)
    pressure = quad(lambda q: q**4 / np.sqrt(q*q + (a*mass)**2) * f(q),
                    0, dilation, order) / (6 * math.pi**2 * a**4)
    number = quad(lambda q: q*q * f(q), 0, dilation, order) / (2*math.pi**2*a**3)
    drho_da = quad(
        lambda q: q*q * (a*mass*mass / np.sqrt(q*q+(a*mass)**2)
                        - 4*np.sqrt(q*q+(a*mass)**2)/a) * f(q),
        0, dilation, order) / (2*math.pi**2*a**4)
    return np.array([rho, pressure, number, drho_da])


def beta_kernel(eps, mass_squared, weights):
    eps = np.asarray(eps)
    return sum(w * eps * np.sqrt(np.maximum(eps*eps-m2, 0))
               for m2, w in zip(mass_squared, weights))


def rational(x):
    return {"fraction": str(x), "decimal": float(x)}


def run():
    m2 = np.array([1., 4., 9.])
    shift = 3.
    mass = np.sqrt(m2)
    newmass = np.sqrt(m2 + shift)
    lam = newmass / mass
    amp = np.array([1/8, 1/4, 3/8])
    weights = np.array([1/2, 1/3, 1/6])
    scales = [1/100, 1/10, 1., 10., 100.]
    rows = []
    residuals = []
    for a in scales:
        original = np.array([moments(a, mass[i], amp[i], order=96) for i in range(3)])
        scaled = np.array([moments(a, newmass[i], amp[i], lam[i], order=113) for i in range(3)])
        rp = float(np.max(np.abs(scaled[:, :2]/original[:, :2] - 1)))
        nn = float(np.max(np.abs(lam*scaled[:, 2]/original[:, 2] - 1)))
        number_exact = amp * 4 / (105*math.pi**2*a**3)
        number_res = float(np.max(np.abs(original[:, 2]/number_exact-1)))
        continuity = float(np.max(np.abs(a*original[:, 3]+3*(original[:, 0]+original[:, 1]))
                                  / (3*(original[:, 0]+original[:, 1]))))
        residuals.extend([rp, nn, number_res, continuity])
        rows.append({"a": a, "rho_P_relative_error": rp,
                     "number_scaled_identity_relative_error": nn,
                     "number_exact_relative_error": number_res,
                     "continuity_relative_residual": continuity,
                     "original_rho_P_n": original[:, :3].tolist(),
                     "scaled_rho_P_n": scaled[:, :3].tolist()})
    assert max(residuals) < 2e-12

    # A source-level local-frame check with nonzero momentum flux; not a
    # non-vacuum Minkowski solution. Curved-spacetime transport is analytic.
    vectors = np.array([[.2, .3, -.1], [.1, -.4, .25], [-.3, .12, .23]])
    integrand_err = 0.
    velocity_err = 0.
    liouville_err = 0.
    for i, p in enumerate(vectors):
        p4 = np.r_[np.sqrt(m2[i]+p@p), p]
        pp4 = np.r_[np.sqrt(m2[i]+shift+lam[i]**2*(p@p)), lam[i]*p]
        measure = 1/p4[0]
        transformed_measure = lam[i]**3/pp4[0]
        f = .12
        A = np.outer(p4, p4)*f*measure
        B = np.outer(pp4, pp4)*f/lam[i]**4*transformed_measure
        integrand_err = max(integrand_err, float(np.max(np.abs(A-B))))
        velocity_err = max(velocity_err, float(np.max(np.abs(p4[1:]/p4[0]-pp4[1:]/pp4[0]))))
        # FLRW Christoffel jet at a=1, H=.3; off-shell derivatives test
        # L' f' = lambda^-3 L f, not numerical proof of PDE solutions.
        dx = np.array([.17, -.2, .3, .4])
        dp = np.array([-.2, .15, .09])
        L = p4@dx - .6*p4[0]*(p@dp)
        Lp = pp4@(dx/lam[i]**4) - .6*pp4[0]*((lam[i]*p)@(dp/lam[i]**5))
        liouville_err = max(liouville_err, abs(Lp-L/lam[i]**3))
    assert max(integrand_err, velocity_err, liouville_err) < 2e-14

    # Nonempty self-consistent FLRW seed. G=1 in fixed units, flat T^3,
    # a(0)=1, Lambda_cosm=0. The positive expanding branch is unique.
    def cosmic_integrand(a_values, scaled=False):
        out = []
        for a in a_values:
            rho = sum(moments(float(a), newmass[i] if scaled else mass[i], amp[i],
                              lam[i] if scaled else 1., order=72 if scaled else 64)[0]
                      for i in range(3))
            out.append(1/(a*math.sqrt(8*math.pi*rho/3)))
        return np.array(out)
    t1 = quad(lambda a: cosmic_integrand(a, False), 1, 2, 32)
    t2 = quad(lambda a: cosmic_integrand(a, True), 1, 2, 37)
    assert abs(t1-t2) < 2e-11

    # One fixed unitary mixing matrix. It is a theoretical certificate,
    # not the fitted PMNS matrix or observed mass splittings.
    row1 = np.sqrt(weights)
    row2 = np.array([math.sqrt(2/5), -math.sqrt(3/5), 0.])
    U = np.array([row1, row2, np.cross(row1, row2)])
    unitarity = float(np.max(np.abs(U@U.T-np.eye(3))))
    beta_old = float(weights@m2)
    beta_new = float(weights@(m2+shift))
    assert abs(beta_new-beta_old-shift) < 1e-14
    ur_bound = F(100)*F(3)*F(8)/(4*F(100)**3)
    phase_rows = []
    for p, t in [(100., 100.), (150., 60.), (300., 100.)]:
        E = np.sqrt(p*p+m2)
        Ep = np.sqrt(p*p+m2+shift)
        common = Ep[0]-E[0]
        V = U@np.diag(np.exp(-1j*t*(E-E[0])))@U.T
        Vp = U@np.diag(np.exp(-1j*t*(Ep-Ep[0])))@U.T
        err = float(np.linalg.norm(Vp-V, 2))
        bound = t*shift*8/(4*p**3)
        probability_error = float(np.max(np.abs(np.abs(Vp)**2-np.abs(V)**2)))
        assert err <= bound + 2e-11
        phase_rows.append({"p": p, "t": t, "reference_energy_shift": float(common),
                           "relative_unitary_difference": err, "analytic_bound": bound,
                           "largest_flavor_probability_difference": probability_error})

    lo, hi = F(5, 4), F(3, 2)
    window_old = quad(lambda e: beta_kernel(e, m2, weights), float(lo), float(hi), 64)
    window_new = quad(lambda e: beta_kernel(e, m2+shift, weights), float(lo), float(hi), 71)
    exact_old = (40*math.sqrt(5)-27)/384
    lower = F(61, 384)  # sqrt(5)>11/5, verified by rational squares.
    per_branch_error = F(1, 32)
    margin = lower-2*per_branch_error
    assert F(11, 5)**2 < 5
    assert abs(window_old-exact_old) < 2e-14 and window_new == 0
    assert window_old > float(lower) and margin > 0

    result = {
        "round": 1062, "status": "author_candidate_passed",
        "count": {"candidate_scientific_calibration_groups": 1,
                  "empirical_groups": 0, "new_cognitive_axioms": 0,
                  "baseline_round": 1061, "baseline_cumulative": 3837,
                  "baseline_stage_groups": 16, "global_count_finalized_by_mainline": False},
        "whole_roadmap_complete": False,
        "units": "hbar=c=1; fixed energy unit mu=1; G=1 only for the fixed FLRW diagnostic",
        "inputs": {"mass_squared": m2.tolist(), "mass_squared_shift": shift,
                   "new_mass_squared": (m2+shift).tolist(), "lambda": lam.tolist(),
                   "occupation_amplitudes": amp.tolist(), "electron_weights": weights.tolist(),
                   "U_real": U.tolist()},
        "rational_certificates": {"radial_number_integral": rational(F(8,105)),
                                  "m_beta_squared_old": rational(F(10,3)),
                                  "m_beta_squared_new": rational(F(19,3)),
                                  "UR_operator_error_bound": rational(ur_bound),
                                  "beta_window_rate_lower": rational(lower),
                                  "per_branch_rate_error_adopted_budget": rational(per_branch_error),
                                  "beta_rate_gap_after_two_errors_lower": rational(margin)},
        "flrw_moments": rows,
        "local_mass_shell_checks": {"stress_integrand_max_error": integrand_err,
                                    "velocity_max_error": velocity_err,
                                    "Liouville_scaling_max_error": liouville_err},
        "flrw_seed": {"space": "flat compact three-torus, fixed coordinate volume 1",
                      "a_initial": 1, "a_final": 2, "G": 1, "cosmological_constant": 0,
                      "elapsed_time_original": t1, "elapsed_time_scaled": t2,
                      "elapsed_time_difference": abs(t1-t2)},
        "mixing": {"unitarity_residual": unitarity,
                   "Delta_m_squared_preserved": np.array_equal(
                       (m2+shift)[:, None]-(m2+shift)[None, :], m2[:, None]-m2[None, :]),
                   "m_beta_squared_original": beta_old, "m_beta_squared_scaled": beta_new,
                   "exact_flavor_diagnostics": phase_rows},
        "beta_window": {"epsilon_over_mu": [float(lo), float(hi)],
                        "leading_reduced_bin_original": window_old,
                        "leading_reduced_bin_scaled": window_new,
                        "closed_form_original": exact_old,
                        "units": "kappa times mu^3 for the stated positive common-prefactor lower bound",
                        "actual_KATRIN_response_or_likelihood_used": False},
        "maximum_moment_relative_residual": max(residuals),
        "physical_boundaries": [
            "classical collisionless one-particle stress, not quantum stress noise",
            "free post-decoupling initial data, not fixed Standard Model thermal history",
            "laboratory beta is another declared preparation, not inserted into the equal-source spacetime",
            "UR flavor equivalence only at leading order; full exact processes need not agree",
            "particle numbers and full preparation resources are not held equal",
            "finite beta response/theory error budget is a condition, not an instrument construction",
            "finite numerical checks do not prove the spacetime or distribution quantifiers"],
        "history_sha256": {p: sha(RESEARCH/p) for p in HISTORY},
        "science_code_sha256": sha(Path(__file__).resolve()),
        "passed": True,
    }
    return result


def compare(old, new, path="root"):
    if isinstance(new, dict):
        assert old.keys() == new.keys(), path
        for key in new:
            compare(old[key], new[key], path+"."+key)
    elif isinstance(new, list):
        assert len(old) == len(new), path
        for i, (x, y) in enumerate(zip(old, new)):
            compare(x, y, path+f"[{i}]")
    elif isinstance(new, float):
        assert math.isclose(old, new, rel_tol=3e-12, abs_tol=3e-13), (path, old, new)
    else:
        assert old == new, (path, old, new)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save-exclusive", action="store_true")
    args = parser.parse_args()
    values = run()
    if args.save_exclusive:
        with RESULT.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(values, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write("\n")
        print("1062 results saved exclusively; author candidate passed.")
    else:
        saved = json.loads(RESULT.read_text(encoding="utf-8"))
        compare(saved, values)
        print("1062 read-only scientific comparison passed.")
