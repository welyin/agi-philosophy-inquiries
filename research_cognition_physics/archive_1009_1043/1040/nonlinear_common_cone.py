"""1040: constitutive linearization, full Maxwell symbol, and common-cone scope.

Default invocation recomputes and compares without writing. --write creates the
result exclusively; no cumulative research count is assigned by this branch.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PHYSICS = HERE.parents[1]
RESULT = HERE / "nonlinear_common_cone_results.json"
HISTORY = [
    "archive_301_341/research_note_305.md",
    "archive_301_341/research_note_312.md",
    "archive_301_341/research_note_317.md",
    "archive_342_369/research_note_347.md",
    "archive_342_369/research_note_352.md",
    "archive_342_369/research_note_369.md",
    "archive_370_428/research_note_382.md",
    "archive_370_428/research_note_383.md",
    "archive_370_428/research_note_384.md",
    "archive_370_428/research_note_386.md",
    "archive_370_428/research_note_425.md",
    "archive_467_530/research_note_522.md",
    "archive_467_530/research_note_523.md",
    "archive_702_741/research_note_720.md",
    "archive_923_934/research_note_923.md",
    "archive_923_934/research_note_930.md",
    "archive_1009_/1035/bridge_ledger_v0_2.json",
    "archive_1009_/research_note_1038.md",
]


def cross_matrix(v):
    x, y, z = v
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)


def invariants(E, B):
    return (E @ E - B @ B) / 2, E @ B


def bi_jet(S, P, T):
    R = np.sqrt(1 - 2*S/T - P*P/(T*T))
    L = T*(1-R)
    LS, LP = 1/R, P/(T*R)
    H = np.array([[1/(T*R**3), P/(T*T*R**3)],
                  [P/(T*T*R**3), 1/(T*R)+P*P/(T**3*R**3)]])
    return L, LS, LP, H


def constitutive(E, B, LS, LP, hessian):
    """Return d(D,H)/d(E,B), using D=L_E and H=-L_B."""
    I, O = np.eye(3), np.zeros((3, 3))
    gradients = np.column_stack([np.r_[E, -B], np.r_[B, E]])
    H6 = LS*np.diag([1, 1, 1, -1, -1, -1])
    H6 += LP*np.block([[O, I], [I, O]])
    H6 += gradients @ hessian @ gradients.T
    return np.diag([1, 1, 1, -1, -1, -1]) @ H6


def canonical_hessian(C):
    eps, xi = C[:3, :3], C[:3, 3:]
    eta, mu_inv = C[3:, :3], C[3:, 3:]
    assert np.allclose(eta, -xi.T, atol=1e-13)
    inv = np.linalg.inv(eps)
    return np.block([[inv, -inv@xi], [-xi.T@inv, mu_inv+xi.T@inv@xi]])


def maxwell_matrix(C, k):
    """State (delta D, delta B); includes both Gauss constraint directions."""
    K, O = cross_matrix(k), np.zeros((3, 3))
    J = np.block([[O, -K], [K, O]])
    HC = canonical_hessian(C)
    N = J @ HC
    # Similar symmetric matrix, for numerically stable frequencies.
    R = np.linalg.cholesky(HC)
    eigenvalues = np.linalg.eigvalsh(R.T @ J @ R)
    return N, HC, eigenvalues


def field_tensor(E, B):
    F4 = np.zeros((4, 4))
    F4[0, 1:], F4[1:, 0] = -E, E
    F4[1:, 1:] = -cross_matrix(B)
    return F4


def stress_from_jet(E, B, L, LS, LP):
    eta = np.diag([-1., 1., 1., 1.])
    F4 = field_tensor(E, B)
    _, P = invariants(E, B)
    return LS*(F4 @ eta @ F4.T) + eta*(L-P*LP)


def lagrangian_density_inverse_metric(invg, E, B, T):
    F4 = field_tensor(E, B)
    volume = np.sqrt(-1/np.linalg.det(invg))
    Fup = invg @ F4 @ invg.T
    S = -np.sum(F4*Fup)/4
    # Covariant F is held fixed: sqrt(-g) P = E dot B in this chart.
    P = E @ B / volume
    L = T*(1-np.sqrt(1-2*S/T-P*P/(T*T)))
    return volume*L


def rational_matrix_witness():
    cases = []
    # Exact all-Hessian necessity on different nonzero invariant backgrounds.
    for e, b in [(F(1, 2), F(2, 3)), (F(0), F(1)), (F(1), F(0)),
                 (F(3, 4), F(-1, 5))]:
        J = [[-b, e], [e, b]]
        H = [[F(2, 7), F(-3, 11)], [F(-3, 11), F(5, 13)]]
        mul = lambda A, B: [[sum(A[i][k]*B[k][j] for k in range(2))
                            for j in range(2)] for i in range(2)]
        JJ = mul(J, J)
        r2 = e*e+b*b
        assert JJ == [[r2, F(0)], [F(0), r2]]
        M = mul(mul(J, H), J)
        recovered = [[x/(r2*r2) for x in row] for row in mul(mul(J, M), J)]
        assert recovered == H
        cases.append({"e": str(e), "b": str(b), "J_squared": str(r2),
                      "hessian_exactly_recovered": True})
    return cases


def build_results():
    rng = np.random.default_rng(1040008)
    max_null_screen = max_constitutive = max_pure_spectrum = 0.
    max_general_spectrum = max_gauss = max_stress = max_energy = 0.
    min_health = 100.
    # Full 6x6 constitutive and first-order equations, rather than only a scalar dispersion formula.
    for _ in range(32):
        e, b = rng.uniform(-.5, .5, 2)
        if e*e+b*b < .01:
            b += .3
        H = rng.normal(size=(2, 2)); H = (H+H.T)/2
        LS, LP = rng.uniform(.4, 2.), rng.uniform(-.3, .3)
        E, B = np.array([0., 0., e]), np.array([0., 0., b])
        C = constitutive(E, B, LS, LP, H)
        M = np.zeros((2, 2))
        k = np.array([1., 0., 0.])
        for j, ep in enumerate([np.array([0., 1., 0.]), np.array([0., 0., 1.])]):
            bp = np.cross(k, ep)
            dp, hp = np.split(C @ np.r_[ep, bp], 2)
            M[:, j] = [hp[2]-dp[1], hp[1]+dp[2]]
        J = np.array([[-b, e], [e, b]])
        predicted = np.diag([-1., 1.]) @ J @ H @ J
        max_null_screen = max(max_null_screen, float(np.max(abs(M-predicted))))
        assert np.allclose(np.linalg.norm(M, 2), (e*e+b*b)*np.linalg.norm(H, 2))
    pure_cases = []
    for T in [1., 2., 4.]:
        for Bmag in [.25, .5, 1., 2.]:
            E, B = np.zeros(3), np.array([0., 0., Bmag])
            S, P = invariants(E, B); L, LS, LP, Hess = bi_jet(S, P, T)
            C = constitutive(E, B, LS, LP, Hess)
            R = np.sqrt(1+Bmag*Bmag/T)
            target = np.diag([1/R, 1/R, R, 1/R, 1/R, 1/R**3])
            max_constitutive = max(max_constitutive, float(np.max(abs(C-target))))
            for k in [np.array([1., 0., 0.]), np.array([0., 0., 1.]), np.array([.6, 0., .8])]:
                N, HC, w = maxwell_matrix(C, k)
                min_health = min(min_health, float(np.linalg.eigvalsh(HC)[0]))
                omega = np.sqrt(k[2]**2+(k[0]**2+k[1]**2)/(R*R))
                expected = np.array([-omega, -omega, 0., 0., omega, omega])
                max_pure_spectrum = max(max_pure_spectrum, float(np.max(abs(w-expected))))
                gauss = np.block([[k.reshape(1, 3), np.zeros((1, 3))],
                                  [np.zeros((1, 3)), k.reshape(1, 3)]])
                max_gauss = max(max_gauss, float(np.max(abs(gauss@N))))
            stress = stress_from_jet(E, B, L, LS, LP)
            u = T*(R-1); pt = T*(1-1/R)
            max_stress = max(max_stress, float(np.max(abs(stress-np.diag([u, pt, pt, -u])))))
            max_energy = max(max_energy, abs(float(stress[0, 0]+L)))
            pure_cases.append({"T": T, "B": Bmag, "R": float(R), "v_perp_squared": 1/(R*R),
                               "energy_density": float(u), "transverse_pressure": float(pt),
                               "parallel_pressure": float(-u)})
    # General weak fields check the double physical characteristic branch, with our q=(-omega,k).
    for _ in range(24):
        E, B = rng.uniform(-.2, .2, 3), rng.uniform(-.4, .4, 3)
        T = rng.uniform(1., 3.)
        k = rng.normal(size=3); k /= np.linalg.norm(k)
        S, P = invariants(E, B); L, LS, LP, Hess = bi_jet(S, P, T)
        C = constitutive(E, B, LS, LP, Hess)
        N, HC, w = maxwell_matrix(C, k)
        min_health = min(min_health, float(np.linalg.eigvalsh(HC)[0]))
        lam = 1/(T-2*S); aa = 1+lam*(E@E); drift = lam*(k@np.cross(E, B))
        rhs = aa*(k@k)-lam*(np.cross(k, E)@np.cross(k, E)+np.cross(k, B)@np.cross(k, B))
        roots = np.sort(np.roots([aa, -2*drift, -rhs]))
        expected = np.sort(np.r_[roots, roots, 0., 0.])
        max_general_spectrum = max(max_general_spectrum, float(np.max(abs(w-expected))))
        # Metric variation of the same L, not an independently fitted stress definition.
        stress = stress_from_jet(E, B, L, LS, LP)
        eta = np.diag([-1., 1., 1., 1.])
        for i in range(4):
            for j in range(i, 4):
                direction = np.zeros((4, 4)); direction[i, j] = 1.; direction[j, i] = 1.
                derivative = np.imag(lagrangian_density_inverse_metric(eta+1j*1e-25*direction, E, B, T))/1e-25
                predicted = -.5*np.sum(stress*direction)
                max_stress = max(max_stress, abs(float(derivative-predicted)))
        D = LS*E+LP*B
        HB = np.sqrt(T*T+T*(D@D+B@B)+np.cross(D, B)@np.cross(D, B))-T
        max_energy = max(max_energy, abs(float(E@D-L-HB)))
    # Finite common-cone tolerance: exact rational squared-speed information.
    bounds = []
    for B2, delta in [(F(1), F(1, 100)), (F(1, 4), F(1, 20)), (F(4), F(1, 5))]:
        Tmin = B2*(1-delta)/delta
        gap = B2/(Tmin+B2)
        assert gap == delta
        bounds.append({"B_squared": str(B2), "squared_speed_budget": str(delta),
                       "necessary_and_sufficient_T_lower_bound": str(Tmin), "boundary_gap": str(gap)})
    # A weak nonlinear anisotropy control: equality of photon polarizations is also an extra requirement.
    a, b, B2 = F(1, 100), F(7, 400), F(1, 4)
    ls = 1-a*B2
    vy2, vz2 = (ls-2*a*B2)/ls, ls/(ls+2*b*B2)
    assert vy2 != vz2 and 0 < vy2 < 1 and 0 < vz2 < 1
    # Two causal BI actions: same zero-field normalization, finite strictly different v^2.
    difference = F(4, 5)-F(1, 2)
    per_description_error = F(1, 100)
    robust = difference-2*per_description_error
    assert robust == F(7, 25)
    thresholds = [max_null_screen, max_constitutive, max_pure_spectrum, max_general_spectrum,
                  max_gauss, max_stress, max_energy]
    assert max(thresholds) < 2e-11 and min_health > 0
    return {
        "round": 1040, "date": "2026-10-08", "scientific_baseline": 1038,
        "new_calibration_groups": 1, "new_adopted_cognitive_axioms": 0,
        "goal_completed": False, "all_checks_passed": True,
        "scope": {"Lorentz_3plus1_and_local_L_SP_adopted": True,
                  "no_birefringence_is_new_adopted_cognitive_axiom": False,
                  "full_quantum_or_UV_completion_proved": False,
                  "background_is_flat_Einstein_solution": False,
                  "actual_instrument_or_EFT_error_certified": False,
                  "general_Boillat_classification_numerically_proved": False},
        "rational_null_screen_cases": rational_matrix_witness(),
        "general_null_screen_hessian_cases": 32,
        "max_null_screen_matrix_residual": max_null_screen,
        "pure_BI_cases": pure_cases,
        "max_constitutive_residual": max_constitutive,
        "max_pure_full_Maxwell_spectrum_residual": max_pure_spectrum,
        "general_BI_full_Maxwell_cases": 24,
        "max_general_full_Maxwell_spectrum_residual": max_general_spectrum,
        "max_Gauss_preservation_residual": max_gauss,
        "minimum_tested_canonical_energy_hessian_eigenvalue": min_health,
        "max_same_action_stress_residual": max_stress,
        "max_Legendre_energy_residual": max_energy,
        "finite_squared_speed_bounds": bounds,
        "two_BI_models": {"B_squared": "1", "T_values": ["1", "4"],
                          "v_squared": ["1/2", "4/5"], "difference": str(difference),
                          "conditional_each_error": str(per_description_error),
                          "conditional_remaining_difference": str(robust)},
        "positive_birefringent_quartic_control": {"L": "S+a*S^2+b*P^2", "a": str(a), "b": str(b),
                                                 "B_squared": str(B2), "vy_squared": str(vy2),
                                                 "vz_squared": str(vz2), "both_subluminal": True},
        "historical_sha256": {p: hashlib.sha256((PHYSICS/p).read_bytes()).hexdigest() for p in HISTORY},
    }


def compare(a, b, path="root"):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a: compare(a[k], b[k], path+"."+k)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)): compare(x, y, path+f"[{i}]")
    elif isinstance(a, float):
        assert abs(a-b) <= 2e-10+1e-10*abs(a), (path, a, b)
    else:
        assert a == b, (path, a, b)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); result = build_results()
    if args.write:
        with RESULT.open("x", encoding="utf8") as f: json.dump(result, f, ensure_ascii=False, indent=2); f.write("\n")
    else:
        compare(result, json.loads(RESULT.read_text(encoding="utf8")))
    print(json.dumps({"round": 1040, "all_checks_passed": True, "mode": "exclusive_write" if args.write else "read_only_compare",
                      "max_full_symbol_residual": max(result["max_pure_full_Maxwell_spectrum_residual"], result["max_general_full_Maxwell_spectrum_residual"]),
                      "max_stress_residual": result["max_same_action_stress_residual"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
