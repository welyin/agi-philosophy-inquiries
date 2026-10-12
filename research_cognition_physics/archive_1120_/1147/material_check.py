"""1147 material-only diagnostics. Default replays and compares saved JSON.
No spectral cutoff is used to estimate the gap, heat trace, or gate error.
The two-mode CAR calculation tests an exact subalgebra, not a full-sea truncation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EPSILON = 0.01
MASSES = (1.0, 1.001)


def j0(x):
    return np.sin(x) / x


def j1(x):
    return np.sin(x) / x**2 - np.cos(x) / x


def radial(mass, order=220):
    lo, hi = 2.0, math.pi
    for _ in range(90):
        q = (lo + hi) / 2.0
        energy = math.hypot(q, mass)
        f = j0(q) - q / (energy + mass) * j1(q)
        if f > 0.0:
            lo = q
        else:
            hi = q
    q = (lo + hi) / 2.0
    energy = math.hypot(q, mass)
    b = q / (energy + mass)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    r = (nodes + 1.0) / 2.0
    w = weights / 2.0
    upper, lower = j0(q*r), b*j1(q*r)
    den = float(np.dot(w, r*r*(upper*upper + lower*lower)))
    scalar = float(np.dot(w, r*r*(upper*upper - lower*lower))/den)
    closed = (energy/2.0 + mass*(energy-1.0)) / (
        energy*(energy-1.0) + mass/2.0
    )
    return {
        "mR": mass,
        "q": q,
        "ER": energy,
        "boundary_residual": float(abs(j0(q)-b*j1(q))),
        "scalar_radial_integral": scalar,
        "scalar_closed_formula": closed,
        "integral_closed_difference": abs(scalar-closed),
    }


def group_radial():
    rows = []
    for mass in MASSES:
        row = radial(mass)
        d = 1e-4
        derivative = (radial(mass+d)["ER"]-radial(mass-d)["ER"])/(2*d)
        row["independent_mass_derivative"] = derivative
        row["derivative_scalar_difference"] = abs(
            derivative-row["scalar_radial_integral"]
        )
        row["j1_numerator_before_exact_heat_normalization"] = (
            math.exp(-2.0*EPSILON*row["ER"]**2)
            * row["scalar_closed_formula"]
        )
        assert row["boundary_residual"] < 2e-15
        assert row["integral_closed_difference"] < 2e-13
        assert row["derivative_scalar_difference"] < 2e-9
        assert 0.0 < row["scalar_closed_formula"] < 1.0
        rows.append(row)
    return {
        "name": "MIT_ground_root_scalar_density_and_mass_derivative",
        "pass": True,
        "fixed_inputs": {"R": 1, "m0": 1, "m1": "1001/1000", "epsilon": "1/100"},
        "rows": rows,
        "strict_scalar_positivity_basis": "monotonic j1(x)/j0(x) on (0,pi), not floating samples",
        "gap": "Delta_D=min(g0,g1,2*e_min)>0, exact spectral definition",
        "heat_trace": "C_epsilon=Tr(exp(-2*epsilon*h_m1^2))<infinity, exact spectral definition",
        "gap_or_heat_trace_truncation_used": False,
        "full_source_norm": "2*C_epsilon/(1+2*C_epsilon)<1",
        "j1": "exp(-2*epsilon*E1^2)*s(m1)/(1+2*C_epsilon)",
        "gate_certificate_claimed": False,
    }


def opnorm(x):
    return float(np.linalg.norm(x, 2))


def group_clock_car():
    energies = [radial(m)["ER"] for m in MASSES]
    e0, e1 = energies
    nu = e1-e0
    ann = np.array([[0.0, 1.0], [0.0, 0.0]], complex)
    I = np.eye(2)
    Z = np.diag([1.0, -1.0])
    b0, b1 = np.kron(ann, I), np.kron(Z, ann)
    n0, n1 = b0.conj().T@b0, b1.conj().T@b1
    vacuum = np.array([1.0, 0.0, 0.0, 0.0], complex)
    u0, u1 = b0.conj().T@vacuum, b1.conj().T@vacuum
    W = np.column_stack((u0, u1))
    P = W@W.conj().T
    X = np.outer(u0, u1.conj()) + np.outer(u1, u0.conj())
    Y = -1j*np.outer(u0, u1.conj()) + 1j*np.outer(u1, u0.conj())
    total_charge = n0+n1
    H = e0*n0 + e1*n1
    effect = (P-Y)/2
    psi = (u0+u1)/math.sqrt(2.0)
    checks = {
        "CAR_different_modes": opnorm(b0@b1.conj().T+b1.conj().T@b0),
        "CAR_same_mode": opnorm(b0@b0.conj().T+b0.conj().T@b0-np.eye(4)),
        "code_isometry": opnorm(W.conj().T@W-np.eye(2)),
        "full_energy_compression": opnorm(
            W.conj().T@H@W-np.diag([e0, e1])
        ),
        "X_total_charge_commutator": opnorm(total_charge@X-X@total_charge),
        "X_flavor_charge_commutator": opnorm(n1@X-X@n1),
    }
    assert all(checks[k] < 2e-15 for k in checks if k != "X_flavor_charge_commutator")
    assert abs(checks["X_flavor_charge_commutator"]-1.0) < 2e-15
    assert 0.0 < nu < MASSES[1]-MASSES[0]
    probabilities = []
    for t in (1.0, 1000.0):
        evolved = np.exp(-1j*np.diag(H)*t)*psi
        observed = float(np.vdot(evolved, effect@evolved).real)
        expected = (1.0+math.sin(nu*t))/2.0
        wrong_mass_only_clock = (1.0+math.sin((MASSES[1]-MASSES[0])*t))/2.0
        assert abs(observed-expected) < 5e-13
        probabilities.append({
            "t": t,
            "actual_bare_CAR_effect_probability": observed,
            "exact_free_clock_formula": expected,
            "formula_difference": abs(observed-expected),
            "incorrect_mass_difference_only_probability": wrong_mass_only_clock,
        })
    return {
        "name": "same_positive_bound_modes_actual_free_clock_and_flavor_readout_algebra",
        "pass": True,
        "E0": e0,
        "E1": e1,
        "free_clock_frequency": nu,
        "matrix_checks": checks,
        "probabilities": probabilities,
        "shared_negative_mode_sea": "fixed spectator in this exact positive-mode subalgebra",
        "normal_full_J_not_replaced_by_two_mode_operator": True,
        "material_pair_processes_simulated_here": False,
        "source_and_interacting_clock_phase_unchanged_claim": False,
        "flavor_U1_is_a_dynamical_symmetry_not_an_imposed_SSR": True,
        "reading_requires_declared_preparation_and_readout_permissions": True,
    }


def run():
    proof = HERE/"material_bound.md"
    source = HERE.parent/"_shared"/"after1146_bound_dirac_feasibility.md"
    groups = [group_radial(), group_clock_car()]
    return {
        "round": 1147,
        "check_groups": 2,
        "all_passed": all(g["pass"] for g in groups),
        "scope": "actual MIT material formulas and free-clock CAR diagnostics; no full gate certificate",
        "source_sha256": {
            "material_bound.md": hashlib.sha256(proof.read_bytes()).hexdigest(),
            "../_shared/after1146_bound_dirac_feasibility.md": hashlib.sha256(source.read_bytes()).hexdigest(),
        },
        "groups": groups,
    }


def compare(saved, actual, location="$"):
    if isinstance(saved, dict):
        if not isinstance(actual, dict) or saved.keys() != actual.keys():
            raise AssertionError(f"{location}: dictionary keys differ")
        for key in saved:
            compare(saved[key], actual[key], location+"."+key)
    elif isinstance(saved, list):
        if not isinstance(actual, list) or len(saved) != len(actual):
            raise AssertionError(f"{location}: list shape differs")
        for i, (left, right) in enumerate(zip(saved, actual)):
            compare(left, right, f"{location}[{i}]")
    elif isinstance(saved, float):
        if not math.isclose(saved, actual, rel_tol=2e-12, abs_tol=5e-13):
            raise AssertionError(f"{location}: numeric diagnostic differs")
    elif saved != actual:
        raise AssertionError(f"{location}: saved value differs")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true",
                        help="Create saved result; refuse to replace a different existing result.")
    args = parser.parse_args()
    actual = run()
    destination = HERE/"material_results.json"
    if args.write and not destination.exists():
        destination.write_text(json.dumps(actual, ensure_ascii=False, indent=2)+"\n",
                               encoding="utf-8")
    if not destination.exists():
        raise FileNotFoundError("Saved results absent. --write is required to create them.")
    saved = json.loads(destination.read_text(encoding="utf-8"))
    compare(saved, actual)
    print(json.dumps({
        "round": 1147,
        "material_check_groups": actual["check_groups"],
        "all_passed": actual["all_passed"],
        "saved_results_match": True,
        "spectral_gap_or_heat_trace_cutoff_used": False,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
