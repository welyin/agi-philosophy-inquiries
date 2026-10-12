"""1121: full CAR momentum defect and independent radial integration.

Default use recomputes and compares results.json without writing.  --write
creates the initial result artifact.  Finite p matrices below check algebra;
they are not a truncation of the physical positive-energy Dirac Hamiltonian.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
TOL = 2e-11


def norm(a):
    return float(np.linalg.norm(a, 2))


def car_annihilators(n):
    out = []
    for j in range(n):
        a = np.zeros((1 << n, 1 << n), dtype=complex)
        for state in range(1 << n):
            if state & (1 << j):
                parity = (state & ((1 << j) - 1)).bit_count()
                a[state ^ (1 << j), state] = (-1) ** parity
        out.append(a)
    return out


CAR = car_annihilators(4)
I16 = np.eye(16, dtype=complex)
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
PARITY = np.diag([(-1) ** k.bit_count() for k in range(16)])


def q(a):
    return sum(
        (a[j, k] * CAR[j].conj().T @ CAR[k]
         for j in range(4) for k in range(4)),
        np.zeros((16, 16), dtype=complex),
    )


def phase(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(1j * t * values)) @ vectors.conj().T


def assert_close(actual, expected, label, atol=TOL):
    err = norm(np.atleast_2d(np.asarray(actual) - np.asarray(expected)))
    assert err <= atol, (label, err)
    return err


def mode_case(sigma, theta, lower_diagonal=2.0):
    # The lower diagonal cancels from V^* p V-p.  This verifies precisely
    # the claimed finite-rank identity without calling p a physical H.
    p4 = np.kron(np.array([[0.0, sigma], [sigma, lower_diagonal]]), I2)
    a4 = np.zeros((4, 4), dtype=complex)
    a4[:2, :2] = X
    u = phase(q(a4), -theta)
    e = math.cos(theta) * I2 - 1j * math.sin(theta) * X
    v = np.eye(4, dtype=complex)
    v[:2, :2] = e
    b_formula = np.block([
        [np.zeros((2, 2)), sigma * (e.conj().T - I2)],
        [sigma * (e - I2), np.zeros((2, 2))],
    ])
    b_direct = v.conj().T @ p4 @ v - p4
    defect_direct = u.conj().T @ q(p4) @ u - q(p4)
    defect = q(b_formula)
    s = abs(math.sin(theta / 2))
    eigen_scale = 2 * sigma * s
    one_expected = eigen_scale * np.array([-1, -1, 1, 1])
    full_expected = eigen_scale * np.array(
        [-2] + [-1] * 4 + [0] * 6 + [1] * 4 + [2]
    )
    residuals = {
        "finite_rank_identity": assert_close(b_direct, b_formula, "B"),
        "full_car_conjugation": assert_close(defect_direct, defect, "qB"),
        "one_particle_spectrum": assert_close(
            np.linalg.eigvalsh(b_formula), one_expected, "one spectrum"),
        "full_car_spectrum": assert_close(
            np.linalg.eigvalsh(defect), full_expected, "full spectrum"),
        "norm_factor_four": assert_close(norm(defect), 4 * sigma * s, "norm"),
        "evenness": assert_close(defect @ PARITY, PARITY @ defect, "even"),
        "reference_identity": assert_close(
            norm(np.kron(defect, np.eye(3))), norm(defect), "reference"),
    }
    even_indices = [k for k in range(16) if k.bit_count() % 2 == 0]
    even_spectrum = np.linalg.eigvalsh(defect[np.ix_(even_indices, even_indices)])
    assert_close(even_spectrum[[0, -1]], [-2 * eigen_scale, 2 * eigen_scale],
                 "even extrema")
    return {
        "sigma": sigma, "theta": theta, "one_particle_norm": norm(b_formula),
        "full_car_norm": norm(defect), "expected_full_norm": 4 * sigma * s,
        "max_residual": max(residuals.values()),
    }, (p4, a4, u, b_formula, defect)


def run():
    car_residual = 0.0
    for j in range(4):
        for k in range(4):
            car_residual = max(car_residual, assert_close(
                CAR[j] @ CAR[k].conj().T + CAR[k].conj().T @ CAR[j],
                I16 if j == k else 0 * I16, "CAR mixed"))
            car_residual = max(car_residual, assert_close(
                CAR[j] @ CAR[k] + CAR[k] @ CAR[j], 0 * I16, "CAR aa"))

    cases = []
    for sigma in (0.0, 1.0, math.sqrt(11 / 3)):
        for theta in (0.0, math.pi / 3, math.pi / 2, math.pi, 2 * math.pi):
            case, _ = mode_case(sigma, theta)
            cases.append(case)

    _, (p4, a4, u, b, defect) = mode_case(1.0, math.pi / 2)
    p_many = q(p4)
    derivative_exact = 1j * (p_many @ u - u @ p_many)
    derivative_errors = []
    for step in (1e-3, 5e-4):
        tp = phase(p_many, step)
        tm = phase(p_many, -step)
        derivative = (tp @ u @ tp.conj().T - tm @ u @ tm.conj().T) / (2 * step)
        derivative_errors.append(norm(derivative - derivative_exact))
    assert derivative_errors[1] < 0.27 * derivative_errors[0]
    derivative_to_defect = assert_close(-1j * u.conj().T @ derivative_exact,
                                       defect, "translation derivative")

    # Normal ordering in a nontrivial pure quasifree covariance.  The filled
    # Slater vector is built independently, not just by applying Tr(SB).
    p2 = np.array([[0.0, 1.0], [1.0, 2.0]])
    vals, vecs = np.linalg.eigh(p2)
    negative = vecs[:, vals < 0]
    s0 = negative @ negative.conj().T
    covariance = np.kron(s0, I2)
    scalar = float(np.trace(covariance @ b).real)
    assert_close(scalar, math.sqrt(2), "normal-order scalar")
    state = np.zeros(16, dtype=complex)
    state[0] = 1.0
    for flavor in np.eye(2):
        occupied = np.kron(negative[:, 0], flavor)
        creation = sum(occupied[k] * CAR[k].conj().T for k in range(4))
        state = creation @ state
    assert_close(np.vdot(state, state), 1.0, "Slater normalization")
    vacuum_expectation = float(np.vdot(state, defect @ state).real)
    assert_close(vacuum_expectation, scalar, "Slater normal-order scalar")
    normal_ordered = defect - scalar * I16
    assert_close(np.vdot(state, normal_ordered @ state), 0.0, "normal order")
    gate_scalar = float(np.trace(covariance @ a4).real)
    u_normal = phase(q(a4) - gate_scalar * I16, -math.pi / 2)
    normal_gate_defect = u_normal.conj().T @ p_many @ u_normal - p_many
    normal_gate_residual = assert_close(normal_gate_defect, defect, "gate phase")

    # Exact beta-polynomial integrals and an independent Gauss-Legendre rule.
    i_norm = sum(Fraction((-1) ** k * math.comb(4, k), 3 + 2 * k)
                 for k in range(5))
    i_grad = 16 * (Fraction(1, 5) - 2 * Fraction(1, 7) + Fraction(1, 9))
    assert i_norm == Fraction(128, 3465)
    assert i_grad == Fraction(128, 315)
    assert i_grad / i_norm == 11
    nodes, weights = leggauss(64)
    t, w = (nodes + 1) / 2, weights / 2
    quadrature_norm = float(np.dot(w, t ** 2 * (1 - t ** 2) ** 4))
    quadrature_grad = float(np.dot(w, 16 * t ** 4 * (1 - t ** 2) ** 2))
    integral_residuals = [
        assert_close(quadrature_norm, float(i_norm), "radial norm"),
        assert_close(quadrature_grad, float(i_grad), "radial gradient"),
    ]
    radial_cases = []
    for radius in (1.0, 2.0):
        c_squared = 3465 / (512 * math.pi * radius ** 3)
        normalized_mass = c_squared * 4 * math.pi * radius ** 3 * quadrature_norm
        grad_sum = c_squared * 4 * math.pi * radius * quadrature_grad
        sigma = math.sqrt(grad_sum / 3)
        assert_close(normalized_mass, 1.0, "radial mass")
        assert_close(grad_sum, 11 / radius ** 2, "radial sigma")
        theta = math.pi / 2
        full = 4 * sigma * abs(math.sin(theta / 2))
        lower = 4 * math.pi / (math.sqrt(3) * radius) * abs(math.sin(theta / 2))
        assert full >= lower
        radial_cases.append({
            "radius": radius, "normalized_mass": normalized_mass,
            "gradient_sum": grad_sum, "sigma_j_squared": grad_sum / 3,
            "theta": theta, "full_car_norm_each_direction": full,
            "at_least_one_direction_lower_bound": lower,
        })

    return {
        "round": 1121,
        "status": "PASS",
        "new_cognitive_axioms": 0,
        "scientific_calibration": 0,
        "scope": "Finite CAR algebra checks and radial quadrature, not a physical device simulation.",
        "groups": {
            "car_and_full_spectrum": {
                "modes": 4, "hilbert_dimension": 16,
                "car_residual": car_residual, "cases": cases,
                "full_spectrum_multiplicities": [1, 4, 6, 4, 1],
                "reference_dimension_checked": 3,
            },
            "translation_derivative": {
                "finite_difference_steps": [1e-3, 5e-4],
                "finite_difference_errors": derivative_errors,
                "derivative_to_defect_residual": derivative_to_defect,
                "interpretation": "Finite algebra identity only; analytic proof uses the actual physical translation automorphism.",
            },
            "normal_ordering": {
                "scalar_trace": scalar,
                "independent_filled_slater_expectation": vacuum_expectation,
                "normal_ordered_norm": norm(normal_ordered),
                "actual_defect_norm": norm(defect),
                "normal_ordered_gate_defect_residual": normal_gate_residual,
                "interpretation": "Nonzero scalar must be retained in the defect; finite covariance is not a Dirac vacuum discretization.",
            },
            "radial_mode": {
                "profile": "C_R(1-r^2/R^2)_+^2 chi; H0^1 intersect H^2, not C-infinity",
                "norm_integral_exact": str(i_norm),
                "gradient_integral_exact": str(i_grad),
                "gradient_over_norm_exact": str(i_grad / i_norm),
                "quadrature_order": 64,
                "max_integral_residual": max(integral_residuals),
                "cases": radial_cases,
            },
        },
        "analytic_only_claims": [
            "Positive-energy Dirac representation and actual unbounded momentum domains.",
            "All smooth compact modes, all radii and angles, arbitrary passive references.",
            "Ball Poincare lower bound and H0^1 density extension.",
            "Conditional controller balance and all-input exact-unitary obstruction with finite moments.",
        ],
        "not_certified": [
            "Actual finite-duration local controller or point Hamiltonian implementation.",
            "Original CO endpoint, physical clock and complete source-family identification.",
            "Attainment within a fixed physical energy budget or old preparation menu.",
            "Unbounded momentum control from ordinary channel distance.",
        ],
    }


def compare(actual, expected, path="root"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], path + "." + key)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i, (left, right) in enumerate(zip(actual, expected)):
            compare(left, right, path + f"[{i}]")
    elif isinstance(actual, float):
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-11), path
    else:
        assert actual == expected, path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the initial results artifact")
    args = parser.parse_args()
    results = run()
    destination = HERE / "results.json"
    if args.write:
        destination.write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        action = "initial results written"
    else:
        compare(results, json.loads(destination.read_text(encoding="utf-8")))
        action = "read-only comparison passed"
    print(f"1121 PASS: 4 groups; 15 CAR cases; {action}; new axioms 0; scientific calibration 0")


if __name__ == "__main__":
    main()
