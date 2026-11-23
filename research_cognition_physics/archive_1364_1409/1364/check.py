"""Exact finite one-down-singlet weak-moment test for round 1364."""

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np


RESULT = Path(__file__).with_name("results.json")
CYCLE = ((0, 1), (1, 2), (2, 0))


def cycle_entries(matrix):
    return np.array([matrix[i, j] for i, j in CYCLE])


def polynomial(left, deficit):
    # Each edge is linear in the candidate heavy squared mass.
    coefficients = np.array([1 + 0j])
    for i, j in CYCLE:
        coefficients = np.polynomial.polynomial.polymul(
            coefficients, [left[i, j], deficit[i, j]]
        )
    return coefficients.imag


def factor_real_parent(hermitian):
    # Input is in a row convention where its upper principal block is real.
    upper = hermitian[:3, :3]
    assert np.max(np.abs(upper.imag)) < 1e-10
    ordinary = np.linalg.cholesky(upper.real)
    coupling = np.linalg.solve(ordinary, hermitian[:3, 3]).conj()
    heavy = np.sqrt(hermitian[3, 3].real - np.vdot(coupling, coupling).real)
    parent = np.zeros((4, 4), dtype=complex)
    parent[:3, :3] = ordinary
    parent[3, :3] = coupling
    parent[3, 3] = heavy
    return parent


def compute():
    ordinary = np.array([[1, 0, 0], [1, 2, 0], [1, 1, 3]], dtype=float)
    expectation = 1 + 1j
    yukawa = np.array([0.5, 0.5, 1.0])
    conjugate_yukawa = np.array([0.5, -0.5, 0.0])
    coupling = yukawa * expectation + conjugate_yukawa * expectation.conjugate()
    parent = np.zeros((4, 4), dtype=complex)
    parent[:3, :3] = ordinary
    parent[3, :3] = coupling
    parent[3, 3] = 4
    hermitian = parent @ parent.conj().T
    squared, eigenvectors = np.linalg.eigh(hermitian)
    light = eigenvectors[:3, :3]
    left = (light * squared[:3]) @ light.conj().T
    deficit = np.eye(3) - light @ light.conj().T
    heavy_column = eigenvectors[:3, 3]
    assert np.allclose(deficit, np.outer(heavy_column, heavy_column.conj()))
    quartets = {
        f"{i}{j};{a}{b}": float(
            (light[i, a] * light[j, b] * light[i, b].conj() * light[j, a].conj()).imag
        ) for i, j in combinations(range(3), 2) for a, b in combinations(range(3), 2)
    }
    assert abs(quartets["01;01"]) > 0.002
    coefficients = polynomial(left, deficit)
    assert abs(coefficients[3]) < 1e-14
    roots = sorted(np.roots(coefficients[:3][::-1]).real.tolist())
    assert roots[0] < 0 and abs(roots[1] - squared[3]) < 1e-10
    relative_phases = np.diag(np.exp(1j * np.array([0.2, -0.87, 1.1])))
    coefficients_after = polynomial(
        relative_phases @ left @ relative_phases.conj().T,
        relative_phases @ deficit @ relative_phases.conj().T,
    )
    assert np.allclose(coefficients, coefficients_after)
    readouts = []
    # Rectangular calibration for reconstructed off-diagonal entries and heavy mass squared.
    delta_left, delta_deficit, delta_mass_sq = 1e-5, 1e-6, 1e-4
    for offset in [0.0, 1.0, 4.0]:
        candidate = squared[3] + offset
        total = left + candidate * deficit
        alternate_spectrum = squared.copy()
        alternate_spectrum[3] = candidate
        alternate = (eigenvectors * alternate_spectrum) @ eigenvectors.conj().T
        assert np.allclose(total, alternate[:3, :3])
        edges = cycle_entries(total)
        value = float(np.prod(edges).imag)
        error_edge = (delta_left + (abs(candidate) + delta_mass_sq) * delta_deficit
                      + np.max(np.abs(cycle_entries(deficit))) * delta_mass_sq)
        error_product = 3 * error_edge * (np.max(np.abs(edges)) + error_edge) ** 2
        excludes_zero = abs(value) > error_product
        assert excludes_zero == bool(offset)
        readouts.append(dict(offset=offset, heavy_sq=float(candidate), invariant=value,
                             max_error=float(error_product), excludes_zero=bool(excludes_zero)))
    reconstructed = factor_real_parent(hermitian)
    assert np.allclose(reconstructed @ reconstructed.conj().T, hermitian)
    phi = 0.2
    rotated = parent @ np.diag([1.0, 1.0, 1.0, np.exp(1j * phi)])
    assert np.allclose(rotated @ rotated.conj().T, hermitian)
    strong_shift = float(np.angle(np.linalg.det(rotated) / np.linalg.det(parent)))
    assert np.isclose(strong_shift, phi)
    leakage = []
    for epsilon in [0.001, 0.01]:
        reverse = 1j * epsilon * ordinary[:, 0]
        perturbed = parent.copy()
        perturbed[:3, 3] = reverse
        z = coupling @ np.linalg.solve(ordinary, reverse) / parent[3, 3]
        ratio = np.linalg.det(perturbed) / np.linalg.det(parent)
        assert np.isclose(ratio, 1 - z) and np.isclose(np.angle(ratio), -np.arctan(epsilon / 4))
        leakage.append(dict(epsilon=epsilon, induced_theta=float(np.angle(ratio))))
    return dict(
        input_units="common illustrative mass unit; tree-level finite matrices",
        squared_masses=squared.tolist(), light_weak_quartets=quartets,
        tree_det_down=float(np.linalg.det(parent).real),
        rank_one_deficit_norm=float(np.linalg.norm(deficit, ord=2)),
        polynomial_constant_first=coefficients[:3].tolist(), candidate_roots=roots,
        calibration=dict(delta_left=delta_left, delta_deficit=delta_deficit, delta_mass_sq=delta_mass_sq),
        heavy_mass_tests=readouts, right_phase_same_left_theta_shift=strong_shift,
        reverse_operator=leakage
    )


def same_result(actual, stored):
    if isinstance(actual, dict):
        return actual.keys() == stored.keys() and all(same_result(v, stored[k]) for k, v in actual.items())
    if isinstance(actual, list):
        return len(actual) == len(stored) and all(same_result(a, b) for a, b in zip(actual, stored))
    if isinstance(actual, float):
        return bool(np.isclose(actual, stored, rtol=2e-9, atol=1e-11))
    return actual == stored


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-new", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = compute()
    if args.write_new:
        with RESULT.open("x", encoding="utf-8") as out:
            json.dump(data, out, indent=2, ensure_ascii=False, allow_nan=False)
            out.write("\n")
    else:
        assert same_result(data, json.loads(RESULT.read_text(encoding="utf-8")))
        print("1364 finite weak-moment, heavy-mass and strong-phase tests passed")


if __name__ == "__main__":
    main()

