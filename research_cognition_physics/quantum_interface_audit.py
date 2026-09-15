"""Exact quantum representations and missing structure of the existing disk model.

These are representations and counterexamples, not a derivation of quantum
theory from cognitive principles. No existing physical operations are changed.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bounded_responses import translation_matrix
from outcome_updates import WINDOW_ATTENUATION as ALPHA, validate_unit_interval
from retention_tradeoff import fresh_branch_matrix


IDENTITY = np.eye(2, dtype=complex)
PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
PAULI_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
PAULI_Z = np.diag([1, -1]).astype(complex)


def density_from_summary(summary):
    w, u, v = np.asarray(summary, dtype=float)
    return (w * IDENTITY + (u * PAULI_Z + v * PAULI_X) / ALPHA) / 2


def summary_from_density(density):
    return np.array([np.trace(density).real, ALPHA * np.trace(PAULI_Z @ density).real,
                     ALPHA * np.trace(PAULI_X @ density).real])


def rotation_unitary(setting):
    return math.cos(setting / 2) * IDENTITY - 1j * math.sin(setting / 2) * PAULI_Y


def noise_kraus(completion="real_rotation"):
    if completion == "real_rotation":
        # -i sigma_y is real. This channel fixes the unobserved y coordinate.
        return (math.sqrt((1 + ALPHA) / 2) * IDENTITY,
                math.sqrt((1 - ALPHA) / 2) * (-1j * PAULI_Y))
    if completion == "isotropic":
        return tuple(math.sqrt(weight) * operator for weight, operator in (
            ((1 + 3 * ALPHA) / 4, IDENTITY), ((1 - ALPHA) / 4, PAULI_X),
            ((1 - ALPHA) / 4, PAULI_Y), ((1 - ALPHA) / 4, PAULI_Z)))
    raise ValueError("Unknown completion.")


def branch_kraus(setting, outcome, contrast, completion="real_rotation"):
    validate_unit_interval(contrast)
    if outcome not in (0, 1):
        raise ValueError("Use outcome 0 or 1.")
    sign = 2 * outcome - 1
    measurement = np.diag(np.sqrt([(1 + sign * contrast) / 2, (1 - sign * contrast) / 2]))
    rotation = rotation_unitary(setting)
    return tuple(rotation @ measurement @ operator @ rotation.conj().T for operator in noise_kraus(completion))


def apply_branch(density, setting, outcome, contrast, completion="real_rotation"):
    return sum(k @ density @ k.conj().T for k in branch_kraus(setting, outcome, contrast, completion))


def effective_readout(setting, outcome, contrast):
    sign = 2 * outcome - 1
    direction = math.cos(setting) * PAULI_Z + math.sin(setting) * PAULI_X
    return (IDENTITY + sign * ALPHA * contrast * direction) / 2


def choi_matrix(kraus):
    vectors = [operator.reshape(4, order="F") for operator in kraus]
    return sum(np.outer(vector, vector.conj()) for vector in vectors)


def two_rebit_states():
    hidden_joint_coordinate = np.kron(PAULI_Y, PAULI_Y)
    identity = np.eye(4)
    return (identity + hidden_joint_coordinate) / 4, (identity - hidden_joint_coordinate) / 4


class QuantumInterfaceAuditTests(unittest.TestCase):
    def test_disk_cone_is_exactly_the_real_two_by_two_positive_cone(self):
        for w, fraction, angle in product((0., .3, 1., 2.), (0., .2, 1.), np.linspace(-2, 5, 9)):
            z = np.array([w, ALPHA * w * fraction * math.cos(angle), ALPHA * w * fraction * math.sin(angle)])
            density = density_from_summary(z)
            np.testing.assert_allclose(density.imag, 0, atol=0)
            np.testing.assert_allclose(summary_from_density(density), z, atol=5e-16)
            np.testing.assert_allclose(np.linalg.eigvalsh(density), [w * (1 - fraction) / 2, w * (1 + fraction) / 2], atol=5e-16)

    def test_both_kraus_completions_match_the_original_entire_branch_map(self):
        # Equality on the linear basis proves equality on every disk preparation.
        for setting, outcome, eta, completion in product((0., .29, -1.7), (0, 1), (0., .144739232095145, .5, 1.), ("real_rotation", "isotropic")):
            for z in np.eye(3):
                actual = apply_branch(density_from_summary(z), setting, outcome, eta, completion)
                expected = density_from_summary(fresh_branch_matrix(setting, outcome, eta) @ z)
                np.testing.assert_allclose(actual, expected, atol=7e-16, rtol=0)

    def test_complete_positivity_trace_preservation_and_actual_readout_effects(self):
        for setting, eta, completion in product((0., .31), (0., .15, .7, 1.), ("real_rotation", "isotropic")):
            total = np.zeros((2, 2), dtype=complex)
            for outcome in (0, 1):
                operators = branch_kraus(setting, outcome, eta, completion)
                self.assertGreaterEqual(np.linalg.eigvalsh(choi_matrix(operators)).min(), -4e-16)
                effect = sum(k.conj().T @ k for k in operators)
                np.testing.assert_allclose(effect, effective_readout(setting, outcome, eta), atol=5e-16)
                total += effect
                if completion == "real_rotation":
                    self.assertLess(max(np.max(np.abs(k.imag)) for k in operators), 1e-16)
            np.testing.assert_allclose(total, IDENTITY, atol=5e-16)

    def test_rotations_and_adaptive_record_probabilities_match_without_extra_states(self):
        z0 = np.array([1., .37, -.28])
        for angle in (-.8, .4, 2.):
            u = rotation_unitary(angle)
            np.testing.assert_allclose(u @ density_from_summary(z0) @ u.conj().T,
                                       density_from_summary(translation_matrix(angle, (1.,)) @ z0), atol=4e-16)
        for record in product((0, 1), repeat=5):
            z = z0.copy()
            density = density_from_summary(z0)
            for step, outcome in enumerate(record):
                setting = .19 * step + (.37 if step and record[step - 1] else -.21)
                eta = (.12, .5, .9, .25, .7)[step]
                z = fresh_branch_matrix(setting, outcome, eta) @ z
                density = apply_branch(density, setting, outcome, eta)
            np.testing.assert_allclose(density, density_from_summary(z), atol=3e-16)
            self.assertAlmostEqual(np.trace(density).real, z[0], places=14)

    def test_two_completions_disagree_on_a_third_coordinate_absent_from_the_project(self):
        outside_disk = (IDENTITY + PAULI_Y) / 2
        outputs = [sum(apply_branch(outside_disk, 0, outcome, 0, completion) for outcome in (0, 1))
                   for completion in ("real_rotation", "isotropic")]
        self.assertAlmostEqual(np.trace(outputs[0] @ PAULI_Y).real, 1, places=14)
        self.assertAlmostEqual(np.trace(outputs[1] @ PAULI_Y).real, ALPHA, places=14)
        self.assertGreater(np.linalg.norm(outputs[0] - outputs[1]), .007)

    def test_zero_strength_fresh_update_still_contracts_the_existing_disk(self):
        density = density_from_summary([1, ALPHA, 0])
        output = sum(apply_branch(density, 0, outcome, 0) for outcome in (0, 1))
        self.assertAlmostEqual(np.trace(output @ output).real, (1 + ALPHA**2) / 2, places=14)
        self.assertLess(np.trace(output @ output).real, 1)
        self.assertAlmostEqual(np.trace(output @ PAULI_Z).real, ALPHA, places=14)

    def test_real_composite_has_a_coordinate_invisible_to_all_local_real_products(self):
        plus, minus = two_rebit_states()
        for density in (plus, minus):
            np.testing.assert_allclose(density.imag, 0, atol=0)
            self.assertGreaterEqual(np.linalg.eigvalsh(density).min(), 0)
            self.assertAlmostEqual(np.trace(density).real, 1)
        for a, b in product((IDENTITY, PAULI_X, PAULI_Z), repeat=2):
            observable = np.kron(a, b)
            self.assertAlmostEqual(np.trace((plus - minus) @ observable).real, 0, places=14)
        effect = (np.eye(4) + np.kron(PAULI_Y, PAULI_Y)) / 2
        self.assertAlmostEqual(np.trace(plus @ effect).real, 1)
        self.assertAlmostEqual(np.trace(minus @ effect).real, 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(QuantumInterfaceAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "single_system_real_quantum_representation": True,
        "scope": "Exact representation of the disk, rotations and the longitudinal=1 fresh family; not all real-quantum effects and transformations are available.",
        "quantum_theory_derived_from_cognition": False,
        "original_readout_visibility": "alpha*eta; its maximum alpha is strictly below one",
        "two_indistinguishable_completions_on_existing_protocols": ["random real rotation noise", "isotropic depolarizing noise"],
        "third_coordinate_expectation_after_zero_strength_readout": [1., ALPHA],
        "zero_strength_output_purity_from_boundary_state": (1 + ALPHA**2) / 2,
        "two_rebit_global_states": "rho_+/- = (I_4 +/- sigma_y tensor sigma_y)/4",
        "local_product_statistics_identical": True,
        "global_effect_probabilities": [1., 0.],
        "next_question": "Which cognitive composition rules select joint states and operations, and do they determine the missing complex-quantum structure?",
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("quantum_interface_audit_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
