"""Local operation equivalence need not survive a real joint-system context.

The real Kraus framework is a declared candidate, not a cognitive deduction.
An extra scalar completes the action of a real 2x2 operation on matrix space.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bipartite_composition import FULL_BASIS, HIDDEN, interaction, visible_summary
from certified_intervals import Interval as I, SCALE, sin_interval
from quantum_interface_audit import (ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                     branch_kraus, effective_readout, noise_kraus,
                                     two_rebit_states)


REAL_BASIS = (IDENTITY, PAULI_X, PAULI_Z)
J = (-1j * PAULI_Y).real


def real_noise_kraus(completion):
    if completion == "real_rotation":
        return tuple(k.real for k in noise_kraus(completion))
    if completion == "isotropic":
        weights = ((1 + 3*ALPHA)/4, (1 - ALPHA)/4, (1 - ALPHA)/4, (1 - ALPHA)/4)
        return tuple(math.sqrt(p) * k for p, k in zip(weights, (IDENTITY.real, PAULI_X.real, J, PAULI_Z.real)))
    raise ValueError("Unknown completion.")


def apply_operation(kraus, matrix):
    return sum(k @ matrix @ k.conj().T for k in kraus)


def local_operation_summary(kraus):
    return np.array([[np.trace(a @ apply_operation(kraus, b)).real / 2 for b in REAL_BASIS]
                     for a in REAL_BASIS])


def determinant_parameter(kraus):
    if any(np.max(np.abs(np.asarray(k).imag)) > 1e-14 for k in kraus):
        raise ValueError("This descriptor requires a real Kraus representation.")
    return float(sum(np.linalg.det(k.real) for k in kraus))


def apply_on_first(kraus, density):
    dimension = density.shape[0] // 2
    return apply_operation(tuple(np.kron(k, np.eye(dimension)) for k in kraus), density)


def completed_action(summary, kappa, matrix):
    coordinates = np.array([np.trace(b @ matrix) / 2 for b in REAL_BASIS])
    output = sum(value * b for value, b in zip(summary @ coordinates, REAL_BASIS))
    return output + kappa * np.trace(PAULI_Y @ matrix) * PAULI_Y / 2


def probe_kappa(kraus):
    density = apply_on_first(kraus, two_rebit_states()[0])
    gate = interaction(math.pi / 2)
    effect = np.kron(IDENTITY, effective_readout(0, 1, 1.))
    rotated = gate @ density @ gate.conj().T
    # General selective maps can have total probability different from one.
    mass = np.trace(density).real
    signed_mean = np.trace((2*effect - np.eye(4)) @ rotated).real
    return {"total_mass": float(mass), "positive_probability": float(np.trace(effect @ rotated).real),
            "signed_record_mean": float(signed_mean), "inferred_kappa": float(signed_mean / ALPHA)}


def strict_gap_certificate():
    alpha = 4 * sin_interval(I.rational(1, 4))
    gap = alpha * (1 - alpha) / 2
    return {"TV_formula": "alpha*(1-alpha)/2", "TV_interval": gap.floats(),
            "strict_lower_exact": str(Fraction(gap.lo, SCALE)),
            "strictly_above_one_two_hundredth": gap.lo * 200 > SCALE}


class AncillaryOperationEquivalenceTests(unittest.TestCase):
    def test_both_completions_have_real_kraus_representations(self):
        for completion in ("real_rotation", "isotropic"):
            kraus = real_noise_kraus(completion)
            np.testing.assert_allclose(sum(k.T @ k for k in kraus), np.eye(2), atol=3e-16)
            for matrix in (*REAL_BASIS, PAULI_Y):
                np.testing.assert_allclose(apply_operation(kraus, matrix), apply_operation(noise_kraus(completion), matrix), atol=3e-16)

    def test_the_entire_local_linear_map_is_identical_but_kappa_differs(self):
        keep, isotropic = [real_noise_kraus(c) for c in ("real_rotation", "isotropic")]
        np.testing.assert_allclose(local_operation_summary(keep), np.diag([1, ALPHA, ALPHA]), atol=3e-16)
        np.testing.assert_allclose(local_operation_summary(keep), local_operation_summary(isotropic), atol=3e-16)
        self.assertAlmostEqual(determinant_parameter(keep), 1.)
        self.assertAlmostEqual(determinant_parameter(isotropic), ALPHA)

    def test_local_equivalence_includes_all_selected_adaptive_records(self):
        initial = (IDENTITY + .3*PAULI_X - .4*PAULI_Z) / 2
        for record in product((0, 1), repeat=5):
            states = [initial.copy(), initial.copy()]
            for step, outcome in enumerate(record):
                setting = .27*step + .19*sum(record[:step])
                eta = (.1, .4, .8, .6, .3)[step]
                for index, completion in enumerate(("real_rotation", "isotropic")):
                    states[index] = apply_operation(branch_kraus(setting, outcome, eta, completion), states[index])
            np.testing.assert_allclose(states[0], states[1], atol=2e-16)

    def test_completions_are_distinguished_on_an_existing_hidden_joint_preparation(self):
        outputs = [apply_on_first(real_noise_kraus(c), two_rebit_states()[0])
                   for c in ("real_rotation", "isotropic")]
        for value, density in zip((1., ALPHA), outputs):
            np.testing.assert_allclose(density, (np.eye(4) + value*HIDDEN)/4, atol=2e-16)
        np.testing.assert_allclose(visible_summary(outputs[0]), visible_summary(outputs[1]), atol=2e-16)
        probes = [probe_kappa(real_noise_kraus(c)) for c in ("real_rotation", "isotropic")]
        gap = probes[0]["positive_probability"] - probes[1]["positive_probability"]
        self.assertAlmostEqual(gap, ALPHA*(1-ALPHA)/2, places=14)
        self.assertTrue(strict_gap_certificate()["strictly_above_one_two_hundredth"])

    def test_reflection_channels_have_continuous_assisted_real_implementations(self):
        rng = np.random.default_rng(52)
        aux_matrix = rng.normal(size=(2, 2))
        aux = aux_matrix @ aux_matrix.T
        aux /= np.trace(aux)
        for reflection in (PAULI_X.real, PAULI_Z.real):
            gate = np.kron(reflection, J)
            np.testing.assert_array_equal(gate.T @ gate, np.eye(4))
            self.assertAlmostEqual(np.linalg.det(gate), 1.)
            matrix = rng.normal(size=(2, 2))
            density = matrix @ matrix.T
            density /= np.trace(density)
            output = gate @ np.kron(density, aux) @ gate.T
            reduced = np.einsum("abcb->ac", output.reshape(2, 2, 2, 2))
            np.testing.assert_allclose(reduced, reflection @ density @ reflection.T, atol=2e-16)

    def test_single_missing_scalar_follows_from_exact_integer_determinant_identity(self):
        for values in product((-2, 0, 3), repeat=4):
            matrix = np.array(values, dtype=int).reshape(2, 2)
            determinant = int(matrix[0, 0]*matrix[1, 1] - matrix[0, 1]*matrix[1, 0])
            np.testing.assert_array_equal(matrix @ J @ matrix.T, determinant * J)

    def test_augmented_summary_predicts_all_matrix_inputs_of_random_real_operations(self):
        rng = np.random.default_rng(152)
        for _ in range(10):
            kraus = tuple(rng.normal(size=(2, 2)) / 4 for _ in range(3))
            summary, kappa = local_operation_summary(kraus), determinant_parameter(kraus)
            for matrix in (*REAL_BASIS, PAULI_Y, *np.eye(4).reshape(4, 2, 2)):
                np.testing.assert_allclose(completed_action(summary, kappa, matrix), apply_operation(kraus, matrix), atol=4e-16)
            self.assertAlmostEqual(probe_kappa(kraus)["inferred_kappa"], kappa, places=14)

    def test_equal_augmented_summaries_survive_kraus_changes_and_arbitrary_real_ancillas(self):
        rng = np.random.default_rng(252)
        kraus = np.array([rng.normal(size=(2, 2)) / 5 for _ in range(3)])
        orthogonal, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        equivalent = np.einsum("ij,jab->iab", orthogonal, kraus)
        np.testing.assert_allclose(local_operation_summary(kraus), local_operation_summary(equivalent), atol=2e-16)
        self.assertAlmostEqual(determinant_parameter(kraus), determinant_parameter(equivalent), places=14)
        for dimension in (2, 3, 4):
            matrix = rng.normal(size=(2*dimension, 2*dimension))
            density = matrix @ matrix.T
            density /= np.trace(density)
            np.testing.assert_allclose(apply_on_first(kraus, density), apply_on_first(equivalent, density), atol=2e-16)

    def test_augmented_descriptor_composes_and_mixes(self):
        a, b = [real_noise_kraus(c) for c in ("real_rotation", "isotropic")]
        composed = tuple(kb @ ka for ka in a for kb in b)
        np.testing.assert_allclose(local_operation_summary(composed), local_operation_summary(b) @ local_operation_summary(a), atol=3e-16)
        self.assertAlmostEqual(determinant_parameter(composed), determinant_parameter(b)*determinant_parameter(a), places=14)
        mixture = tuple(math.sqrt(.3)*k for k in a) + tuple(math.sqrt(.7)*k for k in b)
        self.assertAlmostEqual(determinant_parameter(mixture), .3 + .7*ALPHA, places=14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AncillaryOperationEquivalenceTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"scope": "Real 2x2 Kraus operations with usual real tensor composition and declared joint probes",
              "local_operation_equivalence_implies_ancillary_equivalence": False,
              "local_matrix_basis_order": ["I", "X", "Z"],
              "additional_operation_parameter": "kappa = sum det(K_j), equivalently Phi(Y)=kappa*Y",
              "complete_real_operation_descriptor": "(3x3 real symmetric-sector map, kappa)",
              "kappa_is_a_new_local_state_coordinate": False,
              "probe_reports": {c: probe_kappa(real_noise_kraus(c)) for c in ("real_rotation", "isotropic")},
              "readable_distinction_certificate": strict_gap_certificate(),
              "safe_replacement_criterion_within_this_framework": "Equal symmetric-sector map AND equal kappa",
              "all_real_ancillas_covered_analytically": True,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("ancillary_operation_equivalence_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
