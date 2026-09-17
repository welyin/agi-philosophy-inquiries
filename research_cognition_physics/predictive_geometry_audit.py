"""Rounds 192-193: probability geometry is not a derived quantum phase space.

The exact finite quotient from round 191 is classical. Square-root coordinates
inherit a Fisher metric but have zero pulled-back projective symplectic form.
"""

import argparse
import json
from pathlib import Path
import unittest

import numpy as np


def probability_tangents(size):
    return np.vstack((np.eye(size - 1), -np.ones(size - 1)))


def ray_data(probabilities, phases=None, include_phase_directions=False):
    probabilities = np.asarray(probabilities, dtype=float)
    if probabilities.ndim != 1 or np.any(probabilities <= 0) or not np.isclose(probabilities.sum(), 1):
        raise ValueError("Use an interior normalized probability vector.")
    size = len(probabilities)
    if phases is None:
        phases = np.zeros(size)
    psi = np.sqrt(probabilities) * np.exp(1j * phases)
    derivatives = (np.exp(1j * phases) / (2 * np.sqrt(probabilities)))[:, None] * probability_tangents(size)
    if include_phase_directions:
        phase_derivatives = 1j * np.diag(psi)[:, :size - 1]
        derivatives = np.column_stack((derivatives, phase_derivatives))
    horizontal = derivatives - psi[:, None] * (psi.conj() @ derivatives)[None, :]
    gram = horizontal.conj().T @ horizontal
    return psi, horizontal, gram.real, gram.imag


def root_density(probabilities):
    vector = np.sqrt(np.asarray(probabilities, dtype=float))
    return np.outer(vector, vector)


class PredictiveGeometryAuditTests(unittest.TestCase):
    def test_square_root_probability_metric_is_one_quarter_fisher(self):
        p = np.arange(1., 9.) / 36
        _, _, metric, _ = ray_data(p)
        tangent = probability_tangents(8)
        fisher = tangent.T @ np.diag(1 / p) @ tangent
        np.testing.assert_allclose(metric, fisher / 4, atol=2e-15)
        self.assertGreater(np.linalg.eigvalsh(metric).min(), 0)

    def test_real_probability_slice_has_zero_pulled_back_symplectic_form(self):
        for size in (3, 4, 8):
            p = np.arange(1., size + 1)
            p /= p.sum()
            psi, derivatives, metric, symplectic = ray_data(p)
            np.testing.assert_array_equal(symplectic, np.zeros_like(metric))
            np.testing.assert_allclose(psi.conj() @ derivatives, 0, atol=2e-16)
            self.assertGreater(np.linalg.norm((1j * derivatives).imag), 0)
            np.testing.assert_array_equal(derivatives.imag, np.zeros_like(derivatives.imag))

    def test_adding_independent_phases_completes_the_projective_kahler_tangent(self):
        for size in (2, 3, 8):
            p = np.arange(1., size + 1)
            p /= p.sum()
            _, derivatives, metric, symplectic = ray_data(
                p, np.linspace(0, .7, size), include_phase_directions=True)
            dimension = 2 * (size - 1)
            self.assertEqual(np.linalg.matrix_rank(symplectic), dimension)
            complex_structure = -np.linalg.solve(metric, symplectic)
            np.testing.assert_allclose(complex_structure @ complex_structure,
                                       -np.eye(dimension), atol=2e-14)
            np.testing.assert_allclose(derivatives @ complex_structure, 1j * derivatives, atol=3e-15)
            np.testing.assert_allclose(complex_structure.T @ metric @ complex_structure, metric, atol=1e-14)

    def test_hidden_phase_changes_new_interference_readout_but_not_original_probabilities(self):
        plus = np.array([1., 1.]) / np.sqrt(2)
        minus = np.array([1., -1.]) / np.sqrt(2)
        np.testing.assert_array_equal(abs(plus)**2, abs(minus)**2)
        effect = np.outer(plus, plus)
        self.assertAlmostEqual(plus @ effect @ plus, 1)
        self.assertAlmostEqual(minus @ effect @ minus, 0)

    def test_interpreting_square_roots_as_quantum_pure_states_breaks_preparation_mixing(self):
        for size in (2, 8):
            p = np.ones(size) / size
            actual_randomized_preparation = np.eye(size) / size
            declared_pure_root = root_density(p)
            self.assertAlmostEqual(np.linalg.matrix_rank(declared_pure_root), 1)
            distance = np.abs(np.linalg.eigvalsh(declared_pure_root - actual_randomized_preparation)).sum() / 2
            self.assertAlmostEqual(distance, 1 - 1 / size)
            np.testing.assert_allclose(np.diag(declared_pure_root), np.diag(actual_randomized_preparation), atol=2e-16)
            self.assertAlmostEqual(np.trace(declared_pure_root @ declared_pure_root), 1)
            self.assertAlmostEqual(np.trace(declared_pure_root @ actual_randomized_preparation), 1 / size)

    def test_affine_diagonal_embedding_preserves_mixing_and_independent_composition(self):
        p = np.array([.25, .75])
        q = np.array([.75, .25])
        weight = .25
        np.testing.assert_array_equal(np.diag(weight * p + (1-weight) * q),
                                      weight * np.diag(p) + (1-weight) * np.diag(q))
        np.testing.assert_array_equal(np.diag(np.kron(p, q)), np.kron(np.diag(p), np.diag(q)))
        np.testing.assert_array_equal(np.diag(p) @ np.diag(q), np.diag(q) @ np.diag(p))

    def test_classical_feedback_has_a_quantum_embedding_without_creating_coherence(self):
        transition = np.array([[.75, .25], [.25, .75]])
        kraus = []
        for i, j in np.ndindex(2, 2):
            operator = np.zeros((2, 2))
            operator[i, j] = np.sqrt(transition[i, j])
            kraus.append(operator)
        np.testing.assert_allclose(sum(k.T @ k for k in kraus), np.eye(2), atol=2e-16)
        for p in (np.array([.25, .75]), np.array([1., 0.])):
            output = sum(k @ np.diag(p) @ k.T for k in kraus)
            np.testing.assert_allclose(output, np.diag(transition @ p), atol=2e-16)
        root = root_density([.5, .5])
        output = sum(k @ root @ k.T for k in kraus)
        np.testing.assert_allclose(output, np.eye(2) / 2, atol=2e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(PredictiveGeometryAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 192,
        "also_supports_round": 193,
        "round_191_uncalibrated_simplex_vertices": 8,
        "normalized_belief_dimension": 7,
        "square_root_metric": "one quarter of the Fisher metric in the chosen Fubini-Study normalization",
        "pulled_back_projective_symplectic_form": "zero",
        "independent_relative_phases_needed_for_full_pure_ray_patch": 7,
        "resulting_projective_patch_real_dimension_if_phases_are_added": 14,
        "phase_controls_derived_from_cognitive_tasks": False,
        "uniform_eight_label_root_vs_mixture_trace_distance": "7/8",
        "finite_classical_affine_projection_can_cover_all_nontrivial_quantum_states": False,
        "theorem_scope": "Finite classical latent preparations, all randomizations allowed, affine map into standard quantum density matrices; not arbitrary continuous SoCA or approximate models",
        "odd_belief_dimension_alone_excludes_quantum_theory": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("predictive_geometry_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
