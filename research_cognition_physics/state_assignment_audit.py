"""Round 177: a universal positive consistent affine assignment is a product.

The all-state claim for ordinary real and complex finite quantum states is proved
in the note. These examples audit separate hypotheses, not a general optimizer.
"""

import argparse
import json
from pathlib import Path
import unittest

import numpy as np
from cut_consistency_audit import trace_keep


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])
P0, P1 = np.diag([1., 0.]), np.diag([0., 1.])
PLUS, MINUS = (I + X) / 2, (I - X) / 2


def canonical_pure_assignment(rho):
    values, vectors = np.linalg.eigh(rho)
    root = (vectors * np.sqrt(np.maximum(values, 0))) @ vectors.conj().T
    psi = root.ravel()
    return np.outer(psi, psi.conj())


def correlated_linear_assignment(rho):
    """Consistent extension of classical copying; not positive on all states."""
    off_diagonal = rho - np.diag(np.diag(rho))
    return rho[0, 0] * np.kron(P0, P0) + rho[1, 1] * np.kron(P1, P1) + np.kron(off_diagonal, I / 2)


class StateAssignmentAuditTests(unittest.TestCase):
    def test_fixed_product_assignment_satisfies_all_three_requirements(self):
        rng = np.random.default_rng(177)
        for complex_model in (False, True):
            sigma = np.array([[.4, .1], [.1, .6]])
            raw = rng.normal(size=(3, 3))
            if complex_model:
                raw = raw + 1j * rng.normal(size=raw.shape)
            rho = raw @ raw.conj().T
            rho /= np.trace(rho)
            assigned = np.kron(rho, sigma)
            self.assertGreaterEqual(np.linalg.eigvalsh(assigned).min(), -1e-15)
            np.testing.assert_allclose(trace_keep(assigned, (3, 2), (0,)), rho, atol=2e-16)
            np.testing.assert_allclose(np.kron(.3 * rho + .7 * np.eye(3) / 3, sigma),
                                       .3 * assigned + .7 * np.kron(np.eye(3) / 3, sigma), atol=2e-16)

    def test_two_decompositions_force_same_environment_state(self):
        # For each real symmetric environment basis component, four coefficients
        # s0,s1,s+,s- must satisfy the same block-matrix equality.
        constraint = np.column_stack([P0.ravel(), P1.ravel(), -PLUS.ravel(), -MINUS.ravel()])
        self.assertEqual(np.linalg.matrix_rank(constraint), 3)
        np.testing.assert_array_equal(constraint @ np.ones(4), np.zeros(4))
        y = np.array([[0., -1j], [1j, 0.]])
        complex_constraint = np.column_stack([P0.ravel(), P1.ravel(),
                                             -((I + y) / 2).ravel(), -((I - y) / 2).ravel()])
        stacked = np.vstack([complex_constraint.real, complex_constraint.imag])
        self.assertEqual(np.linalg.matrix_rank(stacked), 3)
        np.testing.assert_array_equal(stacked @ np.ones(4), np.zeros(8))

    def test_classical_copy_is_valid_only_on_its_commuting_domain(self):
        for p in (0., .2, .5, 1.):
            rho = np.diag([p, 1 - p])
            assigned = correlated_linear_assignment(rho)
            self.assertGreaterEqual(np.linalg.eigvalsh(assigned).min(), 0)
            np.testing.assert_array_equal(trace_keep(assigned, (2, 2), (0,)), rho)
        assigned = correlated_linear_assignment(I / 2)
        self.assertGreater(np.linalg.norm(assigned - np.kron(I / 2, I / 2)), 0)

    def test_consistent_linear_correlated_extension_has_exact_negative_witness(self):
        assigned = correlated_linear_assignment(PLUS)
        witness = np.outer([1., 0., -2., 0.], [1., 0., -2., 0.]) / 5
        np.testing.assert_array_equal(trace_keep(assigned, (2, 2), (0,)), PLUS)
        self.assertAlmostEqual(np.trace(witness @ assigned), -.1)
        self.assertAlmostEqual(np.linalg.eigvalsh(assigned).min(), (1 - np.sqrt(2)) / 4)

    def test_positive_coherent_copy_drops_marginal_consistency(self):
        embedding = np.array([[1., 0.], [0., 0.], [0., 0.], [0., 1.]])
        assigned = embedding @ PLUS @ embedding.T
        self.assertEqual(np.trace(assigned @ assigned), 1)
        np.testing.assert_array_equal(trace_keep(assigned, (2, 2), (0,)), I / 2)
        self.assertEqual(np.abs(np.linalg.eigvalsh(PLUS - I / 2)).sum() / 2, .5)

    def test_pure_assignment_exists_but_is_not_affine(self):
        assigned = canonical_pure_assignment(I / 2)
        averaged = (canonical_pure_assignment(P0) + canonical_pure_assignment(P1)) / 2
        for joint in (assigned, averaged):
            np.testing.assert_allclose(trace_keep(joint, (2, 2), (0,)), I / 2, atol=2e-16)
        self.assertAlmostEqual(np.trace(assigned @ assigned), 1)
        self.assertEqual(np.trace(averaged @ averaged), .5)
        self.assertAlmostEqual(np.abs(np.linalg.eigvalsh(assigned - averaged)).sum() / 2, .5)

    def test_purification_consistency_for_complex_and_real_inputs(self):
        for rho in (np.array([[.6, .2], [.2, .4]]), np.array([[.6, .2j], [-.2j, .4]])):
            assigned = canonical_pure_assignment(rho)
            np.testing.assert_allclose(trace_keep(assigned, (2, 2), (0,)), rho, atol=5e-16)
            self.assertAlmostEqual(np.trace(assigned @ assigned).real, 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(StateAssignmentAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 177,
        "assignment_assumptions": ["all local quantum states", "affine", "positive", "local marginal consistent"],
        "real_and_complex_assignment_conclusion": "A(rho)=rho tensor fixed_sigma",
        "classical_conditional_assignments_can_correlate": True,
        "invalid_linear_extension_negative_probability_exact": "-1/10",
        "canonical_purification_nonaffinity_trace_distance": "1/2",
        "state_purification_exists_implies_unknown_input_purification_channel": False,
        "assignment_theorem_is_claimed_original": False,
        "cognitive_principles_force_quantum_structure": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("state_assignment_audit_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
