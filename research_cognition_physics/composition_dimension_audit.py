"""Round 180: conditional dimension selection in a declared matrix family.

Dimension matching is necessary, not a reconstruction of state cones or dynamics.
The quaternionic count is an obstruction to the naive proposed composite, not a
claim that an ordinary quaternionic tensor product defines a physical theory.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np


def state_dimension(n, q):
    return n + q * n * (n - 1) // 2


def matrix_basis(n, complex_model=False):
    basis = []
    for i in range(n):
        entry = np.zeros((n, n), dtype=complex if complex_model else float)
        entry[i, i] = 1
        basis.append(entry)
    for i in range(n):
        for j in range(i + 1, n):
            symmetric = np.zeros((n, n), dtype=complex if complex_model else float)
            symmetric[i, j] = symmetric[j, i] = 1
            basis.append(symmetric)
            if complex_model:
                imaginary = np.zeros((n, n), dtype=complex)
                imaginary[i, j], imaginary[j, i] = 1j, -1j
                basis.append(imaginary)
    return basis


def realification(matrix):
    return np.block([[matrix.real, -matrix.imag], [matrix.imag, matrix.real]])


def product_pairing(n, m, complex_model=False):
    products = [np.kron(a, b) for a in matrix_basis(n, complex_model)
                for b in matrix_basis(m, complex_model)]
    joint = matrix_basis(n * m, complex_model)
    return np.array([[np.trace(a @ b).real for b in joint] for a in products])


class CompositionDimensionAuditTests(unittest.TestCase):
    def test_dimension_defect_polynomial_coefficients_exactly_factor(self):
        # Compare polynomial coefficients in q, not a grid in q.
        for n in range(2, 9):
            for m in range(2, 9):
                a, b = n * (n - 1) // 2, m * (m - 1) // 2
                coefficients = [0, n * m * (n * m - 1) // 2 - n * b - m * a, -a * b]
                scale = Fraction(n * m * (n - 1) * (m - 1), 4)
                self.assertEqual(coefficients, [0, 2 * scale, -scale])

    def test_real_product_measurements_have_the_predicted_missing_directions(self):
        for n, m in ((2, 2), (2, 3), (3, 3)):
            pairing = product_pairing(n, m)
            rank = np.linalg.matrix_rank(pairing)
            self.assertEqual(rank, state_dimension(n, 1) * state_dimension(m, 1))
            self.assertEqual(pairing.shape[1] - rank, n * m * (n - 1) * (m - 1) // 4)

    def test_complex_product_measurements_span_all_joint_directions(self):
        for n, m in ((2, 2), (2, 3), (3, 3)):
            pairing = product_pairing(n, m, True)
            self.assertEqual(pairing.shape, ((n * m) ** 2, (n * m) ** 2))
            self.assertEqual(np.linalg.matrix_rank(pairing), (n * m) ** 2)

    def test_naive_quaternionic_composite_lacks_space_for_independent_products(self):
        for n, m in ((2, 2), (2, 3), (3, 3)):
            required = state_dimension(n, 4) * state_dimension(m, 4)
            proposed = state_dimension(n * m, 4)
            self.assertLess(proposed, required)
        self.assertEqual((state_dimension(4, 4), state_dimension(2, 4) ** 2), (28, 36))

    def test_equivalent_real_coordinates_keep_complex_operational_dimension(self):
        for n in range(1, 5):
            j = np.block([[np.zeros((n, n)), -np.eye(n)], [np.eye(n), np.zeros((n, n))]])
            ambient = matrix_basis(2 * n)
            constraints = np.column_stack([(a @ j - j @ a).ravel() for a in ambient])
            kernel_dimension = len(ambient) - np.linalg.matrix_rank(constraints)
            self.assertEqual(kernel_dimension, n * n)
            embedded = [realification(a) for a in matrix_basis(n, True)]
            for a in embedded:
                np.testing.assert_array_equal(a.T, a)
                np.testing.assert_array_equal(a @ j, j @ a)
            self.assertEqual(np.linalg.matrix_rank(np.array([a.ravel() for a in embedded])), n * n)

    def test_classical_survivor_and_trivial_system_do_not_select_complex(self):
        for n in range(1, 9):
            self.assertEqual(state_dimension(n * 3, 0), state_dimension(n, 0) * 3)
            for q in (0, 1, 2, 4):
                self.assertEqual(state_dimension(n, q), state_dimension(1, q) * state_dimension(n, q))
        # Pure states of the classical 2-by-3 simplex have point marginals.
        for vertex in np.eye(6):
            marginal = vertex.reshape(2, 3).sum(axis=1)
            self.assertEqual(np.count_nonzero(marginal), 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CompositionDimensionAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 180,
        "declared_dimension_family": "K_q(n)=n+q*n*(n-1)/2",
        "assumed_capacity_composition": "n_AB=n_A*n_B",
        "dimension_defect": "q*(2-q)*n*m*(n-1)*(m-1)/4",
        "roots_for_nontrivial_factors": [0, 2],
        "two_capacity_systems": [{"q": q, "local_K": state_dimension(2, q),
                                  "joint_K": state_dimension(4, q),
                                  "product_K": state_dimension(2, q) ** 2} for q in (0, 1, 2, 4)],
        "independent_product_assumptions_are_required": True,
        "dimension_match_alone_proves_full_quantum_theory": False,
        "matrix_family_is_derived_from_cognition": False,
        "equivalent_real_J_representation_changes_operational_K": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("composition_dimension_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
