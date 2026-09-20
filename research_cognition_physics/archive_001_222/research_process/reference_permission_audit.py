"""Round 187: the exact reference sector needed for assisted real tomography.

The proof is in research_note_187.md. Matrix checks cover basis decomposition,
allowed effects, reconstruction, and the failure of real separable resources.
"""

import argparse
import json
from pathlib import Path
import unittest

import numpy as np


J = np.array([[0., -1.], [1., 0.]])
H = np.kron(J, J)


def symmetric_basis(dimension):
    basis = []
    for i in range(dimension):
        diagonal = np.zeros((dimension, dimension))
        diagonal[i, i] = 1
        basis.append(diagonal)
        for j in range(i + 1, dimension):
            off_diagonal = np.zeros((dimension, dimension))
            off_diagonal[i, j] = off_diagonal[j, i] = 1
            basis.append(off_diagonal)
    return basis


def skew_basis(dimension):
    basis = []
    for i in range(dimension):
        for j in range(i + 1, dimension):
            matrix = np.zeros((dimension, dimension))
            matrix[i, j], matrix[j, i] = -1, 1
            basis.append(matrix)
    return basis


def reorder(matrix, dimensions, order):
    count = len(dimensions)
    axes = tuple(order) + tuple(count + index for index in order)
    return matrix.reshape(tuple(dimensions) * 2).transpose(axes).reshape(matrix.shape)


def partial_transpose(matrix, dimensions, sites):
    count = len(dimensions)
    axes = list(range(2 * count))
    for site in sites:
        axes[site], axes[site + count] = axes[site + count], axes[site]
    return matrix.reshape(tuple(dimensions) * 2).transpose(axes).reshape(matrix.shape)


def relation_state(value):
    if not -1 <= value <= 1:
        raise ValueError("The relation parameter must lie in [-1,1].")
    return (np.eye(4) + value * H) / 4


def local_product_in_target_order(left, right, target_dims=(2, 2), ref_dims=(2, 2)):
    a, b = target_dims
    r, s = ref_dims
    return reorder(np.kron(left, right), (a, r, b, s), (0, 2, 1, 3))


def complete_product_moments(joint):
    return np.array([np.trace(joint @ local_product_in_target_order(a, b))
                     for a in symmetric_basis(4) for b in symmetric_basis(4)])


def assisted_reconstruction(rho, dimensions, c):
    if c == 0:
        raise ValueError("A zero reference correlation cannot reveal the hidden sector.")
    a, b = dimensions
    joint = np.kron(rho, relation_state(c))
    recovered = np.zeros_like(rho)
    for bases, hidden in (((symmetric_basis(a), symmetric_basis(b)), False),
                          ((skew_basis(a), skew_basis(b)), True)):
        for left in bases[0]:
            for right in bases[1]:
                target = np.kron(left, right)
                if hidden:
                    observable = local_product_in_target_order(
                        np.kron(left, J), np.kron(right, J), dimensions)
                    moment = np.trace(joint @ observable) / c
                else:
                    observable = local_product_in_target_order(
                        np.kron(left, np.eye(2)), np.kron(right, np.eye(2)), dimensions)
                    moment = np.trace(joint @ observable)
                recovered += moment * target / np.sum(target * target)
    return recovered


class ReferencePermissionAuditTests(unittest.TestCase):
    def test_symmetric_and_skew_tensor_sectors_span_full_real_joint_space(self):
        for a, b in ((2, 2), (3, 2), (3, 3)):
            basis = [np.kron(x, y) for xs, ys in
                     ((symmetric_basis(a), symmetric_basis(b)), (skew_basis(a), skew_basis(b)))
                     for x in xs for y in ys]
            self.assertEqual(len(basis), a * b * (a * b + 1) // 2)
            gram = np.array([[np.sum(x * y) for y in basis] for x in basis])
            np.testing.assert_array_equal(gram, np.diag(np.diag(gram)))
            self.assertTrue(np.all(np.diag(gram) > 0))
            for matrix in basis:
                np.testing.assert_array_equal(matrix, matrix.T)

    def test_reference_projection_selects_exactly_the_skew_skew_sector(self):
        for c in (-1., -.25, 0., .5, 1.):
            tau = relation_state(c)
            projected = (tau - partial_transpose(tau, (2, 2), (1,))) / 2
            np.testing.assert_array_equal(projected, c * H / 4)
            self.assertEqual(np.trace(tau @ H), c)

    def test_full_enlarged_local_readout_has_rank_nine_or_ten_on_target_family(self):
        basis = symmetric_basis(4)
        for c in (0., .25, -1.):
            matrix = np.array([complete_product_moments(np.kron(x, relation_state(c)))
                               for x in basis]).T
            self.assertEqual(np.linalg.matrix_rank(matrix), 9 if c == 0 else 10)

    def test_known_two_rebit_reference_reconstructs_arbitrary_three_by_two_real_state(self):
        rng = np.random.default_rng(187)
        for dimensions in ((2, 2), (3, 2)):
            dimension = int(np.prod(dimensions))
            raw = rng.normal(size=(dimension, dimension))
            rho = raw @ raw.T
            rho /= np.trace(rho)
            for c in (.2, -.5, 1.):
                np.testing.assert_allclose(assisted_reconstruction(rho, dimensions, c), rho, atol=3e-16)

    def test_lifted_skew_observables_are_legal_real_binary_measurements(self):
        for dimension in (2, 3, 4):
            for matrix in skew_basis(dimension):
                observable = np.kron(matrix, J)
                np.testing.assert_array_equal(observable, observable.T)
                effects = [(np.eye(2 * dimension) + sign * observable) / 2 for sign in (-1, 1)]
                np.testing.assert_array_equal(sum(effects), np.eye(2 * dimension))
                for effect in effects:
                    self.assertGreaterEqual(np.linalg.eigvalsh(effect).min(), 0)

    def test_real_separable_reference_and_shared_classical_randomness_do_not_activate(self):
        rng = np.random.default_rng(1871)
        reference = np.zeros((4, 4))
        for _ in range(5):
            x, y = rng.normal(size=(2, 2))
            x /= np.linalg.norm(x)
            y /= np.linalg.norm(y)
            reference += np.kron(np.outer(x, x), np.outer(y, y)) / 5
        np.testing.assert_array_equal(reference, partial_transpose(reference, (2, 2), (1,)))
        hidden = np.kron(H / 2, reference)
        np.testing.assert_allclose(complete_product_moments(hidden), 0, atol=1e-16)

    def test_reference_independent_of_target_gives_product_not_sum_of_relations(self):
        for q, c in ((.5, .25), (-1., .5), (.75, -1.)):
            joint = np.kron(relation_state(q), relation_state(c))
            observable = local_product_in_target_order(H, H)
            self.assertEqual(np.trace(joint @ observable), q * c)
            for a in (-1, 1):
                for b in (-1, 1):
                    effect = local_product_in_target_order((np.eye(4) + a * H) / 2,
                                                           (np.eye(4) + b * H) / 2)
                    self.assertEqual(np.trace(joint @ effect), (1 + a * b * q * c) / 4)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ReferencePermissionAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 187,
        "scope": "Ordinary finite real quantum theory with both target Hilbert dimensions >= 2, all local real effects, known fixed reference independent of target",
        "complete_bipartite_assisted_tomography_iff_reference_skew_skew_sector_nonzero": True,
        "real_separable_reference_is_sufficient": False,
        "known_nonzero_two_rebit_relation_reference_works_for_all_finite_target_dimensions": True,
        "restores_unassisted_local_tomography": False,
        "reference_origin_or_calibration_is_free": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("reference_permission_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
