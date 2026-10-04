"""Round 181: probability duality does not imply self-duality of a state cone."""

import argparse
from itertools import permutations, product
import json
from pathlib import Path
import unittest

import numpy as np


CUBE_RAYS = np.array([(1, *signs) for signs in product((-1, 1), repeat=3)], dtype=int)
DUAL_RAYS = np.array([(1, *(sign * np.eye(3, dtype=int)[i]))
                      for i in range(3) for sign in (-1, 1)])
UNIT = np.array([1, 0, 0, 0])


def cube_effect_bounds(effect):
    return effect[0] - np.abs(effect[1:]).sum(), effect[0] + np.abs(effect[1:]).sum()


def ray_fixing_constraints(rays):
    """Linear equations for A v_i = lambda_i v_i, tangent to ray fixing."""
    count, dimension = rays.shape
    rows = []
    for i, ray in enumerate(rays):
        for j in range(dimension):
            row = np.zeros(dimension * dimension + count)
            row[j * dimension:(j + 1) * dimension] = ray
            row[dimension * dimension + i] = -ray[j]
            rows.append(row)
    return np.array(rows)


class ConeDualityAuditTests(unittest.TestCase):
    def test_dual_positivity_and_effect_bounds_match_all_vertices_exactly(self):
        for coefficients in product(range(-2, 3), repeat=4):
            effect = np.array(coefficients)
            values = CUBE_RAYS @ effect
            self.assertEqual(cube_effect_bounds(effect), (values.min(), values.max()))

    def test_eight_state_rays_and_six_dual_rays_have_exact_facet_incidence(self):
        pairing = CUBE_RAYS @ DUAL_RAYS.T
        self.assertEqual(pairing.shape, (8, 6))
        self.assertTrue(np.all((pairing == 0) | (pairing == 2)))
        for j in range(6):
            facet = CUBE_RAYS[pairing[:, j] == 0]
            self.assertEqual(len(facet), 4)
            self.assertEqual(np.linalg.matrix_rank(facet), 3)
        self.assertEqual(np.linalg.matrix_rank(CUBE_RAYS), 4)
        self.assertEqual(np.linalg.matrix_rank(DUAL_RAYS), 4)

    def test_dual_bipolar_inequalities_recover_cube_cone(self):
        for point in product(range(-2, 3), repeat=4):
            expected = point[0] >= max(abs(x) for x in point[1:])
            self.assertEqual(bool(np.all(DUAL_RAYS @ point >= 0)), expected)

    def test_three_binary_readings_separate_states_with_legal_probabilities(self):
        effects = DUAL_RAYS / 2
        for j in range(0, 6, 2):
            np.testing.assert_array_equal(effects[j] + effects[j + 1], UNIT)
        probabilities = CUBE_RAYS @ effects.T
        self.assertEqual(len(np.unique(probabilities, axis=0)), 8)
        self.assertTrue(np.all((probabilities >= 0) & (probabilities <= 1)))
        self.assertEqual(np.linalg.matrix_rank(effects), 4)

    def test_read_and_reset_is_a_positive_normalized_instrument(self):
        for i in range(3):
            effects = DUAL_RAYS[2 * i:2 * i + 2] / 2
            branches = [np.outer(UNIT, effect) for effect in effects]
            for ray in CUBE_RAYS:
                outputs = [branch @ ray for branch in branches]
                for output in outputs:
                    self.assertTrue(np.all(DUAL_RAYS @ output >= 0))
                np.testing.assert_array_equal(sum(outputs), UNIT)

    def test_pure_transitivity_does_not_remove_the_dual_ray_count_obstruction(self):
        orbit = set()
        for permutation in permutations(range(3)):
            for signs in product((-1, 1), repeat=3):
                rotation = np.diag(signs) @ np.eye(3, dtype=int)[list(permutation)]
                images = CUBE_RAYS[:, 1:] @ rotation.T
                self.assertEqual(set(map(tuple, images)), set(map(tuple, CUBE_RAYS[:, 1:])))
                orbit.add(tuple(rotation @ np.ones(3, dtype=int)))
        self.assertEqual(len(orbit), 8)
        # The connected automorphism group can only scale: one-dimensional
        # ray-fixing solution space, in accord with the face-relation proof.
        equations = ray_fixing_constraints(CUBE_RAYS)
        self.assertEqual(equations.shape[1] - np.linalg.matrix_rank(equations), 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ConeDualityAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = {
        "round": 181,
        "state_cone": "t >= max(abs(x_i)) in R^4",
        "dual_cone": "a >= sum(abs(b_i))",
        "normalized_effect_condition": "sum(abs(b_i)) <= min(a,1-a)",
        "state_extreme_rays": 8,
        "dual_extreme_rays": 6,
        "cone_linearly_isomorphic_to_its_dual": False,
        "pure_state_transitivity": True,
        "ray_fixing_lie_algebra_dimension": 1,
        "all_mathematically_positive_effects_are_physically_available_by_definition": False,
        "prediction_test_duality_implies_self_duality": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("cone_duality_audit_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
