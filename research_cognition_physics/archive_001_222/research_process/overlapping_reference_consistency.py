"""Round 81: exact compatibility of three overlapping YY reference marginals."""

import argparse
import json
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from independent_source_alignment import keep_systems, pair_reference
from role_symmetry_and_swap import pauli_word


SIGNS = tuple(product((-1, 1), repeat=2))
OBSERVABLES = tuple(pauli_word(word) for word in ("YYI", "IYY", "YIY"))


def triangle_weights(q12, q23, q13):
    return {(s, t): (1+s*q12+t*q23+s*t*q13)/4 for s, t in SIGNS}


def triangle_state(q12, q23, q13):
    return (np.eye(8)+sum(float(q)*o for q, o in zip((q12, q23, q13), OBSERVABLES)))/8


def joint_projectors():
    first, second, _ = OBSERVABLES
    return {(s, t): (np.eye(8)+s*first)@(np.eye(8)+t*second)/4 for s, t in SIGNS}


def is_compatible(q12, q23, q13):
    return all(w >= 0 for w in triangle_weights(q12, q23, q13).values())


class OverlappingReferenceConsistencyTests(unittest.TestCase):
    def test_three_joint_observables_commute_and_their_product_is_identity(self):
        for first, second in product(OBSERVABLES, repeat=2):
            np.testing.assert_array_equal(first@second, second@first)
        np.testing.assert_array_equal(OBSERVABLES[0]@OBSERVABLES[1], OBSERVABLES[2])
        np.testing.assert_array_equal(OBSERVABLES[0]@OBSERVABLES[1]@OBSERVABLES[2], np.eye(8))

    def test_joint_projectors_are_complete_real_orthogonal_rank_two_sectors(self):
        projectors = joint_projectors()
        np.testing.assert_array_equal(sum(projectors.values()), np.eye(8))
        for key, p in projectors.items():
            np.testing.assert_array_equal(p.imag, np.zeros((8, 8)))
            np.testing.assert_array_equal(p@p, p)
            self.assertEqual(np.trace(p), 2)
            for other, r in projectors.items():
                if key != other:
                    np.testing.assert_array_equal(p@r, np.zeros((8, 8)))

    def test_weight_formula_is_necessary_for_arbitrary_real_and_complex_wholes(self):
        rng = np.random.default_rng(81)
        for complex_entries in (False, True):
            for _ in range(8):
                raw = rng.normal(size=(8, 8))
                if complex_entries:
                    raw = raw + 1j*rng.normal(size=(8, 8))
                whole = raw@raw.conj().T
                whole /= np.trace(whole)
                q = [np.trace(whole@o).real for o in OBSERVABLES]
                expected = triangle_weights(*q)
                for key, p in joint_projectors().items():
                    self.assertAlmostEqual(expected[key], np.trace(whole@p).real)
                    self.assertGreaterEqual(expected[key], 0)

    def test_every_rational_tetrahedron_point_has_the_required_pair_marginals(self):
        projectors = joint_projectors()
        for numerators in ((1,1,1,1), (0,0,0,1), (1,2,3,4), (0,3,1,2)):
            weights = {key: Fraction(n, sum(numerators)) for key, n in zip(SIGNS, numerators)}
            q = (sum(s*w for (s,t),w in weights.items()),
                 sum(t*w for (s,t),w in weights.items()),
                 sum(s*t*w for (s,t),w in weights.items()))
            self.assertEqual(triangle_weights(*q), weights)
            whole = triangle_state(*q)
            expected = sum(float(weights[key])*p/2 for key, p in projectors.items())
            np.testing.assert_allclose(whole, expected, atol=3e-17)
            for sites, value in zip(((0,1),(1,2),(0,2)), q):
                np.testing.assert_allclose(keep_systems(whole, sites, 3),
                                           pair_reference(float(value)), atol=0)

    def test_valid_agreeing_pair_marginals_can_fail_global_compatibility(self):
        pair = pair_reference(-1)
        self.assertGreaterEqual(np.linalg.eigvalsh(pair).min(), 0)
        for site in ((0,), (1,)):
            np.testing.assert_array_equal(keep_systems(pair, site, 2), np.eye(2)/2)
        weights = triangle_weights(Fraction(-1), Fraction(-1), Fraction(-1))
        self.assertEqual(min(weights.values()), Fraction(-1,2))
        self.assertFalse(is_compatible(-1,-1,-1))
        self.assertAlmostEqual(np.linalg.eigvalsh(triangle_state(-1,-1,-1)).min(), -.25)

    def test_equal_anticorrelation_has_the_exact_one_third_ceiling(self):
        for r in (Fraction(0), Fraction(1,4), Fraction(1,3), Fraction(1,2), Fraction(1)):
            self.assertEqual(is_compatible(-r,-r,-r), r <= Fraction(1,3))
        self.assertEqual(min(triangle_weights(*(-Fraction(1,3),)*3).values()), 0)

    def test_local_real_sign_changes_preserve_compatibility(self):
        q = (Fraction(1,5), Fraction(-2,5), Fraction(1,10))
        initial = sorted(triangle_weights(*q).values())
        for a,b,c in product((-1,1), repeat=3):
            transformed = (a*b*q[0], b*c*q[1], a*c*q[2])
            self.assertEqual(sorted(triangle_weights(*transformed).values()), initial)

    def test_compatible_pair_marginals_need_not_identify_the_whole(self):
        wholes = [(np.eye(8)+sign*pauli_word("XXX"))/8 for sign in (-1,1)]
        for whole in wholes:
            self.assertGreaterEqual(np.linalg.eigvalsh(whole).min(), 0)
            for sites in ((0,1),(1,2),(0,2)):
                np.testing.assert_array_equal(keep_systems(whole, sites, 3), np.eye(4)/4)
        self.assertEqual(np.abs(np.linalg.eigvalsh(wholes[0]-wholes[1])).sum()/2, 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(OverlappingReferenceConsistencyTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 81,
        "specified_pair_marginals": "(I+q_ij YY)/4 for AB, BC and AC",
        "necessary_and_sufficient": "1+s*q_AB+t*q_BC+s*t*q_AC >= 0 for all s,t in {-1,+1}",
        "pair_positivity_and_matching_single_marginals_are_sufficient": False,
        "all_minus_one_counterexample_minimum_weight_exact": "-1/2",
        "equal_anticorrelation_maximum_exact": "1/3",
        "compatible_pair_marginals_determine_a_unique_whole": False,
        "complex_wholes_can_evade_this_fixed_observable_constraint": False,
        "quantum_theory_derived_from_cognition": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}
    }
    if args.write_results:
        Path(__file__).with_name("overlapping_reference_consistency_results.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
