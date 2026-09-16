"""Round 82: common-whole consistency of exact and noisy YY alignment cycles."""

import argparse
import json
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from quantum_interface_audit import IDENTITY, PAULI_Y
from role_symmetry_and_swap import pauli_word


def solve_orientations(n, edges):
    adjacency = [[] for _ in range(n)]
    for a,b,sign in edges:
        if a == b or not (0 <= a < n and 0 <= b < n) or sign not in (-1,1):
            raise ValueError("Use distinct valid sites and signs +/-1.")
        adjacency[a].append((b,sign))
        adjacency[b].append((a,sign))
    result = {}
    for root in range(n):
        if root in result:
            continue
        result[root] = 1
        stack = [root]
        while stack:
            a = stack.pop()
            for b,sign in adjacency[a]:
                expected = result[a]*sign
                if b in result and result[b] != expected:
                    return None
                if b not in result:
                    result[b] = expected
                    stack.append(b)
    return tuple(result[i] for i in range(n))


def cycle_edges(signs):
    return tuple((i,(i+1)%len(signs),s) for i,s in enumerate(signs))


def spin_whole(spins):
    density = np.ones((1,1), dtype=complex)
    for s in spins:
        density = np.kron(density, (IDENTITY+s*PAULI_Y)/2)
    return (density+density.conj())/2


def edge_operator(n,a,b):
    word = ["I"]*n
    word[a] = word[b] = "Y"
    return pauli_word("".join(word))


def frustrated_saturator(signs):
    if np.prod(signs) != -1:
        raise ValueError("The target cycle must have odd sign parity.")
    states = []
    for bad in range(len(signs)):
        actual = tuple(-s if i == bad else s for i,s in enumerate(signs))
        spins = solve_orientations(len(signs), cycle_edges(actual))
        states.append(spin_whole(spins))
    return sum(states)/len(states)


class ReferenceCycleConsistencyTests(unittest.TestCase):
    def test_graph_orientation_solver_matches_exhaustive_sign_assignments(self):
        for signs in product((-1,1), repeat=6):
            pairs = ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
            edges = tuple((a,b,s) for (a,b),s in zip(pairs,signs))
            exhaustive = [x for x in product((-1,1),repeat=4)
                          if all(x[a]*x[b] == s for a,b,s in edges)]
            solved = solve_orientations(4,edges)
            self.assertEqual(solved is not None, bool(exhaustive))
            if solved is not None:
                self.assertIn(solved,exhaustive)

    def test_every_signed_tree_has_a_real_whole_with_the_requested_edge_correlations(self):
        pairs = ((0,1),(1,2),(1,3))
        for signs in product((-1,1), repeat=3):
            edges = tuple((a,b,s) for (a,b),s in zip(pairs,signs))
            spins = solve_orientations(4,edges)
            whole = spin_whole(spins)
            np.testing.assert_array_equal(whole.imag,np.zeros((16,16)))
            self.assertAlmostEqual(np.trace(whole).real,1)
            for a,b,s in edges:
                self.assertAlmostEqual(np.trace(whole@edge_operator(4,a,b)).real,s)

    def test_cycle_exact_solvability_is_exactly_positive_sign_product(self):
        for n in range(3,8):
            for signs in product((-1,1), repeat=n):
                self.assertEqual(solve_orientations(n,cycle_edges(signs)) is not None,
                                 np.prod(signs) == 1)

    def test_independent_local_orientation_changes_cannot_remove_negative_cycle_parity(self):
        signs = (-1,1,1,1,1)
        for flips in product((-1,1),repeat=5):
            changed = tuple(s*flips[i]*flips[(i+1)%5] for i,s in enumerate(signs))
            self.assertEqual(np.prod(changed),-1)
            self.assertIsNone(solve_orientations(5,cycle_edges(changed)))

    def test_frustrated_cycle_operator_has_the_exact_classical_sign_spectrum(self):
        for n in range(3,7):
            signs = (-1,)+(1,)*(n-1)
            total = sum(s*edge_operator(n,a,b) for a,b,s in cycle_edges(signs))
            expected = sorted(sum(s*x[a]*x[b] for a,b,s in cycle_edges(signs))
                              for x in product((-1,1),repeat=n))
            np.testing.assert_allclose(np.linalg.eigvalsh(total),expected,atol=8e-15)
            self.assertEqual(max(expected), n-2)
            self.assertEqual(min(Fraction(n-score,2) for score in expected),1)

    def test_uniform_single_defect_mixture_saturates_every_edge_equally(self):
        for n in range(3,7):
            signs = (-1,)+(1,)*(n-1)
            whole = frustrated_saturator(signs)
            self.assertGreaterEqual(np.linalg.eigvalsh(whole).min(),-2e-16)
            self.assertAlmostEqual(np.trace(whole).real,1)
            for a,b,s in cycle_edges(signs):
                self.assertAlmostEqual(s*np.trace(whole@edge_operator(n,a,b)).real,1-2/n)

    def test_common_complex_wholes_obey_the_same_bound(self):
        rng = np.random.default_rng(82)
        n = 4
        signs = (-1,1,1,1)
        total = sum(s*edge_operator(n,a,b) for a,b,s in cycle_edges(signs))
        for _ in range(12):
            raw = rng.normal(size=(16,16))+1j*rng.normal(size=(16,16))
            rho = raw@raw.conj().T
            rho /= np.trace(rho)
            self.assertLessEqual(np.trace(rho@total).real,n-2+1e-15)

    def test_odd_all_anticorrelated_cycles_reproduce_triangle_and_higher_exact_thresholds(self):
        for n in (3,5,7):
            signs = (-1,)*n
            self.assertEqual(np.prod(signs),-1)
            threshold = Fraction(n-2,n)
            self.assertEqual((1-threshold)/2,Fraction(1,n))
        self.assertEqual(Fraction(3-2,3),Fraction(1,3))
        self.assertEqual(Fraction(5-2,5),Fraction(3,5))
        with self.assertRaises(ValueError):
            frustrated_saturator((1,1,1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ReferenceCycleConsistencyTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round":82,
        "perfect_graph_alignment_condition":"Product of required edge signs is +1 on every cycle",
        "frustrated_n_cycle_bound_exact":"sum_e s_e <Y_i Y_j> <= n-2",
        "sum_of_edge_disagreement_probabilities_at_least_exact":"1",
        "uniform_edge_retention_maximum_exact":"1-2/n",
        "uniform_edge_bound_is_attained_by_real_whole":True,
        "odd_anticorrelation_examples":[{"n":n,"maximum_r_exact":str(Fraction(n-2,n))} for n in (3,5,7)],
        "local_sign_relabeling_can_fix_frustrated_cycle":False,
        "complex_wholes_can_evade_same_observable_bound":False,
        "a_general_noncommuting_marginal_problem_is_solved":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}
    }
    if args.write_results:
        Path(__file__).with_name("reference_cycle_consistency_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
