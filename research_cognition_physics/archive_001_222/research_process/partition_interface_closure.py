"""Round 80: regrouping and recursive interface closure of real matrix models.

Symmetric and antisymmetric matrix sectors form an associative two-sector
description. The latter is relational interface data, not a local real state.
"""

import argparse
import json
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from role_symmetry_and_swap import pauli_word


def sector_bases(d):
    if not isinstance(d, int) or d < 1:
        raise ValueError("Use a positive integer matrix dimension.")
    units = np.eye(d*d).reshape(d*d, d, d)
    symmetric = [units[i*d+i] for i in range(d)]
    skew = []
    for i in range(d):
        for j in range(i+1, d):
            symmetric.append(units[i*d+j] + units[j*d+i])
            skew.append(units[i*d+j] - units[j*d+i])
    return tuple(symmetric), tuple(skew)


def dimensions(d):
    return d*(d+1)//2, d*(d-1)//2


def combine_dimensions(first, second):
    k, l = first
    m, n = second
    return k*m+l*n, k*n+l*m


def apply(kraus, matrix):
    return sum(k @ matrix @ k.T for k in kraus)


def sector_descriptor(kraus):
    d = len(kraus[0])
    if any(np.iscomplexobj(k) and np.any(k.imag) for k in kraus):
        raise ValueError("The declared operations have real Kraus matrices.")
    return tuple(np.array([[np.sum(a * apply(kraus, b))/np.sum(a*a)
                            for b in basis] for a in basis]).reshape(len(basis), len(basis))
                 for basis in sector_bases(d))


def reconstructed_action(descriptor, matrix):
    result = np.zeros_like(matrix, dtype=np.result_type(matrix, float))
    for transform, basis in zip(descriptor, sector_bases(len(matrix))):
        coordinates = np.array([np.sum(b*matrix)/np.sum(b*b) for b in basis])
        for coefficient, b in zip(transform @ coordinates, basis):
            result += coefficient*b
    return result


def full_choi(action, d):
    result = np.zeros((d*d, d*d))
    for i, j in product(range(d), repeat=2):
        unit = np.zeros((d, d))
        unit[i, j] = 1
        result += np.kron(unit, action(unit))
    return result


def opposite_completions(d, sign):
    """Real CP/TP: depolarizing +/- (id-transpose)/(2*d^2)."""
    if sign not in (-1, 1):
        raise ValueError("Use a completion sign.")
    strength = 1/(d*d)
    action = lambda m: np.trace(m)*np.eye(d)/d + sign*strength*(m-m.T)/2
    choi = full_choi(action, d)
    values, vectors = np.linalg.eigh(choi)
    if min(values) <= 0:
        raise AssertionError("Analytic Choi lower bound was violated.")
    return tuple(np.sqrt(value)*vectors[:, i].reshape(d, d).T
                 for i, value in enumerate(values))


def one_rebit_probe(skew):
    d = len(skew)
    j = np.array([[0., 1.], [-1., 0.]])
    norm = np.linalg.norm(skew, 2)
    if norm == 0:
        raise ValueError("A nonzero antisymmetric direction is required.")
    return (np.eye(2*d) + np.kron(skew/norm, j)/2)/(2*d)


class PartitionInterfaceClosureTests(unittest.TestCase):
    def test_integer_sector_bases_cover_every_matrix_and_have_correct_transpose(self):
        for d in range(1, 6):
            sectors = sector_bases(d)
            self.assertEqual(tuple(map(len, sectors)), dimensions(d))
            basis = sectors[0] + sectors[1]
            gram = np.array([[np.sum(a*b) for b in basis] for a in basis])
            np.testing.assert_array_equal(gram, np.diag(np.diag(gram)))
            self.assertTrue(np.all(np.diag(gram) > 0))
            for sign, sector in zip((1, -1), sectors):
                for b in sector:
                    np.testing.assert_array_equal(b.T, sign*b)

    def test_regrouping_dimensions_are_associative_for_arbitrary_block_sizes(self):
        for a, b, c in product((1, 2, 3, 4), repeat=3):
            left = combine_dimensions(combine_dimensions(dimensions(a), dimensions(b)), dimensions(c))
            right = combine_dimensions(dimensions(a), combine_dimensions(dimensions(b), dimensions(c)))
            self.assertEqual(left, right)
            self.assertEqual(left, dimensions(a*b*c))

    def test_tensor_sector_rule_holds_on_all_integer_basis_products(self):
        sectors = sector_bases(2)
        counts = [0, 0]
        for pa, pb, pc in product(range(2), repeat=3):
            for a, b, c in product(sectors[pa], sectors[pb], sectors[pc]):
                left = np.kron(np.kron(a, b), c)
                right = np.kron(a, np.kron(b, c))
                np.testing.assert_array_equal(left, right)
                np.testing.assert_array_equal(left.T, (-1)**(pa+pb+pc)*left)
                counts[(pa+pb+pc) % 2] += 1
        self.assertEqual(tuple(counts), dimensions(8))

    def test_complete_descriptor_reproduces_all_matrix_units_and_selective_mass(self):
        rng = np.random.default_rng(80)
        for d in (2, 3, 4):
            kraus = tuple(rng.normal(size=(d, d))/(4*d) for _ in range(3))
            descriptor = sector_descriptor(kraus)
            for unit in np.eye(d*d).reshape(d*d, d, d):
                np.testing.assert_allclose(reconstructed_action(descriptor, unit),
                                           apply(kraus, unit), atol=2e-17)
            state = np.eye(d)/d
            self.assertAlmostEqual(np.trace(reconstructed_action(descriptor, state)),
                                   np.trace(apply(kraus, state)))

    def test_same_symmetric_action_can_have_distinct_physical_skew_actions_in_larger_blocks(self):
        for d in (2, 3, 4):
            maps = [opposite_completions(d, s) for s in (-1, 1)]
            descriptors = [sector_descriptor(k) for k in maps]
            np.testing.assert_allclose(descriptors[0][0], descriptors[1][0], atol=9e-16)
            for sign, kraus, descriptor in zip((-1, 1), maps, descriptors):
                np.testing.assert_allclose(sum(k.T@k for k in kraus), np.eye(d), atol=2e-15)
                np.testing.assert_allclose(descriptor[1], sign*np.eye(dimensions(d)[1])/d**2, atol=9e-16)

    def test_one_extra_rebit_exposes_every_skew_basis_direction_of_a_block(self):
        for d in (2, 3, 4):
            maps = [opposite_completions(d, s) for s in (-1, 1)]
            for skew in sector_bases(d)[1]:
                rho = one_rebit_probe(skew)
                self.assertGreater(np.linalg.eigvalsh(rho).min(), 0)
                self.assertAlmostEqual(np.trace(rho), 1)
                outputs = [apply(tuple(np.kron(k, np.eye(2)) for k in kraus), rho) for kraus in maps]
                gap = np.linalg.norm(outputs[1]-outputs[0])
                self.assertGreater(gap, 1e-3)
                j = np.array([[0., 1.], [-1., 0.]])
                np.testing.assert_allclose(outputs[1]-outputs[0],
                                           np.kron(skew, j)/(2*d**3), atol=5e-16)

    def test_descriptor_composition_and_mixture_preserve_both_sectors(self):
        a, b = opposite_completions(3, 1), opposite_completions(3, -1)
        da, db = sector_descriptor(a), sector_descriptor(b)
        composed = sector_descriptor(tuple(kb@ka for ka in a for kb in b))
        mixture = sector_descriptor(tuple(np.sqrt(.25)*k for k in a) + tuple(np.sqrt(.75)*k for k in b))
        for sector in (0, 1):
            np.testing.assert_allclose(composed[sector], db[sector]@da[sector], atol=6e-16)
            np.testing.assert_allclose(mixture[sector], .25*da[sector]+.75*db[sector], atol=6e-16)

    def test_real_closure_differs_from_local_interface_completeness(self):
        for a, b in product(range(2, 6), repeat=2):
            ka, la = dimensions(a)
            kb, lb = dimensions(b)
            self.assertEqual(dimensions(a*b)[0]-ka*kb, la*lb)
            self.assertGreater(la*lb, 0)
        # Compare actual matrix pairings, not just the dimension formulas.
        for local_labels, expected_rank in (("IXZ", 9), ("IXYZ", 16)):
            effects = [pauli_word(a+b) for a,b in product(local_labels,repeat=2)]
            whole_basis = [pauli_word(a+b) for a,b in product("IXYZ",repeat=2)]
            pairing = np.array([[np.trace(e@b).real/4 for b in whole_basis] for e in effects])
            self.assertEqual(np.linalg.matrix_rank(pairing),expected_rank)
        classical = [np.diag(np.eye(4)[i]) for i in range(4)]
        classical_pairing = np.array([[np.trace(e@b) for b in classical] for e in classical])
        np.testing.assert_array_equal(classical_pairing,np.eye(4))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(PartitionInterfaceClosureTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 80,
        "real_block_dimensions": [{"d": d, "symmetric": dimensions(d)[0], "skew": dimensions(d)[1]}
                                  for d in (2, 4, 8, 16)],
        "composition_rule": "(K,L)*(M,N)=(KM+LN,KN+LM)",
        "same_rule_after_regrouping_forces_complex_structure": False,
        "complete_real_operation_interface": "Symmetric-sector map plus antisymmetric-sector map",
        "one_additional_rebit_suffices_to_probe_a_missing_skew_direction": True,
        "four_dimensional_block_interface_sizes": {"symmetric_map": [10,10], "skew_map": [6,6]},
        "local_only_dimension_deficit": "L_A*L_B",
        "local_tomography_alone_selects_complex_over_classical": False,
        "scope": "Declared real Kraus framework and usual tensor products; no cognitive derivation of these inputs",
        "quantum_theory_derived_from_cognition": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}
    }
    if args.write_results:
        Path(__file__).with_name("partition_interface_closure_results.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
