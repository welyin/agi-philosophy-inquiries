"""Round 207: finite real-reference growth on a declared tree.

Source preparation, joining and readout use previously implemented real gates.
Density matrices trace out spent systems only for calculation; the resource
account retains them inside a finite whole. This is not autonomous recycling.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import unittest

import numpy as np

from amplified_source_alignment import pointer_branch
from bilocal_record_tomography import embed_operator
from bipartite_composition import interaction
from certified_intervals import Interval as I, SCALE
from independent_source_alignment import keep_systems, pair_reference, real_x_with_reset
from one_bit_real_network import visibility_interval
from operational_effect_closure import majority_error
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_Y, PAULI_Z, branch_kraus, rotation_unitary
from reference_information_cost import target_intervals
from reference_maintenance_budget import maximum_age
from reusable_network_reference import local_record_kraus, direction_vectors


def tensor(items):
    result = np.ones((1, 1), dtype=complex)
    for item in items:
        result = np.kron(result, item)
    return result


def validate_tree(parents):
    if not parents or parents[0] != -1:
        raise ValueError("Root has parent -1.")
    if any(not isinstance(p, int) or not 0 <= p < v for v, p in enumerate(parents[1:], 1)):
        raise ValueError("Each new vertex attaches to a previously present vertex.")


def orientation_distribution(parents, correlations):
    validate_tree(parents)
    if len(correlations) != len(parents)-1 or any(not -1 <= g <= 1 for g in correlations):
        raise ValueError("One allowed correlation per tree edge required.")
    return {s: .5*math.prod((1+g*s[v]*s[parents[v]])/2
                           for v, g in enumerate(correlations, 1))
            for s in product((-1, 1), repeat=len(parents))}


def tree_state(parents, correlations):
    return sum(p*tensor((IDENTITY+s_i*PAULI_Y)/2 for s_i in s)
               for s, p in orientation_distribution(parents, correlations).items())


def path_vertices(parents, a, b):
    validate_tree(parents)
    if not 0 <= a < len(parents) or not 0 <= b < len(parents):
        raise ValueError("Vertices outside tree.")
    first, second = [], []
    while a != -1:
        first.append(a)
        a = parents[a]
    while b not in first:
        second.append(b)
        b = parents[b]
    return first[:first.index(b)+1]+second[::-1]


def path_correlation(parents, correlations, a, b):
    path = path_vertices(parents, a, b)
    return math.prod(correlations[max(u, v)-1] for u, v in zip(path, path[1:]))


def prepared_root_and_link():
    """Keep the purifier in accounting; the returned states are marginals."""
    zero = np.diag([1., 0.])
    bell_gate = interaction(math.pi/2)
    bell = bell_gate @ np.kron(zero, zero) @ bell_gate.T
    mixed = keep_systems(bell, (0,), 2)
    link_gate = interaction(-math.pi/2)
    link = link_gate @ np.kron(mixed, zero) @ link_gate.T
    return mixed, link


def append_records(old, parent, count, contrast=1.):
    """Old systems followed by link half and child; retain all pointer records.

    old may include spectators. Preservation requires its parent-Y commutator
    to vanish. The full old matrix need not be a product or a scalar mixture.
    """
    majority_error(count, .1)
    n = int(round(math.log2(len(old))))
    if old.shape != (2**n, 2**n) or not 0 <= parent < n:
        raise ValueError("Qubit matrix and existing parent required.")
    gate = embed_operator(interaction(math.pi/2), (parent, n), n+2)
    initial = gate @ np.kron(old, pair_reference()) @ gate.T
    states = {(): initial}
    for _ in range(count):
        states = {r+(2*outcome-1,): pointer_branch(rho, n, n+2, outcome, contrast)
                  for r, rho in states.items() for outcome in (0, 1)}
    keep = tuple(range(n))+(n+1,)
    result = {}
    for record, rho in states.items():
        retained = keep_systems(rho, keep, n+2)
        result[record] = retained if sum(record) > 0 else real_x_with_reset(retained, n, n+1)
    return result


def resources(n, count):
    if not isinstance(n, int) or n < 2:
        raise ValueError("At least two agents required.")
    majority_error(count, .1)
    e = n-1
    return {
        "agents": n, "reads_per_join": count, "independent_known_link_pairs": e,
        "source_pure_rebits_including_purifiers": 2+3*e,
        "source_YX_flows": 1+2*e, "retained_source_purifiers": n,
        "active_reference_rebits": n, "spent_link_halves": e,
        "fresh_pure_readout_pointers": count*e,
        "fresh_correction_helpers_upper": e,
        "pure_rebits_total_upper_excluding_detector_and_classical_registers": 2+(count+4)*e,
        "alignment_old_reads": count*e,
        "alignment_YX_flows_upper": (count+2)*e,
        "source_plus_alignment_YX_flows_upper": 1+(count+4)*e,
        "pointer_copy_local_rotations": 3*count*e,
        "majority_sign_messages_bits": e,
        "raw_alignment_record_bits_if_archived": count*e,
        "link_system_endpoint_deliveries": 2*e,
        "star_root_local_YX_flows": (count+1)*e,
        "star_root_reference_shared_wire_YX_flows": e,
        "excluded_from_numerical_prices": ["detector environment and classical registers",
            "majority computation and clock", "delivery distance, latency and storage duration",
            "source/pointer preparation energy and eventual recovery", "subsequent target tasks"]}


@lru_cache(maxsize=None)
def required_reads(power, threshold=F(999, 1000), whole=False):
    """Smallest odd count IN THIS majority family, with certified inequalities."""
    if power < 1 or not 0 < threshold < 1:
        raise ValueError("Positive exponent and interior quality required.")
    needed = I.exact(threshold)
    for count in range(1, 100, 2):
        g = visibility_interval(count)
        quality = ((1+g)/2)**power if whole else g**power
        if quality.lo >= needed.hi:
            return count
        if quality.hi >= needed.lo:
            raise ArithmeticError("Undecided interval boundary.")
    raise ArithmeticError("Count outside finite search range.")


def complex_star_marginal(n, g):
    return tensor([(IDENTITY+PAULI_Y)/2]+[(IDENTITY+g*PAULI_Y)/2]*(n-1))


def storage_channel_all(rho, n, q):
    out = rho.copy()
    for site in range(n):
        z = embed_operator(PAULI_Z, (site,), n)
        out = (1-q)*out+q*z@out@z
    return out


class TreeReferenceGrowthTests(unittest.TestCase):
    def test_perfect_known_sources_have_explicit_real_preparation(self):
        root, link = prepared_root_and_link()
        np.testing.assert_allclose(root, IDENTITY/2, atol=3e-16)
        np.testing.assert_allclose(link, pair_reference(), atol=3e-16)

    def test_actual_pointer_circuit_appends_to_root_and_both_tree_shapes(self):
        for parents in ((-1, 0), (-1, 0, 0, 0), (-1, 0, 1, 2)):
            state = IDENTITY/2
            gs = []
            for child, parent in enumerate(parents[1:], 1):
                count = 3 if child == 1 else 1
                eta = .73
                state = sum(append_records(state, parent, count, eta).values())
                gs.append(1-2*majority_error(count, (1-ALPHA*eta)/2))
                np.testing.assert_allclose(state, tree_state(parents[:child+1], gs), atol=2e-15)

    def test_every_raw_record_preserves_old_correlations_with_spectators(self):
        # A genuinely correlated real state, block diagonal in parent Y.
        p_plus, p_minus = (IDENTITY+PAULI_Y)/2, (IDENTITY-PAULI_Y)/2
        old = .5*(np.kron(p_plus, (IDENTITY+.6*PAULI_Y+.3*PAULI_Z)/2)
                  +np.kron(p_minus, (IDENTITY-.6*PAULI_Y+.3*PAULI_Z)/2))
        for record, rho in append_records(old, 0, 3, .7).items():
            probability = .5*sum(math.prod((1+s*ALPHA*.7*r)/2 for r in record) for s in (-1, 1))
            np.testing.assert_allclose(keep_systems(rho, (0, 1), 3), probability*old, atol=7e-16)

    def test_arbitrary_off_sector_reference_is_not_claimed_preserved(self):
        zero = np.diag([1., 0.])
        joined = sum(append_records(zero, 0, 1).values())
        np.testing.assert_allclose(keep_systems(joined, (0,), 2), IDENTITY/2, atol=4e-16)
        self.assertGreater(np.linalg.norm(keep_systems(joined, (0,), 2)-zero), .7)

    def test_pair_marginals_are_path_products(self):
        parents, gs = (-1, 0, 0, 1, 2), [.7, .6, .8, .9]
        state = tree_state(parents, gs)
        for a in range(5):
            for b in range(a+1, 5):
                np.testing.assert_allclose(keep_systems(state, (a, b), 5),
                                          pair_reference(path_correlation(parents, gs, a, b)), atol=3e-16)

    def test_global_distance_is_edge_error_union_not_worst_pair_distance(self):
        for parents in ((-1, 0, 0, 0), (-1, 0, 1, 2)):
            gs = [.7, .6, .8]
            state, ideal = tree_state(parents, gs), tree_state(parents, [1.]*3)
            actual = np.abs(np.linalg.eigvalsh(state-ideal)).sum()/2
            self.assertAlmostEqual(actual, 1-math.prod((1+g)/2 for g in gs), places=14)
            self.assertGreaterEqual(np.linalg.eigvalsh(state).min(), -3e-16)
            np.testing.assert_allclose(state.imag, 0, atol=1e-16)

    def test_resource_counts_follow_appending_one_leaf(self):
        for n, m in ((2, 1), (4, 5), (1000, 7)):
            row, grown = resources(n, m), resources(n+1, m)
            self.assertEqual(grown["alignment_old_reads"]-row["alignment_old_reads"], m)
            self.assertEqual(grown["source_plus_alignment_YX_flows_upper"]-row["source_plus_alignment_YX_flows_upper"], m+4)
            self.assertEqual(row["pure_rebits_total_upper_excluding_detector_and_classical_registers"],
                             n+row["retained_source_purifiers"]+row["spent_link_halves"]
                             +row["fresh_pure_readout_pointers"]+row["fresh_correction_helpers_upper"])

    def test_certified_minimal_counts_in_declared_majority_family(self):
        for power, whole in ((2, False), (999, False), (999, True), (999999, True)):
            m = required_reads(power, whole=whole)
            threshold = I.exact(F(999, 1000))
            g = visibility_interval(m)
            quality = ((1+g)/2)**power if whole else g**power
            self.assertGreaterEqual(quality.lo, threshold.hi)
            if m > 1:
                prev = visibility_interval(m-2)
                old = ((1+prev)/2)**power if whole else prev**power
                self.assertLess(old.hi, threshold.lo)

    def test_complex_dither_really_uses_old_read_probabilities(self):
        zero = np.diag([1., 0.])
        p_wrong = sum(np.trace(k@zero@k.conj().T).real for k in branch_kraus(0., 0, 1.))
        self.assertAlmostEqual(p_wrong, (1-ALPHA)/2, places=15)
        phase = np.diag([1., 1j])
        state = phase@rotation_unitary(math.pi/2)@zero@rotation_unitary(math.pi/2).T@phase.conj().T
        np.testing.assert_allclose(state, (IDENTITY+PAULI_Y)/2, atol=3e-16)
        error = majority_error(5, p_wrong)
        noisy = (1-error)*state+error*PAULI_Z@state@PAULI_Z
        np.testing.assert_allclose(noisy, (IDENTITY+(1-2*error)*PAULI_Y)/2, atol=3e-16)

    def test_complex_private_star_matches_all_real_effects(self):
        rng = np.random.default_rng(207)
        n, g = 4, .83
        real = tree_state((-1, 0, 0, 0), [g]*3)
        complex_state = complex_star_marginal(n, g)
        np.testing.assert_allclose((complex_state+complex_state.conj())/2, real, atol=3e-16)
        for _ in range(4):
            k = rng.normal(size=(16, 16))
            effect = k.T@k
            self.assertAlmostEqual(np.trace(effect@real).real, np.trace(effect@complex_state).real, places=12)

    def test_matched_storage_decay_and_real_record_equivalence(self):
        n, g, q = 4, .83, .017
        real = tree_state((-1, 0, 0, 0), [g]*3)
        complex_state = complex_star_marginal(n, g)
        for age in range(4):
            np.testing.assert_allclose(complex_state.real, real, atol=5e-16)
            np.testing.assert_allclose(keep_systems(real, (1, 2), n),
                                      pair_reference(g*g*(1-2*q)**(2*age)), atol=5e-16)
            real = storage_channel_all(real, n, q)
            complex_state = storage_channel_all(complex_state, n, q)

    def test_old_task_kraus_stays_in_preserved_reference_sector(self):
        generator = np.kron(PAULI_Y, IDENTITY)
        for direction in direction_vectors(1, 4):
            for k in local_record_kraus(tuple(direction), 1, 5):
                np.testing.assert_allclose(k@generator, generator@k, atol=5e-16)

    def test_invalid_tree_or_counts_are_rejected(self):
        for tree in ((0,), (-1, 1), (-1, 0, 3)):
            with self.assertRaises(ValueError):
                validate_tree(tree)
        with self.assertRaises(ValueError):
            resources(3, 2)


def report():
    gamma = visibility_interval(5)
    g = sum(gamma.floats())/2
    rows = []
    for n in (10, 1000, 1000000):
        rows.append({"n": n, "fixed_m5_star_min_pair_g": g*g,
            "fixed_m5_whole_trace_distance_diagnostic": -math.expm1((n-1)*math.log1p(-(1-g)/2)),
            "star_pair_reads_for_g_at_least_0_999": required_reads(2),
            "path_pair_reads_for_g_at_least_0_999": required_reads(n-1),
            "balanced_tree_pair_reads_sufficient_for_g_at_least_0_999": required_reads(2*(n.bit_length()-1)),
            "whole_reads_for_distance_at_most_0_001": required_reads(n-1, whole=True)})
    return {"round": 207, "scope": "known blank reference attachment; finite internal stocks; ideal gates and logical storage periods",
        "tree_pair_correlation": "product of edge gamma_e along unique path",
        "whole_trace_distance": "1-product_e((1+gamma_e)/2)",
        "old_joint_state_preserved_condition": "commutes with parent Y, including correlated spectators; every raw alignment record",
        "arbitrary_unknown_complex_groups_fused": False,
        "gamma_5_interval": gamma.floats(), "scale_rows": rows,
        "resources_n1000_m5": resources(1000, 5),
        "star_q1e_minus6_worst_pair_maximum_age": maximum_age(F(1, 10**6), gamma**2),
        "old_matched_target_g_interval": target_intervals()[2].floats(),
        "complex_comparison": "declared local phase S; matched private star has identical real-protocol record laws; no total-price optimum claimed",
        "not_yet_closed": ["autonomous reset/control/energy", "continuous-time setup noise", "bounded-degree online growth at unknown final size", "universal complexity lower bound"],
        "quantum_necessity_or_spontaneous_J_derived": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(TreeReferenceGrowthTests))
    if not result.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data["automated_checks"] = {"run": result.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("tree_reference_growth_results.json").write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
