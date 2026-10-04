"""Round 509: current cut-membership records under the unchanged tree flips.

This audits locally readable membership, not all possible macroscopic objects.
Approximate records obey an explicit interaction-budget/generator-error tradeoff.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import itertools
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import coherent_graph_mean_obstruction as old
import branching_tree_distance_audit as trees
from distributed_role_reader import exponential, opnorm

HERE = Path(__file__).resolve().parent
TARGET = HERE/'dynamic_membership_record_audit_results.json'
OBS = {}


def component(tree, n, start, deleted):
    adj = trees.adjacency(tree-set(deleted), n)
    seen, todo = {start}, [start]
    while todo:
        for neighbor in adj[todo.pop()]:
            if neighbor not in seen:
                seen.add(neighbor)
                todo.append(neighbor)
    return seen


def membership(tree, n, b, c):
    assert trees.edge(b, c) in tree
    side = component(tree, n, b, {trees.edge(b, c)})
    return tuple(int(v not in side) for v in range(n))


def central_flip(tree, a, b, c, d):
    changed = {trees.edge(a, b), trees.edge(c, d), trees.edge(a, c), trees.edge(b, d)}
    return frozenset(tree.symmetric_difference(changed))


def witness(i):
    assert i >= 2
    n = 2*i+2
    tree = {trees.edge(j, j+1) for j in range(i-1)}
    tree |= {trees.edge(0, i), trees.edge(i-1, i+1)}
    tree |= {trees.edge(j, i+2+j) for j in range(i)}
    tree = frozenset(tree)
    b, c = i//2-1, i//2
    a = i if b == 0 else b-1
    d = i+1 if c == i-1 else c+1
    flipped = central_flip(tree, a, b, c, d)
    first, second = membership(tree, n, b, c), membership(flipped, n, b, c)
    changed = [v for v in range(n) if first[v] != second[v]]
    return tree, flipped, (a, b, c, d), changed


def embed(local, sites, qubits):
    indices = np.arange(2**qubits, dtype=np.int64)
    selected = np.zeros(2**qubits, dtype=np.int64)
    for site in sites:
        selected = 2*selected+((indices >> (qubits-1-site)) & 1)
    outside = sum(1 << (qubits-1-j) for j in range(qubits) if j not in sites)
    rest = indices & outside
    return np.where(rest[:, None] == rest[None, :], local[selected[:, None], selected[None, :]], 0)


def two_graph_code(changed):
    # Minimal two-configuration witness: graph qubit followed by changed records.
    # Unchanged data/records may be tensor factors; this is not a new tree dynamics.
    w = np.zeros((2**(changed+1), 2), dtype=np.int64)
    w[0, 0] = w[-1, 1] = 1
    return w


class Audit(unittest.TestCase):
    def test_01_unbounded_record_change_for_a_four_owner_flip(self):
        rows = []
        for i in (2, 3, 4, 8, 16, 64):
            tree, flipped, (a, b, c, d), changed = witness(i)
            n = 2*i+2
            left = component(tree, n, a, {trees.edge(a, b)})
            right = component(tree, n, d, {trees.edge(c, d)})
            self.assertEqual(set(changed), left | right)
            self.assertEqual(len(changed), n-4)
            self.assertEqual(len(tree.symmetric_difference(flipped)), 4)
            before, after = trees.adjacency(tree, n), trees.adjacency(flipped, n)
            owners = {v for v in range(n) if before[v] != after[v]}
            self.assertEqual(owners, {a, b, c, d})
            self.assertEqual(trees.tree_flips(tree, n)[flipped], 1)
            self.assertEqual([len(x) for x in before], [3]*i+[1]*(i+2))
            self.assertEqual([len(x) for x in before], [len(x) for x in after])
            rows.append(dict(I=i, N=n, central_edge=[b, c], address_owners_changed=4,
                member_records_changed=len(changed), changed_edges=4,
                unchanged_record_owners=sorted(set(range(n))-set(changed)),
                original_flip_matrix_element=1))
        OBS['extensive_membership_witnesses'] = rows

    def test_02_all_small_central_flips_and_nonpersistent_interface(self):
        rows = []
        for i in (2, 3, 4):
            family = old.previous.ensemble(i)[0]
            n, b, c = 2*i+2, 0, 1
            active = flips = missing_after_other_flip = maximum = 0
            values = set()
            for tree in family:
                if trees.edge(b, c) not in tree:
                    continue
                active += 1
                adj = trees.adjacency(tree, n)
                before = membership(tree, n, b, c)
                for a in sorted(adj[b]-{c}):
                    for d in sorted(adj[c]-{b}):
                        other = central_flip(tree, a, b, c, d)
                        after = membership(other, n, b, c)
                        changed = {v for v in range(n) if before[v] != after[v]}
                        expected = component(tree, n, a, {trees.edge(a, b)})
                        expected |= component(tree, n, d, {trees.edge(c, d)})
                        self.assertEqual(changed, expected)
                        self.assertLessEqual(len(changed), n-4)
                        maximum = max(maximum, len(changed))
                        values.add(len(changed))
                        flips += 1
                missing_after_other_flip += sum(
                    trees.edge(b, c) not in other for other in trees.tree_flips(tree, n))
            self.assertEqual(maximum, n-4)
            if i >= 3:
                self.assertGreater(missing_after_other_flip, 0)
            rows.append(dict(I=i, N=n, all_trees=len(family), interface_present_trees=active,
                central_directed_flips=flips, possible_record_changes=sorted(values),
                maximum_record_changes=maximum,
                other_flips_removing_selected_edge=missing_after_other_flip))
        OBS['exhaustive_tree_certificate'] = rows

    def test_03_exact_support_selection_and_generator_obstruction(self):
        changed, support = 6, 4
        w = two_graph_code(changed)
        rng = np.random.default_rng(509)
        tested, maximum = 0, 0.0
        for count in range(support+1):
            for subset in itertools.combinations(range(changed), count):
                sites = [0]+[v+1 for v in subset]
                matrix = rng.integers(-3, 4, size=(2**len(sites),)*2)
                matrix = matrix+matrix.T
                full = embed(matrix, sites, changed+1)
                value = int((w.T @ full @ w)[1, 0])
                self.assertEqual(value, 0)
                maximum = max(maximum, abs(value))
                tested += 1
        x = np.array([[0, 1], [1, 0]], dtype=np.int64)
        low = sum(embed(x, [j], changed+1) for j in range(changed+1))
        delta = opnorm(low @ w-w @ x)
        self.assertGreaterEqual(delta, 1)
        dressed = np.array([[1]], dtype=np.int64)
        for _ in range(changed+1):
            dressed = np.kron(dressed, x)
        self.assertTrue(np.array_equal(dressed @ w, w @ x))
        OBS['exact_support'] = dict(changed_records=changed, maximum_records_per_term=support,
            tested_integer_terms=tested, maximum_offdiagonal_encoded_element=int(maximum),
            required_target_element=1, low_body_generator_error=old.short(delta),
            high_body_single_witness_intertwining_exact=True,
            arbitrary_auxiliary_factors_covered_by_orthogonality_proof=True)

    def test_04_correlated_approximate_records_and_budget_tradeoff(self):
        changed = 6
        x = np.array([[0, 1], [1, 0]])
        graph_flip = embed(x, [0], changed+1)
        rows = []
        for eps in (Q(1, 10000), Q(1, 100), Q(1, 4)):
            v = np.zeros((2**(changed+1), 2))
            v[0, 0] = v[-1, 1] = np.sqrt(1-float(eps))
            v[2**changed-1, 0] = v[2**changed, 1] = np.sqrt(float(eps))
            self.assertLess(opnorm(v.T @ v-np.eye(2)), 1e-14)
            # Errors on all records are perfectly correlated, not independent.
            for member in range(changed):
                bits = (np.arange(len(v)) >> (changed-1-member)) & 1
                self.assertAlmostEqual(float(np.sum(v[bits == 1, 0]**2)), float(eps))
                self.assertAlmostEqual(float(np.sum(v[bits == 0, 1]**2)), float(eps))
            element = float((v.T @ graph_flip @ v)[1, 0])
            delta = opnorm(graph_flip @ v-v @ x)
            self.assertAlmostEqual(element, 2*np.sqrt(float(eps*(1-eps))))
            self.assertLessEqual(element, 2*np.sqrt(float(eps))+1e-14)
            self.assertGreaterEqual(delta+2*np.sqrt(float(eps)), 1-1e-14)
            # Matching the compressed generator alone still leaves real leakage.
            rescaled = graph_flip/element
            compressed = v.T @ rescaled @ v
            leakage = opnorm(rescaled @ v-v @ compressed)
            self.assertLess(opnorm(compressed-x), 1e-13)
            self.assertGreater(leakage, 0.5)
            rows.append(dict(per_record_error=str(eps), total_term_norm_budget=1,
                encoded_flip_element=old.short(element), generator_error=old.short(delta),
                selection_upper=old.short(2*np.sqrt(float(eps))),
                compressed_matching_requires_norm=old.short(1/element),
                compressed_match_still_has_leakage=old.short(leakage)))
        certificates = []
        for budget, delta in ((1, Q(0)), (10, Q(1, 10)), (100, Q(1, 2))):
            floor = (1-delta)**2/(4*budget**2)
            certificates.append(dict(total_term_norm_budget=budget,
                generator_error_upper=str(delta), necessary_record_error_floor=str(floor)))
        OBS['approximate_records'] = dict(correlated_cases=rows, rational_certificates=certificates,
            bound='abs(kappa) <= delta + 2 sqrt(epsilon) J',
            J_is_declared_sum_of_term_norms=True, generator_error_is_not_channel_error=True)

    def test_05_finite_time_witness_is_not_a_size_uniform_no_go(self):
        changed = 6
        x = np.array([[0, 1], [1, 0]], dtype=np.int64)
        w = two_graph_code(changed)
        low = sum(embed(x, [j], changed+1) for j in range(changed+1))
        duration = Q(1, 128)
        remainder = Q(((changed+1)**2+1), 2)*duration**2
        certificate = duration-remainder
        actual = exponential(float(duration)*low) @ w
        target = w @ exponential(float(duration)*x)
        coordinate_difference = abs((actual-target)[-1, 0])
        self.assertGreaterEqual(coordinate_difference, float(certificate)-1e-13)
        self.assertGreater(certificate, duration/2)
        OBS['finite_time'] = dict(duration=str(duration), low_body_norm_bound=changed+1,
            target_norm_bound=1, vector_coordinate_error=old.short(coordinate_difference),
            certified_amplitude_difference_lower=str(certificate),
            scale_uniform_time_claimed=False, arbitrary_macro_approximation_excluded=False)

    def test_06_parallel_endpoints_do_not_give_a_hamming_latency_bound(self):
        x = np.array([[0, 1], [1, 0]], dtype=np.int64)
        rows = []
        for changed in (2, 4, 6, 8):
            h = sum(embed(x, [j], changed) for j in range(changed))
            state = np.zeros(2**changed, dtype=np.int64)
            state[0] = 1
            for order in range(changed):
                self.assertEqual(int(state[-1]), 0)
                state = h @ state
            self.assertEqual(int(state[-1]), int(np.prod(np.arange(1, changed+1))))
            evolution = exponential(np.pi/2*h)
            self.assertAlmostEqual(abs(evolution[-1, 0]), 1)
            rows.append(dict(changed_records=changed, per_term_record_support=1,
                total_term_norm_budget=changed, first_nonzero_power=changed,
                complete_parallel_flip_time_over_pi='1/2'))
        OBS['finite_delay_boundary'] = dict(cases=rows, endpoint_update_not_full_tree_intertwining=True,
            graph_conditional_routing_or_spatial_locality_not_derived=True,
            hamming_distance_alone_is_not_physical_latency=True)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(round=509, scientific_baseline_round=508,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'coherent_graph_mean_obstruction.py', 'branching_tree_distance_audit.py',
            'distributed_role_reader.py', 'role_covariant_graph_obstruction.py',
            'research_note_444.md', 'research_note_508.md')},
        scope=dict(unchanged_original_NNI_witness=True,
            exact_current_local_membership_has_extensive_support_obstruction=True,
            approximate_record_generator_budget_tradeoff=True,
            arbitrary_reference_at_isometric_interface=True,
            selected_interface_is_declared_input=True,
            edge_present_sector_claimed_invariant=False,
            hamming_distance_used_as_time_lower_bound=False,
            internal_membership_source_generated=False,
            all_macroscopic_membership_or_geometry_excluded=False,
            dimension_three_generated=False, full_GR_goal_completed=False,
            phase_closure_triggered=False), observations=OBS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    data = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(data, ensure_ascii=False, indent=2))
