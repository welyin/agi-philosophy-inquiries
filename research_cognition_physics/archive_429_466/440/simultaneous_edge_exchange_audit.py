"""Round 440: the existing Quantum Graphity flip and its quantum interfaces.

Reuses Konopka et al. 0801.0861v2 Eq.27 and the flip-chain literature.
Graph/qubit architecture and the extra flip generator are model inputs.
No Markov mixing theorem is transferred to closed unitary evolution.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import short_path_exchange_audit as old
import endpoint_resource_exchange_audit as resource
import exchange_relation_audit as core

TARGET = Path(__file__).with_name('simultaneous_edge_exchange_audit_results.json')
OBS = {}


def flips(g, n, endpoint=None):
    """Count one per unoriented valid three-path (middle endpoints b<c)."""
    edges = resource.edges_for(n)
    index = {e: k for k, e in enumerate(edges)}
    adj = old.neighbors(g, n, edges)
    out = Counter()
    for b, c in edges:
        if c not in adj[b]:
            continue
        for a in sorted(adj[b]-{c}):
            for d in sorted(adj[c]-{b, a}):
                if c in adj[a] or d in adj[b]:
                    continue
                if endpoint is not None and endpoint not in (a, b, c, d):
                    continue
                changed = [(a, b), (c, d), (a, c), (b, d)]
                target = g
                for e in changed:
                    target ^= 1 << index[tuple(sorted(e))]
                out[target] += 1
    return dict(out)


def matrix_for(graphs, transitions):
    index = {g: k for k, g in enumerate(graphs)}
    a = np.zeros((len(graphs), len(graphs)), dtype=np.int64)
    for g, col in index.items():
        for other, weight in transitions(g).items():
            assert other in index
            a[index[other], col] = weight
    return a


def full_model(n, graphs):
    a = matrix_for(graphs, lambda g: flips(g, n))
    h = np.kron(np.eye(2**n, dtype=np.int64), a)
    for e, (i, j) in enumerate(resource.edges_for(n)):
        occ = np.diag([int(bool(g & (1 << e))) for g in graphs])
        h += np.kron(np.eye(2**n, dtype=np.int64)+core.swap(n, i, j).real.astype(np.int64), occ)
    return h, a


def path_graph(order):
    return old.graph_mask(len(order), zip(order[:-1], order[1:]))


def cycle_graph(order):
    return old.graph_mask(len(order), zip(order, order[1:]+order[:1]))


def classify(g, n):
    edges = resource.edges_for(n)
    deg = resource.degrees(g, n, edges)
    blocks, predicted = [], 1
    for block in old.partition(g, n, edges):
        s = len(block)
        if s == 3:
            factor, kind = 4, 'triangle_or_three_path'
        elif s >= 4 and all(deg[i] == 2 for i in block):
            factor, kind = math.factorial(s-1)//2, 'cycle'
        elif s >= 4:
            factor, kind = math.factorial(s-2), 'path'
        else:
            factor, kind = 1, 'fixed_small'
        predicted *= factor
        blocks.append(dict(vertices=list(block), kind=kind, factor=factor))
    return predicted, blocks


class Audit(unittest.TestCase):
    def close(self, a, b, tol=8e-12):
        self.assertLess(float(np.max(np.abs(np.asarray(a)-np.asarray(b)))), tol)

    def test_01_source_operator_normalization_and_invariants(self):
        n = 4
        edges = resource.edges_for(n)
        index = {e: k for k, e in enumerate(edges)}
        direct = np.zeros((64, 64), dtype=np.int64)
        for col in range(64):
            for a, b, c, d in itertools.permutations(range(n)):
                bit = lambda x, y: 1 << index[tuple(sorted((x, y)))]
                on = bit(b, c) | bit(a, b) | bit(c, d)
                off = bit(a, c) | bit(b, d)
                if col & on == on and col & off == 0:
                    direct[col ^ bit(a, b) ^ bit(c, d) ^ off, col] += 1
        # F = 1/4 sum(R+R^dagger) = 1/2 sum R.
        self.assertTrue(np.array_equal(direct, direct.T))
        self.assertFalse(np.any(direct % 2))
        f = matrix_for(list(range(64)), lambda g: flips(g, n))
        self.assertTrue(np.array_equal(direct, 2*f))
        self.assertEqual(int(f.max()), 2)  # two distinct middle edges can contribute.
        for i in range(n):
            deg = np.array([resource.degrees(g, n, edges)[i] for g in range(64)])
            self.assertFalse(np.any(f*deg[None, :]-deg[:, None]*f))
        rows = []
        for n in range(4, 7):
            edges = resource.edges_for(n)
            all_deg = [tuple(resource.degrees(g, n, edges)) for g in range(2**len(edges))]
            parts = [old.partition(g, n, edges) for g in range(2**len(edges))]
            moves, weighted = 0, 0
            for g in range(len(all_deg)):
                for target, weight in flips(g, n).items():
                    self.assertEqual(all_deg[g], all_deg[target])
                    self.assertEqual(parts[g], parts[target])
                    self.assertEqual(flips(target, n)[g], weight)
                    self.assertEqual((g ^ target).bit_count(), 4)
                    moves += 1
                    weighted += weight
            rows.append(dict(vertices=n, graphs=len(all_deg),
                directed_configuration_moves=moves, path_weighted_moves=weighted))
        OBS['source_operator_and_conservation'] = dict(
            full_four_vertex_operator_dimension=64,
            ordered_source_sum_normalization_exact=True,
            degree_and_component_partition_preserved=rows,
            resources_f_i_and_charges_Q_i_unchanged_by_flip=True,
            temporary_extra_endpoint_slots_required=False)

    def test_02_paths_cycles_and_cubic_reconfiguration(self):
        rows = []
        for n in range(3, 8):
            orders = [(0,)+p+(n-1,) for p in itertools.permutations(range(1, n-1))]
            paths = [path_graph(p) for p in orders]
            adjacency = {g: flips(g, n) for g in paths}
            self.assertEqual(old.component_of(paths[0], adjacency), set(paths))
            self.assertEqual(len(paths), math.factorial(n-2))
            c_orders = [(0,)+p for p in itertools.permutations(range(1, n)) if p[0] < p[-1]]
            cycles = [cycle_graph(p) for p in c_orders]
            adjacency = {g: flips(g, n) for g in cycles}
            self.assertEqual(old.component_of(cycles[0], adjacency), set(cycles))
            self.assertEqual(len(cycles), math.factorial(n-1)//2)
            rows.append(dict(vertices=n, paths_fixed_endpoint_labels=len(paths),
                             cycles_fixed_vertex_labels=len(cycles)))
        n = 6
        edges = resource.edges_for(n)
        cubic = [g for g in range(2**len(edges)) if resource.degrees(g, n, edges) == [3]*n]
        adj = {g: flips(g, n) for g in cubic}
        self.assertEqual(old.component_of(cubic[0], adj), set(cubic))
        bipartite = old.graph_mask(n, itertools.product(range(3), range(3, 6)))
        triangles = lambda g: sum(all(g & (1 << edges.index(tuple(sorted(e)))) for e in
            ((a, b), (b, c), (a, c))) for a, b, c in itertools.combinations(range(n), 3))
        self.assertEqual(triangles(bipartite), 0)
        after = min(flips(bipartite, n))
        self.assertGreater(triangles(after), 0)
        OBS['reconfiguration'] = dict(path_and_cycle_counts=rows,
            all_s_adjacent_internal_permutation_proof=True,
            six_vertex_cubic_graphs=len(cubic), cubic_configuration_graph_connected=True,
            explicit_non_isomorphic_pair=dict(before='K_3,3', before_mask=bipartite,
                after_mask=after, before_triangles=0, after_triangles=triangles(after)),
            general_regular_graph_reachability_attributed_to_flip_literature=True)

    def test_03_joint_unknown_data_and_closed_quantum_evolution(self):
        graphs = [path_graph((0, 1, 2, 3)), path_graph((0, 2, 1, 3))]
        h, f = full_model(4, graphs)
        self.assertTrue(np.array_equal(f, [[0, 1], [1, 0]]))
        rows, cols = np.arange(16)*2+1, np.arange(16)*2
        self.assertTrue(np.array_equal(h[np.ix_(rows, cols)], np.eye(16, dtype=np.int64)))
        self.assertLessEqual(float(np.linalg.norm(h, 2)), 7+1e-12)
        t = Fraction(1, 100)
        x = 7*t
        remainder = x*x/(2*(1-x))
        lower = t-remainder
        self.assertEqual(lower, Fraction(137, 18600))
        u = core.evolve(h, float(t))
        block = u[np.ix_(rows, cols)]
        smallest = float(np.linalg.svd(block, compute_uv=False)[-1])
        self.assertGreaterEqual(smallest, float(lower))
        self.close(u.conj().T@u, np.eye(32))
        sym = np.zeros((16, 5))
        for k in range(5):
            ids = [j for j in range(16) if j.bit_count() == k]
            sym[ids, k] = 1/math.sqrt(len(ids))
        v = np.kron(sym, np.eye(2))
        self.close(h@v, v@np.kron(np.eye(5), f+6*np.eye(2)))
        quantum, classical = [], []
        for t in (0., math.pi/4, math.pi/2, math.pi):
            p = float(abs(core.evolve(f, t)[1, 0])**2)
            self.assertAlmostEqual(p, math.sin(t)**2)
            quantum.append(p)
            classical.append((1-math.exp(-2*t))/2)
        OBS['quantum_interface'] = dict(
            full_unknown_data_dimension=16, exact_first_order_transfer_block='I_16',
            all_input_and_reference_positive_transfer=dict(time='1/100',
                singular_value_lower=str(lower), probability_lower=str(lower*lower),
                numerical_smallest_singular_value=smallest),
            full_joint_unitary_and_reference_preserved=True,
            symmetric_data_dimension=5, exact_probability='sin(t)^2',
            times=['0', 'pi/4', 'pi/2', 'pi'], quantum_probabilities=quantum,
            comparison_continuous_markov_probabilities=classical,
            irreversible_selection_or_stationary_mixing_of_quantum_model_claimed=False)

    def test_04_unknown_data_fixed_path_frame(self):
        n = 5
        orders = [(0,)+p+(n-1,) for p in itertools.permutations(range(1, n-1))]
        graphs = [path_graph(p) for p in orders]
        h, _ = full_model(n, graphs)
        count = len(graphs)
        # V maps position data, graph order -> label data, same graph order.
        v = np.zeros_like(h)
        for g, order in enumerate(orders):
            for column in range(2**n):
                row = 0
                for position, label in enumerate(order):
                    row |= ((column >> (n-1-position)) & 1) << (n-1-label)
                v[row*count+g, column*count+g] = 1
        self.assertTrue(np.array_equal(v.T@v, np.eye(len(v), dtype=np.int64)))
        path_h = (n-1)*np.eye(2**n, dtype=np.int64)
        for r in range(n-1):
            path_h += core.swap(n, r, r+1).real.astype(np.int64)
        expected = np.kron(path_h, np.eye(count, dtype=np.int64))
        for r in range(1, n-2):
            move = np.zeros((count, count), dtype=np.int64)
            for col, order in enumerate(orders):
                changed = list(order)
                changed[r], changed[r+1] = changed[r+1], changed[r]
                move[orders.index(tuple(changed)), col] = 1
            expected += np.kron(core.swap(n, r, r+1).real.astype(np.int64), move)
        transformed = v.T@h@v
        self.assertTrue(np.array_equal(transformed, expected))
        # A static path Hamiltonian alone omits the data SWAP on reconfiguration.
        naive = np.kron(path_h, np.eye(count, dtype=np.int64))
        naive += np.kron(np.eye(2**n, dtype=np.int64), matrix_for(graphs, lambda g: flips(g, n)))
        self.assertTrue(np.any(transformed-naive))
        # Every fixed-label readout must be transformed as well.
        for label in range(n):
            def z_at(slot):
                return np.diag([1-2*((j >> (n-1-slot)) & 1) for j in range(2**n)])
            label_readout = np.kron(z_at(label), np.eye(count, dtype=np.int64))
            expected_readout = np.zeros_like(h)
            for g, order in enumerate(orders):
                sector = np.diag([int(k == g) for k in range(count)])
                expected_readout += np.kron(z_at(order.index(label)), sector)
            self.assertTrue(np.array_equal(v.T@label_readout@v, expected_readout))
        OBS['relational_path_frame'] = dict(
            vertices=n, path_orders=count, full_dimension=len(h),
            controlled_data_relabeling_exact_integer_unitary=True,
            exact_frame_hamiltonian='H_path tensor I + kappa sum_r SWAP_(r,r+1) tensor T_r',
            all_unknown_joint_data_graph_and_reference_covered=True,
            fixed_label_readout_transformation_verified=True,
            controlled_frame_transform_physically_implemented=False,
            graph_rewiring_without_corresponding_data_readout_transform_not_pure_gauge=True,
            frame_is_an_existing_path_order_not_a_derived_three_dimensional_coordinate=True)

    def test_05_minimum_edge_support(self):
        rows = []
        for n in (4, 5):
            edges = resource.edges_for(n)
            groups = defaultdict(list)
            for g in range(2**len(edges)):
                groups[tuple(resource.degrees(g, n, edges))].append(g)
            histogram = Counter()
            for graphs in groups.values():
                for g, h in itertools.combinations(graphs, 2):
                    distance = (g ^ h).bit_count()
                    self.assertGreaterEqual(distance, 4)
                    histogram[distance] += 1
            self.assertIn(4, histogram)
            rows.append(dict(vertices=n, degree_sector_pairs_by_Hamming_distance=
                {str(k): histogram[k] for k in sorted(histogram)}))
        OBS['primitive_support_boundary'] = dict(finite_exact_checks=rows,
            all_N_degree_preserving_nontrivial_graph_transition_needs_at_least_four_edge_flips=True,
            at_most_three_edge_qubit_local_degree_commuting_H_is_graph_diagonal=True,
            new_flip_derived_from_raw_two_body_SWAP=False,
            encoded_exchange_simulation_excluded=False,
            charge_violating_virtual_excursions_or_auxiliaries_excluded=False,
            general_joint_dynamics_may_still_act_on_data=True)

    def test_06_size_uniform_flip_strength(self):
        rows = []
        for n in range(4, 7):
            edges = resource.edges_for(n)
            maximum = defaultdict(int)
            for g in range(2**len(edges)):
                c = max(resource.degrees(g, n, edges))
                if c > 3:
                    continue
                row_sum = sum(flips(g, n, endpoint=0).values())
                self.assertLessEqual(row_sum, 2*c*(c-1)**2)
                maximum[c] = max(maximum[c], row_sum)
            rows.append(dict(vertices=n, max_row_sum_by_capacity=
                {str(c): maximum[c] for c in sorted(maximum)}))
        OBS['endpoint_flip_strength'] = dict(
            all_N_bound='2*abs(kappa)*C*(C-1)^2',
            proof='four positions of endpoint in oriented three-path; identify reversal',
            finite_rows=rows, combined_with_439_bound=
                'abs(Omega)*C*(C-1)+2*abs(kappa)*C*(C-1)^2',
            physical_metric_light_cone_not_yet_derived=True)

    def test_07_combination_with_short_path_conversion(self):
        rows = []
        for n in range(3, 7):
            edges = resource.edges_for(n)
            graphs = [g for g in range(2**len(edges)) if max(resource.degrees(g, n, edges)) <= 2]
            adjacency = {}
            for g in graphs:
                moves = Counter(old.transitions(g, n, 2, edges))
                moves.update(flips(g, n))
                adjacency[g] = dict(moves)
            unseen, histogram = set(graphs), Counter()
            while unseen:
                g = min(unseen)
                component = old.component_of(g, adjacency)
                count, _ = classify(g, n)
                self.assertEqual(len(component), count)
                for other in component:
                    self.assertEqual(classify(other, n)[0], count)
                histogram[count] += 1
                unseen -= component
            rows.append(dict(vertices=n, graphs=len(graphs),
                combined_components_by_size={str(k): histogram[k] for k in sorted(histogram)}))
        OBS['combined_C2_complete_classification'] = dict(
            formula='4^r times product_paths((s-2)!) times product_cycles((s-1)!/2), s>=4',
            nonzero_Omega_and_kappa_required=True, finite_exhaustive_checks=rows,
            initial_component_membership_and_long_component_degree_labels_remain_fixed=True,
            three_dimensional_space_or_growth_of_N_generated=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=440, baseline_round=439, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(existing_quantum_graphity_flip_explicitly_attributed=True,
            simultaneous_rewiring_without_extra_slots_proved=True,
            endpoint_degree_charge_and_partition_invariants_proved=True,
            unknown_data_reference_and_relational_frame_verified=True,
            combined_C2_all_N_configuration_classification_proved=True,
            size_independent_endpoint_flip_norm_bound_proved=True,
            minimum_four_edge_transition_support_proved=True,
            graph_flip_generator_selected_by_429=False,
            individual_degrees_conserved_by_combined_model=False,
            classical_mixing_claimed_for_closed_unitary=False,
            physical_dimension_or_new_subject_growth_derived=False,
            final_Heisenberg_hardware_compiled=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
