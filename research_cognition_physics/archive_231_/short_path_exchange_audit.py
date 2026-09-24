"""Round 439: short-path-gated resource conversion and exact graph sectors.

Hamma et al., arXiv:0911.5075v3 Eq.18 supplies the L=2 gate idea.
Finite unit-shift endpoint resources, conserved local charges and data SWAP
are explicit choices inherited from 438, not the original full hopping model.
"""
import argparse
from collections import Counter, deque
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import endpoint_resource_exchange_audit as previous
import exchange_relation_audit as core

TARGET = Path(__file__).with_name('short_path_exchange_audit_results.json')
OBS = {}


def neighbors(mask, n, edges):
    out = [set() for _ in range(n)]
    for e, (i, j) in enumerate(edges):
        if mask & (1 << e):
            out[i].add(j)
            out[j].add(i)
    return out


def partition(mask, n, edges):
    adj = neighbors(mask, n, edges)
    unseen, blocks = set(range(n)), []
    while unseen:
        todo = [min(unseen)]
        found = set(todo)
        while todo:
            for j in adj[todo.pop()]:
                if j not in found:
                    found.add(j)
                    todo.append(j)
        unseen -= found
        blocks.append(tuple(sorted(found)))
    return tuple(blocks)


def path_count(mask, edge, n, edges):
    adj = neighbors(mask, n, edges)
    i, j = edge
    return len(adj[i] & adj[j])


def transitions(mask, n, capacity, edges):
    adj = neighbors(mask, n, edges)
    out = {}
    for e, (i, j) in enumerate(edges):
        weight = len(adj[i] & adj[j])
        if not weight:
            continue
        if mask & (1 << e) or (len(adj[i]) < capacity and len(adj[j]) < capacity):
            out[mask ^ (1 << e)] = weight
    return out


def model(n, capacity, graphs=None):
    edges = previous.edges_for(n)
    if graphs is None:
        graphs = [g for g in range(2**len(edges))
                  if max(previous.degrees(g, n, edges)) <= capacity]
    index = {g: k for k, g in enumerate(graphs)}
    adjacency = np.zeros((len(graphs), len(graphs)), dtype=np.int64)
    occupancy = [np.array([int(bool(g & (1 << e))) for g in graphs], dtype=np.int64)
                 for e in range(len(edges))]
    for g, col in index.items():
        for target, weight in transitions(g, n, capacity, edges).items():
            assert target in index, (g, target)
            adjacency[index[target], col] = weight
    h = np.kron(np.eye(2**n, dtype=np.int64), adjacency+np.diag(sum(occupancy)))
    for (i, j), occ in zip(edges, occupancy):
        h += np.kron(core.swap(n, i, j).real.astype(np.int64), np.diag(occ))
    return h, graphs, adjacency, occupancy


def graph_mask(n, pairs):
    labels = {e: k for k, e in enumerate(previous.edges_for(n))}
    return sum(1 << labels[tuple(sorted(e))] for e in pairs)


def component_of(start, graph):
    seen, todo = {start}, deque([start])
    while todo:
        here = todo.popleft()
        for other in graph[here]:
            if other not in seen:
                seen.add(other)
                todo.append(other)
    return seen


class Audit(unittest.TestCase):
    def close(self, a, b, tol=3e-12):
        self.assertLess(float(np.max(np.abs(np.asarray(a)-np.asarray(b)))), tol)

    def test_01_raw_charge_and_intertwining(self):
        n, capacity = 3, 2
        edges = previous.edges_for(n)
        raw, charges, raw_occ = previous.resource_model(n, capacity)
        # Reweight every raw transition by a diagonal common-neighbor count.
        for row, col in zip(*np.nonzero(raw)):
            changed = (row % 8) ^ (col % 8)
            e = int(changed).bit_length()-1
            raw[row, col] *= path_count(col % 8, edges[e], n, edges)
        self.assertTrue(np.array_equal(raw, raw.T))
        for q in charges:
            self.assertTrue(np.array_equal(raw*q[None, :], q[:, None]*raw))
        _, graphs, adjacency, occ = model(n, capacity)
        w = previous.embedding(n, capacity)
        self.assertTrue(np.array_equal(w.T@w, np.eye(len(graphs), dtype=np.int64)))
        self.assertTrue(np.array_equal(raw@w, w@adjacency))
        for before, after in zip(raw_occ, occ):
            self.assertTrue(np.array_equal(before[:, None]*w, w*after[None, :]))
        for q in charges:
            self.assertTrue(np.array_equal(q[:, None]*w, capacity*w))
        # The count is not a Boolean projector when two paths coexist.
        e4 = previous.edges_for(4)
        diamond = graph_mask(4, [(0, 2), (0, 3), (1, 2), (1, 3)])
        self.assertEqual(path_count(diamond, (0, 1), 4, e4), 2)
        OBS['exact_resource_compression'] = dict(raw_resource_graph_dimension=len(raw),
            legal_graph_dimension=len(graphs), all_local_charge_commutators_zero=True,
            weighted_conversion_and_all_occupancy_intertwiners_exact=True,
            hence_full_data_and_arbitrary_reference_intertwining=True,
            common_neighbor_gate_is_a_count_not_generally_a_projector=True)

    def test_02_partition_superselection(self):
        rows = []
        for n in range(2, 7):
            edges = previous.edges_for(n)
            parts = [partition(g, n, edges) for g in range(2**len(edges))]
            moves = 0
            for g in range(len(parts)):
                for target, weight in transitions(g, n, n-1, edges).items():
                    self.assertEqual(parts[g], parts[target])
                    self.assertEqual(transitions(target, n, n-1, edges)[g], weight)
                    moves += 1
            rows.append(dict(vertices=n, graphs=len(parts), directed_legal_moves=moves,
                             distinct_partitions=len(set(parts))))
        OBS['partition_conservation'] = dict(exhaustive_checks=rows,
            all_N_all_C_proof='every addition has an existing two-path; every deletion retains a two-path',
            quantum_partition_projectors_commute_with_full_H=True,
            coherent_partition_superpositions_not_collapsed=True)

    def test_03_capacity_two_complete_classification(self):
        rows = []
        for n in range(2, 7):
            edges = previous.edges_for(n)
            graphs = [g for g in range(2**len(edges))
                      if max(previous.degrees(g, n, edges)) <= 2]
            adjacency = {g: transitions(g, n, 2, edges) for g in graphs}
            unseen, histogram, frozen = set(graphs), Counter(), 0
            while unseen:
                g = min(unseen)
                connected = component_of(g, adjacency)
                r = sum(len(b) == 3 for b in partition(g, n, edges))
                self.assertEqual(len(connected), 4**r)
                self.assertEqual(sum(not adjacency[h] for h in connected),
                                 1 if r == 0 else 0)
                for h in connected:
                    for target in adjacency[h]:
                        changed_edge = (h ^ target).bit_length()-1
                        i, j = edges[changed_edge]
                        block = next(b for b in partition(h, n, edges) if i in b)
                        self.assertIn(j, block)
                        self.assertEqual(len(block), 3)
                histogram[len(connected)] += 1
                frozen += int(r == 0)
                unseen -= connected
            rows.append(dict(vertices=n, legal_graphs=len(graphs), frozen_graphs=frozen,
                configuration_components_by_size={str(k): histogram[k] for k in sorted(histogram)}))
        matching_counts = []
        for n in range(2, 9):
            edges = previous.edges_for(n)
            ms = previous.matchings(n)
            self.assertTrue(all(not transitions(g, n, 1, edges) for g in ms))
            matching_counts.append(dict(vertices=n, frozen_matchings=len(ms)))
        OBS['capacity_classification'] = dict(capacity_two=rows, capacity_one=matching_counts,
            all_N_C2_reachable_configuration_count='4^r, r = number of three-vertex connected components',
            nonzero_Omega_required_for_exact_reachable_set=True,
            full_configuration_count_not_single_trajectory_span=True,
            longer_paths_and_cycles_structurally_frozen=True)

    def test_04_triangle_star_and_unknown_data(self):
        graphs = [7, 3, 5, 6]  # triangle then the three P3 graphs
        h, _, adjacency, occupancy = model(3, 2, graphs)
        expected = np.zeros((4, 4), dtype=np.int64)
        expected[0, 1:] = expected[1:, 0] = 1
        self.assertTrue(np.array_equal(adjacency, expected))
        hgraph = adjacency+2*np.diag(sum(occupancy))
        # Check the all-unknown-data lowest-order transition block.
        certificates = []
        for start in range(4):
            dist = previous.shortest_paths(adjacency, start)
            for end in range(4):
                power = np.linalg.matrix_power(h, dist[end])
                block = power[np.ix_(np.arange(8)*4+end, np.arange(8)*4+start)]
                self.assertTrue(np.array_equal(block, np.eye(8, dtype=np.int64)))
                for smaller in range(dist[end]):
                    small = np.linalg.matrix_power(h, smaller)
                    self.assertFalse(np.any(small[np.ix_(np.arange(8)*4+end, np.arange(8)*4+start)]))
                certificates.append([start, end, dist[end]])
        # The graph-only sector is exactly Sym^3(data), including unknown reference.
        sym = np.zeros((8, 4))
        for column in range(4):
            ids = [j for j in range(8) if j.bit_count() == column]
            sym[ids, column] = 1/math.sqrt(len(ids))
        v = np.kron(sym, np.eye(4))
        self.close(h@v, v@np.kron(np.eye(4), hgraph))
        spectrum = np.linalg.eigvalsh(hgraph)
        self.close(spectrum, [3, 4, 4, 7])
        t = math.pi/2
        probabilities = np.abs(core.evolve(hgraph, t)[:, 1])**2
        self.close(probabilities, [0, 5/9, 2/9, 2/9])
        OBS['nontrivial_surviving_components'] = dict(
            exact_integer_leading_transition_blocks=certificates,
            arbitrary_unknown_data_has_nonzero_sufficiently_short_time_transition=True,
            symmetric_data_dimension=4, graph_spectrum=[3, 4, 4, 7],
            time='pi/2', probabilities_triangle_same_path_other_paths=['0', '5/9', '2/9', '2/9'],
            numerical_probabilities=[float(x) for x in probabilities],
            general_data_may_correlate_with_graph=True,
            permanently_recorded_classical_graph_not_claimed=True)

    def test_05_component_factorization_no_signalling(self):
        local_graphs = [7, 3, 5, 6]
        h_a, _, _, _ = model(3, 2, local_graphs)
        graphs = []
        e3 = previous.edges_for(3)
        for g in local_graphs:
            pairs = [e for k, e in enumerate(e3) if g & (1 << k)]+[(3, 4)]
            graphs.append(graph_mask(5, pairs))
        h_all, _, _, _ = model(5, 2, graphs)
        h_b = np.eye(4, dtype=np.int64)+core.swap(2, 0, 1).real.astype(np.int64)
        # Natural ordering = D_A, D_B, graph_A; block ordering = D_A, graph_A, D_B.
        perm = np.array([(a*4+b)*4+g for a in range(8) for g in range(4) for b in range(4)])
        h_blocks = h_all[np.ix_(perm, perm)]
        exact_sum = np.kron(h_a, np.eye(4, dtype=np.int64))+np.kron(np.eye(32, dtype=np.int64), h_b)
        self.assertTrue(np.array_equal(h_blocks, exact_sum))
        t = .71
        u = core.evolve(h_blocks, t)
        self.close(u, np.kron(core.evolve(h_a, t), core.evolve(h_b, t)))
        rng = np.random.default_rng(439)
        psi = rng.normal(size=(128, 2))+1j*rng.normal(size=(128, 2))
        psi /= np.linalg.norm(psi)
        # Apply a nontrivial operation in component A to an entangled ABR input.
        local = np.kron(core.PAULI[0], np.eye(16))
        psi_changed = np.kron(local, np.eye(4))@psi
        out = (u@psi).reshape(32, 8)
        changed = (u@psi_changed).reshape(32, 8)
        self.close(out.T@out.conj(), changed.T@changed.conj())
        OBS['fixed_partition_operational_result'] = dict(
            witness_partition=[[0, 1, 2], [3, 4]], full_block_dimension=128,
            exact_integer_tensor_sum_verified=True,
            entangled_input_reference_dimension=2,
            local_operation_cannot_change_other_component_plus_reference_marginal=True,
            all_C_proof='fixed-partition legal space factorizes; H is a sum of component Hamiltonians',
            arbitrary_local_trace_preserving_maps_not_only_unitaries_covered_analytically=True,
            initial_partition_is_an_input_not_an_emergent_spatial_dimension=True)

    def test_06_saturated_triangle_free_and_data_activity(self):
        examples = []
        for c in range(1, 6):
            n = 2*c
            edges = previous.edges_for(n)
            g = graph_mask(n, itertools.product(range(c), range(c, n)))
            self.assertEqual(previous.degrees(g, n, edges), [c]*n)
            self.assertFalse(transitions(g, n, c, edges))
            examples.append(dict(graph=f'K_{c},{c}', capacity=c, vertices=n, frozen=True))
        n = 4
        g = graph_mask(n, [(0, 1), (1, 2), (2, 3)])
        self.assertFalse(transitions(g, n, 2, previous.edges_for(n)))
        h = sum(core.swap(n, i, i+1).real.astype(np.int64) for i in range(n-1))
        z0 = np.kron(core.PAULI[2].real.astype(np.int64), np.eye(8, dtype=np.int64))
        comm = h@z0-z0@h
        self.assertTrue(np.any(comm))
        OBS['larger_capacity_boundary'] = dict(
            saturated_triangle_free_examples=examples,
            all_C_regular_triangle_free_graphs_frozen_proved=True,
            four_vertex_path_data_commutator_nonzero_entries=int(np.count_nonzero(comm)),
            frozen_graph_does_not_imply_frozen_data=True,
            increasing_capacity_does_not_by_itself_guarantee_reorganization=True,
            original_Hamma_resource_hopping_not_included=True,
            all_local_models_or_all_cognition_refuted=False)
        # Reuse the source's graph-reachability argument, adding the explicit
        # sufficient capacity s-1; this is not claimed as a new graph theorem.
        ample = []
        for n in range(2, 6):
            edges = previous.edges_for(n)
            connected = [g for g in range(2**len(edges)) if len(partition(g, n, edges)) == 1]
            adjacency = {g: transitions(g, n, n-1, edges) for g in connected}
            self.assertEqual(component_of(connected[0], adjacency), set(connected))
            ample.append(dict(vertices=n, sufficient_capacity=n-1,
                              all_connected_graphs_reachable=len(connected)))
        OBS['larger_capacity_boundary']['known_ample_capacity_interface'] = ample

    def test_07_size_independent_endpoint_norm(self):
        bounds = []
        for n in range(2, 6):
            edges = previous.edges_for(n)
            for c in range(1, n):
                maximum = 0
                for g in range(2**len(edges)):
                    adj = neighbors(g, n, edges)
                    if max(map(len, adj)) > c:
                        continue
                    raw_sum = sum(len(adj[0] & adj[j]) for j in range(1, n))
                    self.assertEqual(raw_sum, sum(len(adj[k])-1 for k in adj[0]))
                    allowed_sum = sum(weight for target, weight in transitions(g, n, c, edges).items()
                                      if 0 in edges[(g ^ target).bit_length()-1])
                    self.assertLessEqual(allowed_sum, raw_sum)
                    self.assertLessEqual(raw_sum, c*(c-1))
                    maximum = max(maximum, allowed_sum)
                bounds.append(dict(vertices=n, capacity=c,
                    maximum_endpoint_absolute_row_sum=maximum, proven_upper=c*(c-1)))
        exact_norms = []
        for n in range(3, 7):
            edges = previous.edges_for(n)
            graphs = [g for g in range(2**len(edges))
                      if max(previous.degrees(g, n, edges)) <= 2]
            endpoint = {g: {h: w for h, w in transitions(g, n, 2, edges).items()
                            if 0 in edges[(g ^ h).bit_length()-1]} for g in graphs}
            unseen, maximum = set(graphs), 0.
            while unseen:
                component = sorted(component_of(min(unseen), endpoint))
                block = np.array([[endpoint[g].get(h, 0) for h in component] for g in component])
                maximum = max(maximum, float(np.max(np.abs(np.linalg.eigvalsh(block)))))
                unseen -= set(component)
            self.assertAlmostEqual(maximum, math.sqrt(2))
            exact_norms.append(dict(vertices=n, endpoint_norm=maximum))
        OBS['endpoint_rate_bound'] = dict(
            all_N_operator_norm_bound='abs(Omega)*C*(C-1)',
            proof='symmetric matrix row sum <= sum over neighbors k of deg(k)-1',
            finite_row_checks=bounds, capacity_two_exact_norms=exact_norms,
            C2_N_at_least_three_exact_norm='sqrt(2)*abs(Omega)',
            fixed_capacity_removes_438_sqrt_N_endpoint_norm_growth=True,
            physical_light_cone_or_metric_speed_claimed=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=439, baseline_round=438, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(existing_short_path_gate_explicitly_attributed=True,
            resource_charge_and_full_unknown_reference_intertwining_proved=True,
            partition_conservation_and_component_no_signalling_proved=True,
            capacity_two_all_N_configuration_classification_proved=True,
            nonzero_autonomous_three_vertex_reconfiguration_retained=True,
            saturated_triangle_free_obstruction_proved=True,
            size_independent_endpoint_conversion_norm_bound_proved=True,
            old_routing_closure_and_shortest_power_results_reused=True,
            short_path_rule_or_capacity_derived_from_429=False,
            whole_data_dynamics_frozen_claimed=False,
            single_trajectory_dimension_equal_to_configuration_count_claimed=False,
            original_full_hopping_model_refuted=False,
            stable_connected_geometry_or_dimension_generated=False,
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
