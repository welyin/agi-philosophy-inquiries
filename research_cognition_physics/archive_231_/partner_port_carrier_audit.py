"""Round 444: reciprocal address carriers for the existing graph-flip model.

Quantum address registers are an established idea (Arrighi et al., 2023).
This audit supplies a bounded-subject-support realization and exact intertwiner
for the particular 440 Hamiltonian.  It does not generate physical locality.
"""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import simultaneous_edge_exchange_audit as graph_model
import branching_tree_distance_audit as trees_model
import exchange_relation_audit as data_model

TARGET = Path(__file__).with_name('partner_port_carrier_audit_results.json')
OBS = {}


def decode(port, n, capacity):
    """Return a simple graph, or None for an invalid raw address configuration."""
    if len(port) != n*capacity:
        return None
    edges = set()
    for u, v in enumerate(port):
        if v == -1:
            continue
        if not 0 <= v < len(port) or u//capacity == v//capacity or port[v] != u:
            return None
        if u < v:
            e = tuple(sorted((u//capacity, v//capacity)))
            if e in edges:
                return None
            edges.add(e)
    return frozenset(edges)


@lru_cache(None)
def lifts(graph, n, capacity):
    neighbors = trees_model.adjacency(graph, n)
    if max(map(len, neighbors), default=0) > capacity:
        return ()
    choices = [list(itertools.permutations(range(capacity), len(a))) for a in neighbors]
    result = []
    for assignments in itertools.product(*choices):
        local = [dict(zip(sorted(neighbors[i]), assignments[i])) for i in range(n)]
        port = [-1]*(n*capacity)
        for i, j in graph:
            u, v = i*capacity+local[i][j], j*capacity+local[j][i]
            port[u], port[v] = v, u
        result.append(tuple(port))
    return tuple(result)


def reciprocal_pairs(port, i, j, capacity):
    return [(u, v) for u in range(i*capacity, (i+1)*capacity)
            for v in [port[u]]
            if j*capacity <= v < (j+1)*capacity and port[v] == u]


def local_rewrite(port, a, b, c, d, capacity, strict=True):
    """One directed partial-bijection kernel; reads only four subject blocks."""
    assert len({a, b, c, d}) == 4
    ab = reciprocal_pairs(port, a, b, capacity)
    cd = reciprocal_pairs(port, c, d, capacity)
    if (not reciprocal_pairs(port, b, c, capacity)
            or reciprocal_pairs(port, a, c, capacity)
            or reciprocal_pairs(port, b, d, capacity)):
        return ()
    if strict and (len(ab) != 1 or len(cd) != 1):
        return ()
    result = []
    for u, v in ab:
        for w, z in cd:
            target = list(port)
            target[u], target[w], target[v], target[z] = w, u, z, v
            result.append(tuple(target))
    return tuple(result)


def port_flips(port, n, capacity):
    """Fast valid-code evaluation of the same four-subject local kernels."""
    graph = decode(port, n, capacity)
    assert graph is not None
    neighbors = trees_model.adjacency(graph, n)
    result = Counter()
    for b, c in sorted(graph):
        for a in sorted(neighbors[b]-{c}):
            for d in sorted(neighbors[c]-{b, a}):
                for target in local_rewrite(port, a, b, c, d, capacity):
                    result[target] += 1
    return result


def graph_mask(graph, n):
    return graph_model.old.graph_mask(n, graph)


def graph_from_mask(mask, n):
    return frozenset(e for k, e in enumerate(graph_model.resource.edges_for(n))
                     if mask & (1 << k))


def all_partial_matchings(items):
    """Every reciprocal matching, allowing same-subject and parallel edges."""
    if not items:
        yield ()
        return
    u, *rest = items
    for pairs in all_partial_matchings(rest):
        yield ((u, -1),)+pairs
    for j, v in enumerate(rest):
        for pairs in all_partial_matchings(rest[:j]+rest[j+1:]):
            yield ((u, v),)+pairs


def raw_matching(pairs, size):
    result = [-1]*size
    for u, v in pairs:
        if v >= 0:
            result[u], result[v] = v, u
    return tuple(result)


def rename_port(port, n, capacity, vertex_permutation, local_permutations):
    transform = {i*capacity+s: vertex_permutation[i]*capacity+local_permutations[i][s]
                 for i in range(n) for s in range(capacity)}
    result = [-1]*len(port)
    for u, v in enumerate(port):
        result[transform[u]] = -1 if v == -1 else transform[v]
    return tuple(result)


class Audit(unittest.TestCase):
    def close(self, a, b, tol=2e-11):
        self.assertLess(float(np.max(np.abs(np.asarray(a)-np.asarray(b)))), tol)

    def test_01_complete_bounded_graph_lifts_and_resource_budget(self):
        n, capacity = 4, 2
        graph_count = total = 0
        all_codes = set()
        for mask in range(64):
            graph = graph_from_mask(mask, n)
            degree = list(map(len, trees_model.adjacency(graph, n)))
            if max(degree) > capacity:
                continue
            family = lifts(graph, n, capacity)
            expected = math.prod(math.factorial(capacity)//math.factorial(capacity-d)
                                 for d in degree)
            self.assertEqual(len(family), expected)
            self.assertEqual(len(set(family)), expected)
            for code in family:
                self.assertEqual(decode(code, n, capacity), graph)
                self.assertNotIn(code, all_codes)
                all_codes.add(code)
            graph_count += 1
            total += len(family)
        raw = [raw_matching(pairs, n*capacity)
               for pairs in all_partial_matchings(list(range(n*capacity)))]
        valid = {p for p in raw if decode(p, n, capacity) is not None}
        self.assertEqual(valid, all_codes)
        rows = [dict(subjects=n, capacity=3, port_qudits=3*n,
                     dimension_per_port=3*n+1,
                     binary_address_budget=3*n*(3*n).bit_length(),
                     old_potential_edge_qubits=n*(n-1)//2)
                for n in (6, 32, 128, 1024)]
        OBS['carrier_and_complete_fibers'] = dict(
            subjects=4, capacity=2, graphs=graph_count, valid_port_basis_states=total,
            reciprocal_matchings_including_invalid=len(raw),
            uniform_lift_isometry_exact=True, binary_storage_examples=rows,
            fixed_local_dimension_claimed=False)

    def test_02_local_extension_validity_and_naive_adjoint_leakage(self):
        n, capacity = 4, 2
        rewrites = 0
        for pairs in all_partial_matchings(list(range(n*capacity))):
            p = raw_matching(pairs, n*capacity)
            valid = decode(p, n, capacity) is not None
            for a, b, c, d in itertools.permutations(range(n)):
                for q in local_rewrite(p, a, b, c, d, capacity):
                    self.assertEqual(valid, decode(q, n, capacity) is not None)
                    self.assertIn(p, local_rewrite(q, a, c, b, d, capacity))
                    changed = {i//capacity for i in range(len(p)) if p[i] != q[i]}
                    self.assertEqual(changed, {a, b, c, d})
                    rewrites += 1
        # Without the "old edges unique" gate, an invalid double edge can be
        # repaired.  The adjoint then leaks OUT of the valid simple-graph code.
        bad = (3, 4, -1, 0, 1, 6, 5, 9, -1, 7, -1, -1)
        self.assertIsNone(decode(bad, 4, 3))
        naive = local_rewrite(bad, 0, 1, 2, 3, 3, strict=False)
        self.assertEqual(len(naive), 2)
        self.assertTrue(all(decode(q, 4, 3) is not None for q in naive))
        self.assertEqual(local_rewrite(bad, 0, 1, 2, 3, 3), ())
        # Unreciprocated external pointers cannot be repaired by the rewrite.
        path = frozenset({(0, 1), (1, 2), (2, 3)})
        p = lifts(path, 5, 2)[0]
        checked_external = 0
        for pointer in range(8):
            bad_external = list(p)
            bad_external[8] = pointer
            bad_external = tuple(bad_external)
            self.assertIsNone(decode(bad_external, 5, 2))
            for q in local_rewrite(bad_external, 0, 1, 2, 3, 2):
                self.assertIsNone(decode(q, 5, 2))
                self.assertEqual(q[8:], bad_external[8:])
                checked_external += 1
        OBS['four_subject_raw_tensor_extension'] = dict(
            reciprocal_matching_directed_rewrites_checked=rewrites,
            global_validity_equivalent_before_and_after=True,
            external_unreciprocated_pointer_examples=checked_external,
            naive_kernel_invalid_to_valid_example=list(bad),
            naive_kernel_valid_targets=len(naive),
            strict_local_gate_blocks_adjoint_leakage=True,
            global_validity_projector_used_to_claim_locality=False)

    def test_03_exact_flip_intertwiner_including_branching_and_multiplicity(self):
        trees, masks, _, _, _ = trees_model.six_vertex_sector()
        families = [
            ('all_four_subject_C2_graphs', 4, 2,
             [graph_from_mask(m, 4) for m in range(64)
              if max(graph_model.resource.degrees(m, 4, graph_model.resource.edges_for(4))) <= 2]),
            ('six_branching_trees', 6, 3, trees)]
        report = []
        for name, n, capacity, graphs in families:
            graph_set = set(graphs)
            directed = weighted = 0
            largest_weight = 0
            for graph in graphs:
                expected = graph_model.flips(graph_mask(graph, n), n)
                received = Counter()
                for p in lifts(graph, n, capacity):
                    transitions = port_flips(p, n, capacity)
                    quotient = Counter()
                    for q, weight in transitions.items():
                        target_graph = decode(q, n, capacity)
                        self.assertIn(target_graph, graph_set)
                        quotient[graph_mask(target_graph, n)] += weight
                        received[q] += weight
                        self.assertEqual(port_flips(q, n, capacity)[p], weight)
                        self.assertEqual(sum(x >= 0 for x in p), sum(x >= 0 for x in q))
                        largest_weight = max(largest_weight, weight)
                    self.assertEqual(dict(quotient), expected)
                    directed += len(transitions)
                    weighted += sum(transitions.values())
                # Integer unnormalised fiber columns Z satisfy F_port Z=ZF.
                wanted = Counter()
                for target_mask, weight in expected.items():
                    target_graph = graph_from_mask(target_mask, n)
                    self.assertEqual(len(lifts(graph, n, capacity)),
                                     len(lifts(target_graph, n, capacity)))
                    for q in lifts(target_graph, n, capacity):
                        wanted[q] = weight
                self.assertEqual(received, wanted)
            report.append(dict(name=name, subjects=n, capacity=capacity,
                graph_basis_states=len(graphs),
                port_basis_states=sum(len(lifts(g, n, capacity)) for g in graphs),
                directed_port_transitions=directed, path_weighted_transitions=weighted,
                maximum_parallel_path_weight=largest_weight))
        self.assertEqual(report[0]['maximum_parallel_path_weight'], 2)
        OBS['exact_graph_intertwining'] = report

    def test_04_full_data_dynamics_and_actual_readout_intertwining(self):
        n, capacity = 4, 2
        graphs = [frozenset({(0, 1), (1, 2), (2, 3)}),
                  frozenset({(0, 2), (1, 2), (1, 3)})]
        codes = [p for g in graphs for p in lifts(g, n, capacity)]
        index = {p: i for i, p in enumerate(codes)}
        f_port = np.zeros((len(codes), len(codes)), dtype=np.int64)
        for p, column in index.items():
            for q, weight in port_flips(p, n, capacity).items():
                f_port[index[q], column] = weight
        z = np.zeros((len(codes), 2), dtype=np.int64)
        for i, p in enumerate(codes):
            z[i, graphs.index(decode(p, n, capacity))] = 1
        self.assertTrue(np.array_equal(z.T @ z, 16*np.eye(2, dtype=int)))
        h_graph, f_graph = graph_model.full_model(n, [graph_mask(g, n) for g in graphs])
        h_port = np.kron(np.eye(16, dtype=np.int64), f_port)
        for i, j in itertools.combinations(range(n), 2):
            occupation = np.diag([int(bool(reciprocal_pairs(p, i, j, capacity))) for p in codes])
            local_data = np.eye(16, dtype=np.int64)+data_model.swap(n, i, j).real.astype(np.int64)
            h_port += np.kron(local_data, occupation)
        total_z = np.kron(np.eye(16, dtype=np.int64), z)
        self.assertTrue(np.array_equal(h_port @ total_z, total_z @ h_graph))
        for b in range(n):
            readout = np.diag([(v >> (n-1-b)) & 1 for v in range(16)])
            old_readout = np.kron(readout, np.eye(2, dtype=int))
            port_readout = np.kron(readout, np.eye(len(codes), dtype=int))
            self.assertTrue(np.array_equal(port_readout @ total_z, total_z @ old_readout))
        w = total_z/4
        u_port, u_graph = data_model.evolve(h_port, .3), data_model.evolve(h_graph, .3)
        error = float(np.linalg.norm(u_port @ w-w @ u_graph, 2))
        self.assertLess(error, 2e-11)
        OBS['complete_data_and_readout'] = dict(
            old_graph_data_dimension=32, port_data_dimension=len(h_port),
            exact_integer_H_intertwiner=True, every_subject_readout_intertwiner=True,
            finite_time_operator_error=error,
            all_encoded_joint_data_graph_reference_inputs_covered=True,
            arbitrary_nonuniform_port_fiber_inputs_equated_to_old_graph=False)

    def test_05_three_subject_support_obstruction_and_four_subject_witness(self):
        trees, *_ = trees_model.six_vertex_sector()
        n, capacity = 6, 3
        tagged_codes = [(p, g) for g in trees for p in lifts(g, n, capacity)]
        comparisons = 0
        for support in itertools.combinations(range(n), 3):
            outside = [i for i in range(n*capacity) if i//capacity not in support]
            classes = {}
            for p, graph in tagged_codes:
                key = tuple(p[i] for i in outside)
                if key in classes:
                    self.assertEqual(classes[key], graph)
                else:
                    classes[key] = graph
                comparisons += 1
        p = tagged_codes[0][0]
        q = next(iter(port_flips(p, n, capacity)))
        changed_blocks = sorted({i//capacity for i in range(len(p)) if p[i] != q[i]})
        self.assertEqual(len(changed_blocks), 4)
        self.assertNotEqual(decode(p, n, capacity), decode(q, n, capacity))
        OBS['support_boundary'] = dict(
            fixed_degree_vector=[1, 3, 3, 1, 1, 1],
            three_subject_subsets=math.comb(n, 3),
            outside_record_checks=comparisons,
            distinct_graphs_never_share_outside_records=True,
            changing_graph_witness_subjects=changed_blocks,
            arbitrary_sum_of_at_most_three_subject_terms_has_no_cross_graph_code_block=True,
            approximate_virtual_processes_or_extra_auxiliaries_excluded=False)

    def test_06_address_relabeling_and_local_port_covariance(self):
        n, capacity = 4, 2
        graph = frozenset({(0, 1), (1, 2), (2, 3), (0, 3)})
        tested = 0
        for p in lifts(graph, n, capacity):
            source_moves = port_flips(p, n, capacity)
            for vertex_perm in itertools.permutations(range(n)):
                # Simultaneously rename the physical owner and every stored target.
                local_perm = [((1, 0) if (i+vertex_perm[i]) % 2 else (0, 1))
                              for i in range(n)]
                renamed = rename_port(p, n, capacity, vertex_perm, local_perm)
                expected = Counter()
                for q, weight in source_moves.items():
                    expected[rename_port(q, n, capacity, vertex_perm, local_perm)] += weight
                self.assertEqual(port_flips(renamed, n, capacity), expected)
                tested += 1
        OBS['address_symmetry'] = dict(
            joint_vertex_and_port_renamings_checked=tested,
            numerical_address_order_not_used_as_distance=True,
            address_renaming_not_identified_with_generated_space=True)


def run():
    OBS.clear()
    buffer = io.StringIO()
    result = unittest.TextTestRunner(stream=buffer, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(buffer.getvalue())
    return dict(round=444, baseline_round=443, date='2026-09-24',
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(
            existing_quantum_address_register_idea_reused=True,
            explicit_subject_tensor_factorization_supplied=True,
            no_independent_qubit_per_potential_edge_required_by_this_realization=True,
            exact_440_graph_data_H_and_subject_readout_intertwining=True,
            unknown_encoded_joint_data_graph_reference_preserved=True,
            four_subject_extension_without_global_code_projector_proved=True,
            exact_three_subject_support_obstruction_in_fixed_degree_code=True,
            universal_fixed_dimension_subject_hardware_constructed=False,
            four_subject_interaction_derived_from_429_two_body_rule=False,
            address_access_and_interaction_architecture_generated=False,
            general_nonuniform_port_states_equivalent_to_old_graph=False,
            all_old_relation_operations_locally_preserved=False,
            three_dimensional_space_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))
