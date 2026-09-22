"""Round 259: boundary-preserving vertex replacement and the cycle budget.

Finite simple undirected channel graphs, supplied degree cap and rewrites.
Numerical labels enumerate vertices, not physical positions or time.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from contact_route_closure_audit import adjacency

EDGE = ((1,), (0,))
TRIANGLE = ((1, 2), (0, 2), (0, 1))


def edges(graph):
    return tuple((u, v) for u, row in enumerate(graph) for v in row if u < v)


def components(graph):
    unseen = set(range(len(graph)))
    count = 0
    while unseen:
        count += 1
        pending = [unseen.pop()]
        while pending:
            for w in graph[pending.pop()]:
                if w in unseen:
                    unseen.remove(w); pending.append(w)
    return count


def statistics(graph):
    n, m, c = len(graph), len(edges(graph)), components(graph)
    return {'vertices': n, 'edges': m, 'components': c,
            'cycle_rank': m-n+c, 'max_degree': max(map(len, graph), default=0)}


def incidence_cycle_rank(graph):
    es = edges(graph)
    incidence = np.zeros((len(graph), len(es)))
    for k, (u, v) in enumerate(es):
        incidence[u, k] = 1; incidence[v, k] = -1
    rank = int(np.linalg.matrix_rank(incidence)) if es else 0
    return len(es)-rank


def replace_vertex(graph, v, gadget, assignment, degree_cap=3):
    """Reassign each old incident edge exactly once; add a connected gadget.

    Child zero temporarily retains index v; the other children get fresh indices.
    These are numerical names only. assignment maps OLD neighbor -> child index.
    """
    if components(gadget) != 1:
        raise ValueError('The replacement gadget must be connected.')
    if set(assignment) != set(graph[v]):
        raise ValueError('Each old boundary edge must be assigned exactly once.')
    if any(i < 0 or i >= len(gadget) for i in assignment.values()):
        raise ValueError('Invalid output port.')
    children = (v,)+tuple(range(len(graph), len(graph)+len(gadget)-1))
    new_edges = [edge for edge in edges(graph) if v not in edge]
    new_edges += [(children[assignment[w]], w) for w in graph[v]]
    new_edges += [(children[a], children[b]) for a, b in edges(gadget)]
    result = adjacency(len(graph)+len(gadget)-1, new_edges)
    if max(map(len, result), default=0) > degree_cap:
        raise ValueError('Output port capacity exceeded.')
    return result


def binary_assignments(graph, v):
    for values in itertools.product((0, 1), repeat=len(graph[v])):
        if max(values.count(0), values.count(1)) <= 2:
            yield dict(zip(graph[v], values))


def subcubic_graphs(n=5):
    possible = tuple(itertools.combinations(range(n), 2))
    for mask in range(1 << len(possible)):
        graph = adjacency(n, [edge for k, edge in enumerate(possible) if mask >> k & 1])
        if max(map(len, graph), default=0) <= 3:
            yield graph


@lru_cache(maxsize=None)
def exhaustive():
    graph_count = split_count = tree_count = 0
    for graph in subcubic_graphs():
        before = statistics(graph)
        graph_count += 1
        tree_count += int(before['components'] == 1 and before['cycle_rank'] == 0)
        for v in range(len(graph)):
            for assignment in binary_assignments(graph, v):
                output = replace_vertex(graph, v, EDGE, assignment)
                after = statistics(output)
                assert after['vertices'] == before['vertices']+1
                assert after['edges'] == before['edges']+1
                assert after['components'] == before['components']
                assert after['cycle_rank'] == before['cycle_rank']
                assert incidence_cycle_rank(output) == before['cycle_rank']
                split_count += 1
    return {'five_vertex_subcubic_graphs': graph_count, 'labeled_trees_among_them': tree_count,
            'admissible_binary_rewrites': split_count, 'incidence_nullity_checks': split_count}


def growing_tree(steps):
    graph = ((),)
    for k in range(steps):
        v = k % len(graph)
        choices = list(binary_assignments(graph, v))
        graph = replace_vertex(graph, v, EDGE, choices[k % len(choices)])
    return graph


def report():
    return {'round': 259, 'degree_cap': 3, 'exhaustive': exhaustive(),
            'single_seed_after_40_binary_splits': statistics(growing_tree(40)),
            'general_replacement_cycle_increment': 'internal_edges - output_vertices + 1',
            'scope': 'Connected boundary-preserving gadgets replace vertices in finite simple undirected channel graphs. Binary edge gadgets preserve cycle rank; a single seed remains a tree. This is not a restriction on all quantum graphs or a spacetime dimension theorem.'}


class Audit(unittest.TestCase):
    def test_01_exhaustive_binary_rewrites(self):
        self.assertGreater(exhaustive()['admissible_binary_rewrites'], 10000)
        self.assertEqual(exhaustive()['labeled_trees_among_them'], 120)

    def test_02_single_seed_remains_tree(self):
        s = statistics(growing_tree(40))
        self.assertEqual((s['vertices'], s['edges'], s['components'], s['cycle_rank']), (41, 40, 1, 0))
        self.assertLessEqual(s['max_degree'], 3)

    def test_03_existing_cycle_not_removed(self):
        graph = TRIANGLE
        for _ in range(8):
            graph = replace_vertex(graph, 0, EDGE, next(binary_assignments(graph, 0)))
        self.assertEqual(statistics(graph)['cycle_rank'], 1)

    def test_04_disconnected_components_preserved(self):
        graph = adjacency(5, [(0, 1), (2, 3)])
        out = replace_vertex(graph, 4, EDGE, {})
        self.assertEqual(components(graph), 3)
        self.assertEqual(components(out), 3)

    def test_05_general_gadget_cycle_budget(self):
        graph = adjacency(4, [(0, 1), (1, 2), (2, 0), (0, 3)])
        square = ((1, 3), (0, 2), (1, 3), (0, 2))
        path3 = ((1,), (0, 2), (1,))
        for gadget in (EDGE, path3, TRIANGLE, square):
            for outputs in itertools.product(range(len(gadget)), repeat=3):
                assignment = dict(zip(graph[0], outputs))
                result = replace_vertex(graph, 0, gadget, assignment, degree_cap=6)
                self.assertEqual(statistics(result)['cycle_rank']-statistics(graph)['cycle_rank'],
                                 len(edges(gadget))-len(gadget)+1)

    def test_06_capacity_and_boundary_rejections(self):
        star = adjacency(4, [(0, 1), (0, 2), (0, 3)])
        with self.assertRaises(ValueError):
            replace_vertex(star, 0, EDGE, {1: 0, 2: 0, 3: 0})
        with self.assertRaises(ValueError):
            replace_vertex(star, 0, EDGE, {1: 0, 2: 1})
        with self.assertRaises(ValueError):
            replace_vertex(star, 0, ((), ()), {1: 0, 2: 1, 3: 1})

    def test_07_relabeling_transports_port_assignment(self):
        graph = adjacency(4, [(0, 1), (0, 2), (0, 3), (2, 3)])
        permutation = (2, 0, 3, 1)
        assignment = {1: 0, 2: 1, 3: 1}
        output = replace_vertex(graph, 0, EDGE, assignment)
        renamed = adjacency(4, [(permutation[u], permutation[v]) for u, v in edges(graph)])
        other = replace_vertex(renamed, permutation[0], EDGE,
                               {permutation[w]: child for w, child in assignment.items()})
        extended = permutation+(4,)
        self.assertEqual(other, adjacency(5, [(extended[u], extended[v]) for u, v in edges(output)]))

    def test_08_contraction_recovers_original_graph(self):
        for graph in subcubic_graphs(4):
            for v in range(4):
                for a in binary_assignments(graph, v):
                    out = replace_vertex(graph, v, EDGE, a)
                    contracted = [(v if x == 4 else x, v if y == 4 else y)
                                  for x, y in edges(out) if set((x, y)) != {v, 4}]
                    self.assertEqual(adjacency(4, contracted), graph)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['runtime'] = {'python': platform.python_version(), 'numpy': np.__version__}
    result['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    target = Path(__file__).with_name('vertex_split_cycle_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
