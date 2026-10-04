"""Round 272: demand, cut capacity and path work on supplied classical graphs.

An undirected edge has shared total service one per logical slot. Jobs are
independent, unicast, lossless and uncompressed. No spatial dimension is fitted.
"""
import argparse
from collections import deque
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np


def cycle(n):
    return tuple(tuple(sorted({(u-1) % n, (u+1) % n})) for u in range(n))


def grid(length, dimension):
    words = list(itertools.product(range(length), repeat=dimension))
    index = {w: i for i, w in enumerate(words)}
    graph = []
    for w in words:
        neighbors = []
        for axis in range(dimension):
            for delta in (-1, 1):
                other = list(w)
                other[axis] += delta
                if 0 <= other[axis] < length:
                    neighbors.append(index[tuple(other)])
        graph.append(tuple(sorted(neighbors)))
    return tuple(graph), words


def distances(graph):
    rows = []
    for root in range(len(graph)):
        row = [-1]*len(graph)
        row[root] = 0
        todo = deque([root])
        while todo:
            u = todo.popleft()
            for v in graph[u]:
                if row[v] < 0:
                    row[v] = row[u]+1
                    todo.append(v)
        assert min(row) >= 0
        rows.append(row)
    return rows


def edge_list(graph):
    return [(u, v) for u, row in enumerate(graph) for v in row if u < v]


def cut(graph, subset):
    return sum((u in subset) != (v in subset) for u, v in edge_list(graph))


def cut_rate_bound(n, s, boundary, mu=F(1)):
    return mu*boundary*F(n-1, 2*s*(n-s))


def cycle_flow(n, rate):
    """All ordered unicast demands, shortest paths, equal split of antipodes."""
    loads = {e: F(0) for e in edge_list(cycle(n))}
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            clockwise = (v-u) % n
            anticlockwise = n-clockwise
            directions = ([1] if clockwise < anticlockwise else
                          [-1] if anticlockwise < clockwise else [1, -1])
            demand = rate/F(n-1)/len(directions)
            for step in directions:
                a = u
                while a != v:
                    b = (a+step) % n
                    loads[tuple(sorted((a, b)))] += demand
                    a = b
    return loads


def benchmark(kind, length, dimension=1):
    if kind == 'cycle':
        n, edges, boundary = length, length, 2
        mean_distance = F(n*n, 4*(n-1))
    else:
        n = length**dimension
        edges = dimension*(length-1)*length**(dimension-1)
        boundary = length**(dimension-1)
        mean_distance = F(dimension*(length*length-1)*n, 3*length*(n-1))
    return {'family': kind, 'linear_size': length, 'dimension_input': dimension,
            'vertices': n, 'edges': edges, 'half_cut': boundary,
            'mean_distance_excluding_self': float(mean_distance),
            'cut_per_node_rate_upper_bound': float(cut_rate_bound(n, n//2, boundary)),
            'path_work_per_node_rate_upper_bound': float(F(edges, n)/mean_distance)}


def report():
    return {
        'round': 272,
        'edge_total_service_per_slot': 1,
        'demand': 'Every node originates rate lambda; each other destination has rate lambda/(N-1).',
        'cycle_rows': [benchmark('cycle', n) for n in (16, 64, 256)],
        'grid_rows': [benchmark('grid', length, d) for d in (1, 2, 3)
                      for length in (4, 8, 16)],
        'cycle_fractional_capacity_certificates': [
            {'vertices': n, 'rate': float(cut_rate_bound(n, n//2, 2)),
             'min_edge_load': float(min(cycle_flow(n, cut_rate_bound(n, n//2, 2)).values())),
             'max_edge_load': float(max(cycle_flow(n, cut_rate_bound(n, n//2, 2)).values()))}
            for n in (4, 8, 16, 32)],
        'three_dimensional_grid_lambda_one_cut_deficit_per_slot': [
            {'length': length, 'deficit': float(F(length**6, 2*(length**3-1))-length**2)}
            for length in (4, 8, 16)],
        'local_load_counterexample': 'Independent adjacent-edge work at rate 1/2 per edge is served immediately at capacity 1; every supplied graph family supports it.',
        'scope': 'Necessary cut and path-work inequalities for supplied classical unicast tasks and fixed service graphs; explicit fractional-flow certificates, not quantum communication, a topology generator, or a derivation of spatial dimension.'
    }


class Audit(unittest.TestCase):
    def test_01_ordered_demand_crossings_for_all_small_subsets(self):
        for n in (4, 6, 8):
            for mask in range(1, (1 << n)-1):
                subset = {u for u in range(n) if mask >> u & 1}
                enumerated = sum(F(1, n-1) for u in range(n) for v in range(n)
                                 if (u in subset) != (v in subset))
                self.assertEqual(enumerated, F(2*len(subset)*(n-len(subset)), n-1))

    def test_02_cycle_mean_distance_from_breadth_first_search(self):
        for n in (4, 6, 8, 16, 32):
            total = sum(map(sum, distances(cycle(n))))
            self.assertEqual(F(total, n*(n-1)), F(n*n, 4*(n-1)))

    def test_03_grid_mean_distance_and_edge_counts_independent_enumeration(self):
        for d in (1, 2, 3):
            for length in (2, 3, 4):
                graph, words = grid(length, d)
                n = len(graph)
                self.assertEqual(len(edge_list(graph)), d*(length-1)*length**(d-1))
                self.assertEqual(F(sum(map(sum, distances(graph))), n*(n-1)),
                                 F(d*(length*length-1)*n, 3*length*(n-1)))

    def test_04_actual_grid_bisections_have_claimed_capacity(self):
        for d in (1, 2, 3):
            for length in (2, 4, 6):
                graph, words = grid(length, d)
                subset = {i for i, w in enumerate(words) if w[0] < length//2}
                self.assertEqual(len(subset)*2, len(graph))
                self.assertEqual(cut(graph, subset), length**(d-1))

    def test_05_cycle_cut_bound_has_fractional_flow_certificate(self):
        for n in (4, 6, 8, 16, 32):
            rate = cut_rate_bound(n, n//2, 2)
            loads = cycle_flow(n, rate)
            self.assertEqual(set(loads.values()), {F(1)})
            hop_work = rate*n*F(n*n, 4*(n-1))
            self.assertEqual(sum(loads.values()), hop_work)

    def test_06_optimistic_cut_queue_accumulates_unavoidable_deficit(self):
        # This oracle instantly routes every crossing job to the cut. A real
        # network cannot complete more cut crossings than the same capacity.
        for length in (4, 8, 16):
            n, capacity = length**3, length**2
            arrival = F(n*n, 2*(n-1))
            q = F(0)
            for t in range(1, 65):
                q = max(F(0), q+arrival-capacity)
                self.assertEqual(q, t*(arrival-capacity))
            self.assertGreater(q, 0)

    def test_07_local_workload_respects_every_cut_without_volume_scaling(self):
        graph = cycle(8)
        for mask in range(1, 255):
            subset = {u for u in range(8) if mask >> u & 1}
            crossing_arrival = sum(F(1, 2) for u, v in edge_list(graph)
                                   if (u in subset) != (v in subset))
            self.assertEqual(crossing_arrival, F(cut(graph, subset), 2))
            self.assertLess(crossing_arrival, cut(graph, subset))

    def test_08_rate_bound_scales_with_length_for_every_fixed_grid_dimension(self):
        for d in (1, 2, 3, 4):
            for length in (4, 8, 16):
                n = length**d
                bound = cut_rate_bound(n, n//2, length**(d-1))
                self.assertEqual(bound, F(2*(n-1), length*n))
                self.assertLess(bound, F(2, length))


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
    target = Path(__file__).with_name('traffic_capacity_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
