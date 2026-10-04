"""Round 263: the same Fisher rule and resource counts need not fix dimension.

A perpetually selected degree-two tip produces a capped ladder, rather than
the uniform-generation Hanoi family. Selection fairness is NOT assumed.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from contact_route_closure_audit import adjacency
from vertex_split_cycle_audit import replace_vertex, TRIANGLE, statistics, edges
from fisher_hanoi_scaling_audit import bfs, level_metrics


@lru_cache(maxsize=None)
def focused_growth(t):
    if t < 1:
        raise ValueError('At least one triangle is required.')
    graph = replace_vertex(((),), 0, TRIANGLE, {})
    for event in range(1, t):
        tip = 2*event
        assert graph[tip] == (tip-2, tip-1)
        graph = replace_vertex(graph, tip, TRIANGLE, {tip-2: 0, tip-1: 1})
    return graph


def capped_ladder(t):
    links = [(2*j, 2*j+1) for j in range(t)]
    links += [(2*j+p, 2*(j+1)+p) for j in range(t-1) for p in (0, 1)]
    links += [(2*(t-1)+p, 2*t) for p in (0, 1)]
    return adjacency(2*t+1, links)


@lru_cache(maxsize=None)
def contrast(n):
    events = (3**n-1)//2
    graph = focused_growth(events)
    diameter = max(max(bfs(graph, v)) for v in range(len(graph)))
    uniform = level_metrics(n)
    return {'uniform_generation': n, 'triangle_events_each': events,
            'shared_counts': statistics(graph),
            'uniform_diameter': uniform['diameter_all_pairs'],
            'focused_diameter_all_pairs': diameter}


def report():
    ds = bfs(focused_growth(128), 0)
    return {'round': 263, 'equal_budget_comparisons': [contrast(n) for n in range(1, 7)],
            'focused_root_balls': [{'radius': r, 'vertices': sum(d <= r for d in ds)}
                                   for r in (1, 2, 4, 8, 16, 32, 64)],
            'focused_infinite_volume_exponent': 1,
            'scope': 'Pure Fisher triangles, the same single seed and identical vertex/edge/cycle counts yield different large-scale growth under complete generations versus persistent tip selection. The tip schedule is intentionally unfair to older vertices. This refutes schedule-independent dimension from local rule and loop density alone, not universality under a specified fair or quantum stochastic law.'}


class Audit(unittest.TestCase):
    def test_01_rewrite_equals_independent_ladder(self):
        for t in range(1, 65):
            self.assertEqual(focused_growth(t), capped_ladder(t))

    def test_02_exact_resource_counts(self):
        for t in range(1, 65):
            row = statistics(focused_growth(t))
            self.assertEqual((row['vertices'], row['edges'], row['cycle_rank']), (2*t+1, 3*t, t))
            self.assertEqual(row['components'], 1)

    def test_03_exact_focused_diameter(self):
        for n in range(1, 7):
            row = contrast(n)
            self.assertEqual(row['focused_diameter_all_pairs'], row['triangle_events_each'])

    def test_04_boundary_capacity_and_corners(self):
        for t in range(1, 65):
            graph = focused_growth(t)
            self.assertEqual([i for i, row in enumerate(graph) if len(row) == 2], [0, 1, 2*t])
            self.assertTrue(all(len(row) <= 3 for row in graph))

    def test_05_equal_counts_different_schedules(self):
        for n in range(1, 7):
            row = contrast(n)
            before = level_metrics(n)
            for key, value in row['shared_counts'].items():
                self.assertEqual(value, before[key])

    def test_06_exact_linear_root_balls(self):
        ds = bfs(focused_growth(128), 0)
        for radius in range(128):
            self.assertEqual(sum(d <= radius for d in ds), 2*radius+1)

    def test_07_stable_interior_exhaustion_is_half_ladder(self):
        for t in range(1, 65):
            old = {(u, v) for u, v in edges(focused_growth(t)) if v < 2*t}
            later = {(u, v) for u, v in edges(focused_growth(t+1)) if v < 2*t}
            self.assertEqual(old, later)
            expected = {(2*j, 2*j+1) for j in range(t)}
            expected |= {(2*j+p, 2*(j+1)+p) for j in range(t-1) for p in (0, 1)}
            self.assertEqual(old, expected)

    def test_08_same_counters_do_not_imply_isomorphism(self):
        for n in range(2, 7):
            row = contrast(n)
            self.assertGreater(row['focused_diameter_all_pairs'], row['uniform_diameter'])


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
    target = Path(__file__).with_name('fisher_update_schedule_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
