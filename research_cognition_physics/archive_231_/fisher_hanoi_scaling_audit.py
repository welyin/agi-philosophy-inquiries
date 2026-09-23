"""Round 261: identify uniform Fisher generations with three-peg Hanoi graphs.

The seed, complete generations, unit edge lengths and rooted embedding are inputs.
This is a deterministic geometric control, not a sampled quantum growth law.
"""
import argparse
from collections import deque
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
from fisher_loop_birth_audit import triangle_all
from vertex_split_cycle_audit import edges, statistics


@lru_cache(maxsize=None)
def hanoi(n):
    """Independent legal-move construction; disk positions smallest first."""
    words = tuple(itertools.product(range(3), repeat=n))
    index = {w: i for i, w in enumerate(words)}
    rows = []
    for w in words:
        adjacent = set()
        for a, b in ((0, 1), (0, 2), (1, 2)):
            for k, peg in enumerate(w):
                if peg in (a, b):
                    new = w[:k]+(a+b-peg,)+w[k+1:]
                    adjacent.add(index[new])
                    break
        rows.append(tuple(sorted(adjacent)))
    return tuple(rows), words


@lru_cache(maxsize=None)
def fisher_with_labels(n):
    if n == 0:
        return ((),), ((),)
    old, labels = fisher_with_labels(n-1)
    output = triangle_all(old)
    new_labels = []
    for v, row in enumerate(old):
        child_pegs = []
        for u in row:
            changed = [(a, b) for a, b in zip(labels[v], labels[u]) if a != b]
            assert len(changed) == 1
            a, b = changed[0]
            child_pegs.append(3-a-b)
        assert len(set(child_pegs)) == len(row)
        child_pegs.extend(sorted(set(range(3))-set(child_pegs)))
        new_labels.extend((peg,)+labels[v] for peg in child_pegs)
    return output, tuple(new_labels)


def bfs(graph, start):
    distances = [-1]*len(graph)
    distances[start] = 0
    queue = deque([start])
    while queue:
        v = queue.popleft()
        for u in graph[v]:
            if distances[u] < 0:
                distances[u] = distances[v]+1
                queue.append(u)
    return distances


@lru_cache(maxsize=None)
def level_metrics(n):
    graph, words = hanoi(n)
    raw, labels = fisher_with_labels(n)
    index = {w: i for i, w in enumerate(words)}
    permutation = tuple(index[w] for w in labels)
    mapped = {tuple(sorted((permutation[u], permutation[v]))) for u, v in edges(raw)}
    assert mapped == set(edges(graph))
    diameter = max(max(bfs(graph, v)) for v in range(len(graph)))
    return {'generation': n, **statistics(graph), 'diameter_all_pairs': diameter,
            'explicit_edge_isomorphism': True, 'mapped_edges': len(mapped)}


@lru_cache(maxsize=None)
def root_metrics():
    # H_8 contains the entire infinite rooted ball of radius 255.
    graph, words = hanoi(8)
    distances = bfs(graph, 0)
    dyadic = [{'n': n, 'radius': 2**n-1,
               'ball_vertices': sum(d <= 2**n-1 for d in distances)}
              for n in range(9)]
    exits = []
    for n in range(7):
        larger, larger_words = hanoi(n+1)
        ds = bfs(larger, 0)
        first = min(ds[i] for i, w in enumerate(larger_words) if w[-1] != 0)
        exits.append({'n': n, 'first_exit_distance': first})
    return {'finite_level_for_root_balls': 8, 'dyadic_balls': dyadic,
            'first_exits_from_nested_copies': exits}


def report():
    return {'round': 261, 'levels': [level_metrics(n) for n in range(1, 7)],
            'rooted_growth': root_metrics(),
            'volume_exponent_log3_over_log2': math.log(3)/math.log(2),
            'scope': 'A single seed, uniform pure triangle generations, unit graph distances and a corner-rooted nested limit identify exactly with simple three-peg Hanoi graphs. This does not establish the dimension of mixed quantum growth, a continuum manifold, or physical space.'}


class Audit(unittest.TestCase):
    def test_01_explicit_isomorphism_not_just_counts(self):
        for n in range(1, 7):
            self.assertTrue(level_metrics(n)['explicit_edge_isomorphism'])

    def test_02_exact_counts_and_cycle_rank(self):
        for n in range(1, 7):
            row = level_metrics(n)
            self.assertEqual((row['vertices'], row['edges'], row['cycle_rank']),
                             (3**n, 3*(3**n-1)//2, (3**n-1)//2))

    def test_03_only_three_corners_have_degree_two(self):
        for n in range(1, 8):
            graph, words = hanoi(n)
            self.assertEqual({words[i] for i, row in enumerate(graph) if len(row) == 2},
                             {(peg,)*n for peg in range(3)})
            self.assertTrue(all(len(row) in (2, 3) for row in graph))

    def test_04_diameter_is_exact(self):
        for n in range(1, 7):
            self.assertEqual(level_metrics(n)['diameter_all_pairs'], 2**n-1)

    def test_05_largest_disk_zero_gives_induced_nested_copy(self):
        for n in range(7):
            old, words = hanoi(n)
            new, new_words = hanoi(n+1)
            index = {w: i for i, w in enumerate(new_words)}
            embedding = {i: index[w+(0,)] for i, w in enumerate(words)}
            image = set(embedding.values())
            mapped = {tuple(sorted((embedding[u], embedding[v]))) for u, v in edges(old)}
            induced = {(u, v) for u, v in edges(new) if u in image and v in image}
            self.assertEqual(mapped, induced)

    def test_06_first_exit_prevents_hidden_ball_shortcuts(self):
        for row in root_metrics()['first_exits_from_nested_copies']:
            self.assertEqual(row['first_exit_distance'], 2**row['n'])

    def test_07_exact_root_balls(self):
        for row in root_metrics()['dyadic_balls']:
            self.assertEqual(row['ball_vertices'], 3**row['n'])

    def test_08_all_integer_radius_bounds(self):
        ds = bfs(hanoi(8)[0], 0)
        exponent = math.log(3)/math.log(2)
        for radius in range(1, 256):
            size = sum(d <= radius for d in ds)
            self.assertGreaterEqual(size+1e-9, radius**exponent/3)
            self.assertLessEqual(size, 3*radius**exponent+1e-9)

    def test_09_peg_names_do_not_change_graph(self):
        graph, words = hanoi(4)
        index = {w: i for i, w in enumerate(words)}
        for pegs in itertools.permutations(range(3)):
            perm = [index[tuple(pegs[x] for x in w)] for w in words]
            renamed = {tuple(sorted((perm[u], perm[v]))) for u, v in edges(graph)}
            self.assertEqual(renamed, set(edges(graph)))


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
    target = Path(__file__).with_name('fisher_hanoi_scaling_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
