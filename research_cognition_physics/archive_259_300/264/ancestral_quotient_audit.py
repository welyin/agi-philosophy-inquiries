"""Round 264: ancestral quotients and distances for asynchronous Fisher growth.

Persistent local port names and ancestry are explicit inputs. Rewrites are
atomic abstract graph operations, not a protocol for messages in flight.
"""
import argparse
from collections import Counter
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import random
import unittest
import numpy as np
from fisher_hanoi_scaling_audit import bfs


def seed():
    return {(): {}}


def copy_state(state):
    return {v: dict(row) for v, row in state.items()}


def rewrite(state, vertex):
    old = state.pop(vertex)
    children = [vertex+(i,) for i in range(3)]
    for child in children:
        state[child] = {}
    for i, j in itertools.combinations(range(3), 2):
        state[children[i]][j] = (children[j], i)
        state[children[j]][i] = (children[i], j)
    for port, (neighbor, other_port) in old.items():
        state[children[port]][port] = (neighbor, other_port)
        state[neighbor][other_port] = (children[port], port)


def complete(state, depth, prefix=()):
    """Finish all required ancestors; deeper leaves are not processed again."""
    while True:
        eligible = [w for w in state if len(w) < depth and w[:len(prefix)] == prefix]
        if not eligible:
            return state
        for w in sorted(eligible):
            rewrite(state, w)


def uniform(n):
    return complete(seed(), n)


def fingerprint(state):
    return tuple((w, tuple(sorted(row.items()))) for w, row in sorted(state.items()))


def graph_data(state):
    words = tuple(sorted(state))
    index = {w: i for i, w in enumerate(words)}
    graph = tuple(tuple(sorted(index[v] for v, _ in state[w].values())) for w in words)
    return graph, words


def canonical_edges(n):
    """Independent closed form for the complete Sierpinski/Hanoi graph."""
    result = set()
    for k in range(n):
        for prefix in itertools.product(range(3), repeat=k):
            for a, b in itertools.combinations(range(3), 2):
                x = prefix+(a,)+(b,)*(n-k-1)
                y = prefix+(b,)+(a,)*(n-k-1)
                result.add(tuple(sorted((x, y))))
    return result


def quotient_edges(state, depth):
    if min(map(len, state)) < depth:
        raise ValueError('The proposed ancestor cut is not complete.')
    counter = Counter()
    for u, row in state.items():
        for v, _ in row.values():
            if u < v and u[:depth] != v[:depth]:
                counter[tuple(sorted((u[:depth], v[:depth])))] += 1
    return counter


def partition(state, depth):
    graph, words = graph_data(state)
    if min(map(len, words)) < depth:
        raise ValueError('Incomplete cut.')
    ancestors = tuple(sorted({w[:depth] for w in words}))
    lookup = {w: i for i, w in enumerate(ancestors)}
    labels = np.array([lookup[w[:depth]] for w in words], dtype=int)
    blocks = [np.flatnonzero(labels == j) for j in range(len(ancestors))]
    return graph, words, ancestors, labels, blocks


def block_diameter(graph, block):
    index = {int(v): j for j, v in enumerate(block)}
    induced = tuple(tuple(index[v] for v in graph[int(u)] if v in index) for u in block)
    ds = [bfs(induced, i) for i in range(len(induced))]
    assert all(min(row) >= 0 for row in ds)
    return max(map(max, ds), default=0)


def metric_check(state, depth, lag):
    graph, _, _, labels, blocks = partition(state, depth)
    coarse = graph_data(uniform(depth))[0]
    distances = np.array([bfs(graph, v) for v in range(len(graph))])
    base = np.array([bfs(coarse, v) for v in range(len(coarse))])
    projected = base[np.ix_(labels, labels)]
    assert np.all(projected <= distances)
    assert np.all(distances <= 2**lag*projected+2**lag-1)
    assert max(map(len, blocks)) <= 3**lag
    assert max(block_diameter(graph, block) for block in blocks) <= 2**lag-1
    return len(graph)**2


def random_snapshot(depth, lag, value, steps=20):
    rng = random.Random(value)
    state = uniform(depth)
    for _ in range(steps):
        eligible = sorted(w for w in state if len(w) < depth+lag)
        if not eligible:
            break
        rewrite(state, rng.choice(eligible))
    return state


@lru_cache(maxsize=None)
def exhaustive_one_lag():
    names = tuple(sorted(uniform(2)))
    pairs = 0
    for mask in range(1 << len(names)):
        state = uniform(2)
        for i, w in enumerate(names):
            if mask >> i & 1:
                rewrite(state, w)
        counter = quotient_edges(state, 2)
        assert set(counter) == canonical_edges(2) and set(counter.values()) == {1}
        pairs += metric_check(state, 2, 1)
    return {'snapshots': 512, 'ordered_distance_pairs': pairs}


@lru_cache(maxsize=None)
def burst_rows():
    rows = []
    for extra in range(1, 7):
        state = uniform(1)
        complete(state, 1+extra, (0,))
        graph, words = graph_data(state)
        start, end = words.index((0,)*(extra+1)), words.index((1,))
        distance = bfs(graph, start)[end]
        assert set(quotient_edges(state, 1)) == canonical_edges(1)
        rows.append({'extra_depth': extra, 'vertices': len(graph),
                     'same_coarse_K3': True, 'coarse_distance': 1,
                     'fine_distance': distance, 'predicted_fine_distance': 2**extra})
    return rows


def fair_demo():
    state = uniform(1)
    depth = 1
    rows = []
    for extra in (1, 2):
        old = set(state)
        complete(state, depth+extra, (0,)*depth)
        burst_size = len(state)
        next_depth = depth+extra+1
        complete(state, next_depth)
        assert old.isdisjoint(state)
        assert all(len(w) == next_depth for w in state)
        rows.append({'start_depth': depth, 'extra_depth': extra,
                     'burst_vertices': burst_size, 'catchup_depth': next_depth,
                     'catchup_vertices': len(state), 'all_old_leaves_processed': True})
        depth = next_depth
    return rows


def report():
    return {'round': 264, 'exhaustive_lag_one': exhaustive_one_lag(),
            'unbounded_burst_counterexamples': burst_rows(), 'fair_finite_prefix': fair_demo(),
            'metric_bound': 'd_Q <= d_G <= 2^b d_Q + (2^b-1)',
            'fiber_size_bound': '3^b', 'fiber_intrinsic_diameter_bound': '2^b-1',
            'scope': 'Pure Fisher growth with persistent local port identities has exact ancestral quotients at complete cuts. Bounded depth overshoot gives uniform metric and volume bounds. Eventual fairness alone allows unbounded finite bursts; neither physical spacetime nor quantum dynamics is derived.'}


class Audit(unittest.TestCase):
    def test_01_persistent_ports_are_reciprocal_and_subcubic(self):
        for value in range(8):
            state = random_snapshot(2, 2, value)
            for u, row in state.items():
                self.assertLessEqual(len(row), 3)
                for p, (v, q) in row.items():
                    self.assertEqual(state[v][q], (u, p))
                    self.assertNotEqual(u, v)

    def test_02_complete_graph_matches_independent_edge_formula(self):
        for n in range(1, 6):
            state = uniform(n)
            actual = {tuple(sorted((u, v))) for u, row in state.items() for v, _ in row.values()}
            self.assertEqual(actual, canonical_edges(n))
            self.assertEqual(len(state), 3**n)

    def test_03_nonuniform_descendants_have_exact_quotient(self):
        for n, value in itertools.product((1, 2, 3), range(8)):
            state = random_snapshot(n, 2, value)
            counter = quotient_edges(state, n)
            self.assertEqual(set(counter), canonical_edges(n))
            self.assertEqual(set(counter.values()), {1})

    def test_04_same_event_set_different_orders(self):
        target = fingerprint(uniform(3))
        for value in range(64):
            rng = random.Random(value)
            state = seed()
            while any(len(w) < 3 for w in state):
                rewrite(state, rng.choice(sorted(w for w in state if len(w) < 3)))
            self.assertEqual(fingerprint(state), target)

    def test_05_adjacent_rewrites_commute_with_transported_ports(self):
        base = uniform(2)
        for u, v in canonical_edges(2):
            first, second = copy_state(base), copy_state(base)
            rewrite(first, u); rewrite(first, v)
            rewrite(second, v); rewrite(second, u)
            self.assertEqual(fingerprint(first), fingerprint(second))

    def test_06_nested_quotients_and_incomplete_cut_rejection(self):
        state = random_snapshot(3, 2, 91)
        fine_quotient = quotient_edges(state, 3)
        for n in (1, 2):
            collapsed = {tuple(sorted((u[:n], v[:n]))) for u, v in fine_quotient if u[:n] != v[:n]}
            self.assertEqual(collapsed, set(quotient_edges(state, n)))
        with self.assertRaises(ValueError):
            quotient_edges(uniform(1), 2)

    def test_07_exhaustive_one_lag_metric_bounds(self):
        self.assertEqual(exhaustive_one_lag()['snapshots'], 512)

    def test_08_mixed_two_lag_metric_bounds(self):
        for n, value in itertools.product((1, 2, 3), range(8)):
            metric_check(random_snapshot(n, 2, value), n, 2)

    def test_09_same_quotient_has_arbitrarily_large_stretch(self):
        for row in burst_rows():
            self.assertEqual(row['fine_distance'], row['predicted_fine_distance'])
            self.assertEqual(row['coarse_distance'], 1)

    def test_10_catchup_stages_are_fair_on_every_old_leaf(self):
        self.assertTrue(all(row['all_old_leaves_processed'] for row in fair_demo()))
        self.assertEqual(fair_demo()[-1]['catchup_vertices'], 729)


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
    target = Path(__file__).with_name('ancestral_quotient_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
