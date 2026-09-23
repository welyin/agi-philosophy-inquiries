"""Round 270: restricting matchings to newest cells freezes older separators."""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import random
import unittest
import numpy as np
from paired_interface_growth_audit import refine, check_ports
from ancestral_quotient_audit import copy_state, graph_data
from vertex_split_cycle_audit import statistics, edges
from fisher_hanoi_scaling_audit import bfs
from coarse_diffusion_audit import laplacian

ROOTS = {(0, 0), (1, 0)}


def starting_state():
    state = {(0,): {1: ((1,), 1)}, (1,): {1: ((0,), 1)}}
    refine(state, (0,), (1,))
    refine(state, (0, 0), (1, 0))
    refine(state, (0, 1), (1, 1))
    return state


def cells(state):
    unseen = set(state)
    result = []
    while unseen:
        u = min(unseen)
        a, b = state[u][0][0], state[u][1][0]
        cell = {u, a, b, state[a][1][0]}
        assert len(cell) == 4 and cell <= unseen
        assert all(state[v][p][0] in cell for v in cell for p in (0, 1))
        result.append(tuple(sorted(cell)))
        unseen -= cell
    return result


def step(state, choices, flips=None):
    groups = cells(state)
    assert len(choices) == len(groups)
    matching = []
    for group, port in zip(groups, choices):
        assert port in (0, 1)
        matching.extend((u, state[u][port][0]) for u in group if u < state[u][port][0])
    assert len({w for edge in matching for w in edge}) == len(state)
    if flips is None:
        flips = [(0, 0)]*len(matching)
    for (u, v), assignment in zip(matching, flips):
        refine(state, u, v, assignment)
    return state


def partition_stats(state):
    selected = {u for u in state if u[:2] in ROOTS}
    boundary = sum((u in selected) != (v in selected)
                   for u, row in state.items() for v, _ in row.values() if u < v)
    return len(selected), boundary


@lru_cache(maxsize=None)
def history(mode, depth, value=270):
    rng = random.Random(value)
    state = starting_state()
    for _ in range(2, depth):
        count = len(cells(state))
        choices = ([int(mode)]*count if mode in ('0', '1')
                   else [rng.randrange(2) for _ in range(count)])
        assignments = ([(rng.randrange(2), rng.randrange(2)) for _ in range(len(state)//2)]
                       if mode == 'mixed_twisted' else None)
        step(state, choices, assignments)
    return state


@lru_cache(maxsize=None)
def comparison():
    rows = []
    for mode in ('1', '0', 'mixed', 'mixed_twisted'):
        for depth in (4, 6, 8):
            state = history(mode, depth)
            graph, _ = graph_data(state)
            size, boundary = partition_stats(state)
            rows.append({'mode': mode, 'depth': depth, **statistics(graph),
                         'diameter': max(max(bfs(graph, v)) for v in range(len(graph))),
                         'protected_half_size': size, 'protected_half_boundary': boundary,
                         'boundary_over_half_volume_power_two_thirds': boundary/size**(2/3),
                         'laplacian_gap': float(np.linalg.eigvalsh(laplacian(graph))[1]),
                         'proved_gap_upper_bound': 8/len(graph)})
    return rows


@lru_cache(maxsize=None)
def exhaustive_two_steps():
    count = 0
    for first in itertools.product((0, 1), repeat=2):
        state = step(starting_state(), first)
        for second in itertools.product((0, 1), repeat=4):
            final = step(copy_state(state), second)
            assert partition_stats(final) == (16, 2)
            check_ports(final)
            count += 1
    return count


def report():
    return {'round': 270, 'exact_two_step_cell_choice_histories': exhaustive_two_steps(),
            'equal_budget_comparison': comparison(),
            'scope': 'After a fixed eight-vertex initial cut, all matchings restricted to ports 0 and 1 stay inside newest cells. The two-edge balanced ancestral separator persists for arbitrary cell choices and endpoint flips. Different diameters do not remove this exact bottleneck; no dimension is inferred from log N/log diameter.'}


class Audit(unittest.TestCase):
    def test_01_cells_are_disjoint_and_cover_graph(self):
        for mode in ('0', '1', 'mixed', 'mixed_twisted'):
            state = history(mode, 6)
            flat = [u for group in cells(state) for u in group]
            self.assertEqual(len(flat), len(set(flat)))
            self.assertEqual(set(flat), set(state))
            self.assertTrue(all(len({u[:2] in ROOTS for u in group}) == 1 for group in cells(state)))

    def test_02_all_small_matching_histories_keep_separator(self):
        self.assertEqual(exhaustive_two_steps(), 64)

    def test_03_randomized_endpoint_assignments_keep_separator(self):
        for value in range(32):
            for depth in (3, 4, 6):
                state = history('mixed_twisted', depth, value)
                self.assertEqual(partition_stats(state), (2**depth, 2))
                check_ports(state)

    def test_04_equal_budget_does_not_determine_diameter(self):
        rows = [r for r in comparison() if r['depth'] == 8]
        self.assertTrue(all((r['vertices'], r['edges']) == (512, 766) for r in rows))
        self.assertGreater(len({r['diameter'] for r in rows}), 1)
        self.assertTrue(all(r['protected_half_boundary'] == 2 for r in rows))

    def test_05_pointwise_rayleigh_bound(self):
        for mode in ('0', '1', 'mixed_twisted'):
            state = history(mode, 5)
            graph, words = graph_data(state)
            f = np.array([1 if w[:2] in ROOTS else -1 for w in words], dtype=float)
            matrix = laplacian(graph)
            self.assertAlmostEqual(float(f@matrix@f/(f@f)), 8/len(graph))
            self.assertLessEqual(float(np.linalg.eigvalsh(matrix)[1]), 8/len(graph)+1e-12)

    def test_06_allowing_an_old_edge_changes_the_protected_cut(self):
        state = starting_state()
        u, v = next((u, v) for u, row in state.items() for v, _ in row.values()
                    if u < v and (u[:2] in ROOTS) != (v[:2] in ROOTS))
        refine(state, u, v)
        self.assertEqual(partition_stats(state), (5, 3))

    def test_07_resource_ledger_agrees_with_operation_counts(self):
        for row in comparison():
            operations = (row['vertices']-2)//2
            self.assertEqual(row['edges'], 1+3*operations)
            self.assertEqual(row['cycle_rank'], operations)
            self.assertLessEqual(row['max_degree'], 3)

    def test_08_known_ladder_control_is_recovered(self):
        for row in comparison():
            if row['mode'] == '1':
                self.assertEqual(row['diameter'], 2**row['depth'])
                expected = 2-2*np.cos(np.pi/2**row['depth'])
                self.assertAlmostEqual(row['laplacian_gap'], expected, places=11)


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
    target = Path(__file__).with_name('recent_cell_matching_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
