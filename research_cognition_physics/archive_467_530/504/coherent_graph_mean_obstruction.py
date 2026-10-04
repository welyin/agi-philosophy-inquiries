"""Round 504: uniform-in-size port error for replacing graph operators by means.

Original 440 edge-controlled exchange, no 501--503 identity tag is added.
The obstruction concerns one coherent mean-Hamiltonian replacement only.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import branching_tree_distance_audit as tree_tools
import role_covariant_graph_obstruction as previous

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'coherent_graph_mean_obstruction_results.json'
OBS = {}


def short(value):
    return float(f'{float(value):.12g}')


@lru_cache(None)
def model(i):
    assert i in (2, 3)
    trees, _, _ = previous.ensemble(i)
    n, g = 2*i+2, len(trees)
    lookup = {tree: j for j, tree in enumerate(trees)}
    f = np.zeros((g, g), dtype=np.int64)
    lap = np.zeros((n*g, n*g), dtype=np.int64)
    adj = []
    for j, tree in enumerate(trees):
        a = np.zeros((n, n), dtype=np.int64)
        for u, v in tree:
            a[u, v] = a[v, u] = 1
        adj.append(a)
        lg = np.diag(a.sum(axis=1))-a
        lap[j::g, j::g] = lg
        for other, weight in tree_tools.tree_flips(tree, n).items():
            f[lookup[other], j] = weight
    hf = np.kron(np.eye(n, dtype=np.int64), f)
    h = hf-lap
    w = np.kron(np.eye(n, dtype=np.int64), np.ones((g, 1), dtype=np.int64))
    bar = -(w.T @ lap @ w)/g
    return trees, f, lap, hf, h, w, bar, adj


def evolve(h, columns, time, order=18):
    result = np.array(columns, dtype=complex)
    term = result.copy()
    for k in range(1, order+1):
        term = (-1j*float(time)/k)*(h @ term)
        result += term
    return result


def certificate():
    t = Q(1, 1024)
    factor = (1-48*t)**2-(Q(3, 4)+18*t)**2
    assert factor == Q(83695, 262144) and factor > Q(1, 4)
    return t, factor, t*t/4


class Audit(unittest.TestCase):
    def test_01_original_full_tensor_and_excitation_interface(self):
        i = 2
        trees, f, _, _, h, _, _, _ = model(i)
        n, g = 6, len(trees)
        raw = np.kron(np.eye(2**n, dtype=np.int64), f)
        for j, tree in enumerate(trees):
            for a, b in tree:
                swap = np.zeros((2**n, 2**n), dtype=np.int64)
                ma, mb = 1 << (n-1-a), 1 << (n-1-b)
                for x in range(2**n):
                    y = x if bool(x & ma) == bool(x & mb) else x ^ ma ^ mb
                    swap[y, x] = 1
                raw[j::g, j::g] += swap
        indices = [g*(1 << (n-1-a))+j for a in range(n) for j in range(g)]
        embed = np.eye(2**n*g, dtype=np.int64)[:, indices]
        self.assertTrue(np.array_equal(raw @ embed, embed @ (h+(n-1)*np.eye(n*g, dtype=np.int64))))
        OBS['full_tensor'] = dict(dimension=384, excitation_dimension=36,
            exact_intertwining=True, graph_and_data_kept=True, identity_tag_added=False)

    def test_02_full_tree_families_and_local_commutator_count(self):
        records = []
        for i in (2, 3, 4):
            trees, _, _ = previous.ensemble(i)
            n, g = 2*i+2, len(trees)
            tree_set = set(trees)
            sums = np.zeros((n, n), dtype=np.int64)
            maximum = 0
            for tree in trees:
                adj = tree_tools.adjacency(tree, n)
                for a, b in tree:
                    sums[a, b] += 1
                    sums[b, a] += 1
                flips = tree_tools.tree_flips(tree, n)
                self.assertEqual(sum(flips.values()), 4*(i-1))
                rows = [0]*n
                for other, weight in flips.items():
                    self.assertIn(other, tree_set)
                    other_adj = tree_tools.adjacency(other, n)
                    for v in range(n):
                        rows[v] += weight*len(adj[v] ^ other_adj[v])
                maximum = max(maximum, max(rows))
            self.assertLessEqual(maximum, 48)
            self.assertTrue(np.all(i*sums[i:, :i] == g))
            records.append(dict(I=i, trees=g, F_row_sum=4*(i-1),
                commutator_absolute_row_sum_max=maximum, leaf_parent_marginals_uniform=True))
        OBS['finite_combinatorial_crosschecks'] = records

    def test_03_source_block_moment_and_graph_leakage(self):
        records = []
        for i in (2, 3):
            trees, f, lap, hf, h, w, bar, _ = model(i)
            n, g = 2*i+2, len(trees)
            r = 4*(i-1)
            k = h-r*np.eye(n*g, dtype=np.int64)
            self.assertTrue(np.array_equal(k @ w, -lap @ w))
            self.assertTrue(np.array_equal(k @ k @ w, -(hf @ lap-lap @ hf) @ w+lap @ lap @ w))
            self.assertLessEqual(np.linalg.norm(k @ k @ w/np.sqrt(g), 2), 84)
            # Mean compression and its source column are checked with integers.
            mean_i = np.rint(i*bar).astype(np.int64)
            self.assertTrue(np.array_equal(i*(w.T @ h @ w-r*g*np.eye(n, dtype=np.int64)), g*mean_i))
            source = i
            vec = k @ w[:, source]
            norm2 = Q(int(vec @ vec), g)
            projected2 = sum(Q(int(x), i)**2 for x in mean_i[:, source])
            self.assertEqual(norm2-projected2, 1-Q(1, i))
            records.append(dict(I=i, uniform_source_second_moment_norm=short(np.linalg.norm(k @ k @ w/np.sqrt(g), 2)),
                graph_departure_t2_coefficient=str(norm2-projected2), compression_exact=True))
        OBS['uniform_source_crosscheck'] = records

    def test_04_size_uniform_rational_error_and_general_concentration(self):
        t, factor, floor = certificate()
        # General collision-probability condition S <= s^2 < 1.
        cases = []
        for s in (Q(1, 4), Q(1, 2), Q(3, 4), Q(99, 100)):
            ts = min(Q(1, 256), (1-s)/132)
            actual = 1-48*ts
            mean = s+18*ts
            self.assertGreaterEqual(actual-mean, (1-s)/2)
            self.assertGreaterEqual(actual+mean, Q(1, 2))
            self.assertGreaterEqual(actual**2-mean**2, (1-s)/4)
            cases.append(dict(s=str(s), time_max=str(ts), gap_coefficient_lower=str((1-s)/4)))
        OBS['uniform_certificate'] = dict(time_max=str(t), rigorous_factor=str(factor),
            conservative_probability_gap=str(floor), ordinary_diamond_gap_lower=str(2*floor),
            all_I_at_least_2=True, graph_source_can_be_unknown_and_entangled=True,
            full_graph_H_norm_used_as_uniform_constant=False,
            general_parent_collision_cases=cases)

    def test_05_actual_port_probabilities_at_fixed_times(self):
        records = []
        for i in (2, 3):
            trees, _, _, _, h, w, bar, _ = model(i)
            n, g = 2*i+2, len(trees)
            for time in (Q(1, 1024), Q(1, 2048)):
                initial = np.zeros(n*g)
                initial[i*g:(i+1)*g] = 1/np.sqrt(g)
                actual = evolve(h, initial, time)
                mean = evolve(bar, np.eye(n)[:, i], time)
                p = float(np.vdot(actual[:i*g], actual[:i*g]).real)
                pm = float(np.vdot(mean[:i], mean[:i]).real)
                x = Q(4*(i-1)+6)*time
                tail = x**19/(math.factorial(19)*(1-x/20))
                gap = time*time/4
                self.assertGreater(p-pm-5*float(tail), float(gap))
                # Full source operator: minimum over all unknown graph states.
                columns = np.zeros((n*g, g))
                columns[i*g:(i+1)*g] = np.eye(g)
                block = evolve(h, columns, time)[:i*g]
                minimum = float(np.linalg.eigvalsh(block.conj().T @ block)[0])
                self.assertGreater(minimum-pm-5*float(tail), float(gap))
                records.append(dict(I=i, time=str(time), actual_uniform_source_probability=short(p),
                    coherent_mean_probability=short(pm), unknown_graph_minimum_probability=short(minimum),
                    probability_gap_lower=str(gap), numerical_evolution_tail_bound=str(tail)))
        OBS['actual_port_diagnostics'] = records

    def test_06_complete_instrument_and_passive_reference(self):
        trees, _, _, _, h, _, _, _ = model(2)
        n, g = 6, len(trees)
        u = evolve(h, np.eye(n*g), Q(1, 1024))
        completeness = np.linalg.norm(u.conj().T @ u-np.eye(n*g), 2)
        self.assertLess(completeness, 1e-12)
        # A fixed leaf and an unknown graph entangled with an internal reference.
        initial = np.zeros((n*g, 2), dtype=complex)
        initial[2*g, 0] = 1/np.sqrt(2)
        initial[2*g+1, 1] = 1j/np.sqrt(2)
        final = u @ initial
        reference_error = np.linalg.norm(final.conj().T @ final-initial.conj().T @ initial)
        self.assertLess(reference_error, 1e-12)
        probabilities = [float(np.linalg.norm(final[b*g:(b+1)*g])**2) for b in range(n)]
        self.assertAlmostEqual(sum(probabilities), 1, places=12)
        self.assertGreater(sum(probabilities[:2]), float(Q(1, 1024)**2*(1-Q(48, 1024))**2))
        OBS['instrument'] = dict(completeness_residual=short(completeness),
            passive_reference_residual=short(reference_error), all_port_outcomes_retained=True,
            internal_event_is_record_postprocessing=True, graph_not_measured=True)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(round=504, scientific_baseline_round=503,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'branching_tree_distance_audit.py', 'role_covariant_graph_obstruction.py', 'research_note_496.md')},
        scope=dict(original_440_generator=True, graph_operator_replacement_by_mean_only=True,
            fixed_time_uniform_size_probability_gap=True, arbitrary_actual_graph_reference=True,
            symmetric_source_preparation_derived=False, other_effective_geometries_excluded=False,
            dimension_three_generated=False, full_GR_goal_completed=False, phase_closure_triggered=False),
        observations=OBS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    data = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(data, ensure_ascii=False, indent=2))
