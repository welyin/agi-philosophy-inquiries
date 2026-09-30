"""Round 506: exact joint graph/port orbit code with coarse role instruments.

The declared symmetric source and coherent role measurement are extra inputs.
Code dimension is not physical spatial dimension. Fine label instruments differ.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import io
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import coherent_graph_mean_obstruction as old

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_orbit_port_interface_results.json'
OBS = {}


def rooted_key(adj, v, parent=-1):
    return '('+''.join(sorted(rooted_key(adj, w, v) for w in adj[v] if w != parent))+')'


@lru_cache(None)
def quotient(i):
    trees = old.previous.ensemble(i)[0]
    n, g = 2*i+2, len(trees)
    lookup = {tree: j for j, tree in enumerate(trees)}
    adjs = [old.tree_tools.adjacency(tree, n) for tree in trees]
    keys = [[rooted_key(adj, v) for j, adj in enumerate(adjs)] for v in range(n)]
    names = sorted({key for row in keys for key in row})
    index = {key: j for j, key in enumerate(names)}
    labels = np.array([[index[key] for key in row] for row in keys], dtype=np.int64)
    q = len(names)
    sizes = np.bincount(labels.ravel(), minlength=q)
    role = np.zeros(q, dtype=np.int64)-1
    for v in range(n):
        for a in labels[v]:
            if role[a] == -1:
                role[a] = int(v < i)
            assert role[a] == int(v < i)
    counts = np.zeros((n*g, q), dtype=np.int64)
    for j, tree in enumerate(trees):
        flips = old.tree_tools.tree_flips(tree, n)
        for v in range(n):
            x = v*g+j
            counts[x, labels[v, j]] -= len(adjs[j][v])
            for w in adjs[j][v]:
                counts[x, labels[w, j]] += 1
            for other, weight in flips.items():
                counts[x, labels[v, lookup[other]]] += weight
    m = np.stack([counts[np.flatnonzero(labels.ravel() == a)[0]] for a in range(q)])
    exact = bool(np.array_equal(counts, m[labels.ravel()]))
    assert exact
    symmetric_counts = sizes[:, None]*m
    assert np.array_equal(symmetric_counts, symmetric_counts.T)
    k = symmetric_counts/np.sqrt(sizes[:, None]*sizes[None, :])
    return dict(I=i, n=n, g=g, q=q, names=names, labels=labels, sizes=sizes,
                role=role, m=m, k=k, counts=counts)


def embedding(data):
    labels = data['labels'].ravel()
    e = np.zeros((len(labels), data['q']))
    e[np.arange(len(labels)), labels] = 1/np.sqrt(data['sizes'][labels])
    return e


def uniform_leaf_source(data):
    return np.sqrt(data['sizes']/(data['g']*(data['I']+2)))*(1-data['role'])


def sequence(h, projectors, times, outcomes):
    value = np.eye(len(h), dtype=complex)
    for t, outcome in zip(times, outcomes):
        value = projectors[outcome] @ old.evolve(h, value, t, order=36)
    return value


class Audit(unittest.TestCase):
    def test_01_complete_joint_orbits_and_integer_intertwining(self):
        records = []
        for i, q in ((2, 2), (3, 4), (4, 7)):
            d = quotient(i)
            self.assertEqual(d['q'], q)
            self.assertEqual(int(d['sizes'].sum()), d['n']*d['g'])
            # h Z = Z M checked on every original basis row, with integer sums.
            self.assertTrue(np.array_equal(d['counts'], d['m'][d['labels'].ravel()]))
            records.append(dict(I=i, original_dimension=d['n']*d['g'], code_dimension=q,
                orbit_sizes=d['sizes'].tolist(), internal_role_flags=d['role'].tolist(),
                rooted_tree_signatures=d['names'], unnormalized_quotient=d['m'].tolist(),
                integer_intertwining_exact=True))
        OBS['joint_codes'] = records

    def test_02_full_hamiltonian_and_coarse_instrument_intertwining(self):
        records = []
        for i in (2, 3):
            d = quotient(i)
            e = embedding(d)
            h = old.model(i)[4]
            qint = np.diag(d['role'])
            pint = np.diag(np.repeat(np.arange(d['n']) < i, d['g']).astype(float))
            self.assertLess(np.linalg.norm(e.T @ e-np.eye(d['q'])), 1e-12)
            self.assertLess(np.linalg.norm(h @ e-e @ d['k']), 1e-12)
            self.assertLess(np.linalg.norm(pint @ e-e @ qint), 1e-12)
            self.assertLess(np.linalg.norm((np.eye(len(h))-pint) @ e-e @ (np.eye(d['q'])-qint)), 1e-12)
            source = e @ uniform_leaf_source(d)
            target = np.zeros(d['n']*d['g'])
            target[i*d['g']:] = 1/np.sqrt((i+2)*d['g'])
            self.assertLess(np.linalg.norm(source-target), 1e-12)
            records.append(dict(I=i, all_code_columns_checked=True, source_inside_code=True,
                joint_graph_information_retained=True, coherent_role_instrument_intertwines=True))
        OBS['full_interface'] = records

    def test_03_unknown_code_reference_and_all_history_branches(self):
        d = quotient(3)
        e = embedding(d)
        h = old.model(3)[4]
        pint = np.diag(np.repeat(np.arange(d['n']) < 3, d['g']).astype(float))
        physical = [np.eye(len(h))-pint, pint]
        qint = np.diag(d['role'])
        coarse = [np.eye(d['q'])-qint, qint]
        times = [Q(1, 32), Q(1, 16)]
        completeness = np.zeros((d['q'], d['q']), dtype=complex)
        maximum = 0.0
        for a in range(2):
            for b in range(2):
                k = sequence(d['k'], coarse, times, [a, b])
                full = old.evolve(h, e, times[0], order=36)
                full = physical[a] @ full
                full = physical[b] @ old.evolve(h, full, times[1], order=36)
                residual = np.linalg.norm(full-e @ k)
                maximum = max(maximum, float(residual))
                self.assertLess(residual, 1e-12)
                completeness += k.conj().T @ k
        self.assertLess(np.linalg.norm(completeness-np.eye(d['q'])), 1e-12)
        # All columns above imply arbitrary unknown input/reference intertwining.
        OBS['histories'] = dict(branches=4, waits=[str(t) for t in times],
            all_unknown_code_columns_preserved=True, max_intertwining_residual=old.short(maximum),
            completeness_residual=old.short(np.linalg.norm(completeness-np.eye(d['q']))),
            passive_reference_included_by_tensor_identity=True)

    def test_04_label_probabilities_agree_but_fine_poststate_leaks(self):
        records = []
        for i in (2, 3):
            d = quotient(i)
            e = embedding(d)
            qint = np.diag(d['role'])
            expected_leak = Q(i+1, i+2)
            source = e @ uniform_leaf_source(d)
            after_support = 0.0
            for v in range(d['n']):
                ev = e[v*d['g']:(v+1)*d['g']]
                compressed = ev.T @ ev
                target = qint/i if v < i else (np.eye(d['q'])-qint)/(i+2)
                self.assertLess(np.linalg.norm(compressed-target), 1e-12)
                projected_source = np.zeros_like(source)
                projected_source[v*d['g']:(v+1)*d['g']] = source[v*d['g']:(v+1)*d['g']]
                after_support += np.linalg.norm(e.T @ projected_source)**2
            self.assertAlmostEqual(1-after_support, float(expected_leak), places=12)
            records.append(dict(I=i, all_single_label_effects_equal_role_uniform_split=True,
                fine_nonselective_instrument_leakage=str(expected_leak), coherent_role_read_leakage='0'))
        OBS['fine_instrument_boundary'] = records

    def test_05_two_step_actual_probability_difference(self):
        t = Q(1, 4096)
        factor = Q(1, 2)-240*t
        self.assertEqual(factor, Q(113, 256))
        self.assertGreater(factor, Q(3, 8))
        rows = []
        for i in (2, 3):
            d = quotient(i)
            e = embedding(d)
            h = old.model(i)[4]
            src = e @ uniform_leaf_source(d)
            full = old.evolve(h, src, t)
            coherent = float(np.linalg.norm(full[:i*d['g']])**2)
            # Read fine ports at time zero, keep every branch and then wait.
            branches = np.zeros((len(h), i+2))
            for c, v in enumerate(range(i, d['n'])):
                branches[v*d['g']:(v+1)*d['g'], c] = src[v*d['g']:(v+1)*d['g']]
            full_branches = old.evolve(h, branches, t)
            fine = float(np.linalg.norm(full_branches[:i*d['g']])**2)
            self.assertGreater(coherent-fine, float(Q(3, 8)*t*t))
            rows.append(dict(I=i, coherent_role_then_internal_probability=old.short(coherent),
                fine_read_merge_then_internal_probability=old.short(fine),
                actual_difference=old.short(coherent-fine)))
        OBS['instrument_operational_certificate'] = dict(time_max=str(t),
            all_I_at_least_2=True, factor_lower=str(factor), conservative_gap=str(Q(3, 8)*t*t),
            initial_role_outcome_certain_for_both=True, diagnostics=rows)

    def test_06_joint_code_exceeds_two_modes_and_keeps_relations(self):
        d = quotient(3)
        source = uniform_leaf_source(d)
        ui = np.sqrt(d['sizes']/(d['g']*3))*d['role']
        p = np.outer(source, source)+np.outer(ui, ui)
        velocity = d['k'] @ source
        leakage = float(np.linalg.norm((np.eye(d['q'])-p) @ velocity)**2)
        self.assertAlmostEqual(leakage, float(Q(2, 15)), places=12)
        # The missing directions live inside the exact joint code, not outside it.
        self.assertGreater(leakage, 0)
        OBS['relation_memory'] = dict(I=3, exact_joint_dimension=4,
            two_role_velocity_leakage_squared='2/15', leakage_stays_inside_joint_code=True,
            code_dimension_identified_with_spatial_dimension=False,
            arbitrary_unknown_graph_sources_covered=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=506, scientific_baseline_round=505,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'coherent_graph_mean_obstruction.py', 'role_covariant_graph_obstruction.py',
            'branching_tree_distance_audit.py', 'research_note_496.md', 'research_note_505.md')},
        scope=dict(joint_graph_port_code=True, all_unknown_code_reference_inputs=True,
            coherent_role_instrument_exact=True, arbitrary_finite_histories_intertwine=True,
            fine_merge_instrument_equivalent=False, arbitrary_unknown_graph_inputs_covered=False,
            symmetric_source_derived=False, coherent_role_reader_implemented=False,
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
