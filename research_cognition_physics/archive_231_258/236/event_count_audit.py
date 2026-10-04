"""Round 236: event refinement, boundary channels and the missing volume input.

A weight defined only on finite-dimensional boundary channels cannot be both
finite, additive for serial composition and nonzero on a monoid containing a
reset channel. This is an algebraic obstruction, not a thermodynamic claim.
"""
import argparse
import json
from pathlib import Path
import platform
import unittest
import numpy as np


def units(d):
    for i in range(d):
        for j in range(d):
            e = np.zeros((d, d), dtype=complex)
            e[i, j] = 1
            yield e


def random_unitary(d, rng):
    q, _ = np.linalg.qr(rng.normal(size=(d, d))+1j*rng.normal(size=(d, d)))
    return q


def factorization_certificate():
    rng = np.random.default_rng(236001)
    cases = []
    for d in [2, 3, 4]:
        target = random_unitary(d, rng)
        for count in [1, 2, 5, 10]:
            factors = []
            product = np.eye(d, dtype=complex)
            for _ in range(count-1):
                v = random_unitary(d, rng)
                factors.append(v)
                product = v@product
            factors.append(target@product.conj().T)
            composed = np.eye(d, dtype=complex)
            for v in factors:
                composed = v@composed
            channel_error = max(np.linalg.norm(composed@e@composed.conj().T-target@e@target.conj().T)
                                for e in units(d))
            cases.append({'dimension': d, 'operation_vertices': count,
                          'same_boundary_channel_error': float(channel_error)})
    return cases


def identity_refinements():
    v = np.eye(2).reshape(-1)/np.sqrt(2)
    rho = np.outer(v, v)
    cases = []
    identity_on_system_and_reference = np.kron(np.eye(2), np.eye(2))
    for n in [1, 2, 5, 10, 50]:
        evolved = rho.copy()
        for _ in range(n):
            evolved = identity_on_system_and_reference@evolved@identity_on_system_and_reference.T
        cases.append({'identity_vertices': n,
                      'boundary_reference_error': float(np.linalg.norm(evolved-rho)),
                      'naive_volume_at_assumed_density_10': n/10})
    return cases


def absorbing_certificate():
    rng = np.random.default_rng(236002)
    errors, tp_errors = [], []
    for d in [2, 3, 4]:
        for _ in range(8):
            iso, _ = np.linalg.qr(rng.normal(size=(3*d, d))+1j*rng.normal(size=(3*d, d)))
            kraus = iso.reshape(3, d, d)
            phi = lambda a:sum(k@a@k.conj().T for k in kraus)
            reset = lambda a:np.trace(a)*np.eye(d)/d
            tp_errors.append(np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(d)))
            errors.append(max(np.linalg.norm(reset(phi(a))-reset(a)) for a in units(d)))
    return {'random_CPTP_examples': len(errors),
            'maximum_trace_preserving_error': float(max(tp_errors)),
            'maximum_reset_after_channel_vs_reset_error': float(max(errors))}


def monoid_certificate():
    # Four deterministic stochastic maps on a bit. All sixteen serial products
    # stay in this monoid. c(g o f)-c(g)-c(f)=0 gives a homogeneous linear system.
    maps = [(0, 1), (1, 0), (0, 0), (1, 1)]
    labels = ['identity', 'flip', 'reset_0', 'reset_1']
    constraints, table = [], []
    for g, go in enumerate(maps):
        row = []
        for f, fo in enumerate(maps):
            h = maps.index(tuple(go[fo[i]] for i in [0, 1]))
            row.append(h)
            equation = np.zeros(4)
            equation[h] += 1
            equation[g] -= 1
            equation[f] -= 1
            constraints.append(equation)
        table.append(row)
    a = np.array(constraints)
    rank = int(np.linalg.matrix_rank(a))
    return {'maps': labels, 'composition_table_g_after_f': table,
            'equations': a.astype(int).tolist(), 'rank': rank,
            'finite_real_additive_weight_nullity': len(maps)-rank,
            'unit_event_count_max_additivity_defect': float(np.max(np.abs(a@np.ones(4))))}


def chain_certificate():
    cases = []
    for interior in [1, 2, 4, 8, 16, 32]:
        n = interior+2
        relation = np.triu(np.ones((n, n), dtype=bool), 1)
        # Each finite stage has endpoints a,b and n-2 intermediate events.
        # Inclusion inserts new events before b; all arrows carry identity M2.
        transitive = not np.any((relation.astype(int)@relation.astype(int)>0) & ~relation)
        between = int(np.sum(relation[0, :] & relation[:, -1]))
        cases.append({'interior_events': between, 'total_events': n,
                      'strict_partial_order': bool(transitive and not np.any(np.diag(relation))),
                      'each_slice_system_dimension': 2})
    return {'finite_stages': cases,
            'analytic_direct_limit': 'a < z_1 < z_2 < ... < b; I(a,b) is countably infinite',
            'finite_samples_do_not_prove_limit_locally_finite': True}


def tagged_positive_control():
    # Added duration is a model input, not a number extracted from the channel.
    x = np.array([[0., 1.], [1., 0.]])
    return {'two_flip_product_identity_error': float(np.linalg.norm(x@x-np.eye(2))),
            'two_flip_added_duration': 2.,
            'empty_identity_duration': 0.,
            'tagged_processes_equal': False,
            'identity_step_durations': [1., 2., 5., 10., 50.],
            'assumption': 'Duration/history is retained in the interface; boundary channel equality no longer identifies the processes.'}


def report():
    return {'round': 236, 'date': '2026-09-22',
            'hypothesis': 'Can a positive additive event-volume be extracted from boundary-channel equivalence alone?',
            'unitary_factorizations': factorization_certificate(),
            'identity_refinements': identity_refinements(),
            'absorbing_reset': absorbing_certificate(),
            'finite_monoid': monoid_certificate(),
            'finite_vs_local_finiteness': chain_certificate(),
            'added_tag_positive_control': tagged_positive_control(),
            'conclusion': 'Boundary channels and finite-dimensionality do not fix an atomic event count, volume density or local finiteness of an infinite history.',
            'not_claimed': ['Physical idle time has zero volume or zero cost',
                            'Different microscopic histories are physically identical',
                            'Causal-set counting is inconsistent',
                            'Identity refinement is a gauge symmetry of nature'],
            'next': 'Choose a conditional history model with retained atomic events and test existing geometric estimators; independently track how event identity and counting could be motivated.',
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__}}


class EventCountTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = report()

    def test_unitary_factorizations_change_counts_preserve_channel(self):
        self.assertEqual(len(self.r['unitary_factorizations']), 12)
        for c in self.r['unitary_factorizations']:
            self.assertLess(c['same_boundary_channel_error'], 3e-13)

    def test_identity_refinement_and_assumed_volume(self):
        cases = self.r['identity_refinements']
        self.assertEqual([c['identity_vertices'] for c in cases], [1, 2, 5, 10, 50])
        self.assertTrue(all(c['boundary_reference_error'] == 0 for c in cases))
        self.assertEqual([c['naive_volume_at_assumed_density_10'] for c in cases], [.1, .2, .5, 1., 5.])

    def test_reset_absorbs_all_test_channels(self):
        c = self.r['absorbing_reset']
        self.assertEqual(c['random_CPTP_examples'], 24)
        self.assertLess(c['maximum_trace_preserving_error'], 3e-13)
        self.assertLess(c['maximum_reset_after_channel_vs_reset_error'], 3e-13)

    def test_exact_finite_composition_constraints(self):
        c = self.r['finite_monoid']
        self.assertEqual(c['rank'], 4)
        self.assertEqual(c['finite_real_additive_weight_nullity'], 0)
        self.assertEqual(c['unit_event_count_max_additivity_defect'], 1.)
        self.assertEqual(c['composition_table_g_after_f'][2], [2]*4)
        self.assertEqual(c['composition_table_g_after_f'][3], [3]*4)

    def test_finite_stages_have_unbounded_intervals(self):
        cases = self.r['finite_vs_local_finiteness']['finite_stages']
        self.assertTrue(all(c['strict_partial_order'] for c in cases))
        self.assertEqual([c['interior_events'] for c in cases], [1, 2, 4, 8, 16, 32])

    def test_added_tag_removes_the_assumed_equivalence(self):
        c = self.r['added_tag_positive_control']
        self.assertEqual(c['two_flip_product_identity_error'], 0.)
        self.assertGreater(c['two_flip_added_duration'], c['empty_identity_duration'])
        self.assertFalse(c['tagged_processes_equal'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(EventCountTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('event_count_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != json.loads(json.dumps(data)):
            raise RuntimeError('Existing result differs; inspect before replacing.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))
