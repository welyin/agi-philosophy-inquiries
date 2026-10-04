"""Round 255: complete local proposals, mutual matching, and pair interactions.

The candidate graph and reliable classical handshake are supplied resources.
Only the logical qubits are enumerated; protected R flags are audited separately.
No maximal-matching, autonomous encounter, or spacetime claim is made.
"""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from record_preserving_feedback_audit import weak

I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])
PLUS = np.ones(2)/np.sqrt(2)
CYCLE = ((1, 3), (0, 2), (1, 3), (0, 2))


def tensor(items):
    result = np.array([1.])
    for item in items:
        result = np.kron(result, item)
    return result


def embed(n, local):
    return tensor([local.get(v, I) for v in range(n)])


def validate(graph):
    for v, neighbors in enumerate(graph):
        if len(neighbors) > 2 or len(set(neighbors)) != len(neighbors):
            raise ValueError('This implementation supports simple graphs of degree <= 2.')
        for w in neighbors:
            if w == v or w < 0 or w >= len(graph) or v not in graph[w]:
                raise ValueError('Candidate edges must be undirected, without self loops.')


def records(graph):
    validate(graph)
    return tuple(itertools.product(*(neighbors or (None,) for neighbors in graph)))


def proposal(graph, v, choice, eta, settings=None):
    neighbors = graph[v]
    if choice not in (neighbors or (None,)):
        raise ValueError('Invalid proposal.')
    if len(neighbors) <= 1:
        return I.copy()
    if not 0 <= eta <= 1:
        raise ValueError('Invalid measurement strength.')
    q = (X if v % 2 == 0 else Z) if settings is None else settings[v]
    return weak(q, 1 if choice == neighbors[0] else -1, eta)


def matching(record):
    return tuple((v, w) for v, w in enumerate(record)
                 if w is not None and v < w and record[w] == v)


def pair_gate(n, edge, theta, axis=Z):
    a, b = edge
    return np.cos(theta)*np.eye(2**n) - 1j*np.sin(theta)*embed(n, {a: axis, b: axis})


def event_data(graph, record, eta=.6, theta=np.pi/4, settings=None):
    n = len(graph)
    operators = {('p', v): embed(n, {v: proposal(graph, v, record[v], eta, settings)})
                 for v in range(n)}
    dependencies = {event: set() for event in operators}
    for a, b in matching(record):
        event = ('u', a, b)
        operators[event] = pair_gate(n, (a, b), theta)
        dependencies[event] = {('p', a), ('p', b)}
    return operators, dependencies


def legal_orders(dependencies, done=()):
    if len(done) == len(dependencies):
        yield done
        return
    seen = set(done)
    for event, parents in dependencies.items():
        if event not in seen and parents <= seen:
            yield from legal_orders(dependencies, done+(event,))


def branch(graph, record, eta=.6, theta=np.pi/4, settings=None):
    result = tensor([proposal(graph, v, choice, eta, settings)
                     for v, choice in enumerate(record)])
    for edge in matching(record):
        result = pair_gate(len(graph), edge, theta) @ result
    return result


def cycle_metrics():
    recs = records(CYCLE)
    hs = [branch(CYCLE, r) for r in recs]
    total = sum(h.conj().T@h for h in hs)
    order_count = 0
    diamond_error = 0.
    for rec, canonical in zip(recs, hs):
        ops, deps = event_data(CYCLE, rec)
        for order in legal_orders(deps):
            value = np.eye(16, dtype=complex)
            for event in order:
                value = ops[event] @ value
            diamond_error = max(diamond_error, float(np.max(np.abs(value-canonical))))
            order_count += 1
    prefix_error = max(float(np.max(np.abs(
        sum((b@a).conj().T@(b@a) for b in hs) - a.conj().T@a))) for a in hs)
    counts = Counter(matching(r) for r in recs)
    distribution = [{'edges': list(key), 'proposal_count': value, 'mixed_state_probability': value/16}
                    for key, value in sorted(counts.items())]
    return {'proposal_records': len(recs), 'distinct_matchings': len(counts),
            'legal_operator_orders': order_count,
            'completeness_max_error': float(np.max(np.abs(total-np.eye(16)))),
            'order_max_error': diamond_error, 'two_epoch_prefix_max_error': prefix_error,
            'maximally_mixed_matching_distribution': distribution}


def report():
    psi = pair_gate(2, (0, 1), np.pi/4) @ tensor([PLUS, PLUS])
    reduced = psi.reshape(2, 2) @ psi.reshape(2, 2).conj().T
    u = pair_gate(3, (0, 1), np.pi/4, Z)
    v = pair_gate(3, (1, 2), np.pi/4, X)
    return {'round': 255, 'eta': .6, 'theta': float(np.pi/4),
            'cycle': cycle_metrics(),
            'selected_pair_reduced_purity': float(np.trace(reduced@reduced).real),
            'overlapping_ZZ_XX_order_max_difference': float(np.max(np.abs(u@v-v@u))),
            'scope': 'Supplied degree-at-most-two encounter graph; local complete proposal instruments; reliable classical handshakes; two-input two-output gates; complete negotiation epochs. No physical geometry is derived.'}


class Audit(unittest.TestCase):
    def test_01_graphs_and_matching_exclusion(self):
        for graph in (CYCLE, ((1,), (0, 2), (1,)), ((),), ((1,), (0,))):
            for r in records(graph):
                endpoints = [v for e in matching(r) for v in e]
                self.assertEqual(len(endpoints), len(set(endpoints)))
        for invalid in (((1,), ()), ((0,),)):
            with self.assertRaises(ValueError):
                records(invalid)

    def test_02_complete_instrument_on_arbitrary_input(self):
        for eta in (0., .6, 1.):
            for graph in (CYCLE, ((1,), (0, 2), (1,)), ((),)):
                hs = [branch(graph, r, eta) for r in records(graph)]
                np.testing.assert_allclose(sum(h.conj().T@h for h in hs),
                                           np.eye(2**len(graph)), atol=2e-14)

    def test_03_all_legal_interleavings(self):
        result = cycle_metrics()
        self.assertEqual(result['legal_operator_orders'], 688)
        self.assertLess(result['order_max_error'], 2e-14)

    def test_04_matching_probabilities_independent_count(self):
        counts = Counter(matching(r) for r in records(CYCLE))
        self.assertEqual(sorted(counts.values()), [1, 1, 2, 3, 3, 3, 3])
        for eta in (0., .6, 1.):
            actual = Counter()
            for rec in records(CYCLE):
                h = branch(CYCLE, rec, eta)
                actual[matching(rec)] += float(np.trace(h.conj().T@h).real/16)
            for m, c in counts.items():
                self.assertAlmostEqual(actual[m], c/16)
        self.assertEqual(counts[()], 2)  # Empty output is retained, not postselected away.

    def test_05_protected_flags_full_local_spaces(self):
        for q in (X, Z):
            ks = [np.kron(I, weak(q, g, .6)) for g in (-1, 1)]
            flag = np.kron(Z, I)
            np.testing.assert_allclose(sum(k.conj().T@flag@k for k in ks), flag, atol=2e-14)
            for k in ks:
                np.testing.assert_allclose(k@flag, flag@k, atol=2e-14)
        gate = pair_gate(4, (1, 3), np.pi/4)  # R0,L0,R1,L1
        for r in (0, 2):
            flag = embed(4, {r: Z})
            np.testing.assert_allclose(gate.conj().T@flag@gate, flag, atol=2e-14)

    def test_06_pair_can_entangle_distinct_branches(self):
        self.assertAlmostEqual(report()['selected_pair_reduced_purity'], .5)

    def test_07_pair_gate_preserves_unknown_reference_capacity(self):
        # Maximal entanglement with a four-dimensional reference: recover with U dagger.
        psi = np.eye(4).reshape(-1)/2
        u = pair_gate(2, (0, 1), .371)
        evolved = np.kron(u, np.eye(4))@psi
        np.testing.assert_allclose(np.kron(u.conj().T, np.eye(4))@evolved, psi, atol=2e-14)
        reduced = evolved.reshape(4, 4)@evolved.reshape(4, 4).conj().T
        np.testing.assert_allclose(reduced, np.eye(4)/4, atol=2e-14)

    def test_08_complete_second_epoch_preserves_first_effect(self):
        self.assertLess(cycle_metrics()['two_epoch_prefix_max_error'], 2e-14)

    def test_09_shared_port_orders_are_not_generally_gauge(self):
        self.assertGreater(report()['overlapping_ZZ_XX_order_max_difference'], .9)

    def test_10_isolated_and_single_edge_resources(self):
        self.assertEqual(records(((),)), ((None,),))
        np.testing.assert_allclose(branch(((),), (None,)), I)
        graph = ((1,), (0,))
        self.assertEqual(records(graph), ((1, 0),))
        self.assertEqual(matching((1, 0)), ((0, 1),))
        np.testing.assert_allclose(branch(graph, (1, 0)), pair_gate(2, (0, 1), np.pi/4))


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
    target = Path(__file__).with_name('mutual_proposal_interaction_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve existing results; use a new result version.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
