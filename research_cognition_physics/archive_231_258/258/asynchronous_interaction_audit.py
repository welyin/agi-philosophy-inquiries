"""Round 258: neighbor completion messages compile pair instruments asynchronously.

No global barrier or physical clock is used. Logical epochs, fixed adjacency,
reliable delivery, quantum memory and joint-gate capability remain inputs.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import random
import unittest
import numpy as np
import mutual_proposal_interaction_audit as pair

PATH4 = ((1,), (0, 2), (1, 3), (2,))


def q_settings(graph, epoch):
    return {v: pair.Z if epoch % 2 == 0 else pair.X for v in range(len(graph))}


def dependencies(graph, history):
    deps = {}
    for k, record in enumerate(history):
        for v, neighbors in enumerate(graph):
            deps['p', k, v] = (set() if k == 0 else
                {('f', k-1, v)} | {('s', k-1, w, v) for w in neighbors})
            for w in neighbors:
                deps['d', k, v, w] = {('p', k, v)}
            deps['c', k, v] = {('p', k, v)} | {('d', k, w, v) for w in neighbors}
        accepted = pair.matching(record)
        for u, v in accepted:
            deps['u', k, u, v] = {('c', k, u), ('c', k, v)}
        for v, neighbors in enumerate(graph):
            deps['f', k, v] = {('c', k, v)} | {('u', k, a, b) for a, b in accepted if v in (a, b)}
            for w in neighbors:
                deps['s', k, v, w] = {('f', k, v)}
    return deps


def ancestors(deps, event):
    result = set()
    pending = list(deps[event])
    while pending:
        current = pending.pop()
        if current not in result:
            result.add(current)
            pending.extend(deps[current])
    return result


def ordered(deps, seed=0, delay=None):
    generator = random.Random(seed)
    done = set()
    result = []
    while len(done) < len(deps):
        ready = [event for event, parents in deps.items() if event not in done and parents <= done]
        if not ready:
            raise ValueError('Deadlocked dependency graph.')
        preferred = [event for event in ready if event != delay]
        event = generator.choice(preferred or ready)
        done.add(event); result.append(event)
    return tuple(result)


def quantum_operators(graph, history, eta=.6, theta=np.pi/4):
    n = len(graph)
    result = {}
    for k, record in enumerate(history):
        for v in range(n):
            result['p', k, v] = pair.embed(n, {v: pair.proposal(graph, v, record[v], eta, q_settings(graph, k))})
        for a, b in pair.matching(record):
            result['u', k, a, b] = pair.pair_gate(n, (a, b), theta)
    return result


def evaluate(graph, history, order=None, theta=np.pi/4):
    ops = quantum_operators(graph, history, theta=theta)
    if order is None:
        order = tuple(event for k in range(len(history))
                      for event in ops if event[1] == k)
    value = np.eye(2**len(graph), dtype=complex)
    for event in order:
        if event in ops:
            value = ops[event]@value
    return value


def footprint(event):
    return set(event[2:]) if event[0] in ('p', 'u') else set()


def legal_order_count(deps):
    events = tuple(deps)
    positions = {event: k for k, event in enumerate(events)}
    parents = [sum(1 << positions[p] for p in deps[event]) for event in events]
    full = (1 << len(events))-1

    @lru_cache(maxsize=None)
    def count(done):
        if done == full:
            return 1
        return sum(count(done | (1 << k)) for k, needed in enumerate(parents)
                   if not done >> k & 1 and needed & done == needed)
    result = count(0)
    return {'events': len(events), 'linear_extensions': result, 'distinct_prefix_masks': count.cache_info().currsize}


@lru_cache(maxsize=None)
def metrics():
    rs = pair.records(PATH4)
    canonical = {}
    worst = 0.
    samples = 0
    for history in itertools.product(rs, repeat=2):
        deps = dependencies(PATH4, history)
        h = evaluate(PATH4, history)
        canonical[history] = h
        for seed in range(32):
            actual = evaluate(PATH4, history, ordered(deps, seed))
            worst = max(worst, float(np.max(np.abs(actual-h))))
            samples += 1
    complete = sum(h.conj().T@h for h in canonical.values())
    prefix = max(float(np.max(np.abs(sum(canonical[a, b].conj().T@canonical[a, b] for b in rs)
                 - (first := evaluate(PATH4, (a,))).conj().T@first))) for a in rs)
    psi = pair.tensor([pair.PLUS]*4)
    with_gate = [float(np.vdot(h@psi, h@psi).real) for h in canonical.values()]
    without_gate = []
    for history in canonical:
        h = evaluate(PATH4, history, theta=0.)
        without_gate.append(float(np.vdot(h@psi, h@psi).real))
    return {'two_epoch_histories': len(canonical), 'sampled_message_schedules': samples,
            'operator_schedule_max_error': worst,
            'completeness_max_error': float(np.max(np.abs(complete-np.eye(16)))),
            'second_epoch_prefix_max_error': prefix,
            'joint_record_total_variation_gate_vs_identity': float(np.abs(np.array(with_gate)-without_gate).sum()/2),
            'proposal_delivery_messages_per_epoch': 6,
            'completion_delivery_messages_per_epoch': 6,
            'both_epoch_message_deliveries': 24}


def report():
    graph = ((1,), (0, 2), (1, 3), (2, 4), (3,))
    history = ((1, 0, 1, 4, 3),)*2
    deps = dependencies(graph, history)
    early = ('p', 1, 0); late = ('p', 0, 4)
    order = ordered(deps, 123, delay=late)
    assert order.index(early) < order.index(late)
    return {'round': 258, 'eta': .6, 'path4': metrics(),
            'two_node_two_epoch_exact_schedule_count': legal_order_count(dependencies(((1,), (0,)), ((1, 0),)*2)),
            'no_global_barrier_witness': {'path_vertices': 5,
                'node0_epoch1_proposal_index': order.index(early),
                'node4_epoch0_proposal_index': order.index(late),
                'node0_can_advance_before_node4_starts': True},
            'scope': 'A finite local synchronizer implementation preserves full quantum instruments on a fixed supplied graph. No global barrier or wall clock; logical epochs, reliable eventual delivery, storage and pair-gate capabilities remain. No autonomous adjacency or physical time is derived.'}


class Audit(unittest.TestCase):
    def test_01_all_shared_quantum_support_is_causally_ordered(self):
        rs = pair.records(PATH4)
        for history in itertools.product(rs, repeat=2):
            deps = dependencies(PATH4, history)
            quantum = list(quantum_operators(PATH4, history))
            past = {e: ancestors(deps, e) for e in quantum}
            for a, b in itertools.combinations(quantum, 2):
                if footprint(a) & footprint(b):
                    self.assertTrue(a in past[b] or b in past[a])

    def test_02_all_pair_gates_have_delivered_mutual_proposals(self):
        for record in pair.records(PATH4):
            deps = dependencies(PATH4, (record,))
            for a, b in pair.matching(record):
                past = ancestors(deps, ('u', 0, a, b))
                self.assertTrue({('d', 0, a, b), ('d', 0, b, a)} <= past)

    def test_03_sampled_schedules_preserve_branch_operators(self):
        self.assertEqual(metrics()['sampled_message_schedules'], 512)
        self.assertLess(metrics()['operator_schedule_max_error'], 3e-14)

    def test_04_complete_two_epoch_instrument(self):
        self.assertLess(metrics()['completeness_max_error'], 3e-14)

    def test_05_complete_later_epoch_preserves_record_prefix(self):
        self.assertLess(metrics()['second_epoch_prefix_max_error'], 3e-14)

    def test_06_no_global_barrier(self):
        self.assertTrue(report()['no_global_barrier_witness']['node0_can_advance_before_node4_starts'])

    def test_07_exact_small_schedule_count(self):
        graph = ((1,), (0,)); history = ((1, 0),)
        deps = dependencies(graph, history)
        enumerated = sum(1 for _ in pair.legal_orders(deps))
        self.assertEqual(legal_order_count(deps)['linear_extensions'], enumerated)
        self.assertGreater(legal_order_count(dependencies(graph, history*2))['linear_extensions'], enumerated)

    def test_08_message_accounting_and_no_empty_queue_deadlock(self):
        for record in pair.records(PATH4):
            deps = dependencies(PATH4, (record,)*2)
            self.assertEqual(sum(e[0] in ('d', 's') for e in deps), 24)
            order = ordered(deps)
            self.assertEqual(len(order), len(deps))
            self.assertEqual(len(set(order)), len(deps))

    def test_09_gate_backreaction_on_next_records_is_retained(self):
        self.assertGreater(metrics()['joint_record_total_variation_gate_vs_identity'], .1)

    def test_10_waiting_for_neighbor_completion_is_explicit(self):
        history = (pair.records(PATH4)[0],)*2
        deps = dependencies(PATH4, history)
        for v, neighbors in enumerate(PATH4):
            past = ancestors(deps, ('p', 1, v))
            self.assertTrue({('f', 0, w) for w in neighbors} <= past)
        # With no fairness assumption an enabled message can be postponed forever;
        # our theorem concerns completed finite histories, not a time upper bound.


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
    target = Path(__file__).with_name('asynchronous_interaction_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
