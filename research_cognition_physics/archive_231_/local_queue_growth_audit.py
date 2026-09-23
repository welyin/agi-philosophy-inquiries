"""Round 273: local fluid backlog can trigger bounded-degree paired growth.

Task splitting, free local migration at an atomic rewrite, primitive service
and construction credits are supplied. This is not an unknown-state protocol.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from paired_interface_growth_audit import refine, check_ports, rail_order
from ancestral_quotient_audit import graph_data
from vertex_split_cycle_audit import statistics
from fisher_hanoi_scaling_audit import bfs


def needed_channels(arrival, service=F(1)):
    count = 1
    while arrival > service*count:
        count *= 2
    return count


@lru_cache(maxsize=None)
def simulate(arrival=F(12), threshold=F(12), service=F(1), ticks=256, credits=None):
    arrival, threshold, service = map(F, (arrival, threshold, service))
    if not (arrival > 0 and threshold >= arrival and service > 0):
        raise ValueError('Require positive rates and threshold at least initial arrival.')
    state = {(0,): {1: ((1,), 1)}, (1,): {1: ((0,), 1)}}
    # One classical fluid queue on the left endpoint of each designated rung.
    queues = {(0,): (F(0), arrival)}
    history, growth_times = [], []
    arrived = delivered = F(0)
    rewrites = 0
    for t in range(1, ticks+1):
        before_count = len(queues)
        new_queues = {}
        proposals = []
        slot_arrival = slot_delivery = F(0)
        for u, (q, rate) in sorted(queues.items()):
            slot_arrival += rate
            sent = min(q+rate, service)
            slot_delivery += sent
            remainder = q+rate-sent
            new_queues[u] = (remainder, rate)
            if remainder > threshold:
                proposals.append((u, state[u][1][0]))
        selected = proposals if credits is None else proposals[:max(0, credits-rewrites)]
        for u, v in selected:
            q, rate = new_queues.pop(u)
            refine(state, u, v)
            new_queues[u+(0,)] = (q/2, rate/2)
            new_queues[u+(1,)] = (q/2, rate/2)
        queues = new_queues
        rewrites += len(selected)
        if selected:
            growth_times.append(t)
        arrived += slot_arrival
        delivered += slot_delivery
        check_ports(state)
        assert sum(rate for q, rate in queues.values()) == arrival
        history.append({'t': t, 'channels_before': before_count, 'channels': len(queues),
                        'splits': len(selected), 'arrived': arrived,
                        'delivered': delivered, 'queued': sum(q for q, rate in queues.values()),
                        'max_local_queue': max(q for q, rate in queues.values()),
                        'slot_delivery': slot_delivery})
    return {'state': state, 'queues': queues, 'history': history,
            'growth_times': growth_times, 'rewrites': rewrites}


def geometry(state):
    first = rail_order(state)
    second = [state[u][1][0] for u in first]
    coordinates = {w: (side, i) for side, row in enumerate((first, second))
                   for i, w in enumerate(row)}
    assert len(coordinates) == len(state)
    actual = {tuple(sorted((coordinates[u], coordinates[v])))
              for u, row in state.items() for v, _ in row.values() if u < v}
    length = len(first)
    expected = {((0, i), (1, i)) for i in range(length)}
    expected |= {((side, i), (side, i+1)) for side in (0, 1) for i in range(length-1)}
    assert actual == expected
    prefix = set(first[:length//2]+second[:length//2])
    boundary = sum((u in prefix) != (v in prefix)
                   for u, row in state.items() for v, _ in row.values() if u < v)
    graph = graph_data(state)[0]
    data_owner = min(w for w in state if w[0] == 0)
    return {**statistics(graph), 'rungs': length,
            'diameter': max(max(bfs(graph, i)) for i in range(len(graph))),
            'longitudinal_half_cut': boundary,
            'original_register_vertex_degree': len(state[data_owner])}


def summary(arrival, threshold, credits=None):
    run = simulate(F(arrival), F(threshold), credits=credits)
    last = run['history'][-1]
    return {'total_arrival_per_slot': float(arrival), 'threshold': float(threshold),
            'construction_credits': credits, 'growth_times': run['growth_times'],
            'rewrites': run['rewrites'], 'new_registers': 2*run['rewrites'],
            'new_primitive_channels_net': 3*run['rewrites'],
            'final_queue': float(last['queued']),
            'largest_local_queue_seen': float(max(row['max_local_queue'] for row in run['history'])),
            'last_slot_delivery': float(last['slot_delivery']),
            'geometry': geometry(run['state'])}


def report():
    return {
        'round': 273,
        'logical_slots': 256,
        'local_rule': 'After arrivals and service, refine a designated rung if its local queue exceeds H; split classical queue and future fluid arrival equally between its children.',
        'unlimited_construction': [summary(F(rate), F(rate)) for rate in (1, 3, 8, 12)],
        'threshold_controls_delay_not_final_channel_count': [summary(F(12), F(h)) for h in (12, 24, 48)],
        'finite_construction': [summary(F(12), F(12), credits=b) for b in (0, 1, 3)],
        'fixed_endpoint_warning': 'At 16 rungs aggregate rung capacity is 16 but the child-zero original data vertex has only 3 incident edges. The divided local workload cannot be reinterpreted as rate-12 transfer from that one fixed source.',
        'scope': 'A conditional classical fluid-work growth controller with exact conservation, bounded queues and explicit construction costs; its aligned output remains a ladder. No spontaneous resources, physical time, quantum task migration, throughput-optimal general routing or selected spatial dimension is claimed.'
    }


class Audit(unittest.TestCase):
    def test_01_all_arrivals_are_delivered_or_remain_queued(self):
        for rate in (F(1, 2), F(1), F(3), F(8), F(12)):
            for row in simulate(rate, rate)['history']:
                self.assertEqual(row['arrived'], rate*row['t'])
                self.assertEqual(row['delivered']+row['queued'], row['arrived'])
                self.assertLessEqual(row['slot_delivery'], row['channels_before'])

    def test_02_local_queue_bound_through_rewrites(self):
        for rate in (F(3, 2), F(3), F(8), F(12)):
            for scale in (1, 2, 4):
                run = simulate(rate, scale*rate)
                self.assertLessEqual(max(r['max_local_queue'] for r in run['history']), scale*rate)

    def test_03_eventual_channel_count_including_exact_service_threshold(self):
        for service in (F(1), F(3, 2)):
            for rate in (F(1, 4), F(1), F(3, 2), F(2), F(3), F(8), F(17, 2), F(12)):
                run = simulate(rate, rate, service)
                count = len(run['queues'])
                self.assertEqual(count, needed_channels(rate, service))
                self.assertTrue(all(r <= service for q, r in run['queues'].values()))
                if rate > service:
                    self.assertGreater(2*rate, count*service)
                    self.assertGreaterEqual(count*service, rate)

    def test_04_actual_graph_counts_match_the_paid_rewrites(self):
        for rate in (F(1), F(3), F(8), F(12)):
            run = simulate(rate, rate)
            data = statistics(graph_data(run['state'])[0])
            self.assertEqual(data['vertices'], 2+2*run['rewrites'])
            self.assertEqual(data['edges'], 1+3*run['rewrites'])
            self.assertEqual(data['cycle_rank'], run['rewrites'])
            self.assertEqual(data['components'], 1)
            self.assertLessEqual(data['max_degree'], 3)

    def test_05_ladder_isomorphism_and_distances_are_recovered_from_edges(self):
        for rate in (F(3), F(8), F(12)):
            data = geometry(simulate(rate, rate)['state'])
            self.assertEqual(data['diameter'], data['rungs'])
            self.assertEqual(data['longitudinal_half_cut'], 2)

    def test_06_threshold_changes_waiting_time_without_encoding_dimension(self):
        runs = [simulate(F(12), F(h)) for h in (12, 24, 48)]
        self.assertEqual({len(r['queues']) for r in runs}, {16})
        times = [r['growth_times'][-1] for r in runs]
        self.assertLess(times[0], times[1])
        self.assertLess(times[1], times[2])

    def test_07_finite_credits_leave_linear_unserved_work(self):
        for credits in (0, 1, 3):
            run = simulate(F(12), F(12), credits=credits)
            self.assertEqual(run['rewrites'], credits)
            self.assertEqual(len(run['queues']), credits+1)
            history = run['history']
            for a, b in zip(history[-32:-1], history[-31:]):
                self.assertEqual(b['queued']-a['queued'], 12-(credits+1))
            self.assertGreaterEqual(history[-1]['queued'], (12-(credits+1))*256)

    def test_08_stable_final_queue_does_not_imply_zero_delay_at_critical_load(self):
        run = simulate(F(8), F(8))
        self.assertGreater(run['history'][-1]['queued'], 0)
        self.assertEqual(run['history'][-1]['queued'], run['history'][-2]['queued'])
        self.assertEqual(run['history'][-1]['slot_delivery'], 8)

    def test_09_fixed_original_source_cut_is_not_the_aggregate_interface(self):
        run = simulate(F(12), F(12))
        state = run['state']
        original = min(w for w in state if w[0] == 0)
        self.assertTrue(all(bit == 0 for bit in original))
        actual_cut = sum(v != original for v, port in state[original].values())
        self.assertEqual(actual_cut, 3)
        self.assertEqual(len(run['queues']), 16)
        self.assertLess(actual_cut, 12)


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
    target = Path(__file__).with_name('local_queue_growth_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
