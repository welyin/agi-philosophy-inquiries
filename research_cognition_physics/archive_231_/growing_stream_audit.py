"""Round 276: exact sample paths for a stream and a retreating receiver.

Slot order: local graph rewrite, input admission, one-hop parallel transmission.
Channels, clock, construction, blank memory and output archives are inputs.
Default execution is read-only.
"""
import argparse
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from contact_route_closure_audit import adjacency
from vertex_split_cycle_audit import EDGE, replace_vertex, statistics
from moving_receiver_audit import growth_at


def retreat_grow(graph, following, owner, frontier):
    predecessor = graph[owner][0]
    new = len(graph)
    graph = replace_vertex(graph, owner, EDGE, {predecessor: 1}, degree_cap=2)
    following[predecessor] = new
    following[new] = owner
    return graph, owner, (owner, new)


def simulate(distance, growth, arrivals, grow_rule=retreat_grow, keep_events=False):
    if distance < 1 or len(growth) != len(arrivals):
        raise ValueError('Positive distance and matching finite schedules required.')
    graph = adjacency(distance+1, [(v,v+1) for v in range(distance)])
    owner = frontier = distance
    following = {v:v+1 for v in range(distance)}
    packets, delivered, history, events = {}, {}, [], []
    admitted = total_growth = hops = 0
    for time, (g,a) in enumerate(zip(growth, arrivals), 1):
        if g not in (0,1) or a not in (0,1):
            raise ValueError('At most one construction and admission per slot.')
        control = None
        if g:
            graph, frontier, control = grow_rule(graph, following, owner, frontier)
            total_growth += 1
            assert control[1] in graph[control[0]]
        if a:
            assert 0 not in packets
            packets[0] = time
            admitted += 1
        after, used, inbound = {}, set(), set()
        before = set(packets)
        moves = []
        for node, task in packets.items():
            target = following[node]
            edge = tuple(sorted((node,target)))
            assert target in graph[node] and edge not in used
            used.add(edge)
            assert target not in inbound
            inbound.add(target)
            moves.append((task,node,target))
            if target == owner:
                assert task not in delivered
                delivered[task] = time
            else:
                assert target not in after
                after[target] = task
        if control:
            # A counted one-bit ready/role handoff uses the newly created edge.
            assert tuple(sorted(control)) not in used
        hops += len(moves)
        packets = after
        assert admitted == len(delivered)+len(packets)
        history.append({'t':time, 'growth':total_growth, 'nodes':len(graph),
                        'admitted':admitted, 'delivered':len(delivered),
                        'pending':len(packets), 'payload_edge_uses':hops,
                        'control_bit_edge_uses':total_growth,
                        'end_slot_max_payload_occupancy':int(bool(packets)),
                        'transient_payload_buffers':max(
                            [int(v in before)+int(v in inbound)
                             for v in before | inbound] or [0]),
                        'blocked_hops':0})
        if keep_events:
            events.append({'t':time, 'moves':moves, 'control':control})
    return {'graph':graph, 'owner':owner, 'frontier':frontier,
            'history':history, 'deliveries':delivered, 'pending':packets,
            'events':events}


def expected_delivery_count(distance, growth, arrivals, time):
    cutoff = max(0, time-sum(growth[:time])-distance+1)
    return sum(arrivals[:cutoff])


def case(policy, distance=8, horizon=256, sparse=False):
    growth = [growth_at(t,policy) for t in range(1,horizon+1)]
    arrivals = [int(not sparse or t % 2 == 1) for t in range(1,horizon+1)]
    run = simulate(distance,growth,arrivals)
    last = run['history'][-1]
    return {'policy':policy, 'distance':distance, 'horizon':horizon,
            'arrival_policy':'odd_slots' if sparse else 'every_slot',
            **last, 'graph':statistics(run['graph']),
            'pending_area':sum(row['pending'] for row in run['history']),
            'output_archive_qubits':last['delivered'],
            'transport_payload_qubits':last['pending'],
            'local_transient_buffer_bound':max(
                row['transient_payload_buffers'] for row in run['history'])}


def report():
    return {'round':276,
            'scope':'Specified expanding path, exact routing and finite resource accounting; no derived spacetime or stationary queue claim.',
            'scenarios':[case(p) for p in ('none','every_slot','pause_every_2',
                                          'pause_every_4','pause_at_squares')]
                         + [case('pause_every_2',sparse=True)],
            'enumeration':{'growth_bits':6,'arrival_bits':6,
                           'distances':[1,2,4],'histories':12288},
            'formulas':{'delivered':'B(T)=A(max(0,T-G(T)-d+1))',
                        'full_stream_pending':'Q(T)=min(T,d+G(T)-1)',
                        'area':'sum_t Q(t)=sum_admitted_s [min(tau_s,T+1)-s]'},
            'added_inputs':['clock and atomic rewrite','one payload or control bit per edge per slot',
                            'up to two transient payload registers per relay',
                            'new nodes/channels and local port handoff',
                            'separate output archive and old receiver register']}


class StreamTests(unittest.TestCase):
    def test_01_all_short_growth_and_input_histories(self):
        for g in itertools.product((0,1),repeat=6):
            for a in itertools.product((0,1),repeat=6):
                for d in (1,2,4):
                    run = simulate(d,g,a)
                    for row in run['history']:
                        t = row['t']
                        self.assertEqual(row['delivered'],expected_delivery_count(d,g,a,t))

    def test_02_exact_finite_area_including_censored_tasks(self):
        for p in ('none','every_slot','pause_every_2','pause_at_squares'):
            g = [growth_at(t,p) for t in range(1,81)]
            a = [int(t % 3 != 0) for t in range(1,81)]
            run = simulate(5,g,a)
            censored = sum(min(run['deliveries'].get(s,81),81)-s
                           for s in range(1,81) if a[s-1])
            self.assertEqual(sum(r['pending'] for r in run['history']),censored)

    def test_03_exact_periodic_task_deadlines(self):
        for k in (2,4):
            run = simulate(3,[int(t % k != 0) for t in range(1,81)],[1]*80)
            for task,time in run['deliveries'].items():
                self.assertEqual(time,k*(3+task-1))

    def test_04_square_pause_eventual_service(self):
        run = simulate(2,[growth_at(t,'pause_at_squares') for t in range(1,101)],[1]*100)
        self.assertEqual(run['deliveries'],{s:(s+1)**2 for s in range(1,10)})

    def test_05_full_stream_node_and_storage_identity(self):
        for p in ('none','every_slot','pause_every_2'):
            run = simulate(8,[growth_at(t,p) for t in range(1,65)],[1]*64)
            for row in run['history']:
                self.assertEqual(row['pending'],min(row['t'],7+row['growth']))
                self.assertEqual(row['nodes'],9+row['growth'])
                self.assertEqual(row['admitted'],row['pending']+row['delivered'])

    def test_06_no_collision_and_two_transient_buffers(self):
        run = simulate(8,[growth_at(t,'pause_every_2') for t in range(1,65)],
                       [1]*64,keep_events=True)
        for row in run['history']:
            self.assertLessEqual(row['transient_payload_buffers'],2)
            self.assertLessEqual(row['end_slot_max_payload_occupancy'],1)
        self.assertEqual(max(r['transient_payload_buffers'] for r in run['history']),2)
        self.assertEqual(sum(len(r['moves']) for r in run['events']),
                         run['history'][-1]['payload_edge_uses'])

    def test_07_never_delivered_stream_is_not_stable(self):
        row = case('every_slot')
        self.assertEqual((row['delivered'],row['pending'],row['blocked_hops']),(0,256,0))
        self.assertEqual(row['control_bit_edge_uses'],256)

    def test_08_known_counts_and_archives(self):
        pairs = {'none':(249,7),'pause_every_2':(121,135),
                 'pause_every_4':(57,199),'pause_at_squares':(9,247)}
        for p,expected in pairs.items():
            row = case(p)
            self.assertEqual((row['delivered'],row['pending']),expected)
            self.assertEqual(row['output_archive_qubits']+row['transport_payload_qubits'],256)
        sparse = case('pause_every_2',sparse=True)
        self.assertEqual((sparse['admitted'],sparse['delivered'],sparse['pending']),(128,61,67))


def main(module_name, result_stem, report_function):
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    tests = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromName(module_name))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    value = report_function()
    value['checks'] = {'run':tests.testsRun,'failures':0,'errors':0}
    value['runtime'] = {'python':platform.python_version(),'numpy':np.__version__}
    if args.write_results:
        path = Path(__file__).with_name(result_stem+'_results.json')
        if path.exists():
            if json.loads(path.read_text(encoding='utf-8')) != value:
                raise RuntimeError('Refusing to replace different saved results.')
        else:
            path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(value,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main(__name__,'growing_stream_audit',report)
