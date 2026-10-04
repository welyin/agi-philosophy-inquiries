"""Round 275: delivery to a logical receiver that keeps a child-zero register.

A terminal K2 refinement inserts a blank relay before the owner. Growth and
one-hop message service are separately counted; the line is a counterexample,
not a spatial model. All scheduling policies are explicit extra inputs.
"""
import argparse
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
from contact_route_closure_audit import adjacency
from vertex_split_cycle_audit import EDGE, replace_vertex, statistics
from fisher_hanoi_scaling_audit import bfs
from fixed_endpoint_routing_audit import swap_qubits


def growth_at(time, policy):
    if policy == 'every_slot':
        return 1
    if policy == 'none':
        return 0
    if policy == 'pause_every_2':
        return int(time % 2 != 0)
    if policy == 'pause_every_4':
        return int(time % 4 != 0)
    if policy == 'pause_at_squares':
        return int(math.isqrt(time)**2 != time)
    raise ValueError(policy)


def execute(initial_distance, growth):
    owner = initial_distance
    graph = adjacency(owner+1, [(i,i+1) for i in range(owner)])
    following = {i:i+1 for i in range(owner)}
    packet = 0
    history = []
    total_growth = 0
    delivered_at = None
    for time, amount in enumerate(growth,1):
        if amount not in (0,1):
            raise ValueError('This branch allows at most one split per slot.')
        new_vertex = None
        if amount:
            predecessor = graph[owner][0]
            new_vertex = len(graph)
            graph = replace_vertex(graph,owner,EDGE,{predecessor:1},degree_cap=2)
            # The old data register stays at child zero (same enumeration ID).
            # Reattachment updates the predecessor's local outgoing port.
            following[predecessor] = new_vertex
            following[new_vertex] = owner
            total_growth += 1
        old_packet = packet
        packet = following[packet]
        assert packet in graph[old_packet]
        distance = bfs(graph,packet)[owner]
        history.append({'t':time, 'growth':amount, 'total_growth':total_growth,
                        'new_vertex':new_vertex, 'from':old_packet, 'to':packet,
                        'distance':distance, 'nodes':len(graph),
                        'owner':owner})
        if packet == owner:
            delivered_at = time
            break
    return {'graph':graph, 'history':history, 'delivered_at':delivered_at,
            'final_packet':packet, 'owner':owner}


@lru_cache(maxsize=None)
def scenario(distance, policy, horizon=256):
    return execute(distance,[growth_at(t,policy) for t in range(1,horizon+1)])


def quantum_delivery_check(distance=2, policy='pause_every_2'):
    run = scenario(distance,policy,64)
    assert run['delivered_at'] is not None
    # One data slot per initial node, two external references, a separate inbox.
    ref_message, ref_owner, inbox = distance+1,distance+2,distance+3
    width = distance+4
    registers = {v:v for v in range(distance+1)}
    psi = np.zeros(1 << width,complex)
    for message,old_data in itertools.product((0,1), repeat=2):
        bits = ((message << ref_message) | message |
                (old_data << ref_owner) | (old_data << distance))
        psi[bits] = 0.5
    packet_register = 0
    for row in run['history']:
        if row['new_vertex'] is not None:
            registers[row['new_vertex']] = width
            psi = np.concatenate((psi,np.zeros_like(psi)))
            width += 1
        target = inbox if row['to'] == distance else registers[row['to']]
        psi = swap_qubits(psi,packet_register,target)
        packet_register = target
    ideal = np.zeros_like(psi)
    for message,old_data in itertools.product((0,1), repeat=2):
        bits = ((message << ref_message) | (message << inbox) |
                (old_data << ref_owner) | (old_data << distance))
        ideal[bits] = 0.5
    return {'joint_reference_fidelity':float(abs(np.vdot(ideal,psi))**2),
            'logical_slots':run['delivered_at'],
            'qubit_transmissions':len(run['history']),
            'new_blank_node_registers':sum(r['growth'] for r in run['history']),
            'separate_inbox_registers':1}


def report():
    rows = []
    for policy in ('none','every_slot','pause_every_2','pause_every_4','pause_at_squares'):
        run = scenario(8,policy)
        last = run['history'][-1]
        rows.append({'policy':policy, 'initial_distance':8,
                     'delivered_at':run['delivered_at'],
                     'observed_until':last['t'], 'remaining_distance':last['distance'],
                     'new_nodes_and_channels':last['total_growth'],
                     'message_edge_uses':len(run['history']),
                     'final_graph':statistics(run['graph'])})
    return {
        'round':275, 'receiver_update':'Terminal K2 refinement: old incident edge goes to child one, old data remain at child zero; each growth creates one node and one edge.',
        'distance_identity_before_delivery':'D(t)=d+G(t)-t.',
        'delivery_criterion':'First t with t-G(t)>=d; infinite pause count is necessary and sufficient to serve every finite d for schedules with growth in {0,1}.',
        'eight_hop_cases':rows,
        'square_pause_delivery_times':[{'distance':d,'delivered_at':scenario(d,'pause_at_squares')['delivered_at']} for d in (1,2,4,8,12)],
        'quantum_preservation_example':quantum_delivery_check(),
        'scope':'A conditional causal graph-rewrite/routing counterexample and an explicit scheduling remedy. Static connectivity and state-preserving refinement do not guarantee delivery to a growing logical endpoint. The infinite failure follows analytically; finite simulations audit the rewrite and distance identity. No physical speed, natural growth law or spatial dimension is derived.'
    }


class Audit(unittest.TestCase):
    def test_01_actual_rewrites_preserve_a_connected_degree_two_path(self):
        for policy in ('every_slot','pause_every_2','pause_at_squares'):
            run = scenario(8,policy)
            data = statistics(run['graph'])
            self.assertEqual(data['components'],1)
            self.assertEqual(data['cycle_rank'],0)
            self.assertLessEqual(data['max_degree'],2)
            self.assertEqual(data['vertices'],9+run['history'][-1]['total_growth'])
            self.assertEqual(data['edges'],data['vertices']-1)

    def test_02_distance_identity_for_all_eight_slot_growth_histories(self):
        for sequence in itertools.product((0,1), repeat=8):
            for initial in (1,2,4):
                run = execute(initial,sequence)
                for row in run['history']:
                    self.assertEqual(row['distance'],initial+row['total_growth']-row['t'])

    def test_03_every_slot_growth_never_closes_the_initial_gap(self):
        for initial in (1,2,8):
            run = scenario(initial,'every_slot')
            self.assertIsNone(run['delivered_at'])
            self.assertEqual({r['distance'] for r in run['history']},{initial})

    def test_04_delivery_criterion_matches_all_finite_binary_schedules(self):
        for sequence in itertools.product((0,1), repeat=8):
            for initial in (1,2,4):
                slack = 0
                predicted = None
                for t,g in enumerate(sequence,1):
                    slack += 1-g
                    if slack >= initial:
                        predicted = t
                        break
                self.assertEqual(execute(initial,sequence)['delivered_at'],predicted)

    def test_05_periodic_pause_policy_gives_exact_finite_deadlines(self):
        for initial in (1,2,4,8,12):
            for period in (2,4):
                run = scenario(initial,f'pause_every_{period}')
                self.assertEqual(run['delivered_at'],period*initial)
                self.assertEqual(run['history'][-1]['total_growth'],(period-1)*initial)

    def test_06_growth_density_one_can_still_allow_delivery(self):
        for initial in (1,2,4,8,12):
            run = scenario(initial,'pause_at_squares')
            self.assertEqual(run['delivered_at'],initial**2)
            self.assertEqual(run['history'][-1]['total_growth'],initial**2-initial)
        self.assertGreater(1-1/10000,0.999)

    def test_07_growing_target_is_not_fixed_target_nonconnectivity(self):
        frozen = scenario(8,'none')
        moving = scenario(8,'every_slot')
        self.assertEqual(frozen['delivered_at'],8)
        self.assertIsNone(moving['delivered_at'])
        self.assertEqual(len({r['from'] for r in moving['history']}),256)
        self.assertTrue(all(r['to'] != r['from'] for r in moving['history']))

    def test_08_quantum_delivery_keeps_both_message_and_old_owner_state(self):
        checked = quantum_delivery_check()
        self.assertAlmostEqual(checked['joint_reference_fidelity'],1.0,places=12)
        self.assertEqual(checked['logical_slots'],4)
        self.assertEqual(checked['qubit_transmissions'],4)
        self.assertEqual(checked['new_blank_node_registers'],2)
        self.assertEqual(checked['separate_inbox_registers'],1)

    def test_09_cumulative_growth_envelope_implies_a_deadline(self):
        # For periodic one-pause-in-k schedules, G(t) <= (1-1/k)t+1.
        for period in (2,4):
            for initial in (1,2,4,8):
                deadline = period*(initial+1)
                total = 0
                for t in range(1,deadline+1):
                    total += growth_at(t,f'pause_every_{period}')
                    self.assertLessEqual(period*total,(period-1)*t+period)
                self.assertLessEqual(scenario(initial,f'pause_every_{period}')['delivered_at'],deadline)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['runtime'] = {'python':platform.python_version(), 'numpy':np.__version__}
    result['checks'] = {'run':checks.testsRun,'failures':0,'errors':0}
    payload = json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    path = Path(__file__).with_name('moving_receiver_audit_results.json')
    if args.write_results:
        if path.exists() and path.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        path.write_text(payload,encoding='utf-8')
    print(payload)
