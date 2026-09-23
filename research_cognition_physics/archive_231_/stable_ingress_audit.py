"""Round 277: a fixed logical inbox and a separate growing tip.

The same unlabelled path history can support different delivery laws.
This is an explicit alternative attachment/role rule, not forced geometry.
"""
import itertools
import unittest
import numpy as np
from vertex_split_cycle_audit import EDGE, replace_vertex
from fisher_hanoi_scaling_audit import bfs
from moving_receiver_audit import growth_at
from growing_stream_audit import simulate, main
from fixed_endpoint_routing_audit import swap_qubits


def anchor_grow(graph, following, owner, frontier):
    predecessor = graph[frontier][0]
    new = len(graph)
    graph = replace_vertex(graph,frontier,EDGE,{predecessor:0},degree_cap=2)
    # Inbox stays on the original owner. Only the construction-tip role moves.
    following[frontier] = new
    return graph, new, (frontier,new)


def path_order(graph):
    order, previous, node = [], None, 0
    while True:
        order.append(node)
        following = [v for v in graph[node] if v != previous]
        if not following:
            return order
        assert len(following) == 1
        previous, node = node, following[0]


def quantum_pipeline():
    run = simulate(2,[1,1,1],[1,1,0],anchor_grow,True)
    # External references 0,1,2; two admitted source inputs 3,4;
    # two distinct relay buffers 5,6; old owner data 7; output archive 8,9.
    # Three newly constructed blank tip registers 10,11,12 remain in |0>.
    width = 13
    psi = np.zeros(1 << width,complex)
    for a,b,c in itertools.product((0,1),repeat=3):
        psi[(a<<0)|(b<<1)|(c<<2)|(a<<3)|(b<<4)|(c<<7)] = 1/np.sqrt(8)
    for row in run['events']:
        touched = set()
        for task,u,v in row['moves']:
            register = lambda node: {0:task+2,1:task+4,2:task+7}[node]
            left,right = register(u),register(v)
            assert left not in touched and right not in touched
            touched.update((left,right))
            psi = swap_qubits(psi,left,right)
    ideal = np.zeros_like(psi)
    for a,b,c in itertools.product((0,1),repeat=3):
        ideal[(a<<0)|(b<<1)|(c<<2)|(a<<8)|(b<<9)|(c<<7)] = 1/np.sqrt(8)
    return {'joint_reference_fidelity':float(abs(np.vdot(ideal,psi))**2),
            'independent_input_qubits':2,'preserved_old_receiver_qubits':1,
            'relay_payload_buffers':2,'output_archive_qubits':2,
            'new_blank_tip_registers':3,'payload_edge_uses':4,
            'control_bit_edge_uses':3,'last_delivery_slot':3,
            'external_reference_qubits':3,'simulated_registers':width}


def report():
    cases = []
    for p in ('none','every_slot','pause_every_2','pause_every_4','pause_at_squares'):
        g = [growth_at(t,p) for t in range(1,257)]
        retreat = simulate(8,g,[1]*256)
        anchor = simulate(8,g,[1]*256,anchor_grow)
        cases.append({'policy':p,'nodes':len(anchor['graph']),
                      'growth':sum(g),'retreat_delivered':len(retreat['deliveries']),
                      'anchor_delivered':len(anchor['deliveries']),
                      'retreat_pending':len(retreat['pending']),
                      'anchor_pending':len(anchor['pending']),
                      'retreat_owner_distance':bfs(retreat['graph'],0)[8],
                      'anchor_owner_distance':bfs(anchor['graph'],0)[8],
                      'anchor_frontier_distance':bfs(anchor['graph'],0)[anchor['frontier']],
                      'anchor_payload_edge_uses':anchor['history'][-1]['payload_edge_uses'],
                      'control_bit_edge_uses':sum(g),
                      'anchor_output_archive_qubits':len(anchor['deliveries'])})
    return {'round':277,
            'scope':'Constructive fixed-ingress protocol with explicit role split, channels, control and memory; not autonomous dimension selection.',
            'comparisons':cases,'quantum_check':quantum_pipeline(),
            'enumeration':{'growth_bits':6,'arrival_bits':6,
                           'distances':[1,2,4],'histories':12288},
            'formula':'B(T)=A(max(0,T-d+1)); tau_s=s+d-1; Q(T)<=d-1',
            'limitation':'Delivery is to the retained logical inbox, not every new cell or the growing frontier.'}


class IngressTests(unittest.TestCase):
    def test_01_all_short_growth_and_input_histories(self):
        for g in itertools.product((0,1),repeat=6):
            for a in itertools.product((0,1),repeat=6):
                for d in (1,2,4):
                    run = simulate(d,g,a,anchor_grow)
                    self.assertEqual(run['deliveries'],
                                     {s:s+d-1 for s in range(1,7) if a[s-1] and s+d-1<=6})
                    self.assertLessEqual(len(run['pending']),d-1)

    def test_02_same_unlabelled_graph_different_marked_owner(self):
        for g in itertools.product((0,1),repeat=6):
            a = [1]*6
            x,y = simulate(3,g,a),simulate(3,g,a,anchor_grow)
            px,py = path_order(x['graph']),path_order(y['graph'])
            self.assertEqual(len(px),len(py))
            self.assertEqual(len(px),len(x['graph']))
            mapping = dict(zip(px,py))
            for u,neighbors in enumerate(x['graph']):
                self.assertEqual({mapping[v] for v in neighbors},set(y['graph'][mapping[u]]))
            self.assertEqual(px.index(3),3+sum(g))
            self.assertEqual(py.index(3),3)

    def test_03_quantum_identity_including_external_references(self):
        q = quantum_pipeline()
        self.assertAlmostEqual(q['joint_reference_fidelity'],1.0,places=12)

    def test_04_actual_edge_uses_and_local_buffers(self):
        run = simulate(4,[1]*32,[1]*32,anchor_grow,True)
        self.assertEqual(sum(len(r['moves']) for r in run['events']),122)
        for row in run['history']:
            self.assertLessEqual(row['transient_payload_buffers'],2)
        self.assertEqual(sum(r['control'] is not None for r in run['events']),32)

    def test_05_delivery_is_not_broadcast_to_growth_tip(self):
        run = simulate(3,[1]*24,[1]*24,anchor_grow)
        self.assertEqual(bfs(run['graph'],run['owner'])[run['frontier']],24)
        self.assertEqual(len(run['deliveries']),22)
        # No payload reaches the tip; this is not a test of collective assimilation.
        self.assertNotIn(run['frontier'],run['pending'])

    def test_06_storage_archive_not_erased_from_budget(self):
        run = simulate(8,[1]*256,[1]*256,anchor_grow)
        self.assertEqual(len(run['deliveries'])+len(run['pending']),256)
        self.assertEqual((len(run['deliveries']),len(run['pending'])),(249,7))
        self.assertEqual(len(run['graph']),265)

    def test_07_growth_pause_or_burst_preserves_deadline(self):
        for g in ([0]*20+[1]*20,[1]*20+[0]*20,[1,0]*20):
            run = simulate(5,g,[1]*40,anchor_grow)
            self.assertEqual(set(time-task+1 for task,time in run['deliveries'].items()),{5})

    def test_08_role_handoff_follows_local_new_edge(self):
        run = simulate(2,[1,0,1,1],[0]*4,anchor_grow,True)
        self.assertEqual([r['control'] for r in run['events']],[(2,3),None,(3,4),(4,5)])
        self.assertEqual(run['owner'],2)
        self.assertEqual(run['frontier'],5)

    def test_09_known_extreme_same_size_comparison(self):
        row = next(r for r in report()['comparisons'] if r['policy']=='every_slot')
        self.assertEqual((row['retreat_delivered'],row['anchor_delivered']),(0,249))
        self.assertEqual((row['retreat_owner_distance'],row['anchor_owner_distance']),(264,8))


if __name__ == '__main__':
    main(__name__,'stable_ingress_audit',report)
