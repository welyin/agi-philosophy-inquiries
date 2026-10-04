"""Round 231: fixed-context influence, possible influence and circuit reachability.

Classical bit routes are diagonal quantum measure/prepare protocols, not a
classical derivation of quantum mechanics. A finite acyclic wiring is an
explicit model input. Matrix examples use a two-dimensional subspace.
"""
import argparse
import json
from pathlib import Path
import platform
import unittest
from collections import deque

import numpy as np


def dagger(x):
    return x.conj().T


def ket(bit):
    return np.eye(2,dtype=complex)[:,bit]


def projector(bit):
    v=ket(bit)
    return np.outer(v,v.conj())


def reset_instrument(output):
    return [np.outer(ket(output),ket(b).conj()) for b in (0,1)]


def chain_distribution(x,y):
    """A prepares x; B records Z then prepares y; C reads Z."""
    rho=projector(x)
    joint=np.zeros((2,2))
    for b,k in enumerate(reset_instrument(y)):
        branch=k@rho@dagger(k)
        for c in (0,1):
            joint[b,c]=np.trace(projector(c)@branch).real
    return joint


def tv(a,b):
    return float(np.sum(np.abs(a-b))/2)


def reachability(adj):
    closure=adj.copy()
    for k in range(len(adj)):
        closure |= closure[:,k,None] & closure[None,k,:]
    return closure


def hops(adj,source,target):
    queue=deque([(source,0)]); visited={source}
    while queue:
        vertex,distance=queue.popleft()
        if vertex==target: return distance
        for nxt in np.flatnonzero(adj[vertex]):
            if int(nxt) not in visited:
                visited.add(int(nxt)); queue.append((int(nxt),distance+1))
    return None


def bit_circuit(adj,source,value):
    """An independent forward evaluation; indices supply a topological ordering.

    Each node reads its incoming orthogonal bit states, computes OR, and
    prepares that orthogonal state on outgoing ports. No unknown state cloning.
    """
    output=np.zeros(len(adj),dtype=int)
    for node in range(len(adj)):
        incoming=np.flatnonzero(adj[:,node])
        output[node]=value if node==source else int(np.any(output[incoming]))
    return output


def enumerate_dag_certificate(n=5):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    reachable=unreachable=0
    for mask in range(1<<len(edges)):
        adj=np.zeros((n,n),dtype=bool)
        for bit,(i,j) in enumerate(edges): adj[i,j]=bool(mask&(1<<bit))
        closure=reachability(adj)
        for source in range(n):
            zero,one=bit_circuit(adj,source,0),bit_circuit(adj,source,1)
            for target in range(n):
                if source==target: continue
                observed=zero[target]!=one[target]
                assert observed==closure[source,target]
                reachable+=int(observed); unreachable+=int(not observed)
    return {'nodes':n,'labelled_dags_in_fixed_topological_order':1<<len(edges),
            'ordered_distinct_pairs_checked':reachable+unreachable,
            'reachable_pairs_with_bit_witness':reachable,
            'unreachable_pairs_without_bit_witness':unreachable}


def trace_first(rho):
    t=rho.reshape(2,2,2,2)
    return np.einsum('abad->bd',t)


def common_cause_certificate():
    bell=np.array([1,0,0,1],dtype=complex)/np.sqrt(2)
    rho=np.outer(bell,bell.conj())
    z=np.diag([1.,-1.]); h=np.array([[1,1],[1,-1]])/np.sqrt(2)
    alice_channels=[[np.eye(2)],[h],reset_instrument(0),reset_instrument(1)]
    bob=[]
    for channel in alice_channels:
        output=sum(np.kron(k,np.eye(2))@rho@dagger(np.kron(k,np.eye(2))) for k in channel)
        bob.append(trace_first(output))
    errors=[float(np.linalg.norm(x-np.eye(2)/2)) for x in bob]
    return {'Z_correlation':float(np.trace(np.kron(z,z)@rho).real),
            'maximum_B_marginal_change_after_A_TP_action':max(errors)}


def ambiguity_certificate():
    chain=np.array([[0,1,0],[0,0,1],[0,0,0]],dtype=bool)
    shortcut=chain.copy(); shortcut[0,2]=True
    return {'same_reachability':bool(np.array_equal(reachability(chain),reachability(shortcut))),
            'chain_A_to_C_hops':hops(chain,0,2),'shortcut_A_to_C_hops':hops(shortcut,0,2),
            'order_compatible_clock_1':[0,1,2],'order_compatible_clock_2':[0,1,4],
            'clock_1_successive_intervals':[1,1],'clock_2_successive_intervals':[1,3]}


class EventOrderTests(unittest.TestCase):
    def test_measure_reset_is_a_complete_quantum_instrument(self):
        for y in (0,1):
            np.testing.assert_allclose(sum(dagger(k)@k for k in reset_instrument(y)),np.eye(2))

    def test_fixed_context_influence_is_not_transitive(self):
        p00,p10,p01=chain_distribution(0,0),chain_distribution(1,0),chain_distribution(0,1)
        self.assertEqual(tv(p00.sum(axis=1),p10.sum(axis=1)),1.)
        self.assertEqual(tv(p00.sum(axis=0),p01.sum(axis=0)),1.)
        self.assertEqual(tv(p00.sum(axis=0),p10.sum(axis=0)),0.)
        # Changing the later preparation cannot change B's already-read result.
        np.testing.assert_allclose(p00.sum(axis=1),p01.sum(axis=1))

    def test_change_of_relay_context_restores_A_to_C(self):
        outputs=[projector(x) for x in (0,1)]  # identity at B
        self.assertEqual(tv(np.diag(outputs[0]),np.diag(outputs[1])),1.)
        # Read and reprepare the observed bit also transmits this orthogonal code.
        for x in (0,1):
            output=sum(projector(b)@projector(x)@projector(b) for b in (0,1))
            np.testing.assert_allclose(output,projector(x))

    def test_all_five_vertex_dags_have_the_claimed_bit_witnesses(self):
        c=enumerate_dag_certificate()
        self.assertEqual(c['ordered_distinct_pairs_checked'],20480)

    def test_correlation_is_not_interventional_signalling(self):
        c=common_cause_certificate()
        self.assertAlmostEqual(c['Z_correlation'],1.)
        self.assertLess(c['maximum_B_marginal_change_after_A_TP_action'],1e-14)

    def test_reachability_does_not_fix_primitive_edges_or_hop_distance(self):
        c=ambiguity_certificate()
        self.assertTrue(c['same_reachability'])
        self.assertEqual((c['chain_A_to_C_hops'],c['shortcut_A_to_C_hops']),(2,1))

    def test_order_does_not_fix_elapsed_time(self):
        t=np.array([0,1,2]); s=t*t
        np.testing.assert_array_equal(t[:,None]<t[None,:],s[:,None]<s[None,:])
        self.assertFalse(np.allclose(np.diff(s),np.diff(t)))

    def test_small_relay_errors_preserve_a_nonzero_signal(self):
        x=np.array([[0.,1.],[1.,0.]])
        for length in (1,2,4,8):
            p=.02
            outputs=[projector(b) for b in (0,1)]
            for _ in range(length): outputs=[(1-p)*r+p*x@r@x for r in outputs]
            gap=tv(np.diag(outputs[0]),np.diag(outputs[1]))
            self.assertAlmostEqual(gap,(1-2*p)**length)
            self.assertGreaterEqual(gap,1-2*length*p-1e-13)


def report():
    p00,p10,p01=chain_distribution(0,0),chain_distribution(1,0),chain_distribution(0,1)
    return {'round':231,'date':'2026-09-21',
            'hypothesis':'Do operational influences themselves provide the transitive relation needed for event order?',
            'analytic_results':['Fixed-context interventional influence need not be transitive.',
                                'In a finite acyclic quantum circuit with freely selectable instruments and nontrivial noiseless wires, possible influence over relay contexts is exactly directed reachability.',
                                'Round228 density is sufficient for a nonzero routing signal; the continuous seed is not required for this weak claim.',
                                'Reachability alone fixes neither primitive adjacency, propagation duration nor a physical metric.'],
            'extra_model_inputs':['finite acyclic typed circuit wiring','independent intervention choices','wire systems carrying two distinguishable states','freely chosen relay contexts; no postselection in signalling tests'],
            'fixed_context_TV':{'A_to_B':tv(p00.sum(axis=1),p10.sum(axis=1)),
                                'B_to_C':tv(p00.sum(axis=0),p01.sum(axis=0)),
                                'A_to_C':tv(p00.sum(axis=0),p10.sum(axis=0))},
            'DAG_certificate':enumerate_dag_certificate(),
            'common_cause':common_cause_certificate(),'order_ambiguity':ambiguity_certificate(),
            'not_derived':['the universe has a fixed acyclic event graph','metric time and spatial dimension','universal limiting speed or Lorentz symmetry','mass, gauge group or gravity'],
            'numerical_scope':'Finite diagonal quantum routing witnesses and small matrix counterexamples; graph input and general analytic proof are separate.',
            'runtime':{'python':platform.python_version(),'numpy':np.__version__},
            'next':'Audit which composition or record-consistency conditions can justify an event order without assuming that order as a hidden primitive.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(EventOrderTests))
    if not result.wasSuccessful(): raise SystemExit(1)
    data=report(); data['checks']={'run':result.testsRun,'failures':0,'errors':0}
    if args.write_results:
        target=Path(__file__).with_name('event_order_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8'))!=data:
            raise RuntimeError('Existing research result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
