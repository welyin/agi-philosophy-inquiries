"""Round233: event-aware grouping of ordinary circuits and memory channels.

Finite wiring and event labels are explicit model inputs. This is a typing
audit, not a derivation of the existence of time or of a physical clock.
"""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
import json
from pathlib import Path
import platform
import unittest

import numpy as np


def partitions(n):
    """Canonical restricted-growth labels: every set partition exactly once."""
    def grow(prefix):
        if len(prefix)==n:
            yield tuple(prefix);return
        for label in range(max(prefix,default=-1)+2):
            yield from grow(prefix+[label])
    yield from grow([])


def quotient(edges,labels):
    return {(labels[a],labels[b]) for a,b in edges if labels[a]!=labels[b]}


def acyclic(n,edges):
    incoming=[0]*n;out=[[] for _ in range(n)]
    for a,b in edges: incoming[b]+=1;out[a].append(b)
    ready=[v for v in range(n) if incoming[v]==0];seen=0
    while ready:
        v=ready.pop();seen+=1
        for nxt in out[v]:
            incoming[nxt]-=1
            if incoming[nxt]==0:ready.append(nxt)
    return seen==n


@lru_cache(None)
def grouping_certificate(n=5):
    possible=[(a,b) for a in range(n) for b in range(a+1,n)]
    orders=list(permutations(range(n)));parts=list(partitions(n))
    # Independent certificate: a topological ordering in which each block
    # occupies one contiguous interval. No quotient algorithm is used here.
    order_masks=[]
    for order in orders:
        pos={v:i for i,v in enumerate(order)}
        order_masks.append(sum(1<<j for j,(a,b) in enumerate(possible) if pos[a]<pos[b]))
    contiguous=[]
    for labels in parts:
        bitset=0
        for j,order in enumerate(orders):
            sequence=[labels[v] for v in order]
            runs=[sequence[0]]+[sequence[i] for i in range(1,n) if sequence[i]!=sequence[i-1]]
            if len(runs)==len(set(runs)):bitset|=1<<j
        contiguous.append(bitset)
    good=bad=0
    for mask in range(1<<len(possible)):
        edges=[e for j,e in enumerate(possible) if mask&(1<<j)]
        topo=sum(1<<j for j,allowed in enumerate(order_masks) if mask&~allowed==0)
        for labels,contig in zip(parts,contiguous):
            by_quotient=acyclic(max(labels)+1,quotient(edges,labels))
            by_schedule=bool(topo&contig)
            assert by_quotient==by_schedule
            good+=int(by_quotient);bad+=int(not by_quotient)
    return {'nodes':n,'fixed_topological_order_dags':1<<len(possible),
            'partitions_per_graph':len(parts),'graph_partition_pairs':good+bad,
            'acyclic_quotients':good,'cyclic_quotients':bad,
            'independent_contiguous_schedule_disagreements':0}


def two_channels_fixed_points(fa,fb,opposite=True):
    fixed=[]
    for a,b in product(range(4),repeat=2):
        oa,ob=fa[a],fb[b]
        # a=x_A+2y_A, b=x_B+2y_B; output indices use the same convention.
        expected=(2*(ob//2),oa%2) if opposite else (0,oa)
        if expected==(a,b):fixed.append((a,b))
    return fixed


@lru_cache(None)
def channel_grouping_certificate():
    controls=list(permutations(range(4)))
    hist={}
    for opposite,name in ((False,'same_order'),(True,'opposite_orders')):
        counts=Counter(len(two_channels_fixed_points(a,b,opposite)) for a,b in product(controls,repeat=2))
        hist[name]={str(k):counts[k] for k in sorted(counts)}
    fa=(1,3,0,2);fb=(0,2,1,3)
    return {'all_two_bit_permutation_pairs':len(controls)**2,'fixed_point_histograms':hist,
            'bad_regrouping_witness':{'A':fa,'B':fb,'fixed_points':two_channels_fixed_points(fa,fb)},
            'wiring_edges':[(0,1),(2,3)],'partition_labels':[0,1,1,0],
            'original_event_graph_acyclic':True,'quotient_edges':[(0,1),(1,0)]}


def partial_trace(matrix,dims,keep):
    keep=tuple(keep);remaining=list(range(len(dims)))
    tensor=matrix.reshape(tuple(dims)*2)
    for subsystem in reversed(range(len(dims))):
        if subsystem not in keep:
            axis=remaining.index(subsystem)
            tensor=np.trace(tensor,axis1=axis,axis2=axis+len(remaining))
            remaining.remove(subsystem)
    size=int(np.prod([dims[k] for k in keep],dtype=int))
    return tensor.reshape(size,size)


def random_isometry(rows,cols,rng):
    z=rng.normal(size=(rows,cols))+1j*rng.normal(size=(rows,cols))
    return np.linalg.qr(z)[0][:,:cols]


def memory_choi(v1,v2,memory1=2,memory2=2):
    """V1: A1 -> B1 M1; V2: A2 M1 -> B2 M2.

    Return R2 ordered (B2,A2,B1,A1), R1 ordered (B1,A1), and full
    Kraus operators ordered outputs (B1,B2), inputs (A1,A2).
    """
    t1=v1.reshape(2,memory1,2)
    t2=v2.reshape(2,memory2,2,memory1)
    kraus=[]
    for m in range(memory2):
        k=np.einsum('bcx,dzc->bdxz',t1,t2[:,m,:,:])
        kraus.append(k.reshape(4,4))
    r2=np.zeros((16,16),complex)
    for k in kraus:
        vec=k.reshape(2,2,2,2).transpose(1,3,0,2).reshape(16)
        r2+=np.outer(vec,vec.conj())
    r1=sum(np.outer(t1[:,m,:].reshape(4),t1[:,m,:].reshape(4).conj()) for m in range(memory1))
    return r2,r1,kraus


def two_slot_residual(r2):
    reduced=partial_trace(r2,[2]*4,[1,2,3])  # remove B2
    r1=partial_trace(reduced,[2]*3,[1,2])/2  # remove A2
    return float(np.linalg.norm(reduced-np.kron(np.eye(2),r1))),r1


def unitary_choi(unitary):
    vec=unitary.reshape(2,2,2,2).transpose(1,3,0,2).reshape(16)
    return np.outer(vec,vec.conj())


@lru_cache(None)
def comb_certificate():
    rng=np.random.default_rng(233);norm_error=comb_error=0.;min_eigen=0.
    for _ in range(48):
        r2,r1,kraus=memory_choi(random_isometry(4,2,rng),random_isometry(4,4,rng))
        norm_error=max(norm_error,float(np.linalg.norm(sum(k.conj().T@k for k in kraus)-np.eye(4))))
        comb_error=max(comb_error,float(np.linalg.norm(partial_trace(r2,[2]*4,[1,2,3])-np.kron(np.eye(2),r1))),
                       float(np.linalg.norm(partial_trace(r1,[2,2],[1])-np.eye(2))))
        min_eigen=min(min_eigen,float(np.linalg.eigvalsh(r2).min()))
    swap=np.eye(4)[[0,2,1,3]]
    r_bad=unitary_choi(swap);bad_residual,_=two_slot_residual(r_bad)
    return {'random_two_slot_memory_channels':48,'maximum_TP_residual':norm_error,
            'maximum_comb_normalization_residual':comb_error,
            'minimum_choi_eigenvalue':min_eigen,
            'SWAP_channel':{'TP_residual':float(np.linalg.norm(swap.conj().T@swap-np.eye(4))),
                            'comb_residual':bad_residual,'later_input_to_earlier_output_TV':1.0}}


class TypedCompositionTests(unittest.TestCase):
    def test_all_five_event_graph_partitions_agree_with_contiguous_schedules(self):
        c=grouping_certificate()
        self.assertEqual(c['graph_partition_pairs'],53248)
        self.assertEqual(c['independent_contiguous_schedule_disagreements'],0)

    def test_separately_antichain_blocks_can_still_make_a_cycle(self):
        edges={(0,1),(2,3)};labels=(0,1,1,0)
        self.assertTrue(acyclic(4,edges))
        self.assertTrue(all(labels[a]!=labels[b] for a,b in edges))
        self.assertFalse(acyclic(2,quotient(edges,labels)))

    def test_same_order_parallel_channels_allow_all_joint_permutations(self):
        self.assertEqual(channel_grouping_certificate()['fixed_point_histograms']['same_order'],{'1':576})

    def test_opposite_order_false_grouping_has_explicit_normalization_failure(self):
        c=channel_grouping_certificate()
        self.assertEqual(c['bad_regrouping_witness']['fixed_points'],[])
        self.assertEqual(c['fixed_point_histograms']['opposite_orders'],{'0':128,'1':320,'2':128})

    def test_reversible_routing_on_a_cycle_can_force_zero_normalization(self):
        for length in range(2,7):
            fixed=0
            for inputs in product(range(2),repeat=length):
                outputs=(1-inputs[0],)+inputs[1:]
                fixed+=all(inputs[(j+1)%length]==outputs[j] for j in range(length))
            self.assertEqual(fixed,0)

    def test_memory_dilations_satisfy_recursive_trace_conditions(self):
        c=comb_certificate()
        self.assertLess(c['maximum_TP_residual'],8e-15)
        self.assertLess(c['maximum_comb_normalization_residual'],8e-15)
        self.assertGreater(c['minimum_choi_eigenvalue'],-8e-15)

    def test_a_channel_can_be_CPTP_but_violate_its_two_slot_order(self):
        c=comb_certificate()['SWAP_channel']
        self.assertEqual(c['TP_residual'],0.)
        self.assertGreater(c['comb_residual'],1.)
        swap=np.eye(4)[[0,2,1,3]]
        early=[]
        for late_bit in (0,1):
            psi=swap@np.eye(4)[:,late_bit]  # first input fixed to zero
            early.append(partial_trace(np.outer(psi,psi),[2,2],[0]))
        self.assertEqual(float(np.linalg.norm(early[0]-early[1],ord='nuc')/2),1.)

    def test_a_valid_memory_channel_can_delay_an_entangled_unknown_input(self):
        v1=np.zeros((4,2));v1[0,0]=v1[1,1]=1.  # B1=0, M1=A1
        swap=np.eye(4)[[0,2,1,3]]                # B2=M1, M2=A2
        r2,r1,kraus=memory_choi(v1,swap)
        self.assertLess(two_slot_residual(r2)[0],1e-14)
        # Bell A1-reference, A2 fixed; preserve its entanglement in B2-reference.
        psi=np.zeros(8);psi[0]=psi[5]=1/np.sqrt(2)  # order A1,A2,R
        rho=np.outer(psi,psi)
        out=sum(np.kron(k,np.eye(2))@rho@np.kron(k.conj().T,np.eye(2)) for k in kraus)
        bell=np.array([1,0,0,1])/np.sqrt(2)
        np.testing.assert_allclose(partial_trace(out,[2,2,2],[1,2]),np.outer(bell,bell),atol=1e-14)
        np.testing.assert_allclose(partial_trace(out,[2,2,2],[0]),np.diag([1.,0.]),atol=1e-14)


def report():
    return {'round':233,'date':'2026-09-22',
            'hypothesis':'When does replacing event blocks by arbitrary joint channels preserve universal normalization?',
            'graph_grouping':grouping_certificate(),'ordinary_channel_controls':channel_grouping_certificate(),
            'quantum_memory_controls':comb_certificate(),
            'analytic_result':'For finite noiseless typed wiring with arbitrary block instruments, quotient acyclicity is equivalent to block-contiguous scheduling and universal normalization.',
            'interpretation':'A persistent agent may require multiple ordered operation slots; coarse lab labels do not grant arbitrary cross-slot channels.',
            'extra_inputs':['Given finite event graph and selected grouping','Two distinguishable states on relevant wires','Unrestricted or dense local replacement permissions','No unlisted communication or postselection'],
            'not_derived':['Existence of a fundamental event graph','Metric time','Spacetime dimension or limiting speed','General characterization of indefinite-order processes'],
            'next':'Check the literal identical-copy condition against a causally controlled order, retaining schedule records and event slots.',
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(TypedCompositionTests))
    if not result.wasSuccessful():raise SystemExit(1)
    data=report();data['checks']={'run':result.testsRun,'failures':0,'errors':0}
    if args.write_results:
        target=Path(__file__).with_name('typed_composition_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8'))!=json.loads(json.dumps(data)):
            raise RuntimeError('Existing research result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
