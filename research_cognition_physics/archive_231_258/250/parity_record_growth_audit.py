"""Round 250: redundant parity records leave an encoded qubit coherent.

Specialize round 249 to two roots and k children sharing those roots.
The supplied CZ interaction records joint parity, not either parent's bit.
Conditional geometry probabilities are reported rather than hidden.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import graph_record_growth_audit as graph
import online_record_growth_audit as online
import random_order_growth_audit as growth

J=np.kron(graph.Z,graph.Z)
PROJECTORS=[(np.eye(4)+s*J)/2 for s in (1,-1)]


def record_map(k):
    """Independent sequential gate construction on arbitrary two-root input."""
    v=np.eye(4)
    for j in range(k):
        n=j+3
        v=graph.phase(n,[(0,n-1),(1,n-1)])[:,None]*np.kron(v,graph.PLUS[:,None])
    return v


def product_record(k,s):
    out=np.ones(1)
    for _ in range(k):
        out=np.kron(out,np.array([1.,s])/math.sqrt(2))
    return out


def sector_formula(k):
    return sum(np.kron(p,product_record(k,s)[:,None]) for p,s in zip(PROJECTORS,(1,-1)))


def root_density_after_loss(psi,k):
    a=(record_map(k)@psi).reshape(4,2**k)
    return a@a.conj().T


def child_density(psi,k,child=0):
    a=(record_map(k)@psi).reshape((4,)+(2,)*k)
    a=np.moveaxis(a,1+child,0).reshape(2,-1)
    return a@a.conj().T


def children_density(psi,k):
    a=(record_map(k)@psi).reshape(4,2**k)
    return a.T@a.conj()


def family(k):
    r=np.zeros((k+2,k+2),dtype=bool)
    r[0,2:]=True; r[1,2:]=True
    return growth.key(r)


def family_probability(k):
    out=F(1)
    for n in range(2,k+2):
        out*=F(2*(n+1),n**3)
    return out


def trace_distance(a,b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2)


def reused_record_child(s):
    root=np.array([1,0,0,1])/math.sqrt(2) if s==1 else np.array([0,1,1,0])/math.sqrt(2)
    before=record_map(1)@root
    after=graph.phase(4,[(2,3)])*np.kron(before,graph.PLUS)
    local=np.moveaxis(after.reshape(4,2,2),1,0).reshape(2,-1)
    recovered=float(after@graph.pauli(4,x=(2,),z=(3,))@after)
    return local@local.T,recovered


def report():
    plus2=np.ones(4)/2
    known0=np.array([1,1,0,0])/math.sqrt(2)
    known1=np.array([0,0,1,1])/math.sqrt(2)
    parity0=np.array([1,0,0,1])/math.sqrt(2)
    parity1=np.array([0,1,1,0])/math.sqrt(2)
    rows=[]
    for k in range(1,5):
        a=record_map(k)@plus2
        xx=np.kron(graph.X,graph.X)
        roots=root_density_after_loss(plus2,k)
        z0=np.kron(graph.Z,graph.I)
        decohered=(roots+z0@roots@z0)/2
        rows.append({'children':k,'geometry_probability_exact':str(family_probability(k)),
                     'parity_probabilities':[float(plus2@p@plus2) for p in PROJECTORS],
                     'single_child_parity_trace_distance':trace_distance(child_density(parity0,k),child_density(parity1,k)),
                     'all_children_individual_root_bit_trace_distance':trace_distance(children_density(known0,k),children_density(known1,k)),
                     'logical_XX_expectation_after_ignoring_children':float(np.trace(roots@xx).real),
                     'logical_XX_after_individual_Z_measurement':float(np.trace(decohered@xx).real),
                     'normalization':float(np.vdot(a,a).real)})
    return {'round':250,'date':'2026-09-22',
            'model':'Two-root cover graph with k common children; specified CZ interaction',
            'recorded_observable':'J=Z_0 Z_1 (parity)',
            'record_isometry':'V_k=P_plus tensor |+>^k + P_minus tensor |->^k',
            'parent_channel_after_children_ignored':'rho -> P_plus rho P_plus + P_minus rho P_minus',
            'quantum_capacity_within_each_record_sector':'One arbitrary encoded qubit, including correlations with a reference',
            'results':rows,
            'record_reuse_counterexample':{
                'geometry_conditional_probability':'16/27 after the three-event star',
                'old_child_parity_trace_distance_after_new_CZ':trace_distance(reused_record_child(1)[0],reused_record_child(-1)[0]),
                'joint_Xold_Znew_expectations':[reused_record_child(s)[1] for s in (1,-1)]},
            'limits':['Conditioned graph family is not assumed typical; its probability is reported',
                      'Blank children, gates and persistent access to root memories are additional resources',
                      'Fixed-fragment records persist only while future interactions protect that fragment or its readout',
                      'Geometry probabilities still do not respond to matter states',
                      'No thermodynamic cost, QCH locality or gravitational dynamics is derived'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_sequential_gates_equal_parity_sector_formula(self):
        for k in range(1,6):
            v=record_map(k)
            np.testing.assert_allclose(v,sector_formula(k),atol=1e-13)
            np.testing.assert_allclose(v.T@v,np.eye(4),atol=1e-13)

    def test_matches_round249_cover_graph_states(self):
        for k in range(1,6):
            np.testing.assert_allclose(record_map(k)@(np.ones(4)/2),graph.graph_state(family(k)),atol=1e-13)

    def test_newborn_X_measurement_is_complete_parity_instrument(self):
        v=record_map(1).reshape(4,2,4)
        for s,p in zip((1,-1),PROJECTORS):
            bra=np.array([1.,s])/math.sqrt(2)
            operator=np.einsum('b,abc->ac',bra,v)
            np.testing.assert_allclose(operator,p,atol=1e-13)
        np.testing.assert_allclose(sum(p.T@p for p in PROJECTORS),np.eye(4),atol=1e-13)

    def test_redundant_children_agree_and_old_records_persist(self):
        psi=np.array([1,2j,3,-1j]);psi=psi/np.linalg.norm(psi)
        expected=[float(np.vdot(psi,p@psi).real) for p in PROJECTORS]
        for k in range(2,6):
            state=record_map(k)@psi
            for j in range(k):
                d=child_density(psi,k,j)
                for s,p in zip((1,-1),expected):
                    self.assertAlmostEqual(float(np.trace(d@(graph.I+s*graph.X)/2).real),p)
            for i in range(k):
                for j in range(i+1,k):
                    observable=graph.pauli(k+2,x=(i+2,j+2))
                    self.assertAlmostEqual(float(np.vdot(state,observable@state).real),1)

    def test_discarding_children_is_only_parity_dephasing(self):
        psi=np.array([1,2j,3,-1j]);psi=psi/np.linalg.norm(psi)
        rho=np.outer(psi,psi.conj())
        expected=sum(p@rho@p for p in PROJECTORS)
        for k in range(1,6):
            np.testing.assert_allclose(root_density_after_loss(psi,k),expected,atol=1e-13)

    def test_encoded_unknown_state_and_reference_survive(self):
        encoded=np.zeros((4,2),dtype=complex)
        encoded[0,0]=1/math.sqrt(2);encoded[3,1]=1j/math.sqrt(2)
        expected=np.outer(encoded.ravel(),encoded.ravel().conj())
        for k in range(1,6):
            out=(record_map(k)@encoded).reshape(4,2**k,2)
            reduced=np.einsum('acr,bcs->arbs',out,out.conj()).reshape(8,8)
            np.testing.assert_allclose(reduced,expected,atol=1e-13)

    def test_redundancy_does_not_reveal_individual_root_bit(self):
        a=np.array([1,1,0,0])/math.sqrt(2)
        b=np.array([0,0,1,1])/math.sqrt(2)
        for k in range(1,6):
            np.testing.assert_allclose(children_density(a,k),children_density(b,k),atol=1e-13)

    def test_each_child_perfectly_distinguishes_parity(self):
        even=np.array([1,0,0,1])/math.sqrt(2)
        odd=np.array([0,1,1,0])/math.sqrt(2)
        for k in range(1,6):
            for j in range(k):
                self.assertAlmostEqual(trace_distance(child_density(even,k,j),child_density(odd,k,j)),1)

    def test_within_sector_interference_survives_coarse_record(self):
        xx=np.kron(graph.X,graph.X); z0=np.kron(graph.Z,graph.I)
        for p in PROJECTORS:
            psi=p@(np.ones(4)/2);psi=psi/np.linalg.norm(psi)
            rho=root_density_after_loss(psi,3)
            dephased=(rho+z0@rho@z0)/2
            self.assertAlmostEqual(float(np.trace(rho@(np.eye(4)+xx)/2)),1)
            self.assertAlmostEqual(float(np.trace(dephased@(np.eye(4)+xx)/2)),.5)

    def test_conditioning_probability_is_not_hidden(self):
        for k in range(1,5):
            self.assertEqual(online.online_level(k+2)[family(k)],family_probability(k))

    def test_reusing_record_child_changes_local_accessibility(self):
        plus,ep=reused_record_child(1);minus,em=reused_record_child(-1)
        np.testing.assert_allclose(plus,np.eye(2)/2,atol=1e-13)
        np.testing.assert_allclose(minus,np.eye(2)/2,atol=1e-13)
        self.assertAlmostEqual(ep,1);self.assertAlmostEqual(em,-1)
        star=family(1)
        r=np.zeros((4,4),dtype=bool);r[:3,:3]=growth.matrix(star);r[:3,3]=True
        self.assertEqual(online.transitions(star)[growth.key(r)],F(16,27))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true'); args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful(): raise SystemExit(1)
    result=report(); result['checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        target=Path(__file__).with_name('parity_record_growth_audit_results.json')
        if target.exists() and target.read_text(encoding='utf-8')!=payload:
            raise SystemExit('Refusing to overwrite a different saved result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
