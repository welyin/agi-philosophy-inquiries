"""Round 252: quantum-sensitive future motifs with protected parity records.

The circuit, weak instrument, motifs and routing are explicit inputs. This is
finite operational feedback, not an autonomous causal-set growth dynamics.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import graph_record_growth_audit as graph
import parity_record_growth_audit as parity

ID=np.eye(4)
S=graph.pauli(2,z=(0,1))
LX=graph.pauli(2,x=(0,1))
LZ=graph.pauli(2,z=(0,))
SIGNS=(1,-1)


def density(v): return np.outer(v,np.conj(v))


def weak(q,g,eta):
    if not 0<=eta<=1 or g not in SIGNS: raise ValueError('Invalid instrument parameters')
    eye=np.eye(len(q))
    return math.sqrt((1+g*eta)/2)*(eye+q)/2+math.sqrt((1-g*eta)/2)*(eye-q)/2


def channel(q,eta,rho):
    return sum((k:=weak(q,g,eta))@rho@k.conj().T for g in SIGNS)


def sequential_probs(rho,eta,x_first=True):
    out=np.zeros((2,2))
    for i,a in enumerate(SIGNS):
        for j,b in enumerate(SIGNS):
            x,z=weak(LX,a,eta),weak(LZ,b,eta)
            k=z@x if x_first else x@z
            out[i,j]=np.trace(k@rho@k.conj().T).real
    return out


def order_formula(rho,eta,x_first=True):
    rx=np.trace(rho@LX).real; rz=np.trace(rho@LZ).real
    gamma=math.sqrt(1-eta**2)
    return np.array([[(1+a*eta*rx+b*eta*gamma*rz)/4 if x_first else
                      (1+a*eta*gamma*rx+b*eta*rz)/4 for b in SIGNS] for a in SIGNS])


def motif_edges(sign,offset=0):
    # U queries the incoming memory. Both subsequent events are in its future.
    return [(offset,offset+1),(offset,offset+2)] if sign==1 else [(offset,offset+1),(offset+1,offset+2)]


def orders(signs):
    edges=motif_edges(signs[0])+motif_edges(signs[1],3)
    return [order for order in itertools.permutations(range(6))
            if all(order.index(a)<order.index(b) for a,b in edges)]


def branch_operator(order,signs,eta):
    a=np.kron(weak(LX,signs[0],eta),ID)
    b=np.kron(ID,weak(LX,signs[1],eta))
    result=np.eye(16)
    for event in order:
        if event==0: result=a@result
        elif event==3: result=b@result
    return result


def graph_order(sign):
    # Two source events, then query U, then a fork or a chain.
    r=np.zeros((5,5),dtype=bool)
    for a,b in [(0,2),(1,2)]+motif_edges(sign,2): r[a,b]=True
    for k in range(5): r|=r[:,k,None]&r[k,None,:]
    return r


def motif_probabilities(rho,eta,q=LX):
    return [float(np.trace(weak(q,g,eta)@rho@weak(q,g,eta)).real) for g in SIGNS]


def report():
    bell_plus=np.array([1,0,0,1])/math.sqrt(2)
    bell_minus=np.array([1,0,0,-1])/math.sqrt(2)
    mixture=np.diag([.5,0,0,.5])
    rows=[]
    for eta in (0.,.2,.6,1.):
        gamma=math.sqrt(1-eta**2)
        p=sequential_probs(density(np.array([1,0,0,0])),eta)
        q=sequential_probs(density(np.array([1,0,0,0])),eta,False)
        rows.append({'eta':eta,'fork_probability_logical_X_plus':motif_probabilities(density(bell_plus),eta)[0],
                     'fork_probability_logical_X_minus':motif_probabilities(density(bell_minus),eta)[0],
                     'fork_probability_dephased_state':motif_probabilities(mixture,eta)[0],
                     'retained_logical_Z':float(np.trace(channel(LX,eta,np.diag([1.,0,0,0]))@LZ).real),
                     'predicted_retained_logical_Z':gamma,
                     'encoded_reference_entanglement_fidelity':(1+gamma)/2,
                     'maximum_order_TV':float(np.sum(abs(p-q))/2)})
    psi=np.zeros(16); psi[0]=psi[-1]=1/math.sqrt(2)
    concurrent=[]
    for signs in itertools.product(SIGNS,repeat=2):
        variants=orders(signs); maps=[branch_operator(o,signs,.6) for o in variants]
        concurrent.append({'outcomes':list(signs),'legal_birth_orders':len(variants),
                           'branch_operator_order_error':max(float(np.max(abs(m-maps[0]))) for m in maps),
                           'probability_for_encoded_Bell_input':float(np.linalg.norm(maps[0]@psi)**2)})
    return {'round':252,'date':'2026-09-22',
            'record_algebra':'S=Z_A Z_B; persistent memory at current boundary',
            'feedback_observables':['X_L=X_A X_B','Z_L=Z_A'],
            'instrument':'K_g(Q)=sqrt((I+g eta Q)/2), fixed eta and setting Q',
            'motif_dictionary':{'+1':'U precedes incomparable V,W','-1':'U precedes V precedes W'},
            'tradeoff_and_order_checks':rows,'two_disjoint_blocks':concurrent,
            'theorem_scope':'Exact identity channel within each record sector forbids informative outputs about its internal state',
            'limits':['Controlled future motifs and interaction parameters are supplied',
                      'Order covariance is verified for the specified finite circuit and its legal interleavings',
                      'No general state-dependent CSG kernel, full QCH classification, continuum or gravity is derived',
                      'A single fixed output distribution is not a witness against all classical hidden-state models'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_positive_square_roots_and_completeness(self):
        for eta in (0.,.2,.6,1.):
            for q in (LX,LZ):
                ks=[weak(q,g,eta) for g in SIGNS]
                for g,k in zip(SIGNS,ks):
                    effect=(ID+g*eta*q)/2
                    values,vectors=np.linalg.eigh(effect)
                    independent=(vectors*np.sqrt(np.maximum(0,values)))@vectors.conj().T
                    np.testing.assert_allclose(k,independent,atol=1e-13)
                np.testing.assert_allclose(sum(k.conj().T@k for k in ks),ID,atol=1e-13)

    def test_each_branch_preserves_record_sectors(self):
        for q in (LX,LZ):
            for eta in (.2,.6,1.):
                for g in SIGNS:
                    k=weak(q,g,eta)
                    np.testing.assert_allclose(k@S,S@k,atol=1e-13)
                for s in SIGNS:
                    p=(ID+s*S)/2
                    np.testing.assert_allclose(sum(weak(q,g,eta)@p@weak(q,g,eta) for g in SIGNS),p,atol=1e-13)

    def test_internal_phase_changes_geometry_with_same_old_record(self):
        for sign in SIGNS:
            v=np.array([1,0,0,sign])/math.sqrt(2)
            self.assertAlmostEqual(float(v@S@v),1)
            self.assertAlmostEqual(motif_probabilities(density(v),.6)[0],(1+sign*.6)/2)
        self.assertAlmostEqual(motif_probabilities(np.diag([.5,0,0,.5]),.6)[0],.5)

    def test_exact_information_disturbance_tradeoff(self):
        for eta in np.linspace(0,1,11):
            gamma=math.sqrt(1-eta**2)
            for q,r in ((LX,LZ),(LZ,LX)):
                rho=(ID+q+.4*S+.4*S@q)/4
                # Dephasing in the r basis damps the q expectation by gamma.
                out=channel(r,eta,rho)
                self.assertAlmostEqual(np.trace(out@q).real,gamma)
            self.assertAlmostEqual(eta**2+gamma**2,1)

    def test_unknown_reference_is_disturbed_by_informative_feedback(self):
        v=np.zeros(8,complex); v[0]=1/math.sqrt(2); v[7]=1j/math.sqrt(2)
        rho=density(v)
        for eta in (0.,.2,.6,1.):
            out=sum((k:=np.kron(weak(LX,g,eta),np.eye(2)))@rho@k.conj().T for g in SIGNS)
            fidelity=float((v.conj()@out@v).real)
            self.assertAlmostEqual(fidelity,(1+math.sqrt(1-eta**2))/2)

    def test_old_redundant_records_and_frozen_outcome(self):
        v=parity.record_map(2)@np.ones(4)/2
        rho=density(v)
        before=np.einsum('arbr->ab',rho.reshape(4,4,4,4))
        output=np.zeros_like(rho)
        for g in SIGNS:
            first=np.kron(weak(LX,g,.6),np.eye(4))
            branch=first@rho@first.conj().T
            frozen_prob=np.trace(branch)
            evolved=sum((k:=np.kron(weak(LZ,b,.6),np.eye(4)))@branch@k.conj().T for b in SIGNS)
            self.assertAlmostEqual(np.trace(evolved),frozen_prob)
            output+=evolved
        for s in SIGNS:
            root=(ID+s*S)/2
            wrong=(np.eye(2)-s*graph.X)/2
            mismatch=np.kron(root,np.kron(wrong,np.eye(2)))
            self.assertAlmostEqual(np.trace(mismatch@output).real,0)
        old_records=np.einsum('aras->rs',rho.reshape(4,4,4,4))
        new_records=np.einsum('aras->rs',output.reshape(4,4,4,4))
        np.testing.assert_allclose(old_records,new_records,atol=1e-13)

    def test_every_legal_concurrent_birth_order(self):
        counts=[]; completeness=np.zeros((16,16))
        for signs in itertools.product(SIGNS,repeat=2):
            variants=orders(signs); expected=np.kron(weak(LX,signs[0],.6),weak(LX,signs[1],.6))
            for order in variants: np.testing.assert_allclose(branch_operator(order,signs,.6),expected,atol=1e-13)
            completeness+=expected.T@expected; counts.append(len(variants))
        self.assertEqual(counts,[80,40,40,20])
        np.testing.assert_allclose(completeness,np.eye(16),atol=1e-13)

    def test_concurrent_entangled_input_and_no_remote_setting_signal(self):
        v=np.zeros(16); v[0]=v[-1]=1/math.sqrt(2); rho=density(v)
        for i,a in enumerate(SIGNS):
            for j,b in enumerate(SIGNS):
                k=np.kron(weak(LX,a,.6),weak(LX,b,.6))
                self.assertAlmostEqual(np.trace(k@rho@k.T), (1+a*b*.36)/4)
            for remote in (LX,LZ):
                prob=sum(np.trace((k:=np.kron(weak(LX,a,.6),weak(remote,b,.6)))@rho@k.T) for b in SIGNS)
                self.assertAlmostEqual(prob,.5)

    def test_noncommuting_order_formula_for_complex_inputs(self):
        for eta in (.2,.6,1.):
            for theta in (.2,.7,1.2):
                for phi in (0,.4,1.5):
                    v=np.array([math.cos(theta),0,0,np.exp(1j*phi)*math.sin(theta)])
                    rho=density(v)
                    for first in (True,False):
                        np.testing.assert_allclose(sequential_probs(rho,eta,first),order_formula(rho,eta,first),atol=1e-13)

    def test_order_total_variation_bound_is_attained(self):
        for eta in (.2,.6,1.):
            rho=np.diag([1.,0,0,0]); first=sequential_probs(rho,eta); second=sequential_probs(rho,eta,False)
            self.assertAlmostEqual(np.sum(abs(first-second))/2,eta*(1-math.sqrt(1-eta**2))/2)
            for v in (np.array([1,0,0,1])/math.sqrt(2),np.array([1,0,0,1j])/math.sqrt(2)):
                rho=density(v)
                rx=abs(np.trace(rho@LX)); rz=abs(np.trace(rho@LZ))
                tv=np.sum(abs(sequential_probs(rho,eta)-sequential_probs(rho,eta,False)))/2
                self.assertAlmostEqual(tv,eta*(1-math.sqrt(1-eta**2))*max(rx,rz)/2)

    def test_generated_motifs_are_distinct_common_future_posets(self):
        counts=[]
        for sign in SIGNS:
            r=graph_order(sign)
            self.assertFalse(np.any(np.diag(r)))
            self.assertTrue(np.all(r[:2,2:]))
            self.assertFalse(np.any(r[2:,:2]))
            counts.append(int(r.sum()))
        self.assertEqual(counts,[8,9])

    def test_noncommuting_with_old_record_fails_protection(self):
        q=graph.pauli(2,x=(0,)); rho=np.diag([1.,0,0,0])
        self.assertAlmostEqual(np.trace(channel(q,.6,rho)@S),.8)
        self.assertAlmostEqual(np.trace(channel(LX,.6,rho)@S),1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true'); args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful(): raise SystemExit(1)
    result=report(); result['checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        target=Path(__file__).with_name('record_preserving_feedback_audit_results.json')
        if target.exists() and target.read_text(encoding='utf-8')!=payload:
            raise SystemExit('Refusing to overwrite a different saved result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
