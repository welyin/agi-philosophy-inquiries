"""Round 253: local quantum instruments rewrite a live port into one or two.

Port identities and fresh ancillas are supplied. Disjoint raw rewrites commute;
uniform global next-port scheduling need not assign equal order probabilities.
"""
import argparse
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
from record_preserving_feedback_audit import weak

I=np.eye(2); X=np.array([[0.,1.],[1.,0.]]); Z=np.diag([1.,-1.])
ROOTS=((0,),(1,))


def observable(address):
    return X if (len(address)-1)%2==0 else Z


def split_embedding(d=2):
    j=np.zeros((d*d,d))
    for x in range(d):
        fresh=0 if d==2 else 2*(x//2)
        j[d*x+fresh,x]=1
    return j


def local_map(sign,eta=.6,q=X,d=2):
    k=weak(q,sign,eta)
    if d==4: k=np.kron(I,k)
    return split_embedding(d)@k if sign==1 else k


def advance(frontier,matrix,address,sign,eta=.6,q=None,d=2):
    """Apply a rectangular local map and canonically order output port factors."""
    frontier=tuple(frontier); at=frontier.index(address)
    q=observable(address) if q is None else q
    k=local_map(sign,eta,q,d)
    count=2 if sign==1 else 1
    children=tuple(address+(i,) for i in range(count))
    others=frontier[:at]+frontier[at+1:]
    old_axes=(at,)+tuple(i for i in range(len(frontier)) if i!=at)+(len(frontier),)
    tensor=matrix.reshape((d,)*len(frontier)+(matrix.shape[1],)).transpose(old_axes)
    tensor=(k@tensor.reshape(d,-1)).reshape((d,)*(count+len(others))+(matrix.shape[1],))
    temporary=children+others; new=tuple(sorted(temporary))
    axes=tuple(temporary.index(a) for a in new)+(len(new),)
    return new,tensor.transpose(axes).reshape(d**len(new),matrix.shape[1])


def pair_order(signs,reverse=False,eta=.6,d=2,global_strength=False):
    frontier=ROOTS; matrix=np.eye(d**2)
    for i in ((1,0) if reverse else (0,1)):
        strength=eta/len(frontier) if global_strength else eta
        frontier,matrix=advance(frontier,matrix,ROOTS[i],signs[i],strength,d=d)
    return frontier,matrix


def schedule_factor(signs,reverse=False):
    # Exactly these two roots are chosen in the first two global steps.
    m=2; factor=F(1)
    for i in ((1,0) if reverse else (0,1)):
        factor/=m
        m+=int(signs[i]==1)
    return factor


def report():
    rows=[]
    for signs in itertools.product((1,-1),repeat=2):
        f,a=pair_order(signs); ff,b=pair_order(signs,True)
        assert f==ff
        rows.append({'outcomes':list(signs),'output_ports':len(f),
                     'raw_operator_diamond_error':float(np.max(abs(a-b))),
                     'uniform_schedule_forward_factor':str(schedule_factor(signs)),
                     'uniform_schedule_reverse_factor':str(schedule_factor(signs,True))})
    _,a=pair_order((1,-1),global_strength=True)
    _,b=pair_order((1,-1),True,global_strength=True)
    return {'round':253,'date':'2026-09-22',
            'rule':'+ consumes one port and emits carrier plus fresh blank; - emits carrier only',
            'local_setting':'X on even local depth, Z on odd local depth; eta=0.6',
            'diamonds':rows,
            'mixed_transcript_probabilities_for_X_plus_inputs':['2/75','1/25'],
            'frontier_dependent_strength_raw_diamond_error':float(np.max(abs(a-b))),
            'coupling_constraint':'For eta in [0,1), eta(m)=eta(m+1) is necessary for the mixed raw diamond in the shared strength ansatz',
            'resources':'One fresh logical qubit for each split; classical record sector copied, unknown logical state not copied',
            'limits':['Marked roots and directed port genealogy are inputs',
                      'Locality is disjoint tensor support, not emergent spatial distance',
                      'One-parent rule generates a forest; no merging, interactions between branches or spacetime limit',
                      'Uniform scheduling is an additional stochastic process, not passive relabelling'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_local_instruments_on_variable_output_spaces(self):
        for d in (2,4):
            for eta in (0,.2,.6,1):
                for angle in (0,.4,1.2):
                    q=math.cos(angle)*X+math.sin(angle)*Z
                    ks=[local_map(g,eta,q,d) for g in (1,-1)]
                    np.testing.assert_allclose(sum(k.T@k for k in ks),np.eye(d),atol=1e-13)

    def test_copied_record_values_are_preserved(self):
        record=np.kron(Z,I)
        for g in (1,-1):
            k=local_map(g,d=4)
            outputs=[record] if g==-1 else [np.kron(record,np.eye(4)),np.kron(np.eye(4),record)]
            for out in outputs: np.testing.assert_allclose(out@k,k@record,atol=1e-13)

    def test_split_retains_unknown_logical_state_and_reference_in_one_sector(self):
        j=split_embedding(4)
        for s in (0,1):
            v=np.zeros(8,complex); v[4*s]=1/math.sqrt(2); v[4*s+3]=1j/math.sqrt(2)
            out=(np.kron(j,I)@v).reshape(4,4,2)
            rho=np.einsum('abr,cbd->arcd',out,out.conj()).reshape(8,8)
            np.testing.assert_allclose(rho,np.outer(v,v.conj()),atol=1e-13)

    def test_all_raw_disjoint_diamonds(self):
        for d in (2,4):
            for eta in (.2,.6,.9):
                for signs in itertools.product((1,-1),repeat=2):
                    f,a=pair_order(signs,eta=eta,d=d); ff,b=pair_order(signs,True,eta,d)
                    self.assertEqual(f,ff); np.testing.assert_allclose(a,b,atol=1e-13)

    def test_locally_different_settings_after_a_split(self):
        f,a=advance(ROOTS,np.eye(4),ROOTS[0],1)
        u,v=f[0],ROOTS[1]
        for signs in itertools.product((1,-1),repeat=2):
            f1,a1=advance(f,a,u,signs[0],.3)
            f1,a1=advance(f1,a1,v,signs[1],.8)
            f2,a2=advance(f,a,v,signs[1],.8)
            f2,a2=advance(f2,a2,u,signs[0],.3)
            self.assertEqual(f1,f2); np.testing.assert_allclose(a1,a2,atol=1e-13)
            self.assertEqual(len(f1),len(set(f1)))

    def test_uniform_next_port_choice_is_complete(self):
        for f in (ROOTS,((0,0),(0,1),(1,))):
            m=len(f); identity=np.eye(2**m); total=np.zeros_like(identity)
            for address in f:
                for sign in (1,-1):
                    _,k=advance(f,identity,address,sign)
                    total+=k.T@k/m
            np.testing.assert_allclose(total,identity,atol=1e-13)

    def test_uniform_scheduler_breaks_equal_transcript_weights(self):
        signs=(1,-1); psi=np.ones(4)/2
        _,a=pair_order(signs); _,b=pair_order(signs,True)
        self.assertEqual(schedule_factor(signs),F(1,6))
        self.assertEqual(schedule_factor(signs,True),F(1,4))
        self.assertAlmostEqual(float(np.linalg.norm(a@psi)**2)*float(schedule_factor(signs)),float(F(2,75)))
        self.assertAlmostEqual(float(np.linalg.norm(b@psi)**2)*float(schedule_factor(signs,True)),float(F(1,25)))

    def test_remote_nonselective_operation_leaves_local_birth_effect(self):
        for eta in (.2,.9):
            local=local_map(1,.6,X)
            expected=np.kron(local.T@local,I)
            effect=sum((k:=np.kron(local,local_map(g,eta,Z))).T@k for g in (1,-1))
            np.testing.assert_allclose(effect,expected,atol=1e-13)

    def test_global_strength_has_an_operator_level_obstruction(self):
        _,a=pair_order((1,-1),global_strength=True)
        _,b=pair_order((1,-1),True,global_strength=True)
        self.assertGreater(np.max(abs(a-b)),.01)
        # Direct tensor expression independently identifies the obstruction.
        expected=np.kron(local_map(1,.3),local_map(-1,.2)-local_map(-1,.3))
        np.testing.assert_allclose(a-b,expected,atol=1e-13)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true'); args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful(): raise SystemExit(1)
    result=report(); result['checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        target=Path(__file__).with_name('local_branching_growth_audit_results.json')
        if target.exists() and target.read_text(encoding='utf-8')!=payload: raise SystemExit('Refusing to overwrite a different result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
