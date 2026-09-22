"""Round 249: graph-state matter on the prescribed classical growth process.

The construction adds one qubit per event and controlled-Z gates on cover
relations. It has readable geometry records and unrecorded matter interference.
It does not select a natural geometry law or establish full QCH locality.
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
import online_record_growth_audit as online
import random_order_growth_audit as growth
import covariant_minima_audit as minima

I=np.eye(2); X=np.array([[0.,1.],[1.,0.]]); Z=np.diag([1.,-1.])
PLUS=np.ones(2)/math.sqrt(2)


def pauli(n,x=(),z=()):
    out=np.ones((1,1))
    for j in range(n):
        out=np.kron(out,(X if j in x else I)@(Z if j in z else I))
    return out


@lru_cache(maxsize=None)
def bits(n):
    return np.array(list(itertools.product([0,1],repeat=n)),dtype=int)


def covers(h):
    r=growth.matrix(h)
    return tuple((i,j) for i in range(len(r)) for j in range(len(r))
                 if r[i,j] and not any(r[i,k] and r[k,j] for k in range(len(r))))


def phase(n,edges):
    b=bits(n)
    count=sum((b[:,i]*b[:,j] for i,j in edges),np.zeros(2**n,dtype=int))
    return 1-2*(count%2)


@lru_cache(maxsize=None)
def graph_state(h):
    return phase(len(h),covers(h))/math.sqrt(2**len(h))


@lru_cache(maxsize=None)
def append_map(child):
    n=len(child)
    edges=[(i,j) for i,j in covers(child) if j==n-1]
    return phase(n,edges)[:,None]*np.kron(np.eye(2**(n-1)),PLUS[:,None])


def stabilizer(h,v):
    neighbors=[j if i==v else i for i,j in covers(h) if v in (i,j)]
    return pauli(len(h),x=(v,),z=neighbors)


def grow_isometry(n):
    parents=tuple(online.online_level(n))
    children=tuple(online.online_level(n+1))
    index={h:j for j,h in enumerate(children)}
    oldq,newq=2**n,2**(n+1)
    v=np.zeros((len(children)*2*newq,len(parents)*2*oldq))
    for i,h in enumerate(parents):
        for child,t in online.transitions(h).items():
            j=index[child]
            flip=minima.minima_count(child)-minima.minima_count(h)
            a=math.sqrt(float(t))*append_map(child)
            for r in range(2):
                v[(2*j+(r^flip))*newq:(2*j+(r^flip)+1)*newq,
                  (2*i+r)*oldq:(2*i+r+1)*oldq]=a
    return v


def grown_state(n):
    state=np.concatenate([np.ones(4)/2,np.zeros(4)])
    for k in range(2,n):
        state=grow_isometry(k)@state
    return state


def frozen_record_probabilities(final_n,snapshot_n=3):
    state=grown_state(snapshot_n)
    q=2**snapshot_n
    saved=np.zeros((len(state),2))
    for i in range(len(state)):
        saved[i,(i//q)%2]=state[i]
    for n in range(snapshot_n,final_n):
        saved=grow_isometry(n)@saved
    return np.sum(saved**2,axis=0)


def fine_classes(h,target=None):
    n=len(h); target=n-1 if target is None else target
    a=append_map(h)
    k=stabilizer(h,target)
    zz=pauli(n,z=(target,))
    classes=[]
    for y in (1,-1):
        for z in (1,-1):
            classes.append((np.eye(2**n)+y*k)/2@(np.eye(2**n)+z*zz)/2@a)
    return classes


def gram(classes,psi):
    vectors=np.column_stack([c@psi for c in classes])
    return vectors.T@vectors


def chsh_star():
    h=next(h for h in online.online_level(3) if minima.minima_count(h)==2)
    psi=graph_state(h)
    # Bipartition: newborn qubit versus the two-parent subsystem.
    b=math.sqrt(2)*(pauli(3,x=(2,),z=(0,1))+pauli(3,x=(0,),z=(2,)))
    return float(psi@b@psi)


def report():
    rows=[]
    for n in range(2,6):
        keys=tuple(online.online_level(n))
        state=grown_state(n).reshape(len(keys),2,2**n)
        geometry=np.sum(state**2,axis=(1,2))
        old=np.array([float(v) for v in online.online_level(n).values()])
        stabilizer_error=max(abs(float(graph_state(h)@stabilizer(h,j)@graph_state(h))-1)
                             for h in keys for j in range(n))
        rows.append({'N':n,'supported_geometries':len(keys),
                     'geometry_probability_error':float(np.max(abs(geometry-old))),
                     'record_R1_probability':float(np.sum(state[:,1,:]**2)),
                     'max_stabilizer_error':stabilizer_error})
    star=next(h for h in online.online_level(3) if minima.minima_count(h)==2)
    d=gram(fine_classes(star),np.ones(4)/2)
    return {'round':249,'date':'2026-09-22','model':'Prescribed online geometry plus cover-graph qubits and geometry latch',
            'finite_growth':rows,'star_geometry_probability':.75,
            'conditional_fine_history_D':d.tolist(),
            'conditional_coherent_stabilizer_plus':float(d[:2,:2].sum()),
            'conditional_plus_after_Z_record':float(np.trace(d[:2,:2])),
            'frozen_N3_record_distribution_at_N5':frozen_record_probabilities(5).tolist(),
            'star_CHSH_newborn_vs_parents':chsh_star(),
            'extra_inputs':['One qubit per event in |+>','CZ interaction on every cover relation',
                            'Round247 geometry probabilities and global controller',
                            'Specified record access and matter measurement'],
            'limits':['Geometry remains the given classical pointer law',
                      'CZ preparation commutation does not prove full QCH axioms for arbitrary interventions',
                      'Only specified record marginals, not an infinite coherent history functional, are extended'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_append_isometry_and_graph_preparation(self):
        for n in range(3,6):
            for h in online.online_level(n):
                a=append_map(h)
                np.testing.assert_allclose(a.T@a,np.eye(2**(n-1)),atol=1e-13)
                parent=tuple(tuple(row[:-1]) for row in h[:-1])
                np.testing.assert_allclose(a@graph_state(parent),graph_state(h),atol=1e-13)

    def test_full_growth_isometry(self):
        for n in range(2,5):
            v=grow_isometry(n)
            np.testing.assert_allclose(v.T@v,np.eye(v.shape[1]),atol=1e-13)

    def test_full_state_matches_closed_form(self):
        for n in range(2,6):
            law=online.online_level(n)
            state=grown_state(n).reshape(len(law),2,2**n)
            for j,(h,p) in enumerate(law.items()):
                r=minima.minima_count(h)-2
                np.testing.assert_allclose(state[j,r],math.sqrt(float(p))*graph_state(h),atol=1e-13)
                np.testing.assert_allclose(state[j,1-r],0,atol=1e-13)

    def test_stabilizer_and_dephasing_operational_contrast(self):
        for h in online.online_level(4):
            psi=graph_state(h)
            for j in range(4):
                k=stabilizer(h,j); z=pauli(4,z=(j,))
                self.assertAlmostEqual(float(psi@k@psi),1)
                self.assertAlmostEqual(float((psi@k@psi+psi@z@k@z@psi)/2),0)

    def test_fine_class_completeness_and_coherent_sum(self):
        for h in online.online_level(3):
            cs=fine_classes(h)
            np.testing.assert_allclose(sum(c.T@c for c in cs),np.eye(4),atol=1e-13)
            np.testing.assert_allclose(sum(cs),append_map(h),atol=1e-13)
            d=gram(cs,np.ones(4)/2)
            self.assertAlmostEqual(float(d.sum()),1)
            self.assertAlmostEqual(float(np.trace(d)),1)
            self.assertAlmostEqual(float(d[:2,:2].sum()),1)
            self.assertAlmostEqual(float(np.trace(d[:2,:2])),.5)
            self.assertAlmostEqual(float(d[0,1]),.25)

    def test_passive_relabelling_of_states_and_marked_queries(self):
        for h in online.online_level(5):
            r=growth.matrix(h)
            for perm in itertools.permutations(range(5)):
                rr=r[np.ix_(perm,perm)]
                if np.any(np.tril(rr)):
                    continue
                hp=growth.key(rr)
                transformed=graph_state(h).reshape((2,)*5).transpose(perm).ravel()
                np.testing.assert_allclose(transformed,graph_state(hp),atol=1e-13)
                self.assertEqual(online.online_level(5)[h],online.online_level(5)[hp])
                v=perm.index(4)
                self.assertAlmostEqual(float(transformed@stabilizer(hp,v)@transformed),1)

    def test_old_snapshot_record_is_preserved(self):
        for n in (3,4,5):
            np.testing.assert_allclose(frozen_record_probabilities(n),[.75,.25],atol=1e-13)

    def test_read_geometry_record_does_not_remove_conditional_matter_fringe(self):
        n=4; law=online.online_level(n)
        state=grown_state(n).reshape(len(law),2,2**n)
        for r in range(2):
            probability=sum(float(state[j,r]@state[j,r]) for j in range(len(law)))
            fringe=sum(float(state[j,r]@stabilizer(h,n-1)@state[j,r]) for j,h in enumerate(law))
            self.assertAlmostEqual(fringe,probability)

    def test_entanglement_is_attached_to_cover_geometry(self):
        self.assertAlmostEqual(chsh_star(),2*math.sqrt(2))
        star=next(h for h in online.online_level(3) if minima.minima_count(h)==2)
        psi=graph_state(star).reshape(4,2)
        np.testing.assert_allclose(psi.T@psi,np.eye(2)/2,atol=1e-13)

    def test_incomparable_preparation_order(self):
        # Two roots, then two incomparable children with the same parents.
        h=((False,False,True,True),(False,False,True,True),
           (False,False,False,False),(False,False,False,False))
        p1=phase(4,[(0,2),(1,2)])
        p2=phase(4,[(0,3),(1,3)])
        np.testing.assert_allclose(p1*p2/4,graph_state(h),atol=1e-13)
        np.testing.assert_allclose(p2*p1/4,graph_state(h),atol=1e-13)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true'); args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful(): raise SystemExit(1)
    result=report(); result['checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        target=Path(__file__).with_name('graph_record_growth_audit_results.json')
        if target.exists() and target.read_text(encoding='utf-8')!=payload:
            raise SystemExit('Refusing to overwrite a different saved result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
