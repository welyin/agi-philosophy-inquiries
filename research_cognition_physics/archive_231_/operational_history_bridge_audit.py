"""Round 248: audit the finite-instrument meaning of growth functionals.

Trace normalization of an exhaustive Kraus-history Gram matrix is necessary.
It differs from D(Omega,Omega)=1. A cited isometric amplitude model supplies
a second counterexample even when both scalar normalizations hold.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import record_sector_growth_audit as sector
import online_record_growth_audit as online


def gram(vectors):
    return vectors.conj().T@vectors


def exhaustive_preparation(n,p=.5):
    # Branch operators C -> C^3 on a fixed one-dimensional input.
    _,v = sector.vectors(n,p)
    q = float(np.trace(gram(v)).real)
    extended = np.zeros((3,v.shape[1]+1))
    extended[:2,:-1] = v
    extended[2,-1] = math.sqrt(max(0,1-q))
    return extended,q


def uta(theta):
    return np.exp(1j*theta)*np.array([math.cos(theta),-1j*math.sin(theta)])


def tree_isometry(depth,theta):
    c = uta(theta)
    v = np.zeros((2**(depth+1),2**depth),dtype=complex)
    for j in range(2**depth):
        v[2*j:2*j+2,j] = c
    return v


def amplitudes(depth,theta=math.pi/4):
    c = uta(theta)
    return np.array([np.prod([c[j] for j in bits]) for bits in itertools.product(range(2),repeat=depth)])


def event_values(a,subset):
    chi = np.zeros(len(a))
    chi[list(subset)] = 1
    f = np.outer(chi,chi)
    mu = float(abs(chi@a)**2)
    return {'cardinality':len(subset),'quantum_measure':mu,
            'diagonal_projector_probability':float(np.sum(abs(a[list(subset)])**2)),
            'unnormalized_effect_norm':float(np.linalg.eigvalsh(f)[-1]),
            'normalized_coherent_projector_probability':mu/len(subset)}


def report():
    rows = []
    for n in range(3,6):
        vectors,q = exhaustive_preparation(n)
        full = gram(vectors)
        rows.append({'N':n,'trace_D_rec':q,
                     'exhaustive_instrument_trace_deficit':1-q,
                     'completed_instrument_trace':float(np.trace(full).real),
                     'orthogonal_failure_completion_total_bare_D':float(full.sum().real),
                     'expected_total_bare_D':2-q})
    a = amplitudes(3)
    middle = [j for j,bits in enumerate(itertools.product(range(2),repeat=3)) if sum(bits) in (1,2)]
    return {'round':248,'date':'2026-09-22',
            'instrument_necessary_condition':'For normalized input and exhaustive Kraus histories: sum_h D(h,h)=1',
            'scope_of_obstruction':'Exact D_rec on the specified exhaustive finite prefix partition; not all generalized quantum measures or all coarse models',
            'round246_completion_audit':rows,
            'cited_isometric_model':{
                'source':'https://arxiv.org/abs/1409.3770',
                'theta':'pi/4','binary_growth_steps':3,
                'sum_amplitudes':[float(a.sum().real),float(a.sum().imag)],
                'squared_norm':float(np.vdot(a,a).real),
                'scalar_D_rank':int(np.linalg.matrix_rank(np.outer(a.conj(),a))),
                'projector_history_D_rank':len(a),
                'six_middle_paths':event_values(a,middle),
                'whole_event_coherent_projection_probability':float(abs(a.sum())**2/len(a))},
            'conclusion':'An isometry does not identify a scalar amplitude-sum event functional with a physical pointer-history functional.',
            'next_interface':'Specify event effects, class operators, actual record algebra and complete interventions before selecting growth dynamics.',
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_round246_fails_exhaustive_trace_normalization(self):
        for p in [.2,.5,.8]:
            for n in range(3,6):
                d = sector.functional(n,p)[1]
                self.assertAlmostEqual(float(d.sum()),1)
                self.assertLess(float(np.trace(d)),1-1e-12)
            self.assertAlmostEqual(float(np.trace(sector.functional(3,p)[1])),(1+p)/2)

    def test_missing_outcome_completes_instrument_but_not_old_total_D(self):
        for n in range(3,6):
            v,q = exhaustive_preparation(n)
            d = gram(v)
            self.assertAlmostEqual(float(np.trace(d)),1)
            np.testing.assert_allclose(d[:-1,:-1],sector.functional(n)[1],atol=1e-13)
            self.assertAlmostEqual(float(d.sum()),2-q)
            self.assertGreater(float(d.sum()),1)
            # Summing all CP branch effects gives I on the preparation input.
            self.assertAlmostEqual(sum(float(np.vdot(v[:,j],v[:,j]).real) for j in range(v.shape[1])),1)

    def test_online_projective_histories_satisfy_both_normalizations(self):
        for n in range(2,6):
            weights=np.array([float(v) for v in online.online_level(n).values()])
            branch_vectors=np.diag(np.sqrt(weights))
            d=gram(branch_vectors)
            self.assertAlmostEqual(float(np.trace(d)),1)
            self.assertAlmostEqual(float(d.sum()),1)

    def test_cited_uta_has_both_local_normalizations(self):
        for theta in [.1,.4,math.pi/4,1.2,2.4]:
            c=uta(theta)
            self.assertAlmostEqual(abs(c.sum()-1),0)
            self.assertAlmostEqual(float(np.vdot(c,c).real),1)
            for depth in range(1,5):
                a=amplitudes(depth,theta)
                self.assertAlmostEqual(abs(a.sum()-1),0)
                self.assertAlmostEqual(float(np.vdot(a,a).real),1)

    def test_isometry_generates_the_same_amplitudes(self):
        state=np.ones(1,dtype=complex)
        for depth in range(4):
            v=tree_isometry(depth,.63)
            np.testing.assert_allclose(v.conj().T@v,np.eye(v.shape[1]),atol=1e-13)
            state=v@state
            np.testing.assert_allclose(state,amplitudes(depth+1,.63),atol=1e-13)

    def test_same_state_different_history_functionals(self):
        a=amplitudes(2)
        scalar=np.outer(a.conj(),a)
        operational=gram(np.diag(a))
        self.assertEqual(np.linalg.matrix_rank(scalar),1)
        self.assertEqual(np.linalg.matrix_rank(operational),4)
        for d in [scalar,operational]:
            self.assertAlmostEqual(float(np.trace(d).real),1)
            self.assertAlmostEqual(float(d.sum().real),1)
        self.assertAlmostEqual(float(scalar[np.ix_([1,2],[1,2])].sum().real),1)
        self.assertAlmostEqual(float(operational[np.ix_([1,2],[1,2])].sum().real),.5)

    def test_quantum_measure_above_one_despite_unit_norm(self):
        a=amplitudes(3)
        middle=[1,2,3,4,5,6]
        values=event_values(a,middle)
        self.assertAlmostEqual(values['quantum_measure'],2.25)
        self.assertAlmostEqual(values['diagonal_projector_probability'],.75)

    def test_unscaled_event_operator_is_not_an_effect(self):
        a=amplitudes(3)
        subset=[1,2,3,4,5,6]
        values=event_values(a,subset)
        self.assertAlmostEqual(values['unnormalized_effect_norm'],6)
        self.assertAlmostEqual(values['normalized_coherent_projector_probability'],.375)
        chi=np.zeros(8);chi[subset]=1
        e=np.outer(chi,chi)/6
        np.testing.assert_allclose(e@e,e,atol=1e-13)
        self.assertGreaterEqual(np.linalg.eigvalsh(np.eye(8)-e).min(),-1e-13)
        self.assertAlmostEqual(float(np.vdot(a,e@a).real),.375)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    check=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not check.wasSuccessful(): raise SystemExit(1)
    result=report()
    result['checks']={'run':check.testsRun,'failures':len(check.failures),'errors':len(check.errors)}
    payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        target=Path(__file__).with_name('operational_history_bridge_audit_results.json')
        if target.exists() and target.read_text(encoding='utf-8')!=payload:
            raise SystemExit('Refusing to overwrite a different saved result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
