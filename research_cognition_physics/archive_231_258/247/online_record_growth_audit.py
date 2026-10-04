"""Round 247: an online realization of round 246's classical comparator.

No final-sector label is sampled or stored initially. The prescribed hazard,
uniform pair choice and geometry memory remain additional model inputs.
The associated pointer-history functional is classical, not D_rec.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import csg_support_audit as csg
import random_order_growth_audit as growth
import covariant_minima_audit as minima
import record_sector_growth_audit as sector


def hazard(n,m,p=F(1,2)):
    if n < 2 or m not in (2,3) or not 0 < p < 1:
        raise ValueError('Audited support is n>=2, m=2 or 3, 0<p<1')
    return (1-p)/(n*(p*(n-1)+1-p)) if m == 2 else F(0)


def pair_closure(r,pair):
    return tuple(j for j in range(len(r)) if j in pair or any(r[j,k] for k in pair))


def transitions(parent,p=F(1,2)):
    r = growth.matrix(parent)
    n = len(r)
    h = hazard(n,minima.minima_count(parent),p)
    out = defaultdict(F)
    if h:
        out[growth.key(csg.append_event(r,()))] += h
    # Enumerate actual pairs; this is independent of the old lambda formula.
    for pair in itertools.combinations(range(n),2):
        child = growth.key(csg.append_event(r,pair_closure(r,pair)))
        out[child] += (1-h)/math.comb(n,2)
    return dict(out)


@lru_cache(maxsize=None)
def online_level(n,p=F(1,2)):
    if n == 2:
        return {growth.key(np.zeros((2,2),dtype=bool)):F(1)}
    out = defaultdict(F)
    for parent,weight in online_level(n-1,p).items():
        for child,t in transitions(parent,p).items():
            out[child] += weight*t
    return dict(out)


def geometry_isometry(n,p=F(1,2)):
    parents = tuple(online_level(n,p))
    children = tuple(online_level(n+1,p))
    index = {h:j for j,h in enumerate(children)}
    w = np.zeros((len(children),len(parents)))
    for i,h in enumerate(parents):
        for child,t in transitions(h,p).items():
            w[index[child],i] = math.sqrt(float(t))
    return parents,children,w


def recorded_isometry(n,p=F(1,2)):
    parents,children,w = geometry_isometry(n,p)
    v = np.zeros((2*len(children),2*len(parents)))
    for i,h in enumerate(parents):
        for j,child in enumerate(children):
            flip = minima.minima_count(child)-minima.minima_count(h)
            if not w[j,i]:
                continue
            for r in range(2):
                v[2*j+(r^flip),2*i+r] = w[j,i]
    return parents,children,v


def recorded_state(n,p=F(1,2)):
    state = np.array([1.,0.])
    for j in range(2,n):
        state = recorded_isometry(j,p)[2]@state
    return state


def record_density(n,p=F(1,2)):
    x = recorded_state(n,p).reshape(-1,2)
    return x.T@x


def failure_of_sector_inference(n,p=F(1,2)):
    # R=1 certifies E3. On R=0, report whichever final sector is likelier.
    return min(p,(1-p)/(n-1))


def report():
    rows = []
    for n in range(2,6):
        law = online_level(n)
        keys,old = sector.comparator(n)
        error = max(abs(float(law.get(h,0))-v) for h,v in zip(keys,old))
        root_mass = sum((v for h,v in law.items() if minima.minima_count(h)==3),F(0))
        rows.append({'N':n,'supported_prefixes':len(law),
                     'normalization_exact':str(sum(law.values(),F(0))),
                     'difference_from_round246_classical_P':error,
                     'record_R1_probability_exact':str(root_mass),
                     'record_R0_probability_exact':str(1-root_mass),
                     'final_sector_inference_error_exact':str(failure_of_sector_inference(n))})
    _,_,v = recorded_isometry(4)
    return {'round':247,'date':'2026-09-22','prior_p':'1/2',
            'model':'Current-prefix hazard and uniform existing-pair closure; no initial final-sector register',
            'hazard_m2':'(1-p)/(n*(p*(n-1)+1-p))',
            'hazard_m3':'0',
            'p_half_hazard':'1/n^2',
            'finite_enumeration':rows,
            'N4_to_N5_recorded_isometry_error':float(np.max(abs(v.T@v-np.eye(v.shape[1])))),
            'terminal_no_third_minimum_probability':'p, by exact telescoping',
            'record_meaning':'R=1 iff a third minimal event has already appeared; R=0 is not an exact E2 certificate',
            'operational_history_D':'diag(P_N); differs from round246 rank-two D_rec',
            'extra_inputs':['Chosen p and time-dependent hazard','Uniform pair selection from the whole finite prefix',
                            'Orthogonal geometry memory and blank record qubits','Finite isometric control; no physical clock or spatial locality supplied'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_pair_closure_multiplicity_matches_lambda(self):
        for n in range(2,5):
            for parent in online_level(n):
                r = growth.matrix(parent)
                counts = defaultdict(int)
                for pair in itertools.combinations(range(n),2):
                    counts[pair_closure(r,pair)] += 1
                for s in csg.ideals(r):
                    if not s:
                        continue
                    k = csg.maximal_count(r,s)
                    coefficient = math.comb(len(s)-k,2-k) if k<=2<=len(s) else 0
                    self.assertEqual(counts[s],coefficient)

    def test_exact_online_probabilities(self):
        for p in [F(1,5),F(1,2),F(4,5)]:
            for n in range(2,6):
                law = online_level(n,p)
                self.assertEqual(sum(law.values(),F(0)),1)
                keys,old = sector.comparator(n,float(p))
                np.testing.assert_allclose([float(law.get(h,0)) for h in keys],old,atol=1e-13)

    def test_normalization_and_prefix_preservation(self):
        for n in range(2,5):
            out = defaultdict(F)
            for child,weight in online_level(n+1).items():
                parent = tuple(tuple(row[:-1]) for row in child[:-1])
                out[parent] += weight
            self.assertEqual(dict(out),online_level(n))

    def test_natural_relabelling_covariance(self):
        representatives = {}
        for h,weight in online_level(5).items():
            key = growth.canonical(h)
            if key in representatives:
                self.assertEqual(weight,representatives[key])
            representatives[key] = weight

    def test_both_isometries_on_full_supported_input_space(self):
        for n in range(2,5):
            for v in [geometry_isometry(n)[2],recorded_isometry(n)[2]]:
                np.testing.assert_allclose(v.T@v,np.eye(v.shape[1]),atol=1e-13)

    def test_record_is_exact_current_event_certificate(self):
        for n in range(2,6):
            state = recorded_state(n).reshape(-1,2)
            for i,(h,weight) in enumerate(online_level(n).items()):
                r = minima.minima_count(h)-2
                self.assertAlmostEqual(state[i,1-r],0)
                self.assertAlmostEqual(state[i,r]**2,float(weight))
            np.testing.assert_allclose(record_density(n),np.diag([.5+.5/(n-1),.5-.5/(n-1)]),atol=1e-13)

    def test_old_records_are_not_overwritten_after_third_minimum(self):
        for n in range(3,5):
            for h in online_level(n):
                if minima.minima_count(h)==3:
                    self.assertTrue(all(minima.minima_count(c)==3 for c in transitions(h)))

    def test_blank_copy_record_and_unknown_reference_preserved(self):
        state = recorded_state(4).reshape(-1,2)
        copied = np.zeros((len(state),2,2))
        for r in range(2):
            copied[:,r,r] = state[:,r]
        rho1 = np.einsum('irs,its->rt',copied,copied)
        rho2 = np.einsum('irs,irt->st',copied,copied)
        np.testing.assert_allclose(rho1,rho2,atol=1e-13)
        self.assertAlmostEqual(float(np.sum(copied[:,0,1]**2+copied[:,1,0]**2)),0)
        bell = np.array([1,0,0,1])/math.sqrt(2)
        joint = np.outer(copied.ravel(),bell)
        np.testing.assert_allclose(joint.T@joint,np.outer(bell,bell),atol=1e-13)

    def test_hazard_telescopes_without_initial_future_label(self):
        for p in [F(1,5),F(1,2),F(4,5)]:
            survival = F(1)
            for n in range(2,100):
                survival *= 1-hazard(n,2,p)
                self.assertEqual(survival,p+(1-p)/n)

    def test_pointer_histories_are_classical_and_not_D_rec(self):
        keys,prob = sector.comparator(4)
        d = np.diag(prob)
        active = [j for j,h in enumerate(keys) if h in online_level(4)]
        self.assertEqual(np.linalg.matrix_rank(d),len(active))
        self.assertGreater(np.max(abs(d-sector.functional(4)[1])),.01)
        self.assertAlmostEqual(float(d.sum()),1)


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
        target=Path(__file__).with_name('online_record_growth_audit_results.json')
        if target.exists() and target.read_text(encoding='utf-8')!=payload:
            raise SystemExit('Refusing to overwrite a different saved result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
