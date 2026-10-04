"""Round 251: a complete parity-controlled growth instrument.

The geometry marginal is exactly an old classical mixture. The marked seed,
global controller and its coupling are inputs, not causal or gravitational laws.
"""
import argparse
from functools import lru_cache
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import online_record_growth_audit as online
import graph_record_growth_audit as graph
import random_order_growth_audit as growth
import covariant_minima_audit as minima

PARAMS=(F(1,5),F(4,5))
S=graph.pauli(2,z=(0,1))
PROJECTORS=((np.eye(4)+S)/2,(np.eye(4)-S)/2)


def root_filter(weights):
    return sum(math.sqrt(float(w))*p for w,p in zip(weights,PROJECTORS))


def transition_maps(parent):
    laws=[online.transitions(parent,p) for p in PARAMS]
    n=len(parent)
    return {c:graph.append_map(c)@np.kron(root_filter([law[c] for law in laws]),np.eye(2**(n-2)))
            for c in laws[0]}


@lru_cache(maxsize=None)
def histories(n):
    if n==2:
        seed=next(iter(online.online_level(2)))
        return {seed:np.eye(4)}
    return {c:k@a for h,a in histories(n-1).items() for c,k in transition_maps(h).items()}


def closed_map(h):
    n=len(h)
    weights=[online.online_level(n,p)[h] for p in PARAMS]
    return graph.phase(n,graph.covers(h))[:,None]*np.kron(
        root_filter(weights),np.ones((2**(n-2),1))/math.sqrt(2**(n-2)))


def geometry_law(n,rho):
    return {h:float(np.trace(a@rho@a.conj().T).real) for h,a in histories(n).items()}


def decoded_parent_channel(n,rho):
    # Conditional inverse CZ, then discard geometry and blank newborns.
    return sum(root_filter([online.online_level(n,p)[h] for p in PARAMS])@rho@
               root_filter([online.online_level(n,p)[h] for p in PARAMS]) for h in histories(n))


def coherence_factor(n):
    return sum(math.sqrt(float(online.online_level(n,PARAMS[0])[h]*
                               online.online_level(n,PARAMS[1])[h])) for h in histories(n))


def analytic_factor(n):
    a,b=[float(p+(1-p)/F(n-1)) for p in PARAMS]
    return math.sqrt(a*b)+math.sqrt((1-a)*(1-b))


def density(v): return np.outer(v,np.conj(v))


def report():
    rows=[]
    for n in range(2,6):
        norms=sum(a.conj().T@a for a in histories(n).values())
        rows.append({'N':n,'geometries':len(histories(n)),
                     'instrument_completeness_error':float(np.max(abs(norms-np.eye(4)))),
                     'third_minimum_probabilities_by_parity':[
                         float((1-p)*F(n-2,n-1)) for p in PARAMS],
                     'decoded_cross_parity_coherence':coherence_factor(n),
                     'closed_coherence_formula':analytic_factor(n)})
    return {'round':251,'date':'2026-09-22','sector_parameters':[str(p) for p in PARAMS],
            'finite_growth':rows,'N3_new_minimum_probabilities':[.4,.1],
            'geometry_marginal':'P_rho(c)=w_plus P_(1/5)(c)+w_minus P_(4/5)(c)=P_p_eff(c)',
            'infinite_decoded_coherence_limit':.8,
            'limits':['Marked initial seed and global growth controller are inputs',
                      'Geometry statistics are exactly a classical mixture; no nonclassical geometry witness',
                      'State-dependent birth of a new minimum is not a local causal implementation',
                      'Preserved information within a parity sector requires access to the conditional graph decoder'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_each_birth_instrument_complete(self):
        for n in range(2,5):
            for h in histories(n):
                ks=transition_maps(h)
                np.testing.assert_allclose(sum(k.conj().T@k for k in ks.values()),np.eye(2**n),atol=1e-13)

    def test_recursive_history_equals_closed_map(self):
        for n in range(2,6):
            for h,a in histories(n).items():
                np.testing.assert_allclose(a,closed_map(h),atol=1e-13)

    def test_geometry_law_is_old_affine_family(self):
        for w in (F(0),F(1,3),F(1)):
            psi=np.array([math.sqrt(float(w)),1j*math.sqrt(float(1-w)),0,0])
            p=w*PARAMS[0]+(1-w)*PARAMS[1]
            for n in (3,4,5):
                law=geometry_law(n,density(psi))
                for h,prob in online.online_level(n,p).items():
                    self.assertAlmostEqual(law[h],float(prob))

    def test_geometry_cannot_detect_cross_sector_phase(self):
        for phase in (0,.7,2.3):
            psi=np.array([1,np.exp(1j*phase),0,0])/math.sqrt(2)
            rho=density(psi); decohered=sum(p@rho@p for p in PROJECTORS)
            a=geometry_law(5,rho); b=geometry_law(5,decohered)
            np.testing.assert_allclose(list(a.values()),list(b.values()),atol=1e-13)

    def test_full_geometry_direct_sum_isometry(self):
        n=3; parents=tuple(histories(n)); children=tuple(histories(n+1))
        v=np.zeros((len(children)*2**(n+1),len(parents)*2**n))
        for i,h in enumerate(parents):
            for c,k in transition_maps(h).items():
                j=children.index(c)
                v[j*2**(n+1):(j+1)*2**(n+1),i*2**n:(i+1)*2**n]=k
        np.testing.assert_allclose(v.T@v,np.eye(v.shape[1]),atol=1e-13)

    def test_operator_prefix_marginals(self):
        for n in range(2,5):
            for h,a in histories(n).items():
                total=sum((k@a).conj().T@(k@a) for k in transition_maps(h).values())
                np.testing.assert_allclose(total,a.conj().T@a,atol=1e-13)

    def test_marked_seed_relabelling(self):
        for h,a in histories(5).items():
            for seed in ((0,1),(1,0)):
                for tail in itertools.permutations(range(2,5)):
                    perm=seed+tail; r=growth.matrix(h)[np.ix_(perm,perm)]
                    if np.any(np.tril(r)): continue
                    hp=growth.key(r)
                    transformed=a.reshape((2,)*7).transpose(perm+tuple(5+i for i in seed)).reshape(32,4)
                    np.testing.assert_allclose(transformed,histories(5)[hp],atol=1e-13)

    def test_protected_parity_and_old_snapshot(self):
        for n in (3,4,5):
            for h,a in histories(n).items():
                np.testing.assert_allclose(graph.pauli(n,z=(0,1))@a,a@S,atol=1e-13)
        rho=np.eye(4)/4
        old=geometry_law(3,rho)
        for n in (4,5):
            recovered={h:0. for h in old}
            for h,p in geometry_law(n,rho).items():
                prefix=tuple(tuple(row[:3]) for row in h[:3]); recovered[prefix]+=p
            np.testing.assert_allclose(list(recovered.values()),list(old.values()),atol=1e-13)

    def test_coherence_loss_formula_after_decoding(self):
        rho=density(np.array([1,1j,0,0])/math.sqrt(2))
        for n in range(2,6):
            self.assertAlmostEqual(coherence_factor(n),analytic_factor(n))
            result=decoded_parent_channel(n,rho)
            self.assertAlmostEqual(abs(result[0,1]),analytic_factor(n)/2)

    def test_unknown_sector_state_and_reference_preserved_after_decode(self):
        for s,indices in enumerate(((0,3),(1,2))):
            psi=np.zeros(8,complex); psi[2*indices[0]]=1/math.sqrt(2); psi[2*indices[1]+1]=1j/math.sqrt(2)
            rho=density(psi)
            for n in (3,5):
                out=np.zeros_like(rho)
                for h in histories(n):
                    b=np.kron(root_filter([online.online_level(n,p)[h] for p in PARAMS]),np.eye(2))
                    out+=b@rho@b.conj().T
                np.testing.assert_allclose(out,rho,atol=1e-13)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true'); args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful(): raise SystemExit(1)
    result=report(); result['checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        target=Path(__file__).with_name('quantum_controlled_birth_audit_results.json')
        if target.exists() and target.read_text(encoding='utf-8')!=payload:
            raise SystemExit('Refusing to overwrite a different saved result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
