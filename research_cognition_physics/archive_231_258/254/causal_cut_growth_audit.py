"""Round 254: complete local-depth cuts remove bookkeeping scheduler bias.

Born record measures on marked port forests are normalized, projectively
consistent and independent of which enabled port is processed first. This is
not equal weighting of natural birth labels or an infinite coherent measure.
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
import local_branching_growth_audit as local

POLICIES=('uniform','left_bias','adaptive','lexicographic')


@lru_cache(maxsize=None)
def one_tree(address,depth):
    if depth==0: return ((),)
    out=[]
    for a in one_tree(address+(0,),depth-1):
        out.append(tuple(sorted(((address,-1),)+a)))
        for b in one_tree(address+(1,),depth-1):
            out.append(tuple(sorted(((address,1),)+a+b)))
    return tuple(out)


@lru_cache(maxsize=None)
def trees(depth,roots=local.ROOTS):
    return tuple(tuple(sorted(sum(parts,()))) for parts in
                 itertools.product(*(one_tree(r,depth) for r in roots)))


def enabled(tree,done):
    return tuple(a for a,g in tree if a not in done and (len(a)==1 or a[:-1] in done))


@lru_cache(maxsize=None)
def linear_orders(tree,done=()):
    if len(done)==len(tree): return (done,)
    return tuple(o for a in enabled(tree,done) for o in linear_orders(tree,done+(a,)))


def selection(eligible,records,policy):
    if policy=='uniform': weights={a:1 for a in eligible}
    elif policy=='left_bias': weights={a:2 if a[0]==0 else 1 for a in eligible}
    elif policy=='adaptive':
        weights={a:1+sum(g==1 and b[0]==a[0] for b,g in records) for a in eligible}
    elif policy=='lexicographic': weights={a:int(a==min(eligible)) for a in eligible}
    elif policy=='roots_first': weights={a:int(len(a)==min(map(len,eligible))) for a in eligible}
    else: raise ValueError(policy)
    total=sum(weights.values())
    return {a:F(w,total) for a,w in weights.items()}


def order_weight(tree,order,policy):
    outcomes=dict(tree); done=(); result=F(1)
    for a in order:
        probabilities=selection(enabled(tree,done),tuple((b,outcomes[b]) for b in done),policy)
        result*=probabilities[a]; done+=(a,)
    return result


def branch(tree,roots=local.ROOTS,order=None):
    order=tuple(a for a,g in tree) if order is None else order
    frontier=roots; matrix=np.eye(2**len(roots)); outcomes=dict(tree)
    for a in order: frontier,matrix=local.advance(frontier,matrix,a,outcomes[a])
    return frontier,matrix


@lru_cache(maxsize=None)
def intrinsic_effects(depth,roots=local.ROOTS):
    return {t:(h:=branch(t,roots)[1]).T@h for t in trees(depth,roots)}


def online_effects(depth,policy,roots=local.ROOTS):
    """Independently enumerate scheduled instruments on current live ports."""
    effects={}; trajectories=0
    def visit(frontier,matrix,records):
        nonlocal trajectories
        active=tuple(a for a in frontier if len(a)-1<depth)
        if not active:
            key=tuple(sorted(records)); effect=matrix.T@matrix
            effects[key]=effects.get(key,np.zeros_like(effect))+effect
            trajectories+=1; return
        for a,weight in selection(active,records,policy).items():
            if not weight: continue
            for g in (1,-1):
                f,h=local.advance(frontier,matrix,a,g)
                visit(f,math.sqrt(float(weight))*h,records+((a,g),))
    visit(roots,np.eye(2**len(roots)),())
    return effects,trajectories


def restrict(tree,depth): return tuple((a,g) for a,g in tree if len(a)-1<depth)


def probability(effect,rho): return float(np.trace(effect@rho).real)


def both_roots_after_two_steps(policy):
    psi=np.ones(4)/2; total=0.
    for a,w1 in selection(local.ROOTS,(),policy).items():
        if not w1: continue
        for g in (1,-1):
            f,h=local.advance(local.ROOTS,np.eye(4),a,g)
            for b,w2 in selection(f,((a,g),),policy).items():
                if not w2 or b not in local.ROOTS or b==a: continue
                for gg in (1,-1):
                    _,hh=local.advance(f,h,b,gg)
                    total+=float(w1*w2)*float(np.linalg.norm(hh@psi)**2)
    return total


def report():
    effects=intrinsic_effects(2); rho=np.zeros((4,4)); rho[0,0]=rho[-1,-1]=rho[0,-1]=rho[-1,0]=.5
    size=defaultdict(float)
    for t,e in effects.items(): size[len(t)]+=probability(e,rho)
    rows=[]
    for policy in POLICIES:
        scheduled,count=online_effects(2,policy)
        rows.append({'policy':policy,'positive_weight_transcripts':count,
                     'largest_effect_difference':max(float(np.max(abs(scheduled[t]-e))) for t,e in effects.items()),
                     'schedule_weight_sum_exact':str(min(sum((order_weight(t,o,policy) for o in linear_orders(t)),F(0)) for t in effects))})
    return {'round':254,'date':'2026-09-22','initial_roots':2,'local_depth_cut':2,
            'terminal_marked_forests':len(effects),
            'legal_transcripts':sum(len(linear_orders(t)) for t in effects),
            'normalization_operator_error':float(np.max(abs(sum(effects.values())-np.eye(4)))),
            'scheduler_comparison':rows,
            'event_count_distribution_for_logical_Bell_input':{str(k):v for k,v in sorted(size.items())},
            'fixed_two_global_steps_counterexample':{
                'uniform_both_roots_processed':both_roots_after_two_steps('uniform'),
                'roots_first_both_roots_processed':both_roots_after_two_steps('roots_first')},
            'extension':'Consistent finite-depth Born record probabilities define a classical probability measure on infinite marked port forests',
            'limits':['Scheduler only chooses enabled events using classical past; it does not alter settings or impose a stopping deadline',
                      'All events within the outcome-dependent depth cut are completed',
                      'No equal probabilities for different birth-number transcripts are assumed',
                      'Only classical pointer records are extended, not an infinite coherent history functional',
                      'Port types, fresh qubits, branch-only topology and local instruments remain additional inputs'],
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


class Audit(unittest.TestCase):
    def test_variable_tree_counts_and_complete_cut_instruments(self):
        self.assertEqual([len(one_tree((0,),d)) for d in range(4)],[1,2,6,42])
        for roots,maximum in (((0,),),3), (local.ROOTS,2):
            for depth in range(maximum+1):
                effects=intrinsic_effects(depth,roots)
                np.testing.assert_allclose(sum(effects.values()),np.eye(2**len(roots)),atol=1e-12)

    def test_every_topological_order_has_the_same_rectangular_operator(self):
        for tree in trees(2):
            frontier,expected=branch(tree)
            for order in linear_orders(tree):
                f,h=branch(tree,order=order)
                self.assertEqual(f,frontier); np.testing.assert_allclose(h,expected,atol=1e-13)

    def test_exact_scheduler_weights_sum_to_one_per_forest(self):
        for tree in trees(2):
            for policy in POLICIES:
                self.assertEqual(sum((order_weight(tree,o,policy) for o in linear_orders(tree)),F(0)),1)

    def test_independent_online_enumeration_all_policies(self):
        expected=intrinsic_effects(2)
        for policy in POLICIES:
            actual,count=online_effects(2,policy)
            self.assertEqual(set(actual),set(expected))
            for t in actual: np.testing.assert_allclose(actual[t],expected[t],atol=1e-13)

    def test_depth_two_to_one_operator_prefix_consistency(self):
        coarse=defaultdict(lambda:np.zeros((4,4)))
        for t,e in intrinsic_effects(2).items(): coarse[restrict(t,1)]+=e
        for t,e in intrinsic_effects(1).items(): np.testing.assert_allclose(coarse[t],e,atol=1e-13)

    def test_one_root_depth_three_to_two_operator_prefix_consistency(self):
        roots=((0,),); coarse=defaultdict(lambda:np.zeros((2,2)))
        for t,e in intrinsic_effects(3,roots).items(): coarse[restrict(t,2)]+=e
        for t,e in intrinsic_effects(2,roots).items(): np.testing.assert_allclose(coarse[t],e,atol=1e-13)

    def test_entangled_input_affects_forest_size_statistics(self):
        for sign in (1,-1):
            v=np.array([1,0,0,sign])/math.sqrt(2); rho=np.outer(v,v)
            size=defaultdict(float)
            for t,e in intrinsic_effects(2).items(): size[len(t)]+=probability(e,rho)
            expected={4:(1+sign*.36)/4,5:(1-sign*.36)/2,6:(1+sign*.36)/4}
            for k in expected: self.assertAlmostEqual(size[k],expected[k])

    def test_internal_phase_remains_visible_without_a_global_schedule(self):
        for sign in (1,-1):
            v=np.array([1,sign])/math.sqrt(2); rho=np.outer(v,v)
            p=sum(probability(e,rho) for t,e in intrinsic_effects(2,((0,),)).items() if dict(t)[(0,)]==1)
            self.assertAlmostEqual(p,(1+.6*sign)/2)

    def test_fixed_global_step_deadline_is_not_a_complete_cut(self):
        self.assertAlmostEqual(both_roots_after_two_steps('uniform'),11/30)
        self.assertAlmostEqual(both_roots_after_two_steps('roots_first'),1)

    def test_all_output_spaces_account_for_fresh_qubits(self):
        for tree in trees(2):
            frontier,h=branch(tree)
            splits=sum(g==1 for a,g in tree)
            self.assertEqual(len(frontier),len(local.ROOTS)+splits)
            self.assertEqual(h.shape,(2**len(frontier),4))
            self.assertTrue(all(len(a)==3 for a in frontier))
            self.assertTrue(all(len(a)==1 or a[:-1] in dict(tree) for a,g in tree))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true'); args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful(): raise SystemExit(1)
    result=report(); result['checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        target=Path(__file__).with_name('causal_cut_growth_audit_results.json')
        if target.exists() and target.read_text(encoding='utf-8')!=payload: raise SystemExit('Refusing to overwrite a different result')
        target.write_text(payload,encoding='utf-8')
    print(payload)
