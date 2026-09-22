"""Round234: identical copies of a classically scheduled ordinary channel.

The random schedule is a past classical record, not an emergent physical clock.
We contrast unrestricted merged laboratories with explicitly ordered memory
interfaces. Classical examples diagnose permissions; no new quantum derivation.
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

SWAP=(0,2,1,3)
FLIP_BOTH_SWAP=tuple(x^3 for x in SWAP)


def single_fixed_points(fa,fb,schedule):
    return [(a,b) for a,b in product(range(2),repeat=2)
            if (a,b)==((fb[b],0) if schedule else (0,fa[a]))]


def branch_weight(fa,fb,s,t):
    count=0
    for a,b in product(range(4),repeat=2):
        oa,ob=fa[a],fb[b]
        x=(ob%2,0) if s else (0,oa%2)
        y=(ob//2,0) if t else (0,oa//2)
        count+=((a,b)==(x[0]+2*y[0],x[1]+2*y[1]))
    return count


def weight(law,fa,fb):
    return float(sum(law[s,t]*branch_weight(fa,fb,s,t) for s,t in product(range(2),repeat=2)))


def independent_law(p):
    return np.outer([p,1-p],[p,1-p])


@lru_cache(None)
def permutation_certificate():
    controls=list(permutations(range(4)))
    hist=Counter(sum(branch_weight(a,b,s,t) for s,t in product(range(2),repeat=2))
                 for a,b in product(controls,repeat=2))
    return {'local_permutation_pairs':576,
            'fair_independent_weight_histogram':{str(k/4):hist[k] for k in sorted(hist)},
            'SWAP_SWAP_branch_weights':[[branch_weight(SWAP,SWAP,s,t) for t in range(2)] for s in range(2)],
            'SWAP_FLIP_branch_weights':[[branch_weight(SWAP,FLIP_BOTH_SWAP,s,t) for t in range(2)] for s in range(2)]}


def scheduled_forward(kernels,s,t):
    """Four single events, two local one-bit memories, all output records kept.

    kernels[agent,copy][o,new_memory,input,old_memory] is stochastic.
    Early events run first, then the two late events. Copy index is a protocol
    label; every schedule branch is an explicitly given DAG.
    """
    first=(s,t)  # 0 means A first, 1 means B first
    order=[(first[k],k,True) for k in range(2)]+[(1-first[k],k,False) for k in range(2)]
    # state = (memory A, memory B, pending wire0, pending wire1, record tuple)
    distribution={(0,0,0,0,()):1.}
    for agent,copy,early in order:
        updated={}
        kernel=kernels[agent,copy]
        for state,prob in distribution.items():
            memory=list(state[:2]);wires=list(state[2:4]);records=state[4]
            incoming=0 if early else wires[copy]
            for out,new_memory in product(range(2),repeat=2):
                p=prob*kernel[out,new_memory,incoming,memory[agent]]
                next_mem=memory.copy();next_mem[agent]=new_memory
                next_wire=wires.copy()
                if early:next_wire[copy]=out
                key=tuple(next_mem+next_wire)+(records+(out,),)
                updated[key]=updated.get(key,0.)+p
        distribution=updated
    record_dist=np.zeros((2,)*4)
    for state,p in distribution.items():record_dist[state[4]]+=p
    return record_dist


@lru_cache(None)
def memory_certificate():
    rng=np.random.default_rng(234)
    norm=early_change=0.
    for _ in range(32):
        kernels=rng.random((2,2,2,2,2,2))
        kernels/=kernels.sum(axis=(2,3),keepdims=True)
        for s,t in product(range(2),repeat=2):
            original=scheduled_forward(kernels,s,t)
            changed=kernels.copy()
            # Change only the late event's operation in each copy.
            for copy,first in enumerate((s,t)):
                agent=1-first
                changed[agent,copy]=0.
                changed[agent,copy,1,0,:,:]=1.
            perturbed=scheduled_forward(changed,s,t)
            norm=max(norm,abs(float(original.sum())-1),abs(float(perturbed.sum())-1))
            early_change=max(early_change,float(np.max(np.abs(original.sum(axis=(2,3))-perturbed.sum(axis=(2,3))))))
    return {'random_memory_protocols':32,'schedule_branches_per_protocol':4,
            'explicit_local_memory_bits':2,'schedule_record_bits':2,
            'maximum_normalization_error':norm,'maximum_early_record_change_after_late_intervention':early_change}


class ScheduledCopyTests(unittest.TestCase):
    def test_each_schedule_and_all_local_single_event_functions_are_valid(self):
        funcs=list(product(range(2),repeat=2))
        for s,fa,fb in product(range(2),funcs,funcs):
            self.assertEqual(len(single_fixed_points(fa,fb,s)),1)

    def test_independent_identical_mixtures_fail_unrestricted_grouping(self):
        for p in np.linspace(0,1,21):
            q=2*p*(1-p)
            self.assertAlmostEqual(weight(independent_law(p),SWAP,SWAP),1+q)
            self.assertAlmostEqual(weight(independent_law(p),SWAP,FLIP_BOTH_SWAP),1-q)
        self.assertEqual(weight(independent_law(.5),SWAP,SWAP),1.5)
        self.assertEqual(weight(independent_law(.5),SWAP,FLIP_BOTH_SWAP),.5)

    def test_all_local_permutations_at_fair_independent_schedule(self):
        self.assertEqual(permutation_certificate()['fair_independent_weight_histogram'],
                         {'0.5':36,'0.75':112,'1.0':280,'1.25':112,'1.5':36})

    def test_general_joint_schedule_law_has_defect_equal_to_disagreement(self):
        rng=np.random.default_rng(235)
        for _ in range(64):
            law=rng.random((2,2));law/=law.sum()
            q=law[0,1]+law[1,0]
            self.assertAlmostEqual(weight(law,SWAP,SWAP),1+q)
            self.assertAlmostEqual(weight(law,SWAP,FLIP_BOTH_SWAP),1-q)

    def test_shared_schedule_is_a_different_correlated_resource(self):
        controls=list(permutations(range(4)))
        for p in (.2,.5,.8):
            correlated=np.diag([p,1-p]);independent=independent_law(p)
            np.testing.assert_allclose(correlated.sum(axis=0),independent.sum(axis=0))
            self.assertAlmostEqual(float(np.abs(correlated-independent).sum()/2),2*p*(1-p))
            for fa,fb in product(controls,repeat=2):
                self.assertAlmostEqual(weight(correlated,fa,fb),1.)

    def test_prefilter_keeps_success_probability_and_failure_record(self):
        for p in (.2,.5,.8):
            law=independent_law(p);success=float(np.trace(law))
            filtered=np.diag(np.diag(law))
            # A physically earlier flag filter: retained success branches have
            # their original subnormalization, and failure remains explicit.
            self.assertAlmostEqual(weight(filtered,SWAP,SWAP),success)
            self.assertAlmostEqual(weight(filtered,SWAP,SWAP)+(1-success),1.)
            self.assertAlmostEqual(weight(filtered/success,SWAP,SWAP),1.)
            self.assertLess(success,1.)

    def test_event_typed_memory_protocols_preserve_normalization_and_past_records(self):
        c=memory_certificate()
        self.assertLess(c['maximum_normalization_error'],3e-15)
        self.assertLess(c['maximum_early_record_change_after_late_intervention'],3e-15)


def report():
    return {'round':234,'date':'2026-09-22',
            'hypothesis':'Does identical independent repeatability imply unrestricted joint control after erasing event slots?',
            'single_use_deterministic_settings':32,
            'identical_independent_copies':permutation_certificate(),
            'analytic_formula':{'disagreement_probability':'q=P(s!=t)',
                                'double_swap_normalization':'1+q',
                                'swap_and_flip_normalization':'1-q',
                                'independent_identical_order_law':'q=2p(1-p)'},
            'fair_case':{'independent_normalization_plus':1.5,'independent_normalization_minus':.5,
                         'shared_schedule_normalization':1.,'prior_flag_filter_success_probability':.5,
                         'correlated_vs_independent_schedule_TV':.5},
            'event_typed_positive_control':memory_certificate(),
            'conclusion':'A causally scheduled ordinary process may fail the untyped identical-copy contract. Keeping event slots permits normalized independent repetition; synchronizing schedules changes the resource.',
            'not_claimed':['Independent copies of quantum channels are impossible',
                           'A shared physical clock is universally necessary',
                           'Time emerges from normalization alone',
                           'The model selects Lorentz symmetry or geometry'],
            'next':'Define event interfaces by operational access and substitution tests, and determine which timing/type information those tests necessarily require.',
            'runtime':{'python':platform.python_version(),'numpy':np.__version__}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ScheduledCopyTests))
    if not result.wasSuccessful():raise SystemExit(1)
    data=report();data['checks']={'run':result.testsRun,'failures':0,'errors':0}
    if args.write_results:
        target=Path(__file__).with_name('scheduled_copy_audit_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8'))!=json.loads(json.dumps(data)):
            raise RuntimeError('Existing research result differs; inspect before replacing.')
        target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,ensure_ascii=False,indent=2))
