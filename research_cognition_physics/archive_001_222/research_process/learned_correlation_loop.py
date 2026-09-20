"""Round 114: a stopped finite learning-query loop and the independence boundary."""
import argparse
import json
import math
import unittest
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, ceil_div
from biased_source_compiler import biased_environment_messages, biased_budget
from fragment_query_compiler import query_branches
from incompatible_relation_writes import (H, calibration_state, written_state,
    forget_memories, old_channel)
from joint_relation_records import interval_fields
from record_correlation_frontier import biased_probe_records
from unknown_correlation_learning import (choose_sign, calibration_visibility,
    source_mean_interval, uniform_risk_certificate, matched_policy_value)


def learned_query_records(old, true_kappa, learned_sign, m=1, q=1):
    """The physical source uses true_kappa; controller receives only learned_sign."""
    if learned_sign not in (-1,0,1):
        raise ValueError("Invalid learned action.")
    c=d=.8
    gates=(math.sqrt(c),)*3
    out={}
    for message,branch in biased_environment_messages(
            written_state(old,c,d),c,(.6,)*3,true_kappa,gates,m).items():
        action=("joint",1) if learned_sign and message==-learned_sign else ("second_only",1)
        for history,value in query_branches(branch,c,d,q,*action).items():
            out[(message,history)]=value
    return out


def stopped_distribution(probability, schedule, action):
    """Directed integer probability propagation, retaining every active count."""
    schedule=tuple(schedule)
    if not schedule or any(not isinstance(n,int) or n<1 for n in schedule) or any(
            x>=y for x,y in zip(schedule,schedule[1:])):
        raise ValueError("Use strictly increasing positive integer looks.")
    p=I.exact(probability)
    if p.lo<0 or p.hi>SCALE:
        raise ValueError("Invalid Bernoulli interval.")
    q=1-p
    lower=[SCALE];upper=[SCALE]
    stopped=[]
    look_of={n:j+1 for j,n in enumerate(schedule)}
    for n in range(1,schedule[-1]+1):
        lo=[0]*(n+1);hi=[0]*(n+1)
        for k in range(n):
            lo[k]+=lower[k]*q.lo//SCALE
            hi[k]+=ceil_div(upper[k]*q.hi,SCALE)
            lo[k+1]+=lower[k]*p.lo//SCALE
            hi[k+1]+=ceil_div(upper[k]*p.hi,SCALE)
        lower,upper=lo,hi
        if n in look_of:
            masses={1:[0,0],-1:[0,0]}
            j=look_of[n]
            for k in range(n+1):
                decision=action(n,k,j)
                if decision not in (-1,0,1):
                    raise ValueError("Invalid action.")
                if decision:
                    masses[decision][0]+=lower[k]
                    masses[decision][1]+=upper[k]
                    lower[k]=upper[k]=0
            stopped.append((n,I(*masses[1]),I(*masses[-1])))
    timeout=I(sum(lower),sum(upper))
    plus=sum((p for n,p,m in stopped),I.exact(0))
    minus=sum((m for n,p,m in stopped),I.exact(0))
    expected=sum((n*(p+m) for n,p,m in stopped),I.exact(0))+schedule[-1]*timeout
    return {"stopped":stopped,"plus":plus,"minus":minus,"timeout":timeout,"expected_samples":expected}


@lru_cache(maxsize=None)
def action_table(n,look,a,target):
    return tuple(choose_sign(n,k,look,F(64,125),a,target) for k in range(n+1))


@lru_cache(maxsize=None)
def loop_distribution(kappa, a=F(1), target=F(1,5), schedule=(64,128,256,512,1024)):
    kappa,a,target=map(F,(kappa,a,target))
    if abs(kappa)>F(64,125):
        raise ValueError("Kappa outside the declared source family.")
    p=(1+source_mean_interval(kappa,a))/2
    return stopped_distribution(p,schedule,lambda n,k,j:action_table(n,j,a,target)[k])


def loop_certificate(kappa, a=F(1), target=F(1,5), schedule=(64,128,256,512,1024)):
    kappa,a,target=map(F,(kappa,a,target))
    result=loop_distribution(kappa,a,target,schedule)
    base,coefficient=matched_policy_value(0,0)
    gain=coefficient*I.exact(kappa)*(result["plus"]-result["minus"])
    wrong=result["minus"] if kappa>0 else result["plus"] if kappa<0 else result["plus"]+result["minus"]
    ncap=schedule[-1]
    query=biased_budget(1,1)
    expected=result["expected_samples"]
    return {"true_kappa_for_evaluation_only_exact":str(kappa),"calibration_response_strength_exact":str(a),
        "activation_margin_exact":str(target),"prespecified_sample_looks":list(schedule),
        "positive_policy_probability":interval_fields(result["plus"]),
        "negative_policy_probability":interval_fields(result["minus"]),
        "timeout_fallback_probability":interval_fields(result["timeout"]),
        "wrong_sign_or_null_activation_probability":interval_fields(wrong),
        "expected_independent_samples":interval_fields(expected),
        "learned_future_H_success":interval_fields(base+gain),
        "fallback_H_success":interval_fields(base),
        "unconditional_H_success_gain":interval_fields(gain),
        "stop_distribution":[{"samples":n,"positive":interval_fields(p),"negative":interval_fields(m)}
                             for n,p,m in result["stopped"]],
        "resources":{
            "per_calibration_sample_pair_rotations":34,"per_calibration_sample_extra_local_rotations":41,
            "per_calibration_sample_pure_slots":7,
            "expected_total_pair_rotations_with_one_query":interval_fields(34*expected+query["pair_rotations"]),
            "expected_total_extra_local_rotations_with_one_query":interval_fields(41*expected+query["extra_local_rotations"]),
            "expected_total_fresh_pure_slots_allocated":interval_fields(7*expected+query["initial_pure_preparations"]),
            "maximum_total_pair_rotations":34*ncap+query["pair_rotations"],
            "maximum_total_extra_local_rotations":41*ncap+query["extra_local_rotations"],
            "maximum_total_fresh_pure_slots_allocated_and_retained":7*ncap+query["initial_pure_preparations"],
            "maximum_actual_readouts":ncap+2,
            "maximum_distributed_calibration_record_bits":ncap,
            "controller_to_query_policy_bits":2,"per_query_environment_return_bits":1,
            "logical_classical_sufficient_statistic":"current count n, positive count k, look j, chosen sign",
            "all_past_quantum_slots_retained_in_closed_whole":True,
            "calibration_and_future_query_use_independent_source_copies":True,
            "transport_and_source_refresh_physical_cost_modelled":False},
        "uniform_confidence_failure_bound":interval_fields(uniform_risk_certificate())}


def protected_variance(kappa, a, count):
    kappa,a=map(F,(kappa,a))
    if not isinstance(count,int) or count<1 or not -1<=kappa<=1 or not 0<=a<=1:
        raise ValueError("Invalid protected block parameters.")
    alpha=calibration_visibility(1)
    theta=I.exact(a*kappa)
    return alpha**2*(1-theta**2)+(1-alpha**2)/count


def finite_sample_activation_ceiling(kappa, count, null_activation_bound):
    # Even arbitrary source-independent quantum controls cannot increase this
    # tensor-product source trace-distance bound; capped stopping can be padded.
    kappa,delta=map(F,(kappa,null_activation_bound))
    if not isinstance(count,int) or count<1 or not 0<=delta<=1 or not -1<=kappa<=1:
        raise ValueError("Invalid bound parameters.")
    return min(F(1),delta+count*abs(kappa)/2)


class LearnedCorrelationLoopTests(unittest.TestCase):
    def test_wrong_and_right_learned_policies_match_complete_original_noisy_trees(self):
        for kappa in (-F(64,125),F(1,5)):
            for learned in (-1,0,1):
                success=0.
                for label in (-1,1):
                    records=learned_query_records(calibration_state(H,label),float(kappa),learned)
                    success+=sum(np.trace(rho).real for guess,rho in records.values() if guess==label)/2
                expected,_=matched_policy_value(kappa,learned)
                self.assertAlmostEqual(success,sum(expected.floats())/2,places=13)

    def test_any_learned_classical_policy_preserves_old_unconditional_reference_channel(self):
        from complex_whole_real_interfaces import random_state
        rho=random_state(np.random.default_rng(114),8)
        for learned in (-1,0,1):
            records=learned_query_records(rho,.512,learned)
            actual=sum(forget_memories(branch) for guess,branch in records.values())
            np.testing.assert_allclose(actual,old_channel(rho,.8,.8),atol=2e-15)

    def test_directed_stopping_recursion_contains_exact_exhaustive_small_tree(self):
        p=F(3,5)
        action=lambda n,k,j:1 if k==n else -1 if k==0 else 0
        output=stopped_distribution(I.exact(p),(2,4),action)
        exact={1:F(0),-1:F(0),0:F(0)};expected=F(0)
        for word in product((0,1),repeat=4):
            mass=p**sum(word)*(1-p)**(4-sum(word))
            sign=0
            for n in (2,4):
                sign=action(n,sum(word[:n]),1)
                if sign: break
            exact[sign]+=mass;expected+=n*mass
        for key,field in ((1,"plus"),(-1,"minus"),(0,"timeout")):
            self.assertLessEqual(F(output[field].lo,SCALE),exact[key])
            self.assertGreaterEqual(F(output[field].hi,SCALE),exact[key])
        self.assertLessEqual(F(output["expected_samples"].lo,SCALE),expected)
        self.assertGreaterEqual(F(output["expected_samples"].hi,SCALE),expected)

    def test_full_stopping_probabilities_normalize_and_learned_gain_is_certified(self):
        row=loop_certificate(F(64,125))
        dist=loop_distribution(F(64,125))
        total=dist["plus"]+dist["minus"]+dist["timeout"]
        self.assertLessEqual(total.lo,SCALE)
        self.assertGreaterEqual(total.hi,SCALE)
        self.assertLess(total.hi-total.lo,10**12)
        self.assertGreater(F(row["unconditional_H_success_gain"]["lower_exact"]),F(15,1000))
        self.assertLess(F(row["expected_independent_samples"]["upper_exact"]),300)

    def test_sign_reflection_and_null_source_have_consistent_decisions(self):
        positive=loop_distribution(F(64,125))
        negative=loop_distribution(-F(64,125))
        self.assertAlmostEqual(positive["plus"].floats()[0],negative["minus"].floats()[0],places=14)
        null=loop_distribution(F(0))
        self.assertLess((null["plus"]+null["minus"]).hi,uniform_risk_certificate().hi)
        self.assertAlmostEqual(null["plus"].floats()[0],null["minus"].floats()[0],places=14)
        row=loop_certificate(F(0))
        self.assertEqual(F(row["unconditional_H_success_gain"]["lower_exact"]),0)

    def test_protected_repeats_share_a_latent_record_and_have_a_variance_floor(self):
        kappa=F(64,125)
        records=biased_probe_records(.8,(.6,)*3,float(kappa),(.9,0.,0.),3)
        mean=second=pair=0.
        for history,rho in records.items():
            weight=np.trace(rho).real
            # Output Y reverses the original Z sign: bit 1 -> -1.
            ys=[1-2*bit for bit in history]
            average=sum(ys)/3
            mean+=weight*average;second+=weight*average**2;pair+=weight*ys[0]*ys[1]
        variance=second-mean**2
        expected=protected_variance(kappa,1,3)
        self.assertAlmostEqual(variance,sum(expected.floats())/2,places=13)
        alpha=calibration_visibility(1).floats()[0]
        self.assertAlmostEqual(pair,alpha**2,places=13)
        self.assertGreater(variance,(1-mean**2)/3+.4)

    def test_finite_cap_cannot_guarantee_detection_of_every_nonzero_correlation(self):
        ceiling=finite_sample_activation_ceiling(F(1,10**6),1024,F(1,100))
        self.assertEqual(ceiling,F(657,62500))
        self.assertLess(ceiling,F(11,1000))
        # Even revealing the single underlying protected bit leaves overlap.
        ideal_single_source_tv=F(64,125)/2
        self.assertLess(ideal_single_source_tv,1)

    def test_budget_counts_every_fresh_purification_and_actual_return_gate(self):
        report=loop_certificate(F(64,125))
        b=report["resources"]
        self.assertEqual(b["maximum_total_fresh_pure_slots_allocated_and_retained"],7176)
        self.assertEqual(b["maximum_total_pair_rotations"],34854)
        self.assertEqual(b["maximum_actual_readouts"],1026)
        self.assertTrue(b["all_past_quantum_slots_retained_in_closed_whole"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(LearnedCorrelationLoopTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":114,
        "controller":"At looks 64,128,256,512,1024 certify kappa>=1/5 or kappa<=-1/5, otherwise continue; timeout uses N only.",
        "examples":[loop_certificate(k) for k in (F(64,125),F(1,5),F(0),-F(64,125))],
        "protected_three_read_variance":interval_fields(protected_variance(F(64,125),1,3)),
        "finite_cap_obstruction":"P_kappa(activate)<=P_0(activate)+N*abs(kappa)/2 for iid family copies and arbitrary source-independent processing",
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("learned_correlation_loop_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
