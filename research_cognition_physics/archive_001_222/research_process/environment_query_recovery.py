"""Round 102: finite noisy environment access restores an erased query conditionally."""

import argparse
import json
import math
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from deferred_relation_queries import memory_channel,h_query_decoder
from incompatible_relation_writes import G,H,calibration_state,written_state,forget_memories,old_channel
from joint_relation_records import (memory_axis,axis_read_memory_tree,read_memory_tree,interval_fields)
from memory_dephasing_threshold import dephasing_unitary,dilate_memory,dephase_memory
from one_bit_real_network import visibility_interval
from quantum_interface_audit import rotation_unitary
from weak_relation_tradeoff import forget_front_pointer,trace_distance


def environment_record_branches(state,c,lam,count):
    """Only the majority sign is exported; the quantum environment is then discarded.

    Environment basis F=lam*Z+sqrt(1-lam^2)*X differs from the G-record basis.
    The original readout noise remains present.
    """
    extended=dilate_memory(state,c,lam)
    raw=axis_read_memory_tree(extended,0,math.acos(lam),count)
    result={r:np.zeros_like(state,dtype=complex) for r in (-1,1)}
    for history,branch in raw.items():
        r=1 if sum(history)>count//2 else -1
        result[r]+=forget_front_pointer(branch)
    return result


def branch_moments(lam,gamma,r):
    return (1+r*lam*gamma)/2,(lam+r*gamma)/2


def analytic_environment_branch(state,c,lam,gamma,r):
    a=np.kron(memory_axis(c),np.eye(len(state)//2))
    plus=(1+lam)*(1+r*gamma)/4
    minus=(1-lam)*(1-r*gamma)/4
    return plus*state+minus*a@state@a


def conditional_query_strategy(c,lam,gamma,r):
    p,a=branch_moments(lam,gamma,r)
    return ("joint",1 if a>=0 else -1) if abs(a)>c*p else ("second_only",1)


def recovery_factor(c,lam,gamma):
    return sum(max(c*p,abs(a)) for r in (-1,1)
               for p,a in [branch_moments(lam,gamma,r)])


def recovery_records(state,c,d,lam,environment_count,query_count):
    gamma=sum(visibility_interval(environment_count).floats())/2
    output={}
    for r,branch in environment_record_branches(written_state(state,c,d),c,lam,environment_count).items():
        strategy,sign=conditional_query_strategy(c,lam,gamma,r)
        if strategy=="joint":
            u=np.kron(h_query_decoder(c,d),np.eye(len(state)))
            leaves=axis_read_memory_tree(u@branch@u.conj().T,0,0.,query_count)
        else:
            leaves=read_memory_tree(branch,1,d,query_count)
        for history,final in leaves.items():
            guess=sign*(1 if sum(history)>query_count//2 else -1)
            output[(r,history)]=(guess,final)
    return output


def interval_abs(x):
    if x.lo>=0: return x
    if x.hi<=0: return -x
    return I(0,max(-x.lo,x.hi))


def interval_max(x,y):
    return I(max(x.lo,y.lo),max(x.hi,y.hi))


def recovery_certificate(c,d,lam,environment_count,query_count):
    c,d,lam=map(Fraction,(c,d,lam))
    if any(not 0<=x<=1 for x in (c,d,lam)): raise ValueError("Use parameters in [0,1].")
    gamma=visibility_interval(environment_count)
    eta=visibility_interval(query_count)
    total=I.exact(0)
    choices=[]
    for r in (-1,1):
        p=(1+r*I.exact(lam)*gamma)/2
        a=(I.exact(lam)+r*gamma)/2
        local=I.exact(c)*p
        joint=interval_abs(a)
        total+=interval_max(local,joint)
        if joint.lo>local.hi:
            decision="joint_positive" if a.lo>=0 else "joint_negative"
        elif local.lo>joint.hi:
            decision="second_only"
        else: decision="tie_or_interval_overlap"
        choices.append({"environment_sign":r,"strategy_certificate":decision,
                        "probability":interval_fields(p)})
    w=(1-I.exact(d*d)).sqrt()
    baseline=I.exact(max(c,lam))
    return {"c_exact":str(c),"d_exact":str(d),"retained_coherence_exact":str(lam),
        "environment_readouts":environment_count,"final_H_readouts":query_count,
        "environment_classical_message_bits":1,
        "remaining_memory_ideal_recovery_factor":interval_fields(total),
        "finite_H_success":interval_fields((1+eta*w*total)/2),
        "without_environment_success":interval_fields((1+eta*w*baseline)/2),
        "strict_gain":interval_fields(eta*w*(total-baseline)/2),
        "choices":choices,"all_environment_histories_optimized":False,
        "optimality_scope":"Fixed two-valued environment message; ideal final memory discrimination",
        "full_quantum_memory_state_restored_by_this_record_only_protocol":False}


class EnvironmentQueryRecoveryTests(unittest.TestCase):
    def test_ideal_eraser_basis_produces_identity_or_G_axis_phase_flip(self):
        rho=random_state(np.random.default_rng(102),2)
        for c,lam in ((0,0),(.8,0),(.8,.7),(.8,1)):
            v=dephasing_unitary(c,lam)@np.kron(np.eye(2)[:,:1],np.eye(2))
            basis=rotation_unitary(math.acos(lam))
            for r,column in ((1,0),(-1,1)):
                k=np.kron(basis[:,column].conj()[None,:],np.eye(2))@v
                op=np.eye(2) if r==1 else memory_axis(c)
                np.testing.assert_allclose(k@rho@k.conj().T,
                    (1+r*lam)*op@rho@op.conj().T/2,atol=5e-16)

    def test_actual_noisy_environment_record_branches_equal_weighted_unitary_channels(self):
        rng=np.random.default_rng(202)
        for size in (4,8):
            state=random_state(rng,size)
            for lam in (0,.7,.95,1):
                for count in (1,3):
                    gamma=sum(visibility_interval(count).floats())/2
                    branches=environment_record_branches(state,.8,lam,count)
                    for r,branch in branches.items():
                        np.testing.assert_allclose(branch,
                            analytic_environment_branch(state,.8,lam,gamma,r),atol=8e-16)
                        self.assertAlmostEqual(np.trace(branch).real,branch_moments(lam,gamma,r)[0],places=13)
                    np.testing.assert_allclose(sum(branches.values()),dephase_memory(state,.8,lam),atol=1e-15)

    def test_all_message_blocks_have_the_claimed_optimal_H_trace_distance(self):
        for c,d,lam in ((0,0,0),(.8,.8,0),(.8,.8,.95),(.3,.7,.6),(1,.4,.2)):
            for count in (1,3):
                gamma=sum(visibility_interval(count).floats())/2
                minus,plus=[environment_record_branches(memory_channel(calibration_state(H,s),c,d),c,lam,count)
                            for s in (-1,1)]
                actual=sum(trace_distance(minus[r],plus[r]) for r in (-1,1))
                self.assertAlmostEqual(actual,math.sqrt(1-d*d)*recovery_factor(c,lam,gamma),places=13)

    def test_actual_adaptive_memory_queries_attain_finite_certificates_without_postselection(self):
        for lam in (Fraction(0),Fraction(1,2),Fraction(19,20)):
            for ecount,qcount in ((1,1),(3,1),(1,3)):
                row=recovery_certificate(Fraction(4,5),Fraction(4,5),lam,ecount,qcount)
                success=0.
                for s in (-1,1):
                    leaves=recovery_records(calibration_state(H,s),.8,.8,float(lam),ecount,qcount)
                    self.assertAlmostEqual(sum(np.trace(b).real for _,b in leaves.values()),1,places=13)
                    success+=sum(np.trace(b).real for guess,b in leaves.values() if guess==s)/2
                self.assertAlmostEqual(success,row["finite_H_success"]["diagnostic"],places=13)

    def test_conditional_query_choice_can_outperform_one_fixed_correction_rule(self):
        gamma=sum(visibility_interval(1).floats())/2
        self.assertEqual(conditional_query_strategy(.8,.95,gamma,1),("joint",1))
        self.assertEqual(conditional_query_strategy(.8,.95,gamma,-1),("second_only",1))
        self.assertGreater(recovery_factor(.8,.95,gamma),max(.8,.95,gamma))
        row=recovery_certificate(Fraction(4,5),Fraction(4,5),Fraction(19,20),1,1)
        self.assertEqual([x["strategy_certificate"] for x in row["choices"]],["second_only","joint_positive"])

    def test_environment_queries_do_not_change_the_old_unconditional_complex_reference_state(self):
        state=random_state(np.random.default_rng(302),8)
        leaves=recovery_records(state,.8,.6,.4,3,1)
        np.testing.assert_allclose(sum(forget_memories(b) for _,b in leaves.values()),
                                   old_channel(state,.8,.6),atol=2e-15)

    def test_noiseless_and_uninformative_environment_limits_and_no_gain_above_original_information(self):
        for c in (0,.4,.8,1):
            for lam in (0,.3,.9,1):
                self.assertAlmostEqual(recovery_factor(c,lam,1),1)
                self.assertAlmostEqual(recovery_factor(c,lam,0),max(c,lam))
                for gamma in (.2,.7,.99):
                    factor=recovery_factor(c,lam,gamma)
                    self.assertGreaterEqual(factor,max(c,lam)-1e-15)
                    self.assertLessEqual(factor,1+1e-15)

    def test_G_query_information_survives_environment_messages_and_H_gain_is_certified(self):
        branches=[environment_record_branches(memory_channel(calibration_state(G,s),.8,.8),.8,0,1)
                  for s in (-1,1)]
        self.assertAlmostEqual(sum(trace_distance(branches[0][r],branches[1][r]) for r in (-1,1)),.6)
        row=recovery_certificate(Fraction(4,5),Fraction(4,5),Fraction(0),1,1)
        self.assertGreater(Fraction(row["strict_gain"]["lower_exact"]),Fraction(56,1000))
        self.assertFalse(row["full_quantum_memory_state_restored_by_this_record_only_protocol"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(EnvironmentQueryRecoveryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":102,"environment_access_assumption":"Retain and control the pure initialized environment in the round 100 dilation",
        "message_probability":"p_r=(1+r*lam*gamma)/2",
        "signed_phase_weight":"a_r=(lam+r*gamma)/2",
        "ideal_H_distinguishability_for_fixed_environment_message":"w*sum_r max(c*p_r,abs(a_r))",
        "query_feedback":"Read N if c*p_r>=abs(a_r); otherwise read joint H axis and multiply by sign(a_r)",
        "examples":[recovery_certificate(Fraction(4,5),Fraction(4,5),lam,e,q)
                    for lam,e,q in ((Fraction(0),1,1),(Fraction(0),3,1),(Fraction(0),3,3),(Fraction(19,20),1,1))],
        "arbitrary_environment_or_feedback_optimized":False,
        "earlier_classical_G_record_can_be_changed_into_eraser_record_after_reading":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("environment_query_recovery_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
