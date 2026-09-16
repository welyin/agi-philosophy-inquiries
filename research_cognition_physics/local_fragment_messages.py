"""Round 105: local fragment messages, useful record compression, and access thresholds."""
import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I
from complex_whole_real_interfaces import random_state
from deferred_relation_queries import memory_channel
from environment_query_recovery import branch_moments,interval_abs,interval_max
from fragment_query_compiler import partial_factor,query_branches
from incompatible_relation_writes import G,H,calibration_state
from independent_source_alignment import keep_systems
from joint_relation_records import axis_read_memory_tree,interval_fields,validated_count
from one_bit_real_network import visibility_interval
from partial_environment_access import accessible_state,partition_overlaps
from memory_dephasing_threshold import dephase_memory
from joint_relation_records import memory_axis
from weak_relation_tradeoff import forget_front_pointer,trace_distance


def local_fragment_records(state,c,overlaps,accessible,counts):
    """Separate chosen eraser measurements; no coherent gates among fragments."""
    accessible=tuple(accessible);counts=tuple(counts)
    partition_overlaps(overlaps,accessible)
    if len(accessible)!=len(counts): raise ValueError("One read count per accessible fragment.")
    leaves={():accessible_state(state,c,overlaps,accessible)}
    for index,count in zip(accessible,counts):
        following={}
        for prefix,branch in leaves.items():
            for history,output in axis_read_memory_tree(branch,0,math.acos(overlaps[index]),count).items():
                r=1 if sum(history)>count//2 else -1
                key=prefix+(r,)
                reduced=forget_front_pointer(output)
                following[key]=following.get(key,np.zeros_like(reduced))+reduced
        leaves=following
    return leaves


def record_moments(overlaps,gammas,record):
    pairs=[branch_moments(lam,gamma,r) for lam,gamma,r in zip(overlaps,gammas,record)]
    return math.prod(p for p,t in pairs),math.prod(t for p,t in pairs)


def analytic_local_branch(state,c,b,overlaps,gammas,record):
    p,t=record_moments(overlaps,gammas,record)
    axis=np.kron(memory_axis(c),np.eye(len(state)//2))
    return (p+b*t)*state/2+(p-b*t)*axis@state@axis/2


def record_action(c,b,p,t):
    return ("joint",1 if t>=0 else -1) if b*abs(t)>c*p else ("second_only",1)


def full_record_factor(c,b,overlaps,gammas):
    return sum(max(c*p,b*abs(t)) for record in product((-1,1),repeat=len(overlaps))
               for p,t in [record_moments(overlaps,gammas,record)])


def homogeneous_certificate(c,d,lam,total_fragments,accessible_count,environment_count,query_count):
    c,d,lam=map(Fraction,(c,d,lam))
    if any(not 0<=x<=1 for x in (c,d,lam)): raise ValueError("Use overlaps in [0,1].")
    if not 0<=accessible_count<=total_fragments: raise ValueError("Invalid fragment counts.")
    validated_count(environment_count);validated_count(query_count)
    n,k,m,q=total_fragments,accessible_count,environment_count,query_count
    gamma=visibility_interval(m);eta=visibility_interval(q)
    a=lam**k;b=lam**(n-k)
    pp=(1+I.exact(lam)*gamma)/2;pm=(1-I.exact(lam)*gamma)/2
    tp=(I.exact(lam)+gamma)/2;tm=interval_abs((I.exact(lam)-gamma)/2)
    def power(x,n):
        out=I.exact(1)
        for _ in range(n): out*=x
        return out
    # Build weighted binomial masses before rounding tiny terms. Multiplying
    # already-rounded powers by C(n,j) would amplify 2^-100 errors up to O(1).
    def binomial_masses(left,right,n):
        masses=[I.exact(1)]
        for _ in range(n):
            following=[I.exact(0) for _ in range(len(masses)+1)]
            for j,mass in enumerate(masses):
                following[j]+=left*mass
                following[j+1]+=right*mass
            masses=following
        return masses
    factor=I.exact(0)
    for p,t in zip(binomial_masses(pp,pm,k),binomial_masses(tp,tm,k)):
        factor+=interval_max(I.exact(c)*p,I.exact(b)*t)
    parity_gamma=power(gamma,k)
    parity=I.exact(0)
    for r in (-1,1):
        p=(1+r*I.exact(a)*parity_gamma)/2
        t=I.exact(b)*(I.exact(a)+r*parity_gamma)/2
        parity+=interval_max(I.exact(c)*p,interval_abs(t))
    witness=I.exact(b)*power(tp,k)-I.exact(c)*power(pp,k)
    w=(1-I.exact(d*d)).sqrt()
    return {"total_fragments":n,"accessible_fragments":k,"overlap_exact":str(lam),
        "per_fragment_readouts":m,"final_readouts":q,"inaccessible_overlap_exact":str(b),
        "full_majority_record_factor":interval_fields(factor),"parity_only_factor":interval_fields(parity),
        "full_majority_record_H_success":interval_fields((1+eta*w*factor)/2),
        "parity_only_H_success":interval_fields((1+eta*w*parity)/2),
        "strict_full_record_over_parity_factor":interval_fields(factor-parity),
        "all_plus_gain_witness":interval_fields(witness),
        "environment_raw_readouts":k*m,"fragment_summary_bits_to_receiver_if_distributed":k,
        "fixed_query_action_bits_after_classical_aggregation":2 if k else 0,
        "action_compression_preserves_full_quantum_instrument":False,
        "raw_histories_within_each_majority_optimized":False}


def ideal_access_threshold(c,lam,n):
    c,lam=Fraction(c),Fraction(lam)
    if not 0<=c<=1 or not 0<=lam<=1 or n<0: raise ValueError("Invalid parameters.")
    return next((k for k in range(n+1) if lam**(n-k)>c),None)


def distributed_budget(n,k,m,q):
    validated_count(m);validated_count(q)
    if not 0<=k<=n: raise ValueError("Invalid fragment counts.")
    em=k*m if m>1 else 0;qm=q if q>1 else 0
    slots=n+2+(k if m>1 else 0)+int(q>1)
    return {"scope":"One protected pointer per accessed site; worst joint query branch; old writes included",
        "pair_rotations":4+n+em+1+qm,
        "extra_local_rotations":6+5*n+(k*(3*m+1) if m>1 else 0)+5+3*qm,
        "global_extra_quantum_slots_peak":slots,"initial_pure_preparations":slots,
        "pointer_resets":k*(m-1)+(q-1),"total_readouts":k*m+q,
        "classical_fragment_to_receiver_bits":k,"coherent_fragment_to_fragment_gates":0}


class LocalFragmentMessageTests(unittest.TestCase):
    def test_actual_separate_measurements_match_every_analytic_message_branch(self):
        rho=random_state(np.random.default_rng(105),4)
        overlaps=(.7,.9,.95);subset=(2,0);counts=(1,3)
        _,b=partition_overlaps(overlaps,subset)
        gamma=[sum(visibility_interval(m).floats())/2 for m in counts]
        actual=local_fragment_records(rho,.8,overlaps,subset,counts)
        for record,branch in actual.items():
            np.testing.assert_allclose(branch,analytic_local_branch(
                rho,.8,b,[overlaps[i] for i in subset],gamma,record),atol=2e-15)
        np.testing.assert_allclose(sum(actual.values()),dephase_memory(rho,.8,math.prod(overlaps)),atol=2e-15)

    def test_H_message_blocks_and_parity_coarsening_have_exact_trace_norms(self):
        overlaps=(.95,.95);gammas=[sum(visibility_interval(1).floats())/2]*2
        pair=[local_fragment_records(memory_channel(calibration_state(H,s),.8,.8),.8,overlaps,(0,1),(1,1))
              for s in (-1,1)]
        full=sum(trace_distance(pair[0][key],pair[1][key]) for key in pair[0])
        parity=sum(trace_distance(
            sum(rho for key,rho in pair[0].items() if math.prod(key)==r),
            sum(rho for key,rho in pair[1].items() if math.prod(key)==r)) for r in (-1,1))
        self.assertAlmostEqual(full,.6*full_record_factor(.8,1,overlaps,gammas),places=13)
        self.assertAlmostEqual(parity,.6*partial_factor(.8,math.prod(overlaps),1,math.prod(gammas)),places=13)
        self.assertGreater(full,parity+1e-4)

    def test_three_action_classical_summary_preserves_optimal_H_value_and_finite_query(self):
        overlaps=(.2,.95);gamma=sum(visibility_interval(1).floats())/2
        grouped=[];success=0.
        for s in (-1,1):
            records=local_fragment_records(memory_channel(calibration_state(H,s),.8,.8),.8,overlaps,(0,1),(1,1))
            groups={}
            for key,branch in records.items():
                action=record_action(.8,1,*record_moments(overlaps,(gamma,gamma),key))
                groups[action]=groups.get(action,np.zeros_like(branch))+branch
            grouped.append(groups)
            for action,branch in groups.items():
                success+=sum(np.trace(rho).real for guess,rho in
                    query_branches(branch,.8,.8,1,*action).values() if guess==s)/2
        self.assertEqual(len(grouped[0]),3)
        actual=sum(trace_distance(grouped[0][a],grouped[1][a]) for a in grouped[0])
        expected=.6*full_record_factor(.8,1,overlaps,(gamma,gamma))
        self.assertAlmostEqual(actual,expected,places=13)
        self.assertAlmostEqual(success,(1+gamma*expected)/2,places=13)

    def test_homogeneous_binomial_certificate_matches_full_small_record_enumeration(self):
        f=Fraction;gamma=sum(visibility_interval(1).floats())/2
        for k in (0,1,2,5):
            row=homogeneous_certificate(f(4,5),f(4,5),f(19,20),5,k,1,1)
            self.assertAlmostEqual(row["full_majority_record_factor"]["diagnostic"],
                full_record_factor(.8,.95**(5-k),[.95]*k,[gamma]*k),places=13)

    def test_96_vs_97_local_access_threshold_has_strict_interval_certificate(self):
        f=Fraction
        self.assertEqual(ideal_access_threshold(f(4,5),f(19,20),100),96)
        previous=homogeneous_certificate(f(4,5),f(4,5),f(19,20),100,96,1,1)
        nextrow=homogeneous_certificate(f(4,5),f(4,5),f(19,20),100,97,1,1)
        self.assertLess(Fraction(previous["all_plus_gain_witness"]["upper_exact"]),0)
        self.assertGreater(Fraction(nextrow["all_plus_gain_witness"]["lower_exact"]),0)
        self.assertGreater(Fraction(nextrow["strict_full_record_over_parity_factor"]["lower_exact"]),0)

    def test_ideal_messages_reach_access_bound_but_perfect_unavailable_record_still_blocks(self):
        for b in (0,.7,1):
            for overlaps in ((0.,0.),(.3,.95),(.9,.9,.9)):
                self.assertAlmostEqual(full_record_factor(.8,b,overlaps,[1]*len(overlaps)),max(.8,b),places=13)
                self.assertLessEqual(full_record_factor(.8,b,overlaps,[.99]*len(overlaps)),max(.8,b)+1e-15)
        self.assertEqual(ideal_access_threshold(Fraction(4,5),0,10),10)

    def test_fragment_alone_stores_G_but_not_H_and_redundant_records_are_not_independent_G_samples(self):
        c=.8;overlaps=(.6,.9)
        for obs,target in ((G,.6*math.sqrt(1-(.6*.9)**2)),(H,0.)):
            memories=[memory_channel(calibration_state(obs,s),c,.8) for s in (-1,1)]
            pair=[keep_systems(accessible_state(rho,c,overlaps,(0,1)),(0,1),4) for rho in memories]
            self.assertAlmostEqual(trace_distance(*pair),target,places=13)

    def test_distributed_budget_counts_separate_pointers_and_messages(self):
        row=distributed_budget(100,100,1,1)
        self.assertEqual(row["pair_rotations"],105)
        self.assertEqual(row["total_readouts"],101)
        self.assertEqual(row["classical_fragment_to_receiver_bits"],100)
        row=distributed_budget(100,100,3,3)
        self.assertEqual(row["global_extra_quantum_slots_peak"],203)
        self.assertEqual(row["pointer_resets"],202)

    def test_large_binomial_intervals_are_informative_and_respect_accessible_bounds(self):
        f=Fraction
        for k,m in ((96,1),(97,1),(100,1),(96,3),(100,3)):
            row=homogeneous_certificate(f(4,5),f(4,5),f(19,20),100,k,m,1)
            value=row["full_majority_record_factor"]
            lo,hi=map(Fraction,(value["lower_exact"],value["upper_exact"]))
            self.assertLess(hi-lo,f(1,10**23))
            self.assertGreaterEqual(lo,f(4,5)-f(1,10**23))
            self.assertLessEqual(hi,max(f(4,5),f(19,20)**(100-k))+f(1,10**23))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(LocalFragmentMessageTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    f=Fraction
    report={"round":105,
        "full_message_H_factor":"sum_record max(c*product_i p_i,b*abs(product_i t_i))",
        "parity_coarsening":"Round104 formula with gamma=product_i gamma_i",
        "query_specific_classical_compression":"At most three actions: N, joint+, joint-; two fixed bits suffice after aggregation",
        "examples":[homogeneous_certificate(f(4,5),f(4,5),f(19,20),100,k,m,1)
                    for k,m in ((96,1),(97,1),(100,1),(96,3),(100,3))],
        "distributed_budgets":[distributed_budget(100,100,m,1) for m in (1,3)],
        "all_environment_measurement_bases_optimized":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("local_fragment_messages_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
