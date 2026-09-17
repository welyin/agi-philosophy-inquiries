"""Round 152: fixed read allocations and full-count label decisions.

Six-label success is the exact classification objective. Four-label gain uses
a proved conditional quantum error bound; optimizing it optimizes a certificate,
not the actual quantum trace distance. All polynomial comparisons are certified.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import unittest

import numpy as np

from centered_history_recovery import GROUPS, centered_decoder, exact_error_interval, purification
from certified_intervals import Interval as I, SCALE, sin_interval
from classical_joining_history import independent_encoding
from coarsened_history_bound import pure_extension_metrics
from flagged_subject_joining import random_state
from joint_relation_records import interval_fields
from minimal_coherent_joining_history import conditional_unitary
from noisy_classical_history import majority_error
from noisy_coherent_history_recovery import label_risk
from noisy_two_bit_history import noisy_error_upper


def allocations(total,width):
    if type(total) is not int or total<0 or type(width) is not int or width<1:
        raise ValueError("Nonnegative integer budget and positive integer width required.")
    if width==1:
        yield (total,)
    else:
        for first in range(total+1):
            for tail in allocations(total-first,width-1): yield (first,)+tail


def add_polynomials(first,second):
    return tuple(a+b for a,b in zip(first,second))


class ReadModel:
    """Gain polynomial sum_k (a_k+b_k*w) p^k q^(N-k)/6, w=2-3*delta4."""
    def __init__(self,total,coarse=False,contrast=F(1,2)):
        if type(total) is not int or total<0: raise ValueError("Nonnegative integer budget required.")
        eta=F(contrast)
        if not 0<eta<=1: raise ValueError("Use 0 < contrast <= 1.")
        self.total=total; self.coarse=bool(coarse); self.width=2 if coarse else 3
        self.labels=4 if coarse else 6; self.contrast=eta
        self.p=(1+4*sin_interval(F(1,4))*eta)/2; self.q=1-self.p
        self.weight=2-3*exact_error_interval()
        self.basis=tuple(self.p**k*self.q**(total-k)/6 for k in range(total+1))
        self.zero=(0,)*(2*(total+1))
        self._values={}; self.comparisons=0

    def value(self,polynomial):
        if polynomial not in self._values:
            n=self.total+1
            self._values[polynomial]=sum(((a+self.weight*b)*power
                for a,b,power in zip(polynomial[:n],polynomial[n:],self.basis)),I.exact(0))
        return self._values[polynomial]

    def choose(self,candidates):
        """Return first maximizer; equal coefficient vectors are exact ties."""
        if not candidates: raise ValueError("At least one candidate required.")
        best=max(range(len(candidates)),key=lambda i:self.value(candidates[i]).lo+self.value(candidates[i]).hi)
        for other in candidates:
            if other==candidates[best]: continue
            difference=tuple(a-b for a,b in zip(candidates[best],other))
            if self.value(difference).lo<=0:
                raise ArithmeticError("100-bit intervals do not resolve a nonidentical polynomial comparison.")
            self.comparisons+=1
        return best

    def terminal(self,signed):
        if (len(signed)!=self.width or any(type(x) is not int for x in signed)
                or sum(abs(x) for x in signed)>self.total
                or (self.total-sum(signed))%2):
            raise ValueError("Signed counts must be reachable at the terminal budget.")
        candidates=[]
        for label in range(self.labels):
            dot=sum((1 if (label>>(self.width-1-j))&1 else -1)*x for j,x in enumerate(signed))
            matches=(self.total+dot)//2
            coefficients=list(self.zero)
            coefficients[matches+(self.total+1 if self.coarse and label<2 else 0)]=1
            candidates.append(tuple(coefficients))
        label=self.choose(candidates)
        return candidates[label],label


def fixed_polynomial(model,counts):
    if (len(counts)!=model.width or any(type(n) is not int or n<0 for n in counts)
            or sum(counts)!=model.total): raise ValueError("Allocation must exhaust the stated budget.")
    result=model.zero
    for ones in product(*(range(n+1) for n in counts)):
        multiplicity=math.prod(math.comb(n,k) for n,k in zip(counts,ones))
        polynomial,_=model.terminal(tuple(2*k-n for n,k in zip(counts,ones)))
        result=add_polynomials(result,tuple(multiplicity*c for c in polynomial))
    return result


@lru_cache(maxsize=None)
def fixed_optimum(total,coarse=False,contrast=F(1,2)):
    model=ReadModel(total,coarse,contrast)
    choices=list(allocations(total,model.width)); scores=[fixed_polynomial(model,ns) for ns in choices]
    best=model.choose(scores)
    return {"reads":total,"allocation":choices[best],"gain_polynomial":scores[best],
            "error_interval":1-model.value(scores[best]),"certified_strict_comparisons":model.comparisons}


def fixed_certificate(total,coarse=False,contrast=F(1,2)):
    result=fixed_optimum(total,coarse,contrast)
    return {"labels":4 if coarse else 6,"contrast_exact":str(F(contrast)),
        "reads":total,"allocation":result["allocation"],
        "external_error_upper":interval_fields(result["error_interval"]),
        "gain_polynomial":result["gain_polynomial"],
        "certified_strict_comparisons":result["certified_strict_comparisons"],
        "fresh_pointer_initializations":total,"raw_readout_bits_retained":total,
        "new_yx_after_first_join":(292 if coarse else 290)+total,
        "new_ry_after_first_join":(113 if coarse else 111)+3*total,
        "new_native_gates_after_first_join":(405 if coarse else 401)+4*total,
        "optimal_fixed_allocation_for_this_certificate":True,
        "optimal_quantum_recovery_claimed":False}


class BudgetedHistoryReadoutTests(unittest.TestCase):
    def test_conditional_pair_error_is_three_halves_of_the_averaged_floor(self):
        rng=np.random.default_rng(152); expected=sum((F(3,2)*exact_error_interval()).floats())/2
        for _ in range(4):
            omega=independent_encoding((random_state(rng,2),random_state(rng,2)))
            for g,labels in enumerate(GROUPS):
                ks=tuple(centered_decoder(g) @ conditional_unitary(r+1)/np.sqrt(len(labels)) for r in labels)
                self.assertAlmostEqual(pure_extension_metrics(ks,purification(omega))["trace_error"],expected if g<2 else 0.)

    def test_allocation_enumeration_includes_zero_even_and_unequal_counts(self):
        for total in range(8):
            values=list(allocations(total,3))
            self.assertEqual(len(values),math.comb(total+2,2))
            self.assertEqual(len(set(values)),len(values))
        self.assertIn((2,1,3),list(allocations(6,3)))

    def test_full_counts_strictly_improve_the_old_nine_read_hard_majority(self):
        new=fixed_optimum(9)["error_interval"]; old=label_risk(majority_error(3))
        self.assertLess(new.hi,old.lo)
        self.assertAlmostEqual(sum(new.floats())/2,.3532187812018871,places=12)
        model=ReadModel(9)
        # Both observed first bits have majority 1, but with different confidence.
        self.assertEqual(model.terminal((3,1,-3))[1],4)
        self.assertEqual(model.terminal((1,3,-3))[1],2)

    def test_weighted_terminal_objective_uses_recovery_cost_not_only_group_prior(self):
        model=ReadModel(1,True,F(1,5))
        self.assertEqual(model.terminal((1,0))[1],2)
        p=sum(model.p.floats())/2
        self.assertGreater((1-p)/3,p/6)  # Ordinary group MAP would instead choose group 0.

    def test_equal_three_reads_reproduce_the_previous_four_label_certificate(self):
        model=ReadModel(6,True); value=1-model.value(fixed_polynomial(model,(3,3)))
        previous=noisy_error_upper(majority_error(3))
        self.assertLessEqual(value.lo,previous.hi); self.assertGreaterEqual(value.hi,previous.lo)

    def test_five_reads_suffice_for_both_fixed_certificate_families(self):
        for coarse in (False,True):
            self.assertLess(fixed_optimum(5,coarse)["error_interval"].hi,SCALE//2)
            self.assertGreater(fixed_optimum(4,coarse)["error_interval"].lo,SCALE//2)
        self.assertEqual(fixed_optimum(5)["allocation"],(1,1,3))

    def test_equal_native_gate_budgets_favor_the_six_label_certificate_on_checked_range(self):
        for four_reads in range(16):
            six=fixed_certificate(four_reads+1); four=fixed_certificate(four_reads,True)
            self.assertEqual(six["new_native_gates_after_first_join"],four["new_native_gates_after_first_join"])
            self.assertLess(fixed_optimum(four_reads+1)["error_interval"].hi,
                            fixed_optimum(four_reads,True)["error_interval"].lo)

    def test_parameter_guards_and_exact_polynomial_ties(self):
        for bad in (-1,True):
            with self.assertRaises(ValueError): ReadModel(bad)
        with self.assertRaises(ValueError): ReadModel(3,contrast=0)
        with self.assertRaises(ValueError): fixed_polynomial(ReadModel(3),(1,1,0))
        with self.assertRaises(ValueError): ReadModel(3).terminal((0,0,0))
        model=ReadModel(1)
        self.assertEqual(model.choose((model.zero,model.zero)),0)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BudgetedHistoryReadoutTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":152,"conditional_pair_error":"3 delta4 / 2",
        "four_label_correct_gain_weights_in_sixths":["2-3delta4","2-3delta4","1","1"],
        "gain_polynomial_basis":"sum_k (a_k+b_k*(2-3delta4))*p^k*(1-p)^(N-k)/6; concatenate a,b",
        "certificate_not_exact_quantum_trace_error":True,
        "frontier":[{"six":fixed_certificate(n),"four":fixed_certificate(n,True)} for n in range(17)],
        "equal_gate_budget_comparison":{"four_reads_range":[0,15],"six_reads":"four_reads+1",
            "six_has_strictly_smaller_certificate_throughout_checked_range":True},
        "classical_decision_computation_and_scheduling_time_in_native_gate_count":False,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("budgeted_history_readout_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"round":152,"checks":report["automated_checks"],
        "selected":[report["frontier"][n] for n in (4,5,6,9)]},indent=2))


if __name__=="__main__": main()
