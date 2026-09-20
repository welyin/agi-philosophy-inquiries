"""Round 156: optimize adaptive reads using sharp per-confusion quantum costs.

The finite Bellman problem is solved exactly at the certified parameter. Its
objective sums separately worst conditional costs, so its optimum is still an
upper-bound certificate for quantum recovery, not a minimax lower bound.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from adaptive_history_readout import AdaptiveSolution, adaptive_solution, updated
from budgeted_history_readout import ReadModel, allocations, fixed_polynomial, add_polynomials
from certified_intervals import Interval as I, SCALE
from history_recovery_loss_matrix import radical_basis, weighted_coefficients, symbols, cost_interval
from joint_relation_records import interval_fields


class DamageModel:
    def __init__(self,total,coarse=False,contrast=F(1,2)):
        old=ReadModel(total,coarse,contrast)
        self.total=old.total; self.coarse=old.coarse; self.width=old.width
        self.labels=old.labels; self.contrast=old.contrast; self.p=old.p; self.q=old.q
        self.weights=weighted_coefficients(coarse)
        powers=tuple(self.p**k*self.q**(total-k)/48 for k in range(total+1))
        self.basis=tuple(root*power for root in radical_basis() for power in powers)
        self.zero=(0,)*(5*(total+1)); self._values={}; self.comparisons=0

    def value(self,polynomial):
        if polynomial not in self._values:
            self._values[polynomial]=sum((n*b for n,b in zip(polynomial,self.basis) if n),I.exact(0))
        return self._values[polynomial]

    def choose(self,candidates):
        if not candidates: raise ValueError("At least one candidate required.")
        best=min(range(len(candidates)),key=lambda i:self.value(candidates[i]).lo+self.value(candidates[i]).hi)
        for other in candidates:
            if other==candidates[best]: continue
            difference=tuple(a-b for a,b in zip(other,candidates[best]))
            if self.value(difference).lo<=0:
                raise ArithmeticError("100-bit intervals do not separate nonidentical cost polynomials.")
            self.comparisons+=1
        return best

    def report_polynomial(self,signed,report):
        if type(report) is not int or report not in range(self.labels): raise ValueError("Invalid report.")
        if (len(signed)!=self.width or any(type(x) is not int for x in signed)
                or sum(abs(x) for x in signed)>self.total or (self.total-sum(signed))%2):
            raise ValueError("Unreachable terminal signed counts.")
        coefficients=list(self.zero)
        for true in range(self.labels):
            dot=sum((1 if (true>>(self.width-1-j))&1 else -1)*x for j,x in enumerate(signed))
            matches=(self.total+dot)//2
            for radical,n in enumerate(self.weights[true][report]):
                coefficients[radical*(self.total+1)+matches]+=n
        return tuple(coefficients)

    def terminal(self,signed):
        candidates=[self.report_polynomial(signed,g) for g in range(self.labels)]
        best=self.choose(candidates)
        return candidates[best],best


class DamageSolution(AdaptiveSolution):
    def __init__(self,total,coarse=False,contrast=F(1,2)):
        self.model=DamageModel(total,coarse,contrast); self.nodes={}
        self.polynomial=self.solve_state(0,(0,)*self.model.width)[0]
        self.error_interval=self.model.value(self.polynomial)


@lru_cache(maxsize=None)
def damage_solution(total,coarse=False,contrast=F(1,2)):
    return DamageSolution(total,coarse,contrast)


def evaluate_existing_policy(policy,model):
    if (policy.model.total,policy.model.coarse)!=(model.total,model.coarse):
        raise ValueError("Policy and cost model must have the same horizon and labels.")
    @lru_cache(maxsize=None)
    def visit(depth,signed):
        choice=policy.nodes[(depth,signed)][1]
        if depth==model.total: return model.report_polynomial(signed,choice)
        return add_polynomials(visit(depth+1,updated(signed,choice,0)),visit(depth+1,updated(signed,choice,1)))
    return visit(0,(0,)*model.width)


def damage_certificate(total,coarse=False):
    solution=damage_solution(total,coarse); model=solution.model
    old=adaptive_solution(total,coarse)
    old_evaluated=model.value(evaluate_existing_policy(old,model))
    return {"labels":model.labels,"reads":total,"contrast_exact":str(model.contrast),
        "old_label_risk_certificate":interval_fields(old.error_interval),
        "old_policy_with_sharp_costs":interval_fields(old_evaluated),
        "new_policy_with_sharp_costs":interval_fields(solution.error_interval),
        "loss_polynomial":solution.polynomial,"bellman_states":len(solution.nodes),
        "certified_strict_comparisons":model.comparisons,
        "reachable_policy_states":len(solution.policy_rows()),
        "new_native_gates_after_first_join":(405 if coarse else 401)+4*total,
        "fresh_pointer_initializations":total,"optimal_for_this_cost_matrix":True,
        "optimal_actual_quantum_trace_distance_claimed":False,
        "independent_encoded_source_promise_required":True}


def brute_cost(weights,patterns,costs,p,remaining):
    if remaining==0: return float(np.min(weights @ costs))
    return min(sum(brute_cost(weights*np.where(patterns[:,j]==bit,p,1-p),patterns,costs,p,remaining-1)
                   for bit in (0,1)) for j in range(patterns.shape[1]))


class DamageAwareHistoryReadoutTests(unittest.TestCase):
    def test_weighted_terminal_polynomial_matches_direct_conditional_cost_sum(self):
        for coarse,signed in ((False,(1,-1,1)),(True,(2,-1))):
            model=DamageModel(3,coarse); p=sum(model.p.floats())/2
            for report in range(model.labels):
                direct=0.
                for true in range(model.labels):
                    matches=(3+sum((1 if (true>>(model.width-1-j))&1 else -1)*x for j,x in enumerate(signed)))//2
                    prior=1/3 if coarse and true<2 else 1/6
                    direct+=prior*p**matches*(1-p)**(3-matches)*sum(cost_interval(symbols(coarse)[true][report]).floats())/2
                self.assertAlmostEqual(sum(model.value(model.report_polynomial(signed,report)).floats())/2,direct)

    def test_bellman_matches_independent_uncompressed_posterior_recursion(self):
        for coarse in (False,True):
            solution=damage_solution(4,coarse); model=solution.model
            costs=np.array([[sum(cost_interval(s).floats())/2 for s in row] for row in symbols(coarse)])
            priors=np.array([1/3,1/3,1/6,1/6]) if coarse else np.ones(6)/6
            patterns=np.array([[(r>>(model.width-1-j))&1 for j in range(model.width)] for r in range(model.labels)])
            self.assertAlmostEqual(brute_cost(priors,patterns,costs,sum(model.p.floats())/2,4),sum(solution.error_interval.floats())/2)

    def test_sharper_analysis_and_changed_policy_are_separately_verified(self):
        solution=damage_solution(3); old=adaptive_solution(3)
        sharper=solution.model.value(evaluate_existing_policy(old,solution.model))
        self.assertLess(sharper.hi,old.error_interval.lo)
        self.assertLess(solution.error_interval.hi,sharper.lo)
        self.assertAlmostEqual(sum(sharper.floats())/2,.486537189718,places=11)
        self.assertAlmostEqual(sum(solution.error_interval.floats())/2,.480595642525,places=11)

    def test_three_fine_reads_and_four_coarse_reads_certify_below_one_half(self):
        self.assertLess(damage_solution(3).error_interval.hi,SCALE//2)
        self.assertGreater(damage_solution(2).error_interval.lo,SCALE//2)
        self.assertLess(damage_solution(4,True).error_interval.hi,SCALE//2)
        self.assertGreater(damage_solution(3,True).error_interval.lo,SCALE//2)

    def test_adaptation_dominates_all_fixed_allocations_for_the_same_loss(self):
        for coarse in (False,True):
            for total in range(7):
                solution=damage_solution(total,coarse); model=solution.model
                for counts in allocations(total,model.width):
                    fixed=model.value(fixed_polynomial(model,counts))
                    self.assertLessEqual(solution.error_interval.lo,fixed.hi)

    def test_new_policy_dominates_reanalyzed_old_policy_through_the_reported_range(self):
        for coarse in (False,True):
            for n in range(13):
                solution=damage_solution(n,coarse)
                previous=solution.model.value(evaluate_existing_policy(adaptive_solution(n,coarse),solution.model))
                self.assertLessEqual(solution.error_interval.lo,previous.hi)
        solution=damage_solution(8,True)
        previous=solution.model.value(evaluate_existing_policy(adaptive_solution(8,True),solution.model))
        self.assertLess(solution.error_interval.hi,previous.lo)

    def test_policy_export_reconstructs_the_certified_root_loss(self):
        for coarse,n in ((False,3),(True,4)):
            solution=damage_solution(n,coarse)
            self.assertEqual(evaluate_existing_policy(solution,solution.model),solution.polynomial)
            self.assertTrue(all(row.get("physical_wire",4) in range(4,4+solution.model.width)
                                for row in solution.policy_rows()))

    def test_parameter_guards_and_exact_ties(self):
        model=DamageModel(3)
        with self.assertRaises(ValueError): model.report_polynomial((0,0,0),0)
        with self.assertRaises(ValueError): model.report_polynomial((1,1,1),6)
        with self.assertRaises(ValueError): evaluate_existing_policy(adaptive_solution(4),model)
        self.assertEqual(model.choose((model.zero,model.zero)),0)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DamageAwareHistoryReadoutTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":156,"objective":"sum of prior-weighted sharp conditional external recovery costs",
        "radical_basis":["1","sqrt(3)","sqrt(33)","sqrt(25-16sqrt(2))","sqrt(25+16sqrt(2))"],
        "polynomial_convention":"five blocks of N+1 integer coefficients; multiply radical*p^k*q^(N-k)/48",
        "frontier":[{"six":damage_certificate(n),"four":damage_certificate(n,True)} for n in range(13)],
        "policies":{"six_3":damage_solution(3).policy_rows(),"six_5":damage_solution(5).policy_rows(),
                    "four_4":damage_solution(4,True).policy_rows()},
        "all_results_accepted":True,"quantum_minimax_error_or_minimum_read_count_claimed":False,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("damage_aware_history_readout_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"round":156,"checks":report["automated_checks"],"selected":[
        {"reads":n,"six":report["frontier"][n]["six"]["new_policy_with_sharp_costs"]["diagnostic"],
         "four":report["frontier"][n]["four"]["new_policy_with_sharp_costs"]["diagnostic"]}
        for n in (2,3,4,5,7,9,12)]},indent=2))


if __name__=="__main__": main()
