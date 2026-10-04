"""Round 153: certified finite-horizon adaptive allocation of native reads.

Bellman optimization covers every classical policy selecting one of the fixed
label wires after each raw result. The objective is a label-risk certificate,
not the quantum recovery minimax problem. Exact polynomial ties are preserved.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest

import numpy as np

from budgeted_history_readout import ReadModel, add_polynomials, fixed_optimum
from certified_intervals import SCALE
from joint_relation_records import interval_fields


def updated(signed,site,bit):
    if type(site) is not int or site not in range(len(signed)) or bit not in (0,1):
        raise ValueError("Use an available wire and a binary result.")
    return tuple(x+(2*bit-1 if j==site else 0) for j,x in enumerate(signed))


class AdaptiveSolution:
    def __init__(self,total,coarse=False,contrast=F(1,2)):
        self.model=ReadModel(total,coarse,contrast); self.nodes={}
        self.polynomial=self.solve_state(0,(0,)*self.model.width)[0]
        self.error_interval=1-self.model.value(self.polynomial)

    def solve_state(self,depth,signed):
        key=(depth,signed)
        if key in self.nodes: return self.nodes[key]
        if (type(depth) is not int or not 0<=depth<=self.model.total
                or len(signed)!=self.model.width or any(type(x) is not int for x in signed)
                or sum(abs(x) for x in signed)>depth or (depth-sum(signed))%2):
            raise ValueError("Unreachable signed-count state.")
        if depth==self.model.total:
            result=self.model.terminal(signed)
        else:
            candidates=[add_polynomials(self.solve_state(depth+1,updated(signed,j,0))[0],
                                        self.solve_state(depth+1,updated(signed,j,1))[0])
                        for j in range(self.model.width)]
            chosen=self.model.choose(candidates); result=(candidates[chosen],chosen)
        self.nodes[key]=result
        return result

    def policy_rows(self):
        reachable={(0,(0,)*self.model.width)}; rows=[]
        for depth in range(self.model.total+1):
            for _,signed in sorted(node for node in reachable if node[0]==depth):
                choice=self.nodes[(depth,signed)][1]
                row={"depth":depth,"signed_counts":signed}
                if depth==self.model.total: row["report"]=choice
                else:
                    row["read_label_wire"]=choice; row["physical_wire"]=4+choice
                    for bit in (0,1): reachable.add((depth+1,updated(signed,choice,bit)))
                rows.append(row)
        return rows


@lru_cache(maxsize=None)
def adaptive_solution(total,coarse=False,contrast=F(1,2)):
    return AdaptiveSolution(total,coarse,contrast)


def adaptive_certificate(total,coarse=False,contrast=F(1,2)):
    solution=adaptive_solution(total,coarse,contrast); model=solution.model
    return {"labels":model.labels,"reads":total,"contrast_exact":str(F(contrast)),
        "external_error_upper":interval_fields(solution.error_interval),
        "gain_polynomial":solution.polynomial,
        "bellman_states":len(solution.nodes),"certified_strict_comparisons":model.comparisons,
        "reachable_policy_states":len(solution.policy_rows()),
        "fresh_pointer_initializations":total,"raw_readout_bits_retained":total,
        "new_native_gates_after_first_join":(405 if coarse else 401)+4*total,
        "optimal_over_all_classical_adaptive_policies_for_this_certificate":True,
        "optimal_quantum_recovery_or_minimum_read_count_claimed":False,
        "all_results_accepted":True,"classical_policy_computation_in_gate_count":False}


def brute_value(weights,patterns,p,remaining):
    """Independent small-depth recursion on full posterior weights, without state merging."""
    if remaining==0: return float(max(weights))
    return max(sum(brute_value(weights*np.where(patterns[:,j]==bit,p,1-p),patterns,p,remaining-1)
                   for bit in (0,1)) for j in range(patterns.shape[1]))


class AdaptiveHistoryReadoutTests(unittest.TestCase):
    def test_signed_counts_are_sufficient_even_when_query_counts_and_order_differ(self):
        p=F(3,4)
        histories=(((0,0),(0,1),(2,1)),((1,1),(1,0),(2,1)))
        for label in range(6):
            probabilities=[]
            for history in histories:
                probability=F(1); signed=(0,0,0)
                for j,bit in history:
                    signed=updated(signed,j,bit)
                    probability*=p if ((label>>(2-j))&1)==bit else 1-p
                self.assertEqual(signed,(0,0,1)); probabilities.append(probability)
            self.assertEqual(*probabilities)

    def test_bellman_matches_independent_full_history_recursion(self):
        for coarse in (False,True):
            solution=adaptive_solution(4,coarse); model=solution.model
            p=sum(model.p.floats())/2; w=sum(model.weight.floats())/2
            weights=np.array([w,w,1,1])/6 if coarse else np.ones(6)/6
            patterns=np.array([[(r>>(model.width-1-j))&1 for j in range(model.width)] for r in range(model.labels)])
            self.assertAlmostEqual(brute_value(weights,patterns,p,4),1-sum(solution.error_interval.floats())/2)

    def test_adaptation_never_worsens_fixed_allocations_and_can_strictly_improve_them(self):
        for coarse in (False,True):
            for n in range(11):
                self.assertLessEqual(adaptive_solution(n,coarse).error_interval.lo,fixed_optimum(n,coarse)["error_interval"].hi)
        self.assertLess(adaptive_solution(5).error_interval.hi,fixed_optimum(5)["error_interval"].lo)
        self.assertLess(adaptive_solution(4,True).error_interval.hi,fixed_optimum(4,True)["error_interval"].lo)

    def test_thresholds_are_certificates_not_claims_about_all_quantum_decoders(self):
        self.assertLess(adaptive_solution(5).error_interval.hi,SCALE//2)
        self.assertGreater(adaptive_solution(4).error_interval.lo,SCALE//2)
        self.assertLess(adaptive_solution(4,True).error_interval.hi,SCALE//2)
        self.assertGreater(adaptive_solution(3,True).error_interval.lo,SCALE//2)
        self.assertAlmostEqual(sum(adaptive_solution(5).error_interval.floats())/2,.442605819302,places=11)
        self.assertAlmostEqual(sum(adaptive_solution(4,True).error_interval.floats())/2,.491774720475,places=11)

    def test_equal_gate_budget_certificate_order_is_certified_over_the_reported_range(self):
        for n in range(16):
            self.assertLess(adaptive_solution(n+1).error_interval.hi,adaptive_solution(n,True).error_interval.lo)

    def test_exported_policy_reconstructs_the_exact_root_polynomial(self):
        for coarse,n in ((False,5),(True,4)):
            solution=adaptive_solution(n,coarse); rows=solution.policy_rows()
            table={(row["depth"],tuple(row["signed_counts"])):row for row in rows}
            def walk(depth,signed):
                row=table[(depth,signed)]
                if depth==n:
                    polynomial,report=solution.model.terminal(signed)
                    self.assertEqual(row["report"],report); return polynomial
                j=row["read_label_wire"]
                self.assertEqual(row["physical_wire"],4+j)
                return add_polynomials(walk(depth+1,updated(signed,j,0)),walk(depth+1,updated(signed,j,1)))
            self.assertEqual(walk(0,(0,)*solution.model.width),solution.polynomial)

    def test_read_policy_really_depends_on_intermediate_results(self):
        for coarse,n in ((False,5),(True,4)):
            rows=adaptive_solution(n,coarse).policy_rows()
            actions_by_depth={d:{r["read_label_wire"] for r in rows if r["depth"]==d} for d in range(n)}
            self.assertTrue(any(len(values)>1 for values in actions_by_depth.values()))

    def test_reachability_and_terminal_boundary_guards(self):
        solution=adaptive_solution(3)
        for state in ((-1,(0,0,0)),(1,(2,0,0)),(2,(1,0,0)),(4,(0,0,0))):
            with self.assertRaises(ValueError): solution.solve_state(*state)
        with self.assertRaises(ValueError): updated((0,0),2,1)
        self.assertEqual(adaptive_solution(0).polynomial,(1,0))
        self.assertEqual(adaptive_solution(0,True).polynomial,(0,1))


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(AdaptiveHistoryReadoutTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":153,"state":"depth and signed counts of ones minus zeros per label wire",
        "bellman":"V_t(x)=max_j [V_(t+1)(x-e_j)+V_(t+1)(x+e_j)]",
        "fixed_number_of_native_reads_all_outcomes_accepted":True,
        "certificate_frontier":[{"six":adaptive_certificate(n),"four":adaptive_certificate(n,True)} for n in range(17)],
        "policies":{"six_5":adaptive_solution(5).policy_rows(),"four_4":adaptive_solution(4,True).policy_rows()},
        "polynomial_basis":"same as budgeted_history_readout_results.json",
        "equal_gate_budget_six_strictly_better_certificate_four_reads_0_through_15":True,
        "source_state_or_closed_quantum_history_used_by_controller":False,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("adaptive_history_readout_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"round":153,"checks":report["automated_checks"],"selected":[
        {"reads":n,"six":report["certificate_frontier"][n]["six"]["external_error_upper"]["diagnostic"],
         "four":report["certificate_frontier"][n]["four"]["external_error_upper"]["diagnostic"]}
        for n in (4,5,6,7,9,12,15)]},indent=2))


if __name__=="__main__": main()
