"""Round 124: a declared storage-noise model and certified periodic refresh.

No memory noise is derived from cognition. The real and reused complex
references have the same score decay when their initial correlation and local
noise are matched. Preparation costs may still differ.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import unittest
import numpy as np

from certified_intervals import Interval as I, SCALE
from independent_source_alignment import pair_reference
from joint_relation_records import interval_fields
from noisy_imaginarity_distillation import resource_state
from one_bit_real_network import visibility_interval, real_feedback_score
from quantum_interface_audit import IDENTITY, PAULI_Z, PAULI_Y
from reference_information_cost import target_intervals


def local_storage_channel(density,probability):
    q=float(probability)
    if not 0<=q<=.5: raise ValueError("Use flip probability in [0,1/2].")
    return (1-q)*density+q*PAULI_Z@density@PAULI_Z


def real_pair_storage(reference,probability):
    q=float(probability)
    if not 0<=q<=.5: raise ValueError("Use flip probability in [0,1/2].")
    out=np.zeros((4,4),dtype=complex)
    for first,p in ((IDENTITY,1-q),(PAULI_Z,q)):
        for second,r in ((IDENTITY,1-q),(PAULI_Z,q)):
            operator=np.kron(first,second)
            out+=p*r*operator@reference@operator.T
    return out


def retained_correlation(age,probability,initial=None):
    q=F(probability)
    if not isinstance(age,int) or age<0 or not 0<=q<=F(1,2):
        raise ValueError("Nonnegative integer age and flip probability in [0,1/2].")
    g0=visibility_interval(5) if initial is None else I.exact(initial)
    return g0*I.exact((1-2*q)**2)**age


def maximum_age(probability,initial=None,target=None):
    q=F(probability)
    g0=visibility_interval(5) if initial is None else I.exact(initial)
    needed=target_intervals()[2] if target is None else I.exact(target)
    if not 0<=q<=F(1,2) or needed.lo<=0 or needed.hi>SCALE:
        raise ValueError("Use allowed noise and a strictly positive target <=1.")
    def classification(age):
        current=retained_correlation(age,q,g0)
        if current.lo>=needed.hi: return True
        if current.hi<needed.lo: return False
        raise ArithmeticError("Intervals cannot decide this boundary; refine precision.")
    if not classification(0): return -1
    if q==0: return None
    lo,hi=0,1
    while classification(hi):
        lo,hi=hi,2*hi
        if hi>2**50: raise ArithmeticError("Age exceeds the chosen certification range.")
    while hi-lo>1:
        mid=(lo+hi)//2
        if classification(mid): lo=mid
        else: hi=mid
    return lo


def refresh_resources(trials,task_span):
    if not isinstance(trials,int) or trials<1 or not isinstance(task_span,int) or task_span<1:
        raise ValueError("Positive integer trial and block lengths required.")
    blocks=(trials+task_span-1)//task_span
    return {"task_trials":trials,"tasks_per_reference_block":task_span,"alignment_blocks":blocks,
            "fresh_target_Bell_pairs":2*trials,"reference_link_pairs":2*blocks,
            "reference_source_rebits_allocated_and_retained":4*blocks,
            "total_original_readouts":20*trials+5*blocks,
            "alignment_message_bits":blocks,
            "pair_rotations_after_source_preparation_upper":23*trials+7*blocks,
            "source_preparation_classical_audit_and_memory_maintenance_costs_separate":True}


@lru_cache(maxsize=1)
def maintenance_certificate():
    target=target_intervals()[2]
    rows=[]
    for q in (F(0),F(1,10**7),F(1,10**6),F(1,10**5)):
        last=maximum_age(q)
        row={"local_Z_flip_probability_per_intertask_cycle":str(q),
             "largest_certified_age":last,
             "task_span":None if last is None else last+1}
        if last is not None:
            row["last_valid_correlation"]=interval_fields(retained_correlation(last,q))
            row["first_invalid_correlation"]=interval_fields(retained_correlation(last+1,q))
        rows.append(row)
    span=maximum_age(F(1,10**6))+1
    return {"target_correlation_for_original_finite_score":interval_fields(target),
            "storage_noise_rows":rows,"example_1000_tasks":refresh_resources(1000,span),
            "noise_and_exact_resets_are_declared_assumptions":True,
            "criterion":"each preplanned task's unconditional network score at least the round 69 value",
            "same_initial_correlation_and_local_noise_give_same_real_complex_score_decay":True,
            "does_not_certify_all_conditional_histories_or_iid_full_records":True}


class ReferenceMaintenanceBudgetTests(unittest.TestCase):
    def test_storage_decay_follows_the_actual_real_CPTP_channel(self):
        for q in (0,.001,.1,.5):
            for g in (0,.6,1):
                np.testing.assert_allclose(real_pair_storage(pair_reference(g),q),
                                           pair_reference(g*(1-2*q)**2),atol=3e-16,rtol=0)

    def test_complex_private_reference_product_has_the_same_noise_decay(self):
        g=.8
        bias=np.sqrt(g)
        for q in (0,.001,.1,.5):
            rho=resource_state(bias)
            reference=pair_reference(g)
            for _ in range(5):
                rho=local_storage_channel(rho,q)
                reference=real_pair_storage(reference,q)
                complex_g=np.trace(rho@PAULI_Y).real**2
                real_g=np.trace(reference@np.kron(PAULI_Y,PAULI_Y)).real
                self.assertAlmostEqual(complex_g,real_g,places=14)
                self.assertAlmostEqual(real_feedback_score(complex_g),real_feedback_score(real_g),places=13)

    def test_zero_storage_noise_requires_no_periodic_refresh(self):
        self.assertIsNone(maximum_age(0))
        initial=visibility_interval(5)
        later=retained_correlation(10**6,0)
        self.assertEqual((initial.lo,initial.hi),(later.lo,later.hi))

    def test_certified_boundary_has_no_off_by_one_task_error(self):
        last=maximum_age(F(1,10**6))
        self.assertEqual(last,13)
        target=target_intervals()[2]
        self.assertGreaterEqual(retained_correlation(last,F(1,10**6)).lo,target.hi)
        self.assertLess(retained_correlation(last+1,F(1,10**6)).hi,target.lo)
        self.assertEqual(maintenance_certificate()["example_1000_tasks"]["tasks_per_reference_block"],14)

    def test_integer_search_matches_exact_rational_bruteforce(self):
        initial,target=F(9,10),F(1,2)
        for q in (F(1,100),F(1,10),F(1,2)):
            last=max(t for t in range(100) if initial*(1-2*q)**(2*t)>=target)
            self.assertEqual(maximum_age(q,initial,target),last)

    def test_refresh_cost_counts_setup_and_keeps_old_reference_slots(self):
        cost=refresh_resources(1000,14)
        self.assertEqual(cost["alignment_blocks"],72)
        self.assertEqual(cost["total_original_readouts"],20360)
        self.assertEqual(cost["alignment_message_bits"],72)
        self.assertEqual(cost["reference_source_rebits_allocated_and_retained"],288)
        self.assertEqual(cost["pair_rotations_after_source_preparation_upper"],23504)

    def test_every_scheduled_age_is_certified_and_larger_blocks_fail(self):
        last=maximum_age(F(1,10**6))
        target=target_intervals()[2]
        for age in range(last+1):
            self.assertGreaterEqual(retained_correlation(age,F(1,10**6)).lo,target.hi)
        self.assertLess(retained_correlation(last+1,F(1,10**6)).hi,target.lo)

    def test_invalid_input_and_initially_insufficient_reference_are_explicit(self):
        self.assertEqual(maximum_age(F(1,10),F(1,4),F(1,2)),-1)
        with self.assertRaises(ValueError): retained_correlation(-1,F(1,10))
        with self.assertRaises(ValueError): maximum_age(F(-1,10))
        with self.assertRaises(ValueError): refresh_resources(10,0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReferenceMaintenanceBudgetTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":124,"maintenance_certificate":maintenance_certificate(),
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0},
            "quantum_theory_derived_from_cognition":False}
    if args.write_results:
        Path(__file__).with_name("reference_maintenance_budget_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
