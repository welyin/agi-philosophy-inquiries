"""Round 70: real two-source Bell statistics with one premeasurement bit.

Each source carries a target Bell pair and a real reference pair. Bob aligns
references, sends the majority sign to Charlie, then decodes target Bell bits.
All outputs are retained. This explicitly changes the no-message hypothesis.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from amplified_source_alignment import alignment_record_states
from bell_network_statistics import SETTINGS, bob_effects, independent_bell_sources, network_score, network_table, endpoint_observables
from certified_intervals import Interval as I, SCALE
from finite_network_separation import compiled_axis_gate, positive_z_effect
from independent_source_alignment import independent_sources, keep_systems, pair_reference, real_x_with_reset
from operational_effect_closure import majority_certificate
from real_network_analytic_bound import bound_certificate
from sourcewise_real_network_optimum import one_real_source, outer_axes


def visibility_interval(count):
    row = majority_certificate(count)
    error = I(I.exact(Fraction(row["error_lower_exact"])).lo,I.exact(Fraction(row["error_upper_exact"])).hi)
    return 1-2*error


@lru_cache(maxsize=4)
def reference_flags_by_circuit(count):
    flags = {-1:np.zeros((4,4),dtype=complex),1:np.zeros((4,4),dtype=complex)}
    for record,state in alignment_record_states(count).items():
        flag = 1 if sum(record)>0 else -1
        flags[flag] += keep_systems(state,(0,2),3)
    return flags


def target_assemblage(bob_visibility=1.):
    tensor = independent_bell_sources().reshape(2,4,2,2,4,2)
    return tuple(np.einsum("kj,ajclkd->acld",effect,tensor).reshape(4,4) for effect in bob_effects(bob_visibility))


def records_for_reference(reference,visibility_a=1.,visibility_c=1.,bob_visibility=1.):
    axes = outer_axes()
    alice = tuple(visibility_a*axis for axis in axes)
    charlie = tuple(visibility_c*(axes[i]+sign*axes[j])/math.sqrt(2) for i,j,sign in SETTINGS)
    # Initially target A,C then reference A,C; regroup as A target,reference;
    # C target,reference. This is a tensor-index permutation, not a physical gate.
    states = [keep_systems(np.kron(target,reference),(0,2,1,3),4) for target in target_assemblage(bob_visibility)]
    effects = {(x,z,a,c):np.kron((np.eye(4)+a*alice[x])/2,(np.eye(4)+c*charlie[z])/2)
               for x,z,a,c in product(range(3),range(6),(-1,1),(-1,1))}
    return {(x,z,a,b,c):float(np.trace(states[b]@effects[x,z,a,c]).real)
            for x,z,a,b,c in product(range(3),range(6),(-1,1),range(4),(-1,1))}


def flagged_network_by_circuit(alignment_count=3,read_count=5,feedforward=True):
    visibility = majority_certificate(read_count)["effective_visibility_diagnostic"]
    result = {}
    for flag,reference in reference_flags_by_circuit(alignment_count).items():
        corrected = real_x_with_reset(reference,1,2) if feedforward and flag<0 else reference
        result.update({(flag,)+key:value for key,value in records_for_reference(corrected,visibility,visibility,visibility).items()})
    return result


def forget_flag(records):
    return {key[1:]:sum(records[(flag,)+key[1:]] for flag in (-1,1)) for key in records if key[0]==1}


def real_feedback_score(correlation,visibility_a=1.,visibility_c=1.,bob_visibility=1.):
    return 2*math.sqrt(2)*visibility_a*visibility_c*(bob_visibility*(1+correlation)+bob_visibility**2)


def feedback_certificate(alignment_count=5,read_count=5):
    alignment,read = visibility_interval(alignment_count),visibility_interval(read_count)
    score = 2*I.exact(2).sqrt()*read**2*(read*(1+alignment)+read**2)
    lower = Fraction(score.lo,SCALE)
    previous_upper = Fraction(bound_certificate()["real_upper_rational"])
    return {"alignment_reads":alignment_count,"reads_per_network_target":read_count,
            "score_interval":score.floats(),"score_lower_exact":str(lower),
            "above_round_68_bound":lower>previous_upper,
            "margin_over_round_68_lower_exact":str(lower-previous_upper),
            "cross_lab_messages":1,"message_bits":1,"message_direction":"Bob to Charlie before outer measurements",
            "independent_cross_lab_sources":2,"real_reference_pairs_per_source":1,
            "old_noisy_reads_per_trial":alignment_count+4*read_count,
            "private_new_y_uses":0,"postselect_outputs":False,
            "violates_no_message_real_theorem":False}


class OneBitRealNetworkTests(unittest.TestCase):
    def test_initial_target_and_reference_factorization_uses_only_two_link_sources(self):
        source = np.kron(one_real_source(),one_real_source())
        regrouped = keep_systems(source,(0,2,4,6,1,3,5,7),8)
        np.testing.assert_allclose(regrouped,np.kron(independent_bell_sources(),independent_sources()),atol=2e-17)
        np.testing.assert_allclose(source.imag,0,atol=0)

    def test_actual_noisy_alignment_flags_have_the_expected_subnormalized_outer_states(self):
        for count in (1,3,5):
            gamma = majority_certificate(count)["effective_visibility_diagnostic"]
            for flag,state in reference_flags_by_circuit(count).items():
                np.testing.assert_allclose(state,pair_reference(flag*gamma)/2,atol=2e-15)

    def test_real_reset_correction_makes_both_flags_identical(self):
        flags = reference_flags_by_circuit(3)
        np.testing.assert_allclose(flags[1],real_x_with_reset(flags[-1],1,2),atol=8e-16)
        np.testing.assert_allclose(sum(flags.values()),np.eye(4)/4,atol=8e-16)

    def test_lifted_outer_measurements_are_actual_old_finite_gates_and_noisy_pointers(self):
        gamma = majority_certificate(3)["effective_visibility_diagnostic"]
        directions = list(np.eye(3))+[(np.eye(3)[i]+s*np.eye(3)[j])/math.sqrt(2) for i,j,s in SETTINGS]
        for direction in directions:
            gate = compiled_axis_gate(direction)
            rt_effect = gate.T@np.kron(np.eye(2),positive_z_effect(3))@gate
            tr_effect = keep_systems(rt_effect,(1,0),2)
            expected = (np.eye(4)+gamma*sum(v*axis for v,axis in zip(direction,outer_axes())))/2
            np.testing.assert_allclose(tr_effect,expected,atol=2e-15)

    def test_all_576_flagged_records_normalize_and_recover_the_score_formula(self):
        records = flagged_network_by_circuit(5,5)
        self.assertEqual(len(records),576)
        self.assertGreaterEqual(min(records.values()),0.)
        for x,z in product(range(3),range(6)):
            self.assertAlmostEqual(sum(records[h,x,z,a,b,c] for h,a,b,c in product((-1,1),(-1,1),range(4),(-1,1))),1.,places=13)
        gamma = majority_certificate(5)["effective_visibility_diagnostic"]
        self.assertAlmostEqual(network_score(forget_flag(records)),real_feedback_score(gamma,gamma,gamma,gamma),places=12)

    def test_reference_correlation_is_exactly_the_complex_bias_product_in_all_records(self):
        for correlation in (0.,.3,1.):
            actual = records_for_reference(pair_reference(correlation),.7,.8,.9)
            a = endpoint_observables(.7,correlation)[0]
            c = endpoint_observables(.8,1.)[1]
            expected = network_table(a,c,bob_effects(.9))
            np.testing.assert_allclose(list(actual.values()),list(expected.values()),atol=2e-16)

    def test_disabling_feedforward_removes_only_the_outer_reference_alignment(self):
        without = forget_flag(flagged_network_by_circuit(3,5,False))
        gamma = majority_certificate(5)["effective_visibility_diagnostic"]
        self.assertAlmostEqual(network_score(without),real_feedback_score(0,gamma,gamma,gamma),places=12)

    def test_finite_readout_certificate_crosses_the_previous_bound_with_no_new_y(self):
        report = feedback_certificate()
        self.assertTrue(report["above_round_68_bound"])
        self.assertGreater(Fraction(report["margin_over_round_68_lower_exact"]),Fraction(4,1000))
        self.assertEqual(report["old_noisy_reads_per_trial"],25)
        self.assertEqual(report["private_new_y_uses"],0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(OneBitRealNetworkTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report = {"round":70,**feedback_certificate(),"scope":"Two independent real link sources, each enlarged by one real reference pair, plus one premeasurement Bob-to-Charlie bit",
              "quantum_theory_derived_from_cognition":False,"automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("one_bit_real_network_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
