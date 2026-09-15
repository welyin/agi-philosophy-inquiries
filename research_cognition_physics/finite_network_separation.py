"""Round 69: a finite noisy protocol crosses the certified real-network bound.

Private weak-Y preparation is an added complex-model ability. No trial output
is postselected. Local resource heralding precedes source distribution/settings.
Exact old real control and independent resets remain declared assumptions.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from bell_network_statistics import SETTINGS, bell_basis, bob_effects, diagonal_score, endpoint_observables, network_score, network_table
from bipartite_composition import interaction
from certified_intervals import Interval as I, SCALE
from complex_control_from_reference import real_lift, trace_reference
from noisy_imaginarity_distillation import distillation_certificate, resource_state, updated_bias
from operational_effect_closure import majority_certificate, majority_effect_by_circuit, reusable_majority_probability
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z, rotation_unitary
from real_network_analytic_bound import bound_certificate


def read_unitary(direction):
    x,y,z = np.asarray(direction,dtype=float)
    if not math.isclose(x*x+y*y+z*z,1.,abs_tol=1e-14):
        raise ValueError("Use a unit Bloch direction.")
    theta,phi = math.acos(float(np.clip(z,-1.,1.))),math.atan2(y,x)
    rz_inverse = math.cos(phi/2)*IDENTITY+1j*math.sin(phi/2)*PAULI_Z
    return rotation_unitary(-theta)@rz_inverse


def compiled_axis_gate(direction):
    x,y,z = direction
    theta,phi = math.acos(float(np.clip(z,-1.,1.))),math.atan2(y,x)
    change = np.kron(IDENTITY,rotation_unitary(-math.pi/2))
    yz = change@interaction(-phi,"YX")@change.conj().T
    return np.kron(IDENTITY,rotation_unitary(-theta))@yz


@lru_cache(maxsize=8)
def positive_z_effect(pointer_count):
    return majority_effect_by_circuit(pointer_count)


def axis_effect_by_circuit(direction,bias,pointer_count):
    gate = compiled_axis_gate(direction)
    lifted_effect = gate.conj().T@np.kron(IDENTITY,positive_z_effect(pointer_count))@gate
    return trace_reference(np.kron(resource_state(bias),IDENTITY)@lifted_effect)


def finite_bias(level):
    bias = .2
    for _ in range(level):
        bias,_ = updated_bias(bias,bias,ALPHA)
    return bias


def finite_protocol_table(level=5,pointer_count=5):
    bias = finite_bias(level)
    directions = np.eye(3)
    alice = tuple(2*axis_effect_by_circuit(n,bias,pointer_count)-IDENTITY for n in directions)
    charlie = tuple(2*axis_effect_by_circuit((directions[i]+sign*directions[j])/math.sqrt(2),bias,pointer_count)-IDENTITY
                    for i,j,sign in SETTINGS)
    z = positive_z_effect(pointer_count)
    basis = bell_basis()
    bob = tuple(basis@np.kron(z if b//2==0 else IDENTITY-z,z if b%2==0 else IDENTITY-z)@basis.conj().T for b in range(4))
    return network_table(alice,charlie,bob)


def finite_separation_certificate(level=5,pointer_count=5):
    reference = distillation_certificate(levels=level)["rows"][level]
    pointer = majority_certificate(pointer_count)
    bias_lower = 1-2*Fraction(reference["trace_distance_upper_exact"])
    visibility_lower = 1-2*Fraction(pointer["error_upper_exact"])
    root = I.exact(2).sqrt()
    # Positive polynomial, monotone in the two biases and three visibilities.
    score_lower_interval = 2*root*I.exact(visibility_lower)**2*(visibility_lower*(1+I.exact(bias_lower)**2)+I.exact(visibility_lower)**2)
    score_lower = Fraction(score_lower_interval.lo,SCALE)
    real_upper = Fraction(bound_certificate()["real_upper_rational"])
    margin = score_lower-real_upper
    y_uses_upper = 2*Fraction(reference["expected_new_y_readout_uses_upper_exact"])
    return {"distillation_level":level,"pointer_reads_per_target":pointer_count,
            "initial_private_y_bias":"1/5","reference_trace_distance_upper_exact":reference["trace_distance_upper_exact"],
            "reference_trace_distance_diagnostic":float(Fraction(reference["trace_distance_upper_exact"])),
            "pointer_error_upper_exact":pointer["error_upper_exact"],"pointer_error_diagnostic":pointer["error_diagnostic"],
            "network_score_lower_exact":str(score_lower),"network_score_lower_diagnostic":float(score_lower),
            "real_bound_upper_exact":str(real_upper),"certified_margin_lower_exact":str(margin),
            "certified_margin_lower_diagnostic":float(margin),"strict_separation":margin>0,
            "network_old_pointer_reads":4*pointer_count,
            "finite_success_tree_raw_states_per_outer_lab":2**level,
            "finite_success_tree_parity_reads_per_outer_lab":2**level-1,
            "fixed_two_lab_tree_including_raw_heralds_success_lower_exact":str(Fraction(1,2**(2*(2**(level+1)-1)))),
            "expected_private_new_y_uses_both_labs_upper_exact":str(y_uses_upper),
            "expected_private_new_y_uses_both_labs_diagnostic":float(y_uses_upper),
            "preparation_cost_is_expected_with_retries":True,"bounded_trial_operation_count":True,
            "preparation_is_before_network_sources_and_settings":True,"postselect_network_outputs":False,
            "extra_cross_party_reference_source":False,"old_control_parameters_are_exact":True}


class FiniteNetworkSeparationTests(unittest.TestCase):
    def test_finite_old_gate_compilation_is_the_required_real_lift(self):
        directions = list(np.eye(3))+[(np.eye(3)[i]+s*np.eye(3)[j])/math.sqrt(2) for i,j,s in SETTINGS]
        for direction in directions:
            gate = compiled_axis_gate(direction)
            np.testing.assert_allclose(gate,real_lift(read_unitary(direction)),atol=4e-16)
            np.testing.assert_allclose(gate.imag,0.,atol=0)
            np.testing.assert_allclose(gate.T@gate,np.eye(4),atol=4e-16)
            n = sum(v*p for v,p in zip(direction,(PAULI_X,PAULI_Y,PAULI_Z)))
            u = read_unitary(direction)
            np.testing.assert_allclose(u.conj().T@PAULI_Z@u,n,atol=4e-16)

    def test_actual_private_reference_and_pointer_circuit_give_anisotropic_effects(self):
        directions = list(np.eye(3))+[(np.eye(3)[i]+s*np.eye(3)[j])/math.sqrt(2) for i,j,s in SETTINGS]
        for m in (1,3,5):
            gamma = majority_certificate(m)["effective_visibility_diagnostic"]
            for direction in directions:
                x,y,z = direction
                actual = axis_effect_by_circuit(direction,.7,m)
                expected = (IDENTITY+gamma*(x*PAULI_X+.7*y*PAULI_Y+z*PAULI_Z))/2
                np.testing.assert_allclose(actual,expected,atol=1e-15)

    def test_reused_reset_pointer_matches_the_many_pointer_effect_on_coherent_inputs(self):
        for density in ((IDENTITY+.4*PAULI_X+.3*PAULI_Y+.2*PAULI_Z)/2,
                        (IDENTITY+.7*PAULI_Y-.6*PAULI_Z)/2):
            actual = reusable_majority_probability(density,5)
            self.assertAlmostEqual(actual,np.trace(density@positive_z_effect(5)).real,places=14)

    def test_actual_four_output_bob_decoder_matches_the_bit_noise_formula(self):
        basis,z = bell_basis(),positive_z_effect(5)
        gamma = majority_certificate(5)["effective_visibility_diagnostic"]
        for b,effect in enumerate(bob_effects(gamma)):
            direct = basis@np.kron(z if b//2==0 else IDENTITY-z,z if b%2==0 else IDENTITY-z)@basis.conj().T
            # The independent 32-record pointer expansion accumulates roundoff.
            # Strict separation below uses dyadic bounds, not this comparison.
            np.testing.assert_allclose(effect,direct,atol=3e-15)

    def test_all_288_finite_probabilities_match_the_anisotropic_model(self):
        actual = finite_protocol_table()
        gamma = majority_certificate(5)["effective_visibility_diagnostic"]
        bias = finite_bias(5)
        a,c = endpoint_observables(gamma,bias)
        expected = network_table(a,c,bob_effects(gamma))
        self.assertEqual(len(actual),288)
        np.testing.assert_allclose(list(actual.values()),list(expected.values()),atol=3e-16)
        self.assertAlmostEqual(network_score(actual),diagonal_score(gamma,gamma,gamma,bias,bias),places=12)

    def test_strict_rational_margin_survives_finite_distillation_and_original_readouts(self):
        report = finite_separation_certificate()
        self.assertTrue(report["strict_separation"])
        self.assertGreater(Fraction(report["certified_margin_lower_exact"]),Fraction(38,10000))
        self.assertGreater(Fraction(report["reference_trace_distance_upper_exact"]),0)
        self.assertGreater(Fraction(report["pointer_error_upper_exact"]),0)
        self.assertAlmostEqual(network_score(finite_protocol_table()),report["network_score_lower_diagnostic"],places=12)

    def test_weaker_resources_do_not_cross_this_conservative_bound(self):
        self.assertFalse(finite_separation_certificate(level=4,pointer_count=5)["strict_separation"])
        self.assertFalse(finite_separation_certificate(level=5,pointer_count=1)["strict_separation"])

    def test_pretrial_heralding_preserves_the_source_product_and_has_positive_finite_success(self):
        # Local readiness events act on private resources only. Conditioning
        # multiplies two independent preparations, not a function of Bob's trial b.
        left,right = resource_state(.4),resource_state(.7)
        pl,pr = Fraction(2,5),Fraction(3,7)
        np.testing.assert_allclose(np.kron(float(pl)*left,float(pr)*right)/float(pl*pr),np.kron(left,right),atol=6e-17)
        report = finite_separation_certificate()
        self.assertEqual(report["network_old_pointer_reads"],20)
        self.assertEqual(report["finite_success_tree_raw_states_per_outer_lab"],32)
        self.assertGreater(Fraction(report["fixed_two_lab_tree_including_raw_heralds_success_lower_exact"]),0)
        self.assertFalse(report["postselect_network_outputs"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FiniteNetworkSeparationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round":69,**finite_separation_certificate(),"actual_circuit_score_diagnostic":network_score(finite_protocol_table()),
              "scope":"Conditional complex extension with added weak-Y preparation, exact old real gates and independent local resets; finite operations after pretrial resource heralding",
              "empirical_experiment_performed":False,"quantum_theory_derived_from_cognition":False,
              "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("finite_network_separation_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
