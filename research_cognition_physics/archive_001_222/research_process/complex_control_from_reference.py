"""Complex qubit control from real SO(4) control and a Y reference.

Ideal reference reuse is catalytic. An imperfect reference retains a shared
conjugation branch, not independent gate noise. Existing assumptions remain
explicit; this is not a derivation of quantum mechanics from cognition.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from ancillary_operation_equivalence import J, apply_operation
from noisy_imaginarity_distillation import distillation_certificate, resource_state
from operational_effect_closure import majority_certificate, majority_effect_by_circuit
from quantum_interface_audit import (ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                     branch_kraus, effective_readout)


def real_lift(matrix):
    matrix = np.asarray(matrix, dtype=complex)
    return np.kron(np.eye(2), matrix.real) - np.kron(J, matrix.imag)


def trace_reference(joint):
    dimension = joint.shape[0]//2
    return np.einsum("abad->bd", joint.reshape(2, dimension, 2, dimension))


def trace_system(joint):
    dimension = joint.shape[0]//2
    return np.einsum("abcb->ac", joint.reshape(2, dimension, 2, dimension))


def lifted_sequence(density, bias, gates):
    joint = np.kron(resource_state(bias), density)
    for gate in gates:
        lifted = real_lift(gate)
        joint = lifted @ joint @ lifted.T
    return joint


def conjugate_branch_output(density, bias, gate):
    first = gate @ density @ gate.conj().T
    second = gate.conj() @ density @ gate.T
    return (1+bias)*first/2 + (1-bias)*second/2


def unitary_to_vector(vector):
    a, b = vector
    return np.array([[a, -b.conjugate()], [b, a.conjugate()]])


def imaginary_rotation(angle):
    return math.cos(angle/2)*IDENTITY - 1j*math.sin(angle/2)*PAULI_X


def effective_y_from_reference(bias, z_effect):
    gate = real_lift(imaginary_rotation(math.pi/2))
    joint_effect = gate.T @ np.kron(IDENTITY, z_effect) @ gate
    return trace_reference(np.kron(resource_state(bias), IDENTITY) @ joint_effect)


def end_to_end_readout_certificate(level=9, pointer_count=15):
    reference = distillation_certificate(levels=level)["rows"][level]
    pointer = majority_certificate(pointer_count)
    reference_error = Fraction(reference["trace_distance_upper_exact"])
    pointer_error = Fraction(pointer["error_upper_exact"])
    error = reference_error + pointer_error - 2*reference_error*pointer_error
    return {"initial_added_y_visibility": "1/5", "distillation_level": level,
            "original_z_pointer_read_count": pointer_count,
            "reference_error_upper_exact": str(reference_error),
            "pointer_error_upper_exact": str(pointer_error),
            "final_y_effect_error_upper_exact": str(error),
            "final_y_effect_error_upper_diagnostic": float(error),
            "strictly_below_1e_minus_12": error < Fraction(1, 10**12),
            "scope": "Exact lifted real control, conditioned reference preparation, old reset/pointer amplification; finite control compilation error is additional"}


class ComplexControlFromReferenceTests(unittest.TestCase):
    def test_real_lift_preserves_products_adjoint_and_identity(self):
        rng = np.random.default_rng(56)
        for dimension in (2, 3, 4):
            a = rng.normal(size=(dimension, dimension)) + 1j*rng.normal(size=(dimension, dimension))
            b = rng.normal(size=(dimension, dimension)) + 1j*rng.normal(size=(dimension, dimension))
            np.testing.assert_allclose(real_lift(a @ b), real_lift(a) @ real_lift(b), atol=4e-15)
            np.testing.assert_array_equal(real_lift(a.conj().T), real_lift(a).T)
            np.testing.assert_array_equal(real_lift(np.eye(dimension)), np.eye(2*dimension))

    def test_complex_unitaries_lift_to_proper_real_orthogonal_matrices(self):
        rng = np.random.default_rng(156)
        for dimension in (2, 3, 4):
            for _ in range(6):
                unitary, _ = np.linalg.qr(rng.normal(size=(dimension, dimension)) + 1j*rng.normal(size=(dimension, dimension)))
                lifted = real_lift(unitary)
                self.assertTrue(np.isrealobj(lifted))
                np.testing.assert_allclose(lifted.T @ lifted, np.eye(2*dimension), atol=1e-15)
                self.assertAlmostEqual(np.linalg.det(lifted), 1., places=13)

    def test_pure_reference_is_returned_as_a_product_even_with_entangled_spectators(self):
        vector = np.array([1, 0, 0, 1j])/math.sqrt(2)
        density = np.outer(vector, vector.conj())
        gate = np.kron(imaginary_rotation(.7), IDENTITY)
        actual = lifted_sequence(density, 1., [gate])
        expected = np.kron(resource_state(1.), gate @ density @ gate.conj().T)
        np.testing.assert_allclose(actual, expected, atol=2e-16)

    def test_mixed_reference_has_two_conjugate_branches_and_unchanged_marginal(self):
        density = (IDENTITY + .4*PAULI_X + .3*PAULI_Y)/2
        for bias in (0., .2, .9, 1.):
            gate = imaginary_rotation(.8)
            actual = lifted_sequence(density, bias, [gate])
            expected = ((1+bias)/2*np.kron(resource_state(1.), gate @ density @ gate.conj().T)
                        +(1-bias)/2*np.kron(resource_state(-1.), gate.conj() @ density @ gate.T))
            np.testing.assert_allclose(actual, expected, atol=3e-16)
            np.testing.assert_allclose(trace_system(actual), resource_state(bias), atol=3e-16)

    def test_complete_complex_single_pure_state_preparation_in_the_reference_limit(self):
        rng = np.random.default_rng(256)
        initial = np.diag([1., 0.])
        for _ in range(12):
            vector = rng.normal(size=2) + 1j*rng.normal(size=2)
            vector /= np.linalg.norm(vector)
            gate = unitary_to_vector(vector)
            for bias in (.2, .99, 1.):
                output = trace_reference(lifted_sequence(initial, bias, [gate]))
                target = np.outer(vector, vector.conj())
                self.assertLessEqual(np.linalg.norm(output-target, ord="nuc")/2, (1-bias)/2 + 4e-16)
                if bias == 1:
                    np.testing.assert_allclose(output, target, atol=3e-16)

    def test_reference_error_bound_does_not_grow_with_sequence_length(self):
        rng = np.random.default_rng(356)
        density = (IDENTITY + .2*PAULI_X + .7*PAULI_Y)/2
        for length, bias in product((1, 2, 10, 50), (.2, .9)):
            gates = []
            combined = IDENTITY.copy()
            for _ in range(length):
                gate, _ = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))
                gates.append(gate)
                combined = gate @ combined
            output = trace_reference(lifted_sequence(density, bias, gates))
            np.testing.assert_allclose(output, conjugate_branch_output(density, bias, combined), atol=2e-14)
            self.assertLessEqual(np.linalg.norm(output-combined @ density @ combined.conj().T, ord="nuc")/2, (1-bias)/2 + 2e-14)

    def test_shared_reference_and_reset_references_have_different_two_gate_statistics(self):
        phase = np.diag([1., 1j])
        initial = (IDENTITY + PAULI_X)/2
        for bias in (.2, .8):
            reused = trace_reference(lifted_sequence(initial, bias, [phase, phase]))
            reset = conjugate_branch_output(conjugate_branch_output(initial, bias, phase), bias, phase)
            np.testing.assert_allclose(reused, (IDENTITY-PAULI_X)/2, atol=2e-16)
            np.testing.assert_allclose(reset, (IDENTITY-bias*bias*PAULI_X)/2, atol=3e-16)

    def test_unchanged_reference_marginal_does_not_mean_no_correlations(self):
        phase = np.diag([1., 1j])
        initial = (IDENTITY + PAULI_X)/2
        for bias in (0., .2, .9):
            joint = lifted_sequence(initial, bias, [phase])
            product_marginals = np.kron(trace_system(joint), trace_reference(joint))
            self.assertAlmostEqual(np.linalg.norm(joint-product_marginals, ord="nuc")/2, (1-bias*bias)/2, places=14)

    def test_full_adaptive_record_distribution_keeps_the_same_global_branch_weight(self):
        bias, initial = .4, (IDENTITY + .5*PAULI_Y + .2*PAULI_Z)/2
        observed, ideal, conjugated = {}, {}, {}
        final_joint = np.zeros((4, 4), dtype=complex)
        for record in product((0, 1), repeat=4):
            joint = np.kron(resource_state(bias), initial)
            states = [initial.copy(), initial.copy()]
            for step, outcome in enumerate(record):
                gate = imaginary_rotation(.3 + .17*step + .23*sum(record[:step]))
                lifted = real_lift(gate)
                joint = lifted @ joint @ lifted.T
                for index, candidate in enumerate((gate, gate.conj())):
                    states[index] = candidate @ states[index] @ candidate.conj().T
                kraus = branch_kraus(.21*step, outcome, (.2, .7, .4, .8)[step])
                joint = apply_operation(tuple(np.kron(IDENTITY, k) for k in kraus), joint)
                states = [apply_operation(kraus, state) for state in states]
            observed[record] = np.trace(joint).real
            ideal[record], conjugated[record] = [np.trace(state).real for state in states]
            self.assertAlmostEqual(observed[record], (1+bias)*ideal[record]/2 + (1-bias)*conjugated[record]/2, places=14)
            final_joint += joint
        self.assertAlmostEqual(sum(observed.values()), 1., places=14)
        distance = sum(abs(observed[r]-ideal[r]) for r in observed)/2
        self.assertLessEqual(distance, (1-bias)/2)
        np.testing.assert_allclose(trace_system(final_joint), resource_state(bias), atol=5e-16)

    def test_real_lifts_of_complete_complex_kraus_maps_remain_complete(self):
        rate = .3
        damp = (np.diag([1., math.sqrt(1-rate)]), np.array([[0, math.sqrt(rate)], [0, 0]]))
        rotation = imaginary_rotation(.7)
        kraus = tuple(k @ rotation for k in damp)
        lifts = tuple(real_lift(k) for k in kraus)
        np.testing.assert_allclose(sum(k.T @ k for k in lifts), np.eye(4), atol=3e-16)
        density = (IDENTITY + .3*PAULI_X)/2
        actual = apply_operation(lifts, np.kron(resource_state(1.), density))
        np.testing.assert_allclose(actual, np.kron(resource_state(1.), apply_operation(kraus, density)), atol=2e-16)

    def test_actual_old_readout_and_pointer_circuit_gain_a_y_direction(self):
        for bias in (.2, .9, 1.):
            effect = effective_y_from_reference(bias, effective_readout(0, 1, 1.))
            np.testing.assert_allclose(effect, (IDENTITY + ALPHA*bias*PAULI_Y)/2, atol=3e-16)
        for count in (1, 3):
            z_effect = majority_effect_by_circuit(count)
            gamma = majority_certificate(count)["effective_visibility_diagnostic"]
            y_effect = effective_y_from_reference(.9, z_effect)
            np.testing.assert_allclose(y_effect, (IDENTITY+.9*gamma*PAULI_Y)/2, atol=7e-16)

    def test_end_to_end_ideal_y_approximation_has_a_strict_finite_resource_certificate(self):
        report = end_to_end_readout_certificate()
        self.assertTrue(report["strictly_below_1e_minus_12"])
        self.assertGreater(Fraction(report["final_y_effect_error_upper_exact"]), 0)
        self.assertLess(Fraction(report["final_y_effect_error_upper_exact"]), Fraction(1, 10**12))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ComplexControlFromReferenceTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"scope": "Y reference from the explicitly added measurement, old SO(4) control closure, usual complex tensor continuation of the enlarged state set",
              "real_lift_formula": "R(U)=I_ref tensor Re(U) - J_ref tensor Im(U), J=-iY",
              "every_complex_qubit_unitary_lifts_to_SO4": True,
              "pure_reference_returned_uncorrelated": True,
              "imperfect_reference_returned_uncorrelated_in_general": False,
              "same_reference_whole_program_trace_or_record_distance_bound": "(1-nu)/2, independent of depth for exact lifts",
              "same_reference_noise_equals_independent_gate_noise": False,
              "finite_gate_compilation_errors_are_additional": True,
              "general_lifted_kraus_formula_implies_all_such_real_maps_available": False,
              "closed_single_state_set_after_weak_y_resource_distillation": "Full complex Bloch ball",
              "end_to_end_y_readout_certificate": end_to_end_readout_certificate(),
              "new_y_readout_derived_from_cognition": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("complex_control_from_reference_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
