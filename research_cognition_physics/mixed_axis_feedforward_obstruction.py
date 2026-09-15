"""Round 65: delayed classical processing cannot replace all feedforward.

An exact minimax gap is proved for a fixed mixed-axis experiment and a real
two-parameter target family. This is not a bound on every real network strategy.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from bipartite_composition import interaction
from certified_intervals import Interval as I, SCALE, sin_interval
from complex_bipartite_closure import local_effect
from independent_source_alignment import conditional_reference, keep_systems
from quantum_interface_audit import ALPHA, IDENTITY
from role_symmetry_and_swap import pauli_word


INPUT_RECORDS = tuple(product((-1, 1), repeat=3))
OUTPUT_RECORDS = tuple(product((-1, 1), repeat=2))


def target_family(x, y):
    if abs(x)+abs(y) > 1+1e-15:
        raise ValueError("Use the physical diamond |x|+|y|<=1.")
    return (np.eye(4)+x*pauli_word("XX")+y*pauli_word("YY"))/4


def mixed_effect(angle, sign, beta=ALPHA):
    # Reference is first, target second. The observable is real symmetric.
    observable = math.cos(angle)*pauli_word("IX")+math.sin(angle)*pauli_word("YY")
    return (np.eye(4)+sign*beta*observable)/2


def mixed_records_by_circuit(density, alignment_visibility=ALPHA):
    result = []
    for flag, a, c in INPUT_RECORDS:
        reference = keep_systems(conditional_reference(flag, alignment_visibility), (0, 2), 3)
        joint = np.kron(density, reference)
        effect_a = embed_operator(mixed_effect(math.pi/4, a), (2, 0), 4)
        effect_c = embed_operator(mixed_effect(math.pi/4, c), (3, 1), 4)
        result.append(np.trace(joint @ effect_a @ effect_c).real/2)
    return np.array(result)


def mixed_records(x, y, gamma=ALPHA, beta=ALPHA):
    return np.array([(1+a*c*beta**2*(x+flag*gamma*y)/2)/8 for flag, a, c in INPUT_RECORDS])


def aligned_target_records(x, y, gamma=ALPHA, beta=ALPHA):
    return np.array([(1+a*c*beta**2*(x+gamma*y)/2)/4 for a, c in OUTPUT_RECORDS])


def optimal_classical_channel(gamma):
    if not 0 <= gamma <= 1:
        raise ValueError("Use alignment visibility in [0,1].")
    # Any output first bit can be uniformly regenerated; parity carries all
    # target information on this family. All records, including r=-1, are kept.
    return np.array([[(1+a_out*c_out*a*c*(1. if flag == 1 else (1-gamma)/(1+gamma)))/4
                      for flag, a, c in INPUT_RECORDS] for a_out, c_out in OUTPUT_RECORDS])


def exact_minimax_gap(gamma, beta):
    return beta*beta*gamma/(4*(1+gamma))


def conjugate_axis_classification(axis):
    axis = np.asarray(axis, dtype=float)
    if axis.shape != (3,) or not np.isclose(np.linalg.norm(axis), 1):
        raise ValueError("Use a unit Bloch direction.")
    if axis[1] == 0:
        return "same_outcomes"
    if axis[0] == 0 and axis[2] == 0:
        return "flipped_outcomes"
    return "cannot_correct_by_binary_outcome_postprocessing"


def gap_certificate():
    alpha = 4*sin_interval(I.rational(1, 4))
    gap = alpha*alpha*alpha/(4*(1+alpha))
    return {"lower_exact": str(Fraction(gap.lo, SCALE)), "upper_exact": str(Fraction(gap.hi, SCALE)),
            "interval": gap.floats(), "strictly_greater_than_0_12": Fraction(gap.lo, SCALE) > Fraction(3, 25)}


class MixedAxisFeedforwardObstructionTests(unittest.TestCase):
    def test_mixed_real_effect_is_implemented_by_an_old_yz_flow_and_noisy_x_read(self):
        for angle, sign in product((0., .2, math.pi/4, math.pi/2), (-1, 1)):
            gate = interaction(-angle, "YZ")
            actual = gate.conj().T @ np.kron(IDENTITY, local_effect("X", sign)) @ gate
            np.testing.assert_allclose(actual, mixed_effect(angle, sign), atol=3e-16)

    def test_mixed_local_effects_are_real_positive_complete_measurements(self):
        plus, minus = [mixed_effect(math.pi/4, sign) for sign in (1, -1)]
        np.testing.assert_array_equal(plus.imag, np.zeros((4, 4)))
        np.testing.assert_allclose(plus+minus, np.eye(4), atol=0)
        np.testing.assert_allclose(np.linalg.eigvalsh(plus), [(1-ALPHA)/2]*2+[(1+ALPHA)/2]*2, atol=3e-16)

    def test_actual_full_record_circuit_matches_the_two_parameter_formula(self):
        for x, y in ((1, 0), (-1, 0), (0, 1), (0, -1), (.3, -.4)):
            rho = target_family(x, y)
            self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(), -2e-16)
            for gamma in (.2, ALPHA, 1.):
                np.testing.assert_allclose(mixed_records_by_circuit(rho, gamma), mixed_records(x, y, gamma), atol=2e-16)

    def test_bad_flag_has_identical_data_for_states_requiring_distinct_aligned_outputs(self):
        for gamma in (.2, ALPHA, 1.):
            first = mixed_records(gamma, 0, gamma)
            second = mixed_records(0, -1, gamma)
            indices = [i for i, record in enumerate(INPUT_RECORDS) if record[0] == -1]
            np.testing.assert_allclose(first[indices], second[indices], atol=0)
            desired_first = aligned_target_records(gamma, 0, gamma)
            desired_second = aligned_target_records(0, -1, gamma)
            self.assertAlmostEqual(np.abs(desired_first-desired_second).sum()/2, ALPHA**2*gamma/2, places=14)

    def test_classical_channel_is_stochastic_and_attains_the_exact_minimax_gap(self):
        for gamma, beta in product((.1, .5, ALPHA, 1.), (.3, ALPHA, 1.)):
            channel = optimal_classical_channel(gamma)
            self.assertGreaterEqual(channel.min(), 0)
            np.testing.assert_allclose(channel.sum(axis=0), 1, atol=2e-16)
            errors = []
            for x, y in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                actual = channel @ mixed_records(x, y, gamma, beta)
                expected = aligned_target_records(x, y, gamma, beta)
                errors.append(np.abs(actual-expected).sum()/2)
            np.testing.assert_allclose(errors, exact_minimax_gap(gamma, beta), atol=2e-16)

    def test_exact_dual_inequality_bounds_every_symmetrized_classical_processor(self):
        for gamma in (Fraction(1, 5), Fraction(4, 5), Fraction(1)):
            for plus, minus in product([Fraction(k, 4) for k in range(-4, 5)], repeat=2):
                dx, dy = 2-plus-minus, 2-plus+minus
                self.assertGreaterEqual(dx+dy, 2)
                self.assertGreaterEqual(max(dx, gamma*dy), gamma*(dx+dy)/(1+gamma))
                self.assertGreaterEqual(max(dx, gamma*dy), 2*gamma/(1+gamma))
            plus, minus = Fraction(1), (1-gamma)/(1+gamma)
            self.assertEqual(max(2-plus-minus, gamma*(2-plus+minus)), 2*gamma/(1+gamma))

    def test_pauli_directions_admit_relabeling_but_mixed_y_directions_do_not(self):
        self.assertEqual(conjugate_axis_classification([1, 0, 0]), "same_outcomes")
        self.assertEqual(conjugate_axis_classification([0, 0, 1]), "same_outcomes")
        self.assertEqual(conjugate_axis_classification([0, 1, 0]), "flipped_outcomes")
        self.assertEqual(conjugate_axis_classification([1/math.sqrt(2), 1/math.sqrt(2), 0]), "cannot_correct_by_binary_outcome_postprocessing")

    def test_postselection_of_good_flag_recovers_target_records_at_half_success(self):
        for x, y in ((1, 0), (0, -1), (.2, .5)):
            records = mixed_records(x, y)
            good = np.array([records[i] for i, record in enumerate(INPUT_RECORDS) if record[0] == 1])
            self.assertAlmostEqual(good.sum(), .5, places=14)
            np.testing.assert_allclose(2*good, aligned_target_records(x, y), atol=0)

    def test_old_noise_gives_a_strictly_positive_certified_gap_above_twelve_percent(self):
        certificate = gap_certificate()
        self.assertTrue(certificate["strictly_greater_than_0_12"])
        self.assertLess(Fraction(certificate["upper_exact"]), Fraction(123, 1000))

    def test_exact_alignment_does_not_remove_the_missing_feedforward_gap(self):
        self.assertEqual(exact_minimax_gap(Fraction(1), Fraction(1)), Fraction(1, 8))
        self.assertEqual(exact_minimax_gap(Fraction(0), Fraction(1)), Fraction(0))
        with self.assertRaises(ValueError):
            target_family(.8, .8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MixedAxisFeedforwardObstructionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 65, "target_family": "rho(x,y)=(II+x XX+y YY)/4, |x|+|y|<=1, all real",
              "scope": "Fixed mixed-axis readout; one full record (r,a,c); any state-independent classical stochastic processor; every trial retained",
              "minimax_record_tv_exact": "beta^2 gamma/[4(1+gamma)]",
              "original_noise_gap": gap_certificate(),
              "optimal_parity_retention_coefficients": {"r_plus": "1", "r_minus": "(1-gamma)/(1+gamma)"},
              "postselection_can_remove_gap_with_success_probability": 0.5,
              "bound_on_all_real_quantum_network_strategies": False,
              "two_source_bell_experiment_claimed": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("mixed_axis_feedforward_obstruction_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
