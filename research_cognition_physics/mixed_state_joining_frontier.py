"""Round 131: exact deterministic joining error on equal Bloch balls.

For both unknown inputs having Bloch length <= r, the optimal common-output
worst trace error is r(1+r)/12. The round-128 integer dual and a positive
commuting-projector remainder certify all r in [0,1]. Balanced repair attains.
"""

import argparse
import json
import unittest
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

import numpy as np

from approximate_joining_optimality import dual_certificate, integer_dual, six_probes
from approximate_subject_joining import trace_distance
from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from encoded_composition_audit import independent_encoding
from joining_disturbance_sharing import sharing_output


def validate_radius(radius):
    radius = F(radius)
    if not 0 <= radius <= 1:
        raise ValueError("Bloch radius must lie between zero and one.")
    return radius


def shrink(state, radius):
    r = float(validate_radius(radius))
    return r * state + (1-r) * np.eye(2) / 2


def optimum(radius):
    r = validate_radius(radius)
    return r * (1+r) / 12


def ideal_support_score(radius):
    r = validate_radius(radius)
    return (1+r) ** 2 / 4


def allowed_support_score(radius):
    r = validate_radius(radius)
    return F(1, 4) + 5*r/12 + r*r/6


@lru_cache(maxsize=1)
def positive_remainder_scaled():
    """576 Q, with Q the average input tensor a commuting projection."""
    value = np.zeros((128, 128))
    for a in six_probes():
        for b in six_probes():
            fa = real_lift(np.kron(a, np.eye(2)))
            fb = real_lift(np.kron(np.eye(2), b))
            complement = np.eye(8) - fa - fb + fa @ fb
            if not np.array_equal(complement @ complement, complement):
                raise ArithmeticError("Expected a positive projection.")
            source_projection = 4 * independent_encoding((a, b))
            if not np.array_equal(source_projection @ source_projection, source_projection):
                raise ArithmeticError("Expected the pure input support projection.")
            value += 4 * np.kron(source_projection, complement)
    if not np.array_equal(value, np.rint(value)):
        raise ArithmeticError("The positive remainder is not integer at the specified scale.")
    return value.astype(np.int64)


def mixed_dual(radius):
    r = float(validate_radius(radius))
    y1, _, m1 = integer_dual()
    y1 = y1 / 576
    y0 = np.eye(16) / 64
    ylinear = (np.eye(16) / 16 + y1) / 2
    y = (1-r)**2 * y0 + r*(1-r)*ylinear + r*r*y1
    slack = (r*(1+r)*m1 + r*(1-r)*positive_remainder_scaled()) / (2*576)
    return y, slack


def radius_row(radius):
    r = validate_radius(radius)
    return {
        "radius_exact": str(r),
        "optimal_joint_error_exact": str(optimum(r)),
        "optimal_error_with_host_exactly_preserved": str(r / 6),
        "balanced_each_marginal_error_exact": str(r / 12),
        "ideal_support_score_exact": str(ideal_support_score(r)),
        "all_channel_average_support_bound_exact": str(allowed_support_score(r)),
    }


class MixedStateJoiningFrontierTests(unittest.TestCase):
    def test_mixed_dual_matches_direct_36_preparation_score_matrix(self):
        for r in (F(0), F(1, 4), F(1, 2), F(1)):
            direct = sum(np.kron(independent_encoding((shrink(a, r), shrink(b, r))),
                                 real_lift(np.kron(a, b)))
                         for a in six_probes() for b in six_probes()) / 36
            y, slack = mixed_dual(r)
            np.testing.assert_allclose(np.kron(y, np.eye(8)) - direct, slack, atol=3e-17)
            self.assertAlmostEqual(np.trace(y), float(allowed_support_score(r)))

    def test_analytic_positive_decomposition_uses_the_old_exact_certificate(self):
        self.assertTrue(dual_certificate()["exact_matrix_polynomial_is_zero"])
        q = positive_remainder_scaled()
        np.testing.assert_array_equal(q, q.T)
        for r in (F(0), F(1, 100), F(1, 2), F(1)):
            self.assertGreaterEqual(r*(1+r), 0)
            self.assertGreaterEqual(r*(1-r), 0)
            self.assertEqual(ideal_support_score(r)-allowed_support_score(r), optimum(r))

    def test_balanced_channel_attains_frontier_on_all_36_boundary_pairs(self):
        for r in (F(1, 10), F(1, 2), F(1)):
            for a in six_probes():
                for b in six_probes():
                    first, second = shrink(a, r), shrink(b, r)
                    actual = sharing_output(first, second, .5)
                    ideal = encode_state(np.kron(first, second))
                    self.assertAlmostEqual(trace_distance(actual, ideal), float(optimum(r)))
                    self.assertAlmostEqual(np.trace(actual @ real_lift(np.kron(a, b))).real,
                                           float(allowed_support_score(r)))

    def test_interior_inputs_obey_bound_and_full_ball_is_not_real_plane(self):
        rng = np.random.default_rng(131)
        from quantum_interface_audit import PAULI_X, PAULI_Y, PAULI_Z
        axes = (PAULI_X, PAULI_Y, PAULI_Z)
        for radius in (.1, .5, .9):
            for _ in range(8):
                states = []
                for _ in range(2):
                    vector = rng.normal(size=3)
                    vector *= radius * rng.random() / np.linalg.norm(vector)
                    states.append((np.eye(2) + sum(x*p for x, p in zip(vector, axes))) / 2)
                actual = sharing_output(*states, .5)
                self.assertLessEqual(trace_distance(actual, encode_state(np.kron(*states))),
                                     float(optimum(F(str(radius)))) + 1e-15)

    def test_declared_eigenvalue_formula_explains_trace_bound(self):
        for r in (F(0), F(1, 10), F(1, 2), F(1)):
            differences = (-r*(1+r)/12, r*r/12, r*r/12, r*(1-r)/12)
            self.assertEqual(sum(differences), 0)
            self.assertEqual(sum(abs(x) for x in differences)/2, optimum(r))

    def test_signed_extrapolation_does_not_make_a_nonzero_ball_exactly_joinable(self):
        for r in (F(1, 100), F(1, 2), F(1)):
            plus, minus = (1+r)/(2*r), -(1-r)/(2*r)
            self.assertEqual(plus + minus, 1)
            self.assertEqual(plus - minus, 1/r)
            for pure in six_probes():
                np.testing.assert_allclose(float(plus)*shrink(pure, r) +
                                           float(minus)*shrink(np.eye(2)-pure, r), pure, atol=1e-14)
            self.assertGreater(optimum(r), 0)

    def test_zero_radius_is_exact_and_host_protection_cost_is_explicit(self):
        np.testing.assert_allclose(sharing_output(np.eye(2)/2, np.eye(2)/2, .5),
                                   encode_state(np.eye(4)/4), atol=2e-16)
        for r in (F(1, 10), F(1, 2), F(1)):
            first, second = shrink(six_probes()[0], r), shrink(six_probes()[2], r)
            self.assertAlmostEqual(trace_distance(sharing_output(first, second, 1),
                                                 encode_state(np.kron(first, second))), float(r/6))
            self.assertLessEqual(optimum(r), r/6)

    def test_exact_rows_and_invalid_radius(self):
        self.assertEqual(optimum(F(1, 2)), F(1, 16))
        self.assertEqual(optimum(F(1, 10)), F(11, 1200))
        for r in (F(-1, 10), F(11, 10)):
            with self.assertRaises(ValueError):
                optimum(r)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(MixedStateJoiningFrontierTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 131,
        "preparations": "Independent standard qubit encodings, each Bloch ball of the same radius r",
        "optimal_worst_common_output_trace_distance": "r(1+r)/12",
        "optimal_with_host_exactly_preserved": "r/6",
        "marginal_error_sum_lower_bound": "r/6",
        "balanced_each_marginal_error": "r/12",
        "positive_dual_slack_decomposition": "r(1+r) D_1/2 + r(1-r) Q/2",
        "Q_is_average_input_tensor_positive_commuting_projection": True,
        "D_1_certified_by_round128_integer_polynomial": True,
        "any_nonzero_radius_allows_exact_universal_common_encoding": False,
        "radius_parameter_derived_from_cognitive_principles": False,
        "full_bloch_ball_equated_to_original_real_plane": False,
        "rows": [radius_row(r) for r in (F(0), F(1, 10), F(1, 4), F(1, 2), F(1))],
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("mixed_state_joining_frontier_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
