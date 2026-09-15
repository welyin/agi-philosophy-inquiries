"""Finite stopping protocols and rigorous integer probability certificates."""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from outcome_updates import WINDOW_ATTENUATION as ALPHA
from retention_tradeoff import fresh_branch_matrix


def hidden_branches(contrast, alpha=ALPHA):
    flip = (1 - alpha) / 2
    transition = np.array([[1 - flip, flip], [flip, 1 - flip]])
    return [transition @ np.diag([(1 + sign * contrast) / 2, (1 - sign * contrast) / 2]) for sign in (-1, 1)]


def stopped_probabilities(contrast, boundary, cap, alpha=ALPHA):
    """Start from the positive bit. Timeout outputs a fair coin."""
    if boundary < 1 or cap < 1:
        raise ValueError("Use positive count boundary and finite cap.")
    width = 2 * boundary - 1
    active = np.zeros((width, 2))
    active[boundary - 1, 0] = 1
    hit_minus = hit_plus = 0.0
    branches = hidden_branches(contrast, alpha)
    for _ in range(cap):
        minus, plus = [active @ branch.T for branch in branches]
        hit_minus += minus[0].sum()
        hit_plus += plus[-1].sum()
        updated = np.zeros_like(active)
        updated[:-1] += minus[1:]
        updated[1:] += plus[:-1]
        active = updated
    return float(hit_plus), float(hit_minus), float(active.sum())


def alpha_interval():
    lower = 1 - Fraction(1, 96) + Fraction(1, 30720) - Fraction(1, 4096 * 5040)
    upper = lower + Fraction(1, 65536 * 362880)
    return lower, upper


def certified_gain(contrast, boundary, cap, bits=96):
    """Lower integer masses for a rational alpha, then a coupling error bound.

    Every nonnegative probability multiplication rounds downward. The missing
    total mass bounds all lost probability, including intermediate rounding.
    """
    contrast = Fraction(contrast)
    if not 0 <= contrast <= 1 or boundary < 1 or cap < 1:
        raise ValueError("Use a rational contrast, positive boundary and finite cap.")
    alpha_low, alpha_high = alpha_interval()
    denominator = 1 << bits
    floor_scaled = lambda value: value.numerator * denominator // value.denominator
    flip = (1 - alpha_low) / 2
    transition = [[1 - flip, flip], [flip, 1 - flip]]
    emissions = [[(1 + sign * contrast) / 2, (1 - sign * contrast) / 2] for sign in (-1, 1)]
    emission_lower = [[floor_scaled(value) for value in row] for row in emissions]
    branch_lower = [[[floor_scaled(transition[new][old] * emissions[mark][old]) for old in range(2)]
                     for new in range(2)] for mark in range(2)]
    width = 2 * boundary - 1
    active = [[0, 0] for _ in range(width)]
    active[boundary - 1][0] = denominator
    hits = [0, 0]
    for _ in range(cap):
        updated = [[0, 0] for _ in range(width)]
        for index, masses in enumerate(active):
            for mark, step in enumerate((-1, 1)):
                destination = index + step
                if destination < 0 or destination >= width:
                    hits[mark] += sum(masses[old] * emission_lower[mark][old] // denominator for old in range(2))
                else:
                    for new in range(2):
                        updated[destination][new] += sum(masses[old] * branch_lower[mark][new][old] // denominator for old in range(2))
        active = updated
    active_mass = sum(map(sum, active))
    lost_mass = denominator - sum(hits) - active_mass
    if lost_mass < 0:
        raise ArithmeticError("Lower probability masses cannot exceed one.")
    rational_gain_lower = Fraction(hits[1] - hits[0] - lost_mass, denominator)
    # At most cap bit-flips can differ under a maximal coupling. The signed
    # output expectation changes by at most cap*|alpha-alpha_low|.
    coupling_error = cap * (alpha_high - alpha_low)
    actual_gain_lower = rational_gain_lower - coupling_error
    return {"gain_lower": actual_gain_lower, "visibility_lower": alpha_low * actual_gain_lower,
            "rounding_mass_bound": Fraction(lost_mass, denominator), "coupling_gain_error_bound": coupling_error}


def pi_lower_bound():
    def arctangent_bounds(inverse, terms):
        total = sum(((-1)**index * Fraction(1, (2 * index + 1) * inverse**(2 * index + 1)) for index in range(terms)), Fraction(0))
        remainder = Fraction(1, (2 * terms + 1) * inverse**(2 * terms + 1))
        return (total, total + remainder) if terms % 2 == 0 else (total - remainder, total)
    lower_five, _ = arctangent_bounds(5, 12)
    _, upper_239 = arctangent_bounds(239, 4)
    return 16 * lower_five - 4 * upper_239


def certified_witnesses():
    rows = []
    for eta, boundary in ((Fraction(17, 100), 8), (Fraction(11, 50), 7)):
        certificate = certified_gain(eta, boundary, 512)
        visibility = certificate["visibility_lower"]
        plus, minus, timeout = stopped_probabilities(float(eta), boundary, 512)
        rows.append({"contrast": float(eta), "boundary": boundary, "cap": 512,
                     "gain_lower": float(certificate["gain_lower"]), "visibility_lower": float(visibility),
                     "visibility_lower_exact": str(visibility), "two_over_pi_rational_upper": str(2 / pi_lower_bound()),
                     "exceeds_entire_disk_bound_certified": bool(visibility > 2 / pi_lower_bound()),
                     "exceeds_four_preparation_bound_certified": bool(visibility > 0 and visibility**2 > Fraction(1, 2)),
                     "rounding_mass_bound": float(certificate["rounding_mass_bound"]),
                     "coupling_gain_error_bound": float(certificate["coupling_gain_error_bound"]),
                     "floating_score": (1 + ALPHA * (plus - minus) / math.sqrt(2)) / 2,
                     "timeout": timeout})
    return rows


def scan():
    rows = []
    for eta in (0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3):
        candidates = []
        for boundary in range(1, 21):
            plus, minus, timeout = stopped_probabilities(eta, boundary, 512)
            gain = plus - minus
            candidates.append((gain, boundary, timeout))
        gain, boundary, timeout = max(candidates)
        rows.append({"contrast": eta, "boundary": boundary, "cap": 512, "gain": gain,
                     "visibility": ALPHA * gain, "four_preparation_score": (1 + ALPHA * gain / math.sqrt(2)) / 2,
                     "timeout": timeout})
    return rows


class StoppedReadoutWitnessTests(unittest.TestCase):
    def test_finite_protocol_at_point_seventeen_excludes_entire_disk_models(self):
        certificate = certified_gain(Fraction(17, 100), 8, 512)
        self.assertGreater(certificate["visibility_lower"], 2 / pi_lower_bound())

    def test_point_twenty_two_has_a_rational_four_preparation_violation(self):
        certificate = certified_gain(Fraction(11, 50), 7, 512)
        self.assertGreater(certificate["visibility_lower"]**2, Fraction(1, 2))

    def test_machin_bound_is_consistent_with_reference_pi(self):
        bound = pi_lower_bound()
        self.assertGreater(bound, Fraction(314159, 100000))
        self.assertLessEqual(float(bound), math.pi)

    def test_hidden_branches_equal_existing_predictive_block(self):
        transform = np.array([[1, 1], [1, -1]])
        for eta in (0.1, 0.2, 0.5):
            for mark, branch in enumerate(hidden_branches(eta)):
                np.testing.assert_allclose(transform @ branch @ transform / 2,
                                           fresh_branch_matrix(0, mark, eta)[:2, :2], atol=1e-14)

    def test_stopping_and_timeout_probabilities_normalize(self):
        for eta, boundary, cap in ((0, 2, 7), (0.2, 4, 128), (0.1, 20, 512)):
            probabilities = stopped_probabilities(eta, boundary, cap)
            self.assertGreaterEqual(min(probabilities), 0)
            self.assertAlmostEqual(sum(probabilities), 1)

    def test_boundary_one_is_exactly_one_read(self):
        for eta in (0, 0.2, 0.5):
            plus, minus, timeout = stopped_probabilities(eta, 1, 8)
            self.assertAlmostEqual(plus - minus, eta)
            self.assertAlmostEqual(timeout, 0)

    def test_finite_stopping_rule_matches_word_enumeration(self):
        eta, boundary, cap = 0.2, 2, 6
        probabilities = np.zeros(3)
        matrices = hidden_branches(eta)
        for marks in product((0, 1), repeat=cap):
            state = np.array([1.0, 0.0])
            count, decided = 0, None
            for mark in marks:
                state = matrices[mark] @ state
                count += 2 * mark - 1
                if decided is None and abs(count) == boundary:
                    decided = int(count > 0)
            probabilities[2 if decided is None else decided] += state.sum()
        plus, minus, timeout = stopped_probabilities(eta, boundary, cap)
        np.testing.assert_allclose(probabilities, [minus, plus, timeout], atol=1e-14)

    def test_alpha_series_encloses_the_original_parameter(self):
        lower, upper = alpha_interval()
        self.assertLess(float(lower), ALPHA)
        self.assertGreater(float(upper), ALPHA)

    def test_downward_rounding_certificate_is_conservative(self):
        eta, boundary, cap = Fraction(1, 5), 4, 128
        certificate = certified_gain(eta, boundary, cap)
        plus, minus, _ = stopped_probabilities(float(eta), boundary, cap)
        self.assertLess(float(certificate["gain_lower"]), plus - minus)
        self.assertLess(plus - minus - float(certificate["gain_lower"]), 1e-7)
        self.assertLess(float(certificate["rounding_mass_bound"]), 1e-20)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(StoppedReadoutWitnessTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    rows = scan()
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "certified_witnesses": certified_witnesses(),
               "bounded_stopping_rule_scan_not_global_optimization": rows,
               "scope": "Finite stopping at absolute accumulated count boundary, at most 512 reads and a fair timeout output. Contrast 0.17 excludes preparation-noncontextual models of the entire disk by the 2/pi bound; contrast 0.22 also violates the finite four-preparation inequality. Integer downward rounding, a rational sinc interval and a coupling bound certify the exclusions. This is not optimal over all sequential strategies and does not exclude contextual classical models."}
    if args.write_results:
        Path(__file__).with_name("stopped_readout_witness_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
