"""A continuous two-tail martingale construction at the exact PN endpoint.

The stationary scalar part is a continuum of point masses, so this is not a
finite intermediate-label factorization. See research_note_45.md for the proof.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from adaptive_arc_bridge import endpoint_bracket
from certified_intervals import Interval as I, SCALE, pi_interval, sin_interval, cos_interval
from compensated_angle_kernel import preconvolution_density
from outcome_updates import WINDOW_ATTENUATION as ALPHA


CROSSING_UPPER = Fraction(521, 1000)
UPPER_TAIL_MEAN_X_LOWER = Fraction(1993, 2000)  # .9965
UPPER_TAIL_MEAN_Y_LOWER = Fraction(719, 10000)  # .0719
LOWER_TAIL_MEAN_Y_LOWER = Fraction(1707, 2000)  # .8535
UPPER_TAIL_ANGLE_UPPER = Fraction(1443, 10000)  # .1443


@lru_cache(maxsize=1)
def interval_parameters():
    lower, upper = endpoint_bracket()
    eta = I(I.exact(lower).lo, I.exact(upper).hi)
    alpha = 4 * sin_interval(I.rational(1, 4))
    return alpha, eta, (1 - eta**2).sqrt()


def density_factor_interval(s, parameters=None):
    alpha, eta, q = interval_parameters() if parameters is None else parameters
    v2 = (eta * s / alpha)**2
    return q**3 * (1 + v2) / (1 - v2)**2


@lru_cache(maxsize=1)
def constants_certificate():
    alpha, eta, q = interval_parameters()
    c = I.exact(CROSSING_UPPER)
    half_axis_factor = q**3 * (1 + eta**2 / 2) / (1 - eta**2 / 2)**2
    margins = {
        "alpha_positive": alpha,
        "alpha_below_one": 1 - alpha,
        "eta_positive": eta,
        "eta_below_one": 1 - eta,
        "source_density_smaller_at_zero": alpha - q**3,
        "crossing_before_0.521": density_factor_interval(c) * (1 - c**2).sqrt() - (alpha**2 - c**2).sqrt(),
        "crossing_bound_below_alpha_over_sqrt_two": alpha / I.exact(2).sqrt() - c,
        "upper_tail_angle_below_0.1443": alpha - cos_interval(I.exact(UPPER_TAIL_ANGLE_UPPER)),
        "upper_tail_mean_x_above_0.9965": 1 - I.exact(UPPER_TAIL_ANGLE_UPPER)**2 / 6 - UPPER_TAIL_MEAN_X_LOWER,
        "upper_tail_mean_x_bound_above_alpha": UPPER_TAIL_MEAN_X_LOWER - alpha,
        "upper_tail_mean_y_above_0.0719": (1 - alpha) / UPPER_TAIL_ANGLE_UPPER - UPPER_TAIL_MEAN_Y_LOWER,
        "lower_tail_mean_y_above_0.8535": 1 - c**2 - LOWER_TAIL_MEAN_Y_LOWER**2,
        "factor_below_one_beyond_quarter_angle": 1 - half_axis_factor,
    }
    return {"verified": all(value.lo > 0 for value in margins.values()),
            "strict_margin_lower_bounds_exact": {name: str(Fraction(value.lo, SCALE)) for name, value in margins.items()},
            "strict_margin_diagnostics": {name: value.lo / SCALE for name, value in margins.items()},
            "crossing_upper": str(CROSSING_UPPER),
            "upper_tail_mean_x_lower": str(UPPER_TAIL_MEAN_X_LOWER),
            "upper_tail_mean_y_lower": str(UPPER_TAIL_MEAN_Y_LOWER),
            "lower_tail_mean_y_lower": str(LOWER_TAIL_MEAN_Y_LOWER)}


def line_lower_interval(s):
    u, v, ell = UPPER_TAIL_MEAN_X_LOWER, UPPER_TAIL_MEAN_Y_LOWER, LOWER_TAIL_MEAN_Y_LOWER
    return v + (ell - v) * (u - s) / u


def angular_cell_margin(index, cells):
    """Enclose the complete angle interval, not only its center or endpoints."""
    if type(cells) is not int or cells < 1 or type(index) is not int or not 0 <= index < cells:
        raise ValueError("Invalid angular interval.")
    pi = pi_interval()
    lower, upper = pi * Fraction(index, 4 * cells), pi * Fraction(index + 1, 4 * cells)
    theta = I(lower.lo, upper.hi)
    alpha, eta, q = interval_parameters()
    cosine, sine = cos_interval(theta), sin_interval(theta)
    s, height = alpha * cosine, alpha * sine
    outer_height = (1 - s**2).sqrt()
    v2 = eta**2 * cosine**2
    factor = q**3 * (1 + v2) / (1 - v2)**2
    ratio = height / outer_height
    return height * (1 - factor) + (factor - ratio) * line_lower_interval(s)


@lru_cache(maxsize=2)
def verify_endpoint(cells=4096):
    constants = constants_certificate()
    margins = [angular_cell_margin(index, cells) for index in range(cells)]
    index = min(range(cells), key=lambda j: margins[j].lo)
    minimum = margins[index].lo
    lower, upper = endpoint_bracket()
    return {"verified": constants["verified"] and minimum > 0,
            "endpoint_definition": "Unique positive root of alpha*(sqrt(1-eta^2)+eta*asin(eta))=1",
            "endpoint_bracket_exact": [str(lower), str(upper)],
            "constants": constants, "angular_cells": cells, "certified_angle_range": "[0, pi/4]",
            "minimum_angle_margin_lower_bound_exact": str(Fraction(minimum, SCALE)),
            "minimum_angle_margin_lower_bound_diagnostic": minimum / SCALE,
            "weakest_cell_index": index,
            "rest_of_source": "At theta>=pi/4 the density factor is below one; the overlap region is automatic.",
            "scope": "The coupling uses the exact endpoint root and its exact first-moment equality, not the entire bracket. The stationary scalar component is not a finite-label factorization."}


def density_factor_float(s, eta):
    v2 = (eta * np.asarray(s) / ALPHA)**2
    return (1 - eta**2)**1.5 * (1 + v2) / (1 - v2)**2


@lru_cache(maxsize=1)
def tail_model_float():
    eta = float(sum(endpoint_bracket()) / 2)
    def density_difference(s):
        outer = 2 / (math.pi * np.sqrt(1 - s**2))
        inner = 2 * density_factor_float(s, eta) / (math.pi * np.sqrt(ALPHA**2 - s**2))
        return outer - inner
    lower, upper = 0., ALPHA
    for _ in range(60):
        middle = (lower + upper) / 2
        if density_difference(middle) > 0:
            lower = middle
        else:
            upper = middle
    crossing = (lower + upper) / 2
    nodes, weights = np.polynomial.legendre.leggauss(128)
    s = (nodes + 1) * crossing / 2
    lower_weights = weights * crossing / 2 * density_difference(s)
    lower_mass = lower_weights.sum()
    lower_x = lower_weights @ s / lower_mass
    lower_y = lower_weights @ np.sqrt(1 - s**2) / lower_mass
    angle = math.acos(ALPHA)
    upper_mass = 2 * angle / math.pi
    upper_x = math.sin(angle) / angle
    upper_y = (1 - ALPHA) / angle
    return {"eta": eta, "crossing": crossing, "lower_mass": float(lower_mass), "upper_mass": upper_mass,
            "lower_mean_x": float(lower_x), "lower_mean_y": float(lower_y),
            "upper_mean_x": upper_x, "upper_mean_y": upper_y,
            "stationary_mass": 1 - float(lower_mass) - upper_mass,
            "lower_nodes": s, "lower_weights": lower_weights / lower_mass,
            "upper_nodes": np.cos((nodes + 1) * angle / 2), "upper_weights": weights / 2}


def scalar_conditional_float(s, inner_height=None):
    model = tail_model_float()
    s = np.asarray(s)
    height = np.sqrt(np.maximum(0, ALPHA**2 - s**2)) if inner_height is None else np.asarray(inner_height)
    g = np.sqrt(1 - s**2)
    stay = np.minimum(1, height / (density_factor_float(s, model["eta"]) * g))
    excess = s > model["crossing"]
    lower = np.where(excess, (model["upper_mean_x"] - s) / (model["upper_mean_x"] - model["lower_mean_x"]), 0)
    upper = np.where(excess, 1 - lower, 0)
    # Outside the excess region both tail coefficients are irrelevant and zero.
    stay = np.where(excess, stay, 1)
    expected_x = stay * s + (1 - stay) * (lower * model["lower_mean_x"] + upper * model["upper_mean_x"])
    expected_height = stay * g + (1 - stay) * (lower * model["lower_mean_y"] + upper * model["upper_mean_y"])
    bias = height / expected_height
    return {"stay": stay, "lower_tail": (1 - stay) * lower, "upper_tail": (1 - stay) * upper,
            "expected_x": expected_x, "expected_height": expected_height, "vertical_bias": bias,
            "height_margin": expected_height - height}


def continuous_conditional_mean(angles):
    angles = np.asarray(angles)
    x, y = ALPHA * np.cos(angles), ALPHA * np.sin(angles)
    row = scalar_conditional_float(np.abs(x), np.abs(y))
    return np.column_stack((np.sign(x) * row["expected_x"], np.sign(y) * row["vertical_bias"] * row["expected_height"]))


def source_quadrature():
    model = tail_model_float()
    nodes, weights = np.polynomial.legendre.leggauss(160)
    crossing_angle = math.acos(model["crossing"] / ALPHA)
    angles, mass = [], []
    for lower, upper in ((0., crossing_angle), (crossing_angle, math.pi / 2)):
        theta = lower + (nodes + 1) * (upper - lower) / 2
        angles.extend(theta)
        mass.extend(weights * (upper - lower) / math.pi * preconvolution_density(theta, model["eta"]))
    return np.array(angles), np.array(mass)


class EndpointTailBridgeTests(unittest.TestCase):
    def test_exact_endpoint_constants_and_entire_angular_intervals_are_positive(self):
        report = verify_endpoint()
        self.assertTrue(report["verified"])
        self.assertGreater(Fraction(report["minimum_angle_margin_lower_bound_exact"]), 0)

    def test_source_density_crossing_and_tail_parameters_have_independent_integral_checks(self):
        model = tail_model_float()
        self.assertTrue(0 < model["lower_mean_x"] < model["crossing"] < float(CROSSING_UPPER) < ALPHA < model["upper_mean_x"] < 1)
        self.assertGreater(model["lower_mean_y"], float(LOWER_TAIL_MEAN_Y_LOWER))
        self.assertGreater(model["upper_mean_y"], float(UPPER_TAIL_MEAN_Y_LOWER))
        self.assertGreater(model["upper_mean_x"], float(UPPER_TAIL_MEAN_X_LOWER))
        s = model["crossing"]
        self.assertAlmostEqual(density_factor_float(s, model["eta"]) * math.sqrt(1 - s*s), math.sqrt(ALPHA**2 - s*s), places=14)
        theta, weights = source_quadrature()
        self.assertAlmostEqual(weights.sum(), 1, places=14)
        self.assertAlmostEqual(float(weights @ (ALPHA * np.cos(theta))), 2 / math.pi, places=14)

    def test_conditional_probabilities_are_positive_and_preserve_both_coordinates(self):
        angles = np.r_[np.linspace(-4, 9, 20001), 0, math.pi / 2, math.pi, 3 * math.pi / 2]
        x, y = ALPHA * np.cos(angles), ALPHA * np.sin(angles)
        row = scalar_conditional_float(np.abs(x), np.abs(y))
        for name in ("stay", "lower_tail", "upper_tail", "vertical_bias"):
            self.assertGreaterEqual(row[name].min(), 0)
            self.assertLessEqual(row[name].max(), 1 + 1e-15)
        self.assertGreater(row["height_margin"].min(), .004)
        np.testing.assert_allclose(row["stay"] + row["lower_tail"] + row["upper_tail"], 1, atol=3e-16)
        np.testing.assert_allclose(continuous_conditional_mean(angles), np.column_stack((x, y)), atol=5e-16)

    def test_continuous_target_marginal_matches_independent_high_harmonic_integrals(self):
        model = tail_model_float()
        theta, weights = source_quadrature()
        s = ALPHA * np.cos(theta)
        row = scalar_conditional_float(s, ALPHA * np.sin(theta))
        for frequency in list(range(13)) + [24, 32]:
            phi = lambda t: np.cos(frequency * np.arccos(t))
            lower = model["lower_weights"] @ phi(model["lower_nodes"])
            upper = model["upper_weights"] @ phi(model["upper_nodes"])
            actual = weights @ (row["stay"] * phi(s) + row["lower_tail"] * lower + row["upper_tail"] * upper)
            expected = 1 if frequency == 0 else 2 * math.sin(frequency * math.pi / 2) / (frequency * math.pi)
            self.assertAlmostEqual(float(actual), expected, places=12)

    def test_tail_mass_balance_and_nonzero_continuum_component(self):
        model = tail_model_float()
        theta, weights = source_quadrature()
        row = scalar_conditional_float(ALPHA * np.cos(theta), ALPHA * np.sin(theta))
        for name, expected in (("lower_tail", model["lower_mass"]), ("upper_tail", model["upper_mass"]), ("stay", model["stationary_mass"])):
            self.assertAlmostEqual(float(weights @ row[name]), expected, places=13)
        self.assertGreater(model["stationary_mass"], .9)

    def test_full_selective_density_at_the_endpoint_for_both_marks_and_rotations(self):
        from compensated_angle_kernel import warped_angle
        from noncontextual_update_kernel import summary_density
        from retention_tradeoff import fresh_branch_matrix
        eta = tail_model_float()["eta"]
        outgoing = np.linspace(-3, 9, 173)
        for setting in (0, .29, -1.7):
            rotation = np.array([[math.cos(setting), -math.sin(setting)], [math.sin(setting), math.cos(setting)]])
            for mark in (0, 1):
                relative = outgoing - setting
                mean = continuous_conditional_mean(warped_angle(relative, mark, eta)) @ rotation.T
                response = (1 + (2 * mark - 1) * eta * np.cos(relative)) / 2
                for summary in ([1, ALPHA, 0], [1, ALPHA / math.sqrt(2), ALPHA / math.sqrt(2)], [.4, .1, -.2]):
                    actual = response / (2 * math.pi) * (summary[0] + mean @ summary[1:] / ALPHA)
                    expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, outgoing)
                    np.testing.assert_allclose(actual, expected, atol=3e-15, rtol=0)


def results_report(checks):
    model = tail_model_float()
    diagnostics = {key: value for key, value in model.items() if not isinstance(value, np.ndarray)}
    s = np.linspace(model["crossing"], ALPHA, 100001)
    margin = scalar_conditional_float(s)["height_margin"]
    diagnostics.update(minimum_height_margin_sampled=float(margin.min()), minimum_at_source_x_sampled=float(s[margin.argmin()]))
    return {"certificate": verify_endpoint(), "floating_diagnostics": diagnostics,
            "endpoint_feasibility_proven": True, "exact_family_threshold": "eta_* = eta_B, and the endpoint is attained",
            "finite_bridge_at_endpoint": False, "logarithmic_finite_label_sufficiency_proven": False,
            "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(EndpointTailBridgeTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = results_report(checks)
    if not report["certificate"]["verified"]:
        raise SystemExit(1)
    if args.write_results:
        Path(__file__).with_name("endpoint_tail_bridge_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
