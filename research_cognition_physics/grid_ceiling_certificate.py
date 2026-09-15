"""Certify an exact implicit grid strength ceiling using an existing half-block candidate."""

import argparse
import hashlib
import json
import math
import unittest
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, pi_interval, sin_interval, cos_interval
from polygon_limit_analysis import factor_interval, grid_certificate
from polygon_martingale_certificate import certified_geometry, quantize_source
from half_face_certificate import (CASES, RECERTIFIED_CONTRAST, candidate_path, read_candidate,
                                   half_interval_geometry, _verify_half_geometry, repaired_half_float)
from half_polygon_reduction import expand_half_coupling


@lru_cache(maxsize=8)
def grid_root_bracket(size):
    if type(size) is not int or size < 32 or size % 4:
        raise ValueError("Use N>=32 divisible by four.")
    alpha = 4 * sin_interval(I.rational(1, 4))
    cosine = cos_interval(pi_interval() / size)
    lower, upper = Fraction(0), Fraction(1, 5)

    def equation(eta):
        return alpha * factor_interval(eta) - cosine

    if equation(lower).hi >= 0 or equation(upper).lo <= 0:
        raise ArithmeticError("The starting interval does not bracket a positive root.")
    for _ in range(64):
        middle = (lower + upper) / 2
        value = equation(middle)
        if value.hi < 0:
            lower = middle
        elif value.lo > 0:
            upper = middle
        else:
            break
    if upper - lower > Fraction(1, 10**18):
        raise ArithmeticError("Increase arithmetic precision before claiming a tight bracket.")
    return lower, upper, equation(lower), equation(upper)


@lru_cache(maxsize=2)
def ceiling_geometry(case):
    size, _, radius_float = CASES[case]
    radius = Fraction(str(radius_float))
    lower, upper, left_value, right_value = grid_root_bracket(size)
    if left_value.hi >= 0 or right_value.lo <= 0:
        raise ArithmeticError("Root endpoint signs are not strict.")
    eta = I(I.exact(lower).lo, I.exact(upper).hi)
    geometry = dict(half_interval_geometry(size, eta, radius))
    full = certified_geometry(size, eta, radius)
    # At the exact root, alpha*G(eta)/cos(pi/N)=1 identically. These substitutions
    # apply only to that root; the complete enclosing interval is NOT feasible.
    geometry["epsilon"] = I.exact(0)
    geometry["axis_radius_ratio"] = I.exact(radius) / full["beta"]
    quarter = size // 4
    geometry["sources"] = [full["sources"][i % size] for i in range(-quarter + 1, quarter)]
    return geometry


def verify_grid_ceiling(case="2048", data=None):
    data = read_candidate(case) if data is None else data
    size, _, radius = CASES[case]
    actual = tuple(data["report"][key] for key in ("size", "contrast", "source_radius"))
    if actual != CASES[case]:
        raise ValueError("Candidate parameters do not match the selected grid.")
    lower, upper, left_value, right_value = grid_root_bracket(size)
    report = _verify_half_geometry(data, ceiling_geometry(case), "eta_projection_" + str(size), str(Fraction(str(radius))))
    report.update({"root_definition": "alpha*G(eta)=cos(pi/N), unique positive root",
                   "root_bracket_exact": [str(lower), str(upper)],
                   "root_bracket_float_diagnostic": [float(lower), float(upper)],
                   "equation_at_lower_upper_bound_exact": str(Fraction(left_value.hi, SCALE)),
                   "equation_at_upper_lower_bound_exact": str(Fraction(right_value.lo, SCALE)),
                   "root_bracket_width_exact": str(upper - lower),
                   "scope": ("The exact root alone is certified feasible. It is the maximum strength for the N fixed equal output arcs centered at 2*pi*j/N, not the continuous endpoint and not a feasibility claim for the whole bracket."
                             if report["verified"] else
                             "This candidate did not certify the exact grid root. Failed positivity certification is not an infeasibility proof.")})
    return report


@lru_cache(maxsize=1)
def root_float_model():
    lower, upper, _, _ = grid_root_bracket(2048)
    eta = (lower + upper) / 2
    geometry, half = repaired_half_float("2048", eta)
    geometry = dict(geometry, epsilon=0, axis_radius_ratio=geometry["radius"] / geometry["beta"])
    return geometry, expand_half_coupling(geometry, half)


class GridCeilingCertificateTests(unittest.TestCase):
    def test_strict_root_bracket_and_independent_projection_exclusion(self):
        for size in (128, 512, 2048):
            lower, upper, left, right = grid_root_bracket(size)
            self.assertLess(left.hi, 0)
            self.assertGreater(right.lo, 0)
            self.assertLess(upper - lower, Fraction(1, 10**18))
            self.assertFalse(grid_certificate(size, lower)["grid_excluded"])
            self.assertTrue(grid_certificate(size, upper)["grid_excluded"])
        self.assertGreater(grid_root_bracket(2048)[0], RECERTIFIED_CONTRAST)

    def test_saved_candidate_is_strict_inside_root_half_block(self):
        report = verify_grid_ceiling()
        self.assertTrue(report["verified"])
        self.assertEqual(report["uniform_reserve_cap_interval"], [0.0, 0.0])
        self.assertGreater(report["repaired_half_entry_lower_bound"], 0)
        self.assertEqual(hashlib.sha256(candidate_path("2048").read_bytes()).hexdigest(),
                         "a5c6087afee1546fc69416fcf2c0172c7e668b3ad637536394e4a941977a34f0")

    def test_corruption_is_rejected_even_when_epsilon_is_exactly_zero(self):
        candidate = read_candidate("2048")
        candidate["exceptions"][0][2] += 1 << 65
        self.assertFalse(verify_grid_ceiling(data=candidate)["verified"])

    def test_full_core_has_nonnegative_entries_exact_zeros_and_correct_moments(self):
        geometry, matrix = root_float_model()
        size = geometry["size"]
        angles = np.arange(size) * 2 * math.pi / size
        unit = np.stack((np.cos(angles), np.sin(angles)), axis=1)
        sources = np.vstack((geometry["radius"] * unit, np.zeros(2)))
        self.assertGreaterEqual(float(matrix.min()), 0)
        self.assertEqual(np.count_nonzero(matrix[:-1][unit[:, 0, None] * unit[:, 0] < -1e-12]), 0)
        np.testing.assert_allclose(matrix.sum(axis=1), geometry["weights"], atol=1e-15)
        np.testing.assert_allclose(matrix.sum(axis=0), 1 / size, atol=1e-15)
        np.testing.assert_allclose(matrix @ (geometry["beta"] * unit), geometry["weights"][:, None] * sources, atol=1e-15)

    def test_lifted_root_instrument_matches_entire_selective_density(self):
        from compensated_angle_kernel import warped_angle
        from noncontextual_update_kernel import summary_density
        from retention_tradeoff import fresh_branch_matrix
        from outcome_updates import WINDOW_ATTENUATION as alpha
        geometry, matrix = root_float_model()
        size, eta = geometry["size"], geometry["contrast"]
        conditional = matrix / geometry["weights"][:, None]
        outgoing, setting = np.linspace(-2, 7, 41), .29
        relative = outgoing - setting
        centers = np.arange(size) * 2 * math.pi / size + setting
        for mark in (0, 1):
            left, lower, upper, center = quantize_source(warped_angle(relative, mark, eta), size, geometry["radius"])
            probabilities = (lower[:, None] * conditional[left] + upper[:, None] * conditional[(left + 1) % size]
                             + center[:, None] * conditional[-1])
            for summary in ([1, alpha / math.sqrt(2), alpha / math.sqrt(2)], [.4, .1, -.2]):
                arc_mass = (summary[0] + geometry["beta"] / alpha * (summary[1] * np.cos(centers)
                            + summary[2] * np.sin(centers))) / size
                response = (1 + (2 * mark - 1) * eta * np.cos(relative)) / 2
                observed = response * size / (2 * math.pi) * (probabilities @ arc_mass)
                expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, outgoing)
                np.testing.assert_allclose(observed, expected, atol=2e-13)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(GridCeilingCertificateTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = {"certificate": verify_grid_ceiling(),
            "candidate_sha256": hashlib.sha256(candidate_path("2048").read_bytes()).hexdigest(),
            "candidate_generated_at_contrast": "0.14473", "new_solver_run_required": False,
            "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("grid_ceiling_certificate_results.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(data, indent=2))


if __name__ == "__main__":
    main()
