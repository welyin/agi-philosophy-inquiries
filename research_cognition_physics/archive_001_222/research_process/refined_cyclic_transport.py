"""Sharper harmonic bounds for an exact twenty-term continuous circle kernel."""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction as F
from pathlib import Path

import numpy as np

from atomic_reset_kernel import atom_target
from compensated_angle_kernel import preconvolution_density as density_r, warped_angle
from noncontextual_update_kernel import summary_density
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from retention_tradeoff import fresh_branch_matrix
from stopped_readout_witness import alpha_interval


ORDER = 20
STEP = math.pi / ORDER
COSINE = math.cos(STEP)
CEILING = F(591, 5000)
SHIFTS = (0.0, STEP, -STEP, math.pi)


def sqrt_interval(value, digits=30):
    """Exact rational enclosure; integer square comparisons are the certificate."""
    scale = 10**digits
    numerator = math.isqrt(value.numerator * scale**2 // value.denominator)
    return F(numerator, scale), F(numerator + 1, scale)


def cosine_interval():
    low5, high5 = sqrt_interval(F(5))
    low18 = sqrt_interval(10 + 2 * low5)[0] / 4
    high18 = sqrt_interval(10 + 2 * high5)[1] / 4
    return sqrt_interval((1 + low18) / 2)[0], sqrt_interval((1 + high18) / 2)[1]


def parameters(contrast):
    if not 0 <= contrast <= 0.15:
        raise ValueError("Use [0,0.15] for construction diagnostics; only <=0.1182 is certified.")
    q = math.sqrt(1 - contrast**2)
    return q, (contrast / (1 + q))**2


def green_part(angles, contrast):
    indices = np.arange(ORDER)
    green = (ORDER**2 - 1) / (12 * ORDER) - indices * (ORDER - indices) / (2 * ORDER)
    return np.sum(green * (density_r(np.asarray(angles)[..., None] + indices * STEP, contrast) - 1), axis=-1)


def resonance(angles, contrast):
    q, r = parameters(contrast)
    a = r**ORDER
    z = a * np.exp(2j * ORDER * np.asarray(angles))
    deviation = 2 * np.real(z / (1 - z) + 2 * ORDER * q * z / (1 - z)**2)
    average = 2 * (a / (1 - a) + 2 * ORDER * q * a / (1 - a)**2)
    return deviation, average


def reverse_masses(angles, contrast):
    """The offset equals the exact sum of nonresonant Fourier amplitudes."""
    deviation, average = resonance(angles, contrast)
    v = average + deviation
    u = green_part(0.0, contrast) + green_part(angles, contrast)
    r = density_r(np.asarray(angles), contrast)
    a = ((1 + ALPHA) * r - v) / 2 - (1 + COSINE) * u
    c = ((1 - ALPHA) * r - v) / 2 - (1 - COSINE) * u
    return a, u, c, v


def rational_certificate(ceiling=CEILING):
    alpha_low, alpha_high = alpha_interval()
    cosine_low, cosine_high = cosine_interval()
    # A rational square comparison proves this lower bound at the ceiling.
    q_low = 1 - ceiling**2 / 2 - ceiling**4 / 8 - ceiling**6 / 8
    assert q_low > 0 and q_low**2 < 1 - ceiling**2
    r = (1 - q_low) / (1 + q_low)
    # c1(q) decreases in q. For the remaining coefficients use q <= 1.
    c1 = r * (1 + 2 * q_low)
    c2 = 5 * r**2
    tail2 = 2 * r**2 * (5 / (1 - r) + 2 * r / (1 - r)**2)
    tail3 = 2 * r**3 * (7 / (1 - r) + 2 * r / (1 - r)**2)
    power = r**ORDER
    v = 2 * (power / (1 - power) + 2 * ORDER * power / (1 - power)**2)
    low_delta, high_delta = 1 - alpha_high, 1 - alpha_low
    low_sum, high_sum = 1 + alpha_low, 1 + alpha_high
    c_margin = (low_delta - (2 / (1 + cosine_low) - 2 * low_delta) * c1
                - high_delta * tail2 - c2 / (2 * (1 + cosine_low) * cosine_low**2)
                - tail3 / (1 + cosine_low) - 2 * v)
    a_margin = (low_sum - (2 / (1 - cosine_high) - 2 * low_sum) * c1
                - high_sum * tail2 - c2 / (2 * (1 - cosine_high) * cosine_low**2)
                - tail3 / (1 - cosine_high) - 2 * v)
    return {"twice_antipodal_lower_bound": c_margin, "twice_stay_lower_bound": a_margin}


def rational_failure_at_point_12():
    eta = F(3, 25)
    alpha_low, _ = alpha_interval()
    _, cosine_high = cosine_interval()
    q_low, q_high = sqrt_interval(1 - eta**2)
    c1_low = (1 - q_high) * (1 + 2 * q_high) / (1 + q_high)
    # U_p(0) >= c1/(2 sin^2 d); dropping the positive reset gives an upper bound.
    return (1 - alpha_low) * (1 + eta**2) / q_low - 2 * c1_low / (1 + cosine_high)


def forward_atoms(incoming, outcome, contrast, setting=0.0):
    q, _ = parameters(contrast)
    sign = 2 * outcome - 1
    result = []
    for index, shift in enumerate(SHIFTS):
        psi = np.asarray(incoming) - setting - shift
        a, u, c, _ = reverse_masses(psi, contrast)
        weight = q**3 * (a, u, u, c)[index] / (2 * (1 - sign * contrast * np.cos(psi))**2 * density_r(psi, contrast))
        result.append((weight, atom_target(psi, outcome, contrast) + setting))
    return result


def reset_density(outgoing, outcome, contrast, setting=0.0):
    relative = np.asarray(outgoing) - setting
    psi = warped_angle(relative, outcome, contrast)
    v = reverse_masses(psi, contrast)[3]
    return (1 + (2 * outcome - 1) * contrast * np.cos(relative)) * v / (4 * math.pi * density_r(psi, contrast))


def pushed_density(summary, outgoing, outcome, contrast, setting=0.0):
    outgoing = np.asarray(outgoing)
    relative = outgoing - setting
    psi = warped_angle(relative, outcome, contrast)
    jacobian = math.sqrt(1 - contrast**2) / (1 + (2 * outcome - 1) * contrast * np.cos(relative))
    result = summary[0] * reset_density(outgoing, outcome, contrast, setting)
    for index, shift in enumerate(SHIFTS):
        incoming = psi + setting + shift
        weight, _ = forward_atoms(incoming, outcome, contrast, setting)[index]
        result = result + weight * summary_density(summary, incoming) * jacobian
    return result


def sampled_diagnostics():
    grid = np.linspace(0, math.pi, 4097)
    rows = []
    for eta in (0.11, 0.117, float(CEILING), 0.12):
        a, u, c, v = reverse_masses(grid, eta)
        rows.append({"contrast": eta, "minimum_stay_sampled": float(a.min()),
                     "minimum_shift_sampled": float(u.min()), "minimum_antipodal_sampled": float(c.min()),
                     "antipodal_at_zero": float(reverse_masses(0.0, eta)[2]),
                     "reset_total_mass": float(resonance(0.0, eta)[1])})
    return rows


class RefinedCyclicTransportTests(unittest.TestCase):
    def test_nested_radical_intervals_are_exact_enclosures(self):
        for value in (F(5), F(2), F(123, 100)):
            lower, upper = sqrt_interval(value)
            self.assertLessEqual(lower**2, value)
            self.assertGreater(upper**2, value)
        lower, upper = cosine_interval()
        self.assertAlmostEqual(float(lower), COSINE, places=14)
        self.assertLess(upper - lower, F(1, 10**28))

    def test_global_rational_certificate_is_positive(self):
        for margin in rational_certificate().values():
            self.assertGreater(margin, 0)

    def test_fourier_offset_and_tail_bounds(self):
        eta = float(CEILING)
        q, r = parameters(eta)
        k = np.arange(1, 100)
        k = k[k % ORDER != 0]
        amplitudes = r**k * (1 + 2 * k * q) / (2 * np.sin(k * STEP)**2)
        self.assertAlmostEqual(float(amplitudes.sum()), float(green_part(0, eta)), places=12)
        grid = np.linspace(0, math.pi, 501)
        self.assertLessEqual(float(np.max(np.abs(green_part(grid, eta)))), float(amplitudes.sum()) + 1e-12)
        old_offset = ((1 + eta**2) / q - 1) / (4 * math.sin(STEP)**2)
        self.assertLess(float(green_part(0, eta)), old_offset)

    def test_reverse_mass_mean_and_uniform_marginal(self):
        grid = np.linspace(-3, 4, 211)
        for eta in (0, 0.05, float(CEILING)):
            a, u, c, v = reverse_masses(grid, eta)
            for mass in (a, u, c):
                self.assertGreaterEqual(float(mass.min()), -1e-12)
            self.assertGreaterEqual(float(v.min()), -1e-60)
            np.testing.assert_allclose(a + 2 * u + c + v, density_r(grid, eta), atol=1e-12)
            np.testing.assert_allclose(a + 2 * COSINE * u - c, ALPHA * density_r(grid, eta), atol=1e-12)
            left = reverse_masses(grid - STEP, eta)[1]
            right = reverse_masses(grid + STEP, eta)[1]
            opposite = reverse_masses(grid - math.pi, eta)[2]
            np.testing.assert_allclose(a + left + right + opposite + resonance(0, eta)[1], 1, atol=1e-12)

    def test_forward_kernel_normalizes_for_each_input(self):
        grid = np.linspace(-3, 6, 61)
        for eta in (0, 0.07, float(CEILING)):
            for setting in (0, 0.4):
                total = resonance(0, eta)[1] + sum(w for mark in (0, 1) for w, _ in forward_atoms(grid, mark, eta, setting))
                np.testing.assert_allclose(total, 1, atol=1e-12)

    def test_entire_selective_density_matches_original_fresh_family(self):
        grid = np.linspace(-2, 8, 29)
        for summary in ([1, 0, 0], [1, ALPHA, 0], [1, 0, ALPHA], [0.4, 0.1, -0.2]):
            for eta in (0, 0.05, float(CEILING)):
                for setting in (0, 0.7):
                    for mark in (0, 1):
                        expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, grid)
                        np.testing.assert_allclose(pushed_density(summary, grid, mark, eta, setting), expected, atol=1e-12)

    def test_independent_measure_integrals_and_reset_relative_mass(self):
        grid = np.arange(2048) * 2 * math.pi / 2048
        summary, eta, setting = np.array([1, 0.3, -0.4]), float(CEILING), 0.3
        reset_mass = 0.0
        for mark in (0, 1):
            expected = fresh_branch_matrix(setting, mark, eta) @ summary
            atoms = forward_atoms(grid, mark, eta, setting)
            reset = reset_density(grid, mark, eta, setting)
            reset_mass += np.mean(reset) * 2 * math.pi
            for frequency in (0, 1, 2, 3, 6):
                for function in (np.cos, np.sin):
                    actual = (np.mean(summary_density(summary, grid) * sum(w * function(frequency * target) for w, target in atoms))
                              + summary[0] * np.mean(reset * function(frequency * grid))) * 2 * math.pi
                    goal = 0.0
                    if frequency == 0 and function is np.cos:
                        goal = expected[0]
                    if frequency == 1:
                        goal = expected[1 if function is np.cos else 2] / (2 * ALPHA)
                    self.assertAlmostEqual(float(actual), float(goal), places=11)
        self.assertAlmostEqual(float(reset_mass / resonance(0, eta)[1]), 1, places=11)

    def test_actual_preparation_mixtures_remain_identical(self):
        grid = np.linspace(0, 6, 37)
        first, second = np.array([1, ALPHA, 0]), np.array([1, 0, ALPHA])
        mixed = (pushed_density(first, grid, 1, float(CEILING)) + pushed_density(second, grid, 1, float(CEILING))) / 2
        np.testing.assert_allclose(mixed, pushed_density((first + second) / 2, grid, 1, float(CEILING)), atol=1e-12)

    def test_point_12_is_a_certified_failure_of_this_specific_kernel(self):
        self.assertLess(rational_failure_at_point_12(), 0)
        self.assertLess(float(reverse_masses(0, 0.12)[2]), 0)

    def test_zero_strength_read_keeps_the_original_disturbance(self):
        a, u, c, v = reverse_masses(np.array([0.2, 0.7]), 0)
        np.testing.assert_allclose(a, (1 + ALPHA) / 2)
        np.testing.assert_allclose(c, (1 - ALPHA) / 2)
        np.testing.assert_allclose(u, 0)
        np.testing.assert_allclose(v, 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RefinedCyclicTransportTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    cert = rational_certificate()
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "certified_all_contrasts_up_to": float(CEILING), "cyclic_order": ORDER,
               "rational_margins_exact": {key: str(value) for key, value in cert.items()},
               "rational_margins_decimal": {key: float(value) for key, value in cert.items()},
               "specific_kernel_failure_twice_antipodal_upper_bound_at_point_12": float(rational_failure_at_point_12()),
               "specific_kernel_failure_exact_upper_bound": str(rational_failure_at_point_12()),
               "sampled_diagnostics": sampled_diagnostics(),
               "scope": "All eta<=0.1182 have an explicit full preparation-noncontextual continuous-angle instrument. The sharper Fourier offset and correlated angular positivity bounds are certified rationally. Failure at eta=0.12 concerns only this formula, not all couplings. Twenty Green terms do not mean twenty hidden labels."}
    if args.write_results:
        Path(__file__).with_name("refined_cyclic_transport_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in results.items() if "exact" not in key}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
