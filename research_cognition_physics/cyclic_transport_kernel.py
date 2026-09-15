"""Exact circle transport with symmetric small shifts and a resonant reset.

Twelve shifts define a continuous-angle Green function, not a finite hidden
state approximation. A rational certificate proves all contrasts <= 0.11.
"""

import argparse
import json
import math
import platform
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from atomic_reset_kernel import atom_target, minimum_reset_weight
from compensated_angle_kernel import preconvolution_density, warped_angle
from noncontextual_update_kernel import summary_density
from outcome_updates import WINDOW_ATTENUATION as ALPHA
from retention_tradeoff import fresh_branch_matrix
from stopped_readout_witness import alpha_interval


ORDER = 12
STEP = math.pi / ORDER
COSINE = math.cos(STEP)
SHIFTS = (0.0, STEP, -STEP, math.pi)


def cyclic_green():
    j = np.arange(ORDER)
    return (ORDER**2 - 1) / (12 * ORDER) - j * (ORDER - j) / (2 * ORDER)


def resonance(angles, contrast):
    """Return H-1 and its positive Fourier amplitude bound V0.

    Directly averaging R and subtracting one loses this tiny but essential
    correction. The closed form retains it even when it is about 1e-29.
    """
    if not 0 <= contrast <= 0.11:
        raise ValueError("This implementation is certified for 0 <= contrast <= 0.11.")
    q = math.sqrt(1 - contrast**2)
    power = (contrast / (1 + q)) ** (2 * ORDER)
    z = power * np.exp(2j * ORDER * np.asarray(angles))
    deviation = 2 * np.real(z / (1 - z) + 2 * ORDER * q * z / (1 - z)**2)
    ceiling = 2 * (power / (1 - power) + 2 * ORDER * q * power / (1 - power)**2)
    return deviation, ceiling


def shift_mass(angles, contrast):
    angles = np.asarray(angles)
    q = math.sqrt(1 - contrast**2)
    maximum_r = (1 + contrast**2) / q
    offset = (maximum_r - 1) / (4 * math.sin(STEP)**2)
    particular = np.sum(cyclic_green() * (preconvolution_density(angles[..., None] + np.arange(ORDER) * STEP, contrast) - 1), axis=-1)
    return offset + particular, offset


def reverse_masses(angles, contrast):
    """Unnormalized masses A,U,C,V relative to uniform target angle.

    Divide by R to get conditional probabilities for x=psi, psi+-STEP,
    psi+pi, and a uniform independent input angle, respectively.
    """
    deviation, reset_average = resonance(angles, contrast)
    r = preconvolution_density(np.asarray(angles), contrast)
    u, _ = shift_mass(angles, contrast)
    v = reset_average + deviation
    a = ((1 + ALPHA) * r - v) / 2 - (1 + COSINE) * u
    c = ((1 - ALPHA) * r - v) / 2 - (1 - COSINE) * u
    return a, u, c, v


def floating_margins(contrast):
    q = math.sqrt(1 - contrast**2)
    _, v0 = resonance(0, contrast)
    _, u0 = shift_mass(0, contrast)
    return {"twice_antipodal_lower_bound": (1 - ALPHA) * q**3 - 2 * v0 - 4 * (1 - COSINE) * u0,
            "twice_stay_lower_bound": (1 + ALPHA) * q**3 - 2 * v0 - 4 * (1 + COSINE) * u0,
            "shift_offset": float(u0), "reset_total_mass": float(v0)}


def rational_certificate():
    eta = Fraction(11, 100)
    alpha_low, alpha_high = alpha_interval()
    q_low = 1 - eta**2 / 2 - eta**4 / 2
    # sin^2(pi/12)=(2-sqrt(3))/4; cos(pi/12)=(sqrt(6)+sqrt(2))/4.
    sine_square_low = Fraction(6698, 100000)
    cosine_low, cosine_high = Fraction(9659, 10000), Fraction(966, 1000)
    u_high = ((1 + eta**2) / q_low - 1) / (4 * sine_square_low)
    power_high = (eta / (1 + q_low)) ** (2 * ORDER)
    v_high = 2 * (power_high / (1 - power_high) + 2 * ORDER * power_high / (1 - power_high)**2)
    return {"twice_antipodal_lower_bound": (1 - alpha_high) * q_low**3 - 2 * v_high - 4 * (1 - cosine_low) * u_high,
            "twice_stay_lower_bound": (1 + alpha_low) * q_low**3 - 2 * v_high - 4 * (1 + cosine_high) * u_high}


def forward_atoms(incoming, outcome, contrast, setting=0.0):
    """Four (weight, output angle) pairs for one observed outcome."""
    sign = 2 * outcome - 1
    q = math.sqrt(1 - contrast**2)
    result = []
    for index, shift in enumerate(SHIFTS):
        psi = np.asarray(incoming) - setting - shift
        a, u, c, _ = reverse_masses(psi, contrast)
        mass = (a, u, u, c)[index]
        density = q**3 / (2 * (1 - sign * contrast * np.cos(psi))**2)
        weight = density * mass / preconvolution_density(psi, contrast)
        target = atom_target(psi, outcome, contrast) + setting
        result.append((weight, target))
    return result


def forward_reset_density(outgoing, outcome, contrast, setting=0.0):
    relative = np.asarray(outgoing) - setting
    psi = warped_angle(relative, outcome, contrast)
    _, _, _, v = reverse_masses(psi, contrast)
    f = (1 + (2 * outcome - 1) * contrast * np.cos(relative)) / 2
    return f * v / (2 * math.pi * preconvolution_density(psi, contrast))


def pushed_output_density(summary, outgoing, outcome, contrast, setting=0.0):
    """Change variables separately in each forward atom, then add reset."""
    outgoing = np.asarray(outgoing)
    relative = outgoing - setting
    psi = warped_angle(relative, outcome, contrast)
    jacobian = math.sqrt(1 - contrast**2) / (1 + (2 * outcome - 1) * contrast * np.cos(relative))
    result = summary[0] * forward_reset_density(outgoing, outcome, contrast, setting)
    for index, shift in enumerate(SHIFTS):
        incoming = psi + setting + shift
        weight, _ = forward_atoms(incoming, outcome, contrast, setting)[index]
        result = result + weight * summary_density(summary, incoming) * jacobian
    return result


def weak_output_moment(summary, outcome, contrast, setting, frequency, sine=False, nodes=2048):
    grid = np.arange(nodes) * 2 * math.pi / nodes
    function = np.sin if sine else np.cos
    response = sum(weight * function(frequency * target) for weight, target in forward_atoms(grid, outcome, contrast, setting))
    atoms = np.mean(summary_density(summary, grid) * response) * 2 * math.pi
    reset = summary[0] * np.mean(forward_reset_density(grid, outcome, contrast, setting) * function(frequency * grid)) * 2 * math.pi
    return float(atoms + reset)


class CyclicTransportKernelTests(unittest.TestCase):
    def test_cyclic_green_solves_the_discrete_poisson_equation(self):
        green = cyclic_green()
        expected = np.full(ORDER, -1 / ORDER)
        expected[0] += 1
        np.testing.assert_allclose(2 * green - np.roll(green, 1) - np.roll(green, -1), expected, atol=1e-14)
        self.assertAlmostEqual(float(sum(green)), 0, places=14)

    def test_resonance_closed_form_keeps_tiny_nonzero_mass(self):
        grid = np.arange(480) * 2 * math.pi / 480
        deviation, ceiling = resonance(grid, 0.11)
        self.assertGreater(ceiling, 3e-29)
        self.assertLess(ceiling, 4e-29)
        self.assertLess(abs(float(np.mean(deviation))), 1e-42)
        self.assertLessEqual(float(np.max(np.abs(deviation))), ceiling * (1 + 1e-14))
        # Numerical quadrature confirms the reset normalization relative to its tiny scale.
        mass = sum(np.mean(forward_reset_density(grid, mark, 0.11)) * 2 * math.pi for mark in (0, 1))
        self.assertAlmostEqual(float(mass / ceiling), 1, places=12)

    def test_green_solution_removes_all_nonresonant_harmonics(self):
        angles = np.linspace(-1, 5, 101)
        for eta in (0, 0.03, 0.11):
            u, _ = shift_mass(angles, eta)
            left, _ = shift_mass(angles - STEP, eta)
            right, _ = shift_mass(angles + STEP, eta)
            deviation, _ = resonance(angles, eta)
            np.testing.assert_allclose(2 * u - left - right, preconvolution_density(angles, eta) - 1 - deviation, atol=1e-14)

    def test_rational_bounds_used_in_certificate(self):
        self.assertGreater(Fraction(173208, 100000)**2, 3)
        self.assertLess(Fraction(244948, 100000)**2, 6)
        self.assertGreater(Fraction(244950, 100000)**2, 6)
        self.assertLess(Fraction(141421, 100000)**2, 2)
        self.assertGreater(Fraction(141422, 100000)**2, 2)
        self.assertGreater((Fraction(244948, 100000) + Fraction(141421, 100000)) / 4, Fraction(9659, 10000))
        self.assertLess((Fraction(244950, 100000) + Fraction(141422, 100000)) / 4, Fraction(966, 1000))
        for margin in rational_certificate().values():
            self.assertGreater(margin, 0)

    def test_all_reverse_masses_are_nonnegative_and_normalized(self):
        grid = np.linspace(-math.pi, math.pi, 1001)
        for eta in (0, 0.03, 0.08, 0.11):
            a, u, c, v = reverse_masses(grid, eta)
            self.assertGreaterEqual(float(np.min(a)), 0)
            self.assertGreaterEqual(float(np.min(u)), -1e-14)
            self.assertGreaterEqual(float(np.min(c)), 0)
            self.assertGreaterEqual(float(np.min(v)), -1e-42)
            np.testing.assert_allclose(a + 2 * u + c + v, preconvolution_density(grid, eta), atol=1e-14)

    def test_conditional_vector_mean_is_exactly_alpha_times_direction(self):
        grid = np.linspace(-3, 3, 61)
        a, u, c, _ = reverse_masses(grid, 0.11)
        actual = a * np.exp(1j * grid) + u * (np.exp(1j * (grid + STEP)) + np.exp(1j * (grid - STEP))) + c * np.exp(1j * (grid + math.pi))
        expected = ALPHA * preconvolution_density(grid, 0.11) * np.exp(1j * grid)
        np.testing.assert_allclose(actual, expected, atol=1e-14)

    def test_reverse_coupling_has_exact_uniform_input_marginal(self):
        grid = np.linspace(-2, 8, 201)
        for eta in (0, 0.07, 0.11):
            a, _, _, _ = reverse_masses(grid, eta)
            _, left, _, _ = reverse_masses(grid - STEP, eta)
            _, right, _, _ = reverse_masses(grid + STEP, eta)
            _, _, opposite, _ = reverse_masses(grid - math.pi, eta)
            _, reset = resonance(0, eta)
            np.testing.assert_allclose(a + left + right + opposite + reset, 1, atol=1e-14)

    def test_forward_kernel_normalizes_at_every_input(self):
        grid = np.linspace(-3, 6, 71)
        for eta in (0, 0.04, 0.11):
            for setting in (0, 0.7):
                _, reset = resonance(0, eta)
                total = reset + sum(weight for mark in (0, 1) for weight, _ in forward_atoms(grid, mark, eta, setting))
                np.testing.assert_allclose(total, 1, atol=1e-14)

    def test_full_selective_output_density_matches_original_instrument(self):
        outgoing = np.linspace(-2, 8, 37)
        for summary in ([1, 0, 0], [1, ALPHA, 0], [1, 0, ALPHA], [0.4, 0.1, -0.2]):
            for eta in (0, 0.04, 0.11):
                for setting in (0, 0.7):
                    for mark in (0, 1):
                        expected = summary_density(fresh_branch_matrix(setting, mark, eta) @ summary, outgoing)
                        np.testing.assert_allclose(pushed_output_density(summary, outgoing, mark, eta, setting), expected, atol=1e-14)

    def test_direct_forward_measure_integrals_include_higher_moments(self):
        summary = np.array([1, 0.3, -0.4])
        for setting in (0, 0.7):
            for mark in (0, 1):
                output = fresh_branch_matrix(setting, mark, 0.11) @ summary
                self.assertAlmostEqual(weak_output_moment(summary, mark, 0.11, setting, 0), output[0], places=13)
                self.assertAlmostEqual(weak_output_moment(summary, mark, 0.11, setting, 1), output[1] / (2 * ALPHA), places=13)
                self.assertAlmostEqual(weak_output_moment(summary, mark, 0.11, setting, 1, True), output[2] / (2 * ALPHA), places=13)
                for frequency in (2, 3, 6):
                    self.assertAlmostEqual(weak_output_moment(summary, mark, 0.11, setting, frequency), 0, places=13)

    def test_actual_preparation_mixtures_keep_the_same_entire_distribution(self):
        grid = np.linspace(0, 6, 53)
        first, second = np.array([1, ALPHA, 0]), np.array([1, 0, ALPHA])
        mixed = (pushed_output_density(first, grid, 1, 0.11) + pushed_output_density(second, grid, 1, 0.11)) / 2
        np.testing.assert_allclose(mixed, pushed_output_density((first + second) / 2, grid, 1, 0.11), atol=1e-14)

    def test_zero_read_preserves_original_alpha_disturbance_and_old_formula_fails_at_011(self):
        grid = np.array([0.2, 0.7, 2.3])
        a, u, c, v = reverse_masses(grid, 0)
        np.testing.assert_allclose(a, (1 + ALPHA) / 2)
        np.testing.assert_allclose(c, (1 - ALPHA) / 2)
        np.testing.assert_allclose(u, 0)
        np.testing.assert_allclose(v, 0)
        self.assertLess(minimum_reset_weight(0.11), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CyclicTransportKernelTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    exact = rational_certificate()
    results = {"environment": {"python": platform.python_version(), "numpy": np.__version__},
               "automated_checks": {"run": checks.testsRun, "failures": len(checks.failures), "errors": len(checks.errors)},
               "certified_all_contrasts_up_to": 0.11, "cyclic_order": ORDER,
               "floating_margins_at_ceiling": floating_margins(0.11),
               "rational_lower_bounds_decimal": {key: float(value) for key, value in exact.items()},
               "rational_lower_bounds_exact": {key: str(value) for key, value in exact.items()},
               "scope": "A continuous hidden angle with four atomic shifts per observed outcome and a tiny but essential continuous reset exactly implements all full selective affine preparation densities for every eta<=0.11. The cyclic Green function has twelve terms; this is not a twelve-state hidden model. Rational global positivity certificates and analytic full-distribution identities cover all directions and finite adaptive compositions. No optimality is claimed."}
    if args.write_results:
        Path(__file__).with_name("cyclic_transport_kernel_results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    # The exact rational certificates are persisted; keep terminal output compact.
    print(json.dumps({key: value for key, value in results.items() if key != "rational_lower_bounds_exact"}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
