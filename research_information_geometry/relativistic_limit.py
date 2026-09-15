"""Test a conditional-shift quantum walk against its controlled Dirac limit."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np


IDENTITY = np.eye(2, dtype=complex)
SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA_Z = np.diag([1.0, -1.0]).astype(complex)
SPEED = 1.0
MASS_FREQUENCY = 0.7
STEP_SIZES = (0.2, 0.1, 0.05, 0.025)
FINAL_TIME = 2.0


def coin_rotation(angle):
    return math.cos(angle) * IDENTITY - 1j * math.sin(angle) * SIGMA_X


def conditional_shift_mode(wave_number, spacing):
    return np.diag([np.exp(-1j * wave_number * spacing), np.exp(1j * wave_number * spacing)])


def walk_mode(wave_number, time_step, speed=SPEED, mass_frequency=MASS_FREQUENCY, symmetric=True):
    if time_step <= 0 or speed <= 0:
        raise ValueError("Time step and propagation speed must be positive.")
    shift = conditional_shift_mode(wave_number, speed * time_step)
    if symmetric:
        half_coin = coin_rotation(mass_frequency * time_step / 2)
        return half_coin @ shift @ half_coin
    return shift @ coin_rotation(mass_frequency * time_step)


def dirac_hamiltonian(wave_number, speed=SPEED, mass_frequency=MASS_FREQUENCY):
    return speed * wave_number * SIGMA_Z + mass_frequency * SIGMA_X


def dirac_propagator(wave_number, duration, speed=SPEED, mass_frequency=MASS_FREQUENCY):
    energy = math.hypot(speed * wave_number, mass_frequency)
    if energy == 0:
        return IDENTITY.copy()
    return math.cos(energy * duration) * IDENTITY - 1j * math.sin(energy * duration) / energy * dirac_hamiltonian(wave_number, speed, mass_frequency)


def quasienergy(wave_number, time_step, speed=SPEED, mass_frequency=MASS_FREQUENCY):
    phase_cosine = math.cos(speed * wave_number * time_step) * math.cos(mass_frequency * time_step)
    return math.acos(max(-1.0, min(1.0, phase_cosine))) / time_step


def group_velocity(wave_number, time_step, speed=SPEED, mass_frequency=MASS_FREQUENCY):
    momentum_phase = speed * wave_number * time_step
    coin_angle = mass_frequency * time_step
    phase_cosine = math.cos(momentum_phase) * math.cos(coin_angle)
    denominator = math.sqrt(max(0, 1 - phase_cosine**2))
    if denominator < 1e-14:
        raise ValueError("The selected eigenphase branch is not differentiable at this point.")
    return speed * math.cos(coin_angle) * math.sin(momentum_phase) / denominator


def real_space_step(state, time_step, mass_frequency=MASS_FREQUENCY, symmetric=True):
    if state.ndim != 2 or state.shape[1] != 2 or time_step <= 0:
        raise ValueError("Use a two-component position state and a positive time step.")
    angle = mass_frequency * time_step
    mixed = state @ coin_rotation(angle / 2 if symmetric else angle).T
    shifted = np.empty_like(mixed)
    shifted[:, 0] = np.roll(mixed[:, 0], 1)
    shifted[:, 1] = np.roll(mixed[:, 1], -1)
    return shifted @ coin_rotation(angle / 2).T if symmetric else shifted


def evolve_walk(state, time_step, steps, mass_frequency=MASS_FREQUENCY, symmetric=True):
    evolved = np.array(state, dtype=complex, copy=True)
    for _ in range(steps):
        evolved = real_space_step(evolved, time_step, mass_frequency, symmetric)
    return evolved


def mode_convergence():
    wave_numbers = (-2.0, -0.8, 0.0, 0.8, 2.0)
    rows = []
    for time_step in STEP_SIZES:
        steps = round(FINAL_TIME / time_step)
        mode_errors = {}
        for name, symmetric in (("ordered", False), ("symmetric", True)):
            errors = [
                np.linalg.norm(np.linalg.matrix_power(walk_mode(wave_number, time_step, symmetric=symmetric), steps)
                               - dirac_propagator(wave_number, FINAL_TIME), ord=2)
                for wave_number in wave_numbers
            ]
            mode_errors[name] = float(max(errors))
        energy_errors = [abs(quasienergy(wave_number, time_step) - math.hypot(SPEED * wave_number, MASS_FREQUENCY))
                         for wave_number in wave_numbers]
        rows.append({
            "time_step": time_step,
            "spacing": SPEED * time_step,
            "steps": steps,
            "ordered_operator_error": mode_errors["ordered"],
            "symmetric_operator_error": mode_errors["symmetric"],
            "max_positive_quasienergy_error": max(energy_errors),
        })
    return {
        "wave_numbers": wave_numbers,
        "final_time": FINAL_TIME,
        "rows": rows,
        "ordered_fitted_log_error_log_step_slope": float(np.polyfit(np.log(STEP_SIZES), np.log([row["ordered_operator_error"] for row in rows]), 1)[0]),
        "symmetric_fitted_log_error_log_step_slope": float(np.polyfit(np.log(STEP_SIZES), np.log([row["symmetric_operator_error"] for row in rows]), 1)[0]),
    }


def support_check():
    site_count, steps, time_step = 101, 12, 0.1
    origin = site_count // 2
    initial = np.zeros((site_count, 2), dtype=complex)
    initial[origin] = np.array([1, 1j]) / math.sqrt(2)
    evolved = evolve_walk(initial, time_step, steps)
    probabilities = np.sum(abs(evolved)**2, axis=1)
    distances = abs(np.arange(site_count) - origin)
    return {
        "sites": site_count,
        "steps": steps,
        "time_step": time_step,
        "spacing": SPEED * time_step,
        "total_probability": float(probabilities.sum()),
        "probability_outside_step_cone": float(probabilities[distances > steps].sum()),
        "max_site_distance_with_nonzero_probability": int(distances[probabilities > 0].max()),
        "caveat": "Finite-step support is imposed by the nearest-neighbor shift; the simulation ends before periodic wraparound.",
    }


def gaussian_packet(time_step, interval_length=40.0, center=-4.0, width=1.5, wave_number=0.8):
    spacing = SPEED * time_step
    site_count = round(interval_length / spacing)
    if not math.isclose(site_count * spacing, interval_length):
        raise ValueError("The interval must contain an integer number of lattice cells.")
    positions = (np.arange(site_count) - site_count // 2) * spacing
    _, spinors = np.linalg.eigh(dirac_hamiltonian(wave_number))
    spinor = spinors[:, 1]
    envelope = np.exp(-(positions - center)**2 / (4 * width**2) + 1j * wave_number * positions)
    packet = envelope[:, None] * spinor[None, :]
    packet /= np.linalg.norm(packet)
    return positions, packet


def continuum_grid_solution(initial, spacing, duration, speed=SPEED, mass_frequency=MASS_FREQUENCY):
    frequencies = 2 * math.pi * np.fft.fftfreq(initial.shape[0], d=spacing)
    transformed = np.fft.fft(initial, axis=0)
    energies = np.sqrt((speed * frequencies)**2 + mass_frequency**2)
    applied_hamiltonian = speed * frequencies[:, None] * (transformed @ SIGMA_Z.T) + mass_frequency * (transformed @ SIGMA_X.T)
    evolved_modes = np.cos(energies * duration)[:, None] * transformed
    evolved_modes -= 1j * (duration * np.sinc(energies * duration / math.pi))[:, None] * applied_hamiltonian
    return np.fft.ifft(evolved_modes, axis=0)


def packet_convergence():
    rows = []
    for time_step in STEP_SIZES:
        positions, initial = gaussian_packet(time_step)
        steps = round(FINAL_TIME / time_step)
        reference = continuum_grid_solution(initial, SPEED * time_step, FINAL_TIME)
        reference_density = np.sum(abs(reference)**2, axis=1)
        row = {"time_step": time_step, "sites": len(positions), "steps": steps}
        for name, symmetric in (("ordered", False), ("symmetric", True)):
            evolved = evolve_walk(initial, time_step, steps, symmetric=symmetric)
            density = np.sum(abs(evolved)**2, axis=1)
            row[name] = {
                "state_l2_error": float(np.linalg.norm(evolved - reference)),
                "density_total_variation": float(0.5 * np.sum(abs(density - reference_density))),
                "norm_squared": float(np.sum(density)),
                "position_mean_error": float(abs(np.sum(positions * (density - reference_density)))),
            }
        rows.append(row)
    return {
        "interval_length": 40.0,
        "packet_center": -4.0,
        "packet_width": 1.5,
        "central_wave_number": 0.8,
        "spinor": "Positive-energy eigenvector of the continuum Hamiltonian at the central wave number.",
        "reference": "Independent Fourier propagation by the continuum Dirac generator on the same periodic interval and sampled initial state.",
        "norm_convention": "Discrete sum of two-component absolute squares is one, corresponding to sqrt(spacing)-weighted continuum samples.",
        "rows": rows,
        "ordered_fitted_log_error_log_step_slope": float(np.polyfit(np.log(STEP_SIZES), np.log([row["ordered"]["state_l2_error"] for row in rows]), 1)[0]),
        "symmetric_fitted_log_error_log_step_slope": float(np.polyfit(np.log(STEP_SIZES), np.log([row["symmetric"]["state_l2_error"] for row in rows]), 1)[0]),
        "caveat": "A smooth low-momentum packet and finite duration, not a uniform-in-momentum or infinite-time convergence claim. Periodic boundary conditions and the spatial coordinate are inputs.",
    }


def zone_edge_control():
    return {
        "rows": [
            {
                "time_step": time_step,
                "wave_number_pi_over_spacing": math.pi / (SPEED * time_step),
                "lattice_positive_quasienergy": quasienergy(math.pi / (SPEED * time_step), time_step),
                "continuum_positive_energy": math.hypot(math.pi / time_step, MASS_FREQUENCY),
            }
            for time_step in STEP_SIZES
        ],
        "meaning": "Keeping lattice phase k*spacing=pi does not keep physical momentum finite; the fixed-low-momentum continuum approximation is not uniform over the Brillouin zone.",
    }


class RelativisticLimitTests(unittest.TestCase):
    def test_packet_convergence_against_continuum_reference(self):
        results = packet_convergence()
        self.assertAlmostEqual(results["ordered_fitted_log_error_log_step_slope"], 1, delta=0.08)
        self.assertAlmostEqual(results["symmetric_fitted_log_error_log_step_slope"], 2, delta=0.08)
        for row in results["rows"]:
            for name in ("ordered", "symmetric"):
                self.assertAlmostEqual(row[name]["norm_squared"], 1, places=11)
                self.assertLessEqual(row[name]["density_total_variation"], row[name]["state_l2_error"] + 1e-12)

    def test_fourier_reference_agrees_with_independent_mode_solution(self):
        site_count, spacing, mode_index = 60, 0.1, 2
        wave_number = 2 * math.pi * mode_index / (site_count * spacing)
        spinor = np.array([0.3j, 1], dtype=complex)
        spinor /= np.linalg.norm(spinor)
        profile = np.exp(1j * wave_number * spacing * np.arange(site_count)) / math.sqrt(site_count)
        initial = profile[:, None] * spinor[None, :]
        actual = continuum_grid_solution(initial, spacing, 0.73)
        expected = profile[:, None] * (dirac_propagator(wave_number, 0.73) @ spinor)[None, :]
        np.testing.assert_allclose(actual, expected, atol=1e-12)

    def test_zone_edge_is_not_a_uniform_continuum_limit(self):
        for time_step in STEP_SIZES:
            wave_number = math.pi / (SPEED * time_step)
            expected = math.pi / time_step - MASS_FREQUENCY
            self.assertAlmostEqual(quasienergy(wave_number, time_step), expected, places=10)
            difference = math.hypot(SPEED * wave_number, MASS_FREQUENCY) - quasienergy(wave_number, time_step)
            self.assertGreater(difference, 0.69)

    def test_zone_shift_changes_walk_sign(self):
        for symmetric in (False, True):
            shifted = walk_mode(0.8 + math.pi / (SPEED * 0.1), 0.1, symmetric=symmetric)
            np.testing.assert_allclose(shifted, -walk_mode(0.8, 0.1, symmetric=symmetric), atol=1e-12)

    def test_walk_modes_are_unitary_and_reflection_symmetric(self):
        for time_step in STEP_SIZES:
            for wave_number in (-3.2, -0.8, 0, 1.3):
                for symmetric in (False, True):
                    update = walk_mode(wave_number, time_step, symmetric=symmetric)
                    np.testing.assert_allclose(update.conj().T @ update, IDENTITY, atol=1e-12)
                    self.assertAlmostEqual(abs(np.linalg.det(update) - 1), 0, places=12)
                    reflected = SIGMA_X @ walk_mode(-wave_number, time_step, symmetric=symmetric) @ SIGMA_X
                    np.testing.assert_allclose(update, reflected, atol=1e-12)

    def test_eigenvalues_match_independent_trace_dispersion(self):
        for time_step in STEP_SIZES:
            for wave_number in (-2, 0, 0.8, 2):
                phase = quasienergy(wave_number, time_step) * time_step
                for symmetric in (False, True):
                    eigenvalues = np.linalg.eigvals(walk_mode(wave_number, time_step, symmetric=symmetric))
                    np.testing.assert_allclose(np.sort(np.angle(eigenvalues)), [-phase, phase], atol=1e-12)

    def test_massless_walk_is_exact_conditional_translation(self):
        for wave_number in (-2.1, 0, 1.4):
            update = walk_mode(wave_number, 0.1, mass_frequency=0)
            exact = dirac_propagator(wave_number, 0.1, mass_frequency=0)
            np.testing.assert_allclose(update, exact, atol=1e-12)

    def test_dirac_generator_and_squared_dispersion(self):
        time_step, wave_number = 1e-5, 0.8
        expected = dirac_hamiltonian(wave_number)
        update = walk_mode(wave_number, time_step)
        estimate = 1j * (update - update.conj().T) / (2 * time_step)
        np.testing.assert_allclose(estimate, expected, atol=1e-9)
        np.testing.assert_allclose(expected @ expected, ((SPEED * wave_number)**2 + MASS_FREQUENCY**2) * IDENTITY, atol=1e-12)

    def test_global_mode_convergence_orders(self):
        results = mode_convergence()
        ordered_errors = [row["ordered_operator_error"] for row in results["rows"]]
        symmetric_errors = [row["symmetric_operator_error"] for row in results["rows"]]
        self.assertTrue(all(later < earlier for earlier, later in zip(ordered_errors, ordered_errors[1:])))
        self.assertTrue(all(later < earlier for earlier, later in zip(symmetric_errors, symmetric_errors[1:])))
        self.assertAlmostEqual(results["ordered_fitted_log_error_log_step_slope"], 1, delta=0.08)
        self.assertAlmostEqual(results["symmetric_fitted_log_error_log_step_slope"], 2, delta=0.08)

    def test_real_space_step_agrees_with_mode_update(self):
        site_count, time_step, mode_index = 48, 0.1, 3
        momentum_phase = 2 * math.pi * mode_index / site_count
        wave_number = momentum_phase / (SPEED * time_step)
        spinor = np.array([1, 0.3j])
        spinor /= np.linalg.norm(spinor)
        profile = np.exp(1j * momentum_phase * np.arange(site_count)) / math.sqrt(site_count)
        initial = profile[:, None] * spinor[None, :]
        for symmetric in (False, True):
            actual = real_space_step(initial, time_step, symmetric=symmetric)
            expected = profile[:, None] * (walk_mode(wave_number, time_step, symmetric=symmetric) @ spinor)[None, :]
            np.testing.assert_allclose(actual, expected, atol=1e-12)

    def test_compact_support_stays_in_imposed_step_cone(self):
        result = support_check()
        self.assertAlmostEqual(result["total_probability"], 1, places=12)
        self.assertEqual(result["probability_outside_step_cone"], 0)
        self.assertLessEqual(result["max_site_distance_with_nonzero_probability"], result["steps"])

    def test_group_velocity_is_bounded_by_lattice_speed(self):
        time_step = 0.1
        for momentum_phase in np.linspace(-math.pi, math.pi, 101):
            velocity = group_velocity(momentum_phase / (SPEED * time_step), time_step)
            self.assertLessEqual(abs(velocity), SPEED + 1e-12)

    def test_dispersion_correction_has_predicted_second_order_coefficient(self):
        wave_number, time_step = 0.8, 0.02
        expected = -(SPEED * wave_number)**2 * MASS_FREQUENCY**2 / 3
        actual = (quasienergy(wave_number, time_step)**2 - (SPEED * wave_number)**2 - MASS_FREQUENCY**2) / time_step**2
        self.assertAlmostEqual(actual, expected, delta=2e-5)

    def test_fixed_coin_does_not_give_fixed_finite_mass_in_this_scaling(self):
        angle = 0.3
        gaps = [quasienergy(0, time_step, mass_frequency=angle / time_step) for time_step in STEP_SIZES]
        np.testing.assert_allclose(gaps, [angle / time_step for time_step in STEP_SIZES], atol=1e-12)
        self.assertAlmostEqual(gaps[-1] / gaps[0], 8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RelativisticLimitTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {
        "scope": "A standard two-channel unitary quantum-walk construction and its 1+1 dimensional Dirac limit; not emergent gravity or a novelty claim.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "units": "hbar=1, speed=spacing/time_step=1, mass_frequency=m*c^2/hbar=0.7.",
        "fixed_inputs": ["One-dimensional lattice and external clock", "Two internal channels", "Opposite conditional shifts", "Reflection-compatible onsite mixing", "Coin angle scaled as mass_frequency*time_step"],
        "criteria_fixed_before_run": ["Exact unitarity and one-cell support per update", "Mode eigenphases agree with cos(phase)=cos(c*k*time_step)*cos(mass_frequency*time_step)", "At fixed low momentum and total time, the walk approaches the independently computed Dirac propagator", "Keeping the coin angle fixed is not a fixed-finite-mass continuum limit under spacing=c*time_step"],
        "mode_convergence": mode_convergence(),
        "packet_convergence": packet_convergence(),
        "finite_support": support_check(),
        "fixed_coin_control": {"angle": 0.3, "time_steps": STEP_SIZES, "positive_gaps": [0.3 / time_step for time_step in STEP_SIZES]},
        "zone_edge_control": zone_edge_control(),
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    encoded = json.dumps(results, indent=2)
    if arguments.write_results:
        output = Path(__file__).resolve().parent / "relativistic_results.json"
        output.write_text(encoded + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(encoded)


if __name__ == "__main__":
    main()