"""Check variable-speed unitary walks and the limits of their metric interpretation."""

import argparse
import json
import math
import platform
import unittest
from pathlib import Path

import numpy as np

from relativistic_limit import SIGMA_X, SIGMA_Z


SIGMA_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
BASE_SPEED = 1.0
INTERVAL_LENGTH = 20.0
VELOCITY_SCALE = 0.6
PROFILE_AMPLITUDE = 0.25
TIME_STEPS = (0.4, 0.2, 0.1, 0.05)
FINAL_TIME = 4.0


def profile(positions):
    spatial_frequency = 2 * math.pi / INTERVAL_LENGTH
    denominator = 1 + PROFILE_AMPLITUDE * np.cos(spatial_frequency * positions)
    velocity = VELOCITY_SCALE / denominator
    first_derivative = VELOCITY_SCALE * PROFILE_AMPLITUDE * spatial_frequency * np.sin(spatial_frequency * positions) / denominator**2
    second_derivative = (
        VELOCITY_SCALE * PROFILE_AMPLITUDE * spatial_frequency**2 * np.cos(spatial_frequency * positions) / denominator**2
        + 2 * VELOCITY_SCALE * PROFILE_AMPLITUDE**2 * spatial_frequency**2 * np.sin(spatial_frequency * positions)**2 / denominator**3
    )
    return velocity, first_derivative, second_derivative


def optical_coordinate(positions):
    spatial_frequency = 2 * math.pi / INTERVAL_LENGTH
    return (positions + PROFILE_AMPLITUDE * np.sin(spatial_frequency * positions) / spatial_frequency) / VELOCITY_SCALE


def rotation_angle(velocity):
    ratio = np.asarray(velocity) / BASE_SPEED
    if np.any(ratio <= 0) or np.any(ratio >= 1):
        raise ValueError("This chart uses a smooth speed strictly between zero and the base speed.")
    return 2 * np.arccos(ratio)


def rotate(state, angles):
    cosine = np.cos(np.asarray(angles) / 2)
    sine = np.sin(np.asarray(angles) / 2)
    return np.stack((cosine * state[:, 0] - sine * state[:, 1],
                     sine * state[:, 0] + cosine * state[:, 1]), axis=1)


def shift(state, cells):
    return np.stack((np.roll(state[:, 0], cells), np.roll(state[:, 1], -cells)), axis=1)


def conjugated_shift(state, angles, cells):
    return rotate(shift(rotate(state, -angles), cells), angles)


def walk_step(state, angles, time_step, mass_frequency=0.0, symmetric=True, inverse=False):
    if state.ndim != 2 or state.shape[1] != 2 or time_step <= 0:
        raise ValueError("Use a two-component lattice state and a positive time step.")
    sign = -1 if inverse else 1
    evolved = rotate(state, sign * np.asarray(mass_frequency) * time_step)
    if symmetric:
        evolved = shift(evolved, sign)
        evolved = conjugated_shift(evolved, angles, 2 * sign)
        evolved = shift(evolved, sign)
    elif inverse:
        evolved = conjugated_shift(evolved, angles, -2)
        evolved = shift(evolved, -2)
    else:
        evolved = shift(evolved, 2)
        evolved = conjugated_shift(evolved, angles, 2)
    return rotate(evolved, sign * np.asarray(mass_frequency) * time_step)


def evolve_walk(state, angles, time_step, duration, mass_frequency=0.0, symmetric=True):
    steps = round(duration / time_step)
    if not math.isclose(steps * time_step, duration):
        raise ValueError("Duration must contain an integer number of macro steps.")
    evolved = state.copy()
    for _ in range(steps):
        evolved = walk_step(evolved, angles, time_step, mass_frequency, symmetric)
    return evolved


def spatial_grid(time_step):
    spacing = BASE_SPEED * time_step / 4
    site_count = round(INTERVAL_LENGTH / spacing)
    if not math.isclose(site_count * spacing, INTERVAL_LENGTH):
        raise ValueError("The periodic interval must contain an integer number of cells.")
    return (np.arange(site_count) - site_count // 2) * spacing, spacing


def principal_matrix(angles):
    return BASE_SPEED * (
        (1 + np.cos(angles))[:, None, None] * SIGMA_Z
        + np.sin(angles)[:, None, None] * SIGMA_X
    ) / 2


def principal_derivative(positions):
    velocity, velocity_derivative, _ = profile(positions)
    angles = rotation_angle(velocity)
    angle_derivative = -2 * velocity_derivative / (BASE_SPEED * np.sqrt(1 - (velocity / BASE_SPEED)**2))
    return BASE_SPEED * angle_derivative[:, None, None] * (
        -np.sin(angles)[:, None, None] * SIGMA_Z
        + np.cos(angles)[:, None, None] * SIGMA_X
    ) / 2


def spectral_derivative(state, spacing):
    wave_numbers = 2 * math.pi * np.fft.fftfreq(len(state), d=spacing)
    return np.fft.ifft(1j * wave_numbers[:, None] * np.fft.fft(state, axis=0), axis=0)


def target_generator(state, positions, spacing, mass_frequency=0.0):
    angles = rotation_angle(profile(positions)[0])
    spatial_part = -np.einsum("nab,nb->na", principal_matrix(angles), spectral_derivative(state, spacing))
    connection_part = -0.5 * np.einsum("nab,nb->na", principal_derivative(positions), state)
    mass_part = -1j * np.asarray(mass_frequency)[..., None] * (state @ SIGMA_Y.T)
    return spatial_part + connection_part + mass_part


def periodic_envelope(coordinates):
    optical_length = INTERVAL_LENGTH / VELOCITY_SCALE
    wave_numbers = 2 * math.pi * np.arange(-50, 51) / optical_length
    coefficients = np.exp(-1.5**2 * (wave_numbers - 0.8)**2 / 2)
    normalization = math.sqrt(optical_length * float(np.sum(abs(coefficients)**2)))
    center = float(optical_coordinate(np.array([-4.0]))[0])
    return np.exp(1j * np.outer(coordinates - center, wave_numbers)) @ coefficients / normalization


def exact_massless_packet(positions, spacing, duration):
    velocity = profile(positions)[0]
    angles = rotation_angle(velocity)
    coordinates = optical_coordinate(positions)
    optical_state = np.stack((math.sqrt(0.7) * periodic_envelope(coordinates - duration),
                              1j * math.sqrt(0.3) * periodic_envelope(coordinates + duration)), axis=1)
    return math.sqrt(spacing) * rotate(optical_state, angles / 2) / np.sqrt(velocity)[:, None]


def packet_convergence():
    rows = []
    for time_step in TIME_STEPS:
        positions, spacing = spatial_grid(time_step)
        angles = rotation_angle(profile(positions)[0])
        initial = exact_massless_packet(positions, spacing, 0)
        normalization = np.linalg.norm(initial)
        initial /= normalization
        reference = exact_massless_packet(positions, spacing, FINAL_TIME) / normalization
        reference_density = np.sum(abs(reference)**2, axis=1)
        row = {"time_step": time_step, "spacing": spacing, "sites": len(positions), "macro_steps": round(FINAL_TIME / time_step)}
        for name, symmetric in (("ordered", False), ("symmetric", True)):
            evolved = evolve_walk(initial, angles, time_step, FINAL_TIME, symmetric=symmetric)
            density = np.sum(abs(evolved)**2, axis=1)
            row[name] = {
                "state_l2_error": float(np.linalg.norm(evolved - reference)),
                "density_total_variation": float(0.5 * np.sum(abs(density - reference_density))),
                "norm_squared": float(density.sum()),
            }
        row["reference_discrete_norm_squared"] = float(reference_density.sum())
        rows.append(row)
    return {
        "duration": FINAL_TIME,
        "initial_state": "Periodic Fourier packet in y(x), width parameter 1.5, central wave number 0.8; chiral probabilities 0.7 and 0.3.",
        "reference": "Exact massless characteristics in optical coordinate y with v^(-1/2) density weighting and local internal rotation.",
        "rows": rows,
        "ordered_error_step_slope": float(np.polyfit(np.log(TIME_STEPS), np.log([row["ordered"]["state_l2_error"] for row in rows]), 1)[0]),
        "symmetric_error_step_slope": float(np.polyfit(np.log(TIME_STEPS), np.log([row["symmetric"]["state_l2_error"] for row in rows]), 1)[0]),
    }


def generator_convergence(with_mass=False):
    rows = []
    for time_step in TIME_STEPS:
        positions, spacing = spatial_grid(time_step)
        angles = rotation_angle(profile(positions)[0])
        state = exact_massless_packet(positions, spacing, 0)
        state /= np.linalg.norm(state)
        masses = 0.7 * profile(positions)[0] if with_mass else 0.0
        forward = walk_step(state, angles, time_step, masses)
        backward = walk_step(state, angles, time_step, masses, inverse=True)
        estimate = (forward - backward) / (2 * time_step)
        expected = target_generator(state, positions, spacing, masses)
        rows.append({"time_step": time_step, "generator_l2_error": float(np.linalg.norm(estimate - expected))})
    return {
        "mass_frequency": "0.7*v(x)" if with_mass else "0",
        "rows": rows,
        "error_step_slope": float(np.polyfit(np.log(TIME_STEPS), np.log([row["generator_l2_error"] for row in rows]), 1)[0]),
    }


def metric_data(position, representative):
    velocity, velocity_derivative, _ = profile(np.array([position]))
    velocity, velocity_derivative = float(velocity[0]), float(velocity_derivative[0])
    if representative == "optical_flat":
        metric = np.diag([1.0, -1 / velocity**2])
        spatial_derivative = np.diag([0.0, 2 * velocity_derivative / velocity**3])
    elif representative == "unit_spatial":
        metric = np.diag([velocity**2, -1.0])
        spatial_derivative = np.diag([2 * velocity * velocity_derivative, 0.0])
    else:
        raise ValueError("Unknown conformal metric representative.")
    derivatives = np.zeros((2, 2, 2))
    derivatives[1] = spatial_derivative
    return metric, derivatives


def christoffel_symbols(position, representative):
    metric, derivatives = metric_data(position, representative)
    inverse_metric = np.linalg.inv(metric)
    connection = np.zeros((2, 2, 2))
    for upper in range(2):
        for first_lower in range(2):
            for second_lower in range(2):
                connection[upper, first_lower, second_lower] = 0.5 * sum(
                    inverse_metric[upper, contracted] * (
                        derivatives[first_lower, contracted, second_lower]
                        + derivatives[second_lower, contracted, first_lower]
                        - derivatives[contracted, first_lower, second_lower]
                    ) for contracted in range(2)
                )
    return connection


def ricci_from_connection(position, representative, difference_step=1e-4):
    connection = christoffel_symbols(position, representative)
    connection_derivative = (
        christoffel_symbols(position + difference_step, representative)
        - christoffel_symbols(position - difference_step, representative)
    ) / (2 * difference_step)
    ricci = np.zeros((2, 2))
    for first_index in range(2):
        for second_index in range(2):
            value = connection_derivative[1, first_index, second_index]
            if second_index == 1:
                value -= sum(connection_derivative[contracted, first_index, contracted] for contracted in range(2))
            value += sum(
                connection[outer, outer, inner] * connection[inner, first_index, second_index]
                - connection[outer, second_index, inner] * connection[inner, first_index, outer]
                for outer in range(2) for inner in range(2)
            )
            ricci[first_index, second_index] = value
    return ricci


def curvature_from_connection(position, representative, difference_step=1e-4):
    metric, _ = metric_data(position, representative)
    ricci = ricci_from_connection(position, representative, difference_step)
    return float(np.sum(np.linalg.inv(metric) * ricci))


def geometry_diagnostics():
    positions = np.array([-5.0, 0.0, 5.0, 10.0])
    velocity, _, second_derivative = profile(positions)
    rows = []
    for index, position in enumerate(positions):
        einstein_residuals = {}
        for representative in ("optical_flat", "unit_spatial"):
            metric, _ = metric_data(position, representative)
            ricci = ricci_from_connection(position, representative)
            scalar_curvature = float(np.sum(np.linalg.inv(metric) * ricci))
            einstein_residuals[representative] = float(np.max(abs(ricci - 0.5 * scalar_curvature * metric)))
        rows.append({
            "position": float(position),
            "coordinate_characteristic_speed": float(velocity[index]),
            "flat_representative_curvature": curvature_from_connection(position, "optical_flat"),
            "unit_spatial_representative_curvature": curvature_from_connection(position, "unit_spatial"),
            "unit_spatial_analytic_curvature": float(2 * second_derivative[index] / velocity[index]),
            "proper_clock_rate_flat": 1.0,
            "proper_clock_rate_unit_spatial": float(velocity[index]),
            "mass_coefficient_flat_for_fixed_proper_mass_0_7": 0.7,
            "mass_coefficient_unit_spatial_for_fixed_proper_mass_0_7": float(0.7 * velocity[index]),
            "two_dimensional_einstein_tensor_max_residual": einstein_residuals,
        })
    sample_positions, spacing = spatial_grid(0.1)
    state = exact_massless_packet(sample_positions, spacing, 0)
    state /= np.linalg.norm(state)
    angles = rotation_angle(profile(sample_positions)[0])
    naive = -np.einsum("nab,nb->na", principal_matrix(angles), spectral_derivative(state, spacing))
    corrected = target_generator(state, sample_positions, spacing)
    return {
        "metric_representatives": {
            "optical_flat": "ds0^2=dt^2-dx^2/v(x)^2=dt^2-dy^2, dy=dx/v(x)",
            "unit_spatial": "ds1^2=v(x)^2*dt^2-dx^2=v(x)^2*ds0^2",
            "shared_characteristics": "dx/dt=+v(x) and -v(x)",
        },
        "curvature_convention": "R_mn=d_l Gamma^l_mn-d_n Gamma^l_ml+Gamma^l_ls Gamma^s_mn-Gamma^l_ns Gamma^s_ml, signature (+,-).",
        "rows": rows,
        "naive_generator_probability_rate": float(2 * np.vdot(state, naive).real),
        "corrected_generator_probability_rate": float(2 * np.vdot(state, corrected).real),
        "ray_travel_time_minus4_to_plus4": float(optical_coordinate(np.array([4.0]))[0] - optical_coordinate(np.array([-4.0]))[0]),
        "interpretation": "Two externally chosen conformal representatives share the same massless densitized equation but have different curvature and proper-clock scales. A cone alone does not select a full metric.",
    }


class VariableGeometryTests(unittest.TestCase):
    def test_massive_step_converges_to_prescribed_mass_coefficient(self):
        result = generator_convergence(with_mass=True)
        self.assertAlmostEqual(result["error_step_slope"], 2, delta=0.15)
        errors = [row["generator_l2_error"] for row in result["rows"]]
        self.assertTrue(all(later < earlier for earlier, later in zip(errors, errors[1:])))

    def test_two_dimensional_einstein_tensor_vanishes_even_for_curved_representative(self):
        for position in (-5.0, 0.0, 3.0, 10.0):
            for representative in ("optical_flat", "unit_spatial"):
                metric, _ = metric_data(position, representative)
                ricci = ricci_from_connection(position, representative)
                scalar_curvature = float(np.sum(np.linalg.inv(metric) * ricci))
                np.testing.assert_allclose(ricci - 0.5 * scalar_curvature * metric, 0, atol=1e-8)

    def test_rotating_to_characteristic_basis_has_no_extra_massless_term(self):
        positions, spacing = spatial_grid(0.1)
        velocity, velocity_derivative, _ = profile(positions)
        angles = rotation_angle(velocity)
        state = exact_massless_packet(positions, spacing, 0.3)
        rotated_state = rotate(state, -angles / 2)
        expected = -velocity[:, None] * (spectral_derivative(rotated_state, spacing) @ SIGMA_Z.T)
        expected -= 0.5 * velocity_derivative[:, None] * (rotated_state @ SIGMA_Z.T)
        actual = rotate(target_generator(state, positions, spacing), -angles / 2)
        np.testing.assert_allclose(actual, expected, atol=1e-10)

    def test_connection_term_is_needed_for_norm_conservation(self):
        diagnostics = geometry_diagnostics()
        self.assertGreater(abs(diagnostics["naive_generator_probability_rate"]), 1e-4)
        self.assertLess(abs(diagnostics["corrected_generator_probability_rate"]), 1e-10)

    def test_metric_representatives_are_conformal_and_share_null_directions(self):
        for position in (-5.0, 0.0, 3.0, 10.0):
            velocity = float(profile(np.array([position]))[0][0])
            flat, _ = metric_data(position, "optical_flat")
            curved, _ = metric_data(position, "unit_spatial")
            np.testing.assert_allclose(curved, velocity**2 * flat, atol=1e-12)
            for direction in (np.array([1, velocity]), np.array([1, -velocity])):
                self.assertAlmostEqual(float(direction @ flat @ direction), 0, places=12)
                self.assertAlmostEqual(float(direction @ curved @ direction), 0, places=12)

    def test_curvature_calculation_distinguishes_shared_cones(self):
        for position in (-5.0, 0.0, 3.0, 10.0):
            velocity, _, second_derivative = profile(np.array([position]))
            expected = float(2 * second_derivative[0] / velocity[0])
            for step in (1e-4, 5e-5):
                self.assertAlmostEqual(curvature_from_connection(position, "optical_flat", step), 0, delta=1e-8)
                self.assertAlmostEqual(curvature_from_connection(position, "unit_spatial", step), expected, delta=1e-8)
        self.assertGreater(abs(curvature_from_connection(0, "unit_spatial")), 0.03)

    def test_fixed_proper_mass_can_distinguish_the_two_representatives(self):
        positions, spacing = spatial_grid(0.1)
        state = exact_massless_packet(positions, spacing, 0)
        velocity = profile(positions)[0]
        flat = target_generator(state, positions, spacing, 0.7)
        unit_spatial = target_generator(state, positions, spacing, 0.7 * velocity)
        expected_difference = -1j * (0.7 * (velocity - 1))[:, None] * (state @ SIGMA_Y.T)
        np.testing.assert_allclose(unit_spatial - flat, expected_difference, atol=1e-12)
        self.assertGreater(float(np.linalg.norm(unit_spatial - flat)), 0.1)

    def test_variable_rotations_preserve_norm_and_have_inverse(self):
        generator = np.random.default_rng(401)
        state = generator.normal(size=(31, 2)) + 1j * generator.normal(size=(31, 2))
        state /= np.linalg.norm(state)
        angles = generator.uniform(0.5, 2.3, len(state))
        masses = generator.uniform(0.1, 0.8, len(state))
        for symmetric in (False, True):
            updated = walk_step(state, angles, 0.1, masses, symmetric)
            restored = walk_step(updated, angles, 0.1, masses, symmetric, inverse=True)
            self.assertAlmostEqual(float(np.linalg.norm(updated)), 1, places=12)
            np.testing.assert_allclose(restored, state, atol=1e-12)

    def test_one_macro_step_has_four_cell_support(self):
        site_count, steps = 161, 10
        origin = site_count // 2
        state = np.zeros((site_count, 2), dtype=complex)
        state[origin] = [1 / math.sqrt(2), 1j / math.sqrt(2)]
        angles = 1.2 + 0.2 * np.sin(2 * math.pi * np.arange(site_count) / site_count)
        for _ in range(steps):
            state = walk_step(state, angles, 0.1)
        probabilities = np.sum(abs(state)**2, axis=1)
        self.assertEqual(float(probabilities[abs(np.arange(site_count) - origin) > 4 * steps].sum()), 0)
        self.assertAlmostEqual(float(probabilities.sum()), 1, places=12)

    def test_characteristic_speeds_are_plus_and_minus_profile(self):
        positions, _ = spatial_grid(0.2)
        velocity = profile(positions)[0]
        matrices = principal_matrix(rotation_angle(velocity))
        np.testing.assert_allclose(np.linalg.eigvalsh(matrices), np.stack((-velocity, velocity), axis=1), atol=1e-12)

    def test_optical_coordinate_derivative_matches_inverse_velocity(self):
        positions = np.linspace(-8, 8, 17)
        step = 1e-5
        derivative = (optical_coordinate(positions + step) - optical_coordinate(positions - step)) / (2 * step)
        np.testing.assert_allclose(derivative, 1 / profile(positions)[0], atol=1e-9)

    def test_generator_matches_principal_and_half_derivative_connection(self):
        result = generator_convergence()
        self.assertAlmostEqual(result["error_step_slope"], 2, delta=0.12)
        errors = [row["generator_l2_error"] for row in result["rows"]]
        self.assertTrue(all(later < earlier for earlier, later in zip(errors, errors[1:])))

    def test_characteristic_solution_satisfies_continuum_equation(self):
        positions, spacing = spatial_grid(0.1)
        step = 1e-5
        state = exact_massless_packet(positions, spacing, 0.4)
        time_derivative = (exact_massless_packet(positions, spacing, 0.4 + step)
                           - exact_massless_packet(positions, spacing, 0.4 - step)) / (2 * step)
        np.testing.assert_allclose(time_derivative, target_generator(state, positions, spacing), atol=1e-8)

    def test_packet_converges_to_independent_characteristics(self):
        result = packet_convergence()
        self.assertAlmostEqual(result["ordered_error_step_slope"], 1, delta=0.12)
        self.assertAlmostEqual(result["symmetric_error_step_slope"], 2, delta=0.12)
        for row in result["rows"]:
            self.assertAlmostEqual(row["reference_discrete_norm_squared"], 1, places=10)
            for name in ("ordered", "symmetric"):
                self.assertAlmostEqual(row[name]["norm_squared"], 1, places=10)

    def test_variable_mass_changes_no_principal_eigenvalues(self):
        positions, spacing = spatial_grid(0.05)
        velocity = profile(positions)[0]
        angles = rotation_angle(velocity)
        state = exact_massless_packet(positions, spacing, 0)
        masses = 0.4 + 0.2 * np.cos(2 * math.pi * positions / INTERVAL_LENGTH)
        difference = target_generator(state, positions, spacing, masses) - target_generator(state, positions, spacing)
        expected = -1j * masses[:, None] * (state @ SIGMA_Y.T)
        np.testing.assert_allclose(difference, expected, atol=1e-12)
        np.testing.assert_allclose(np.linalg.eigvalsh(principal_matrix(angles)), np.stack((-velocity, velocity), axis=1), atol=1e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(VariableGeometryTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {
        "scope": "Externally prescribed smooth variable-speed quantum walk; background kinematics, not self-consistent gravity.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "units_and_inputs": {"hbar": 1, "base_speed": BASE_SPEED, "interval_length": INTERVAL_LENGTH,
                             "spacing_rule": "a=base_speed*time_step/4", "velocity_scale": VELOCITY_SCALE,
                             "profile_amplitude": PROFILE_AMPLITUDE,
                             "profile": "v(x)=velocity_scale/(1+profile_amplitude*cos(2*pi*x/interval_length))"},
        "criteria_before_first_run": [
            "Local rotations and shifts must be exactly unitary even for position-dependent angles.",
            "The continuum generator must include A*d_x+(d_x A)/2, not just a variable multiplier of the derivative.",
            "Characteristic speeds must equal eigenvalues plus/minus v(x), independent of an onsite mass term.",
            "Massless packet evolution must converge to an independent characteristic solution, not only another matrix update.",
            "A position-dependent coordinate speed will not be treated as proof of nonzero curvature or of a sourced gravity equation.",
        ],
        "generator_convergence": generator_convergence(),
        "prescribed_mass_generator_convergence": generator_convergence(with_mass=True),
        "packet_convergence": packet_convergence(),
        "geometry_diagnostics": geometry_diagnostics(),
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    encoded = json.dumps(results, indent=2)
    if arguments.write_results:
        output = Path(__file__).resolve().parent / "variable_geometry_results.json"
        output.write_text(encoded + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    print(encoded)


if __name__ == "__main__":
    main()