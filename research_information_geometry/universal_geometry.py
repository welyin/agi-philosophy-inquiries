"""Test when two Dirac-type probes admit one minimally coupled static metric."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp

from backreaction import centered_derivative
from relativistic_limit import SIGMA_X, SIGMA_Z


LENGTH = 20.0
BACKGROUND_AMPLITUDE = 0.2
BASE_SPEED = 0.7


@dataclass(frozen=True)
class Probe:
    name: str
    proper_mass: float
    lapse_charge: float
    spatial_charge: float

    def __post_init__(self):
        if self.proper_mass <= 0:
            raise ValueError("A calibrated nonzero proper mass is required for this metric reconstruction.")
        if not all(math.isfinite(value) for value in (
            self.proper_mass, self.lapse_charge, self.spatial_charge
        )):
            raise ValueError("Probe parameters must be finite.")


REFERENCE = Probe("A", 0.6, 1.0, 0.0)
CASES = (
    ("universal", Probe("B", 1.1, 1.0, 0.0)),
    ("same_cones_different_clocks", Probe("B", 1.1, 1.5, 0.5)),
    ("same_clocks_different_cones", Probe("B", 1.1, 1.0, 0.4)),
)


def background(positions):
    wave_number = 2 * math.pi / LENGTH
    positions = np.asarray(positions)
    return (
        BACKGROUND_AMPLITUDE * np.cos(wave_number * positions),
        -BACKGROUND_AMPLITUDE * wave_number * np.sin(wave_number * positions),
    )


def coefficients(field, probe):
    lapse = np.exp(probe.lapse_charge * np.asarray(field))
    spatial_scale = np.exp(probe.spatial_charge * np.asarray(field)) / BASE_SPEED
    speed = lapse / spatial_scale
    mass_coefficient = probe.proper_mass * lapse
    return speed, mass_coefficient


def reconstructed_metric(field, probe):
    speed, mass_coefficient = coefficients(field, probe)
    lapse = mass_coefficient / probe.proper_mass
    spatial_scale = lapse / speed
    metric = np.zeros(np.shape(lapse) + (2, 2))
    metric[..., 0, 0] = lapse**2
    metric[..., 1, 1] = -spatial_scale**2
    return metric


def comparison_residuals(field, first, second):
    first_speed, first_mass = coefficients(field, first)
    second_speed, second_mass = coefficients(field, second)
    cone_residual = np.log(second_speed / first_speed)
    clock_residual = np.log((second_mass / second.proper_mass) / (first_mass / first.proper_mass))
    return {
        "max_absolute_log_cone_ratio": float(np.max(abs(cone_residual))),
        "max_absolute_log_normalized_mass_ratio": float(np.max(abs(clock_residual))),
        "interpretation": "Coefficient agreement in a fixed static diagonal chart with constant calibrated proper masses and no additional forces. Not a test of all possible metric-plus-matter theories.",
    }


def apply_hamiltonian(state, field, probe, spacing):
    speed, mass_coefficient = coefficients(field, probe)
    derivative = speed[:, None] * centered_derivative(state, spacing)
    derivative += centered_derivative(speed[:, None] * state, spacing)
    return -0.5j * (derivative @ SIGMA_Z.T) + mass_coefficient[:, None] * (state @ SIGMA_X.T)


def local_dispersion(position, momentum, probe):
    field, _ = background(position)
    speed, mass_coefficient = coefficients(field, probe)
    return float(math.hypot(float(speed) * momentum, float(mass_coefficient)))


def ray_derivative(position, momentum, probe):
    field, field_derivative = background(position)
    speed, mass_coefficient = coefficients(field, probe)
    speed_derivative = speed * (probe.lapse_charge - probe.spatial_charge) * field_derivative
    mass_derivative = mass_coefficient * probe.lapse_charge * field_derivative
    energy = math.hypot(float(speed) * momentum, float(mass_coefficient))
    position_rate = speed**2 * momentum / energy
    momentum_rate = -(speed * speed_derivative * momentum**2 + mass_coefficient * mass_derivative) / energy
    return float(position_rate), float(momentum_rate)


def rest_acceleration(position, probe):
    field, field_derivative = background(position)
    speed, _ = coefficients(field, probe)
    return float(-speed**2 * probe.lapse_charge * field_derivative)


def initial_momentum_for_velocity(position, coordinate_velocity, probe):
    field, _ = background(position)
    speed, mass_coefficient = coefficients(field, probe)
    if abs(coordinate_velocity) >= speed:
        raise ValueError("The initial ray must be timelike for this probe.")
    return float(mass_coefficient * coordinate_velocity / (speed * math.sqrt(float(speed**2 - coordinate_velocity**2))))


def integrate_ray(probe, duration=4.0, initial_position=LENGTH / 4, initial_velocity=0.0):
    momentum = initial_momentum_for_velocity(initial_position, initial_velocity, probe)

    def right_hand_side(time, variables):
        return ray_derivative(variables[0], variables[1], probe)

    solution = solve_ivp(right_hand_side, (0.0, duration), (initial_position, momentum),
                         method="DOP853", rtol=1e-11, atol=1e-13,
                         t_eval=np.linspace(0, duration, 81))
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def geodesic_acceleration(position, coordinate_velocity, probe):
    field, field_derivative = background(position)
    speed, _ = coefficients(field, probe)
    gamma_x_tt = speed**2 * probe.lapse_charge * field_derivative
    gamma_x_xx = probe.spatial_charge * field_derivative
    gamma_t_tx = probe.lapse_charge * field_derivative
    return float(-gamma_x_tt + (2 * gamma_t_tx - gamma_x_xx) * coordinate_velocity**2)


def integrate_geodesic(probe, duration=4.0, initial_position=LENGTH / 4, initial_velocity=0.0):
    def right_hand_side(time, variables):
        return variables[1], geodesic_acceleration(variables[0], variables[1], probe)

    solution = solve_ivp(right_hand_side, (0.0, duration), (initial_position, initial_velocity),
                         method="RK45", rtol=1e-11, atol=1e-13,
                         t_eval=np.linspace(0, duration, 81))
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def null_travel_time(probe, source=0.0, target=LENGTH / 4):
    def integrand(position):
        speed, _ = coefficients(background(position)[0], probe)
        return 1 / float(speed)

    duration, error = quad(integrand, source, target, epsabs=1e-12, epsrel=1e-12)
    return {"duration": float(duration), "quadrature_error_estimate": float(error)}


def trajectory_diagnostics():
    first = integrate_ray(REFERENCE)
    rows = []
    for name, second in CASES:
        ray = integrate_ray(second)
        geometric = integrate_geodesic(second)
        energy = np.array([local_dispersion(position, momentum, second) for position, momentum in ray.y.T])
        rows.append({
            "case": name,
            "initial_position": LENGTH / 4,
            "initial_coordinate_velocity": 0.0,
            "duration": 4.0,
            "final_position_A": float(first.y[0, -1]),
            "final_position_B": float(ray.y[0, -1]),
            "final_position_difference_B_minus_A": float(ray.y[0, -1] - first.y[0, -1]),
            "max_sampled_position_difference": float(np.max(abs(ray.y[0] - first.y[0]))),
            "max_sampled_hamiltonian_ray_vs_geodesic_position_difference": float(np.max(abs(ray.y[0] - geometric.y[0]))),
            "max_sampled_absolute_ray_energy_drift": float(np.max(abs(energy - energy[0]))),
            "null_travel_time_A_0_to_5": null_travel_time(REFERENCE),
            "null_travel_time_B_0_to_5": null_travel_time(second),
        })
    return {
        "rows": rows,
        "methods": "Positive local mass-shell Hamilton equations (DOP853), independently checked against static-metric coordinate-time geodesics (RK45).",
        "scope": "Leading geometric-optics/point-particle diagnostics. No exact finite-width quantum-packet free-fall simulation or laboratory data.",
    }


def static_unitarity_diagnostics():
    sites, spacing = 24, LENGTH / 24
    field, _ = background(np.arange(sites) * spacing)
    dimension = 2 * sites
    basis = np.eye(dimension, dtype=complex)
    generator = np.random.default_rng(7071)
    state = generator.normal(size=dimension) + 1j * generator.normal(size=dimension)
    state /= np.linalg.norm(state)
    rows = []
    for name, probe in (("reference_A", REFERENCE),) + CASES:
        hamiltonian = np.column_stack([
            apply_hamiltonian(basis[:, index].reshape(sites, 2), field, probe, spacing).ravel()
            for index in range(dimension)
        ])
        spectrum, eigenvectors = np.linalg.eigh(hamiltonian)
        evolved = eigenvectors @ (np.exp(-2.3j * spectrum) * (eigenvectors.conj().T @ state))
        rows.append({
            "case": name,
            "probe": asdict(probe),
            "hermiticity_residual": float(np.linalg.norm(hamiltonian - hamiltonian.conj().T)),
            "norm_squared_drift": float(abs(np.vdot(evolved, evolved) - np.vdot(state, state))),
            "static_energy_drift": float(abs(np.vdot(evolved, hamiltonian @ evolved) - np.vdot(state, hamiltonian @ state))),
        })
    return {"sites": sites, "duration": 2.3, "rows": rows,
            "meaning": "Every probe has ordinary unitary, static-energy-conserving dynamics, including nonuniversal cases. Coarse finite matrices check algebra only, not a continuum error bound."}


def sensitivity_diagnostics():
    positions = np.array([0.0, LENGTH / 4, LENGTH / 2])
    values, _ = background(positions)
    cone_design = np.array([[field, -field] for field in values])
    clock_design = np.array([[field, 0.0] for field in values])
    combined = np.concatenate((cone_design, clock_design), axis=0)
    error = 0.001
    max_field = float(np.max(abs(values)))
    return {
        "unknown_differences": ["delta_alpha=alpha_B-alpha_A", "delta_beta=beta_B-beta_A"],
        "data": "log(v_B/v_A)=(delta_alpha-delta_beta)*q; log[(mu_B/m_B)/(mu_A/m_A)]=delta_alpha*q.",
        "background_values": values.tolist(),
        "cone_only_rank": int(np.linalg.matrix_rank(cone_design, tol=1e-12)),
        "clock_only_rank": int(np.linalg.matrix_rank(clock_design, tol=1e-12)),
        "combined_rank": int(np.linalg.matrix_rank(combined, tol=1e-12)),
        "combined_singular_values": np.linalg.svd(combined, compute_uv=False).tolist(),
        "illustration_not_experimental_limits": {
            "assumed_absolute_log_ratio_bounds": error,
            "known_max_absolute_background": max_field,
            "absolute_delta_alpha_bound": error / max_field,
            "absolute_delta_alpha_minus_delta_beta_bound": error / max_field,
            "conservative_absolute_delta_beta_bound": 2 * error / max_field,
        },
        "calibration_caveat": "Charge bounds require an independently normalized and measured q. With q unknown, experiments constrain products and probe ratios, not alpha and beta separately. The numerical bound is illustrative and is not taken from MICROSCOPE or clock experiments.",
    }


def coefficient_diagnostics():
    positions = np.array([0.0, LENGTH / 4, LENGTH / 2])
    field, _ = background(positions)
    outputs = {}
    for name, second in CASES:
        first_speed, first_mass = coefficients(field, REFERENCE)
        second_speed, second_mass = coefficients(field, second)
        first_acceleration = rest_acceleration(LENGTH / 4, REFERENCE)
        second_acceleration = rest_acceleration(LENGTH / 4, second)
        outputs[name] = {
            "probe_A": asdict(REFERENCE),
            "probe_B": asdict(second),
            "coefficient_residuals": comparison_residuals(field, REFERENCE, second),
            "positions": positions.tolist(),
            "background_values": field.tolist(),
            "speed_A": first_speed.tolist(),
            "speed_B": second_speed.tolist(),
            "normalized_mass_A": (first_mass / REFERENCE.proper_mass).tolist(),
            "normalized_mass_B": (second_mass / second.proper_mass).tolist(),
            "clock_ratio_B_over_A": ((second_mass / second.proper_mass) / (first_mass / REFERENCE.proper_mass)).tolist(),
            "initial_rest_acceleration_A": first_acceleration,
            "initial_rest_acceleration_B": second_acceleration,
            "coordinate_acceleration_contrast_at_release": 2 * abs(second_acceleration - first_acceleration) / (abs(first_acceleration) + abs(second_acceleration)),
            "probe_metric_A": reconstructed_metric(field, REFERENCE).tolist(),
            "probe_metric_B": reconstructed_metric(field, second).tolist(),
        }
    return {
        "background": {"length": LENGTH, "amplitude": BACKGROUND_AMPLITUDE,
                       "formula": "q(x)=0.2*cos(2*pi*x/20)", "external_and_static": True},
        "reference_calibration": "x=L/4, where q=0, normalized mass coefficient=1, and all coordinate speeds=0.7.",
        "cases": outputs,
    }


class UniversalGeometryTests(unittest.TestCase):
    def test_shared_metric_rays_ignore_proper_mass_at_matched_initial_velocity(self):
        heavier = replace(REFERENCE, proper_mass=1.7)
        for initial_velocity in (0.0, 0.15):
            first = integrate_ray(REFERENCE, initial_velocity=initial_velocity)
            second = integrate_ray(heavier, initial_velocity=initial_velocity)
            np.testing.assert_allclose(first.y[0], second.y[0], atol=1e-9, rtol=0)
            np.testing.assert_allclose(first.y[1] / REFERENCE.proper_mass, second.y[1] / heavier.proper_mass, atol=1e-9, rtol=0)

    def test_hamiltonian_rays_match_independent_geodesic_integration(self):
        for _, probe in CASES:
            first = integrate_ray(probe, initial_velocity=0.1)
            second = integrate_geodesic(probe, initial_velocity=0.1)
            np.testing.assert_allclose(first.y[0], second.y[0], atol=1e-8, rtol=0)
            first_velocities = np.array([ray_derivative(position, momentum, probe)[0] for position, momentum in first.y.T])
            np.testing.assert_allclose(first_velocities, second.y[1], atol=1e-8, rtol=0)

    def test_same_cones_can_have_different_massive_release_trajectories(self):
        first = integrate_ray(REFERENCE)
        second = integrate_ray(CASES[1][1])
        self.assertGreater(float(abs(first.y[0, -1] - second.y[0, -1])), 0.05)
        self.assertAlmostEqual(null_travel_time(REFERENCE)["duration"], null_travel_time(CASES[1][1])["duration"], places=11)

    def test_null_and_rest_tests_probe_distinct_coefficients(self):
        second = CASES[2][1]
        self.assertAlmostEqual(rest_acceleration(LENGTH / 4, REFERENCE), rest_acceleration(LENGTH / 4, second))
        self.assertGreater(abs(null_travel_time(REFERENCE)["duration"] - null_travel_time(second)["duration"]), 0.1)

    def test_static_energy_and_norm_do_not_enforce_universality(self):
        diagnostics = static_unitarity_diagnostics()
        for row in diagnostics["rows"]:
            self.assertLess(row["hermiticity_residual"], 1e-12)
            self.assertLess(row["norm_squared_drift"], 1e-12)
            self.assertLess(row["static_energy_drift"], 1e-11)

    def test_two_normalized_ratios_identify_two_charge_differences(self):
        diagnostics = sensitivity_diagnostics()
        self.assertEqual(diagnostics["cone_only_rank"], 1)
        self.assertEqual(diagnostics["clock_only_rank"], 1)
        self.assertEqual(diagnostics["combined_rank"], 2)
        bound = diagnostics["illustration_not_experimental_limits"]
        self.assertAlmostEqual(bound["absolute_delta_alpha_bound"], 0.005)
        self.assertAlmostEqual(bound["conservative_absolute_delta_beta_bound"], 0.01)

    def test_nonuniversal_conformal_case_can_be_rewritten_with_variable_proper_mass(self):
        second = CASES[1][1]
        fields = np.linspace(-0.2, 0.2, 17)
        first_speed, first_mass = coefficients(fields, REFERENCE)
        second_speed, second_mass = coefficients(fields, second)
        common_lapse = first_mass / REFERENCE.proper_mass
        effective_mass = second.proper_mass * np.exp((second.lapse_charge - REFERENCE.lapse_charge) * fields)
        np.testing.assert_allclose(second_speed, first_speed, atol=1e-12)
        np.testing.assert_allclose(second_mass, effective_mass * common_lapse, atol=1e-12)
        self.assertGreater(float(np.ptp(effective_mass)), 0.1)

    def test_different_proper_masses_allow_same_metric(self):
        field = np.linspace(-0.2, 0.2, 17)
        second = CASES[0][1]
        np.testing.assert_allclose(reconstructed_metric(field, REFERENCE), reconstructed_metric(field, second), atol=1e-12)
        first_speed, first_mass = coefficients(field, REFERENCE)
        second_speed, second_mass = coefficients(field, second)
        np.testing.assert_allclose(first_speed, second_speed, atol=1e-12)
        self.assertGreater(float(np.max(abs(first_mass - second_mass))), 0.4)

    def test_equal_cones_do_not_force_equal_clock_normalization(self):
        second = CASES[1][1]
        field = np.linspace(-0.2, 0.2, 17)
        residuals = comparison_residuals(field, REFERENCE, second)
        self.assertLess(residuals["max_absolute_log_cone_ratio"], 1e-12)
        self.assertAlmostEqual(residuals["max_absolute_log_normalized_mass_ratio"], 0.1)
        first_metric = reconstructed_metric(field, REFERENCE)
        second_metric = reconstructed_metric(field, second)
        np.testing.assert_allclose(second_metric, np.exp(field)[:, None, None] * first_metric, atol=1e-12)

    def test_equal_clocks_do_not_force_equal_cones(self):
        residuals = comparison_residuals(np.linspace(-0.2, 0.2, 17), REFERENCE, CASES[2][1])
        self.assertLess(residuals["max_absolute_log_normalized_mass_ratio"], 1e-12)
        self.assertAlmostEqual(residuals["max_absolute_log_cone_ratio"], 0.08)

    def test_reference_calibration_cannot_hide_spatially_varying_mismatch(self):
        field_at_reference, _ = background(LENGTH / 4)
        for _, probe in CASES:
            speed, mass_coefficient = coefficients(field_at_reference, probe)
            self.assertAlmostEqual(float(speed), BASE_SPEED)
            self.assertAlmostEqual(float(mass_coefficient / probe.proper_mass), 1)
        _, reference_mass = coefficients(0.2, REFERENCE)
        _, different_mass = coefficients(0.2, CASES[1][1])
        self.assertGreater(abs(float(different_mass / CASES[1][1].proper_mass - reference_mass / REFERENCE.proper_mass)), 0.1)

    def test_every_case_is_hermitian_despite_metric_mismatch(self):
        generator = np.random.default_rng(707)
        sites, spacing = 40, LENGTH / 40
        field, _ = background(np.arange(sites) * spacing)
        first = generator.normal(size=(sites, 2)) + 1j * generator.normal(size=(sites, 2))
        second = generator.normal(size=(sites, 2)) + 1j * generator.normal(size=(sites, 2))
        for probe in (REFERENCE,) + tuple(probe for _, probe in CASES):
            first_pair = np.vdot(first, apply_hamiltonian(second, field, probe, spacing))
            second_pair = np.vdot(apply_hamiltonian(first, field, probe, spacing), second)
            self.assertAlmostEqual(abs(first_pair - second_pair), 0, delta=1e-11)

    def test_geometric_mass_shell_matches_local_dirac_spectrum(self):
        for probe in (REFERENCE,) + tuple(probe for _, probe in CASES):
            for position in (0.0, 3.0, 7.0):
                field, _ = background(position)
                speed, mass_coefficient = coefficients(field, probe)
                momentum = 0.4
                energy = local_dispersion(position, momentum, probe)
                hamiltonian = speed * momentum * SIGMA_Z + mass_coefficient * SIGMA_X
                np.testing.assert_allclose(np.linalg.eigvalsh(hamiltonian), [-energy, energy], atol=1e-12)
                inverse_metric = np.linalg.inv(reconstructed_metric(field, probe))
                covector = np.array([energy, momentum])
                self.assertAlmostEqual(float(covector @ inverse_metric @ covector), probe.proper_mass**2, places=11)

    def test_ray_derivative_matches_hamiltonian_finite_differences(self):
        step = 1e-5
        for _, probe in CASES:
            position, momentum = 3.0, 0.4
            position_rate, momentum_rate = ray_derivative(position, momentum, probe)
            position_difference = (local_dispersion(position, momentum + step, probe)
                                   - local_dispersion(position, momentum - step, probe)) / (2 * step)
            momentum_difference = -(local_dispersion(position + step, momentum, probe)
                                    - local_dispersion(position - step, momentum, probe)) / (2 * step)
            self.assertAlmostEqual(position_rate, position_difference, delta=1e-9)
            self.assertAlmostEqual(momentum_rate, momentum_difference, delta=1e-9)

    def test_rest_acceleration_is_mass_independent_only_with_fixed_couplings(self):
        position = LENGTH / 4
        other_mass = replace(REFERENCE, proper_mass=1.7)
        self.assertAlmostEqual(rest_acceleration(position, REFERENCE), rest_acceleration(position, other_mass))
        first = rest_acceleration(position, REFERENCE)
        second = rest_acceleration(position, CASES[1][1])
        self.assertAlmostEqual(second / first, 1.5)
        self.assertAlmostEqual(2 * abs(second - first) / (abs(first) + abs(second)), 0.4)
        for probe in (REFERENCE, other_mass, CASES[1][1]):
            field, _ = background(position)
            speed, mass_coefficient = coefficients(field, probe)
            momentum_rate = ray_derivative(position, 0.0, probe)[1]
            self.assertAlmostEqual(float(speed**2 / mass_coefficient * momentum_rate), rest_acceleration(position, probe))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(UniversalGeometryTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = {
        "scope": "Two prescribed Dirac-type probe species on one external static scalar background; a restricted common-metric compatibility test, not a derivation of universal coupling or gravity.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "units": "hbar=c_reference=1; common external coordinates and clocks; proper masses are independent fixed calibration inputs.",
        "criteria_before_first_run": [
            "A single static diagonal minimally coupled metric requires equal v_s and equal mu_s/m_s for all massive probes in the same calibration.",
            "The same scalar background and Hermitian probe generators are allowed to fail this common-metric requirement.",
            "Equal characteristic cones alone must not be called equal proper-clock or massive-probe geometry.",
            "Ray equations are leading Hamilton-Jacobi diagnostics, not an exact quantum wave-packet/free-fall experiment.",
        ],
        "coefficient_diagnostics": coefficient_diagnostics(),
        "trajectory_diagnostics": trajectory_diagnostics(),
        "static_unitarity_diagnostics": static_unitarity_diagnostics(),
        "sensitivity_diagnostics": sensitivity_diagnostics(),
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    encoded = json.dumps(results, indent=2)
    if arguments.write_results:
        output = Path(__file__).resolve().parent / "universal_geometry_results.json"
        output.write_text(encoded + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    compact = {
        "automated_checks": results["automated_checks"],
        "coefficient_cases": {
            name: {key: case[key] for key in ("coefficient_residuals", "clock_ratio_B_over_A",
                                            "initial_rest_acceleration_A", "initial_rest_acceleration_B",
                                            "coordinate_acceleration_contrast_at_release")}
            for name, case in results["coefficient_diagnostics"]["cases"].items()
        },
        "trajectory_diagnostics": results["trajectory_diagnostics"],
        "sensitivity_diagnostics": results["sensitivity_diagnostics"],
        "scope": results["scope"],
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()