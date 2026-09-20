"""Construction 09: a declared four-dimensional weak tensor metric contract.

The quadratic metric action, dimension, coupling and prescribed source are
inputs. Linearized constraints and probe predictions are conditional outputs,
not a derivation of gravity from cognition or a full nonlinear completion.
"""

import argparse
import json
import math
import unittest
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from baseline import IDENTITY, real_lift
from clock_calibration import clock_effect, clock_gate
from common_orientation_structure import encode_state
from massive_envelope import massive_hopping, postprocessing_table
from quantum_chain_background import initial_data
from background_field_limit import harmonic_flow, ground_covariance


ETA = np.diag([1., -1., -1., -1.])


def trace_reverse(tensor):
    return tensor - ETA * np.einsum("ab,ab->", ETA, tensor) / 2


def linear_einstein_symbol(covector, tensor):
    covector, tensor = np.asarray(covector), np.asarray(tensor)
    raised = ETA @ covector
    invariant = covector @ raised
    trace = np.einsum("ab,ab->", ETA, tensor)
    divergence = raised @ tensor
    double = raised @ tensor @ raised
    return (invariant * tensor + np.outer(covector, covector) * trace
            - np.outer(covector, divergence) - np.outer(divergence, covector)
            + ETA * (double - invariant * trace)) / 2


def harmonic_response(covector, source, gravitational_constant):
    covector, source = np.asarray(covector), np.asarray(source)
    invariant = covector @ ETA @ covector
    if (source.shape != (4, 4) or not np.allclose(source, source.T, rtol=0, atol=1e-12)
            or np.linalg.norm((ETA @ covector) @ source) > 1e-12 or abs(invariant) < 1e-12
            or not np.isfinite(gravitational_constant) or gravitational_constant <= 0):
        raise ValueError("Use a symmetric conserved source, positive G and a non-null Fourier mode.")
    return trace_reverse(16 * math.pi * gravitational_constant * source / invariant)


def quadratic_action_symbol(covector, tensor, source, gravitational_constant):
    raised = ETA @ tensor @ ETA
    return (np.sum(raised * linear_einstein_symbol(covector, tensor)) / (32 * math.pi * gravitational_constant)
            - np.sum(raised * source) / 2)


@dataclass(frozen=True)
class GaussianSource:
    mass: float = 1.0
    width: float = 1.5
    gravitational_constant: float = 0.02

    def __post_init__(self):
        if (not all(np.isfinite(value) and value > 0 for value in (self.mass, self.width, self.gravitational_constant))
                or self.peak_potential >= 0.05):
            raise ValueError("Use positive finite parameters and the declared weak-field peak below 0.05.")

    @property
    def peak_potential(self):
        return self.gravitational_constant * self.mass * math.sqrt(2 / math.pi) / self.width

    def density(self, radius):
        radius = np.asarray(radius)
        return self.mass * np.exp(-radius ** 2 / (2 * self.width ** 2)) / (2 * math.pi * self.width ** 2) ** 1.5

    def radial_jet(self, radius):
        if not np.isfinite(radius) or radius < 0:
            raise ValueError("Use a finite nonnegative radius.")
        strength = self.gravitational_constant * self.mass
        coefficient = math.sqrt(2 / math.pi) / self.width
        if radius < 1e-3 * self.width:
            scaled = radius / self.width
            value = -strength * coefficient * (1 - scaled ** 2 / 6 + scaled ** 4 / 40 - scaled ** 6 / 336)
            first = strength * coefficient / self.width * (scaled / 3 - scaled ** 3 / 10 + scaled ** 5 / 56)
            second = strength * coefficient / self.width ** 2 * (1 / 3 - 3 * scaled ** 2 / 10 + 5 * scaled ** 4 / 56)
            return value, first, second
        enclosed = math.erf(radius / (math.sqrt(2) * self.width))
        derivative = coefficient * math.exp(-radius ** 2 / (2 * self.width ** 2))
        value = -strength * enclosed / radius
        first = strength * (enclosed / radius ** 2 - derivative / radius)
        second = strength * (derivative / self.width ** 2 + 2 * derivative / radius ** 2 - 2 * enclosed / radius ** 3)
        return value, first, second

    def cartesian_jet(self, position):
        position = np.asarray(position, dtype=float)
        if position.shape != (3,) or not np.all(np.isfinite(position)):
            raise ValueError("Use a finite three-dimensional position.")
        radius = np.linalg.norm(position)
        value, first, second = self.radial_jet(radius)
        if radius == 0:
            return value, np.zeros(3), second * np.eye(3)
        direction = position / radius
        hessian = first / radius * np.eye(3) + (second - first / radius) * np.outer(direction, direction)
        return value, first * direction, hessian


def linear_einstein_from_hessian(second):
    ricci = np.zeros((4, 4))
    for first in range(4):
        for second_index in range(4):
            ricci[first, second_index] = sum(
                ETA[axis, axis] * (second[axis, first, axis, second_index]
                                  + second[axis, second_index, axis, first]
                                  - second[axis, axis, first, second_index]
                                  - second[first, second_index, axis, axis]) / 2
                for axis in range(4))
    return ricci - ETA * np.einsum("ab,ab->", ETA, ricci) / 2


def weak_metric_jet(position, source=GaussianSource()):
    potential, gradient, hessian = source.cartesian_jet(position)
    metric = ETA + 2 * potential * np.eye(4)
    first = np.zeros((4, 4, 4))
    second = np.zeros((4, 4, 4, 4))
    for axis in range(3):
        first[axis + 1] = 2 * gradient[axis] * np.eye(4)
        for other in range(3):
            second[axis + 1, other + 1] = 2 * hessian[axis, other] * np.eye(4)
    return metric, first, second


def nonlinear_einstein_from_jet(metric, first, second):
    inverse = np.linalg.inv(metric)
    inverse_first = np.array([-inverse @ derivative @ inverse for derivative in first])
    bracket = np.zeros((4, 4, 4))
    bracket_first = np.zeros((4, 4, 4, 4))
    for lower in range(4):
        for first_index in range(4):
            for second_index in range(4):
                bracket[lower, first_index, second_index] = (first[first_index, lower, second_index]
                    + first[second_index, lower, first_index] - first[lower, first_index, second_index])
                for derivative in range(4):
                    bracket_first[derivative, lower, first_index, second_index] = (
                        second[derivative, first_index, lower, second_index]
                        + second[derivative, second_index, lower, first_index]
                        - second[derivative, lower, first_index, second_index])
    connection = np.einsum("rl,lmn->rmn", inverse, bracket) / 2
    connection_first = (np.einsum("drl,lmn->drmn", inverse_first, bracket)
                        + np.einsum("rl,dlmn->drmn", inverse, bracket_first)) / 2
    ricci = (np.einsum("aamn->mn", connection_first) - np.einsum("nama->mn", connection_first)
             + np.einsum("amn,bab->mn", connection, connection)
             - np.einsum("bma,anb->mn", connection, connection))
    return ricci - metric * np.einsum("ab,ab->", inverse, ricci) / 2


def curvature_symbol(covector, tensor):
    curvature = np.zeros((4, 4, 4, 4))
    for first, second, third, fourth in np.ndindex(curvature.shape):
        curvature[first, second, third, fourth] = (
            -covector[third] * covector[second] * tensor[first, fourth]
            -covector[fourth] * covector[first] * tensor[second, third]
            +covector[fourth] * covector[second] * tensor[first, third]
            +covector[third] * covector[first] * tensor[second, fourth]) / 2
    return curvature


@dataclass(frozen=True)
class RingMetric:
    source: GaussianSource | None = GaussianSource()
    length: float = 20.0
    center_offset: float = 4.0

    def __post_init__(self):
        if not np.isfinite(self.length) or self.length <= 0 or not np.isfinite(self.center_offset) or self.center_offset < 0:
            raise ValueError("Use a finite positive circumference and nonnegative circle offset.")

    def potential(self, positions):
        positions = np.asarray(positions, dtype=float)
        radius = self.length / (2 * math.pi)
        distances = np.sqrt(np.maximum(0., self.center_offset ** 2 + radius ** 2
                    + 2 * self.center_offset * radius * np.cos(2 * math.pi * positions / self.length)))
        if self.source is None:
            return np.zeros_like(positions)
        return np.array([self.source.radial_jet(float(distance))[0] for distance in distances.flat]).reshape(positions.shape)

    def lapse(self, positions):
        return np.sqrt(1 + 2 * self.potential(positions))

    def spatial_scale(self, positions):
        return np.sqrt(1 - 2 * self.potential(positions))

    def speed(self, positions):
        return self.lapse(positions) / self.spatial_scale(positions)

    @property
    def conservative_speed_floor(self):
        peak = 0. if self.source is None else self.source.peak_potential
        return math.sqrt((1 - 2 * peak) / (1 + 2 * peak))


def probe_unitary(sites, mass, duration, metric=RingMetric()):
    if not np.isfinite(duration) or duration < 0:
        raise ValueError("Use a finite nonnegative propagation interval.")
    hamiltonian = massive_hopping(sites, mass, metric)
    values, vectors = np.linalg.eigh(hamiltonian)
    data_gate = (vectors * np.exp(-1j * duration * values)) @ vectors.T
    rates = np.array([metric.lapse(0.), 1.2 * metric.lapse(metric.length / 2)])
    clocks = np.kron(clock_gate(rates[0], duration), clock_gate(rates[1], duration))
    return np.kron(data_gate, clocks), hamiltonian, rates


def probe_record_probabilities(state, sites):
    thresholds, denominator = postprocessing_table(sites, center_fraction=0.5)
    amplitudes = state.reshape(sites, 4)
    probabilities = np.zeros((2, 2, 2))
    for data_bit, first_clock, second_clock in np.ndindex(probabilities.shape):
        weights = thresholds / denominator if data_bit == 0 else 1 - thresholds / denominator
        effect = np.kron(clock_effect("X", first_clock), clock_effect("Y", second_clock))
        probabilities[data_bit, first_clock, second_clock] = np.einsum(
            "ja,ab,jb,j->", amplitudes.conj(), effect, amplitudes, weights).real
    if probabilities.min() < -1e-12 or abs(probabilities.sum() - 1) > 1e-12:
        raise ValueError("Retain a normalized complete data-clock record law.")
    return probabilities


def probe_schedule(sites, duration, metric=RingMetric()):
    massive_hopping(sites, 0.35, metric)
    if not np.isfinite(duration) or duration < 0:
        raise ValueError("Use a finite nonnegative read parameter.")
    delay = metric.length / sites / metric.conservative_speed_floor
    sources = [("position", site) for site in range(sites)] + [("clock", 0), ("clock", sites // 2)]
    receipt, schedule = duration, []
    for kind, source in sources:
        hops = min(source, sites - source)
        receipt += hops * delay
        schedule.append({"kind": kind, "source_site": source, "read_at": duration, "hops": hops,
                         "allocated_hop_interval": delay, "received_by_collector_at": receipt})
    return schedule


def release_probe_record(positions, clock_bits, random_integer, receipt, duration, metric=RingMetric()):
    schedule = probe_schedule(len(positions), duration, metric)
    if (any(bit not in (0, 1) for bit in positions) or sum(positions) != 1
            or len(clock_bits) != 2 or any(bit not in (0, 1) for bit in clock_bits)
            or not isinstance(random_integer, int) or not 0 <= random_integer < 4096):
        raise ValueError("Use all position and clock bits plus twelve independent fair bits.")
    if not np.isfinite(receipt) or receipt < schedule[-1]["received_by_collector_at"]:
        raise ValueError("Wait for every conservatively scheduled raw record.")
    thresholds, _ = postprocessing_table(len(positions), center_fraction=0.5)
    return int(random_integer >= thresholds[list(positions).index(1)]), *clock_bits


def nonlinear_residual_report():
    rows = []
    points = ([0., 0., 0.], [1., 0.5, 0.3], [4., 0., 0.])
    for coupling in (0.02, 0.01, 0.005):
        source = GaussianSource(gravitational_constant=coupling)
        absolute, linear_norm = [], []
        for position in points:
            metric, first, second = weak_metric_jet(position, source)
            linear = linear_einstein_from_hessian(second)
            residual = nonlinear_einstein_from_jet(metric, first, second) - linear
            absolute.append(float(np.linalg.norm(residual)))
            linear_norm.append(float(np.linalg.norm(linear)))
        rows.append({"G": coupling, "peak_potential": source.peak_potential,
                     "maximum_sampled_nonlinear_residual": max(absolute), "maximum_sampled_linear_norm": max(linear_norm)})
    return {"sample_points": points, "rows": rows,
            "scope": "residual of the declared eta+h metric against the zeroth-order source; not an error bound to a solved nonlinear spacetime"}


def run_report():
    source, metric = GaussianSource(), RingMetric()
    sites, duration = 32, 2.
    prepared = np.kron(initial_data(sites), np.ones(4) / 2)
    rows = []
    for mass in (0.35, 0.55):
        gate, hamiltonian, rates = probe_unitary(sites, mass, duration, metric)
        state = gate @ prepared
        flat_gate, _, _ = probe_unitary(sites, mass, duration, RingMetric(source=None))
        records = probe_record_probabilities(state, sites)
        flat_records = probe_record_probabilities(flat_gate @ prepared, sites)
        rows.append({"proper_mass": mass, "joint_record_probabilities": records.tolist(),
                     "flat_reference_joint_record_probabilities": flat_records.tolist(),
                     "record_TV_from_flat_reference": float(np.abs(records - flat_records).sum() / 2),
                     "data_only_probability_change_from_flat": float(records[0].sum() - flat_records[0].sum()),
                     "smooth_data_probability": float(records[0].sum()),
                     "clock_A_X_plus_probability": float(records[:, 0, :].sum()),
                     "clock_C_Y_plus_probability": float(records[:, :, 0].sum()),
                     "static_data_energy_initial": float(np.vdot(initial_data(sites), hamiltonian @ initial_data(sites)).real),
                     "static_data_energy_final": float(np.vdot(state, np.kron(hamiltonian, np.eye(4)) @ state).real),
                     "clock_rates": rates.tolist()})
    points = np.linspace(0, metric.length, 65)
    lapse, spatial, speed = metric.lapse(points), metric.spatial_scale(points), metric.speed(points)
    positions = [0., metric.length / 2]
    far, near = metric.potential(positions)
    peak = source.peak_potential
    linear_ratio = 1 + near - far
    ratio_remainder = ((near ** 2 + 2 * abs(near * far)) / (1 - 2 * peak) ** 2
                       + 3 * far ** 2 / (1 - 2 * peak) ** 2.5) / 2
    schedule = probe_schedule(sites, duration, metric)
    return {"model_version": "C9", "route_round": "09", "research_scope": "sufficiency only",
            "declared_spacetime_dimension": "3+1", "signature": "+---", "units": "hbar=c_reference=1, no SI calibration",
            "action_input": "I2=(1/(32 pi G)) integral h^{mu nu} G1_{mu nu}[h] - (1/2) integral h^{mu nu} T_{mu nu}",
            "derived_linear_equation": "G1[h]=8 pi G T; harmonic gauge gives box(trace_reverse(h))=-16 pi G T",
            "source": vars(source), "peak_potential": source.peak_potential,
            "isolated_boundary_condition": "Phi tends to zero at spatial infinity; positive Gaussian mass on R3, no periodic source subtraction",
            "ring_probe": {"length": metric.length, "center_offset": metric.center_offset,
                           "circle_radius": metric.length / (2 * math.pi), "metric_potential_at_clocks": metric.potential(positions).tolist(),
                           "clock_lapses": metric.lapse(positions).tolist(), "lapse_ratio_near_to_far": float(metric.lapse(positions)[1] / metric.lapse(positions)[0]),
                           "first_order_lapse_ratio": float(linear_ratio),
                           "chosen_metric_ratio_Taylor_remainder_bound": float(ratio_remainder),
                           "ratio_bound_scope": "square-root evaluation of the chosen eta+h metric only, not an error bound to nonlinear Einstein dynamics",
                           "lapse_min_max": [float(lapse.min()), float(lapse.max())],
                           "spatial_factor_min_max": [float(spatial.min()), float(spatial.max())],
                           "speed_min_max": [float(speed.min()), float(speed.max())],
                           "conservative_speed_floor": metric.conservative_speed_floor},
            "probe_trials": rows, "nonlinear_residual": nonlinear_residual_report(),
            "free_tensor_benchmark": {"polarizations": ["plus", "cross"], "nonzero_wave_number": 0.7,
                                      "polarization_norm_squared": 2, "canonical_amplitude_scale": "h_polarization/sqrt(16 pi G) for one normalized real spatial mode",
                                      "joint_two_oscillator_vacuum_energy": 0.7,
                                      "scope": "two selected free TT polarizations only, not an interacting quantum metric or sourced quantum-gravity completion"},
            "resources_per_fresh_trial": {"chain_qubits": sites, "clock_qubits": 2, "raw_outcome_bits": sites + 2,
                                          "coarse_joint_outcomes": 8, "independent_fair_output_bits": 12,
                                          "threshold_table_bits": sites * 13, "bit_hops": sum(event["hops"] for event in schedule),
                                          "read_parameter": duration, "allocated_hop_interval": schedule[0]["allocated_hop_interval"],
                                          "all_records_ready_by": schedule[-1]["received_by_collector_at"], "resets": 0},
            "fresh_probe_trials": 2, "cumulative_logical_qubit_initializations": 2 * (sites + 2),
            "probe_schedule": schedule,
            "limits": ["quadratic metric action, four dimensions, G, static source and universal minimal probe coupling are inputs",
                       "prescribed classical source is not computed from the unknown probe-state ensemble; no self-consistent quantum metric-matter backreaction here",
                       "source is conserved in the flat zeroth-order description, not a self-supported nonlinear static dust solution",
                       "a constrained one-dimensional ring uses the induced metric; no complete three-dimensional Dirac reduction or confining apparatus is modeled",
                       "sqrt lapse and spatial factors define exact finite probe dynamics on the chosen weak metric, not controlled higher-order GR predictions",
                       "C5 fixed-profile continuum error bounds do not automatically apply to this new source-generated profile",
                       "nonlinear residual is a diagnostic, not a distance bound to an exact Einstein solution",
                       "source preparation, support, confinement, synchronization, detectors and measurement energy remain unclosed",
                       "no lab data, sample counts, full nonlinear gravity or derivation of the adopted action from cognition"]}


class WeakMetricResponseTests(unittest.TestCase):
    def test_linear_tensor_operator_has_gauge_and_conservation_identities(self):
        generator = np.random.default_rng(9)
        for _ in range(12):
            covector = generator.integers(-3, 4, 4)
            direction = generator.integers(-3, 4, 4)
            tensor = generator.integers(-3, 4, (4, 4))
            tensor = tensor + tensor.T
            gauge = np.outer(covector, direction) + np.outer(direction, covector)
            np.testing.assert_array_equal(linear_einstein_symbol(covector, gauge), np.zeros((4, 4)))
            np.testing.assert_array_equal((ETA @ covector) @ linear_einstein_symbol(covector, tensor), np.zeros(4))

    def test_conserved_static_mode_solves_the_harmonic_source_equation(self):
        covector = np.array([0., 1., 2., -1.])
        source = np.zeros((4, 4))
        source[0, 0] = 0.3
        perturbation = harmonic_response(covector, source, 0.02)
        np.testing.assert_allclose((ETA @ covector) @ trace_reverse(perturbation), 0, atol=1e-15)
        np.testing.assert_allclose(linear_einstein_symbol(covector, perturbation), 8 * math.pi * 0.02 * source, atol=1e-15)

    def test_isolated_positive_gaussian_source_solves_poisson_without_zero_mode_removal(self):
        source = GaussianSource()
        for position in ([0., 0., 0.], [1e-5, 0., 0.], [0.8, 1.2, -0.5], [5., 0., 0.]):
            _, _, hessian = source.cartesian_jet(position)
            density = source.density(np.linalg.norm(position))
            self.assertAlmostEqual(np.trace(hessian), 4 * math.pi * source.gravitational_constant * density)
            _, _, second = weak_metric_jet(position, source)
            expected = np.zeros((4, 4))
            expected[0, 0] = 8 * math.pi * source.gravitational_constant * density
            np.testing.assert_allclose(linear_einstein_from_hessian(second), expected, atol=2e-14)

    def test_two_declared_transverse_tensor_waves_satisfy_the_vacuum_equation(self):
        covector = np.array([-1., 0., 0., 1.])
        plus = np.diag([0., 1., -1., 0.])
        cross = np.zeros((4, 4))
        cross[1, 2] = cross[2, 1] = 1
        for polarization in (plus, cross):
            np.testing.assert_array_equal((ETA @ covector) @ trace_reverse(polarization), np.zeros(4))
            np.testing.assert_array_equal(linear_einstein_symbol(covector, polarization), np.zeros((4, 4)))

    def test_quadratic_action_operator_is_self_adjoint_and_has_the_stated_variation(self):
        generator = np.random.default_rng(18)
        covector = np.array([2., 1., -1., 3.])
        first, second = (generator.normal(size=(4, 4)) for _ in range(2))
        first, second = first + first.T, second + second.T
        left = np.sum((ETA @ first @ ETA) * linear_einstein_symbol(covector, second))
        right = np.sum((ETA @ second @ ETA) * linear_einstein_symbol(covector, first))
        self.assertAlmostEqual(left, right)
        source = np.zeros((4, 4))
        source[0, 0] = 0.1
        step = 0.125
        finite_variation = (quadratic_action_symbol(covector, first + step * second, source, 0.02)
                    - quadratic_action_symbol(covector, first - step * second, source, 0.02)) / (2 * step)
        expected = np.sum((ETA @ second @ ETA) * (linear_einstein_symbol(covector, first) / (16 * math.pi * 0.02) - source / 2))
        self.assertAlmostEqual(finite_variation, expected)

    def test_tensor_tidal_curvature_is_nonzero_and_pure_gauge_has_none(self):
        covector = np.array([-1., 0., 0., 1.])
        polarization = np.diag([0., 1., -1., 0.])
        curvature = curvature_symbol(covector, polarization)
        self.assertAlmostEqual(curvature[0, 1, 0, 1], 0.5)
        self.assertAlmostEqual(curvature[0, 2, 0, 2], -0.5)
        direction = np.array([1., 2., 3., -1.])
        gauge = np.outer(covector, direction) + np.outer(direction, covector)
        np.testing.assert_array_equal(curvature_symbol(covector, gauge), np.zeros((4, 4, 4, 4)))

    def test_chosen_weak_metric_has_second_order_nonlinear_residual_not_zero(self):
        rows = nonlinear_residual_report()["rows"]
        for first, second in zip(rows, rows[1:]):
            ratio = first["maximum_sampled_nonlinear_residual"] / second["maximum_sampled_nonlinear_residual"]
            self.assertGreater(ratio, 3.8)
            self.assertLess(ratio, 4.3)

    def test_same_ring_metric_sets_both_species_speeds_and_mass_clock_scales(self):
        metric = RingMetric()
        positions = np.linspace(0, metric.length, 33)
        first_mass, second_mass = 0.35, 0.55
        np.testing.assert_allclose(first_mass * metric.lapse(positions) / first_mass,
                                   second_mass * metric.lapse(positions) / second_mass)
        np.testing.assert_allclose(metric.speed(positions), metric.lapse(positions) / metric.spatial_scale(positions))
        self.assertAlmostEqual(float(metric.lapse(0.)), float(metric.lapse(metric.length)))
        self.assertGreaterEqual(metric.speed(positions).min(), metric.conservative_speed_floor)
        self.assertLess(float(metric.lapse(metric.length / 2)), float(metric.lapse(0.)))
        first_hamiltonian = massive_hopping(32, first_mass, metric)
        second_hamiltonian = massive_hopping(32, second_mass, metric)
        np.testing.assert_allclose(first_hamiltonian - np.diag(np.diag(first_hamiltonian)),
                       second_hamiltonian - np.diag(np.diag(second_hamiltonian)))
        np.testing.assert_allclose(np.diag(first_hamiltonian) / first_mass, np.diag(second_hamiltonian) / second_mass)

    def test_probe_dynamics_preserves_norm_energy_and_common_J_representation(self):
        gate, hamiltonian, _ = probe_unitary(8, 0.35, 2.)
        initial = np.kron(initial_data(8), np.ones(4) / 2)
        final = gate @ initial
        np.testing.assert_allclose(gate.conj().T @ gate, np.eye(32), atol=3e-14)
        energy = np.kron(hamiltonian, np.eye(4))
        self.assertAlmostEqual(np.vdot(initial, energy @ initial).real, np.vdot(final, energy @ final).real)
        real_gate = real_lift(gate)
        np.testing.assert_allclose(real_gate @ encode_state(np.outer(initial, initial.conj())) @ real_gate.T,
                                   encode_state(np.outer(final, final.conj())), atol=3e-14)

    def test_fixed_source_preserves_complete_probe_preparation_mixtures(self):
        gate, _, _ = probe_unitary(8, 0.35, 1.)
        first, second = np.eye(8)[:, 0], np.eye(8)[:, 2]
        def evolved_density(data):
            state = gate @ np.kron(data, np.ones(4) / 2)
            return np.outer(state, state.conj())
        direct = (evolved_density(first) + evolved_density(second)) / 2
        other = sum(evolved_density((first + sign * second) / math.sqrt(2)) for sign in (1, -1)) / 2
        np.testing.assert_allclose(direct, other, atol=3e-14)

    def test_clock_and_data_records_keep_all_raw_bits_and_cannot_arrive_early(self):
        gate, _, rates = probe_unitary(8, 0.35, 2.)
        state = gate @ np.kron(initial_data(8), np.ones(4) / 2)
        records = probe_record_probabilities(state, 8)
        self.assertAlmostEqual(records.sum(), 1.)
        self.assertAlmostEqual(records[:, 0, :].sum(), (1 + math.cos(2 * rates[0])) / 2)
        schedule = probe_schedule(8, 2.)
        positions = [0] * 8
        positions[3] = 1
        receipt = schedule[-1]["received_by_collector_at"]
        with self.assertRaises(ValueError):
            release_probe_record(positions, [0, 1], 0, receipt - 0.01, 2.)
        count = sum(release_probe_record(positions, [0, 1], random, receipt, 2.)[0] == 0 for random in range(4096))
        thresholds, _ = postprocessing_table(8, center_fraction=0.5)
        self.assertEqual(count, thresholds[3])

    def test_free_tensor_polarizations_have_positive_canonical_oscillator_dynamics(self):
        wave_number = 0.7
        stiffness_matrix = wave_number ** 2 * np.eye(2)
        flow = harmonic_flow(stiffness_matrix, 3.)
        covariance = ground_covariance(stiffness_matrix)
        np.testing.assert_allclose(flow @ covariance @ flow.T, covariance, atol=3e-14)
        self.assertAlmostEqual(np.trace(covariance[2:, 2:] + stiffness_matrix @ covariance[:2, :2]) / 2, wave_number)

    def test_invalid_source_and_unresolved_null_response_are_rejected(self):
        with self.assertRaises(ValueError):
            GaussianSource(gravitational_constant=1.)
        with self.assertRaises(ValueError):
            harmonic_response(np.array([1., 0., 0., 1.]), np.zeros((4, 4)), 0.02)
        source = np.ones((4, 4))
        with self.assertRaises(ValueError):
            harmonic_response(np.array([0., 1., 0., 0.]), source, 0.02)

    def test_nonlinear_origin_tensor_matches_an_independent_closed_expression(self):
        source = GaussianSource()
        potential, _, curvature = source.radial_jet(0.)
        temporal, spatial = 1 + 2 * potential, 1 - 2 * potential
        expected = np.diag([6 * temporal * curvature / spatial ** 2,
                            *([2 * curvature * (1 / temporal - 1 / spatial)] * 3)])
        np.testing.assert_allclose(nonlinear_einstein_from_jet(*weak_metric_jet([0., 0., 0.], source)), expected, atol=2e-14)

    def test_positive_source_mass_and_isolated_far_field_normalization(self):
        source = GaussianSource()
        nodes, weights = np.polynomial.legendre.leggauss(96)
        radius = 6 * source.width * (nodes + 1)
        mass = 6 * source.width * np.dot(weights, 4 * math.pi * radius ** 2 * source.density(radius))
        self.assertAlmostEqual(mass, source.mass)
        far = 100 * source.width
        potential, force, _ = source.radial_jet(far)
        self.assertAlmostEqual(-potential * far, source.gravitational_constant * source.mass)
        self.assertAlmostEqual(force * far ** 2, source.gravitational_constant * source.mass)

    def test_first_order_clock_ratio_has_a_declared_convention_remainder_bound(self):
        metric = RingMetric()
        far, near = metric.potential([0., metric.length / 2])
        peak = metric.source.peak_potential
        bound = ((near ** 2 + 2 * abs(near * far)) / (1 - 2 * peak) ** 2
                 + 3 * far ** 2 / (1 - 2 * peak) ** 2.5) / 2
        ratio = float(metric.lapse(metric.length / 2) / metric.lapse(0.))
        self.assertLess(abs(ratio - (1 + near - far)), bound)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    outcome = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(WeakMetricResponseTests))
    if not outcome.wasSuccessful():
        raise SystemExit(1)
    if arguments.write_results:
        report = run_report()
        Path(__file__).with_name("weak_metric_response_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"tests": outcome.testsRun, "ring": report["ring_probe"], "probes": report["probe_trials"],
                          "nonlinear_residual": report["nonlinear_residual"], "resources": report["resources_per_fresh_trial"]}, indent=2))