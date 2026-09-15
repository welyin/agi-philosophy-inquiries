"""Compare joint quantum evolution with factorized source-background trajectories."""

import argparse
import json
import math
import platform
import unittest
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np

from relativistic_limit import IDENTITY, SIGMA_X, SIGMA_Z


SIGMA_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
PLUS_Z = np.array([1, 0], dtype=complex)
MINUS_Z = np.array([0, 1], dtype=complex)
PLUS_X = np.array([1, 1], dtype=complex) / math.sqrt(2)
MINUS_X = np.array([1, -1], dtype=complex) / math.sqrt(2)


@dataclass(frozen=True)
class Parameters:
    cutoff: int = 32
    frequency: float = 1.0
    splitting: float = 0.8
    coupling: float = 0.7
    transverse: float = 0.0

    def __post_init__(self):
        if self.cutoff < 3 or self.frequency <= 0:
            raise ValueError("Use at least three oscillator levels and a positive frequency.")
        if not all(math.isfinite(value) for value in (self.frequency, self.splitting, self.coupling, self.transverse)):
            raise ValueError("Hamiltonian parameters must be finite.")


def density(state):
    return np.outer(state, state.conj())


def expectation(state, operator):
    return float(np.trace(state @ operator).real)


def entropy_nats(state):
    eigenvalues = np.linalg.eigvalsh(state)
    if eigenvalues.min() < -1e-11:
        raise ValueError("A density matrix has a negative eigenvalue beyond rounding error.")
    positive = eigenvalues[eigenvalues > 1e-14]
    return float(-np.sum(positive * np.log(positive)))


def trace_distance(first, second):
    difference = first - second
    return float(np.sum(abs(np.linalg.eigvalsh((difference + difference.conj().T) / 2))) / 2)


def coherent_projection(amplitude, cutoff):
    state = np.empty(cutoff, dtype=complex)
    state[0] = math.exp(-abs(amplitude)**2 / 2)
    for level in range(1, cutoff):
        state[level] = state[level - 1] * amplitude / math.sqrt(level)
    return state


class JointModel:
    def __init__(self, parameters):
        self.parameters = parameters
        cutoff = parameters.cutoff
        oscillator_identity = np.eye(cutoff)
        lowering = np.diag(np.sqrt(np.arange(1, cutoff)), 1).astype(complex)
        self.number = np.diag(np.arange(cutoff, dtype=float))
        self.position = (lowering + lowering.conj().T) / math.sqrt(2)
        self.momentum = -1j * (lowering - lowering.conj().T) / math.sqrt(2)
        self.matter_energy = np.kron((parameters.splitting * SIGMA_Z + parameters.transverse * SIGMA_X) / 2,
                                    oscillator_identity)
        self.background_energy = np.kron(IDENTITY, parameters.frequency * (self.number + oscillator_identity / 2))
        self.interaction_energy = parameters.coupling * np.kron(SIGMA_Z, self.position)
        self.hamiltonian = self.matter_energy + self.background_energy + self.interaction_energy
        self.eigenvalues, self.eigenvectors = np.linalg.eigh(self.hamiltonian)

    def initial_vector(self, spinor):
        vacuum = np.zeros(self.parameters.cutoff, dtype=complex)
        vacuum[0] = 1
        return np.kron(spinor, vacuum)

    def evolve_vector(self, state, duration):
        coefficients = self.eigenvectors.conj().T @ state
        return self.eigenvectors @ (np.exp(-1j * duration * self.eigenvalues) * coefficients)

    def evolve_density(self, state, duration):
        transformed = self.eigenvectors.conj().T @ state @ self.eigenvectors
        phases = np.exp(-1j * duration * self.eigenvalues)
        return self.eigenvectors @ (phases[:, None] * transformed * phases.conj()[None, :]) @ self.eigenvectors.conj().T

    def ensemble_density(self, spinors, duration):
        return sum(density(self.evolve_vector(self.initial_vector(spinor), duration)) for spinor in spinors) / len(spinors)

    def reduced_states(self, joint_density):
        tensor = joint_density.reshape(2, self.parameters.cutoff, 2, self.parameters.cutoff)
        return np.trace(tensor, axis1=1, axis2=3), np.trace(tensor, axis1=0, axis2=2)

    def diagnostics(self, joint_density):
        matter, background = self.reduced_states(joint_density)
        mean_position = expectation(background, self.position)
        mean_momentum = expectation(background, self.momentum)
        source_mean = expectation(matter, SIGMA_Z)
        energies = {name: expectation(joint_density, operator) for name, operator in (
            ("matter", self.matter_energy), ("background", self.background_energy),
            ("interaction", self.interaction_energy), ("total", self.hamiltonian))}
        return {
            "trace": float(np.trace(joint_density).real),
            "joint_purity": float(np.trace(joint_density @ joint_density).real),
            "matter_purity": float(np.trace(matter @ matter).real),
            "background_purity": float(np.trace(background @ background).real),
            "matter_entropy_nats": entropy_nats(matter),
            "background_entropy_nats": entropy_nats(background),
            "matter_visibility": float(2 * abs(matter[0, 1])),
            "position_mean": mean_position,
            "position_variance": expectation(background, self.position @ self.position) - mean_position**2,
            "momentum_mean": mean_momentum,
            "momentum_variance": expectation(background, self.momentum @ self.momentum) - mean_momentum**2,
            "source_position_covariance": expectation(joint_density, np.kron(SIGMA_Z, self.position)) - source_mean * mean_position,
            "mean_oscillator_number": expectation(background, self.number),
            "highest_level_population": float(background[-1, -1].real),
            "top_four_levels_population": float(np.diag(background).real[-4:].sum()),
            "energies": energies,
        }


def analytic_displacement(parameters, duration):
    if parameters.transverse != 0:
        raise ValueError("The displaced-oscillator solution here requires a conserved sigma_z source.")
    return parameters.coupling / (math.sqrt(2) * parameters.frequency) * (np.exp(-1j * parameters.frequency * duration) - 1)


def analytic_joint_vector(parameters, spinor, duration):
    amplitude = analytic_displacement(parameters, duration)
    global_phase = np.exp(-1j * parameters.frequency * duration / 2 + 1j * parameters.coupling**2
                          * (parameters.frequency * duration - math.sin(parameters.frequency * duration))
                          / (2 * parameters.frequency**2))
    first = spinor[0] * np.exp(-1j * parameters.splitting * duration / 2) * coherent_projection(amplitude, parameters.cutoff)
    second = spinor[1] * np.exp(1j * parameters.splitting * duration / 2) * coherent_projection(-amplitude, parameters.cutoff)
    return global_phase * np.concatenate((first, second))


def analytic_observables(parameters, duration):
    amplitude = analytic_displacement(parameters, duration)
    displacement_position = math.sqrt(2) * amplitude.real
    displacement_momentum = math.sqrt(2) * amplitude.imag
    visibility = math.exp(-2 * abs(amplitude)**2)
    return {
        "conditional_position_plus_z": displacement_position,
        "conditional_momentum_plus_z": displacement_momentum,
        "balanced_source_position_variance": 0.5 + displacement_position**2,
        "balanced_source_momentum_variance": 0.5 + displacement_momentum**2,
        "coherent_spin_visibility": visibility,
        "coherent_spin_purity": (1 + visibility**2) / 2,
        "mean_oscillator_number": abs(amplitude)**2,
        "background_energy_gain": parameters.frequency * abs(amplitude)**2,
        "interaction_energy": parameters.coupling * displacement_position,
    }


def factorized_joint_density(parameters, spinor, duration):
    amplitude = analytic_displacement(parameters, duration)
    source_mean = float(np.vdot(spinor, SIGMA_Z @ spinor).real)
    force_phase = source_mean * parameters.coupling**2 / parameters.frequency**2 * (
        math.sin(parameters.frequency * duration) - parameters.frequency * duration)
    angle = parameters.splitting * duration / 2 + force_phase
    evolved_spinor = spinor * np.array([np.exp(-1j * angle), np.exp(1j * angle)])
    background = coherent_projection(source_mean * amplitude, parameters.cutoff)
    return density(np.kron(evolved_spinor, background))


def preparation_comparison(parameters):
    model = JointModel(parameters)
    duration = math.pi / parameters.frequency
    eigen_ensemble = model.ensemble_density((PLUS_Z, MINUS_Z), duration)
    coherent_ensemble = model.ensemble_density((PLUS_X, MINUS_X), duration)
    direct = model.evolve_density(model.ensemble_density((PLUS_Z, MINUS_Z), 0), duration)
    factorized_eigen = sum(factorized_joint_density(parameters, state, duration) for state in (PLUS_Z, MINUS_Z)) / 2
    factorized_coherent = sum(factorized_joint_density(parameters, state, duration) for state in (PLUS_X, MINUS_X)) / 2
    _, exact_background = model.reduced_states(eigen_ensemble)
    _, factorized_background = model.reduced_states(factorized_coherent)
    coherent_pure = density(model.evolve_vector(model.initial_vector(PLUS_X), duration))
    coherent_product = factorized_joint_density(parameters, PLUS_X, duration)
    return {
        "duration": duration,
        "initial_preparation": "I_spin/2 tensor oscillator vacuum, decomposed into equal z eigenstates or equal x eigenstates.",
        "exact_ensemble_joint_trace_distance": trace_distance(eigen_ensemble, coherent_ensemble),
        "direct_mixed_vs_ensemble_trace_distance": trace_distance(direct, eigen_ensemble),
        "exact_mixed_output": model.diagnostics(eigen_ensemble),
        "factorized_z_ensemble": model.diagnostics(factorized_eigen),
        "factorized_x_ensemble": model.diagnostics(factorized_coherent),
        "factorized_ensemble_joint_trace_distance": trace_distance(factorized_eigen, factorized_coherent),
        "exact_vs_factorized_x_background_trace_distance": trace_distance(exact_background, factorized_background),
        "coherent_pure_input_exact": model.diagnostics(coherent_pure),
        "coherent_pure_input_factorized": model.diagnostics(coherent_product),
        "analytic_reference": analytic_observables(parameters, duration),
        "entanglement_caveat": "The mixed I/2 input evolves into a separable mixture with classical source-background correlations. Entanglement entropy is claimed only for the separate pure +x input.",
        "factorization_definition": "A product-state mean field retaining oscillator vacuum variance 1/2, not an approximation that sets all background uncertainty to zero.",
    }


def cutoff_study():
    cases = []
    for coupling in (0.7, 1.4):
        rows = []
        for cutoff in (6, 8, 12, 20, 32):
            parameters = Parameters(cutoff=cutoff, coupling=coupling)
            model = JointModel(parameters)
            errors = []
            variance_errors = []
            omitted_norms = []
            boundary_populations = []
            for duration in np.linspace(0, 2 * math.pi, 25):
                exact_projection = analytic_joint_vector(parameters, PLUS_X, float(duration))
                actual = model.evolve_vector(model.initial_vector(PLUS_X), float(duration))
                errors.append(float(np.linalg.norm(actual - exact_projection)))
                omitted_norms.append(max(0.0, 1 - float(np.vdot(exact_projection, exact_projection).real)))
                _, background = model.reduced_states(density(actual))
                position_mean = expectation(background, model.position)
                position_variance = expectation(background, model.position @ model.position) - position_mean**2
                variance_errors.append(abs(position_variance - analytic_observables(parameters, float(duration))["balanced_source_position_variance"]))
                boundary_populations.append(float(background[-1, -1].real))
            rows.append({
                "cutoff": cutoff,
                "max_sampled_vector_error_against_unnormalized_analytic_projection": max(errors),
                "max_sampled_position_variance_error": max(variance_errors),
                "max_sampled_analytic_omitted_norm_squared": max(omitted_norms),
                "max_sampled_highest_level_population": max(boundary_populations),
            })
        cases.append({"coupling": coupling, "rows": rows})
    return {
        "time_interval": [0.0, 2 * math.pi],
        "time_samples": 25,
        "cases": cases,
        "interpretation": "Spectral evolution with an oscillator cutoff is compared with the independently displaced infinite-oscillator state projected without renormalization. Small boundary occupation alone is not used as a proof of accuracy.",
    }


def approximation_study():
    rows = []
    for coupling in (0.0, 0.05, 0.1, 0.2, 0.4, 0.7):
        parameters = Parameters(coupling=coupling)
        model = JointModel(parameters)
        duration = math.pi / parameters.frequency
        actual = density(model.evolve_vector(model.initial_vector(PLUS_X), duration))
        product = factorized_joint_density(parameters, PLUS_X, duration)
        actual_spin, actual_background = model.reduced_states(actual)
        product_spin, product_background = model.reduced_states(product)
        amplitude = analytic_displacement(parameters, duration)
        rows.append({
            "coupling_over_frequency": coupling / parameters.frequency,
            "branch_displacement_squared": abs(amplitude)**2,
            "matter_trace_distance": trace_distance(actual_spin, product_spin),
            "matter_trace_distance_analytic": -math.expm1(-2 * abs(amplitude)**2) / 2,
            "background_trace_distance": trace_distance(actual_background, product_background),
            "joint_trace_distance": trace_distance(actual, product),
            "joint_trace_distance_analytic": math.sqrt(-math.expm1(-abs(amplitude)**2)),
            "excess_position_variance_missed": analytic_observables(parameters, duration)["balanced_source_position_variance"] - 0.5,
        })
    error_budget = 0.01
    return {
        "initial_state": "Pure +x source with oscillator vacuum; factorized comparison uses the same vacuum fluctuations.",
        "duration": math.pi,
        "rows": rows,
        "all_time_analytic_error_budget_for_this_state": {
            "trace_distance_budget": error_budget,
            "coupling_ratio_bound_for_matter_only": math.sqrt(-math.log1p(-2 * error_budget)) / 2,
            "coupling_ratio_bound_for_joint_state": math.sqrt(-math.log1p(-error_budget**2) / 2),
        },
        "caveat": "Bounds are specific to this longitudinal single-mode model, balanced pure source, and stated trace-distance observable. They do not establish a general mean-field criterion or classical-gravity limit.",
    }


def joint_time_trace(parameters):
    model = JointModel(parameters)
    rows = []
    for duration in (0, math.pi / 2, math.pi, 3 * math.pi / 2, 2 * math.pi):
        state = density(model.evolve_vector(model.initial_vector(PLUS_X), duration))
        diagnostics = model.diagnostics(state)
        rows.append({"duration": duration, **{key: diagnostics[key] for key in (
            "matter_visibility", "matter_entropy_nats", "position_variance", "momentum_variance",
            "source_position_covariance", "mean_oscillator_number", "energies")}})
    return {
        "rows": rows,
        "initial_background_position_variance": 0.5,
        "probability_interpretation": "Loss of reduced visibility is reversible correlation with one oscillator, not outcome selection or an irreversible collapse law.",
    }


def transverse_control():
    parameters = Parameters(transverse=0.6)
    model = JointModel(parameters)
    duration = 4.0
    exact_z = model.ensemble_density((PLUS_Z, MINUS_Z), duration)
    exact_x = model.ensemble_density((PLUS_X, MINUS_X), duration)
    higher_model = JointModel(replace(parameters, cutoff=48))
    higher = higher_model.evolve_vector(higher_model.initial_vector(PLUS_X), duration).reshape(2, 48)
    actual = model.evolve_vector(model.initial_vector(PLUS_X), duration)
    projected = higher[:, :parameters.cutoff].ravel()
    initial = model.diagnostics(density(model.initial_vector(PLUS_X)))
    final = model.diagnostics(density(actual))
    return {
        "parameters": asdict(parameters),
        "duration": duration,
        "exact_ensemble_joint_trace_distance": trace_distance(exact_z, exact_x),
        "source_commutator_norm": float(np.linalg.norm(model.hamiltonian @ np.kron(SIGMA_Z, np.eye(parameters.cutoff))
                                                       - np.kron(SIGMA_Z, np.eye(parameters.cutoff)) @ model.hamiltonian, ord=2)),
        "initial_pure_state_energy": initial["energies"],
        "final_pure_state_energy": final["energies"],
        "final_pure_state_matter_purity": final["matter_purity"],
        "final_source_position_covariance": final["source_position_covariance"],
        "cutoff_32_vs_48_state_l2_difference": float(np.linalg.norm(actual - projected)),
        "higher_cutoff_omitted_norm_squared": max(0.0, 1 - float(np.vdot(projected, projected).real)),
        "scope": "This noncommuting control has no displaced-oscillator analytic solution used here. It tests exact finite-dimensional linear evolution and cutoff agreement, not a general interacting-field theorem; no mean-field performance is claimed for it.",
    }


class QuantumBackgroundTests(unittest.TestCase):
    def setUp(self):
        self.parameters = Parameters()
        self.model = JointModel(self.parameters)

    def test_hamiltonian_is_hermitian_and_source_is_conserved(self):
        np.testing.assert_allclose(self.model.hamiltonian, self.model.hamiltonian.conj().T, atol=1e-13)
        source = np.kron(SIGMA_Z, np.eye(self.parameters.cutoff))
        np.testing.assert_allclose(self.model.hamiltonian @ source, source @ self.model.hamiltonian, atol=1e-13)

    def test_cutoff_commutator_has_explicit_boundary_term(self):
        commutator = self.model.position @ self.model.momentum - self.model.momentum @ self.model.position
        expected = np.eye(self.parameters.cutoff, dtype=complex)
        expected[-1, -1] -= self.parameters.cutoff
        np.testing.assert_allclose(commutator, 1j * expected, atol=2e-14)

    def test_analytic_conditional_vectors_match_spectral_evolution(self):
        for spinor in (PLUS_Z, MINUS_Z, PLUS_X, MINUS_X):
            for duration in (0, 0.3, 1.0, math.pi, 2 * math.pi):
                actual = self.model.evolve_vector(self.model.initial_vector(spinor), duration)
                expected = analytic_joint_vector(self.parameters, spinor, duration)
                np.testing.assert_allclose(actual, expected, atol=1e-11)

    def test_equivalent_preparations_produce_identical_joint_states(self):
        for duration in (0, 0.7, math.pi, 5.1):
            first = self.model.ensemble_density((PLUS_Z, MINUS_Z), duration)
            second = self.model.ensemble_density((PLUS_X, MINUS_X), duration)
            self.assertLess(trace_distance(first, second), 1e-12)
            direct = self.model.evolve_density(self.model.ensemble_density((PLUS_Z, MINUS_Z), 0), duration)
            self.assertLess(trace_distance(first, direct), 1e-12)

    def test_background_variance_and_coherence_match_analytic_predictions(self):
        for duration in (0, 0.5, math.pi, 2 * math.pi):
            evolved = density(self.model.evolve_vector(self.model.initial_vector(PLUS_X), duration))
            actual = self.model.diagnostics(evolved)
            expected = analytic_observables(self.parameters, duration)
            self.assertAlmostEqual(actual["position_variance"], expected["balanced_source_position_variance"], places=10)
            self.assertAlmostEqual(actual["matter_visibility"], expected["coherent_spin_visibility"], places=10)
            self.assertAlmostEqual(actual["matter_purity"], expected["coherent_spin_purity"], places=10)

    def test_energy_exchange_includes_interaction_energy(self):
        initial = density(self.model.initial_vector(PLUS_X))
        initial_energy = expectation(initial, self.model.hamiltonian)
        for duration in (0.4, math.pi, 4.0):
            evolved = self.model.evolve_density(initial, duration)
            diagnostics = self.model.diagnostics(evolved)
            self.assertAlmostEqual(diagnostics["energies"]["total"], initial_energy, places=11)
            self.assertAlmostEqual(diagnostics["energies"]["background"] - self.parameters.frequency / 2,
                                   -diagnostics["energies"]["interaction"], places=11)
            self.assertAlmostEqual(diagnostics["energies"]["matter"], 0, places=11)

    def test_pure_entanglement_and_mixed_classical_correlations_are_distinct(self):
        duration = math.pi
        pure = density(self.model.evolve_vector(self.model.initial_vector(PLUS_X), duration))
        mixed = self.model.ensemble_density((PLUS_Z, MINUS_Z), duration)
        pure_spin, pure_background = self.model.reduced_states(pure)
        _, mixed_background = self.model.reduced_states(mixed)
        self.assertAlmostEqual(entropy_nats(pure), 0, places=11)
        self.assertAlmostEqual(entropy_nats(pure_spin), entropy_nats(pure_background), places=11)
        self.assertGreater(entropy_nats(pure_spin), 0.6)
        self.assertAlmostEqual(entropy_nats(mixed), math.log(2), places=11)
        self.assertLess(trace_distance(pure_background, mixed_background), 1e-12)

    def test_factorized_branch_extension_retains_preparation_dependence(self):
        comparison = preparation_comparison(self.parameters)
        self.assertAlmostEqual(comparison["factorized_z_ensemble"]["position_variance"], 2.46, places=10)
        self.assertAlmostEqual(comparison["factorized_x_ensemble"]["position_variance"], 0.5, places=10)
        self.assertAlmostEqual(comparison["coherent_pure_input_factorized"]["matter_visibility"], 1, places=10)
        self.assertLess(comparison["coherent_pure_input_exact"]["matter_visibility"], 0.15)

    def test_source_eigenstate_product_ansatz_is_exact_here(self):
        for spinor in (PLUS_Z, MINUS_Z):
            for duration in (0.3, math.pi):
                exact = density(self.model.evolve_vector(self.model.initial_vector(spinor), duration))
                self.assertLess(trace_distance(exact, factorized_joint_density(self.parameters, spinor, duration)), 1e-11)

    def test_background_returns_and_visibility_revives(self):
        duration = 2 * math.pi / self.parameters.frequency
        evolved = density(self.model.evolve_vector(self.model.initial_vector(PLUS_X), duration))
        result = self.model.diagnostics(evolved)
        self.assertAlmostEqual(result["matter_visibility"], 1, places=10)
        self.assertAlmostEqual(result["mean_oscillator_number"], 0, places=10)
        self.assertAlmostEqual(result["position_variance"], 0.5, places=10)

    def test_low_lying_spectrum_matches_displaced_oscillator_levels(self):
        expected = sorted(self.parameters.frequency * (level + 0.5)
                          - self.parameters.coupling**2 / (2 * self.parameters.frequency)
                          + sign * self.parameters.splitting / 2
                          for level in range(4) for sign in (-1, 1))
        np.testing.assert_allclose(self.model.eigenvalues[:6], expected[:6], atol=1e-11)

    def test_reduced_channel_is_linear_positive_and_trace_preserving(self):
        generator = np.random.default_rng(606)
        factor = generator.normal(size=(2, 2)) + 1j * generator.normal(size=(2, 2))
        spin = factor @ factor.conj().T
        spin /= np.trace(spin)
        vacuum = density(self.model.initial_vector(PLUS_Z)[:self.parameters.cutoff])
        evolved = self.model.evolve_density(np.kron(spin, vacuum), 1.2)
        reduced, _ = self.model.reduced_states(evolved)
        expected = spin.copy()
        multiplier = analytic_observables(self.parameters, 1.2)["coherent_spin_visibility"]
        expected[0, 1] *= multiplier * np.exp(-1j * self.parameters.splitting * 1.2)
        expected[1, 0] = expected[0, 1].conj()
        np.testing.assert_allclose(reduced, expected, atol=1e-11)
        self.assertAlmostEqual(float(np.trace(reduced).real), 1, places=11)
        self.assertGreaterEqual(float(np.linalg.eigvalsh(reduced).min()), -1e-12)

    def test_cutoff_refinement_matches_analytic_infinite_oscillator(self):
        coarse_parameters = Parameters(cutoff=8, coupling=1.4)
        fine_parameters = replace(coarse_parameters, cutoff=32)
        errors = []
        for parameters in (coarse_parameters, fine_parameters):
            model = JointModel(parameters)
            actual = model.evolve_vector(model.initial_vector(PLUS_X), math.pi)
            expected = analytic_joint_vector(parameters, PLUS_X, math.pi)
            errors.append(float(np.linalg.norm(actual - expected)))
        self.assertGreater(errors[0], 1e-3)
        self.assertLess(errors[1], 1e-7)
        self.assertLess(errors[1], errors[0] / 1000)

    def test_weak_coupling_error_formulas_and_zero_coupling_control(self):
        study = approximation_study()
        for row in study["rows"]:
            self.assertAlmostEqual(row["matter_trace_distance"], row["matter_trace_distance_analytic"], places=10)
            self.assertAlmostEqual(row["joint_trace_distance"], row["joint_trace_distance_analytic"], places=10)
        self.assertLess(study["rows"][0]["joint_trace_distance"], 1e-12)
        self.assertLess(study["rows"][1]["matter_trace_distance"], 0.01)
        self.assertGreater(study["rows"][2]["matter_trace_distance"], 0.01)
        bounds = study["all_time_analytic_error_budget_for_this_state"]
        self.assertLess(bounds["coupling_ratio_bound_for_joint_state"], bounds["coupling_ratio_bound_for_matter_only"])

    def test_joint_energy_has_a_cutoff_independent_lower_bound(self):
        for transverse in (0.0, 0.6):
            for cutoff in (8, 16, 32):
                parameters = Parameters(cutoff=cutoff, transverse=transverse)
                model = JointModel(parameters)
                lower_bound = parameters.frequency / 2 - math.hypot(parameters.splitting, transverse) / 2
                lower_bound -= parameters.coupling**2 / (2 * parameters.frequency)
                self.assertGreaterEqual(float(model.eigenvalues.min()), lower_bound - 1e-12)

    def test_noncommuting_source_control_preserves_preparation_equivalence(self):
        results = transverse_control()
        self.assertGreater(results["source_commutator_norm"], 0.5)
        self.assertLess(results["exact_ensemble_joint_trace_distance"], 1e-12)
        self.assertLess(results["cutoff_32_vs_48_state_l2_difference"], 1e-9)
        self.assertAlmostEqual(results["initial_pure_state_energy"]["total"], results["final_pure_state_energy"]["total"], places=11)
        self.assertGreater(abs(results["initial_pure_state_energy"]["matter"] - results["final_pure_state_energy"]["matter"]), 0.01)

    def test_coherence_can_revive_without_markov_semigroup(self):
        full = analytic_observables(self.parameters, 2 * math.pi)["coherent_spin_visibility"]
        independent_halves = analytic_observables(self.parameters, math.pi)["coherent_spin_visibility"]**2
        self.assertAlmostEqual(full, 1, places=12)
        self.assertLess(independent_halves, 0.03)

    def test_remote_marginal_unchanged_by_reduced_local_channel(self):
        duration = math.pi
        visibility = analytic_observables(self.parameters, duration)["coherent_spin_visibility"]
        rotation = np.diag([np.exp(-1j * self.parameters.splitting * duration / 2),
                            np.exp(1j * self.parameters.splitting * duration / 2)])
        kraus = (math.sqrt((1 + visibility) / 2) * rotation,
                 math.sqrt((1 - visibility) / 2) * rotation @ SIGMA_Z)
        bell = density(np.array([1, 0, 0, 1], dtype=complex) / math.sqrt(2))
        evolved = sum(np.kron(operator, IDENTITY) @ bell @ np.kron(operator.conj().T, IDENTITY) for operator in kraus)
        remote = np.trace(evolved.reshape(2, 2, 2, 2), axis1=0, axis2=2)
        np.testing.assert_allclose(remote, IDENTITY / 2, atol=1e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(QuantumBackgroundTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    parameters = Parameters()
    results = {
        "scope": "A standard longitudinal two-level system plus one quantum oscillator; not a quantization of the full earlier field model or a theory of spacetime.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "parameters": asdict(parameters),
        "units": "hbar=1; Q=(a+a_dagger)/sqrt(2), P=(a-a_dagger)/(i*sqrt(2)); oscillator energy=omega*(n+1/2).",
        "criteria_before_first_run": [
            "Joint linear unitary evolution must agree for different ensembles of the same initial joint density matrix.",
            "Numerical vectors and oscillator moments must match an independently derived conditional-displacement solution.",
            "The comparison retains the same initial oscillator vacuum uncertainty in exact and factorized models.",
            "Entanglement claims refer only to initially pure joint states; separable correlated mixed outputs are distinguished.",
        ],
        "preparation_comparison": preparation_comparison(parameters),
        "joint_time_trace": joint_time_trace(parameters),
        "cutoff_study": cutoff_study(),
        "approximation_study": approximation_study(),
        "noncommuting_source_control": transverse_control(),
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    encoded = json.dumps(results, indent=2)
    if arguments.write_results:
        output = Path(__file__).resolve().parent / "quantum_background_results.json"
        output.write_text(encoded + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    comparison = results["preparation_comparison"]
    compact = {
        "automated_checks": results["automated_checks"],
        "exact_ensemble_joint_trace_distance": comparison["exact_ensemble_joint_trace_distance"],
        "position_variances": {name: comparison[name]["position_variance"] for name in (
            "exact_mixed_output", "factorized_z_ensemble", "factorized_x_ensemble")},
        "pure_source_visibility": {"exact": comparison["coherent_pure_input_exact"]["matter_visibility"],
                                   "factorized": comparison["coherent_pure_input_factorized"]["matter_visibility"]},
        "joint_time_trace": results["joint_time_trace"],
        "cutoff_study": results["cutoff_study"],
        "approximation_study": results["approximation_study"],
        "noncommuting_source_control": results["noncommuting_source_control"],
        "scope": results["scope"],
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()