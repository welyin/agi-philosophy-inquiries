"""Reproduce the finite-dimensional examples in the hypothesis manuscript."""

import json
import math
import unittest

import numpy as np


TOLERANCE = 1e-12


def entropy(state):
    eigenvalues = np.linalg.eigvalsh(state)
    if eigenvalues.min() < -TOLERANCE:
        raise ValueError("The state must be positive semidefinite.")
    positive = eigenvalues[eigenvalues > 0]
    return float(-np.sum(positive * np.log(positive)))


def kraus_operators(strength):
    if strength < 0:
        raise ValueError("The dimensionless strength must be nonnegative.")
    success = np.diag([1.0, math.exp(-strength / 2)])
    failure = np.diag([0.0, math.sqrt(-math.expm1(-strength))])
    return success, failure


def measurement_branches(state, strength):
    branches = []
    for operator in kraus_operators(strength):
        unnormalized = operator @ state @ operator.conj().T
        probability = float(np.trace(unnormalized).real)
        conditional = unnormalized / probability if probability > 0 else None
        branches.append((probability, conditional))
    return branches


def channel(state, strength):
    return sum(
        operator @ state @ operator.conj().T
        for operator in kraus_operators(strength)
    )


def local_channel(joint_state, strength):
    operators = [np.kron(operator, np.eye(2)) for operator in kraus_operators(strength)]
    return sum(operator @ joint_state @ operator.conj().T for operator in operators)


def remote_marginal(joint_state):
    return np.trace(joint_state.reshape(2, 2, 2, 2), axis1=0, axis2=2)


def exact_evolution(state, rate, duration, frequency):
    if rate < 0 or duration < 0:
        raise ValueError("The rate and duration must be nonnegative.")
    unitary = np.diag([1.0, np.exp(-1j * frequency * duration)])
    return unitary @ channel(state, 2 * rate * duration) @ unitary.conj().T


def random_density(generator, dimension):
    factor = generator.normal(size=(dimension, dimension)) + 1j * generator.normal(
        size=(dimension, dimension)
    )
    state = factor @ factor.conj().T
    return state / np.trace(state)


def reproducible_results():
    strength = math.log(9)
    initial = np.diag([0.1, 0.9])
    probability, filtered = measurement_branches(initial, strength)[0]
    old_derivative = filtered - initial
    plus = np.full((2, 2), 0.5)
    nonselective = channel(plus, strength)
    bell_vector = np.array([1, 0, 0, 1], dtype=complex) / math.sqrt(2)
    bell_state = np.outer(bell_vector, bell_vector.conj())
    success = np.kron(kraus_operators(strength)[0], np.eye(2))
    selected_bell = success @ bell_state @ success.conj().T
    selected_bell /= np.trace(selected_bell)
    return {
        "interpretation": "Analytic toy-model checks; no experimental data or gravity derivation.",
        "strength_beta_gap": strength,
        "counterexample_initial_populations": np.diag(initial).tolist(),
        "counterexample_success_probability": probability,
        "counterexample_filtered_populations": np.diag(filtered).tolist(),
        "counterexample_initial_entropy_over_kB": entropy(initial),
        "counterexample_filtered_entropy_over_kB": entropy(filtered),
        "old_entropy_rate_over_kB_kappa": float(
            -np.sum(np.diag(old_derivative) * np.log(np.diag(initial)))
        ),
        "nonselective_plus_populations": np.diag(nonselective).tolist(),
        "nonselective_plus_eigenvalues": np.linalg.eigvalsh(nonselective).tolist(),
        "nonselective_plus_entropy_over_kB": entropy(nonselective),
        "nonselective_plus_visibility": float(2 * abs(nonselective[0, 1])),
        "bell_remote_unconditional_populations": np.diag(
            remote_marginal(local_channel(bell_state, strength))
        ).real.tolist(),
        "bell_remote_postselected_populations": np.diag(
            remote_marginal(selected_bell)
        ).real.tolist(),
        "illustrative_rate_bound_not_data": {
            "visibility_ratio_lower_bound": 0.99,
            "duration_seconds": 1.0,
            "activity_difference_dimensionless": 0.5,
            "alpha_upper_bound_per_second": -math.log(0.99) / 0.5,
        },
    }


class HypothesisModelTests(unittest.TestCase):
    def assert_density(self, state):
        np.testing.assert_allclose(state, state.conj().T, atol=TOLERANCE)
        self.assertAlmostEqual(float(np.trace(state).real), 1.0, places=12)
        self.assertGreaterEqual(float(np.linalg.eigvalsh(state).min()), -TOLERANCE)

    def test_entropy_counterexample_and_old_continuous_equation(self):
        results = reproducible_results()
        self.assertAlmostEqual(results["counterexample_initial_entropy_over_kB"], 0.3250829733914482)
        self.assertAlmostEqual(results["counterexample_filtered_entropy_over_kB"], math.log(2))
        self.assertGreater(results["old_entropy_rate_over_kB_kappa"], 0)

    def test_completeness_including_endpoints(self):
        for strength in (0.0, 1e-10, math.log(9), 10.0, 1000.0):
            operators = kraus_operators(strength)
            completeness = sum(operator.conj().T @ operator for operator in operators)
            np.testing.assert_allclose(completeness, np.eye(2), atol=TOLERANCE)

    def test_nonselective_map_preserves_density_and_populations(self):
        generator = np.random.default_rng(20260914)
        for _ in range(32):
            state = random_density(generator, 2)
            updated = channel(state, math.log(9))
            self.assert_density(updated)
            np.testing.assert_allclose(np.diag(updated), np.diag(state), atol=TOLERANCE)

    def test_branch_average_equals_channel(self):
        generator = np.random.default_rng(8)
        for _ in range(16):
            state = random_density(generator, 2)
            branches = measurement_branches(state, math.log(9))
            self.assertAlmostEqual(sum(probability for probability, _ in branches), 1)
            averaged = np.zeros((2, 2), dtype=complex)
            for probability, conditional in branches:
                if conditional is not None:
                    self.assert_density(conditional)
                    averaged += probability * conditional
            np.testing.assert_allclose(averaged, channel(state, math.log(9)), atol=TOLERANCE)

    def test_zero_probability_branch_is_not_normalized(self):
        branches = measurement_branches(np.diag([1.0, 0.0]), math.log(9))
        self.assertEqual(branches[1][0], 0)
        self.assertIsNone(branches[1][1])

    def test_no_signalling_on_entangled_and_random_inputs(self):
        generator = np.random.default_rng(42)
        bell_vector = np.array([1, 0, 0, 1], dtype=complex) / math.sqrt(2)
        states = [np.outer(bell_vector, bell_vector.conj())]
        states.extend(random_density(generator, 4) for _ in range(32))
        for state in states:
            updated = local_channel(state, math.log(9))
            self.assert_density(updated)
            np.testing.assert_allclose(remote_marginal(updated), remote_marginal(state), atol=TOLERANCE)

    def test_postselection_changes_remote_conditional_state_only(self):
        results = reproducible_results()
        np.testing.assert_allclose(results["bell_remote_unconditional_populations"], [0.5, 0.5])
        np.testing.assert_allclose(results["bell_remote_postselected_populations"], [0.9, 0.1])

    def test_coherence_and_entropy_of_plus_state(self):
        results = reproducible_results()
        np.testing.assert_allclose(results["nonselective_plus_eigenvalues"], [1 / 3, 2 / 3])
        self.assertAlmostEqual(results["nonselective_plus_visibility"], 1 / 3)
        self.assertAlmostEqual(results["nonselective_plus_entropy_over_kB"], 0.6365141682948128)

    def test_semigroup_and_hamiltonian_phase(self):
        initial = np.full((2, 2), 0.5)
        rate, frequency = 0.4, 1.7
        direct = exact_evolution(initial, rate, 0.8, frequency)
        two_steps = exact_evolution(exact_evolution(initial, rate, 0.3, frequency), rate, 0.5, frequency)
        np.testing.assert_allclose(direct, two_steps, atol=TOLERANCE)
        expected_coherence = 0.5 * np.exp((-rate + 1j * frequency) * 0.8)
        self.assertAlmostEqual(abs(direct[0, 1] - expected_coherence), 0, places=12)

    def test_generator_and_energy_balance(self):
        state = np.array([[0.4, 0.2j], [-0.2j, 0.6]])
        rate, frequency, step = 0.4, 1.7, 1e-7
        hamiltonian_over_hbar = np.diag([0.0, frequency])
        sigma_z = np.diag([1.0, -1.0])
        derivative = -1j * (hamiltonian_over_hbar @ state - state @ hamiltonian_over_hbar)
        derivative += (rate / 2) * (sigma_z @ state @ sigma_z - state)
        numeric = (exact_evolution(state, rate, step, frequency) - state) / step
        np.testing.assert_allclose(numeric, derivative, atol=1e-7)
        self.assertAlmostEqual(abs(np.trace(hamiltonian_over_hbar @ derivative)), 0, places=12)

    def test_positive_metric_pullback_samples(self):
        generator = np.random.default_rng(95)
        for _ in range(16):
            factor = generator.normal(size=(6, 6))
            information_metric = factor.T @ factor
            projection = generator.normal(size=(6, 4))
            pulled_back = projection.T @ information_metric @ projection
            self.assertGreaterEqual(float(np.linalg.eigvalsh(pulled_back).min()), -1e-10)

    def test_illustrative_exclusion_bound(self):
        bound = reproducible_results()["illustrative_rate_bound_not_data"]["alpha_upper_bound_per_second"]
        self.assertAlmostEqual(bound, 0.0201006717070029)
        self.assertAlmostEqual(math.exp(-bound * 0.5), 0.99)


if __name__ == "__main__":
    print(json.dumps(reproducible_results(), indent=2), flush=True)
    unittest.main(verbosity=2)