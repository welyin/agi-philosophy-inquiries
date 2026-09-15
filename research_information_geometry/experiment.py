"""Test state-only reconstruction and phase-response learning on four qubits."""

import argparse
import json
import math
import platform
import unittest
from collections import deque
from pathlib import Path

import numpy as np


SITE_COUNT = 4
CHAIN_EDGES = ((0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0))
STAR_EDGES = ((0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0))
TOLERANCE = 1e-11


def graph_laplacian(site_count, edges):
    hamiltonian = np.zeros((site_count, site_count), dtype=float)
    seen_edges = set()
    for source, target, weight in edges:
        if not (0 <= source < site_count and 0 <= target < site_count):
            raise ValueError("An endpoint is outside the site set.")
        if source == target or not np.isfinite(weight) or weight <= 0:
            raise ValueError("Use distinct endpoints and positive finite weights.")
        edge_key = tuple(sorted((source, target)))
        if edge_key in seen_edges:
            raise ValueError("Duplicate undirected edge.")
        seen_edges.add(edge_key)
        hamiltonian[source, source] += weight
        hamiltonian[target, target] += weight
        hamiltonian[source, target] -= weight
        hamiltonian[target, source] -= weight
    return hamiltonian


def uniform_one_excitation_state(site_count):
    state = np.zeros(2**site_count, dtype=complex)
    for site in range(site_count):
        state[1 << (site_count - site - 1)] = 1 / math.sqrt(site_count)
    return state


def local_operator(operator, site, site_count):
    embedded = np.ones((1, 1), dtype=complex)
    for current in range(site_count):
        embedded = np.kron(embedded, operator if current == site else np.eye(2))
    return embedded


def full_spin_hamiltonian(site_count, edges):
    raising = np.array([[0, 0], [1, 0]], dtype=complex)
    number = np.diag([0.0, 1.0])
    raising_operators = [local_operator(raising, site, site_count) for site in range(site_count)]
    number_operators = [local_operator(number, site, site_count) for site in range(site_count)]
    hamiltonian = np.zeros((2**site_count, 2**site_count), dtype=complex)
    for source, target, weight in edges:
        hopping = raising_operators[source] @ raising_operators[target].conj().T
        hamiltonian += weight * (
            number_operators[source] + number_operators[target] - hopping - hopping.conj().T
        )
    return hamiltonian


def observable_covariance(state_vector, operators):
    centered = []
    for operator in operators:
        applied = operator @ state_vector
        mean = np.vdot(state_vector, applied)
        centered.append(applied - mean * state_vector)
    centered = np.stack(centered, axis=1)
    return (centered.conj().T @ centered).real


def local_phase_fisher_matrix(state_vector, site_count):
    generators = [local_operator(np.diag([0.0, 1.0]), site, site_count) for site in range(site_count)]
    return 4 * observable_covariance(state_vector, generators)


def thermal_state(hamiltonian, inverse_temperature):
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    weights = np.exp(-inverse_temperature * (eigenvalues - eigenvalues.min()))
    weights /= weights.sum()
    return (eigenvectors * weights) @ eigenvectors.conj().T


def recover_hamiltonian_from_thermal_state(state, inverse_temperature):
    if inverse_temperature <= 0:
        raise ValueError("Use a known positive inverse temperature.")
    eigenvalues, eigenvectors = np.linalg.eigh(state)
    if eigenvalues.min() <= 0:
        raise ValueError("The Gibbs state must have full rank in the chosen space.")
    recovered = (eigenvectors * (-np.log(eigenvalues) / inverse_temperature)) @ eigenvectors.conj().T
    return recovered - np.eye(state.shape[0]) * np.trace(recovered) / state.shape[0]


def reduced_density(state_vector, sites, site_count):
    selected = tuple(sites)
    remainder = tuple(site for site in range(site_count) if site not in selected)
    amplitudes = state_vector.reshape((2,) * site_count).transpose(selected + remainder)
    amplitudes = amplitudes.reshape(2 ** len(selected), -1)
    return amplitudes @ amplitudes.conj().T


def entropy_nats(state):
    eigenvalues = np.linalg.eigvalsh(state)
    if eigenvalues.min() < -TOLERANCE:
        raise ValueError("The state is not positive semidefinite.")
    positive = eigenvalues[eigenvalues > TOLERANCE]
    return float(-np.sum(positive * np.log(positive)))


def mutual_information_matrix(state_vector, site_count):
    single_entropies = [
        entropy_nats(reduced_density(state_vector, (site,), site_count))
        for site in range(site_count)
    ]
    mutual_information = np.zeros((site_count, site_count))
    for source in range(site_count):
        for target in range(source + 1, site_count):
            joint_entropy = entropy_nats(reduced_density(state_vector, (source, target), site_count))
            value = single_entropies[source] + single_entropies[target] - joint_entropy
            mutual_information[source, target] = value
            mutual_information[target, source] = value
    return mutual_information


def propagator(hamiltonian, time):
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    return (eigenvectors * np.exp(-1j * time * eigenvalues)) @ eigenvectors.conj().T


def phase_kick_probabilities(hamiltonian, source, phase, time):
    site_count = hamiltonian.shape[0]
    initial = np.full(site_count, 1 / math.sqrt(site_count), dtype=complex)
    initial[source] *= np.exp(1j * phase)
    evolved = propagator(hamiltonian, time) @ initial
    return np.abs(evolved) ** 2


def phase_response(hamiltonian, source, time):
    positive = phase_kick_probabilities(hamiltonian, source, math.pi / 2, time)
    negative = phase_kick_probabilities(hamiltonian, source, -math.pi / 2, time)
    return positive - negative


def reconstructed_weights(hamiltonian, time_step):
    if time_step <= 0:
        raise ValueError("The response time step must be positive.")
    site_count = hamiltonian.shape[0]
    weights = np.zeros_like(hamiltonian)
    for source in range(site_count):
        weights[:, source] = -site_count * phase_response(hamiltonian, source, time_step) / (4 * time_step)
    np.fill_diagonal(weights, 0.0)
    return weights


def graph_distances(site_count, edges):
    adjacency = [[] for _ in range(site_count)]
    for source, target, _ in edges:
        adjacency[source].append(target)
        adjacency[target].append(source)
    distances = np.full((site_count, site_count), np.inf)
    for source in range(site_count):
        distances[source, source] = 0
        frontier = deque([source])
        while frontier:
            current = frontier.popleft()
            for target in adjacency[current]:
                if np.isinf(distances[source, target]):
                    distances[source, target] = distances[source, current] + 1
                    frontier.append(target)
    return distances


def analytic_star_probability(time):
    return math.sin(2 * time) ** 2 / 4


def analytic_chain_probability(time):
    root_two = math.sqrt(2)
    amplitude = (
        0.25
        - (2 + root_two) * np.exp(-1j * (2 - root_two) * time) / 8
        + 0.25 * np.exp(-2j * time)
        - (2 - root_two) * np.exp(-1j * (2 + root_two) * time) / 8
    )
    return float(abs(amplitude) ** 2)


def hypothesis_protocol():
    return {
        "question": "Does one static state plus fixed site factorization determine transport, without a generator?",
        "scope": "Connected positive weighted graph Laplacians, restricted to the conserved one-excitation sector.",
        "units": "hbar=1; edge weights are energies in a fixed energy unit; time is in its inverse unit.",
        "not_tested": "General state-plus-update reconstruction, gravitational dynamics, consciousness, or a Lorentzian continuum.",
        "static_disconfirmation": "Equal full ground state and pairwise mutual information, unequal labeled transition probabilities.",
        "dynamic_check": "At short times J_ab = -(N/4) d[p_b(+pi/2)-p_b(-pi/2)]/dt at zero for distinct sites.",
        "fixed_comparison": "N=4, connected nonisomorphic chain and star, three unit-weight edges and equal Hamiltonian trace.",
        "failure_conditions": [
            "A graph does not have the claimed unique zero-energy ground state within the one-excitation sector.",
            "The two Hamiltonians give identical transport in the selected labeled experiment.",
            "The independently evolved phase-kick states disagree with the analytic response identity.",
            "The reconstructed weights do not converge to the actual weights as the time step decreases.",
        ],
        "evidence_type": "Exact finite-model identities and numerical checks, not laboratory observations or a novelty claim.",
    }


def experiment_results():
    state = uniform_one_excitation_state(SITE_COUNT)
    common_mutual_information = mutual_information_matrix(state, SITE_COUNT)
    fisher = local_phase_fisher_matrix(state, SITE_COUNT)
    edge_basis = [
        full_spin_hamiltonian(SITE_COUNT, ((source, target, 1.0),))
        for source in range(SITE_COUNT)
        for target in range(source + 1, SITE_COUNT)
    ]
    covariance = observable_covariance(state, edge_basis)
    models = {}
    for name, edges in (("chain", CHAIN_EDGES), ("star", STAR_EDGES)):
        hamiltonian = graph_laplacian(SITE_COUNT, edges)
        spectrum, eigenvectors = np.linalg.eigh(hamiltonian)
        true_weights = np.diag(np.diag(hamiltonian)) - hamiltonian
        selected_times = (0.025, 0.05, 0.1, math.pi / 4, 1.0)
        probabilities = [float(abs(propagator(hamiltonian, time)[3, 0]) ** 2) for time in selected_times]
        estimated_power = math.log(probabilities[1] / probabilities[0]) / math.log(2)
        response_step = 1e-4
        learned = reconstructed_weights(hamiltonian, response_step)
        learned_edges = [
            (source, target, 1.0)
            for source in range(SITE_COUNT)
            for target in range(source + 1, SITE_COUNT)
            if learned[source, target] > 0.5
        ]
        models[name] = {
            "edges": edges,
            "hamiltonian": hamiltonian.tolist(),
            "spectrum": spectrum.tolist(),
            "trace": float(np.trace(hamiltonian)),
            "ground_state_fidelity_with_uniform_sector_state": float(abs(np.sum(eigenvectors[:, 0]) / math.sqrt(SITE_COUNT)) ** 2),
            "common_pairwise_mutual_information_nats": common_mutual_information.tolist(),
            "graph_distances_input": graph_distances(SITE_COUNT, edges).astype(int).tolist(),
            "transport_source": 0,
            "transport_target": 3,
            "times": selected_times,
            "transport_probabilities": probabilities,
            "estimated_short_time_probability_power": estimated_power,
            "exact_leading_probability_power": 2 * int(graph_distances(SITE_COUNT, edges)[0, 3]),
            "local_phase_response_at_time_0_05": float(phase_response(hamiltonian, 0, 0.05)[3]),
            "response_step": response_step,
            "reconstructed_weights": learned.tolist(),
            "max_weight_error": float(np.max(np.abs(learned - true_weights))),
            "recovered_graph_distances": graph_distances(SITE_COUNT, learned_edges).astype(int).tolist(),
            "response_convergence": [
                {"time_step": time_step, "max_weight_error": float(np.max(np.abs(reconstructed_weights(hamiltonian, time_step) - true_weights)))}
                for time_step in (0.04, 0.02, 0.01, 0.005)
            ],
        }
    return {
        "protocol": hypothesis_protocol(),
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "static_identifiability": {
            "local_phase_fisher_matrix": fisher.tolist(),
            "local_phase_fisher_eigenvalues": np.linalg.eigvalsh(fisher).tolist(),
            "allowed_edge_operator_count": len(edge_basis),
            "edge_operator_covariance": covariance.tolist(),
            "edge_operator_covariance_rank": int(np.linalg.matrix_rank(covariance, tol=TOLERANCE)),
            "edge_operator_covariance_nullity": len(edge_basis) - int(np.linalg.matrix_rank(covariance, tol=TOLERANCE)),
            "caveat": "A special non-identifiable eigenstate family; not a theorem that static state reconstruction always fails.",
        },
        "models": models,
        "shot_noise_illustration_not_data": {
            "shots_per_phase": 1_000_000,
            "time_step": 0.001,
            "approximate_weight_standard_error": SITE_COUNT / (4 * 0.001) * math.sqrt(2 * 0.25 * 0.75 / 1_000_000),
            "assumptions": "Independent binomial counts, short time p~1/N, no drift, perfect preparation and phase control.",
        },
    }


class GeometryResearchTests(unittest.TestCase):
    def test_full_spin_model_conserves_number_and_restricts_to_laplacian(self):
        number = sum(local_operator(np.diag([0.0, 1.0]), site, SITE_COUNT) for site in range(SITE_COUNT))
        indices = [1 << (SITE_COUNT - site - 1) for site in range(SITE_COUNT)]
        ground = uniform_one_excitation_state(SITE_COUNT)
        vacuum = np.zeros(2**SITE_COUNT)
        vacuum[0] = 1
        for edges in (CHAIN_EDGES, STAR_EDGES):
            full_hamiltonian = full_spin_hamiltonian(SITE_COUNT, edges)
            projected = full_hamiltonian[np.ix_(indices, indices)]
            np.testing.assert_allclose(projected, graph_laplacian(SITE_COUNT, edges), atol=TOLERANCE)
            np.testing.assert_allclose(full_hamiltonian @ number, number @ full_hamiltonian, atol=TOLERANCE)
            np.testing.assert_allclose(full_hamiltonian @ ground, 0, atol=TOLERANCE)
            np.testing.assert_allclose(full_hamiltonian @ vacuum, 0, atol=TOLERANCE)
            self.assertGreaterEqual(float(np.linalg.eigvalsh(full_hamiltonian).min()), -TOLERANCE)

    def test_full_tensor_phase_experiment_matches_restricted_calculation(self):
        ground = uniform_one_excitation_state(SITE_COUNT)
        phase = math.pi / 2
        for edges in (CHAIN_EDGES, STAR_EDGES):
            full_hamiltonian = full_spin_hamiltonian(SITE_COUNT, edges)
            for source in range(SITE_COUNT):
                occupation = local_operator(np.diag([0.0, 1.0]), source, SITE_COUNT)
                kick = np.eye(2**SITE_COUNT) + (np.exp(1j * phase) - 1) * occupation
                evolved = propagator(full_hamiltonian, 0.37) @ kick @ ground
                measured = [
                    float(np.vdot(evolved, local_operator(np.diag([0.0, 1.0]), target, SITE_COUNT) @ evolved).real)
                    for target in range(SITE_COUNT)
                ]
                expected = phase_kick_probabilities(graph_laplacian(SITE_COUNT, edges), source, phase, 0.37)
                np.testing.assert_allclose(measured, expected, atol=TOLERANCE)

    def test_local_phase_fisher_matches_formula_and_fidelity(self):
        ground = uniform_one_excitation_state(SITE_COUNT)
        fisher = local_phase_fisher_matrix(ground, SITE_COUNT)
        expected = 4 * (np.eye(SITE_COUNT) / SITE_COUNT - np.ones((SITE_COUNT, SITE_COUNT)) / SITE_COUNT**2)
        np.testing.assert_allclose(fisher, expected, atol=TOLERANCE)
        np.testing.assert_allclose(np.linalg.eigvalsh(fisher), [0, 1, 1, 1], atol=TOLERANCE)
        for direction in (np.array([1, 0, 0, 0]), np.arange(SITE_COUNT), np.ones(SITE_COUNT)):
            step = 1e-4
            overlap = np.mean(np.exp(1j * step * direction))
            fidelity_estimate = 4 * (1 - abs(overlap) ** 2) / step**2
            self.assertAlmostEqual(fidelity_estimate, float(direction @ fisher @ direction), delta=1e-6)

    def test_static_edge_covariance_has_no_identifying_rank(self):
        state = uniform_one_excitation_state(SITE_COUNT)
        operators = [
            full_spin_hamiltonian(SITE_COUNT, ((source, target, 1.0),))
            for source in range(SITE_COUNT)
            for target in range(source + 1, SITE_COUNT)
        ]
        covariance = observable_covariance(state, operators)
        np.testing.assert_allclose(covariance, np.zeros((6, 6)), atol=TOLERANCE)
        self.assertEqual(int(np.linalg.matrix_rank(covariance, tol=TOLERANCE)), 0)

    def test_known_full_rank_gibbs_state_is_a_positive_static_control(self):
        states = []
        for edges in (CHAIN_EDGES, STAR_EDGES):
            hamiltonian = graph_laplacian(SITE_COUNT, edges)
            state = thermal_state(hamiltonian, 0.7)
            states.append(state)
            recovered = recover_hamiltonian_from_thermal_state(state, 0.7)
            traceless = hamiltonian - np.eye(SITE_COUNT) * np.trace(hamiltonian) / SITE_COUNT
            np.testing.assert_allclose(recovered, traceless, atol=TOLERANCE)
        self.assertGreater(float(np.linalg.norm(states[0] - states[1])), 0.1)

    def test_periodic_chain_dispersion_is_quadratic_not_relativistic(self):
        site_count = 12
        edges = tuple((site, (site + 1) % site_count, 1.0) for site in range(site_count))
        spectrum = np.linalg.eigvalsh(graph_laplacian(site_count, edges))
        momenta = 2 * math.pi * np.arange(site_count) / site_count
        expected = np.sort(2 * (1 - np.cos(momenta)))
        np.testing.assert_allclose(spectrum, expected, atol=TOLERANCE)
        momentum = 0.001
        self.assertAlmostEqual(2 * (1 - math.cos(momentum)) / momentum**2, 1, delta=1e-6)

    def test_connected_graphs_share_unique_sector_ground_state(self):
        uniform = np.full(SITE_COUNT, 1 / math.sqrt(SITE_COUNT))
        for edges in (CHAIN_EDGES, STAR_EDGES):
            hamiltonian = graph_laplacian(SITE_COUNT, edges)
            spectrum = np.linalg.eigvalsh(hamiltonian)
            np.testing.assert_allclose(hamiltonian @ uniform, 0, atol=TOLERANCE)
            self.assertEqual(int(np.sum(abs(spectrum) < TOLERANCE)), 1)
            self.assertGreater(spectrum[1], 0)
            self.assertAlmostEqual(float(np.trace(hamiltonian)), 6)

    def test_reduced_states_and_mutual_information_are_not_product_trick(self):
        state = uniform_one_excitation_state(SITE_COUNT)
        self.assertAlmostEqual(float(np.vdot(state, state).real), 1)
        expected_pair = np.array([[0.5, 0, 0, 0], [0, 0.25, 0.25, 0], [0, 0.25, 0.25, 0], [0, 0, 0, 0]])
        for source in range(SITE_COUNT):
            for target in range(source + 1, SITE_COUNT):
                np.testing.assert_allclose(reduced_density(state, (source, target), SITE_COUNT), expected_pair, atol=TOLERANCE)
        mutual_information = mutual_information_matrix(state, SITE_COUNT)
        binary_entropy_quarter = -0.25 * math.log(0.25) - 0.75 * math.log(0.75)
        expected = 2 * binary_entropy_quarter - math.log(2)
        np.testing.assert_allclose(mutual_information[~np.eye(SITE_COUNT, dtype=bool)], expected)
        self.assertGreater(expected, 0.4)

    def test_spectra_match_closed_forms(self):
        chain = np.linalg.eigvalsh(graph_laplacian(SITE_COUNT, CHAIN_EDGES))
        star = np.linalg.eigvalsh(graph_laplacian(SITE_COUNT, STAR_EDGES))
        np.testing.assert_allclose(chain, [0, 2 - math.sqrt(2), 2, 2 + math.sqrt(2)], atol=TOLERANCE)
        np.testing.assert_allclose(star, [0, 1, 1, 4], atol=TOLERANCE)

    def test_propagator_unitary_and_stationary_ground(self):
        uniform = np.full(SITE_COUNT, 1 / math.sqrt(SITE_COUNT))
        for edges in (CHAIN_EDGES, STAR_EDGES):
            evolution = propagator(graph_laplacian(SITE_COUNT, edges), 0.73)
            np.testing.assert_allclose(evolution.conj().T @ evolution, np.eye(SITE_COUNT), atol=TOLERANCE)
            np.testing.assert_allclose(evolution @ uniform, uniform, atol=TOLERANCE)

    def test_transport_matches_independent_analytic_expressions(self):
        for edges, analytic in ((CHAIN_EDGES, analytic_chain_probability), (STAR_EDGES, analytic_star_probability)):
            hamiltonian = graph_laplacian(SITE_COUNT, edges)
            for time in (0.05, 0.2, math.pi / 4, 1, 2):
                probability = float(abs(propagator(hamiltonian, time)[3, 0]) ** 2)
                self.assertAlmostEqual(probability, analytic(time), places=12)

    def test_distinct_transport_for_identical_static_state(self):
        chain_probability = analytic_chain_probability(math.pi / 4)
        star_probability = analytic_star_probability(math.pi / 4)
        self.assertAlmostEqual(star_probability, 0.25)
        self.assertGreater(abs(star_probability - chain_probability), 0.2)

    def test_leading_powers_follow_graph_distance(self):
        for edges, expected_distance in ((CHAIN_EDGES, 3), (STAR_EDGES, 1)):
            hamiltonian = graph_laplacian(SITE_COUNT, edges)
            for power in range(expected_distance):
                self.assertAlmostEqual(float(np.linalg.matrix_power(hamiltonian, power)[3, 0]), 0)
            self.assertAlmostEqual(float(np.linalg.matrix_power(hamiltonian, expected_distance)[3, 0]), (-1) ** expected_distance)
            time = 0.025
            probability = float(abs(propagator(hamiltonian, time)[3, 0]) ** 2)
            doubled = float(abs(propagator(hamiltonian, 2 * time)[3, 0]) ** 2)
            self.assertAlmostEqual(math.log(doubled / probability) / math.log(2), 2 * expected_distance, delta=0.01)

    def test_phase_response_identity_uses_independent_kicked_states(self):
        for edges in (CHAIN_EDGES, STAR_EDGES):
            hamiltonian = graph_laplacian(SITE_COUNT, edges)
            for time in (0, 0.02, 0.3, 1.7):
                evolution = propagator(hamiltonian, time)
                for source in range(SITE_COUNT):
                    response = phase_response(hamiltonian, source, time)
                    expected = -4 * evolution[:, source].imag / SITE_COUNT
                    np.testing.assert_allclose(response, expected, atol=TOLERANCE)

    def test_response_recovers_weights_and_graph(self):
        for edges in (CHAIN_EDGES, STAR_EDGES):
            hamiltonian = graph_laplacian(SITE_COUNT, edges)
            expected_weights = np.diag(np.diag(hamiltonian)) - hamiltonian
            learned = reconstructed_weights(hamiltonian, 1e-4)
            np.testing.assert_allclose(learned, expected_weights, atol=1e-6)
            np.testing.assert_array_equal(learned > 0.5, expected_weights > 0.5)

    def test_response_has_second_order_bias(self):
        for edges in (CHAIN_EDGES, STAR_EDGES):
            hamiltonian = graph_laplacian(SITE_COUNT, edges)
            expected = np.diag(np.diag(hamiltonian)) - hamiltonian
            coarse_error = np.max(abs(reconstructed_weights(hamiltonian, 0.02) - expected))
            fine_error = np.max(abs(reconstructed_weights(hamiltonian, 0.01) - expected))
            self.assertAlmostEqual(float(coarse_error / fine_error), 4, delta=0.01)

    def test_positive_weighted_graph_response(self):
        edges = ((0, 1, 0.7), (1, 2, 1.3), (2, 3, 0.4), (3, 0, 0.9), (0, 2, 0.2))
        hamiltonian = graph_laplacian(SITE_COUNT, edges)
        expected = np.diag(np.diag(hamiltonian)) - hamiltonian
        np.testing.assert_allclose(reconstructed_weights(hamiltonian, 1e-4), expected, atol=1e-6)

    def test_zero_time_no_remote_response(self):
        for edges in (CHAIN_EDGES, STAR_EDGES):
            hamiltonian = graph_laplacian(SITE_COUNT, edges)
            for source in range(SITE_COUNT):
                np.testing.assert_allclose(phase_response(hamiltonian, source, 0), 0, atol=TOLERANCE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    arguments = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(GeometryResearchTests)
    checks = unittest.TextTestRunner(verbosity=2).run(suite)
    if not checks.wasSuccessful():
        raise SystemExit(1)
    results = experiment_results()
    results["automated_checks"] = {"run": checks.testsRun, "failures": 0, "errors": 0}
    encoded = json.dumps(results, indent=2)
    if arguments.write_results:
        output = Path(__file__).resolve().parent / "results.json"
        output.write_text(encoded + "\n", encoding="utf-8")
        print(f"Results written to {output}")
    compact = {
        "tests_passed": checks.testsRun,
        "shared_pairwise_mutual_information_nats": results["models"]["chain"]["common_pairwise_mutual_information_nats"][0][1],
        "static_identifiability": results["static_identifiability"],
        "models": {
            name: {key: model[key] for key in ("spectrum", "trace", "ground_state_fidelity_with_uniform_sector_state", "transport_probabilities", "estimated_short_time_probability_power", "max_weight_error", "local_phase_response_at_time_0_05")}
            for name, model in results["models"].items()
        },
        "shot_noise_illustration_not_data": results["shot_noise_illustration_not_data"],
        "interpretation": results["protocol"]["evidence_type"],
    }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()