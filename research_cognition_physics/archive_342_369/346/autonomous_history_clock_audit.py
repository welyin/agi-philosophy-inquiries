"""Round 346: finite programmed Feynman clock, transfer and resource bounds."""
import functools
import itertools
import math
import unittest
import numpy as np
from growing_stream_audit import main


def catalogue(name):
    identity = np.eye(2, dtype=complex)
    hadamard = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
    flip = np.array([[0, 1], [1, 0]], complex)
    if name == "real":
        return (identity, hadamard, flip)
    if name == "binary":
        return (identity, flip)
    phase = np.diag([1, np.exp(1j*np.pi/4)])
    cnot = np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                    [0, 0, 0, 1], [0, 0, 1, 0]], complex)
    return (np.kron(hadamard, identity), np.kron(phase, identity), cnot)


def exponential(hamiltonian, time):
    values, vectors = np.linalg.eigh(hamiltonian)
    return (vectors * np.exp(-1j*time*values)) @ vectors.conj().T


def hop_weights(length, omega=1.):
    return omega/2 * np.sqrt(np.arange(1, length+1) *
                              np.arange(length, 0, -1))


def bare_clock(weights):
    return np.diag(weights, 1) + np.diag(weights, -1)


def capped_omega(length, cap=1.):
    peak_factor = math.sqrt(((length+1)**2)//4)
    return 2*cap/peak_factor


@functools.lru_cache(maxsize=None)
def processor(length=2, name="complex", omega=1.):
    gates = catalogue(name)
    q, d = len(gates), len(gates[0])
    programs = tuple(itertools.product(range(q), repeat=length))
    internal = len(programs)*d
    controlled = []
    for slot in range(length):
        operator = np.zeros((internal, internal), complex)
        for index, program in enumerate(programs):
            block = slice(index*d, (index+1)*d)
            operator[block, block] = gates[program[slot]]
        controlled.append(operator)
    weights = hop_weights(length, omega)
    dimension = (length+1)*internal
    hamiltonian = np.zeros((dimension, dimension), complex)
    dressing = np.zeros_like(hamiltonian)
    cumulative = [np.eye(internal, dtype=complex)]
    for index, operator in enumerate(controlled):
        left = slice(index*internal, (index+1)*internal)
        right = slice((index+1)*internal, (index+2)*internal)
        hamiltonian[right, left] = weights[index]*operator
        hamiltonian[left, right] = weights[index]*operator.conj().T
        cumulative.append(operator @ cumulative[-1])
    for clock, unitary in enumerate(cumulative):
        block = slice(clock*internal, (clock+1)*internal)
        dressing[block, block] = unitary
    return {"H": hamiltonian, "W": dressing, "V": controlled,
            "U": cumulative, "programs": programs, "data_dimension": d,
            "internal_dimension": internal, "weights": weights}


def basis(index, dimension):
    result = np.zeros(dimension, complex)
    result[index] = 1
    return result


def circuit(gates, program):
    result = np.eye(len(gates[0]), dtype=complex)
    for symbol in program:
        result = gates[symbol] @ result
    return result


def initial_isometry(clock, program, model):
    nclock = len(model["weights"])+1
    labels = np.kron(basis(clock, nclock),
                     basis(program, len(model["programs"])))
    return np.kron(labels.reshape(-1, 1),
                   np.eye(model["data_dimension"]))


def clock_probabilities(vector, length, internal):
    return np.sum(np.abs(vector.reshape(length+1, internal))**2, axis=1)


def binomial_clock(length, omega, time):
    sine, cosine = np.sin(omega*time/2), np.cos(omega*time/2)
    return np.array([math.comb(length, j)*sine**(2*j) *
                     cosine**(2*(length-j)) for j in range(length+1)])


def reference_error(model, program, time):
    d = model["data_dimension"]
    unitary = exponential(model["H"], time)
    incoming = initial_isometry(0, program, model)
    outgoing = initial_isometry(len(model["weights"]), program, model)
    target = model["U"][-1][program*d:(program+1)*d,
                             program*d:(program+1)*d]
    # Maximally entangled input tests the complete input channel, not only
    # one product input. Output phase is physically irrelevant.
    bell = np.eye(d).reshape(-1)/np.sqrt(d)
    actual = np.kron(unitary @ incoming, np.eye(d)) @ bell
    ideal = np.kron(outgoing @ target, np.eye(d)) @ bell
    phase = (-1j)**len(model["weights"])
    return float(np.linalg.norm(actual-phase*ideal))


class Checks(unittest.TestCase):
    def test_01_fixed_catalogue_and_controlled_unitarity(self):
        model = processor()
        for operator in (*catalogue("complex"), *model["V"], model["W"]):
            self.assertLess(np.linalg.norm(operator.conj().T @ operator -
                                            np.eye(len(operator))), 1e-12)

    def test_02_real_matrix_direct_exponential_against_dressed_chain(self):
        model = processor(2, "real")
        self.assertLess(np.max(np.abs(model["H"].imag)), 1e-15)
        bare = bare_clock(model["weights"])
        for time in (.31, 1.2, np.pi):
            expected = model["W"] @ np.kron(exponential(bare, time),
                np.eye(model["internal_dimension"])) @ model["W"].conj().T
            self.assertLess(np.linalg.norm(exponential(model["H"], time) -
                                           expected), 3e-13)

    def test_03_complex_matrix_independent_hamiltonian_conjugacy(self):
        model = processor()
        self.assertLess(np.linalg.norm(model["H"]-model["H"].conj().T), 1e-13)
        expected = model["W"] @ np.kron(bare_clock(model["weights"]),
            np.eye(model["internal_dimension"])) @ model["W"].conj().T
        self.assertLess(np.linalg.norm(model["H"]-expected), 1e-12)

    def test_04_full_spectrum_and_multiplicities(self):
        for length, name in ((1, "complex"), (2, "complex"), (3, "binary")):
            model = processor(length, name, .7)
            expected = np.repeat(.7*np.arange(-length/2, length/2+1),
                                 model["internal_dimension"])
            self.assertLess(np.max(np.abs(np.linalg.eigvalsh(model["H"]) -
                                         expected)), 1e-12)

    def test_05_every_program_same_hamiltonian_exact_operator_output(self):
        model = processor()
        unitary = exponential(model["H"], np.pi)
        for index, program in enumerate(model["programs"]):
            actual = unitary @ initial_isometry(0, index, model)
            expected = -initial_isometry(2, index, model) @ circuit(
                catalogue("complex"), program)
            self.assertLess(np.linalg.norm(actual-expected), 2e-13)

    def test_06_complete_channel_with_entangled_reference(self):
        model = processor()
        for index in (0, 1, 2, 5, 8):
            self.assertLess(reference_error(model, index, np.pi), 2e-13)

    def test_07_mixed_data_reference_preserves_target_joint_action(self):
        model = processor(1, "binary")
        rng = np.random.default_rng(34607)
        raw = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
        rho = raw @ raw.conj().T
        rho /= np.trace(rho)
        incoming = initial_isometry(0, 1, model)
        actual_map = np.kron(exponential(model["H"], np.pi) @ incoming,
                            np.eye(2))
        expected_map = np.kron(initial_isometry(1, 1, model) @
                              catalogue("binary")[1], np.eye(2))
        self.assertLess(np.linalg.norm(actual_map @ rho @ actual_map.conj().T -
                         expected_map @ rho @ expected_map.conj().T), 2e-13)

    def test_08_coherent_program_is_not_an_unchanged_catalyst(self):
        model = processor(1, "binary")
        vector = np.kron(np.kron(basis(0, 2),
            np.array([1, 1])/np.sqrt(2)), basis(0, 2))
        final = exponential(model["H"], np.pi) @ vector
        clock_program_data = final.reshape(2, 2, 2)
        program_rho = np.einsum("cpd,cqd->pq", clock_program_data,
                               clock_program_data.conj())
        self.assertLess(np.linalg.norm(program_rho-np.eye(2)/2), 2e-13)
        self.assertAlmostEqual(float(np.trace(program_rho@program_rho).real),
                               .5, places=13)

    def test_09_program_projectors_conserved_at_all_times(self):
        model = processor()
        d = model["data_dimension"]
        for index in range(len(model["programs"])):
            labels = np.outer(basis(index, len(model["programs"])),
                              basis(index, len(model["programs"])).conj())
            projector = np.kron(np.eye(3), np.kron(labels, np.eye(d)))
            self.assertLess(np.linalg.norm(projector@model["H"] -
                                           model["H"]@projector), 1e-13)

    def test_10_clock_distribution_matches_binomial_formula(self):
        for length in (1, 2, 5, 12):
            omega = capped_omega(length)
            for fraction in (.2, .5, .9, 1., 1.3):
                time = fraction*np.pi/omega
                clock = exponential(bare_clock(hop_weights(length, omega)),
                                    time) @ basis(0, length+1)
                self.assertLess(np.linalg.norm(abs(clock)**2 -
                    binomial_clock(length, omega, time)), 2e-13)

    def test_11_per_hop_cap_imposes_linear_completion_time(self):
        for length in range(1, 65):
            omega = capped_omega(length, .8)
            self.assertAlmostEqual(float(np.max(hop_weights(length, omega))),
                                   .8, places=13)
            time = np.pi/omega
            self.assertGreaterEqual(time, np.pi*length/(4*.8)-1e-13)
            self.assertLessEqual(time, np.pi*(length+1)/(4*.8)+1e-13)

    def test_12_pointer_velocity_bound_for_nonengineered_chains(self):
        rng = np.random.default_rng(34612)
        for length in (3, 6, 10):
            weights = rng.uniform(.1, 1., size=length)
            h = bare_clock(weights)
            pointer = np.diag(np.arange(length+1))
            velocity = 1j*(h@pointer-pointer@h)
            cap = float(np.max(weights))
            self.assertLessEqual(np.linalg.norm(velocity, 2), 2*cap+1e-13)
            for time in (.2, .7, 1.3):
                vector = exponential(h, time) @ basis(0, length+1)
                mean_pointer = float(np.vdot(vector, pointer@vector).real)
                self.assertLessEqual(mean_pointer, 2*cap*time+1e-13)

    def test_13_total_energy_moments_are_conserved(self):
        model = processor()
        vector = initial_isometry(0, 2, model) @ (np.ones(4)/2)
        h = model["H"]
        for time in (.3, 1.5, np.pi, 2*np.pi):
            evolved = exponential(h, time) @ vector
            self.assertLess(abs(np.vdot(evolved, h@evolved)), 2e-13)
            self.assertAlmostEqual(float(np.vdot(h@evolved, h@evolved).real),
                                   float(model["weights"][0]**2), places=12)

    def test_14_exact_recurrence_undoes_the_computation(self):
        model = processor()
        unitary = exponential(model["H"], 2*np.pi)
        self.assertLess(np.linalg.norm(unitary-np.eye(len(unitary))), 3e-13)
        incoming = initial_isometry(0, 2, model) @ (np.ones(4)/2)
        final = unitary @ incoming
        distribution = clock_probabilities(final, 2, model["internal_dimension"])
        self.assertAlmostEqual(float(distribution[0]), 1., places=13)
        self.assertLess(float(distribution[-1]), 1e-24)

    def test_15_engineering_is_required_not_any_path_has_exact_arrival(self):
        weights = hop_weights(6)
        weights[2] *= .8
        amplitude = exponential(bare_clock(weights), np.pi)[-1, 0]
        self.assertLess(float(abs(amplitude)**2), .95)

    def test_16_gate_order_reversal_is_detectable(self):
        gates = catalogue("complex")
        program = (0, 1)
        correct = circuit(gates, program)
        reversed_order = circuit(gates, tuple(reversed(program)))
        self.assertGreater(np.linalg.norm(correct-reversed_order), .5)
        model = processor()
        index = model["programs"].index(program)
        actual = exponential(model["H"], np.pi) @ initial_isometry(0, index, model)
        wrong = -initial_isometry(2, index, model) @ reversed_order
        self.assertGreater(np.linalg.norm(actual-wrong), .5)

    def test_17_readout_timing_bound_and_preterminal_channel(self):
        for length in (2, 8, 32):
            omega = capped_omega(length)
            epsilon = .01
            delta = 2*np.sqrt(epsilon/length)/omega
            probability = np.cos(omega*delta/2)**(2*length)
            self.assertGreaterEqual(probability, 1-epsilon-1e-14)
        model = processor()
        time = .47*np.pi
        actual = exponential(model["H"], time) @ initial_isometry(0, 1, model)
        data_vector = np.array([1, 1j, .3, -.7j], complex)
        data_vector /= np.linalg.norm(data_vector)
        input_state = np.outer(data_vector, data_vector.conj())
        joint = actual @ input_state @ actual.conj().T
        diagonal = joint.reshape(3, 9, 4, 3, 9, 4)
        data = np.einsum("cpdcpe->de", diagonal)
        probabilities = binomial_clock(2, 1., time)
        gates = catalogue("complex")
        program = model["programs"][1]
        expected = sum(probabilities[j]*circuit(gates, program[:j]) @ input_state
            @ circuit(gates, program[:j]).conj().T for j in range(3))
        self.assertLess(np.linalg.norm(data-expected), 2e-13)


def report():
    model = processor()
    transfers = [{"program": list(program),
                  "reference_vector_error": reference_error(model, index, np.pi)}
                 for index, program in enumerate(model["programs"])]
    budgets = []
    for length in (1, 2, 4, 8, 16, 32, 64):
        omega = capped_omega(length)
        budgets.append({"gates": length, "max_hop_norm": 1.,
                        "omega": omega, "exact_arrival_time": float(np.pi/omega),
                        "all_path_pointer_lower_bound": length/2,
                        "hamiltonian_norm": omega*length/2,
                        "spectral_width": omega*length,
                        "dimension_for_three_symbols_and_four_data_states":
                            (length+1)*3**length*4,
                        "readout_delta_for_failure_bound_0.01":
                            float(2*np.sqrt(.01/length)/omega)})
    return {"round": 346,
            "scope": "Engineered finite autonomous processor with internal "
                     "orthogonal program and exact arrival; no derived "
                     "spatial locality or GR and no autonomous preparation/readout.",
            "full_demonstrator_dimension": len(model["H"]),
            "catalogue": ["H_on_first_qubit", "T_on_first_qubit", "CNOT"],
            "program_count": len(model["programs"]),
            "arrival_time_omega_one": float(np.pi),
            "program_transfers": transfers,
            "resource_budgets": budgets,
            "clock_return_time_omega_one": float(2*np.pi),
            "classical_program_label_preserved": True,
            "arbitrary_coherent_program_state_preserved": False,
            "added_inputs": ["chosen unitary catalogue",
                            "orthogonal length-L program memory",
                            "L+1 clock levels initially localized at zero",
                            "engineered square-root hopping couplings",
                            "prepared unknown data and any external reference",
                            "calibrated readout time and access",
                            "no claim of bounded geometric interaction range"]}


if __name__ == "__main__":
    main(__name__, "autonomous_history_clock_audit", report)
