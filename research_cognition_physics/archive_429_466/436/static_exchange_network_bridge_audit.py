"""Round 436: apply static Hamiltonian simulation to the full round-433 network.

The full Heisenberg-exchange existence claim uses the cited CMP theorem chain.
This code constructs and checks its first, six-mediator, no-Y layer, certifies
the target response, and audits the complete-input dynamical error contract.
It does not claim to output or diagonalize a final Heisenberg coupling list.
"""
import argparse
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import exchange_relation_audit as core
import quantum_participation_audit as network
from encoded_exchange_response_audit import norm

TARGET = Path(__file__).with_name('static_exchange_network_bridge_audit_results.json')
OBS = {}
PAULI = dict(zip('IXYZ', [np.eye(2), *core.PAULI]))


def matrix(word):
    out = np.array([[1.]])
    for letter in word:
        out = np.kron(out, PAULI[letter])
    return out


def target_terms():
    terms = {'IIIIII': Fraction(9, 4)}
    for e, (i, j) in enumerate(network.EDGES):
        for letter, coefficient in [('X', Fraction(1)), ('Z', -Fraction(3, 4))]:
            word = ['I']*6
            word[3+e] = letter
            terms[''.join(word)] = coefficient
        for letter in 'XYZ':
            word = ['I']*6
            word[i] = word[j] = letter
            terms[''.join(word)] = Fraction(1, 4)
            word[3+e] = 'Z'
            terms[''.join(word)] = -Fraction(1, 4)
    return terms


class FirstLayer:
    """Explicit 12-qubit intermediate; sparse Pauli actions on code columns."""
    def __init__(self):
        self.h = network.model()[0].real
        terms = target_terms()
        self.y_terms = [(word, c) for word, c in terms.items() if 'Y' in word]
        self.b_terms = {word+'I'*6: float(c) for word, c in terms.items() if 'Y' not in word}
        self.b_terms['I'*12] += sum(float(abs(c)) for _, c in self.y_terms)
        self.a_terms = {}
        self.k = []
        for q, (word, coefficient) in enumerate(self.y_terms):
            positions = [i for i, letter in enumerate(word) if letter == 'Y']
            assert len(positions) == 2
            left = ['I']*6
            right = list(word)
            for i in positions:
                left[i] = 'X'
                right[i] = 'Z'
            left, right = ''.join(left), ''.join(right)
            amplitude = math.sqrt(float(abs(coefficient))/2)
            sign = 1 if coefficient > 0 else -1
            self.k.append(amplitude*(matrix(left).real+sign*matrix(right).real))
            ancilla = ['I']*6
            ancilla[q] = 'X'
            self.a_terms[left+''.join(ancilla)] = amplitude
            self.a_terms[right+''.join(ancilla)] = sign*amplitude
        self.dimension = 4096
        self.indices = np.arange(self.dimension, dtype=np.int64)
        self.heavy = np.array([(j & 63).bit_count() for j in range(self.dimension)], float)
        self.v = np.zeros((4096, 64))
        self.v[64*np.arange(64), np.arange(64)] = 1
        self.actions = {}

    def apply(self, terms, columns):
        out = np.zeros_like(columns, dtype=np.result_type(columns.dtype, float))
        for word, coefficient in terms.items():
            if word not in self.actions:
                assert set(word) <= {'I', 'X', 'Z'}
                xmask = sum(1 << (11-j) for j, letter in enumerate(word) if letter == 'X')
                zpositions = [11-j for j, letter in enumerate(word) if letter == 'Z']
                parity = np.zeros(self.dimension, np.int64)
                for position in zpositions:
                    parity ^= (self.indices >> position) & 1
                self.actions[word] = (self.indices ^ xmask, 1-2*parity)
            permutation, phase = self.actions[word]
            out += coefficient*phase[:, None]*columns[permutation]
        return out

    def a(self, columns):
        return self.apply(self.a_terms, columns)

    def b(self, columns):
        return self.apply(self.b_terms, columns)

    def h0(self, columns):
        return self.heavy[:, None]*columns

    def inverse(self, columns):
        out = np.zeros_like(columns)
        np.divide(columns, self.heavy[:, None], out=out, where=self.heavy[:, None] > 0)
        return out


def first_layer_bound(mu, time):
    e, t = abs(mu), abs(time)
    return math.sqrt(3)*(2+39*t)*e+math.sqrt(6)*(6+63*t)*e*e


def exact_target_certificate(order=60):
    h, swaps, occupancies, blank = network.model()
    h_int = np.rint(h.real).astype(int)
    o = np.kron(np.eye(8), occupancies[0]-occupancies[1]).real.astype(int)
    rho_scaled = np.kron(np.eye(8)-swaps[0], blank).real.astype(int)
    nonzero_h = [(i, j, int(h_int[i, j])) for i, j in zip(*np.nonzero(h_int))]
    rho_entries = [(i, j, int(rho_scaled[i, j])) for i, j in zip(*np.nonzero(rho_scaled))]
    ad = o.astype(object)
    total, traces = Fraction(0), []
    for k in range(order+1):
        trace = sum(c*int(ad[j, i]) for i, j, c in rho_entries)
        traces.append(trace)
        if k % 2:
            assert trace == 0
        else:
            total += Fraction((-1)**(k//2)*trace, 4*2**k*math.factorial(k))
        next_ad = np.zeros((64, 64), dtype=object)
        for i, j, c in nonzero_h:
            next_ad[i] += c*ad[j]
            next_ad[:, j] -= c*ad[:, i]
        ad = next_ad
    # ||H|| <= 9, ||O||=1, t=1/2; exp(9) < 3^9.
    tail = Fraction(3**9*9**(order+1), math.factorial(order+1))
    lower, upper = total-tail, total+tail
    assert lower > Fraction(7, 500)
    return dict(order=order, time='1/2', scaled_low_order_traces=traces[:7],
                contrast_interval=[float(lower), float(upper)],
                lower_rational=str(lower), upper_rational=str(upper),
                tail_upper_rational=str(tail), tail_upper_float=float(tail),
                certified_greater_than='7/500')


def contract_example():
    """Small abstract contract check; NOT a Heisenberg simulator construction."""
    h = .4*core.PAULI[0]+.3*core.PAULI[2]
    theta, epsilon = .03, .02
    v = np.vstack([np.eye(2), np.zeros((2, 2))]).astype(complex)
    rotation = np.kron(np.array([[math.cos(theta), -math.sin(theta)],
                                [math.sin(theta), math.cos(theta)]]), np.eye(2))
    w = rotation@v
    diagonal = np.zeros((4, 4), complex)
    diagonal[:2, :2] = h+epsilon*core.PAULI[2]
    diagonal[2:, 2:] = 10*np.eye(2)
    hs = rotation@diagonal@rotation.conj().T
    return h, hs, v, w, norm(w-v), epsilon


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=5e-11):
        self.assertLess(norm(a-b), tolerance)

    def test_01_whole_network_is_eligible(self):
        h = network.model()[0]
        terms = target_terms()
        rebuilt = sum(float(c)*matrix(word) for word, c in terms.items())
        self.close(rebuilt, h)
        self.close(h, h.conj().T)
        self.assertEqual(float(np.max(np.abs(h.imag))), 0.)
        self.assertEqual(len(terms), 25)
        self.assertEqual(sum('Y' in word for word in terms), 6)
        self.assertEqual(max(sum(letter != 'I' for letter in word) for word in terms), 3)
        self.assertTrue(all(word.count('Y') % 2 == 0 for word in terms))
        heis = sum(np.kron(p, p) for p in core.PAULI)
        self.close(heis, 2*core.swap(2, 0, 1)-np.eye(4))
        OBS['target_contract'] = dict(original_system_qubits=6, dimension=64,
            pauli_terms=25, terms_with_two_Y=6, maximum_term_locality=3,
            real_symmetric=True, operator_norm=norm(h), analytic_norm_upper=9,
            full_shared_data_and_three_relation_registers_included=True,
            all_complex_input_states_still_allowed=True,
            exact_Pauli_manifest={word: str(c) for word, c in terms.items()},
            Heisenberg_equals_twice_SWAP_minus_identity=True)

    def test_02_simultaneous_six_mediator_identity(self):
        layer = FirstLayer()
        no_y = sum(float(c)*matrix(word).real for word, c in target_terms().items() if 'Y' not in word)
        counter = sum(float(abs(c)) for _, c in layer.y_terms)
        correction = sum(k@k for k in layer.k)
        effective = no_y+counter*np.eye(64)-correction
        self.close(effective, layer.h)
        local_errors = []
        for (word, coefficient), k in zip(layer.y_terms, layer.k):
            error = norm(-k@k+abs(float(coefficient))*np.eye(64)-float(coefficient)*matrix(word))
            self.assertLess(error, 1e-12)
            local_errors.append(error)
        image = layer.a(layer.v)
        self.close(layer.h0(image), image)
        self.close(layer.v.T@layer.b(layer.v)-image.T@image, layer.h)
        self.assertTrue(all(set(word) <= set('IXZ') for word in [*layer.a_terms, *layer.b_terms]))
        commutator = max(norm(a@b-b@a) for a, b in itertools.combinations(layer.k, 2))
        self.assertGreater(commutator, .2)
        OBS['explicit_first_layer'] = dict(original_qubits=6, mediator_qubits=6,
            intermediate_qubits=12, first_layer_only_not_final_Heisenberg_size=True,
            intermediate_dimension=4096, largest_Pauli_support=4,
            mediator_ground_state='|000000>', scalar_counterterm=counter,
            maximum_single_term_error=max(local_errors), whole_effective_error=norm(effective-layer.h),
            maximum_shared_data_K_commutator=commutator,
            intermediate_A_terms=layer.a_terms, intermediate_B_terms=layer.b_terms,
            noncommuting_data_terms_run_simultaneously=True,
            full_final_Heisenberg_coupling_list_generated=False)

    def test_03_first_layer_unknown_input_error_bound(self):
        layer = FirstLayer()
        v, h = layer.v, layer.h
        w1 = -layer.inverse(layer.a(v))
        w2 = -layer.inverse(layer.a(w1))
        self.close(layer.h0(w1)+layer.a(v), np.zeros_like(v))
        self.close(layer.h0(w2)+layer.a(w1)+layer.b(v), v@h)
        a3 = layer.a(w2)+layer.b(w1)-w1@h
        a4 = layer.b(w2)-w2@h
        constants = [('w1', w1, math.sqrt(3)), ('w2', w2, 3*math.sqrt(6)),
                     ('a3', a3, 39*math.sqrt(3)), ('a4', a4, 63*math.sqrt(6))]
        for _, value, upper in constants:
            self.assertLessEqual(norm(value), upper+1e-10)
        mu = 1/13
        w = v+mu*w1+mu*mu*w2
        residual = layer.h0(w)+mu*layer.a(w)+mu*mu*layer.b(w)-w@(mu*mu*h)
        error = norm(residual-(mu**3*a3+mu**4*a4))
        self.assertLess(error, 1e-11)
        OBS['first_layer_uniform_bound'] = dict(
            Hamiltonian='Delta*H0 + sqrt(Delta)*A + B; mu=1/sqrt(Delta)',
            operator_error_bound='sqrt(3)*(2+39|t|)*mu + sqrt(6)*(6+63|t|)*mu^2',
            constants=[dict(name=name, numerical_norm=norm(value), analytic_upper=upper) for name, value, upper in constants],
            full_code_columns_checked=64, exact_residual_identity_error=error,
            witness_mu='1/100000', witness_Delta=10**10,
            bound_at_time_501_over_1000=first_layer_bound(1e-5, .501),
            all_inputs_and_arbitrary_reference_included=True,
            large_energy_scale_recorded_as_cost=True,
            full_spectral_simulator_run_claimed=False)

    def test_04_certified_signal_survives_static_simulation(self):
        certificate = exact_target_certificate()
        self.assertEqual(certificate['scaled_low_order_traces'][:5], [0, 0, 0, 0, 24])
        h, swaps, occ, blank = network.model()
        rho = np.kron((np.eye(8)-swaps[0])/4, blank)
        o = np.kron(np.eye(8), occ[0]-occ[1])
        u = core.evolve(h, .5)
        gap = float(np.trace(o@u@rho@u.conj().T).real)
        self.assertLess(abs(gap-sum(certificate['contrast_interval'])/2), 2e-13)
        self.assertAlmostEqual(norm(1j*(h@o-o@h)), 2)
        eta, epsilon, horizon = Fraction(1, 10000), Fraction(1, 10000), Fraction(501, 1000)
        b = 2*eta+horizon*epsilon
        window_target_gap = Fraction(7, 500)-2*Fraction(1, 1000)
        simulator_gap = window_target_gap-2*b
        self.assertGreater(simulator_gap, Fraction(114, 10000))
        OBS['whole_network_signal'] = dict(certificate=certificate, spectral_cross_check=gap,
            response_observable='n_12 - n_13', derivative_norm=2,
            window=['499/1000', '501/1000'], target_window_lower=str(window_target_gap),
            static_simulation_contract=dict(eta=str(eta), epsilon=str(epsilon), horizon=str(horizon)),
            reference_complete_trace_distance_upper=str(b),
            simulator_window_gap_lower=str(simulator_gap), simulator_window_gap_lower_float=float(simulator_gap),
            applies_to_any_simulator_meeting_the_full_contract=True,
            final_Heisenberg_hardware_numerically_simulated=False)

    def test_05_reference_complete_contract_check(self):
        h, hs, v, w, eta, epsilon = contract_example()
        self.close(w.conj().T@w, np.eye(2))
        self.close(hs@w, w@(w.conj().T@hs@w))
        self.assertAlmostEqual(norm(w.conj().T@hs@w-h), epsilon)
        time = .7
        actual_map = core.evolve(hs, time)@v
        ideal_map = v@core.evolve(h, time)
        bound = 2*eta+time*epsilon
        self.assertLess(norm(actual_map-ideal_map), bound)
        psi = np.diag([1., 1j])/math.sqrt(2)
        actual, ideal = actual_map@psi, ideal_map@psi
        trace_distance = math.sqrt(max(0., 1-abs(np.vdot(actual, ideal))**2))
        self.assertLess(trace_distance, bound)
        self.close(actual.conj().T@actual, psi.conj().T@psi)
        OBS['reference_contract_check'] = dict(test_type='abstract 4D contract example, not Heisenberg construction',
            eta=eta, epsilon=epsilon, time=time, operator_error=norm(actual_map-ideal_map),
            entangled_complex_input_trace_distance=trace_distance, analytic_upper=bound,
            encoding_is_complex_linear_isometry=True,
            anti_linear_state_conjugation_or_input_estimation_used=False)

    def test_06_readout_keeps_leakage_and_all_outcomes(self):
        h, hs, v, w, eta, epsilon = contract_example()
        psi = np.array([1., 1j])/math.sqrt(2)
        actual = core.evolve(hs, .7)@v@psi
        ideal = core.evolve(h, .7)@psi
        p = v@v.conj().T
        effects = [np.outer(v[:, j], v[:, j].conj()) for j in range(2)]+[np.eye(4)-p]
        self.close(sum(effects), np.eye(4))
        probabilities = [float(np.vdot(actual, effect@actual).real) for effect in effects]
        self.assertGreater(probabilities[2], 1e-5)
        self.assertAlmostEqual(sum(probabilities), 1.)
        self.assertLess(abs(probabilities[1]-abs(ideal[1])**2), 2*eta+.7*epsilon)
        OBS['readout_contract'] = dict(test_type='same abstract contract example',
            complete_outcome_probabilities=probabilities,
            leakage_retained_as_outcome=True, probability_renormalization_or_postselection=False,
            induced_effects_and_initial_encoded_preparation_are_resources=True,
            arbitrary_future_interventions_automatically_implemented=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=436, baseline_round=435, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(full_finite_network_static_exchange_simulation_exists_by_cited_theorems=True,
            explicit_first_no_Y_layer_checked_on_all_64_input_columns=True,
            complete_unknown_reference_dynamical_contract_proved=True,
            target_signal_and_simulation_budget_certified=True,
            known_universal_Hamiltonian_theorems_not_claimed_as_new=True,
            final_Heisenberg_coupling_list_or_hardware_run_completed=False,
            target_dependent_resources_and_nonuniform_weights_still_inputs=True,
            exact_finite_resource_universal_implementation_claimed=False,
            natural_architecture_or_actual_positions_generated=False,
            three_dimensional_space_unconditionally_derived=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    else:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
