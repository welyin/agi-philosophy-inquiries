"""Round 434: a continuous conditional response from pair exchanges alone.

Four-qubit singlet codes and effective interactions are established tools.
This explicit finite construction audits a missing primitive of round 433;
the block partition, coupling hierarchy and prepared code states are inputs.
"""
import argparse
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import exchange_relation_audit as core

TARGET = Path(__file__).with_name('encoded_exchange_response_audit_results.json')
OBS = {}


def norm(a):
    return float(np.linalg.norm(a, 2))


def code_block():
    basis = np.eye(16)
    zero = (basis[5]-basis[6]-basis[9]+basis[10])/2
    one = (2*basis[3]+2*basis[12]-basis[5]-basis[6]-basis[9]-basis[10])/math.sqrt(12)
    return np.stack([zero, one], axis=1).astype(complex)


def model():
    block = code_block()
    code = np.kron(block, block)
    c4 = sum(core.swap(4, a, b) for a, b in itertools.combinations(range(4), 2))
    h0 = np.kron(c4, np.eye(16))+np.kron(np.eye(16), c4)
    za = -core.swap(8, 0, 1)
    zb = -core.swap(8, 4, 5)
    xb = (core.swap(8, 5, 6)-core.swap(8, 4, 6))/math.sqrt(3)
    t = core.swap(8, 0, 4)+core.swap(8, 1, 5)-np.eye(256)
    local = xb/6+za/12+zb/4
    # Spectral inverse on the complement of the four-dimensional code.
    energies, vectors = np.linalg.eigh(h0)
    positive = energies > 1e-9
    inverse = (vectors[:, positive]/energies[positive])@vectors[:, positive].conj().T
    w1 = -inverse@t@code
    w2 = -inverse@t@w1
    k2 = code.conj().T@t@w1
    logical_local = code.conj().T@local@code
    effective = k2+logical_local
    x, _, z = core.PAULI
    expected_k2 = (-5*np.eye(4)/12-(np.kron(z, np.eye(2))+np.kron(np.eye(2), z))/12
                   -np.kron(z, z)/6)
    expected = -5*np.eye(4)/12+(np.kron(np.eye(2), x)+np.kron(np.eye(2), z)-np.kron(z, z))/6
    return dict(block=block, code=code, c4=c4, h0=h0, t=t, local=local,
                inverse=inverse, w1=w1, w2=w2, k2=k2, effective=effective,
                expected_k2=expected_k2, expected=expected, za=za, zb=zb, xb=xb)


def evolution_columns(h, columns, time):
    values, vectors = np.linalg.eigh(h)
    return vectors@(np.exp(-1j*time*values)[:, None]*(vectors.conj().T@columns))


def error_bound(epsilon, tau):
    e, s = abs(epsilon), abs(tau)
    return math.sqrt(3)*((.5+s)*e+(3/8+21*s/64)*e*e)


class Audit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = model()

    def close(self, a, b, tolerance=1e-10):
        self.assertLess(norm(a-b), tolerance)

    def test_01_singlet_code_and_exact_internal_directions(self):
        m, x, z = self.m, core.PAULI[0], core.PAULI[2]
        v, v4 = m['code'], m['block']
        self.close(v.conj().T@v, np.eye(4))
        self.close(m['h0']@v, np.zeros_like(v))
        eigenvalues = np.linalg.eigvalsh(m['c4'])
        distinct, counts = np.unique(np.rint(eigenvalues).astype(int), return_counts=True)
        self.assertEqual(list(zip(distinct.tolist(), counts.tolist())), [(0, 2), (2, 9), (6, 5)])
        for name, logical in [('za', np.kron(z, np.eye(2))),
                              ('zb', np.kron(np.eye(2), z)),
                              ('xb', np.kron(np.eye(2), x))]:
            self.close(m[name]@v, v@logical)
            self.assertLessEqual(norm(m[name]), 1+1e-12)
        for pauli in core.PAULI:
            total = np.zeros((16, 16), complex)
            for j in range(4):
                term = np.array([[1.]])
                for k in range(4):
                    term = np.kron(term, pauli if j == k else np.eye(2))
                total += term
            self.close(total@v4, np.zeros_like(v4))
        OBS['code'] = dict(physical_qubits=8, logical_qubits=2,
            one_block_spectrum_and_multiplicity=[[0, 2], [2, 9], [6, 5]],
            one_block_gap=2, zero_total_spin_excitation_gap=4,
            code_dimension=4, physical_three_dimensional_space_assumed=False)

    def test_02_virtual_exchange_and_nonzero_conditional_term(self):
        m, v, t = self.m, self.m['code'], self.m['t']
        self.close(v.conj().T@t@v, np.zeros((4, 4)))
        self.close(m['h0']@t@v, 4*t@v)
        self.close(m['k2'], m['expected_k2'])
        self.close(v.conj().T@t@t@v, -4*m['expected_k2'])
        one_link = core.swap(8, 0, 4)-np.eye(256)/2
        one_link_k2 = -v.conj().T@one_link@m['inverse']@one_link@v
        self.close(one_link_k2, -3*np.eye(4)/16)
        self.assertAlmostEqual(norm(t@v), math.sqrt(3))
        self.assertAlmostEqual(norm(m['k2']), 3/4)
        OBS['second_order'] = dict(k2='-5 I/12 - (Z_D+Z_R)/12 - Z_D Z_R/6',
            conditional_ZZ_coefficient=-1/6, first_order_logical_interaction=0,
            one_cross_link_gives_only_scalar=-3/16,
            two_cross_link_paths_needed_in_this_construction=True,
            k2_operator_formula_error=norm(m['k2']-m['expected_k2']))

    def test_03_fixed_pair_exchange_hamiltonian(self):
        m, v = self.m, self.m['code']
        self.close(m['local']@v, v@(v.conj().T@m['local']@v))
        self.close(m['effective'], m['expected'])
        self.assertLessEqual(norm(m['local']), .5+1e-12)
        e = 1/512
        h = m['h0']+e*m['t']+e*e*m['local']
        weights = {pair: 1. for group in (range(4), range(4, 8))
                   for pair in itertools.combinations(group, 2)}
        weights[(0, 4)] = e
        weights[(1, 5)] = e
        weights[(0, 1)] -= e*e/12
        weights[(4, 5)] -= e*e/4
        weights[(5, 6)] += e*e/(6*math.sqrt(3))
        weights[(4, 6)] -= e*e/(6*math.sqrt(3))
        reconstructed = sum(weight*core.swap(8, *pair) for pair, weight in weights.items())-e*np.eye(256)
        self.close(h, reconstructed)
        self.assertGreater(min(weights.values()), 0)
        OBS['primitive_implementation'] = dict(total_nonzero_pair_couplings=len(weights),
            all_physical_pair_weights_positive_at_witness=True,
            time_independent=True, pulse_schedule_used=False,
            new_fundamental_flip_or_controlled_swap_used=False,
            effective_hamiltonian='-5 I/12 + (X_R + Z_R - Z_D Z_R)/6',
            weak_couplings_and_block_partition_are_prescribed=True)

    def test_04_uniform_dynamical_error_for_unknown_inputs(self):
        m, v = self.m, self.m['code']
        h0, t, local, k = (m[key] for key in ('h0', 't', 'local', 'effective'))
        w1, w2 = m['w1'], m['w2']
        self.close(h0@w1+t@v, np.zeros_like(v))
        self.close(h0@w2+t@w1+local@v, v@k)
        a3 = t@w2+local@w1-w1@k
        a4 = local@w2-w2@k
        self.assertLessEqual(norm(w1), math.sqrt(3)/4+1e-12)
        self.assertLessEqual(norm(w2), 3*math.sqrt(3)/16+1e-12)
        self.assertLessEqual(norm(a3), math.sqrt(3)+1e-12)
        self.assertLessEqual(norm(a4), 21*math.sqrt(3)/64+1e-12)
        # An exact finite polynomial identity, independent of spectral propagation.
        e = 1/7
        w = v+e*w1+e*e*w2
        h = h0+e*t+e*e*local
        self.close(h@w-w@(e*e*k), e**3*a3+e**4*a4)
        rows = []
        tau = 3*math.pi
        for denominator in (64, 128, 256, 512):
            e = 1/denominator
            physical_time = tau/(e*e)
            actual = evolution_columns(h0+e*t+e*e*local, v, physical_time)
            ideal = v@core.evolve(k, tau)
            error = norm(actual-ideal)
            leakage = norm(actual-v@(v.conj().T@actual))**2
            bound = error_bound(e, tau)
            self.assertLess(error, bound)
            rows.append(dict(epsilon_denominator=denominator, tau=tau,
                physical_time=physical_time, operator_error=error,
                analytic_operator_error_bound=bound, maximum_code_leakage_probability=leakage))
        OBS['uniform_error'] = dict(
            proof='polynomial intertwiner plus Duhamel; no normalized dressing or adiabatic switch required',
            bound='sqrt(3)*[(1/2+|tau|)*|epsilon| + (3/8+21|tau|/64)*epsilon^2]',
            covers_arbitrary_encoded_state_and_internal_reference=True,
            finite_tau_with_physical_time_proportional_to_epsilon_inverse_squared=True,
            numerical_scan_is_validation_not_proof=True, rows=rows)

    def test_05_activation_on_a_finite_time_window(self):
        m, v, e = self.m, self.m['code'], 1/512
        tau = 3*math.pi
        n = (np.eye(256)-m['zb'])/2
        n_logical = np.kron(np.eye(2), np.diag([0., 1.]))
        self.close(n@v, v@n_logical)
        actual = evolution_columns(m['h0']+e*m['t']+e*e*m['local'], v, tau/(e*e))
        probabilities = [float(np.vdot(actual[:, j], n@actual[:, j]).real) for j in (0, 2)]
        ideal = [math.sin(tau/6)**2, math.sin(math.sqrt(5)*tau/6)**2/5]
        bound = error_bound(e, tau)
        for p, p0 in zip(probabilities, ideal):
            self.assertLess(abs(p-p0), bound)
        window_bound = error_bound(e, tau+.25)
        contrast_lower_bound = math.cos(1/24)**2-.2-2*window_bound
        self.assertGreater(contrast_lower_bound, .72)
        self.assertGreater(probabilities[0]-probabilities[1], .97)
        # Ideal evolution entangles an initially independent logical |+>|0>.
        logical_input = np.array([1., 0., 1., 0.])/math.sqrt(2)
        ideal_output = core.evolve(m['effective'], tau)@logical_input
        reduced_d = ideal_output.reshape(2, 2)@ideal_output.reshape(2, 2).conj().T
        concurrence = 2*math.sqrt(max(0., float(np.linalg.det(reduced_d).real)))
        self.assertGreater(concurrence, .98)
        OBS['activation'] = dict(initial_relation_state='logical |0>',
            data_inputs='logical Z_D eigenvalues +1 and -1',
            ideal_probabilities_at_tau_3pi=ideal, actual_probabilities_at_epsilon_inverse_512=probabilities,
            actual_contrast=probabilities[0]-probabilities[1],
            guaranteed_tau_window=[tau-.25, tau+.25],
            analytic_window_contrast_lower_bound=contrast_lower_bound,
            analytic_window_equal_prior_success_lower_bound=(1+contrast_lower_bound)/2,
            ideal_output_concurrence_for_product_input=concurrence,
            exact_stopping_time_required=False,
            persistent_record_or_spatial_neighbor_generated=False)

    def test_06_linearity_reference_and_naive_projection_boundary(self):
        m, v, e = self.m, self.m['code'], 1/512
        tau = 3*math.pi
        h = m['h0']+e*m['t']+e*e*m['local']
        actual = evolution_columns(h, v, tau/(e*e))
        # Two logical inputs and a mixed preparation share exactly the same H.
        f = np.array([1., 0., 1j, 0.])/math.sqrt(2)
        g = np.array([0., 1., 1., 0.])/math.sqrt(2)
        rf, rg = np.outer(f, f.conj()), np.outer(g, g.conj())
        mixture_error = norm(actual@(.3*rf+.7*rg)@actual.conj().T-
                             (.3*np.outer(actual@f, (actual@f).conj())+.7*np.outer(actual@g, (actual@g).conj())))
        self.assertLess(mixture_error, 1e-12)
        # A two-dimensional reference entangled with unknown data; R is blank.
        initial = np.zeros((4, 2), complex)
        initial[0, 0] = initial[2, 1] = 1/math.sqrt(2)
        output = actual@initial
        reference_error = norm(output.conj().T@output-initial.conj().T@initial)
        recovered = evolution_columns(h, output, -tau/(e*e))
        inverse_error = norm(recovered-v@initial)
        self.assertLess(reference_error, 2e-12)
        self.assertLess(inverse_error, 3e-12)
        # PH P misses virtual excursions: its logical generator is only local.
        projected_generator = v.conj().T@m['local']@v
        projected = core.evolve(projected_generator, tau)
        n = np.kron(np.eye(2), np.diag([0., 1.]))
        projected_probabilities = [float(np.vdot(projected[:, j], n@projected[:, j]).real) for j in (0, 2)]
        self.assertAlmostEqual(*projected_probabilities)
        OBS['quantum_consistency'] = dict(affinity_error=mixture_error,
            reference_marginal_error=reference_error, whole_inverse_recovery_error=inverse_error,
            naive_PHP_activation_contrast=projected_probabilities[0]-projected_probabilities[1],
            virtual_excursions_are_essential=True,
            unknown_state_preserved_globally_not_unchanged_data_marginal=True)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=434, baseline_round=433, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(continuous_exchange_only_conditional_response_constructed=True,
            uniform_finite_time_error_and_reference_contract_proved=True,
            known_singlet_code_and_effective_interaction_tools_attributed=True,
            new_fundamental_flip_or_controlled_exchange_required_in_this_example=False,
            exact_raw_433_hamiltonian_implemented=False,
            complete_433_three_edge_model_implemented=False,
            partition_couplings_gap_and_code_preparation_are_inputs=True,
            derived_from_429_exchange_alone=False,
            full_cognitive_implementation_completed=False,
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
