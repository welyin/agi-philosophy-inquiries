"""Round 397: interaction supports from calibrated short-time responses.

Known local matrix algebras and a fixed closed finite Hamiltonian are inputs.
No graph, geometric coordinates, or tensor-Markov channel decomposition is
supplied to the estimator. Small examples use exhaustive product tomography;
this is not a scalable or autonomous discovery claim.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('instantaneous_coupling_audit_results.json')
PAULI = (np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex), np.diag([1, -1]).astype(complex))


def tensor(items):
    result = np.ones((1, 1), complex)
    for item in items:
        result = np.kron(result, item)
    return result


@lru_cache(None)
def word(labels):
    return tensor(PAULI[label] for label in labels)


def local(n, site, axis):
    return word(tuple(axis if i == site else 0 for i in range(n)))


def comm(a, b):
    return a@b-b@a


def hs(a):
    return float(np.linalg.norm(a)/np.sqrt(a.shape[0]))


def unitary(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values))@vectors.conj().T


def subsets(n):
    return [tuple(i for i in range(n) if mask & (1 << i)) for mask in range(1, 2**n)]


def exact_derivatives(h, n):
    return {(i, axis): 1j*comm(h, local(n, i, axis))
            for i in range(n) for axis in (1, 2, 3)}


def cumulative_weights(derivatives, n):
    # Input consists solely of local Heisenberg derivatives, not H or its graph.
    result = {}
    for sites in subsets(n):
        value = 0.
        for axes in itertools.product((1, 2, 3), repeat=len(sites)):
            operator = derivatives[(sites[0], axes[0])]
            for site, axis in zip(sites[1:], axes[1:]):
                operator = comm(operator, local(n, site, axis))
            value += hs(operator)**2
        result[sites] = value/(8**len(sites))
    return result


def exact_support_weights(cumulative):
    return {sites: sum((-1)**(len(more)-len(sites))*value
                       for more, value in cumulative.items() if set(sites) <= set(more))
            for sites in cumulative}


def pauli_support_weights(h, n):
    result = {sites: 0. for sites in subsets(n)}
    for labels in itertools.product(range(4), repeat=n):
        sites = tuple(i for i, label in enumerate(labels) if label)
        if sites:
            coefficient = np.trace(word(labels)@h)/(2**n)
            result[sites] += float(abs(coefficient)**2)
    return result


def canonical_audit():
    rng = np.random.default_rng(397)
    cases = []
    for n in (2, 3, 4):
        raw = rng.normal(size=(2**n, 2**n))+1j*rng.normal(size=(2**n, 2**n))
        h = (raw+raw.conj().T)/(2*np.sqrt(2**n))
        truth = pauli_support_weights(h, n)
        cumulative = cumulative_weights(exact_derivatives(h, n), n)
        recovered = exact_support_weights(cumulative)
        error = max(abs(truth[s]-recovered[s]) for s in truth)
        identity_error = max(abs(cumulative[s]-sum(v for t, v in truth.items() if set(s) <= set(t)))
                             for s in truth)
        cases.append(dict(qubits=n, nonempty_supports=2**n-1,
                          support_reconstruction_error=error, commutator_identity_error=identity_error))
    return cases


def frame_audit():
    n = 3
    h = .7*word((3, 3, 0))+.5*word((0, 1, 1))+.3*word((1, 2, 3))+.2*word((2, 0, 0))
    rotations = [unitary(.17*PAULI[1]+(.2+.13*i)*PAULI[2]-.31*PAULI[3], 1.) for i in range(n)]
    u = tensor(rotations)
    before = cumulative_weights(exact_derivatives(h, n), n)
    after = cumulative_weights(exact_derivatives(u@h@u.conj().T, n), n)
    return dict(independent_local_frame_error=max(abs(before[s]-after[s]) for s in before),
                common_spatial_axis_assumed=False)


def finite_pulse_audit():
    n = 3
    h = .7*word((3, 3, 0))+.31*word((1, 2, 3))
    a, b = local(n, 0, 1), local(n, 1, 2)
    rho = tensor((PAULI[0]+PAULI[axis])/2 for axis in (2, 1, 3))
    target = float(np.trace(rho@comm(comm(h, a), b)).real)
    pplus, pminus = unitary(b, np.pi/4), unitary(b, -np.pi/4)
    cases = []
    bound_m = 1.01
    for t in (.01, .005, .0025):
        u = unitary(h, t)
        evolved_a = u.conj().T@a@u
        fplus = np.trace(rho@pplus.conj().T@evolved_a@pplus).real
        fminus = np.trace(rho@pminus.conj().T@evolved_a@pminus).real
        estimator = float((fplus-fminus)/t)
        operator_estimator = comm(evolved_a, b)/(1j*t)
        cases.append(dict(time=t, measured_slope=estimator, exact_derivative=target,
                          finite_pulse_identity_error=float(abs(estimator-np.trace(rho@operator_estimator).real)),
                          truncation_error=abs(estimator-target), proven_bias_bound=4*bound_m**2*t))
    return cases


def hyperedge_audit():
    strength = .6
    triple = strength*word((3, 3, 3))
    triangle = strength*(word((3, 3, 0))+word((3, 0, 3))+word((0, 3, 3)))
    q_triple = cumulative_weights(exact_derivatives(triple, 3), 3)
    q_triangle = cumulative_weights(exact_derivatives(triangle, 3), 3)
    return dict(pair_weights_triple=[q_triple[s] for s in ((0, 1), (0, 2), (1, 2))],
                pair_weights_triangle=[q_triangle[s] for s in ((0, 1), (0, 2), (1, 2))],
                triple_support_weight_triple=q_triple[(0, 1, 2)],
                triple_support_weight_triangle=q_triangle[(0, 1, 2)],
                exact_pair_support_in_triple=exact_support_weights(q_triple)[(0, 2)])


def relay_audit():
    g, k = .7, .9
    h = g*word((1, 1, 0))+k*word((0, 3, 3))
    a, b = word((3, 0, 0)), word((0, 0, 1))
    first = comm(comm(h, a), b)
    second = comm(comm(h, comm(h, a)), b)
    q = cumulative_weights(exact_derivatives(h, 3), 3)
    cases = []
    r = np.sqrt(g*g+k*k)
    for t in (.04, .02, .01):
        u = unitary(h, t)
        value = hs(comm(u.conj().T@a@u, b))
        predicted = 4*g*k*np.sin(r*t)**2/(r*r)
        cases.append(dict(time=t, measured_commutator=value, exact_formula=float(predicted),
                          quadratic_coefficient=value/(t*t)))
    return dict(pair_strengths={','.join(map(str, s)): float(np.sqrt(q[s]))
                                for s in ((0, 1), (0, 2), (1, 2))},
                first_order_ac=hs(first), second_order_coefficient_error=hs(second-8j*g*k*word((2, 2, 2))),
                expected_quadratic_limit=4*g*k, cases=cases)


def condition_middle(operator, state):
    # Indices are A,M,C,A',M',C'; Tr_M(state_M operator).
    return np.einsum('abcdef,eb->acdf', operator.reshape((2,)*6), state).reshape(4, 4)


def hidden_memory_audit():
    h = .6*word((3, 3, 3))
    mixed_sum = prepared_sum = 0.
    for alpha, beta in itertools.product((1, 2, 3), repeat=2):
        d = comm(comm(h, local(3, 0, alpha)), local(3, 2, beta))
        mixed_sum += hs(condition_middle(d, PAULI[0]/2))**2
        prepared_sum += hs(condition_middle(d, (PAULI[0]+PAULI[3])/2))**2
    return dict(full_ac_strength=.6, frozen_mixed_memory_strength=float(np.sqrt(mixed_sum/64)),
                prepared_memory_strength=float(np.sqrt(prepared_sum/64)),
                unobserved_internal_memory_can_hide_first_order_signal=True)


def product_tomography(h, t, n, outcome_error):
    rng = np.random.default_rng(2397)
    u = unitary(h, t)
    # Store actual product-preparation / local-readout expectation data.
    data = {}
    axes_list = list(itertools.product((1, 2, 3), repeat=n))
    signs_list = list(itertools.product((-1, 1), repeat=n))
    observables = {(i, axis): local(n, i, axis) for i in range(n) for axis in (1, 2, 3)}
    for axes in axes_list:
        for signs in signs_list:
            rho = tensor((PAULI[0]+sign*PAULI[axis])/2 for axis, sign in zip(axes, signs))
            evolved = u@rho@u.conj().T
            for key, observable in observables.items():
                value = float(np.trace(evolved@observable).real)+rng.uniform(-outcome_error, outcome_error)
                data[(axes, signs, key)] = float(np.clip(value, -1, 1))
    derivatives = {}
    maximum_error = 0.
    for key, observable in observables.items():
        estimate = np.zeros_like(observable)
        for labels in itertools.product(range(4), repeat=n):
            axes = tuple(label or 3 for label in labels)
            coefficient = sum(np.prod([sign for label, sign in zip(labels, signs) if label])
                              *data[(axes, signs, key)] for signs in signs_list)/(2**n)
            estimate += coefficient*word(labels)
        maximum_error = max(maximum_error, hs(estimate-u.conj().T@observable@u))
        derivatives[key] = (estimate-observable)/t
    return derivatives, dict(unique_product_preparations=6**n,
                             local_expectation_settings=len(data),
                             maximum_heisenberg_hs_error=maximum_error,
                             proven_heisenberg_hs_error=2**n*outcome_error,
                             quantum_shot_sampling_performed=False)


def finite_error_audit():
    n = 3
    h = .7*word((1, 1, 0))+.9*word((0, 3, 3))
    m, epsilon = 1.6, 1e-7
    eta = 2**n*epsilon
    t = np.sqrt(eta/(2*m*m))
    derivatives, tomography = product_tomography(h, t, n, epsilon)
    estimated = cumulative_weights(derivatives, n)
    truth = cumulative_weights(exact_derivatives(h, n), n)
    radius = .75*(2*m*m*t+eta/t)
    pairs = []
    for s in ((0, 1), (0, 2), (1, 2)):
        measured, actual = np.sqrt(estimated[s]), np.sqrt(truth[s])
        pairs.append(dict(sites=list(s), measured_strength=float(measured), true_strength=float(actual),
                          absolute_error=float(abs(measured-actual)), certified_present=bool(measured>radius)))
    return dict(norm_bound=m, single_expectation_error=epsilon, chosen_time=float(t),
                pair_strength_error_radius=float(radius), tomography=tomography, pairs=pairs,
                recovery_uses_nonzero_gap_promise=True,
                small_signal_certifies_exact_absence_without_gap=False)


@lru_cache(None)
def report():
    return dict(round=397,
                scope='Conditional recovery of intrinsic Hamiltonian supports from a known finite tensor interface and calibrated short-time data; no input interaction graph or disjoint-channel tensor-Markov assumption. Interface access, closed fixed H, repeatable preparation and timing remain inputs; no spatial dimension or autonomous cognitive completion is derived.',
                canonical_cases=canonical_audit(), local_frames=frame_audit(),
                finite_pulses=finite_pulse_audit(), hyperedge=hyperedge_audit(),
                coherent_relay=relay_audit(), hidden_memory=hidden_memory_audit(),
                finite_error=finite_error_audit(),
                input_interaction_graph_required=False,
                disjoint_channel_tensor_markov_required=False,
                closed_hamiltonian_time_evolution_assumed=True,
                tensor_interface_and_repeatable_calibration_assumed=True,
                recovered_couplings_identified_with_spatial_contact=False,
                full_position_generation_completed=False)


class Audit(unittest.TestCase):
    def test_canonical_supports_recovered_from_commutators(self):
        for case in report()['canonical_cases']:
            self.assertLess(case['support_reconstruction_error'], 2e-12)
            self.assertLess(case['commutator_identity_error'], 2e-12)

    def test_independent_local_frame_invariance(self):
        self.assertLess(report()['local_frames']['independent_local_frame_error'], 2e-12)

    def test_two_finite_control_pulses_realize_response_witness(self):
        for case in report()['finite_pulses']:
            self.assertAlmostEqual(case['exact_derivative'], 2.8)
            self.assertLess(case['finite_pulse_identity_error'], 2e-12)
            self.assertLessEqual(case['truncation_error'], case['proven_bias_bound'])

    def test_pair_graph_does_not_replace_many_body_support(self):
        case = report()['hyperedge']
        self.assertTrue(np.allclose(case['pair_weights_triple'], case['pair_weights_triangle']))
        self.assertAlmostEqual(case['triple_support_weight_triple'], .36)
        self.assertAlmostEqual(case['triple_support_weight_triangle'], 0)
        self.assertAlmostEqual(case['exact_pair_support_in_triple'], 0)

    def test_coherent_relay_is_second_order_without_direct_coupling(self):
        case = report()['coherent_relay']
        self.assertAlmostEqual(case['pair_strengths']['0,2'], 0)
        self.assertAlmostEqual(case['first_order_ac'], 0)
        self.assertLess(case['second_order_coefficient_error'], 2e-12)
        for point in case['cases']:
            self.assertAlmostEqual(point['measured_commutator'], point['exact_formula'])
        self.assertLess(abs(case['cases'][-1]['quadratic_coefficient']-case['expected_quadratic_limit']), .001)

    def test_uncontrolled_internal_memory_hides_coupling(self):
        case = report()['hidden_memory']
        self.assertAlmostEqual(case['frozen_mixed_memory_strength'], 0)
        self.assertAlmostEqual(case['prepared_memory_strength'], case['full_ac_strength'])

    def test_product_measurements_and_finite_error_recovery(self):
        case = report()['finite_error']
        data = case['tomography']
        self.assertLessEqual(data['maximum_heisenberg_hs_error'], data['proven_heisenberg_hs_error'])
        for pair in case['pairs']:
            self.assertLessEqual(pair['absolute_error'], case['pair_strength_error_radius'])
            self.assertEqual(pair['certified_present'], pair['true_strength']>0)
            if pair['true_strength']>0:
                self.assertGreater(pair['true_strength'], 2*case['pair_strength_error_radius'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result['checks'] = dict(run=checked.testsRun, failures=len(checked.failures), errors=len(checked.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
