"""Round 447: internal swaps and a recorded-partner conditional phase.
The reciprocal penalty and the matching-dependent logical encoding are inputs.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import io
import itertools as it
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import reciprocal_exchange_dynamics_audit as old

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'internal_partner_response_audit_results.json'
OBS = {}
S_LOCAL = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=np.int64)
X_LOCAL = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=np.int64)


def matrix_unit(size, row, col):
    a = np.zeros((size, size), dtype=np.int64)
    a[row, col] = 1
    return a


def physical_word(matching, local):
    # -1 is the shared blank; a marker always survives in at least one register.
    return tuple(v for j, state in zip(matching, local)
                 for v in ((j, j) if state == 0 else
                           ((j, -1) if state == 1 else (-1, j))))


def raw_energy(word):
    address = word[::2]
    return len(address)-2*sum(address[i] == j and address[j] == i
                              for i, j in it.combinations(range(len(address)), 2))


def build(n):
    matches = old.matching_words(n)
    local = tuple(it.product(range(3), repeat=n))
    words = tuple(physical_word(m, q) for m in matches for q in local)
    lookup = {w: i for i, w in enumerate(words)}
    assert len(lookup) == len(words)
    v = np.zeros((len(words), len(words)), dtype=np.int64)
    for col, w in enumerate(words):
        for site in range(n):
            nxt = old.swap(w, (2*site, 2*site+1))
            v[lookup[nxt], col] += 1
    h0 = np.array([raw_energy(w) for w in words], dtype=np.int64)
    bits = tuple(it.product(range(2), repeat=n))
    z = np.zeros((len(words), len(matches)*len(bits)), dtype=np.int64)
    for k, m in enumerate(matches):
        for b, q in enumerate(bits):
            z[lookup[physical_word(m, q)], k*len(bits)+b] = 1
    return matches, local, words, h0, v, z


def energies(epsilon):
    x = float(epsilon)
    e_single = -x*x/(1+np.sqrt(1+x*x))
    h_double = np.array([[0, np.sqrt(2)*x, 0],
                         [np.sqrt(2)*x, 2, np.sqrt(2)*x],
                         [0, np.sqrt(2)*x, 2]])
    eigenvalues, eigenvectors = np.linalg.eigh(h_double)
    e11 = eigenvalues[0]
    return np.array([2*x, x+e_single, x+e_single, e11]), e11-2*e_single, eigenvectors[:, 0]


def logical_energies(n, epsilon):
    pair, _, _ = energies(epsilon)
    return np.array([sum(pair[2*bits[i]+bits[j]]
                        for i, j in enumerate(m) if i < j)
                     for m in old.matching_words(n)
                     for bits in it.product(range(2), repeat=n)])


def reduced(psi, site, local_dimension, n, reference=1):
    tensor = psi.reshape((local_dimension,)*n+(reference,))
    other = [i for i in range(n) if i != site]
    matrix = tensor.transpose([site, n]+other).reshape(local_dimension*reference, -1)
    return matrix@matrix.conj().T


def trace_distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def poly_mul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out[i+j] = out.get(i+j, F(0))+x*y
    return out


class Audit(unittest.TestCase):
    def test_01_raw_invariant_marker_code_and_pair_factorization(self):
        reports = []
        for n in (2, 4):
            matches, local, words, h0, v, z = build(n)
            block = 3**n
            for k, m in enumerate(matches):
                expected = [2*sum(q[i] == 2 or q[j] == 2
                                  for i, j in enumerate(m) if i < j) for q in local]
                np.testing.assert_array_equal(h0[k*block:(k+1)*block], expected)
                self.assertEqual(np.count_nonzero(v[k*block:(k+1)*block, :k*block]), 0)
                self.assertEqual(np.count_nonzero(v[k*block:(k+1)*block, (k+1)*block:]), 0)
                for q in local:
                    w = physical_word(m, q)
                    recovered = tuple(next(value for value in w[2*i:2*i+2]
                                           if value != -1) for i in range(n))
                    self.assertEqual(recovered, m)
            np.testing.assert_array_equal(z.T@z, np.eye(z.shape[1], dtype=int))
            reports.append(dict(subjects=n, matching_sectors=len(matches),
                complete_invariant_dimension=len(words), unknown_logical_code_dimension=z.shape[1]))
        np.testing.assert_array_equal(S_LOCAL@S_LOCAL, np.eye(3, dtype=int))
        self.assertGreater(np.linalg.norm(S_LOCAL@X_LOCAL-X_LOCAL@S_LOCAL), 1)
        OBS['invariant_model'] = dict(exact_raw_swap_checks=reports,
            matching_recoverable_from_both_registers_at_all_times=True,
            marker_and_logical_state_not_identified=True,
            tag_carried_data_conservation_of_446_broken=True,
            conditional_logical_encoding_is_explicit_input=True)

    def test_02_exact_pair_spectrum_connected_shift(self):
        x = F(1, 16)
        levels, chi, ground = energies(x)
        self.assertLess(chi, 0)
        e = levels[-1]
        self.assertLess(abs(e**3-4*e**2+(4-4*float(x*x))*e+4*float(x*x)), 1e-14)
        # Exact series substitution e(z)=-z+z^3/4, z=epsilon^2.
        p = {1: F(-1), 3: F(1, 4)}
        p2, p3 = poly_mul(p, p), poly_mul(poly_mul(p, p), p)
        residual = {k: p3.get(k, 0)-4*p2.get(k, 0)+4*p.get(k, 0)
                    -4*p.get(k-1, 0)+(4 if k == 1 else 0) for k in range(1, 4)}
        self.assertEqual(residual, {1: 0, 2: 0, 3: 0})
        # h_double = h_one tensor I + I tensor h_one - 2|ee><ee|.
        single = np.array([[0, float(x)], [float(x), 2]])
        ev, vv = np.linalg.eigh(single)
        product_ground = np.kron(vv[:, 0], vv[:, 0])
        double = np.kron(single, np.eye(2))+np.kron(np.eye(2), single)
        double[3, 3] -= 2
        trial = float(product_ground@double@product_ground)
        self.assertLess(trial, 2*ev[0])
        self.assertLessEqual(e, trial+1e-14)
        _, _, _, h0, v, z = build(2)
        h = np.diag(h0)+float(x)*v
        # Each computational bit pattern is a separate occupancy sector.
        for col in range(4):
            indices = np.flatnonzero(np.abs(z[:, col])+np.abs(v@z[:, col])
                                     +np.abs(v@v@z[:, col]))
            self.assertAlmostEqual(float(np.linalg.eigvalsh(h[np.ix_(indices, indices)])[0]),
                                   levels[col], places=13)
        OBS['spectrum'] = dict(epsilon='1/16', dimensionless_logical_energies=list(map(old.short, levels)),
            dimensionless_connected_shift=old.short(chi),
            exact_connected_shift_series='-epsilon^4/4 + 3*epsilon^6/8 + O(epsilon^8)',
            strict_negative_shift_for_all_nonzero_lambda_proved_by_variation=True,
            fourth_order_truncation_used_for_long_time_propagation=False)

    def test_03_bare_code_all_time_bound(self):
        x = F(1, 16)
        levels, chi, ground = energies(x)
        sin_double = float(np.sqrt(np.sum(np.abs(ground[1:])**2)))
        resolvent_bound = np.sqrt(2)*float(x)/(2-np.sqrt(2)*float(x))
        self.assertLess(sin_double, resolvent_bound)
        self.assertLess(resolvent_bound, float(x))
        # sqrt(2)/(2-sqrt(2)/4)<1 follows from 50<64, for all 0<x<=1/4.
        self.assertLess(50, 64)
        _, _, _, h0, v, z = build(2)
        h = np.diag(h0)+float(x)*v
        reports = []
        for phase_offset in (-F(1, 4), F(0), F(1, 4)):
            time = (np.pi+float(phase_offset))/abs(chi)
            u = old.unitary(h, time)
            target = z*np.exp(-1j*time*levels)
            error = old.opnorm(u@z-target)
            leakage = old.opnorm(u@z-z@(z.T@u@z))**2
            self.assertLess(error, 2*float(x))
            self.assertLess(leakage, 4*float(x*x))
            reports.append(dict(phase_offset=str(phase_offset), dimensionless_time=old.short(time),
                all_input_operator_error=old.short(error),
                maximum_outside_code_probability=old.short(leakage)))
        OBS['all_time_bound'] = dict(maximum_double_ground_sine=old.short(sin_double),
            analytic_per_pair_operator_bound='2*abs(epsilon)',
            analytic_N_subject_operator_bound='N*abs(epsilon)',
            tensor_reference_uniform=True, exact_eigenvalue_phases_used=True,
            no_adiabatic_preparation_or_dressed_code_required=True,
            long_time_numerical_cross_checks=reports)

    def test_04_unknown_coherent_matching_logic_and_reference(self):
        n, x, time = 4, F(1, 16), 31
        matches, local, words, h0, v, z = build(n)
        u = old.unitary(np.diag(h0)+float(x)*v, time)
        el = logical_energies(n, x)
        error = old.opnorm(u@z-z*np.exp(-1j*time*el))
        self.assertLess(error, n*float(x))
        rng = np.random.default_rng(447)
        psi = rng.normal(size=(48, 3))+1j*rng.normal(size=(48, 3))
        psi /= np.linalg.norm(psi)
        actual = u@z@psi
        target = z@(np.exp(-1j*time*el)[:, None]*psi)
        distance = np.sqrt(max(0, 1-abs(np.vdot(actual, target))**2))
        self.assertLess(distance, n*float(x))
        # Full subject logical X is a sum over internal marker sectors; no graph measurement.
        readout_checks = 0
        for site in range(n):
            physical_x = np.zeros((len(words), len(words)), dtype=int)
            index = {w: i for i, w in enumerate(words)}
            for col, w in enumerate(words):
                a, d = w[2*site:2*site+2]
                if a == d and a >= 0:
                    nxt = list(w); nxt[2*site+1] = -1
                elif a >= 0 and d == -1:
                    nxt = list(w); nxt[2*site+1] = a
                else:
                    continue
                physical_x[index[tuple(nxt)], col] = 1
            expected = np.eye(1, dtype=int)
            for i in range(n):
                expected = np.kron(expected, np.array([[0, 1], [1, 0]]) if i == site else np.eye(2, dtype=int))
            np.testing.assert_array_equal(z.T@physical_x@z, np.kron(np.eye(len(matches), dtype=int), expected))
            readout_checks += 1
        OBS['unknown_interface'] = dict(coherent_matching_sectors=3, logical_qubits=4,
            reference_dimension=3, all_input_operator_error=old.short(error),
            sampled_complete_reference_trace_distance=old.short(distance),
            exact_local_readout_intertwiners=readout_checks,
            general_unknown_raw_N_plus_1_data_claimed=False,
            logical_encoding_is_isometry_without_cloning_unknown_states=True)

    def test_05_partner_signal_and_fixed_matching_no_cross_signal(self):
        x = F(1, 16)
        levels, chi, _ = energies(x)
        _, _, _, h0, v, z = build(2)
        h = np.diag(h0)+float(x)*v
        lower = F(1)-F(1, 128)-4*x
        self.assertEqual(lower, F(95, 128))
        reports = []
        for offset in (-F(1, 4), F(0), F(1, 4)):
            time = (np.pi+float(offset))/abs(chi)
            u = old.unitary(h, time)
            outputs = []
            for bit in (0, 1):
                logical = np.zeros(4); logical[2*bit:2*bit+2] = 1/np.sqrt(2)
                outputs.append(reduced(u@z@logical, 1, 3, 2))
            d = trace_distance(*outputs)
            self.assertGreater(d, float(lower))
            dephased = []
            for output in outputs:
                restricted = output.copy()
                restricted[0, 1:] = 0
                restricted[1:, 0] = 0
                dephased.append(restricted)
            restricted_d = trace_distance(*dephased)
            self.assertLess(restricted_d, 4*float(x))
            reports.append(dict(phase_offset=str(offset), actual_receiver_trace_distance=old.short(d),
                conserved_sector_readout_trace_distance=old.short(restricted_d)))
        # Fixed matching (01)(23) gives an exact tensor product of pair evolutions.
        matches, local, words, h04, v4, z4 = build(4)
        h_pair_0 = np.diag(h0)
        np.testing.assert_array_equal(np.diag(h04[:81]),
            np.kron(h_pair_0, np.eye(9, dtype=int))+np.kron(np.eye(9, dtype=int), h_pair_0))
        np.testing.assert_array_equal(v4[:81, :81],
            np.kron(v, np.eye(9, dtype=int))+np.kron(np.eye(9, dtype=int), v))
        OBS['physical_signal'] = dict(epsilon='1/16',
            exact_window='abs(chi)*t in [pi-1/4, pi+1/4]',
            strict_signal_lower=old.fraction_data(lower),
            all_time_per_pair_outside_code_upper=old.fraction_data(4*x*x),
            actual_receiver_includes_excited_record_state=True,
            no_postselection=True, window_samples=reports,
            fixed_matching_cross_pair_no_signalling_from_exact_tensor_factorization=True,
            interventions_must_preserve_marker_three_state_interface=True,
            matching_reconfiguration_and_this_processor_combined=False)
        OBS['physical_signal']['large_phase_signal_needs_cross_sector_preparation_and_readout'] = True
        OBS['physical_signal']['uniform_trace_distance_does_not_supply_fixed_phase_clock'] = True

    def test_06_unrestricted_raw_data_fourth_order_counterexample(self):
        m = (1, 0, 3, 2)
        other = (2, 3, 0, 1)
        initial = tuple(v for pair in zip(m, other) for v in pair)
        final_expected = tuple(v for pair in zip(other, m) for v in pair)
        count = Counter()
        coefficient = F(0)
        for order in it.permutations(range(4)):
            word = initial
            denominators = []
            for k, site in enumerate(order):
                word = old.swap(word, (2*site, 2*site+1))
                if k < 3:
                    denominators.append(raw_energy(word))
            self.assertEqual(word, final_expected)
            self.assertTrue(all(denominators))
            count[tuple(denominators)] += 1
            coefficient -= F(1, int(np.prod(denominators)))
        self.assertEqual(count, Counter({(2, 2, 2): 16, (2, 4, 2): 8}))
        self.assertEqual(coefficient, F(-5, 2))
        OBS['raw_data_boundary'] = dict(initial_address=list(m), initial_raw_data=list(other),
            final_address=list(other), final_raw_data=list(m),
            fourth_order_path_denominators={str(k): v for k, v in sorted(count.items())},
            leading_off_diagonal_coefficient_in_lambda4_over_Delta3=str(coefficient),
            unrestricted_data_would_reconfigure_matching_at_same_order_as_phase=True,
            blank_marker_code_needed_for_the_claimed_fixed_partner_interface=True)

    def test_07_single_marker_exact_logical_qubits(self):
        n = 4
        matches, local, words, h0, v, _ = build(n)
        selected = [j for j, q in enumerate(local) if all(value in (1, 2) for value in q)]
        indices = [81*k+j for k in range(len(matches)) for j in selected]
        embedding = np.eye(243, dtype=int)[:, indices]
        qh = embedding.T@(np.diag(h0)+v)@embedding
        np.testing.assert_array_equal((np.diag(h0)+v)@embedding, embedding@qh)
        expected = np.zeros((48, 48), dtype=int)
        bits = tuple(it.product(range(2), repeat=n))
        bit_index = {b: i for i, b in enumerate(bits)}
        for k, m in enumerate(matches):
            for col, b in enumerate(bits):
                j = k*16+col
                expected[j, j] = 2*sum(b[i] or b[partner] for i, partner in enumerate(m) if i < partner)
                for site in range(n):
                    nxt = list(b); nxt[site] = 1-nxt[site]
                    expected[k*16+bit_index[tuple(nxt)], j] += 1
        np.testing.assert_array_equal(qh, expected)
        # Address blank is a single fixed physical observable, covering unknown M.
        readouts = 0
        for site in range(n):
            raw = np.diag([int(w[2*site] == -1) for w in words])
            want = np.diag([b[site] for m in matches for b in bits])
            np.testing.assert_array_equal(raw@embedding, embedding@want)
            readouts += 1
        time = 0.7
        physical = old.unitary(np.diag(h0)+v, time)@embedding
        target = embedding@old.unitary(qh, time)
        err = old.opnorm(physical-target)
        self.assertLess(err, 1e-12)
        OBS['single_marker_exact_code'] = dict(subjects=n, coherent_matching_sectors=3,
            physical_invariant_dimension=48, logical_code_dimension=48,
            two_states=['|partner,blank>', '|blank,partner>'],
            exact_unknown_joint_graph_logic_reference_intertwining=True,
            graph_identity_recoverable_from_both_registers=True,
            instantaneous_address_registers_always_reciprocal=False,
            address_blank_readout_intertwiners=readouts,
            finite_time_full_operator_error=old.short(err),
            per_pair_exact_H='2*Delta*(I-|00><00|)+lambda*(X_i+X_j)',
            no_low_energy_or_weak_coupling_approximation=True,
            internal_swap_is_logical_X=True)

    def test_08_fixed_population_readout_signal_certificate(self):
        x = np.array([[0, 1], [1, 0]], dtype=int)
        h = np.diag([0, 2, 2, 2])+np.kron(x, np.eye(2, dtype=int))+np.kron(np.eye(2, dtype=int), x)
        effect = np.diag([0, 1, 0, 1])
        commutator = h@effect-effect@h
        np.testing.assert_array_equal(commutator.T@commutator, np.eye(4, dtype=int))
        # Rational Taylor certificate at t=1, with ||H||<=4 and exp(4)<81.
        real = np.zeros((4, 4), dtype=object)
        imag = np.zeros((4, 4), dtype=object)
        power = np.eye(4, dtype=object)
        for degree in range(33):
            factor = F(1, math.factorial(degree))*(-1 if degree % 4 in (1, 2) else 1)
            if degree % 2:
                imag += factor*power
            else:
                real += factor*power
            power = power@h
        probabilities = [sum(real[row, col]**2+imag[row, col]**2 for row in (1, 3))
                         for col in (0, 2)]
        tail = F(81*4**33, math.factorial(33))
        at_one = probabilities[1]-probabilities[0]-4*tail-2*tail**2
        window = at_one-F(1, 100)
        self.assertGreater(at_one, F(11, 100))
        self.assertGreater(window, F(1, 10))
        numerical = []
        for time in (F(199, 200), F(1), F(201, 200)):
            u = old.unitary(h, float(time))
            ps = [float(sum(abs(u[row, col])**2 for row in (1, 3))) for col in (0, 2)]
            gap = ps[1]-ps[0]
            self.assertGreater(gap, float(window))
            numerical.append(dict(time=str(time), probabilities=list(map(old.short, ps)),
                                  signed_fixed_effect_gap=old.short(gap)))
        OBS['fixed_readout_certificate'] = dict(Delta=1, lambda_value=1,
            initial_messages=['|00>', '|10>'],
            same_receiver_effect='address register is blank',
            rational_taylor_degree=32, exponential_operator_tail=old.fraction_data(tail),
            strict_at_t1_lower=old.fraction_data(at_one),
            interval=['199/200', '201/200'], derivative_absolute_bound=2,
            strict_window_lower=old.fraction_data(window), simple_lower='1/10',
            time_dependent_optimal_measurement_used=False,
            postselection_used=False, small_epsilon_assumption_used=False,
            preparation_and_actual_readout_apparatus_still_inputs=True,
            numerical_cross_checks=numerical)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=447, baseline_round=446, date='2026-09-24',
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(internal_substructure_swap_breaks_446_tagged_data_conservation=True,
            exact_matching_marker_code_and_pair_factorization=True,
            arbitrary_encoded_matching_logic_reference_covered=True,
            no_separate_cross_subject_data_exchange_added=True,
            strict_conditional_phase_for_all_nonzero_lambda=True,
            all_time_bare_code_bound_uses_exact_eigenvalue_phases=True,
            actual_recorded_partner_signal_and_fixed_sector_no_cross_signal=True,
            unrestricted_raw_data_fourth_order_scope_counterexample=True,
            single_marker_exact_qubit_code_and_record_selected_pair_H=True,
            same_fixed_population_effect_signal_above_one_tenth_on_window=True,
            low_energy_phase_readout_limitation_explicit=True,
            reciprocal_consistency_penalty_derived_from_429=False,
            initial_matching_dependent_encoding_generated_autonomously=False,
            arbitrary_raw_data_preserves_fixed_matching=False,
            dynamic_partner_reassignment_combined_with_this_processor=False,
            potential_contact_or_spatial_dimension_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as output:
            output.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))
