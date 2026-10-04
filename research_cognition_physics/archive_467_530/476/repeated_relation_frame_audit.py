"""Round 476: repeated relational preparations deform the reference shape.

Scientific baseline 474. Two fixed-H waits are separated by an explicitly
supplied internal storage swap and fresh independent data preparation. Tracing
old data computes the graph marginal; it is not physical information deletion.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import relational_direction_distance_response as previous
import sequential_tree_distance_readout_audit as reader

TARGET = Path(__file__).with_name('repeated_relation_frame_audit_results.json')
OBS = {}
RX = (1, 0, 0)
RY = (0, 1, 0)


def short(x):
    return float(f'{float(x):.11g}')


def apply(channel, rho, ref=1):
    return reader.apply_map(channel, rho, g=6, r=ref)


def graph_observables():
    trees = previous.system()[0]
    return np.array([[len(previous.old.path_between(tree, 6, a, b))-1 for tree in trees]
                     for a, b in itertools.combinations(previous.LEAVES, 2)], dtype=np.int64)


def channel_from_eta(eta, t):
    values, vectors = np.linalg.eigh(eta)
    assert values.min() > -2e-12
    mask = values > 1e-13
    columns = vectors[:, mask]*np.sqrt(values[mask])
    u = previous.unitary(t).reshape(64, 6, 64, 6)
    kraus = np.einsum('dgfi,fb->bdgi', u, columns, optimize=True)
    channel = np.einsum('bdgi,bdhj->ghij', kraus, kraus.conj(), optimize=True)
    return channel, kraus


@lru_cache(None)
def channel(r, t):
    return channel_from_eta(previous.data_state(r), t)


def repeated(t, u):
    first = apply(channel(RX, t)[0], np.ones((6, 6))/6)
    return apply(channel(RY, u)[0], first)


def integer_gaussian_dot(a, b):
    # Python integers avoid any product/sum overflow in the mixed coefficients.
    real, imag = 0, 0
    for i, j in itertools.product(range(6), repeat=2):
        ar, ai = int(a[i, j].real), int(a[i, j].imag)
        br, bi = int(b[j, i].real), int(b[j, i].imag)
        real += ar*br-ai*bi
        imag += ar*bi+ai*br
    return real, imag


@lru_cache(None)
def certificate():
    # A priori Gaussian-integer bounds, before any float-based recurrence.
    assert 512*18**8 < 2**53
    assert 6*128*18**8*2**8 < 2**53
    h = previous.system()[1]
    raw = previous.data_basis_raw()
    eta_x, eta_y = raw[0]+raw[1], raw[0]+raw[2]
    o = 2*graph_observables()[0]-5
    state = np.kron(eta_x, np.ones((6, 6), dtype=np.int64))
    observable = np.diag(np.tile(o, 64)).astype(complex)
    r, t = [], []
    for k in range(9):
        assert np.array_equal(state, np.round(state.real)+1j*np.round(state.imag))
        assert np.array_equal(observable, np.round(observable.real)+1j*np.round(observable.imag))
        assert np.max(np.abs(state.real)) <= 2*18**k
        assert np.max(np.abs(state.imag)) <= 2*18**k
        assert np.max(np.abs(observable.real)) <= 18**k
        assert np.max(np.abs(observable.imag)) <= 18**k
        r.append(state.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2))
        t.append(np.einsum('df,fgdh->gh', eta_y,
                          observable.reshape(64, 6, 64, 6), optimize=True))
        if k < 8:
            state = h@state-state@h
            observable = h@observable-observable@h
    mixed, coefficients, continuous = {}, [], []
    for n in range(9):
        total = Q(0)
        for k in range(n+1):
            ell = n-k
            real, imag = integer_gaussian_dot(t[ell], r[k])
            phase = (-1j)**k*(1j)**ell
            value = int(phase.real)*real-int(phase.imag)*imag
            assert int(phase.real)*imag+int(phase.imag)*real == 0
            coefficient = Q(value, 49152*math.factorial(k)*math.factorial(ell))
            mixed[k, ell] = coefficient
            total += coefficient
        coefficients.append(total)
        z = sum(int(o[i])*r[n][i, i] for i in range(6))*(-1j)**n
        assert z.imag == 0
        continuous.append(Q(int(z.real)*2**n, 768*math.factorial(n)))
    return mixed, coefficients, continuous


def remainder(t):
    x = 36*t
    return x**9/(2*math.factorial(9)*(1-x/10))


class Audit(unittest.TestCase):
    def close(self, a, b, tol=3e-11):
        self.assertLess(float(np.linalg.norm(np.asarray(a)-np.asarray(b))), tol)

    def test_01_storage_handoff_and_complete_graph_reference_channels(self):
        for r in (RX, RY):
            e, k = channel(r, Q(1, 5))
            complete = np.einsum('bdgi,bdgj->ij', k.conj(), k, optimize=True)
            self.close(complete, np.eye(6))
        raw = np.arange(1, 13)+1j*np.arange(12, 0, -1)
        psi = raw/np.linalg.norm(raw)
        rho = np.outer(psi, psi.conj())
        out = apply(channel(RY, Q(1, 5))[0], apply(channel(RX, Q(1, 5))[0], rho, 2), 2)
        ref = lambda x: x.reshape(6, 2, 6, 2).trace(axis1=0, axis2=2)
        self.close(ref(out), ref(rho))
        self.assertGreaterEqual(float(np.linalg.eigvalsh(out).min()), -1e-12)
        # Arbitrary old D-G-R amplitudes move to storage, without deletion.
        old = (np.arange(12)+1j*np.arange(12)[::-1]).reshape(2, 3, 2)
        fresh = np.array([1, 1j])/math.sqrt(2)
        original = np.einsum('dgr,p->dgpr', old, fresh)
        handed = original.transpose(2, 1, 0, 3)
        target = np.einsum('p,sgr->pgsr', fresh, old)
        self.assertTrue(np.array_equal(handed, target))
        OBS['complete_handoff_contract'] = dict(
            natural_channel='E_s(u) E_r(t)', actual_U_Kraus_completeness_checked=True,
            arbitrary_graph_reference_channel_extension_checked=True,
            arbitrary_old_D_G_R_preserved_by_register_swap=True,
            old_data_are_isolated_internal_storage_not_physically_traced_away=True,
            fresh_eta_independent_of_old_data_graph_reference=True,
            same_port_labels_but_new_reference_carriers=True,
            storage_swap_source_preparation_clock_and_isolation_are_extra_inputs=True)

    def test_02_exact_mixed_Taylor_operator_certificate(self):
        h = previous.system()[1]
        self.assertLessEqual(int(np.abs(h).sum(axis=1).max()), 9)
        eta_y = previous.data_basis_raw()[0]+previous.data_basis_raw()[2]
        absolute_gaussian_sum = int(np.abs(eta_y.real).sum()+np.abs(eta_y.imag).sum())
        self.assertEqual(absolute_gaussian_sum, 512)
        self.assertLess(1024*18**8, 2**53)
        mixed, coefficients, continuous = certificate()
        expected = [Q(1, 6), Q(0), Q(0), Q(0), Q(4, 9), Q(43, 192),
                    -Q(25483, 5760), -Q(29327, 17280), Q(9063731, 483840)]
        self.assertEqual(coefficients, expected)
        self.assertEqual(continuous, [Q(1, 6)]+[Q(0)]*8)
        leading = {(k, ell): c for (k, ell), c in mixed.items()
                   if k+ell == 4 and c}
        self.assertEqual(leading, {(1, 3): Q(5, 18), (2, 2): Q(1, 6)})
        OBS['exact_mixed_coefficients'] = dict(
            distance='D_03', centered_observable='O=2*D_03-5*I, norm=1',
            data_state_raw_denominator=64, graph_state_raw_denominator=6,
            mixed_contraction_denominator_before_factorials=49152,
            all_nonzero_mixed=[dict(first_degree=k, second_degree=ell, coefficient=str(c))
                for (k, ell), c in mixed.items() if c],
            equal_wait_coefficients=[str(c) for c in coefficients],
            continuous_2tau_coefficients=[str(c) for c in continuous],
            constants_in_lists_are_D_minus_5_over_2=True,
            rigorous_Gaussian_integer_intermediate_bound=1024*18**8,
            below_exact_double_integer_limit=True,
            final_mixed_products_use_unbounded_Python_integers=True,
            time_even_symmetry_not_assumed=True)

    def test_03_strict_finite_window_and_continuous_comparator(self):
        coefficients = certificate()[1]
        tmin, tmax = Q(1, 256), Q(1, 128)
        correction = sum(abs(coefficients[n])*tmax**(n-4) for n in range(5, 9))
        normalized_tail = remainder(tmax)/tmax**4
        self.assertLess(correction+normalized_tail, Q(2, 45))
        self.assertLess(Q(4, 9)+correction+normalized_tail, Q(1, 2))
        self.assertLess(normalized_tail, Q(1, 40))
        # Repeated D minus initial >2/5 tau^4; continuous deviation <=tail.
        # Thus repeated minus continuous >3/8 tau^4 throughout the window.
        self.assertEqual(Q(2, 5)-Q(1, 40), Q(3, 8))
        OBS['strict_window'] = dict(
            equal_wait_window=[str(tmin), str(tmax)],
            repeated_distance_minus_initial_lower='(2/5)*tau^4',
            repeated_distance_minus_initial_upper='tau^4/2',
            whole_two_stage_Taylor_tail='(36*tau)^9/[2*9!*(1-36*tau/10)]',
            normalized_high_order_polynomial_bound=str(correction),
            normalized_tail_bound=str(normalized_tail),
            continuous_without_handoff_distance_deviation_at_most_same_tail=True,
            repeated_minus_continuous_lower='(3/8)*tau^4',
            uniform_repeated_initial_gap=str(Q(2, 5)*tmin**4),
            uniform_protocol_comparison_gap=str(Q(3, 8)*tmin**4),
            no_all_time_single_wait_theorem_imported=True)

    def test_04_full_U_shape_population_change(self):
        distances = graph_observables()
        rows = []
        for tau in (Q(1, 128), Q(1, 5)):
            before = apply(channel(RX, tau)[0], np.ones((6, 6))/6)
            after = repeated(tau, tau)
            continuous = apply(channel(RX, 2*tau)[0], np.ones((6, 6))/6)
            means = distances@np.diag(after).real
            first = distances@np.diag(before).real
            comparator = distances@np.diag(continuous).real
            difference = float(means[0]-comparator[0])
            self.assertGreater(difference, 0)
            if tau == Q(1, 128):
                self.assertGreater(difference, float(Q(3, 8)*tau**4))
                prediction = float(Q(5, 2)+sum(c*tau**n for n, c in enumerate(certificate()[1])))
                self.assertLess(abs(float(means[0])-prediction), float(remainder(tau))+5e-14)
            matching_probs = [float(np.diag(after).real[distances[k] == 2].sum()) for k in range(3)]
            self.close(sum(matching_probs), 1.)
            rows.append(dict(wait=str(tau), first_leaf_means=[short(v) for v in first],
                repeated_leaf_means=[short(v) for v in means],
                uninterrupted_leaf_means=[short(v) for v in comparator],
                repeated_minus_uninterrupted_D03=short(difference),
                three_split_class_probabilities=[short(v) for v in matching_probs]))
        OBS['actual_population_witnesses'] = dict(
            leaf_pair_order=['03', '04', '05', '34', '35', '45'], rows=rows,
            changed_shape_is_graph_population_not_only_phase=True,
            large_time_example_not_used_for_short_window_proof=True)

    def test_05_each_fresh_batch_can_be_fully_relational(self):
        changes = []
        for r in (RX, RY):
            eta = previous.data_state(r)
            twirled = previous.restricted_twirl(eta)
            transformed, _ = channel_from_eta(twirled, Q(1, 5))
            original = channel(r, Q(1, 5))[0]
            changes.append(float(np.linalg.norm(transformed-original)))
            self.close(transformed, original, 3e-13)
        OBS['no_external_absolute_axis'] = dict(
            full_graph_channel_twirl_difference_norms=[short(v) for v in changes],
            each_batch_collective_SU2_twirl_preserves_E_analytically=True,
            independent_twirl_of_fresh_batches_keeps_composed_channel=True,
            relative_source_reference_preparations_and_label_identification_still_inputs=True,
            twirl_not_claimed_lossless_for_arbitrary_unknown_quantum_source=True)

    def test_06_actual_shape_readout_and_finite_error_budget(self):
        rho = repeated(Q(1, 5), Q(1, 5))
        probe = Q(1, 65536)
        values = []
        # D03=3-P_013-P_023 on this known six-tree department.
        for path in ((0, 1, 3), (0, 2, 3)):
            work = rho.copy()
            for edge in reader.path_edges(path):
                minus, plus, _ = reader.instrument(edge, probe)
                work = apply((plus-minus)/float(probe), work)
            values.append(float(np.trace(work).real))
        actual = float(graph_observables()[0]@np.diag(rho).real)
        estimate = 3-sum(values)
        h = Q(25)
        delta = (math.expm1(float(2*h*probe))-float(2*h*probe))/float(probe)
        self.assertLess(abs(estimate-actual), 2*((1+delta)**2-1))
        gap = Q(3, 8)*Q(1, 256)**4
        error, m, paths = gap/4, 2, 2
        read_error = gap/8
        preparation_process_budget = gap/4
        self.assertEqual(preparation_process_budget/2+read_error, error)
        a = read_error/paths
        u = a/(48*m*h*h)
        gamma = a*u/(8*m)
        x = 2*h*u
        rbound = x*x/(2*(1-x/3))
        bias = (1+(rbound+gamma)/u)**m-1
        self.assertLess(paths*bias, read_error/2)
        self.assertGreater(sum(Q(7)**k/math.factorial(k) for k in range(25)), 800)
        copies = math.ceil(56/(a*a*u**4))
        self.assertEqual(gap-2*error, gap/2)
        OBS['actual_record_and_resource_budget'] = dict(
            actual_D03_at_two_waits_one_fifth=short(actual),
            actual_two_path_instrument_estimate=short(estimate),
            diagnostic_probe_wait=str(probe),
            exact_distance_identity='D03=3-P013-P023',
            uniform_protocol_gap=str(gap), per_protocol_error=str(error),
            source_first_handoff_and_natural_wait_full_process_budget=str(preparation_process_budget),
            source_process_distance_error_at_most=str(preparation_process_budget/2),
            per_protocol_readout_error=str(read_error),
            per_path_total_error=str(a), probe_wait=str(u),
            full_fixed_instrument_diamond_budget=str(gamma),
            control_only_budget=str(gamma/2), total_control_duration_budget=str(gamma/(4*h)),
            copies_per_path_per_protocol=str(copies), total_trials=str(4*copies),
            remaining_protocol_gap=str(gap/2), joint_failure_probability_at_most='1/100',
            max_data_carriers_per_repeat_trial_before_other_ancillas=24,
            handoff_and_readout_errors_are_full_reference_contracts=True,
            no_feedback_from_past_random_records_or_residual_controller_memory=True,
            old_carriers_records_and_controllers_remain_inside_total_system=True,
            natural_preparation_first_handoff_and_wait_errors_require_complete_process_contract=True,
            no_specific_autonomous_handoff_hardware_simulated=True,
            no_huge_repetition_experiment_or_autonomous_readout_implementation_claim=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=476, baseline_round=474, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(specific_repeated_fresh_reference_protocol_changes_leaf_shape=True,
            original_fixed_H_during_each_natural_wait=True,
            ideal_storage_handoff_and_independent_fresh_batch_are_added_contracts=True,
            complete_old_quantum_information_retained_in_isolated_internal_memory=True,
            graph_partial_trace_is_only_reduced_description=True,
            not_a_counterexample_to_all_possible_reference_protocols=True,
            not_an_axiomatic_refutation_of_cognition_or_three_dimensional_space=True,
            position_not_required_to_predict_all_future_internal_information=True,
            physical_three_space_or_GR_derived=False, full_GR_goal_completed=False,
            phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert result == json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
