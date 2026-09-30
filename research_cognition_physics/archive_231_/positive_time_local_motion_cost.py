"""Round 485: finite endpoint reaccess does not supply vanishing two-sided action.

Scientific baselines 480 and 481. The actual finite-pulse center is retained.
Allowed continuation is H+u(t)K for positive time and finite integral |u|.
The explicit action budget 9*T+integral|u| is an input, not a unique cognitive
cost or physical dissipation. No conclusion against space or other controls.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import correlated_reference_local_access as access
import positive_time_coordinate_reaccess as replay
import rigid_leaf_reference_audit as frame
import symmetric_control_reachability_audit as words
import sequential_tree_distance_readout_audit as reader

TARGET = Path(__file__).with_name('positive_time_local_motion_cost_results.json')
OBS = {}
HNORM = 9
QNORM = 1
AXNORM = 2*HNORM*QNORM
VLOWER = F(1, 8)
ACTION_CUTOFF = F(1, 576)
SPEED_LOWER = F(1, 16)
# Durations and signed pulse areas; all generators are dt*H+area*K.
SCHEDULE = ((F(1, 200000), F(1, 40000)),
            (F(1, 400000), -F(1, 80000)),
            (F(1, 400000), F(1, 80000)))


def short(value):
    return float(f'{float(value):.12g}')


def action_budget(schedule):
    return sum((HNORM*duration+abs(area) for duration, area in schedule), F(0))


def graph_from_columns(columns, denominator=384):
    state = columns@columns.conj().T/denominator
    return state.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)


def q_from_columns(columns):
    weights = np.sum(abs(columns)**2, axis=1)/384
    return np.tile(words.qgraph(), (1, 64))@weights


@lru_cache(None)
def actual_continuation():
    h = frame.system()[1]
    k = words.s12()
    p = access.numerical_control()[0]
    center = p@access.pure_columns()
    output = center.copy()
    continuation = np.eye(384, dtype=complex)
    points = [q_from_columns(output)]
    for duration, area in SCHEDULE:
        step = words.unitary(float(duration)*h+float(area)*k)
        output = step@output
        continuation = step@continuation
        points.append(q_from_columns(output))
    return p, center, output, continuation, points


class Audit(unittest.TestCase):
    def test_01_actual_observable_and_control_identities(self):
        h = frame.system()[1]
        k0 = words.s12()
        self.assertTrue(np.array_equal(k0, k0.astype(np.int64)))
        k = k0.astype(np.int64)
        self.assertTrue(np.issubdtype(h.dtype, np.integer))
        q = np.tile(words.qgraph(), (1, 64)).astype(np.int64)
        self.assertTrue(np.array_equal(h, h.T))
        self.assertEqual(int(np.max(np.sum(abs(h), axis=1))), HNORM)
        self.assertTrue(np.array_equal(k@k, np.eye(384, dtype=np.int64)))
        commutators = []
        for diagonal in q:
            self.assertLessEqual(int(np.max(abs(diagonal))), QNORM)
            k_comm = k*diagonal[None, :]-diagonal[:, None]*k
            self.assertEqual(int(np.max(abs(k_comm))), 0)
            # [H,Q] is real antisymmetric; i[H,Q] is Hermitian.
            h_comm = h*diagonal[None, :]-diagonal[:, None]*h
            self.assertTrue(np.array_equal(h_comm, -h_comm.T))
            commutators.append(int(np.max(np.sum(abs(h_comm), axis=1))))
        OBS['exact_control_and_observables'] = dict(
            active_dimension=384, H_norm_bound=HNORM, K='SWAP_12 tensor I_G',
            K_norm=1, Q_order=['D13-D10', 'D14-D10', 'D15-D10'],
            Q_norm_bound=QNORM, K_commutes_with_each_Q_exactly=True,
            commutator_row_norms=commutators,
            general_velocity_operator_bound=AXNORM,
            original_H_remains_on_and_negative_u_is_allowed=True,
            no_negative_H_or_reset_or_instantaneous_pulse_permission=True,
            integer_matrix_equalities_not_floating_tolerance_tests=True)

    def test_02_strict_velocity_at_actual_finite_pulse_center(self):
        certificate = access.certificate()
        rational_velocity = certificate['jacobian'][0][1]
        self.assertLess(certificate['entry_error'], access.ENTRY_ERROR)
        error = access.ENTRY_ERROR+324*access.TAU
        lower, upper = rational_velocity-error, rational_velocity+error
        self.assertGreater(lower, VLOWER)
        self.assertLess(upper, F(3, 20))
        _, center, _, _, _ = actual_continuation()
        h = frame.system()[1]
        q = np.tile(words.qgraph()[0], 64)
        hcomm = h*q[None, :]-q[:, None]*h
        velocity = float(np.sum(center.conj()*(1j*hcomm@center)).real/384)
        self.assertLess(abs(velocity-float(rational_velocity)), float(error)+3e-12)
        OBS['actual_center_velocity_certificate'] = dict(
            preparation_and_correlated_branch_record_round=480,
            W_center='exp(-i*(2/5)*H) exp(-i*(tau*H+(pi/2)*K)) exp(-i*(1/5)*H)',
            tau=str(access.TAU), center_is_actual_finite_H_on_output=True,
            rational_ideal_velocity=str(rational_velocity),
            finite_pulse_and_integer_entry_error_upper=str(error),
            rational_certified_lower=str(lower), rational_certified_upper=str(upper),
            strict_simple_velocity_interval=['1/8', '3/20'],
            actual_center_numerical_velocity=short(velocity),
            exact_frozen_480_Python_integer_certificate_recomputed=True,
            numerical_velocity_not_used_as_strict_proof=True)

    def test_03_all_finite_controls_action_cone_and_diagnostic(self):
        self.assertEqual(VLOWER-2*AXNORM*ACTION_CUTOFF, SPEED_LOWER)
        self.assertEqual(VLOWER-AXNORM*F(1, 576)-2*AXNORM*F(1, 1152), SPEED_LOWER)
        duration = sum((t for t, _ in SCHEDULE), F(0))
        action = action_budget(SCHEDULE)
        self.assertEqual(duration, F(1, 100000))
        self.assertEqual(action, F(7, 50000))
        self.assertLess(action, ACTION_CUTOFF)
        _, _, _, _, points = actual_continuation()
        delta = points[-1][0]-points[0][0]
        strong_speed = VLOWER-2*AXNORM*action
        self.assertGreaterEqual(delta+3e-12, float(strong_speed*duration))
        for index, (dt, _) in enumerate(SCHEDULE):
            self.assertGreaterEqual(points[index+1][0]-points[index][0]+3e-12,
                                    float(SPEED_LOWER*dt))
        OBS['uniform_action_cone'] = dict(
            admitted_controls='real L1 control u on [0,T], T>0; in particular every finite piecewise-constant word',
            declared_budget='A=9*T+integral_0^T |u(t)| dt',
            all_prefix_trace_norm_distance_upper='2*A',
            all_prefix_velocity_change_upper='36*A',
            cutoff=str(ACTION_CUTOFF), uniform_positive_velocity_lower=str(SPEED_LOWER),
            endpoint_inequality='Qx(T)-Qx(0) >= T/16 > 0 whenever A <= 1/576',
            derivation_uses_unitary_integral_and_commutator_norm_not_time_grid=True,
            arbitrary_large_finite_strength_is_covered_by_integrated_budget=True,
            initial_trace_error_eta_velocity_lower='1/8-18*eta-36*A',
            robust_example='eta<=1/576 and A<=1/1152 still give velocity>=1/16',
            diagnostic_schedule=[dict(duration=str(t), signed_area=str(a),
                                      strength=str(a/t)) for t, a in SCHEDULE],
            diagnostic_elapsed_time=str(duration), diagnostic_action=str(action),
            diagnostic_stronger_lower_gap=str(strong_speed*duration),
            diagnostic_simple_lower_gap=str(SPEED_LOWER*duration),
            diagnostic_Q_at_segment_endpoints=[[short(v) for v in p] for p in points],
            diagnostic_actual_delta_Qx=short(delta))

    def test_04_finite_reaccess_and_nonvanishing_reverse_action(self):
        frozen = json.loads(replay.TARGET.read_text(encoding='utf-8'))
        self.assertEqual((frozen['round'], frozen['tests_run'], frozen['failures'], frozen['errors']),
                         (481, 6, 0, 0))
        radius, _, _, target, _ = replay.budgets()
        self.assertEqual(radius, F(1, 246240000))
        self.assertEqual(target, F(1, 59097600000000))
        examples = [target/F(2**n) for n in (1, 2, 4, 8)]
        self.assertTrue(all(0 < epsilon < target for epsilon in examples))
        self.assertEqual(2*(target/4), target/2)
        OBS['endpoint_reaccess_versus_cost'] = dict(
            frozen_positive_time_reaccess_round=481,
            current_word_P='the finite H-on center W_tau(1/5,2/5,pi/2)',
            actual_current_Q_target_sup_radius=str(target),
            backward_target_family='y_n=Q_current-(delta/2**n)*e_x, n>=1',
            example_negative_Qx_offsets=[str(v) for v in examples],
            all_these_exact_endpoints_have_finite_positive_time_implementations_by_481=True,
            every_exact_implementation_in_current_H_plus_uK_class_has_action_greater_than=str(ACTION_CUTOFF),
            infimum_reverse_action_is_at_least=str(ACTION_CUTOFF),
            approximate_reverse_endpoint_with_Qx_error_less_than_epsilon_also_obstructed=True,
            no_reverse_action_sequence_tends_to_zero_as_target_tends_to_current_Q=True,
            no_uniform_positive_lower_bound_on_elapsed_time_claimed=True,
            not_a_repeat_of_415_torus_exact_cost_or_416_charged_turn_counterexample=True,
            no_condition_that_Q_determine_all_future_states=True)

    def test_05_retained_reference_record_and_unknown_information(self):
        p, center, output, continuation, _ = actual_continuation()
        errors, weights = [], []
        for columns in (center, output):
            for subset in (slice(0, 4), slice(4, 28)):
                branch = graph_from_columns(columns[:, subset], 192)
                errors.append(max(abs(np.diag(branch).real@frame.distance(a, b)-8/3)
                                  for a, b in itertools.combinations(frame.LEAVES, 2)))
                weights.append(np.linalg.norm(columns[:, subset])**2/384)
        self.assertLess(max(errors), 5e-12)
        self.assertLess(max(abs(v-.5) for v in weights), 4e-12)
        rng = np.random.default_rng(485)
        arbitrary = rng.normal(size=(384, 3))+1j*rng.normal(size=(384, 3))
        arbitrary /= np.linalg.norm(arbitrary)
        reverse_error = np.linalg.norm(continuation.conj().T@(continuation@arbitrary)-arbitrary)
        self.assertLess(reverse_error, 5e-12)
        OBS['reference_and_information'] = dict(
            all_admitted_H_plus_uK_controls_share_475_leaf_S4_and_collective_SU2_symmetries=True,
            correlated_seed_each_branch_has_zero_22_sector=True,
            full_time_leaf_matching_probabilities_one_third_preserved=True,
            conditional_on_retained_classical_branch_record_also_preserved=True,
            diagnostic_maximum_leaf_mean_error=short(max(errors)),
            diagnostic_branch_weights=[short(x) for x in weights],
            arbitrary_unknown_reference_information_kept_by_full_control_unitarity=True,
            diagnostic_inverse_identity_residual=short(reverse_error),
            mathematical_inverse_not_admitted_as_physical_negative_time=True,
            no_arbitrary_R_conditional_leaf_moment_invariance_from_unitarity_alone=True,
            original_data_graph_correlations_never_reset_during_continuation=True,
            classical_control_history_may_be_retained_in_internal_storage=True,
            autonomous_generation_of_controllers_not_proved=True)

    def test_06_actual_CP_readout_and_full_error_budget(self):
        _, center, output, _, _ = actual_continuation()
        probe = F(1, 65536)
        records, measured, truth = [], [], []
        for columns in (center, output):
            graph = graph_from_columns(columns)
            true = float(words.qgraph()[0]@np.diag(graph).real)
            means, local = [], []
            for leaf in (0, 3):
                minus, plus, _ = reader.instrument((1, leaf), probe)
                probabilities = [float(np.trace(np.einsum('ghij,ij->gh', op, graph)).real)
                                 for op in (minus, plus)]
                self.assertLess(abs(sum(probabilities)-1), 5e-12)
                self.assertGreaterEqual(min(probabilities), -5e-12)
                means.append((probabilities[1]-probabilities[0])/float(probe))
                local.append(dict(ordered_edge=[1, leaf], measured_port=1,
                                  probabilities=[short(x) for x in probabilities]))
            truth.append(true)
            measured.append(means[0]-means[1])
            records.append(local)
        x = 18*probe
        bias = x*x/(2*(1-x/3)*probe)
        self.assertLess(np.max(abs(np.array(truth)-measured)), float(2*bias)+3e-9)
        duration = sum((t for t, _ in SCHEDULE), F(0))
        gap = SPEED_LOWER*duration
        source_error = gap/8
        a = gap/16
        h = a/(16*81)
        gamma = a*h/4
        x = 18*h
        self.assertLessEqual((x*x/(2*(1-x/3))+gamma)/h, a/2)
        copies = replay.ceil_fraction(80/(h*h*a*a))
        self.assertGreaterEqual(copies*h*h*a*a/8, 10)
        self.assertGreater(sum(F(10**j, math.factorial(j)) for j in range(10)), 800)
        remaining = gap-2*source_error-4*a
        self.assertEqual(remaining, gap/2)
        OBS['actual_readout_and_cost_ledger'] = dict(
            actual_CP_instrument_round=472, diagnostic_probe_wait=str(probe),
            diagnostic_records=records,
            true_Qx_before_after=[short(v) for v in truth],
            instrument_Qx_before_after=[short(v) for v in measured],
            diagnostic_each_Qx_bias_upper=str(2*bias),
            diagnostic_probe_not_claimed_to_certify_the_tiny_gap=True,
            rigorous_ideal_forward_gap=str(gap),
            each_setting_complete_preparation_and_control_diamond_budget=str(source_error),
            same_actual_prepared_state_for_both_edge_readouts_in_each_setting=True,
            edge_choice_only_after_common_source_preparation_no_prior_record_feedback=True,
            independent_edge_specific_source_errors_would_need_a_different_budget=True,
            each_adjacency_mean_error=str(a), certified_probe_wait=str(h),
            complete_storage_handoff_plus_instrument_diamond_budget=str(gamma),
            separate_settings=2, independently_estimated_means=4,
            independent_repeated_known_histories_per_mean=str(copies),
            total_history_preparations=str(4*copies),
            simultaneous_failure_probability_less_than='1/100',
            retained_positive_measured_difference_lower=str(remaining),
            all_finite_preparation_clock_control_probe_and_storage_resources_counted=True,
            no_enormous_sampling_run_or_cloning_of_unknown_input=True,
            old_D_moved_to_internal_isolated_S_and_G_remains_active=True,
            old_S_G_C_R_correlations_retained_at_handoff=True,
            final_CP_readout_can_disturb_G_and_is_not_claimed_reversible=True,
            theorem_concerns_pre_readout_endpoint_and_probe_wait_not_added_to_motion=True,
            action_budget_not_identified_with_dissipation_or_unique_cognitive_resource=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=485, baseline_round=481, scientific_baselines=[480,481],
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(actual_finite_pulse_center_has_strict_small_action_one_sided_motion=True,
            exact_finite_backward_Q_reaccess_still_holds=True,
            chosen_action_cannot_supply_continuous_two_sided_displacement_cost=True,
            result_limited_to_declared_H_plus_uK_class_and_budget=True,
            no_physical_space_or_other_control_implementation_no_go=True,
            no_qubit_or_three_control_parameter_dimension_selection=True,
            no_complete_cognitive_axiom_countermodel_claim=True,
            full_GR_goal_completed=False, phase_closure_triggered=False))


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
