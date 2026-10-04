"""Round 481: positive-time reaccess to the same current distance-readout chart.

Scientific baseline 480. Repetition, calibrated controls, the correlated source,
and the actual data-readout contract remain explicit model inputs. The theorem
uses recurrence of a known finite unitary word and a quantitative covering map;
it neither implements exact inverse evolution nor identifies physical space.
"""
import argparse
from fractions import Fraction as F
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import correlated_reference_local_access as access
import rigid_leaf_reference_audit as frame
import symmetric_control_reachability_audit as words
import sequential_tree_distance_readout_audit as reader

TARGET = Path(__file__).with_name('positive_time_coordinate_reaccess_results.json')
OBS = {}
DIM = 384
B = 2*access.BINV
L = 2052


def ceil_fraction(x):
    return (x.numerator+x.denominator-1)//x.denominator


def short(x):
    return float(f'{float(x):.11g}')


def budgets():
    _, old_radius, _, _, _ = access.exact_budgets()
    radius = old_radius/2
    epsilon = radius/(16*B)
    target_radius = radius/(4*B)
    m = ceil_fraction(32*DIM/epsilon)
    return old_radius, radius, epsilon, target_radius, m


def rational_rotation_return(tolerance=F(1, 100), limit=1000):
    # P=[[3,-4],[4,3]]/5. Every numerator here is a Python integer.
    real, imag, denominator = 1, 0, 1
    for n in range(1, limit+1):
        real, imag = 3*real-4*imag, 4*real+3*imag
        denominator *= 5
        assert real*real+imag*imag == denominator*denominator
        frobenius_squared = F(4*(denominator-real), denominator)
        if frobenius_squared < tolerance*tolerance/4:
            return n, real, imag, denominator, frobenius_squared
    raise RuntimeError('The diagnostic search limit was insufficient.')


def graph_from_columns(columns, denominator=384):
    state = columns@columns.conj().T/denominator
    return state.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)


def movement_data():
    # Actual finite H-on word from 480; no ideal negative-time factor is used.
    if not hasattr(movement_data, 'cached'):
        p = access.numerical_control()[0]
        c = access.pure_columns()
        p2 = p@p
        p3 = p@p2
        movement_data.cached = (p, p2, p3, c)
    return movement_data.cached


class Audit(unittest.TestCase):
    def test_01_frozen_chart_and_current_center_cover(self):
        frozen = json.loads(access.TARGET.read_text(encoding='utf-8'))
        self.assertEqual((frozen['round'], frozen['tests_run'], frozen['failures'], frozen['errors']),
                         (480, 6, 0, 0))
        old_radius, radius, eps, target, m = budgets()
        self.assertEqual(old_radius, F(1, 246240000))
        self.assertEqual(B*L*old_radius, F(1, 4))
        self.assertLess(3*old_radius/4, F(1, 5))
        self.assertEqual(B*target+radius/4+2*B*eps, 5*radius/8)
        self.assertEqual(old_radius/4+radius, 3*old_radius/4)
        self.assertEqual(target, F(1, 59097600000000))
        self.assertEqual(eps, F(1, 236390400000000))
        self.assertLessEqual(F(8*DIM, m), eps/4)
        bound = m**DIM
        OBS['current_center_cover'] = dict(
            frozen_baseline_round=480, original_parameter_cube_radius=str(old_radius),
            inverse_Jacobian_sup_norm_upper=B, Jacobian_Lipschitz_sup_norm_upper=L,
            admitted_basepoint_distance_from_old_center=str(old_radius/4),
            local_parameter_search_radius=str(radius),
            remaining_old_domain_margin=str(old_radius/4),
            actual_current_Q_target_ball_sup_radius=str(target),
            recurrence_operator_error_budget=str(eps),
            Brouwer_self_map_radius_ratio='5/8',
            basepoint_may_be_any_point_in_inner_quarter_cube=True,
            current_state_may_be_any_known_positive_word_P_of_same_seed=True,
            current_Q_must_equal_f_of_admitted_basepoint=True,
            no_assumption_that_three_Q_predict_future_state=True)
        OBS['finite_recurrence_bound'] = dict(
            Hilbert_dimension=DIM, Dirichlet_integer_mesh=m,
            maximum_positive_repetition_count=f'{m}**{DIM}',
            maximum_counter_bit_length=bound.bit_length(),
            exact_integer_bound_computed_not_executed=True,
            recurrence_search_can_depend_on_entire_current_word_P=True,
            no_negative_H_or_compact_generated_subgroup_assumed=True,
            each_atomic_rational_exponential_error='eps/(64*d*n*m_P)',
            rational_products_may_be_evaluated_exactly_after_atomic_approximation=True,
            full_Frobenius_matrix_error_upper='eps/32 < eps/8')

    def test_02_exact_finite_recurrence_certificate(self):
        n, a, b, denominator, norm2 = rational_rotation_return()
        self.assertEqual(n, 393)
        self.assertNotEqual((a, b), (denominator, 0))
        self.assertGreater(norm2, 0)
        eps = F(1, 100)
        m = ceil_fraction(32*2/eps)
        self.assertLess(n, m*m)
        self.assertLess(norm2, eps*eps/4)
        self.assertLess(eps/2+eps/8, eps)
        self.assertLess(eps/4+eps/8, eps/2)
        # ||P^(n-1)-P^*|| = ||P^n-I||, by unitary invariance.
        p = np.array([[.6, -.8], [.8, .6]])
        measured = np.linalg.norm(np.linalg.matrix_power(p, n)-np.eye(2), 2)
        self.assertLess(abs(measured*measured-float(norm2/2)), 1e-13)
        OBS['exact_small_return_certificate'] = dict(
            word='[[3,-4],[4,3]]/5', dimension=2, positive_repetitions=n,
            Frobenius_residual_squared=short(norm2),
            residual_operator_norm=short(math.sqrt(float(norm2/2))),
            exact_nonzero_residual_proves_no_exact_rollback_in_this_example=True,
            real_numerator_bits=a.bit_length(), imaginary_numerator_bits=abs(b).bit_length(),
            denominator_bits=denominator.bit_length(), all_powers_exact_gaussian_integers=True,
            universal_Dirichlet_bound=m*m,
            certificate_search='rational matrix error <= eps/8; accept squared Frobenius residual < (eps/2)^2',
            existence_margin_for_some_candidate_within_Dirichlet_bound='true Frobenius residual < eps/4',
            diagnostic_393rd_candidate_only_claims_the_stated_acceptance_test=True,
            certified_operator_error_less_than='5*eps/8',
            example_not_a_model_of_spatial_dimension=True)

    def test_03_center_shift_requires_parameter_adjustment(self):
        # A finite example of the covering lemma; not a substitute for the 480 chart.
        offset = np.array([1/8, -1/8, 1/16])
        def nonlinear(v):
            return np.roll(v, -1)**2/16
        def g(v):
            return v+nonlinear(v)+offset
        targets = [np.zeros(3)]+[np.array(s)/4 for s in itertools.product((-1, 1), repeat=3)]
        errors, adjusted = [], []
        for target in targets:
            v = np.zeros(3)
            for _ in range(16):
                v = target-nonlinear(v)-offset
            errors.append(np.max(abs(g(v)-target)))
            adjusted.append(v)
            self.assertLessEqual(np.max(abs(v)), 7/16)
        self.assertGreater(np.linalg.norm(adjusted[0]), .1)
        self.assertLess(max(errors), 1e-13)
        self.assertEqual(F(1, 4)+F(1, 16)+F(1, 8), F(7, 16))
        OBS['covering_lemma_finite_crosscheck'] = dict(
            map='g_j(v)=v_j+v_(j+1)^2/16+offset_j',
            offset=['1/8', '-1/8', '1/16'], cube_radius=1,
            derivative_remainder_sup_norm_bound='1/8',
            tested_targets=9, fixed_point_iterations=16,
            self_map_radius_bound='7/16', maximum_endpoint_residual=short(max(errors)),
            parameter_adjustment_for_original_current_center=[short(x) for x in adjusted[0]],
            original_current_center_not_assumed_equal_shifted_map_center=True,
            generic_quantitative_lemma_crosscheck_not_a_second_spatial_model=True)

    def test_04_actual_positive_words_keep_same_system_and_flag(self):
        p, p2, p3, c = movement_data()
        readouts, leaf_errors, flag_weights = [], [], []
        for w in (p, p2, p3):
            y = w@c
            graph = graph_from_columns(y)
            readouts.append([short(x) for x in words.qgraph()@np.diag(graph).real])
            weights = []
            for subset in (slice(0, 4), slice(4, 28)):
                branch = graph_from_columns(y[:, subset], 192)
                leaf_errors.append(max(abs(np.diag(branch).real@frame.distance(a, b)-8/3)
                                       for a, b in itertools.combinations(frame.LEAVES, 2)))
                weights.append(short(np.linalg.norm(y[:, subset])**2/384))
            flag_weights.append(weights)
        self.assertLess(max(leaf_errors), 5e-12)
        self.assertGreater(np.max(abs(np.array(readouts[0])-np.array(readouts[1]))), .001)
        rng = np.random.default_rng(481)
        unknown = rng.normal(size=(384, 3))+1j*rng.normal(size=(384, 3))
        unknown /= np.linalg.norm(unknown)
        inverse_residual = np.linalg.norm(p3.conj().T@(p3@unknown)-unknown)
        self.assertLess(inverse_residual, 5e-12)
        OBS['actual_six_tree_word_check'] = dict(
            finite_H_on_pulse_duration=str(access.TAU),
            repeated_word_parameters=['1/5', '2/5', 'pi/2'],
            positive_word_counts=[1, 2, 3], Q_readouts=readouts,
            same_initial_columns_evolved_without_any_repreparation=True,
            branch_record_probabilities=flag_weights,
            maximum_conditional_branch_leaf_mean_error=short(max(leaf_errors)),
            unknown_data_graph_reference_inverse_identity_residual=short(inverse_residual),
            inverse_used_only_to_check_unitarity_not_executed_as_control=True,
            three_short_words_not_claimed_to_meet_the_tiny_recurrence_budget=True,
            full_finite_Dirichlet_repetition_protocol_not_numerically_executed=True)

    def test_05_finite_control_and_error_costs(self):
        _, radius, _, target, m = budgets()
        tolerance = target/16
        n = 393  # A diagnostic count; the general formulas retain arbitrary n.
        old_error = tolerance/12
        each_repeat_error = tolerance/(12*max(1, n-1))
        last_word_error = tolerance/12
        total = old_error+(n-1)*each_repeat_error+last_word_error
        self.assertEqual(total, tolerance/4)
        grid_spacing = tolerance/(16*38)
        function_calibration = tolerance/16
        self.assertEqual(38*grid_spacing+2*function_calibration, 3*tolerance/16)
        self.assertLess(3*tolerance/16, tolerance/4)
        steps = ceil_fraction(2*radius/grid_spacing)
        self.assertGreater(steps, 0)
        # Rectangular pulse and waits: full channel error <= 2*(9*dt+9*du+9*dtau+dtheta).
        parameter_error = last_word_error/56
        self.assertEqual(2*(9+9+9+1)*parameter_error, last_word_error)
        OBS['finite_control_cost_ledger'] = dict(
            requested_final_Q_sup_error=str(tolerance),
            whole_old_process_error_budget='zeta/12',
            per_repeated_P_word_diamond_budget='zeta/(12*max(1,n-1))',
            last_W_word_diamond_budget='zeta/12',
            accumulated_control_process_diamond_budget='zeta/4',
            three_parameter_Q_Lipschitz_sup_norm_upper=38,
            finite_parameter_grid_mesh=str(grid_spacing),
            rigorous_Q_minus_target_evaluation_error=str(function_calibration),
            finite_grid_points_upper=str((steps+1)**3),
            endpoint_inclusive_grid='z_i-r+2*r*j/k, j=0..k, k=ceil(2*r/h)',
            certified_grid_acceptance_target='zeta/4',
            elapsed_time='(n-1)*T_P + t + u + tau',
            complete_new_history_time='n*T_P + t + u + tau',
            primitive_segments='(n-1)*m_P + 3',
            peak_exchange_strength='maximum of past |theta_j|/tau_j and |theta|/tau',
            rectangular_pulse_unitary_error_upper='9*|delta_tau|+|delta_theta|',
            drift_model_error_adds_full_channel_bound='2*T_total*||delta_H||',
            counter_program_storage_time_calibration_and_isolation_are_resources=True,
            complete_history_word_P_must_remain_known_and_repeatable=True,
            known_computable_controls_required_for_terminating_numerical_search=True,
            no_claim_of_autonomous_finite_controller_for_all_real_requests=True,
            previous_unknown_noise_cannot_be_reduced_by_recurrence=True,
            finite_precision_guarantee_is_for_each_preassigned_finite_protocol=True,
            arbitrary_R_error_quantifier_is_complete_process_diamond=True,
            isolated_gate_calibration_not_substituted_for_memoryful_step_contract=True,
            old_control_ancillas_must_be_isolated_or_included_in_complete_step_map=True)

    def test_06_actual_readout_and_finite_statistics(self):
        _, p2, _, c = movement_data()
        graph = graph_from_columns(p2@c)
        true = words.qgraph()@np.diag(graph).real
        probe = F(1, 65536)
        means, records = [], []
        for leaf in (0, 3, 4, 5):
            minus, plus, _ = reader.instrument((1, leaf), probe)
            probabilities = [float(np.trace(np.einsum('ghij,ij->gh', op, graph)).real)
                             for op in (minus, plus)]
            self.assertLess(abs(sum(probabilities)-1), 4e-12)
            self.assertGreaterEqual(min(probabilities), -4e-12)
            means.append((probabilities[1]-probabilities[0])/float(probe))
            records.append([short(x) for x in probabilities])
        measured = means[0]-np.array(means[1:])
        x = 18*probe
        bias = x*x/(2*(1-x/3)*probe)
        self.assertLess(np.max(abs(measured-true)), float(2*bias)+2e-9)
        _, _, _, target, _ = budgets()
        tolerance = target/16
        a = tolerance/8
        h = a/(16*81)
        gamma = a*h/4
        x = 18*h
        self.assertLessEqual((x*x/(2*(1-x/3))+gamma)/h, a/2)
        copies = ceil_fraction(80/(h*h*a*a))
        self.assertGreaterEqual(copies*h*h*a*a/8, 10)
        self.assertGreater(sum(F(10**j, math.factorial(j)) for j in range(10)), 800)
        self.assertEqual(tolerance/4+tolerance/4+2*a, 3*tolerance/4)
        OBS['actual_readout_and_statistics'] = dict(
            actual_CP_instrument_round=472, ordered_edges=[[1,j] for j in (0,3,4,5)],
            diagnostic_probe_wait=str(probe), two_outcome_probabilities=records,
            true_Q=[short(x) for x in true], measured_diagnostic_Q=[short(x) for x in measured],
            diagnostic_Q_bias_bound=str(2*bias),
            certified_Q_error_tolerance=str(tolerance),
            each_adjacency_mean_error=str(a), certified_probe_wait=str(h),
            whole_storage_handoff_plus_instrument_diamond_budget=str(gamma),
            independent_copies_per_edge=str(copies), total_repeated_histories=str(4*copies),
            simultaneous_failure_probability_less_than='1/100',
            final_Q_error_upper='3*zeta/4 < zeta',
            old_data_graph_reference_correlations_saved_not_physically_discarded=True,
            final_readout_requires_fresh_independent_probe_and_storage_contract=True,
            tomography_or_mean_estimation_not_free_nondisturbing_readout_of_one_system=True,
            copies_replay_known_preparation_and_whole_history_not_clone_unknown_R=True,
            declared_tiny_error_resource_protocol_not_actually_sampled=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=481, baseline_round=480, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(positive_time_local_current_Q_reaccess_in_fixed_model=True,
            exact_covering_is_a_mathematical_continuous_control_statement=True,
            finite_precision_full_process_and_measurement_costs_explicit=True,
            current_Q_does_not_have_to_predict_all_future_internal_states=True,
            no_exact_full_state_rollback_or_free_negative_H_claim=True,
            no_generated_control_group_closedness_assumption=True,
            no_source_or_controller_generation_from_cognitive_principles=True,
            no_spatial_displacement_composition_or_dimension_selection=True,
            endpoint_cover_does_not_keep_intermediate_trajectory_in_local_chart=True,
            no_small_time_local_controllability_or_vanishing_cost_claim=True,
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
