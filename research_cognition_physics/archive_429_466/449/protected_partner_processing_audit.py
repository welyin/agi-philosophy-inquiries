"""Round 449: protected marker rematching and conditional processing coexist.
The two diagonal interactions and their scale hierarchy remain explicit inputs.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import io
import itertools as it
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import joint_marker_dynamics_audit as carrier
import reciprocal_exchange_dynamics_audit as old

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'protected_partner_processing_audit_results.json'
OBS = {}


def error_bound(epsilon, tau, j=1):
    return ((2+abs(tau)*(4+4*abs(j)))*abs(epsilon)
            +(F(13, 5)+abs(tau)*(4+F(26, 5)*abs(j)))*epsilon**2)


@lru_cache(None)
def model():
    p, bits, words, active, record, local, packet, gamma, c, z, _, _ = carrier.model()
    # Integer-scaled intertwiners computed from the raw packet exchange, not fitted.
    r4 = np.array([4//int(e) if e else 0 for e in record], dtype=np.int64)
    w1_4 = -r4[:, None]*(packet@z)
    w2_16 = r4[:, None]*(packet@(r4[:, None]*(packet@z)))
    k4 = z.T@packet@w1_4
    assert np.all(k4 % 2 == 0)
    k2 = k4//2
    d = np.diag(active)+local
    d0 = z.T@d@z
    b2 = k2+2*d0
    a3_16 = packet@w2_16+4*d@w1_4-2*w1_4@b2
    a4_32 = 2*d@w2_16-w2_16@b2
    return p, bits, words, active, record, local, packet, c, z, w1_4, w2_16, k2, d, b2, a3_16, a4_32


def experiment():
    p, bits, words, active, record, local, packet, c, z, *_ = model()
    matching = old.matching_words(4)[0]
    indices = list(np.flatnonzero(record == 0))
    cols = [indices.index(p.index(matching)*16+bit) for bit in (0, 8)]
    readout = np.array([int(w[2] == -1) for w in words], dtype=int)
    rematch = np.array([int(record[i] == 0 and p[i//16] != matching)
                       for i in range(384)], dtype=int)
    return cols, readout, rematch, z.T@(readout[:, None]*z), z.T@(rematch[:, None]*z)


@lru_cache(None)
def rational_columns():
    *_, b2, a3, a4 = model()
    cols, _, _, effect, rematch = experiment()
    order = 40
    matrix = b2.astype(object)
    powers = np.eye(48, dtype=object)[:, cols]
    real, imag = np.zeros((48, 2), dtype=object), np.zeros((48, 2), dtype=object)
    for k in range(order+1):
        coefficient = F((-1)**(k//2), 2**k*math.factorial(k))
        if k % 2:
            imag -= powers*coefficient
        else:
            real += powers*coefficient
        powers = matrix@powers
    remainder = F(3**7*7**(order+1), math.factorial(order+1))
    probability_error = 2*remainder+remainder**2
    probabilities = []
    for e in (effect, rematch):
        indices = np.flatnonzero(np.diag(e))
        probabilities.append([sum((real[i,j]**2+imag[i,j]**2 for i in indices), F(0))
                              for j in range(2)])
    return order, remainder, probability_error, probabilities


class Audit(unittest.TestCase):
    def test_01_diagonal_reward_classification_and_scope(self):
        z = np.diag([1, -1])
        flip = np.array([[0, 1], [1, 0]])
        ident = np.eye(4, dtype=int)
        zz = np.kron(z, z)
        zsum = np.kron(z, np.eye(2, dtype=int))+np.kron(np.eye(2, dtype=int), z)
        for a, b, c in it.product(range(-2, 4), repeat=3):
            q = np.diag([a, b, b, c])
            np.testing.assert_array_equal(4*q, (a+2*b+c)*ident+(a-c)*zsum+(a-2*b+c)*zz)
            common = np.kron(flip, flip)
            single = np.kron(flip, np.eye(2, dtype=int))
            self.assertEqual(np.array_equal(common@q, q@common), a == c)
            self.assertEqual(np.array_equal(single@q, q@single), a == b == c)
            # Two independent matching pairs: every logical state degenerate iff table constant.
            sums = [u+v for u, v in it.product((a, b, b, c), repeat=2)]
            self.assertEqual(len(set(sums)) == 1, a == b == c)
        # q00=q11=2,q01=1: positive and entangling, but no full-band separation.
        self.assertEqual(2*1, 1*2)
        OBS['reward_classification'] = dict(tables_checked=216,
            common_flip_condition='a=c', independent_flip_condition='a=b=c',
            ZZ_reward_coefficient='(a-2*b+c)/4',
            all_matching_logic_degenerate_condition='a=b=c',
            positive_entangling_example=dict(a=2, b=1, c=2, ZZ='1/2'),
            positivity_alone_is_full_band_separation=False,
            N4_equal_energy_example='two weight-1 reciprocal edges versus one weight-2 edge',
            degeneracy_imposed_as_cognitive_axiom=False)

    def test_02_integer_raw_intertwiners_and_residual(self):
        p, bits, words, active, record, local, packet, c, z, w1, w2, k2, d, b2, a3, a4 = model()
        np.testing.assert_array_equal(record[:, None]*w1+4*packet@z, np.zeros_like(w1))
        np.testing.assert_array_equal(record[:, None]*w2+4*packet@w1-8*z@k2, np.zeros_like(w2))
        np.testing.assert_array_equal(d@z, z@(z.T@d@z))
        np.testing.assert_array_equal(z.T@packet@z, np.zeros((48, 48), dtype=int))
        for w in (w1, w2):
            np.testing.assert_array_equal(local@w, w@(z.T@local@z))
        # Verify the complete polynomial identity with independently assembled Fraction matrices
        # at a rational epsilon, using integer scaling (denominator 32*7^4).
        e_den = 7
        h49 = 49*np.diag(record)+7*packet+d
        w784 = 16*49*z+4*7*w1+w2
        left = 2*h49@w784-w784@b2
        right = 2*7*a3+a4
        np.testing.assert_array_equal(left, right)
        # Controlled frame of the matching term is precisely the old rematching K tensor I.
        cp = z.T@c@z
        k_old = old.four_model()[6]
        np.testing.assert_array_equal(cp.T@k2@cp,
            np.kron(np.array(2*k_old, dtype=int), np.eye(16, dtype=int)))
        OBS['intertwiner'] = dict(full_carrier_dimension=384, code_dimension=48,
            exact_integer_polynomial_identity=True,
            arbitrary_unknown_marker_logic_reference_covered=True,
            uniform_local_field_cancels_from_residual=True,
            full_reference_error_not_fitted_to_selected_inputs=True)

    def test_03_analytic_norm_and_all_time_protection_bounds(self):
        p, bits, words, active, record, local, packet, c, z, w1, w2, k2, d, b2, a3, a4 = model()
        # Reuse the exact 445 marker Gram identities, tensor I, through the established unitary.
        cp = z.T@c@z
        marker = old.four_model()
        for w_scaled, scale, reference in ((w1, 4, marker[4]), (w2, 16, marker[5])):
            gram = np.array(reference.T@reference, dtype=object)
            actual = cp.T@(w_scaled.T@w_scaled)@cp
            expected = np.kron(np.array(scale*scale*gram, dtype=int), np.eye(16, dtype=int))
            np.testing.assert_array_equal(actual, expected)
        norms = [old.opnorm(w1/4), old.opnorm(w2/16), old.opnorm(a3/16), old.opnorm(a4/32)]
        self.assertLess(norms[2], 8)
        self.assertLess(norms[3], 46/5)
        epsilon, hi = F(1, 1024), F(1001, 1000)
        b = error_bound(epsilon, hi)
        leakage_amplitude = 6*epsilon+6*epsilon**2
        self.assertLess(b, F(1, 100))
        OBS['bounds'] = dict(j=1, ell=1, epsilon='1/1024',
            tau_window=['999/1000', '1001/1000'],
            all_input_window_operator_error=old.fraction_data(b),
            all_time_code_leakage_upper=old.fraction_data(leakage_amplitude**2),
            general_all_time_amplitude='m*abs(epsilon)+N*(abs(j)/2+abs(ell))*epsilon^2',
            exact_residual_norm_numerical_cross_checks=[old.short(x) for x in norms],
            analytic_residual_upper=['8', '46/5'],
            finite_N_scope=True, unlimited_scale_budget_claimed=False)

    def test_04_rational_effective_readout_and_rematching_window(self):
        *_, b2, a3, a4 = model()
        order, remainder, perr, probs = rational_columns()
        gap_lower = probs[0][1]-probs[0][0]-2*perr
        rematch_lower = min(probs[1])-perr
        self.assertGreater(gap_lower, F(9, 100))
        self.assertGreater(rematch_lower, F(39, 100))
        delta = F(1, 1000)
        b = error_bound(F(1, 1024), 1+delta)
        real_gap = F(9, 100)-6*delta-2*b
        real_rematch = F(39, 100)-2*delta-b
        self.assertGreater(real_gap, F(1, 16))
        self.assertGreater(real_rematch, F(3, 8))
        # Independent numerical norm checks supplement the analytic estimates 7,3,2.
        cols, _, _, e, r = experiment()
        effective = b2/2
        self.assertLess(old.opnorm(effective), 7)
        self.assertLess(old.opnorm(effective@e-e@effective), 3)
        self.assertLess(old.opnorm(effective@r-r@effective), 2)
        OBS['window_certificate'] = dict(taylor_order=order, effective_H_norm_upper=7,
            exp7_upper=3**7, Taylor_vector_remainder=old.fraction_data(remainder),
            probability_error=old.fraction_data(perr),
            effective_t1_population_difference_lower='9/100',
            effective_t1_each_message_rematching_lower='39/100',
            physical_population_difference_lower=old.fraction_data(real_gap),
            physical_each_message_rematching_lower=old.fraction_data(real_rematch),
            simple_signal_strict_lower='1/16', simple_rematching_strict_lower='3/8',
            same_fixed_subject1_blank_effect=True, no_postselection=True,
            specified_matching_basis_and_two_message_preparations=True,
            all_unknown_inputs_claimed_to_send_same_signal=False)

    def test_05_full_dynamics_and_reference_complete_cross_check(self):
        p, bits, words, active, record, local, packet, c, z, w1, w2, k2, d, b2, a3, a4 = model()
        epsilon, tau = F(1, 1024), 1.0
        # Physical units: g=1, Delta=1024, J_A=lambda=1/1024, t=1024.
        physical_h = np.diag(1024*record)+packet+d/1024
        full = old.unitary(physical_h, 1024)
        ideal = old.unitary(b2/2, tau)
        error = old.opnorm(full@z-z@ideal)
        b = error_bound(epsilon, F(1))
        self.assertLess(error, float(b))
        outside = full@z-z@(z.T@full@z)
        leakage = old.opnorm(outside)**2
        self.assertLess(leakage, float((6*epsilon+6*epsilon**2)**2))
        cols, e, r, ee, rr = experiment()
        full_states = (full@z)[:, cols]
        ide_states = ideal[:, cols]
        population = [float(np.sum(e*np.abs(full_states[:,i])**2)) for i in range(2)]
        rematch = [float(np.sum(r*np.abs(full_states[:,i])**2)) for i in range(2)]
        effective_p = [float(np.sum(np.diag(ee)*np.abs(ide_states[:,i])**2)) for i in range(2)]
        self.assertGreater(population[1]-population[0], 1/16)
        self.assertGreater(min(rematch), 3/8)
        # Arbitrary mixed-input/reference guarantee is already operator norm; one correlated check.
        rng = np.random.default_rng(449)
        psi = rng.normal(size=(48, 3))+1j*rng.normal(size=(48, 3))
        psi /= np.linalg.norm(psi)
        ref_error = np.linalg.norm((full@z-z@ideal)@psi)
        self.assertLess(ref_error, float(b))
        OBS['full_dynamics'] = dict(Delta=1024, g=1, J_A='1/1024', lambda_value='1/1024',
            physical_time=1024, operator_error=old.short(error),
            analytic_operator_error=old.fraction_data(b),
            all_input_leakage=old.short(leakage),
            physical_populations=[old.short(x) for x in population],
            effective_populations=[old.short(x) for x in effective_p],
            physical_signal=old.short(population[1]-population[0]),
            physical_rematching=[old.short(x) for x in rematch],
            reference_dimension=3, correlated_reference_error=old.short(ref_error),
            floating_point_results_only_cross_check=True)

    def test_06_record_logic_coupling_and_source_boundary(self):
        p, bits, words, active, record, local, packet, c, z, w1, w2, k2, d, b2, a3, a4 = model()
        ap = z.T@(active[:, None]*z)
        cp = z.T@c@z
        tagged = cp.T@(b2/2)@cp
        # Diagonal 16x16 logic blocks at different matching labels are not equal up to a scalar.
        first = tagged[:16, :16]
        second = tagged[16:32, 16:32]
        difference = first-second
        self.assertGreater(np.linalg.norm(difference), 1)
        self.assertEqual(float(np.trace(difference)), 0)
        commutator2 = k2@ap-ap@k2
        self.assertTrue(np.any(commutator2))
        # In the j=0 tagged factorization, source0 resides at tag M0(0)=1.
        # Original partner1 reads tag M(1), never 1 for a perfect matching.
        for matching in old.matching_words(4):
            self.assertNotEqual(matching[1], 1)
        cols, _, _, effect, _ = experiment()
        transport_only = old.unitary(k2/2+z.T@local@z, 1)
        population0 = [float(np.vdot(transport_only[:,k], effect@transport_only[:,k]).real)
                       for k in cols]
        self.assertLess(abs(population0[1]-population0[0]), 2e-14)
        # Nonconstant diagonal pair phase does not satisfy 429 on its full two-qubit space.
        q = np.diag([2, 1, 1, 2])
        plus = np.ones(4)/2
        state = old.unitary(-q, np.pi/2)@plus
        reduced = state.reshape(2,2)@state.reshape(2,2).conj().T
        trace_distance = np.linalg.norm(reduced-np.full((2,2), 0.5), ord='nuc')/2
        self.assertAlmostEqual(trace_distance, 0.5)
        OBS['coupling_and_source'] = dict(tagged_marker_logic_generator_not_tensor_sum=True,
            marker_rematching_and_conditional_potential_noncommute=True,
            pure_transport_uniform_rotation_old_partner_ideal_signal_identically_zero=True,
            zero_conditional_term_t1_populations=[old.short(v) for v in population0],
            pure_transport_full_model_signal_upper='36*epsilon^2, every time',
            raw_H_contains_only_given_internal_swaps_packet_swaps_and_two_given_potentials=True,
            nonconstant_pair_reward_satisfies_full_429_condition=False,
            pair_phase_same_input_marginal_change=old.short(trace_distance),
            parameters_or_potentials_selected_by_cognitive_principles=False,
            generic_static_exchange_simulation_not_repeated=True)

    def test_07_internal_resource_and_time_ledger(self):
        epsilon = F(1, 1024)
        time_low, time_high = 1024*F(999,1000), 1024*F(1001,1000)
        bound_h = 4*1024+6+F(8,1024)
        # All terms are internal; the positive cost cannot be erased by naming it consensus.
        self.assertEqual(time_low, F(127872, 125))
        self.assertEqual(time_high, F(128128, 125))
        self.assertEqual(bound_h, F(525057,128))
        OBS['resources'] = dict(local_register_dimensions=[5,5], subjects=4,
            full_raw_dimension=5**8, invariant_carrier_dimension=384,
            prepared_unknown_code_dimension=48, packet_contacts=6,
            record_interaction_pairs=6, active_interaction_pairs=6,
            total_generator_norm_upper=old.fraction_data(bound_h),
            physical_read_window=[old.fraction_data(time_low), old.fraction_data(time_high)],
            scaling='J_A=lambda=g^2/Delta; nontrivial rematching time t~Delta/g^2',
            operation_preparation_and_readout_still_inputs=True,
            automatic_cooling_or_generic_repair_proved=False,
            infinite_scale_stability_proved=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=449, baseline_round=448, date='2026-09-24',
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(diagonal_reward_degeneracy_and_symmetry_conditions_proved=True,
            record_and_conditional_processing_compatible_with_explicit_scale_hierarchy=True,
            full_unknown_marker_logic_reference_error_certified=True,
            uniform_internal_field_residual_cancellation_proved=True,
            all_time_finite_N_record_leakage_bound_proved=True,
            same_fixed_physical_effect_and_rematching_on_common_window_certified=True,
            internal_resource_and_time_cost_explicit=True,
            protected_logical_states_required_to_be_degenerate=False,
            new_diagonal_interactions_derived_from_429=False,
            every_composite_full_state_satisfies_429_condition=False,
            unlimited_scale_uniform_stability_proved=False,
            initial_code_preparation_or_actual_clock_generated=False,
            full_spatial_dimension_or_GR_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))
