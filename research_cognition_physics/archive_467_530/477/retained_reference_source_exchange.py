"""Round 477: retaining data/graph correlations during symmetric source exchange.

Scientific baseline 475. General symmetry protection and an actual graph-
distance response under a counted internal exchange control; not a displacement
or a derivation of the controller, three-dimensional space, or GR.
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
import rigid_leaf_reference_audit as ref
import sequential_tree_distance_readout_audit as reader

TARGET = Path(__file__).with_name('retained_reference_source_exchange_results.json')
OBS = {}


@lru_cache(None)
def swap_data():
    mapping = []
    for x in range(64):
        bits = [(x >> (5-j)) & 1 for j in range(6)]
        bits[1], bits[2] = bits[2], bits[1]
        mapping.append(sum(bits[j] << (5-j) for j in range(6)))
    return np.eye(64, dtype=np.int64)[mapping]


def full_swap():
    return np.kron(swap_data(), np.eye(6, dtype=np.int64))


def exp_unitary(generator):
    values, vectors = np.linalg.eigh(generator)
    return (vectors*np.exp(-1j*values))@vectors.conj().T


def qoperators():
    return [np.kron(np.eye(64, dtype=np.int64),
            np.diag(ref.distance(1, j)-ref.distance(1, 0))) for j in (3, 4, 5)]


def comm(h, x):
    return h@x-x@h


@lru_cache(None)
def exact_coefficients():
    _, h, _ = ref.system()
    s = full_swap()
    indices = np.argmax(s, axis=1)
    raw = [np.kron(a, np.ones((6, 6))) for a in ref.source_numerators()]
    result = {'continuous': [], 'swap': [], 'difference': []}
    # Prior Gaussian-integer bound includes trace, affine input, binomials,
    # and subtraction. Every component is exactly represented below 2**53.
    prior_bound = 4*384**2*4*(2*18)**3
    assert prior_bound < 2**53
    for degree in range(4):
        normal_rows, swap_rows, difference_rows = [], [], []
        for q in qoperators():
            derivatives = [q]
            for _ in range(degree):
                derivatives.append(comm(h, derivatives[-1]))
            normal = (2**degree)*derivatives[degree]
            controlled = np.zeros((384, 384), dtype=np.int64)
            for a in range(degree+1):
                value = derivatives[degree-a][np.ix_(indices, indices)]
                for _ in range(a):
                    value = comm(h, value)
                controlled += math.comb(degree, a)*value
            rows = []
            for operator in (normal, controlled, controlled-normal):
                row = []
                for numerator in raw:
                    z = np.einsum('ij,ji->', numerator, operator)*(1j**degree)
                    assert z.imag == 0 and z.real == round(z.real)
                    row.append(F(int(z.real), 384*math.factorial(degree)))
                rows.append(row)
            normal_rows.append(rows[0]); swap_rows.append(rows[1]); difference_rows.append(rows[2])
        result['continuous'].append(normal_rows)
        result['swap'].append(swap_rows)
        result['difference'].append(difference_rows)
    return result, prior_bound


def ceil_fraction(x):
    return (x.numerator+x.denominator-1)//x.denominator


def constants():
    t0, t1 = F(1, 1024), F(1, 512)
    gap = F(2, 3)*t0*t0
    tau = gap/288
    area_error = gap/16
    return t0, t1, gap, tau, area_error


def record_moments(unitary):
    basis = np.eye(2)
    plus = (basis[:, 0]+basis[:, 1])/math.sqrt(2)
    plus_y = (basis[:, 0]+1j*basis[:, 1])/math.sqrt(2)
    moments = [np.zeros((2, 2), dtype=complex) for _ in range(3)]
    masks = [(ref.distance(0, j) == 2) for j in (3, 4, 5)]
    for a, b in itertools.product(range(2), repeat=2):
        columns = [np.kron(ref.tensor([basis[:, a], basis[:, q], basis[:, b],
                   plus, plus_y, basis[:, 0]]), np.ones(6)/math.sqrt(6))/math.sqrt(2)
                   for q in range(2)]
        evolved = (unitary@np.column_stack(columns)).reshape(64, 6, 2)
        for j, mask in enumerate(masks):
            data = evolved[:, mask, :].reshape(-1, 2)
            moments[j] += data.T@data.conj()/4
    return moments


class Audit(unittest.TestCase):
    def test_01_exact_control_symmetries(self):
        _, h, _ = ref.system()
        s = full_swap()
        self.assertTrue(np.array_equal(s@s, np.eye(384)))
        for _, dm, gm, _ in ref.permutations():
            total = np.array([6*d+g for d in dm for g in gm])
            self.assertTrue(np.array_equal(s[np.ix_(total, total)], s))
        for a in range(3):
            rotation = np.kron(sum(ref.pauli_at(v, a) for v in range(6)), np.eye(6))
            self.assertTrue(np.array_equal(s@rotation, rotation@s))
        for q in qoperators():
            self.assertTrue(np.array_equal(s@q, q@s))
        self.assertFalse(np.array_equal(s@h, h@s))
        self.assertTrue(all((1, 2) in t for t in ref.system()[0]))
        OBS['exact_symmetry'] = dict(control_is_existing_always_present_edge_S12=True,
            separate_S4_and_collective_SU2_covariance=True,
            pulse_alone_does_not_change_any_graph_observable=True,
            control_does_not_commute_with_full_H=True,
            all_finite_products_and_integrable_real_schedules_protected=True,
            preservation_applies_to_full_correlated_joint_state=True,
            source_and_jointly_invariant_graph_reference_extensions_from_475=True,
            no_data_reset_or_new_triad_between_segments=True)

    def test_02_exact_actual_distance_response(self):
        coeff, bound = exact_coefficients()
        zero = [[F(0)]*4 for _ in range(3)]
        self.assertEqual(coeff['continuous'][0], zero)
        self.assertEqual(coeff['continuous'][1], zero)
        self.assertEqual(coeff['swap'][0], zero)
        self.assertEqual(coeff['swap'][1], zero)
        self.assertEqual(coeff['swap'][3], zero)
        for j in range(3):
            for mu in range(4):
                diagonal = (mu == j+1)
                self.assertEqual(coeff['continuous'][2][j][mu], F(-8, 3) if diagonal else 0)
                self.assertEqual(coeff['swap'][2][j][mu], F(-4, 3) if diagonal else 0)
                self.assertEqual(coeff['difference'][2][j][mu], F(4, 3) if diagonal else 0)
        self.assertEqual(coeff['difference'][3][0][1], 0)
        self.assertNotEqual(coeff['continuous'][3], zero)
        OBS['exact_response'] = dict(
            affine_source_order=['constant', 'x', 'y', 'z'],
            observable_order=['D13-D10', 'D14-D10', 'D15-D10'],
            Taylor_coefficients={key: [[[str(v) for v in row] for row in degree]
                for degree in array] for key, array in coeff.items()},
            Gaussian_integer_prior_bound=bound,
            specific_Qx_source_x_cubic_zero=True,
            no_general_time_evenness_claim=True)

    def test_03_rational_window_and_finite_pulse(self):
        t0, t1, gap, tau, area_error = constants()
        # ||H||<=9, both free legs total Taylor superoperator norm<=36.
        # After the exactly vanishing cubic witness coefficient, two tails:
        # |R|<=2*(36t)^4/[24*(1-36t/5)].
        coefficient_tail = F(2*36**4, 24)*t1*t1/(1-F(36, 5)*t1)
        self.assertLess(coefficient_tail, F(2, 3))
        self.assertEqual(gap, F(1, 1572864))
        pulse_error = 36*tau+2*area_error
        self.assertEqual(pulse_error, gap/4)
        self.assertEqual(F(2, 3)*t0*t0, gap)
        OBS['rational_certificate'] = dict(free_leg_window=[str(t0), str(t1)],
            normalized_fourth_order_tail_upper=str(coefficient_tail),
            ideal_response_bounds=['2*t^2/3', '2*t^2'],
            uniform_ideal_gap=str(gap), finite_pulse_duration_upper=str(tau),
            pulse_exchange_area_error_upper=str(area_error),
            pulse_area_target='pi/2', constant_pulse_amplitude='(pi/2+delta_theta)/tau',
            original_H_remains_on_during_pulse=True,
            matched_total_time_reference='U(2*t+tau)',
            finite_pulse_and_matched_time_error=str(pulse_error),
            finite_pulse_response_lower=str(3*gap/4),
            finite_control_not_derived_autonomously=True)

    def test_04_full_unitary_control_witness(self):
        _, h, _ = ref.system()
        s = full_swap()
        sg = np.ones((6, 6))/6
        rho = np.kron(ref.eta((1., 0., 0.)), sg)
        q = qoperators()[0]
        t0, t1, gap, tau, area_error = constants()
        records = []
        worst_leaf_error = 0.
        for t_exact in (t0, t1):
            t = float(t_exact)
            u = exp_unitary(t*h)
            ideal = u@s@u
            pulse = exp_unitary(float(tau)*h+(math.pi/2+float(area_error))*s)
            actual = u@pulse@u
            reference = exp_unitary((2*t+float(tau))*h)
            def state(w):
                return w@rho@w.conj().T
            ideal_difference = np.trace(q@(state(ideal)-state(u@u))).real
            actual_difference = np.trace(q@(state(actual)-state(reference))).real
            self.assertGreater(ideal_difference, float(F(2, 3)*t_exact*t_exact)-2e-13)
            self.assertLess(ideal_difference, float(2*t_exact*t_exact)+2e-13)
            self.assertGreater(actual_difference, float(3*gap/4)-2e-13)
            population = np.diag(state(actual)).real.reshape(64, 6).sum(axis=0)
            error = max(abs(population@ref.distance(a, b)-8/3)
                        for a, b in itertools.combinations(ref.LEAVES, 2))
            worst_leaf_error = max(worst_leaf_error, error)
            records.append(dict(free_leg_time=str(t_exact),
                ideal_distance_difference=float(ideal_difference),
                finite_pulse_distance_difference=float(actual_difference),
                leaf_mean_error=float(error)))
        self.assertLess(worst_leaf_error, 2e-12)
        OBS['full_unitary_witness'] = dict(samples=records,
            actual_384_dimensional_distance_observable=True,
            numerical_check_not_replacement_for_rational_window=True)

    def test_05_unknown_information_and_multistep_protection(self):
        _, h, _ = ref.system()
        s = full_swap()
        total = np.eye(384, dtype=complex)
        schedule = [(.13, .27), (.21, -.19), (.08, .35)]
        for duration, theta in schedule:
            total = exp_unitary(duration*h+theta*s)@total
        self.assertLess(np.max(abs(total.conj().T@total-np.eye(384))), 2e-12)
        rng = np.random.default_rng(477)
        unknown = rng.normal(size=(384, 3))+1j*rng.normal(size=(384, 3))
        unknown /= np.linalg.norm(unknown)
        outgoing = total@unknown
        recovered = total.conj().T@outgoing
        error = np.linalg.norm(recovered-unknown)
        self.assertLess(error, 2e-12)
        # Actual arbitrary source-reference entanglement test in protected prep.
        moments = record_moments(total)
        reference_error = max(np.max(abs(m-np.eye(2)/6)) for m in moments)
        self.assertLess(reference_error, 2e-12)
        OBS['complete_information'] = dict(schedule_duration_and_exchange_area=schedule,
            fully_unknown_data_graph_reference_dimension=1152,
            inverse_reconstruction_error=float(error),
            Bell_source_reference_matching_moment_error=float(reference_error),
            mathematical_inverse_not_claimed_to_be_physically_available=True,
            no_subsystem_discarded_between_controls=True,
            no_claim_six_graph_populations_or_graph_coherences_constant=True)

    def test_06_counted_data_only_readout(self):
        _, _, gap, _, _ = constants()
        source_process_diamond = gap/16
        epsilon = gap/32
        h = 9
        probe = epsilon/(16*h*h)
        gamma = epsilon*probe/4
        x = 2*h*probe
        remainder = x*x/(2*(1-x/3))
        self.assertLessEqual(remainder/probe, epsilon/4)
        self.assertLessEqual((remainder+gamma)/probe, epsilon/2)
        copies = ceil_fraction(64/(probe*probe*epsilon*epsilon))
        exponent = copies*probe*probe*epsilon*epsilon/8
        self.assertGreaterEqual(exponent, 8)
        e8_lower = sum(F(8**k, math.factorial(k)) for k in range(9))
        self.assertGreater(e8_lower, 800)
        # Four means, each bounded [1,2], have centered norm 1/2.
        source_error = 4*source_process_diamond/2
        record_error = 4*epsilon
        retained_gap = 3*gap/4-source_error-record_error
        self.assertEqual(retained_gap, gap/2)
        # Diagnostic evaluation of the actual original-U CP instrument.
        # This moderate probe time is not the tiny statistical proof budget.
        _, h_original, _ = ref.system()
        t0, _, _, tau, area_error = constants()
        free = exp_unitary(float(t0)*h_original)
        pulse = exp_unitary(float(tau)*h_original+(math.pi/2+float(area_error))*full_swap())
        controlled = free@pulse@free
        baseline = exp_unitary(float(2*t0+tau)*h_original)
        initial = np.kron(ref.eta((1., 0., 0.)), np.ones((6, 6))/6)
        diagnostic_probe = F(1, 65536)
        diagnostic_x = 18*diagnostic_probe
        diagnostic_bound = diagnostic_x**2/(2*(1-diagnostic_x/3)*diagnostic_probe)
        diagnostic = []
        for name, w in (('finite_pulse', controlled), ('matched_time_continuous', baseline)):
            joint = w@initial@w.conj().T
            # This trace computes statistics after a counted storage handoff;
            # it does not authorize physical deletion of the correlated data.
            graph = joint.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)
            estimates, truths, edge_rows = [], [], []
            for leaf in (3, 0):
                edge = (1, leaf)  # Preserve the measured first port 1.
                minus, plus, _ = reader.instrument(edge, diagnostic_probe)
                probabilities = [float(np.trace(np.einsum('ghij,ij->gh', a, graph)).real)
                                 for a in (minus, plus)]
                self.assertLess(abs(sum(probabilities)-1), 2e-12)
                self.assertGreaterEqual(min(probabilities), -2e-12)
                adjacency_estimate = (probabilities[1]-probabilities[0])/float(diagnostic_probe)
                estimate = 2-adjacency_estimate
                truth = float(np.diag(graph).real@ref.distance(1, leaf))
                self.assertLessEqual(abs(estimate-truth), float(diagnostic_bound)+2e-9)
                estimates.append(estimate); truths.append(truth)
                edge_rows.append(dict(ordered_probe_edge=list(edge), measured_port=1,
                    record_probabilities_minus_plus=probabilities,
                    actual_instrument_distance=estimate, graph_mean_distance=truth,
                    absolute_bias=abs(estimate-truth)))
            diagnostic.append(dict(protocol=name, edge_readings=edge_rows,
                Qx_from_actual_data_instruments=estimates[0]-estimates[1],
                true_Qx=truths[0]-truths[1]))
        OBS['data_only_readout'] = dict(protocol_round=472,
            actual_CP_diagnostic_probe=str(diagnostic_probe),
            actual_CP_diagnostic_each_distance_bias_bound=str(diagnostic_bound),
            actual_CP_instrument_diagnostics=diagnostic,
            endpoint_identity='fixed internal label1 and fixed leaf labels0,3',
            graph_identity='D1j=2I-n1j', measured_edges=[[1, 3], [1, 0]],
            source_preparation_and_remaining_timing_process_diamond_budget=str(source_process_diamond),
            per_mean_total_estimation_error=str(epsilon),
            probe_wait=str(probe), complete_instrument_diamond_budget=str(gamma),
            copies_per_mean=str(copies), total_independent_source_trials=str(4*copies),
            union_failure_probability_upper='1/100',
            finite_pulse_then_source_process_then_readout_gap=str(retained_gap),
            complete_source_process_error_includes_preparation_natural_clock_and_unlisted_controls=True,
            old_full_data_graph_correlations_stored_before_final_readout=True,
            handoff_duration_and_error_in_complete_first_instrument_budget=True,
            fresh_probe_clock_isolation_and_fixed_instruments_counted=True,
            no_actual_enormous_sampling_run=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=477, baseline_round=475, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(retained_correlations_support_arbitrary_finite_symmetric_control_sequences=True,
            same_original_interaction_type_but_added_control_contract=True,
            no_autonomous_controller_generation_claim=True,
            actual_source_leaf_distance_response=True,
            source_label_is_not_relabelled_when_internal_data_are_swapped=True,
            conditional_finite_time_data_only_readout_with_resource_budget=True,
            reference_shape_is_equal_time_distribution_not_classical_rigid_frame=True,
            no_replayable_endpoint_displacement_group_derived=True,
            three_dimensional_space_derived=False,
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
