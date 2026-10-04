"""Round 480: correlated preparation gives a strictly certified local 3-readout chart.

Scientific baseline 479. The correlated source/graph preparation, classical
branch record and independently controlled internal exchange remain inputs.
No claim of preparation from the older product seed, a displacement group,
physical three-dimensional space from cognition, or GR.
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
import rigid_leaf_reference_audit as frame
import symmetric_control_reachability_audit as prior
import singlet_supply_capacity_audit as fixed
import sequential_tree_distance_readout_audit as reader

TARGET = Path(__file__).with_name('correlated_reference_local_access_results.json')
OBS = {}
BINV = 15000
ENTRY_ERROR = F(1, 10**14)
TAU = F(1, 10**10)


def short(x):
    return float(f'{float(x):.13g}')


def pure_columns():
    """Gaussian-integer columns C, with rho=C C*/384 and 4+24 columns."""
    columns = []
    for source, graph_modes in ((0, [None]), (1, list(range(6)))):
        for a, b in itertools.product(range(2), repeat=2):
            data = np.zeros(64, dtype=complex)
            for v1, v3, v4 in itertools.product(range(2), repeat=3):
                bits = [a, v1, b, v3, v4, 0]
                index = sum(z << (5-i) for i, z in enumerate(bits))
                data[index] = (1j**v4)*(1j**v1 if source else 1)
            for g in graph_modes:
                graph = np.ones(6) if g is None else np.eye(6)[g]
                columns.append(np.kron(data, graph))
    return np.column_stack(columns)


def integer_columns():
    c = pure_columns()
    # astype(object) alone leaves Python floats; cast through int64 FIRST.
    real = c.real.astype(np.int64).astype(object)
    imag = c.imag.astype(np.int64).astype(object)
    assert all(isinstance(v, int) for v in real.flat)
    assert all(isinstance(v, int) for v in imag.flat)
    return real, imag


def determinant(a):
    return sum(((-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3)))
               *math.prod(a[i][p[i]] for i in range(3))
               for p in itertools.permutations(range(3)))


def inverse3(a):
    d = determinant(a)
    out = []
    for i in range(3):
        row = []
        for j in range(3):
            rs = [v for v in range(3) if v != j]
            cs = [v for v in range(3) if v != i]
            row.append(((-1)**(i+j))*(a[rs[0]][cs[0]]*a[rs[1]][cs[1]]
                       -a[rs[0]][cs[1]]*a[rs[1]][cs[0]])/d)
        out.append(row)
    return out


@lru_cache(None)
def certificate():
    h = frame.system()[1]
    cr, ci = integer_columns()
    q = np.array([np.tile(v, 64) for v in prior.qgraph()], dtype=np.int64)
    swap_indices = np.argmax(prior.s12(), axis=1)
    result = [[0]*3 for _ in range(3)]
    error = F(0)
    scale = 1 << 80
    blocks = []
    for excitation in range(7):
        indices = np.array([6*x+g for x in range(64) if x.bit_count() == excitation
                            for g in range(6)])
        subh = h[np.ix_(indices, indices)]
        hi = subh.astype(object)
        local = {v: j for j, v in enumerate(indices)}
        permutation = np.array([local[swap_indices[v]] for v in indices])
        ur, ui, b, err, x = fixed.dyadic_unitary_certificate(
            subh, time=1, squarings=2, degree=14, bits=80)
        assert b == scale
        assert all(isinstance(v, int) for v in ur.flat)
        assert all(isinstance(v, int) for v in ui.flat)
        error = max(error, err)
        def apply(v):
            return ur@v[0]-ui@v[1], ur@v[1]+ui@v[0]
        def ham(v):
            return hi@v[0], hi@v[1]
        def minus_i(v):
            return v[1], -v[0]
        def swap(v):
            return v[0][permutation], v[1][permutation]
        start = apply((cr[indices], ci[indices]))
        output = apply(apply(swap(start)))
        derivative_t = minus_i(apply(apply(swap(ham(start)))))
        derivative_u = minus_i(ham(output))
        derivative_theta = minus_i(apply(apply(start)))
        for row in range(3):
            weights = q[row, indices, None].astype(object)
            for col, derivative in enumerate((derivative_t, derivative_u, derivative_theta)):
                numerator = np.sum(weights*(derivative[0]*output[0]+derivative[1]*output[1]))
                assert isinstance(numerator, int)
                result[row][col] += 2*numerator
        blocks.append(dict(excitation=excitation, dimension=len(indices),
                           unitary_norm_error=str(err), scaled_generator_norm=str(x)))
    jacobian = [[F(v, 384*scale**6) for v in row] for row in result]
    det = determinant(jacobian)
    inv = inverse3(jacobian)
    inverse_norm = max(sum(abs(x) for x in row) for row in inv)
    w_error = (1+error)**3-1
    entry_error = 18*w_error*(2+w_error)
    return dict(jacobian=jacobian, determinant=det, inverse_inf_norm=inverse_norm,
                unitary_error=error, W_error=w_error, entry_error=entry_error, blocks=blocks)


def exact_budgets():
    beta = 3*BINV*(ENTRY_ERROR+324*TAU)
    radius = F(1, 8*2052*BINV)
    step = radius/2
    lower_lipschitz = F(3, 8*BINV)
    separation = step*lower_lipschitz
    return beta, radius, step, lower_lipschitz, separation


@lru_cache(None)
def numerical_control():
    h = frame.system()[1]
    s = prior.s12()
    t, u, theta = .2, .4, math.pi/2
    right, left = prior.unitary(t*h), prior.unitary(u*h)
    generator = float(TAU)*h+theta*s
    values, vectors = np.linalg.eigh(generator)
    phases = np.exp(-1j*values)
    pulse = (vectors*phases)@vectors.conj().T
    # Stable exact Frechet derivative formula for the Hermitian exponential.
    mid = (values[:, None]+values[None, :])/2
    difference = values[:, None]-values[None, :]
    derivative_kernel = -1j*np.exp(-1j*mid)*np.sinc(difference/(2*math.pi))
    inner_s = vectors.conj().T@s@vectors
    derivative_pulse = vectors@(derivative_kernel*inner_s)@vectors.conj().T
    w = left@pulse@right
    derivatives = [w@(-1j*h), (-1j*h)@w, left@derivative_pulse@right]
    c = pure_columns()
    y = w@c
    q = np.tile(prior.qgraph(), (1, 64))
    jacobian = np.array([[2*np.sum(q[row, :, None]*(dw@c)*y.conj()).real/384
                         for dw in derivatives] for row in range(3)])
    output = y@y.conj().T/384
    graph = output.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)
    return w, jacobian, graph


class Audit(unittest.TestCase):
    def test_01_correlated_seed_and_record(self):
        cr, ci = integer_columns()
        self.assertEqual(cr.shape, (384, 28))
        self.assertEqual(int(np.sum(cr*cr+ci*ci)), 384)
        # Exact Gaussian integer equality with the intended two-branch density.
        real = cr@cr.T+ci@ci.T
        imag = ci@cr.T-cr@ci.T
        nums = frame.source_numerators()
        target = np.kron(nums[0]+nums[1], np.ones((6, 6)))+np.kron(nums[0]+nums[2], np.eye(6))
        self.assertTrue(np.array_equal(2*real, target.real))
        self.assertTrue(np.array_equal(2*imag, target.imag))
        self.assertEqual(int(np.sum(cr[:, :4]**2+ci[:, :4]**2)), 192)
        self.assertEqual(int(np.sum(cr[:, 4:]**2+ci[:, 4:]**2)), 192)
        h = frame.system()[1]
        charges = np.array([x.bit_count() for x in range(64) for _ in range(6)])
        self.assertTrue(np.array_equal(charges[:, None]*h, h*charges[None, :]))
        s = prior.s12()
        self.assertTrue(np.array_equal(charges[:, None]*s, s*charges[None, :]))
        self.assertTrue(np.array_equal(prior.K@np.array([1, 1, 0]), np.array([1, -1, 0])))
        OBS['preparation_contract'] = dict(
            state='one half eta(ex) x |s><s| plus one half eta(ey) x I_G/6',
            internal_branch_record_dimension=2, branch_probabilities=['1/2', '1/2'],
            pure_Gaussian_column_count=28, density_denominator=384,
            source_marginal_Bloch=['1/2', '1/2', '0'],
            graph_marginal='(|s><s|+I_G/6)/2',
            same_control_applies_to_both_branches_no_flag_feedback=True,
            graph_mixed_source_entropy_and_triad_preparation_are_resources=True,
            source_graph_correlations_not_assumed_generated_from_old_product_seed=True,
            conserved_excitation_blocks_used_only_for_computation=True,
            no_physical_dephasing_or_discarding_interblock_coherences=True)

    def test_02_strict_fixed_point_Jacobian(self):
        c = certificate()
        self.assertGreater(c['determinant'], F(1, 2_000_000))
        self.assertLess(c['determinant'], F(3, 5_000_000))
        self.assertLess(c['inverse_inf_norm'], BINV)
        self.assertLess(c['unitary_error'], F(4, 10**17))
        self.assertLess(c['entry_error'], ENTRY_ERROR)
        self.assertLess(3*BINV*ENTRY_ERROR, F(1, 4))
        OBS['strict_Jacobian_certificate'] = dict(
            point=['1/5', '2/5', 'pi/2'], parameter_order=['first_wait_t', 'second_wait_u', 'exchange_area_theta'],
            observable_order=['D13-D10', 'D14-D10', 'D15-D10'],
            rational_Jacobian=[[str(x) for x in row] for row in c['jacobian']],
            decimal_Jacobian=[[short(x) for x in row] for row in c['jacobian']],
            exact_rational_determinant_interval=['1/2000000', '3/5000000'],
            determinant_decimal=short(c['determinant']),
            rational_inverse_inf_norm_less_than=BINV,
            unitary_error=str(c['unitary_error']), entry_error_upper=str(ENTRY_ERROR),
            exact_computed_entry_error=str(c['entry_error']),
            ideal_Neumann_upper=str(3*BINV*ENTRY_ERROR),
            blocks=c['blocks'], fixed_point_bits=80, Taylor_degree=14, squarings=2,
            all_object_arithmetic_inputs_checked_as_Python_integers=True,
            no_floating_SVD_or_determinant_used_for_nonzero_proof=True)

    def test_03_finite_pulse_and_uniform_local_chart(self):
        beta, radius, step, lipschitz, separation = exact_budgets()
        self.assertLess(beta, F(1, 2))
        self.assertEqual(2*BINV*2052*radius, F(1, 4))
        self.assertLess(radius, F(1, 5))
        self.assertEqual(lipschitz, F(1, 40000))
        self.assertEqual(separation, F(1, 19699200000000))
        OBS['finite_control_and_local_chart'] = dict(
            pulse_duration=str(TAU), peak_exchange_strength_near_point='pi/(2*tau)',
            original_H_remains_on=True,
            pulse_unitary_and_theta_derivative_error_upper='9*tau',
            extra_Jacobian_entry_error_upper='324*tau',
            finite_pulse_Neumann_upper=str(beta),
            true_inverse_inf_norm_upper=2*BINV,
            Jacobian_local_Lipschitz_inf_bound=2052,
            parameter_cube_sup_radius=str(radius),
            normalized_Jacobian_variation_upper='1/4',
            readout_lower_Lipschitz_in_sup_norm=str(lipschitz),
            four_settings='p0 and p0+step*e_j for j=1,2,3',
            step=str(step), any_two_settings_readout_separation_lower=str(separation),
            inverse_function_theorem_gives_open_three_readout_neighborhood=True,
            chart_is_local_readout_access_not_a_spatial_displacement_group=True)

    def test_04_full_unitary_and_preserved_branch_record(self):
        w, jacobian, graph = numerical_control()
        c = certificate()
        approximate = np.array([[float(x) for x in row] for row in c['jacobian']])
        self.assertLess(np.max(abs(jacobian-approximate)), float(ENTRY_ERROR+324*TAU)+2e-12)
        self.assertLess(np.max(abs(w.conj().T@w-np.eye(384))), 3e-12)
        columns = pure_columns()
        branch_errors = []
        for subset in (slice(0, 4), slice(4, 28)):
            y = w@columns[:, subset]
            branch_graph = (y@y.conj().T/192).reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)
            branch_errors.append(max(abs(np.diag(branch_graph).real@frame.distance(a, b)-8/3)
                                     for a, b in itertools.combinations(frame.LEAVES, 2)))
        self.assertLess(max(branch_errors), 3e-12)
        rng = np.random.default_rng(480)
        unknown = rng.normal(size=(384, 3))+1j*rng.normal(size=(384, 3))
        unknown /= np.linalg.norm(unknown)
        reverse_error = np.linalg.norm(w.conj().T@(w@unknown)-unknown)
        self.assertLess(reverse_error, 3e-12)
        q = prior.qgraph()@np.diag(graph).real
        self.assertGreater(abs(q[2]), 1e-5)
        OBS['numerical_and_reference_checks'] = dict(
            finite_pulse_Jacobian=[[short(x) for x in row] for row in jacobian],
            finite_pulse_determinant=short(np.linalg.det(jacobian)),
            source_leaf_readouts=[short(x) for x in q],
            fixed_product_marginals_would_require_Qz_zero=True,
            maximum_branch_leaf_mean_errors=[short(x) for x in branch_errors],
            untouched_classical_branch_flag_retains_matching_probabilities_one_third=True,
            unknown_full_data_graph_reference_inverse_error=short(reverse_error),
            mathematical_inverse_not_claimed_as_available_time_reversal=True,
            arbitrary_R_information_preserved_by_unitarity=True,
            conditional_reference_moment_protection_requires_each_branch_contract_475=True,
            no_arbitrary_purification_decorrelation_claim=True)

    def test_05_actual_CP_data_readout(self):
        _, _, graph = numerical_control()
        true = prior.qgraph()@np.diag(graph).real
        probe = F(1, 65536)
        means, rows = [], []
        for leaf in (0, 3, 4, 5):
            minus, plus, _ = reader.instrument((1, leaf), probe)
            probabilities = [float(np.trace(np.einsum('ghij,ij->gh', op, graph)).real)
                             for op in (minus, plus)]
            self.assertLess(abs(sum(probabilities)-1), 3e-12)
            self.assertGreaterEqual(min(probabilities), -3e-12)
            means.append((probabilities[1]-probabilities[0])/float(probe))
            rows.append(dict(ordered_edge=[1, leaf], measured_port=1,
                             record_probabilities=[short(x) for x in probabilities]))
        measured = means[0]-np.array(means[1:])
        x = 18*probe
        bias = x*x/(2*(1-x/3)*probe)
        self.assertLess(np.max(abs(measured-true)), float(2*bias)+2e-9)
        OBS['actual_data_instruments'] = dict(
            protocol_round=472, actual_U_CP_instruments=True, probe_wait=str(probe),
            records=rows, true_Q=[short(x) for x in true],
            instrument_Q=[short(x) for x in measured],
            each_Q_diagnostic_bias_bound=str(2*bias),
            old_correlations_kept_in_internal_isolated_storage_for_final_readout=True,
            diagnostic_not_the_tiny_certified_16_mean_statistics_run=True)

    def test_06_finite_precision_distinguishability_budget(self):
        _, _, _, _, separation = exact_budgets()
        source_error = separation/8
        mean_error = separation/16
        h = 9
        probe = mean_error/(16*h*h)
        gamma = mean_error*probe/4
        x = 2*h*probe
        self.assertLessEqual((x*x/(2*(1-x/3))+gamma)/probe, mean_error/2)
        copies = (80/(probe*probe*mean_error*mean_error))
        copies = (copies.numerator+copies.denominator-1)//copies.denominator
        exponent = copies*probe*probe*mean_error*mean_error/8
        self.assertGreaterEqual(exponent, 10)
        self.assertGreater(sum(F(10**j, math.factorial(j)) for j in range(10)), 3200)
        remaining = separation-2*source_error-4*mean_error
        self.assertEqual(remaining, separation/2)
        OBS['resource_ledger'] = dict(
            each_setting_complete_source_and_control_process_diamond_budget=str(source_error),
            each_distance_mean_estimation_error=str(mean_error),
            probe_wait=str(probe), complete_first_instrument_diamond_budget=str(gamma),
            settings=4, separately_estimated_means=16,
            independent_copies_per_mean=str(copies), total_independent_seed_copies=str(16*copies),
            joint_failure_probability_at_most='1/100',
            retained_pairwise_three_readout_separation=str(remaining),
            full_correlated_seed_and_branch_record_preparation_counted=True,
            preparation_timing_control_precision_storage_and_handoff_errors_counted=True,
            fixed_instruments_no_old_record_feedback=True,
            no_cloning_of_single_unknown_input_or_actual_enormous_sampling_run=True,
            finite_point_experiment_does_not_certify_an_arbitrary_noisy_calibration_derivative=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=480, baseline_round=479, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(strict_local_three_readout_access_from_one_fixed_correlated_seed=True,
            full_global_unitary_and_classical_branch_record_retained=True,
            old_product_preparation_two_coefficient_theorem_not_contradicted=True,
            correlated_seed_not_derived_from_old_product_state_or_cognitive_axioms=True,
            extra_reference_control_and_measurement_contracts_counted=True,
            displacement_group_or_physical_space_not_derived=True,
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
