"""Round 456: one always-on converter and reader, with an exact compatibility boundary.

Reuse frozen 430, 431 and 455. No gate scheduling or free auxiliary reset.
Integer coefficient Gram and rational Taylor certificates precede numerical checks.
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
import exchange_relation_audit as old
import relational_interface_conversion_audit as converter

TARGET = Path(__file__).with_name('continuous_interface_join_audit_results.json')
OBS = {}
A, CARRIER, OUT = (0, 1, 2), 4, (3, 5, 6)
PAIRS = list(itertools.combinations(range(7), 2))
EDGES = [(0, 4), (1, 4), (2, 4), (0, 1), (1, 2), (2, 5), (5, 6)]


def short(x):
    return float(f'{float(x):.12g}')


@lru_cache(None)
def model():
    swaps = {pair: old.swap(7, *pair) for pair in PAIRS}
    c = sum(swaps[a, CARRIER] for a in A)
    h = sum(swaps[pair] for pair in EDGES)
    vin = np.kron(converter.model()['vin'], converter.model()['s'][:, None])
    effect = (np.eye(128)-swaps[5, 6])/2
    pout = np.kron(converter.model()['pout'], np.eye(4))
    return swaps, c, h, vin, effect, pout


def coefficient_constraints():
    rows = []
    for a, b in itertools.combinations(A, 2):
        row = np.zeros(21, dtype=np.int64)
        row[PAIRS.index((a, CARRIER))] = 1
        row[PAIRS.index((b, CARRIER))] = -1
        rows.append(row)
    for r in OUT:
        for a in A:
            row = np.zeros(21, dtype=np.int64)
            row[PAIRS.index(tuple(sorted((a, r))))] = 1
            row[PAIRS.index(tuple(sorted((CARRIER, r))))] = -1
            rows.append(row)
    return np.array(rows)


def compatible_basis():
    groups = [[pair] for pair in itertools.combinations(A, 2)]
    groups += [[(a, CARRIER) for a in A]]
    groups += [[tuple(sorted((a, r))) for a in (*A, CARRIER)] for r in OUT]
    groups += [[pair] for pair in itertools.combinations(OUT, 2)]
    return np.array([[int(pair in group) for pair in PAIRS] for group in groups]).T


@lru_cache(None)
def certificate(order=100):
    # rho_plus-rho_minus = i D/(8 sqrt(3)); E = E2/2.
    a, b, c, *_ = old.operators()
    s = old.swap(2, 0, 1).real
    difference = np.kron(np.kron((a@b-b@a).real, np.eye(4)-s), np.eye(4)-s).astype(int).astype(object)
    rho48 = np.kron(np.kron(3*np.eye(8)-(a+b+c).real, np.eye(4)-s), np.eye(4)-s).astype(int).astype(object)
    swaps, _, _, _, effect, pout = model()
    effect2 = (2*effect.real).astype(int).astype(object)
    # Four-spin singlet projector = (J^2-2 I)(J^2-6 I)/12, J^2=sum_{i<j} Sij.
    j2 = sum(old.swap(4, i, j).real for i, j in itertools.combinations(range(4), 2))
    code12 = np.kron((j2-2*np.eye(16))@(j2-6*np.eye(16)), np.eye(8)).astype(int).astype(object)
    assert np.linalg.norm(np.array(code12, float)/12-pout) < 1e-12
    permutations = [np.argmax(swaps[p].real, axis=1) for p in EDGES]
    read_sum, code_sum, traces = F(0), F(0), []
    for k in range(order+1):
        n = int(np.sum(difference.T*effect2))
        m = int(np.sum(rho48.T*code12))
        traces.append(n)
        if k % 2:
            read_sum += F((-1)**((k+1)//2)*n, 16*math.factorial(k))*F(3, 2)**k
            assert m == 0
        else:
            assert n == 0
            code_sum += F((-1)**(k//2)*m, 576*math.factorial(k))*F(3, 2)**k
        effect2 = sum(effect2[p, :]-effect2[:, p] for p in permutations)
        code12 = sum(code12[p, :]-code12[:, p] for p in permutations)
    tail = F(3**21*21**(order+1), math.factorial(order+1))
    lo, hi = F(1732050807568877, 10**15), F(1732050807568878, 10**15)
    assert lo**2 < 3 < hi**2 and read_sum < 0
    lower, upper = -read_sum/hi-2*tail, -read_sum/lo+2*tail
    assert lower > F(1, 3)
    assert code_sum+tail < F(3, 5)
    return dict(order=order, time='3/2', low_order_integer_reader_traces=traces[:10],
        signed_reader_sum_before_sqrt3=str(read_sum), code_probability_sum=str(code_sum),
        single_probability_tail=str(tail), single_probability_tail_float=float(tail),
        contrast_interval=[str(lower), str(upper)],
        contrast_interval_float=[float(lower), float(upper)],
        maximally_mixed_GL_code_probability_interval=[str(code_sum-tail), str(code_sum+tail)],
        code_probability_interval_float=[float(code_sum-tail), float(code_sum+tail)],
        exact_contrast_exceeds_one_third=True, exact_code_probability_below_three_fifths=True)


class Audit(unittest.TestCase):
    def close(self, a, b, tolerance=2e-11):
        self.assertLess(float(np.linalg.norm(a-b)), tolerance)

    def test_01_all_weight_exact_commutator_gram(self):
        swaps, c, *_ = model()
        columns = np.array([(c@s-s@c).real.astype(np.int64).reshape(-1) for s in swaps.values()]).T
        gram = columns.T@columns
        constraints = coefficient_constraints()
        self.assertTrue(np.array_equal(gram, 192*constraints.T@constraints))
        basis = compatible_basis()
        self.assertTrue(np.array_equal(columns@basis, np.zeros((128**2, 10), dtype=int)))
        self.assertEqual(np.linalg.matrix_rank(constraints), 11)
        self.assertEqual(np.linalg.matrix_rank(basis), 10)
        # Exact rank: three star-difference rows have rank 2; each of the three
        # disjoint outside four-edge blocks has three independent differences.
        self.assertEqual(2+3*3, 11)
        OBS['commutation_classification'] = dict(raw_spins=7, independent_pair_weights=21,
            exact_integer_gram=gram.tolist(), norm_factor=192,
            constraint_rank=11, kernel_dimension=10,
            exact_gram_equals_192_BtB=True, kernel_basis=basis.tolist(),
            analytical_general_N_factor='3*2^(N-1)', full_raw_commutation_is_extra_contract=True)

    def test_02_entire_compatible_family_cannot_export_L(self):
        swaps, c, *_ = model()
        e = np.kron(old.encoding(), np.eye(16))  # G,L,outside raw 3,4,5,6
        for coefficients in compatible_basis().T:
            h = c+sum(w*swaps[p] for w, p in zip(coefficients, PAIRS))
            reduced = e.conj().T@h@e
            self.close(h@e, e@reduced)
            # Reorder to L,(G,3,4,5,6), and verify additive factorization.
            k = reduced.reshape(2, 2, 16, 2, 2, 16).transpose(1, 0, 2, 4, 3, 5).reshape(64, 64)
            l = old.partial(k, [2, 32], (0,))/32
            rest = old.partial(k, [2, 32], (1,))/2-np.trace(k)*np.eye(32)/64
            self.close(k, np.kron(l, np.eye(32))+np.kron(np.eye(2), rest))
        OBS['compatible_family'] = dict(all_ten_kernel_generators_checked=True,
            original_A_code_invariant=True,
            all_real_combinations_factorize_as='h_L tensor I + I_L tensor h_G_rest',
            outside_output_depends_only_on_initial_reduced_G_rest_R=True,
            new_L_information_cannot_reach_external_reader=True,
            arbitrary_other_noncommuting_dynamics_excluded=False)

    def test_03_one_fixed_joint_generator_and_reference_safe_effect(self):
        swaps, c, h, vin, effect, _ = model()
        self.close(vin.conj().T@vin, np.eye(4))
        self.assertGreater(np.linalg.norm(c@(h-c)-(h-c)@c), 1)
        generator_residuals, effect_residuals = [], []
        for mu in old.PAULI:
            joint = sum(converter.local(7, i, mu) for i in range(7))
            self.close(joint@vin, vin@np.kron(mu, np.eye(2)))
            self.close(joint@h, h@joint)
            self.close(joint@effect, effect@joint)
            generator_residuals.append(short(np.linalg.norm(joint@vin-vin@np.kron(mu, np.eye(2)))))
        for t in (0., .75, 1.5):
            v = old.evolve(h, t)@vin
            f = v.conj().T@effect@v
            f_l = old.partial(f, [2, 2], (1,))/2
            self.close(f, np.kron(np.eye(2), f_l))
            effect_residuals.append(short(np.linalg.norm(f-np.kron(np.eye(2), f_l))))
        OBS['always_on_joint_model'] = dict(edges=[list(p) for p in EDGES], weights=[1]*7,
            old_converter_kept_on=True, old_reader_chain_reused_on_sites=[0, 1, 2, 5, 6],
            initial_independent_singlets=[[3, 4], [5, 6]],
            full_raw_dimension=128, input_order='G,L', all_unknown_GL_R_inputs_allowed=True,
            collective_intertwiner_residuals=generator_residuals,
            compressed_effect_factorization_residuals=effect_residuals,
            Schur_identity='Vin^* U^* E U Vin = I_G tensor F_L for every t',
            classical_outcome_with_reference_depends_only_on_rho_LR=True,
            entire_reader_density_independent_of_G_claimed=False)

    def test_04_exact_signal_and_finite_window(self):
        cert = certificate()
        self.assertEqual(cert['low_order_integer_reader_traces'], [0, 0, 0, 0, 0, 1440, 0, 59136, 0, 2024352])
        _, _, h, vin, effect, _ = model()
        y = old.PAULI[1]
        rho = [vin@np.kron(np.eye(2)/2, (np.eye(2)+sign*y)/2)@vin.conj().T for sign in (1, -1)]
        probabilities = []
        u = old.evolve(h, 1.5)
        outputs = [u@r@u.conj().T for r in rho]
        for r in outputs:
            probabilities.append(float(np.trace(effect@r).real))
        gap = probabilities[1]-probabilities[0]
        self.assertLess(abs(gap-sum(cert['contrast_interval_float'])/2), 1e-12)
        readers = [old.partial(r, [2]*7, (5, 6)) for r in outputs]
        self.assertAlmostEqual(old.distance(*readers), gap)
        derivative_norm = float(np.linalg.norm(h@effect-effect@h, 2))
        self.assertAlmostEqual(derivative_norm, math.sqrt(3)/2)
        self.assertEqual(F(1, 3)-2*F(1, 24), F(1, 4))
        OBS['exact_certificate'] = cert
        OBS['actual_read_window'] = dict(centre_singlet_probabilities=list(map(short, probabilities)),
            centre_contrast=short(gap), interval=['35/24', '37/24'],
            contrast_strict_lower_bound='1/4', equal_prior_success_strict_lower_bound='5/8',
            single_effect_derivative_norm=short(derivative_norm), exact_stop_needed=False,
            postselection=False, permanent_readable_record_claimed=False)

    def test_05_old_conversion_guarantee_is_not_inherited(self):
        _, c, h, vin, _, pout = model()
        u = old.evolve(h, 1.5)
        f = vin.conj().T@u.conj().T@pout@u@vin
        prob = float(np.trace(f).real/4)
        self.assertLess(abs(prob-sum(certificate()['code_probability_interval_float'])/2), 1e-12)
        old_u = old.evolve(c, 1.5)
        old_prob = float(np.trace(vin.conj().T@old_u.conj().T@pout@old_u@vin).real/4)
        self.assertAlmostEqual(old_prob, 1-3*math.cos(1.5)**2/4)
        # 3 < pi < 22/7 => 0 < pi/2-3/2 < 1/14.
        self.assertGreater(1-F(3, 4)*F(1, 14)**2, F(99, 100))
        OBS['conversion_boundary'] = dict(joint_target_code_probability=short(prob),
            joint_target_code_effect_eigenvalues=list(map(short, np.linalg.eigvalsh(f))),
            isolated_converter_code_probability=short(old_prob),
            exact_joint_probability_upper='3/5', exact_isolated_probability_lower='99/100',
            half_diamond_distance_to_455_ideal_strict_lower='2/5',
            half_diamond_distance_to_isolated_at_same_time_strict_lower='39/100',
            input_witness='rho_GL=I4/4 and two independent singlets',
            this_generator_not_all_noncommuting_generators=True)

    def test_06_global_information_and_internal_resource_accounts(self):
        _, _, h, vin, _, _ = model()
        v = old.evolve(h, 1.5)@vin
        self.close(v.conj().T@v, np.eye(4))
        rho = old.density(12, np.random.default_rng(456))
        joint = np.kron(v, np.eye(3))
        output = joint@rho@joint.conj().T
        recovered = joint.conj().T@output@joint
        self.close(recovered, rho)
        self.close(old.partial(output, [128, 3], (1,)), old.partial(rho, [4, 3], (1,)))
        OBS['information_and_resources'] = dict(full_output_isometric_for_all_unknown_references=True,
            dimension_three_reference_recovery_residual=short(np.linalg.norm(recovered-rho)),
            raw_input_subject_qubits=3, initially_prepared_auxiliary_qubits=4,
            independent_pure_singlet_pairs_supplied=2,
            old_carrier_assignment_or_old_local_L_state_preserved_claimed=False,
            resources_returned_or_free_reset_claimed=False,
            preparations_connectivity_effect_access_and_time_units_are_inputs=True)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=456, baseline_round=455, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors), python=platform.python_version(),
        numpy=np.__version__, observations=OBS,
        scope=dict(all_static_pair_commuting_addons_classified=True,
            no_new_L_export_under_extra_commutation_contract=True,
            single_always_on_noncommuting_reader_positive_example=True,
            unknown_GL_reference_preserved_globally=True,
            reference_safe_classical_effect_on_L=True,
            exact_rational_signal_and_conversion_failure_certificates=True,
            finite_read_window_proved=True, singlets_and_contact_table_are_inputs=True,
            old_455_conversion_guarantee_preserved=False, autonomous_handoff_completed=False,
            full_joint_macro_subject_implemented=False, universal_noncommuting_no_go=False,
            physical_positions_or_three_dimensions_derived=False,
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
    print(json.dumps(dict(round=result['round'], tests=result['tests_run'], failures=result['failures'], errors=result['errors'])))
