"""Round 461: independently unknown subjects do not close in one spin-half port.

Baseline 459; independent of round 460. The specified target includes all
helper degrees of freedom in the new subject. The result is not a ban on
composition: a larger coherent representation interface is explicitly kept.
"""
import argparse
from fractions import Fraction
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
import recursive_exchange_interface_audit as recursive

TARGET = Path(__file__).with_name('independent_subject_admission_audit_results.json')
OBS = {}
TRIPLET = np.array([[1, 0, 0], [0, 1/math.sqrt(2), 0],
                    [0, 1/math.sqrt(2), 0], [0, 0, 1]], complex)
SINGLET = (np.eye(4)-old.swap(2, 0, 1))/2
SPIN_ONE = [TRIPLET.conj().T@recursive.collective(2, mu)@TRIPLET for mu in range(3)]


def short(value):
    return float(f'{float(value):.12g}')


@lru_cache(None)
def integer_casimir(n):
    # 4 J^2 = 4 sum_ij SWAP_ij + n(4-n) I, entirely integer.
    return 4*sum(old.swap(n, a, b).real.astype(np.int64)
                 for a in range(n) for b in range(a+1, n))+n*(4-n)*np.eye(2**n, dtype=np.int64)


@lru_cache(None)
def exact_sector(n, twice_j):
    if twice_j not in range(n % 2, n+1, 2):
        return np.zeros((2**n, 2**n), dtype=np.int64), 1
    identity = np.eye(2**n, dtype=np.int64)
    target = twice_j*(twice_j+2)
    numerator = identity.copy()
    denominator = 1
    for twice_k in range(n % 2, n+1, 2):
        if twice_k != twice_j:
            other = twice_k*(twice_k+2)
            numerator = numerator@(integer_casimir(n)-other*identity)
            denominator *= target-other
    return numerator, denominator


def sector(n, twice_j):
    numerator, denominator = exact_sector(n, twice_j)
    return numerator.astype(complex)/denominator


def effect_direct(n, tau):
    p = sector(n+2, 1).reshape(4, 2**n, 4, 2**n)
    return np.einsum('iajb,ba->ij', p, tau)


def effect_moments(n, tau):
    p_half, p_three = sector(n, 1), sector(n, 3)
    j = [recursive.collective(n, mu) for mu in range(3)]
    q_half = float(np.trace(p_half@tau).real)
    q_three = float(np.trace(p_three@tau).real)
    mean_half = [np.trace(p_half@j[mu]@tau) for mu in range(3)]
    mean_three = [np.trace(p_three@j[mu]@tau) for mu in range(3)]
    second_three = [[np.trace(p_three@j[mu]@j[nu]@tau)
                     for nu in range(3)] for mu in range(3)]
    triplet = q_half*np.eye(3)/3-2*sum(SPIN_ONE[mu]*mean_half[mu] for mu in range(3))/3
    triplet = triplet+(sum(SPIN_ONE[mu]@SPIN_ONE[nu]*second_three[mu][nu]
                           for mu in range(3) for nu in range(3))
                       -.5*sum(SPIN_ONE[mu]*mean_three[mu] for mu in range(3))
                       -1.5*q_three*np.eye(3))/6
    return q_half*SINGLET+TRIPLET@triplet@TRIPLET.conj().T, q_half, q_three


def raw_old_embedding():
    # Order G_A,L_A,G_B,L_B,helper; all 32 columns retained.
    return np.kron(np.kron(old.encoding(), old.encoding()), np.eye(2))


def running_hamiltonian():
    # Old A processing 01,12 and old B processing 34,45 remain present.
    edges = [(i, i+1) for i in range(6)]+[(2, 6)]
    return sum(old.swap(7, a, b) for a, b in edges), edges


class Audit(unittest.TestCase):
    def close(self, left, right, tolerance=8e-11):
        self.assertLess(float(np.linalg.norm(left-right)), tolerance)

    def test_01_general_auxiliary_effect_including_sector_coherence(self):
        rows = []
        for n in (1, 3, 5):
            rng = np.random.default_rng(46100+n)
            tau = old.density(2**n, rng)
            computed, q_half, q_three = effect_moments(n, tau)
            direct = effect_direct(n, tau)
            self.close(computed, direct)
            dephased = sum(sector(n, j)@tau@sector(n, j)
                           for j in range(n % 2, n+1, 2))
            self.close(effect_direct(n, dephased), direct)
            coherence = float(np.linalg.norm(tau-dephased))
            if n >= 3:
                self.assertGreater(coherence, .01)
            self.assertGreater(np.linalg.eigvalsh(direct).min(), -1e-11)
            self.assertLess(np.linalg.eigvalsh(direct).max(), 1+1e-11)
            rows.append(dict(auxiliary_qubits=n, spin_half_weight=short(q_half),
                spin_three_halves_weight=short(q_three),
                cross_sector_coherence_norm=short(coherence)))
        tau_even = old.density(4, np.random.default_rng(46102))
        self.close(effect_direct(2, tau_even), np.zeros((4, 4)))
        OBS['general_effect'] = dict(numerical_auxiliaries=rows,
            moment_formula_covers_arbitrary_finite_SU2_representations=True,
            arbitrary_cross_j_and_multiplicity_coherences_allowed=True,
            acceptance_effect_insensitive_to_cross_auxiliary_j_coherence=True,
            parity_mismatch_gives_zero=True,
            raw_helper_preparation_not_assumed_invariant=True)

    def test_02_exact_minimax_and_attaining_auxiliary(self):
        aligned = []
        for pauli in old.PAULI:
            for sign in (-1, 1):
                rho = (np.eye(2)+sign*pauli)/2
                aligned.append(np.kron(rho, rho))
        self.close(sum(aligned)/6, (np.eye(4)-SINGLET)/3)
        rows = []
        for n in (1, 3, 5):
            tau = old.density(2**n, np.random.default_rng(46120+n))
            effect, q_half, q_three = effect_moments(n, tau)
            probabilities = [float(np.trace(effect@rho).real) for rho in aligned]
            average = q_half/3+q_three/6
            self.assertAlmostEqual(sum(probabilities)/6, average)
            self.assertLessEqual(min(probabilities), average+1e-12)
            self.assertLessEqual(average, 1/3+1e-12)
            rows.append(dict(auxiliary_qubits=n,
                exact_design_average=short(average),
                one_of_six_independent_aligned_inputs_probability=short(min(probabilities))))
        optimal = effect_direct(1, np.eye(2)/2)
        self.close(optimal, SINGLET+(np.eye(4)-SINGLET)/3)
        self.close(np.linalg.eigvalsh(optimal), [1/3, 1/3, 1/3, 1])
        # A fixed pure helper can attain the same bound without leaving a
        # purification carrier outside: its G is Bell-correlated with its
        # own multiplicity L inside the three actual auxiliary qubits.
        pure_helper = (old.encoding()[:, 0]+old.encoding()[:, 3])/math.sqrt(2)
        pure_tau = np.outer(pure_helper, pure_helper.conj())
        self.close(effect_direct(3, pure_tau), optimal)
        self.assertAlmostEqual(float(np.trace(pure_tau@pure_tau).real), 1)
        self.assertEqual(Fraction(1, 3)-Fraction(1, 6), Fraction(1, 6))
        # Pure auxiliary references perform worse for their aligned input;
        # for Bloch length r the exact worst product-input value is (1-r)/3.
        rho_up = np.diag([1, 0])
        for radius in (0, Fraction(1, 2), 1):
            tau = (np.eye(2)+float(radius)*old.PAULI[2])/2
            probability = np.trace(effect_direct(1, tau)@np.kron(rho_up, rho_up)).real
            self.assertAlmostEqual(probability, float((1-radius)/3))
        OBS['minimax'] = dict(value='1/3', fixed_auxiliary_may_have_arbitrary_finite_size=True,
            witness_inputs_are_independent_and_parallel=True,
            finite_six_direction_design_is_exact=True, sampled_auxiliary_checks=rows,
            attaining_auxiliary='one qubit I/2',
            pure_attaining_auxiliary='three qubits V3 (|G0,L0>+|G1,L1>)/sqrt(2)',
            pure_attainment_retains_its_purification_inside_helper=True,
            attainment_is_sector_weight_not_unknown_state_preserving_success=True,
            attaining_effect='P_s + P_t/3',
            independent_input_probability='1/2 - dot(r_A,r_B)/6',
            independent_probability_range=['1/3', '2/3'],
            arbitrary_joint_input_probability_range=['1/3', '1'],
            auxiliary_qubit_bloch_radius_worst_case='(1-r)/3',
            upper_bound_is_general_analytic_not_a_state_scan=True)

    def test_03_closed_exchange_conservation_and_support_error(self):
        numerator, denominator = exact_sector(7, 1)
        self.assertEqual(denominator, -23040)
        self.assertEqual(int(np.trace(numerator)), 28*denominator)
        for a in range(7):
            for b in range(a+1, 7):
                perm = recursive.permutation(7, a, b)
                self.assertTrue(np.array_equal(numerator[perm, :], numerator[:, perm]))
        # Two independently prepared old subjects G=0,L=0, plus I/2 helper.
        a = old.encoding()[:, 0]
        vectors = [np.kron(np.kron(a, a), np.eye(2)[:, c]) for c in range(2)]
        integers = [np.rint(2*v.real).astype(np.int64) for v in vectors]
        for v, integer in zip(vectors, integers):
            self.close(v, integer/2)
        weight = sum(Fraction(int(v@numerator@v), 8*denominator) for v in integers)
        self.assertEqual(weight, Fraction(1, 3))
        h, edges = running_hamiltonian()
        p = numerator/denominator
        u = old.evolve(h, .73)
        final_weight = sum(np.vdot(u@v, p@u@v).real/2 for v in vectors)
        self.assertAlmostEqual(final_weight, 1/3)
        self.assertEqual(1-weight, Fraction(2, 3))
        OBS['conservation_and_error'] = dict(raw_qubits=7,
            exact_projector_denominator=denominator, spin_half_rank=28,
            every_pair_exchange_commutes_by_integer_certificate=True,
            permanently_present_edges=[list(edge) for edge in edges],
            independent_input_spin_half_weight='1/3',
            weight_at_all_times='1/3',
            worst_input_trace_distance_to_every_target_supported_state_at_least='2/3',
            half_diamond_distance_to_every_target_supported_channel_at_least='2/3',
            at_least_one_input_assertion_not_every_input=True,
            all_helper_degrees_in_target_subject=True)

    def test_04_postselection_changes_reference_and_retry_cannot_help(self):
        # Effective order G_A,G_B,helper,R. G_A is Bell-correlated with R;
        # G_B=0 and helper=I/2 are independent of that joint input.
        rho4 = np.zeros((16, 16), dtype=np.int64)
        for helper in range(2):
            indices = [((a*2+0)*2+helper)*2+a for a in range(2)]
            for i in indices:
                for j in indices:
                    rho4[i, j] += 1
        c = sum(old.swap(3, a, b).real.astype(np.int64)
                for a, b in ((0, 1), (0, 2), (1, 2)))
        q = 3*np.eye(8, dtype=np.int64)-c  # 3 P_{j=1/2}
        projected = np.kron(q, np.eye(2, dtype=np.int64))@rho4@np.kron(q, np.eye(2, dtype=np.int64))
        trace = int(np.trace(projected))
        self.assertEqual(trace, 18)
        probability = Fraction(trace, 36)
        self.assertEqual(probability, Fraction(1, 2))
        reference_numerator = np.trace(projected.reshape(8, 2, 8, 2), axis1=0, axis2=2)
        self.assertTrue(np.array_equal(reference_numerator, np.diag([6, 12])))
        self.assertEqual(Fraction(1, 2)-Fraction(1, 3), Fraction(1, 6))
        h = old.swap(3, 0, 2)+2*old.swap(3, 1, 2)+old.swap(3, 0, 1)
        p = q/3
        failure = np.eye(8)-p
        self.close(p@old.evolve(h, .57)@failure, np.zeros((8, 8)))
        OBS['postselection'] = dict(reference_initial_state=['1/2', '1/2'],
            success_probability='1/2',
            reference_success_numerator=[6, 12], reference_success_denominator=18,
            reference_success_state=['1/3', '2/3'],
            reference_trace_distance_change='1/6',
            unknown_state_isometry_on_success_proved=False,
            failed_branch_cannot_reenter_half_sector_under_same_closed_exchange=True,
            measurement_implementation_not_assumed_free=True)

    def test_05_larger_coherent_interface_keeps_unknown_information(self):
        p_half, p_three = sector(7, 1), sector(7, 3)
        q = p_half+p_three
        self.close(q@q, q)
        self.assertAlmostEqual(np.trace(p_half).real, 28)
        self.assertAlmostEqual(np.trace(p_three).real, 56)
        w = raw_old_embedding()
        self.close(w.conj().T@w, np.eye(32))
        self.close(q@w, w)
        h, _ = running_hamiltonian()
        self.close(h@q, q@h)
        values, vectors = np.linalg.eigh(q)
        v = vectors[:, values > .5]
        self.assertEqual(v.shape, (128, 84))
        effective = v.conj().T@h@v
        self.close(h@v, v@effective)
        rng = np.random.default_rng(46105)
        psi = rng.normal(size=(32, 3))+1j*rng.normal(size=(32, 3))
        psi /= np.linalg.norm(psi)
        initial = w@psi
        u = old.evolve(h, .73)
        final = u@initial
        self.close(final, v@old.evolve(effective, .73)@v.conj().T@initial)
        self.close(u.conj().T@final, initial)
        self.close(final.conj().T@final, initial.conj().T@initial)
        self.close(q@final, final)
        # The already established 431 reader is a member of this allowed
        # constant-exchange family, not a new numeric signal experiment.
        inherited_h = sum(old.swap(7, a, a+1) for a in range(4))
        self.close(inherited_h@q, q@inherited_h)
        OBS['alternative_interface'] = dict(
            physical_qubits=7, initial_full_code_columns=32,
            spin_sectors=['1/2', '3/2'], multiplicities=[14, 14],
            coherent_total_dimension=84, arbitrary_unknown_reference_dimension_checked=3,
            old_private_L_and_G_information_not_discarded=True,
            all_fixed_pair_exchange_Hamiltonians_preserve_interface=True,
            old_three_qubit_codes_not_required_preserved=True,
            all_j_coherences_are_kept_as_quantum_state=True,
            actual_private_record_function_reused_from_round=431,
            arbitrary_control_and_full_tomography_not_derived=True,
            special_preparation_into_global_half_sector_not_required=True)

    def test_06_sector_label_dephasing_changes_future_partner_readout(self):
        # Effective old G_A=G_B=0, helper=I/2, exactly an allowed product
        # admission input. The later partner is independently prepared in 0.
        c = sum(old.swap(3, a, b).real.astype(np.int64)
                for a, b in ((0, 1), (0, 2), (1, 2)))
        half3 = 3*np.eye(8, dtype=np.int64)-c
        initial2 = np.diag([1, 1, 0, 0, 0, 0, 0, 0])
        dephased18 = half3@initial2@half3+c@initial2@c
        zero = np.diag([1, 0])
        effect2 = np.eye(16, dtype=np.int64)-old.swap(4, 0, 3).real.astype(np.int64)
        numerator = int(np.trace(effect2@np.kron(dephased18, zero)))
        self.assertEqual(numerator, 2)
        self.assertEqual(Fraction(numerator, 36), Fraction(1, 18))
        self.assertEqual(int(np.trace(effect2@np.kron(initial2, zero))), 0)
        # Full physical realization: A and B are the unchanged round-430
        # three-qubit encodings with private L=0, followed by one helper.
        w = np.kron(np.kron(old.encoding()[:, [0, 2]], old.encoding()[:, [0, 2]]), np.eye(2))
        raw_initial = w@(initial2/2)@w.conj().T
        raw_dephased = (sector(7, 1)@raw_initial@sector(7, 1)
                       +sector(7, 3)@raw_initial@sector(7, 3))
        self.close(raw_dephased, w@(dephased18/18)@w.conj().T)
        # Total singlet on the three actual A spins and the new physical
        # qubit realizes the effective A/new-partner singlet effect.
        a_new = np.kron(old.encoding()[:, [0, 2]], np.eye(2))
        self.close(sector(4, 0)@a_new, a_new@SINGLET)
        OBS['coherence_and_new_partner'] = dict(
            old_input='A=V3|G0,L0>, B=V3|G0,L0>, helper=I/2',
            new_partner='one independently prepared |0> qubit',
            alignment_is_an_explicit_preparation_relation=True,
            unknown_qubit_cloning_not_used=True,
            actual_effect='total singlet of old A three spins and new qubit',
            probability_without_sector_dephasing='0',
            probability_after_sector_dephasing='1/18',
            exact_integer_numerator=2, exact_integer_denominator=36,
            arbitrary_unreachable_84_dimensional_input_used=False,
            not_just_classical_j_weights_in_alternative_interface=True)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=461, baseline_round=459, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(arbitrary_finite_fixed_independent_auxiliary_covered=True,
            exact_minimax_for_global_spin_half_target='1/3',
            all_helpers_retained_inside_specified_target_subject=True,
            unknown_subjects_cannot_compose_claimed=False,
            larger_coherent_continuously_usable_interface_constructed=True,
            initial_contacts_and_preparations_are_inputs=True,
            independent_of_round_460=True,
            complete_recursive_cognitive_architecture_proved=False,
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
