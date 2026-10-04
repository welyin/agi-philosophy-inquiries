"""Round 511: extracting an anchored reader with unknown quantum membership.

GHZ deletion is an established primitive. This audit retains coherent membership,
all classical outcomes, the continuing original drift, and the source contract.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import json
from pathlib import Path
import platform
import unittest

import numpy as np
import historical_membership_interface as history
import heralded_dynamic_cat_source as source
from distributed_role_reader import exponential, operators, opnorm

HERE = Path(__file__).resolve().parent
TARGET = HERE/'coherent_member_relay_results.json'
OBS = {}
HAD = np.array([[1., 1.], [1., -1.]])/np.sqrt(2)
Z = np.diag([1., -1.])


def local_k(member, result):
    if member:
        return np.linalg.matrix_power(Z, result)/np.sqrt(2)
    return np.outer([1., 0.], HAD[result])


def local_unitary():
    # Tensor order M,A,R, with R initialized in zero.
    swap = np.eye(4)[[0, 2, 1, 3]]
    u0 = np.kron(np.eye(2), HAD) @ swap
    u1 = np.diag([1., 1., 1., -1.]) @ np.kron(np.eye(2), HAD)
    result = np.zeros((8, 8))
    result[:4, :4], result[4:, 4:] = u0, u1
    return result


def apply_local(vec, gate, vertex):
    """A first; arbitrary spectator tensor follows. Little-endian A bits."""
    zero = np.flatnonzero((np.arange(len(vec)) & (1 << vertex)) == 0)
    one = zero | (1 << vertex)
    out = np.empty_like(vec, dtype=complex)
    out[zero] = gate[0, 0]*vec[zero]+gate[0, 1]*vec[one]
    out[one] = gate[1, 0]*vec[zero]+gate[1, 1]*vec[one]
    return out


def cat(n):
    out = np.zeros(2**n, dtype=complex)
    out[0] = out[-1] = 1/np.sqrt(2)
    return out


def selected_cat(n, members, anchor):
    out = np.zeros(2**n, dtype=complex)
    out[0] = out[members | (1 << anchor)] = 1/np.sqrt(2)
    return out


def branch(vec, n, members, anchor, record, correct=True):
    others = [v for v in range(n) if v != anchor]
    out = vec.copy()
    for j, v in enumerate(others):
        out = apply_local(out, local_k((members >> v) & 1, (record >> j) & 1), v)
    if correct and record.bit_count() % 2:
        out = apply_local(out, Z, anchor)
    return out


def pure_difference(a, b):
    """Raw trace norm of |a><a|-|b><b|, for unnormalized vectors."""
    aa, bb = np.vdot(a, a).real, np.vdot(b, b).real
    return float(np.sqrt(max(0., (aa+bb)**2-4*abs(np.vdot(a, b))**2)))


def short(value):
    return float(f'{float(value):.12g}')


class Audit(unittest.TestCase):
    def test_01_complete_local_unitary_and_instrument(self):
        u = local_unitary()
        self.assertLess(opnorm(u.T @ u-np.eye(8)), 1e-14)
        worst = 0.
        complete = np.zeros((4, 4), complex)
        for r in range(2):
            actual = u[r::2, ::2]
            expected = np.zeros((4, 4))
            expected[:2, :2], expected[2:, 2:] = local_k(0, r), local_k(1, r)
            worst = max(worst, opnorm(actual-expected))
            complete += actual.conj().T @ actual
        self.assertLess(worst, 1e-14)
        self.assertLess(opnorm(complete-np.eye(4)), 1e-14)
        OBS['local_instrument'] = dict(unitary_dimension=8, instrument_input_dimension=4,
            full_unitary_to_kraus_error=short(worst), completeness_error=short(opnorm(complete-np.eye(4))),
            local_blank_record_qubits=1, membership_is_not_measured=True)

    def test_02_all_members_all_records_and_reference_isometry(self):
        rows = []
        for n in (2, 3, 4, 6):
            worst = 0.
            norm_error = 0.
            exit_error = 0.
            for anchor in (0, n-1):
                for m in range(2**n):
                    target = selected_cat(n, m, anchor)
                    for z in range(2**(n-1)):
                        actual = branch(cat(n), n, m, anchor, z)
                        worst = max(worst, float(np.linalg.norm(actual-target/np.sqrt(2**(n-1)))))
                        norm_error = max(norm_error, abs(np.vdot(actual, actual).real-2.**(1-n)))
                    for v in range(n):
                        if v != anchor and not ((m >> v) & 1):
                            exit_error = max(exit_error, float(np.linalg.norm(target[(np.arange(2**n) >> v) & 1 == 1])))
            # Distinct M columns have orthogonal M labels even when cat overlaps.
            w = np.zeros((2**(2*n), 2**n), complex)
            for m in range(2**n):
                w[m*2**n:(m+1)*2**n, m] = selected_cat(n, m, 0)
            gram_error = opnorm(w.conj().T @ w-np.eye(2**n))
            self.assertLess(max(worst, norm_error, exit_error, gram_error), 1e-13)
            rows.append(dict(N=n, anchors_tested=[0,n-1], membership_patterns=2**n,
                records_per_pattern=2**(n-1), full_branch_column_error=short(worst),
                record_effect_error=short(norm_error), conditional_relay_exit_error=short(exit_error),
                complete_W_isometry_error=short(gram_error)))
        OBS['complete_coherent_contract'] = rows

    def test_03_record_leakage_marginal_disturbance_and_pending_parity(self):
        n, anchor = 3, 0
        c0, c1 = selected_cat(n, 0, anchor), selected_cat(n, 2, anchor)
        overlap = float(np.vdot(c0, c1).real)
        reduced = np.array([[1., overlap], [overlap, 1.]])/2
        initial = np.ones((2, 2))/2
        marginal_difference = float(np.abs(np.linalg.eigvalsh(reduced-initial)).sum())
        self.assertAlmostEqual(overlap, .5)
        self.assertAlmostEqual(marginal_difference, .5)
        # The full map has an exact left inverse on its image, including a reference.
        w = np.zeros((2*2**n, 2), complex)
        w[:2**n, 0], w[2**n:, 1] = c0, c1
        bell = np.eye(2)/np.sqrt(2)
        joint = w @ bell
        recovery = np.linalg.norm(w.conj().T @ joint-bell)
        self.assertLess(recovery, 1e-14)
        self.assertLess(np.linalg.norm(joint.conj().T @ joint-np.eye(2)/2), 1e-14)
        # Without randomized retained-member results, a transcript leaks membership.
        naive_all_erased = np.ones(2**(n-1))/2**(n-1)
        naive_all_retained = np.eye(1, 2**(n-1), 0).ravel()
        tv = .5*np.abs(naive_all_erased-naive_all_retained).sum()
        self.assertAlmostEqual(tv, .75)
        # Ignoring the parity before delivery dephases the anchored cat.
        m = 6
        rho = sum(np.outer(v, v.conj()) for v in (
            branch(cat(n), n, m, anchor, z, correct=False) for z in range(4)))
        target = selected_cat(n, m, anchor)
        lost = float(np.abs(np.linalg.eigvalsh(rho-np.outer(target, target.conj()))).sum())
        self.assertAlmostEqual(lost, 1.)
        OBS['information_scope'] = dict(different_membership_reader_overlap=overlap,
            original_M_raw_trace_norm_change=marginal_difference,
            complete_MR_left_inverse_error=short(recovery),
            naive_transcript_total_variation=tv, missing_parity_raw_trace_norm_error=lost,
            old_M_marginal_preservation_claimed=False, physical_local_inverse_claimed=False)

    def test_04_actual_history_and_continuing_drift(self):
        a = history.model()
        n, d = a['n'], a['d']
        old = np.zeros((d, 2), complex)
        old[2*a['G'], 0] = old[3*a['G']+1, 1] = 1/np.sqrt(2)
        generated = np.einsum('mij,jr->mir', history.history(Q(1,2)), old)
        u = exponential(float(Q(1,7))*a['h'])
        advanced = np.einsum('ij,mjr->mir', u, generated)
        target = np.array([selected_cat(n,m,0)[:,None,None]*advanced[m][None,:,:]
                           for m in range(2**n)])
        first = np.array([selected_cat(n,m,0)[:,None,None]*generated[m][None,:,:]
                          for m in range(2**n)])
        evolve_after = np.einsum('ij,majr->mair', u, first)
        commuting = float(np.linalg.norm(target-evolve_after))
        z = 19
        actual = np.array([branch(cat(n)[:,None,None]*advanced[m][None,:,:], n,m,0,z)
                           for m in range(2**n)])
        error = float(np.linalg.norm(actual-target/np.sqrt(2**(n-1))))
        self.assertLess(max(error, commuting), 1e-13)
        # Tracing M,A leaves all old D/G and reference marginals at their free drift.
        old_dr = sum(x.reshape(-1,1) @ x.reshape(1,-1).conj() for x in advanced)
        after_dr = sum(x.reshape(-1,1) @ x.reshape(1,-1).conj()
                       for m in target for x in m)
        spectator = float(np.linalg.norm(old_dr-after_dr))
        self.assertLess(spectator, 1e-13)
        OBS['actual_history_bridge'] = dict(history_time='1/2', free_wait='1/7',
            arbitrary_reference_example_dimension=2, original_DG_dimension=d,
            history_branch_to_W_error=short(error), continuous_h_commutation_error=short(commuting),
            DGR_marginal_preservation_error=short(spectator),
            M_included_in_preserved_marginal=False)

    def test_05_leaf_members_need_outside_assistance(self):
        a = history.model()
        leaf, other = 2, 3
        graph = next(j for j,t in enumerate(a['trees']) if (0,leaf) in t and (0,other) in t)
        self.assertTrue(all((leaf,other) not in t for t in a['trees']))
        start = np.zeros((64,36,1),complex)
        start[0,leaf*a['G']+graph,0] = 1
        state = start
        coefficients = []
        m = (1 << leaf) | (1 << other)
        for order in range(1,5):
            state = history.apply_historical(state)
            coefficients.append(int(round(state[m,other*a['G']+graph,0].real)))
        self.assertEqual(coefficients, [0,0,0,1])
        t = Q(1,131072)
        remainder = Q(11**5,120)*t**5
        self.assertLessEqual(remainder, t**4/48)
        amplitude_floor = t**4/48
        # With anchor one of the two selected leaves, the extracted reader is Bell.
        target = selected_cat(6,m,leaf)
        basis = [0, (1 << leaf), (1 << other), m]
        bell = target[basis]
        self.assertTrue(np.allclose(bell, [1/np.sqrt(2),0,0,1/np.sqrt(2)]))
        OBS['leaf_counterexample_and_delivery'] = dict(initial_graph_index=graph,
            members=[leaf,other], anchor=leaf, H_powers_to_selected_output=coefficients,
            fourth_amplitude_coefficient='1/24', certified_time=str(t),
            amplitude_floor=str(amplitude_floor), probability_floor=str(amplitude_floor**2),
            all_member_internal_detectors_zero=True, extracted_pair_is_Bell=True,
            outside_shared_resource_and_parity_delivery_required=True)

    def test_06_heralded_source_composition_with_old_correlated_members(self):
        a = source.finite_model()
        n, d, anchor, k = a['n'], a['d'], 0, 12
        members = [12, 48]
        old = np.zeros((2,d,2),complex)
        old[0,12,0] = old[1,19,1] = 1/np.sqrt(2)
        source_out = np.empty((2,2**n,d,2),complex)
        for j in range(2):
            for z in range(2**n):
                source_out[j,z] = np.linalg.matrix_power(a['layers'][z],k) @ old[j]/np.sqrt(2**n)
        success = float(np.vdot(source_out,source_out).real)
        free = exponential(float(k*a['T'])*a['h'])
        target = np.array([selected_cat(n,m,anchor)[:,None,None]*(free @ old[j])[None,:,:]
                           for j,m in enumerate(members)])
        source_ideal = np.array([cat(n)[:,None,None]*(free @ old[j])[None,:,:] for j in range(2)])
        before = pure_difference(source_out/np.sqrt(success),source_ideal)
        full_record_difference = 0.
        probabilities = []
        total = 0.
        for z in range(2**(n-1)):
            actual = np.array([branch(source_out[j],n,m,anchor,z) for j,m in enumerate(members)])
            probability = float(np.vdot(actual,actual).real)
            total += probability
            probabilities.append(probability/success)
            full_record_difference += pure_difference(actual/np.sqrt(success),target/np.sqrt(2**(n-1)))
        p0 = Q(1,2**(n-1))
        r2 = (1-p0)*Q(3,4)**(2*k)
        bound = 2*np.sqrt(float(r2/p0))
        self.assertAlmostEqual(total,success,places=13)
        self.assertGreaterEqual(success,float(p0)-1e-12)
        self.assertLessEqual(success,float(p0+r2)+1e-12)
        self.assertLessEqual(full_record_difference,before+1e-10)
        self.assertLessEqual(full_record_difference,bound+1e-10)
        OBS['finite_source_composition'] = dict(layers=k, success_probability=short(success),
            ideal_success_floor=str(p0), conditional_source_raw_trace_norm_error=short(before),
            output_with_all_records_raw_trace_norm_error=short(full_record_difference),
            analytic_conditional_bound=short(bound),
            actual_record_probability_min=short(min(probabilities)),
            actual_record_probability_max=short(max(probabilities)),
            total_success_retained_error=short(abs(total-success)),
            conditional_normalization_is_not_a_diamond_channel=True)

    def test_07_selected_reader_is_not_the_original_global_parity_interface(self):
        a = history.model()
        n, d, m, anchor = a['n'], a['d'], 12, 2
        support = m | (1 << anchor)
        q = np.diag(a['q'][m])
        u0, u1, good, ideal = operators(a['h'], q, 1/64)
        odd_effect = np.zeros((d,d),complex)
        complete = np.zeros_like(odd_effect)
        correct_effect = np.zeros_like(odd_effect)
        for z in range(2**n):
            sign = (-1)**((z & support).bit_count())
            kz = (u0+sign*u1)/np.sqrt(2**(n+1))
            effect = kz.conj().T @ kz
            complete += effect
            if z.bit_count()%2:
                odd_effect += effect
            if (z & support).bit_count()%2:
                correct_effect += effect
        self.assertLess(opnorm(complete-np.eye(d)),1e-12)
        self.assertLess(opnorm(odd_effect-np.eye(d)/2),1e-12)
        self.assertLess(opnorm(correct_effect-good[1].conj().T @ good[1]),1e-12)
        psi = np.zeros(d,complex)
        psi[anchor*a['G']:(anchor+1)*a['G']] = 1/np.sqrt(a['G'])
        wanted = float(np.linalg.norm(ideal[1] @ psi)**2)
        selected = float(np.vdot(psi,correct_effect @ psi).real)
        self.assertGreater(wanted,.99)
        OBS['readout_non_substitution'] = dict(members=[2,3],anchor=anchor,read_time='1/64',
            ordinary_full_parity_effect_error_from_half_I=short(opnorm(odd_effect-np.eye(d)/2)),
            naive_full_parity_probability=.5, intended_Lueders_probability=short(wanted),
            known_selected_parity_probability=short(selected),
            unknown_membership_parity_processor_implemented=False)

    def test_08_recorded_membership_readout_and_all_input_error_bound(self):
        a = history.model()
        n, d, anchor, tau = a['n'], a['d'], 0, 1/64
        worst_fine = worst_complete = worst_dilation = 0.
        largest_commutator = 0.
        for m in range(2**n):
            support = m | (1 << anchor)
            q = np.diag(a['q'][m])
            u0,u1,ks,ideals = operators(a['h'],q,tau)
            complete = np.zeros((d,d),complex)
            counts = [0,0]
            for x in range(2**n):
                r = (x & support).bit_count()%2
                counts[r] += 1
                actual = (u0+(-1)**r*u1)/np.sqrt(2**(n+1))
                worst_fine = max(worst_fine,float(np.max(np.abs(actual-ks[r]/np.sqrt(2**(n-1))))))
                complete += actual.conj().T @ actual
            self.assertEqual(counts,[2**(n-1)]*2)
            worst_complete = max(worst_complete,opnorm(complete-np.eye(d)))
            # Common outcome dilation gives an ordinary diamond upper bound.
            dilation = 2*opnorm(np.vstack([ks[r]-ideals[r] for r in range(2)]))
            worst_dilation = max(worst_dilation,dilation)
            largest_commutator = max(largest_commutator,opnorm(a['h'] @ q-q @ a['h']))
        bound = 3*np.pi*tau/np.sqrt(2)
        self.assertLess(max(worst_fine,worst_complete),1e-12)
        self.assertLessEqual(worst_dilation,bound+1e-12)
        self.assertLessEqual(largest_commutator,3+1e-12)
        OBS['recorded_membership_readout'] = dict(all_membership_patterns=2**n,
            complete_records_per_pattern=2**n, uniform_records_per_selected_parity=2**(n-1),
            fine_instrument_factor_error=short(worst_fine),completeness_error=short(worst_complete),
            direct_sum_dilation_diamond_upper_bound=short(worst_dilation),
            size_uniform_ordinary_diamond_bound=short(bound),
            member_record_is_retained=True,comparison_target_also_records_membership=True,
            arbitrary_initial_input_and_reference_covered_by_proof=True,
            old_statistics_equivalence_requires_M_block_diagonal_future_menu=True,
            autonomous_record_transport_implemented=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output,verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(round=511, scientific_baseline_round=510, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors), python=platform.python_version(),
        numpy=np.__version__, dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ['historical_membership_interface.py','heralded_dynamic_cat_source.py',
                         'research_note_508.md','research_note_510.md']},
        scope=dict(coherent_unknown_membership_supported=True,
            full_transcript_independent_of_membership_for_ideal_cat=True,
            complete_old_information_preserved_in_joint_isometry=True,
            original_h_continues=True, heralded_source_error_composed=True,
            recorded_membership_readout_closed_with_declared_controls=True,
            old_membership_reduced_state_unchanged=False,
            physical_removal_by_unknown_member_list=False,
            autonomous_classical_record_routing_derived=False,
            shared_resource_and_anchor_eliminated=False,
            independent_fixed_cat_readout_formula_inherited=False,
            pure_original_swap_only_implementation=False,
            multiple_complete_spatial_endpoints_generated=False,
            dimension_three_generated=False, full_GR_goal_completed=False,
            phase_closure_triggered=False), observations=OBS)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
