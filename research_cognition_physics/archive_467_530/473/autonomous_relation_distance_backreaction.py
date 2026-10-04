"""Round 473: internal SU(2)-invariant data changes a readable tree distance.

Baseline 472. Same six-tree Hamiltonian, initially maximally mixed graph.
No graph preparation phase, time-dependent force, postselection or axis rotation.
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
import branching_tree_distance_audit as old

TARGET = Path(__file__).with_name('autonomous_relation_distance_backreaction_results.json')
OBS = {}


def short(x):
    return float(f'{float(x):.12g}')


def trnorm(a):
    assert np.linalg.norm(a-a.conj().T) < 1e-9
    return float(np.abs(np.linalg.eigvalsh((a+a.conj().T)/2)).sum())


def comm(a, b):
    return a@b-b@a


@lru_cache(None)
def system():
    trees, _, h, f, d = old.six_vertex_sector()
    h = np.rint(h.real).astype(np.int64)
    f = np.rint(f.real).astype(np.int64)
    s = old.core.swap(6, 0, 5).real.astype(np.int64)
    identity = np.eye(64, dtype=np.int64)
    singlet = (identity-s)/32
    triplet = (identity+s)/96
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    return trees, h, f, d, s, singlet, triplet, eigenvalues, eigenvectors


def graph_trace(a):
    return a.reshape(64, 6, 64, 6).trace(axis1=1, axis2=3)


def data_trace(a):
    return a.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)


@lru_cache(None)
def derivatives():
    _, h, _, d, _, _, _, _, _ = system()
    # Covers every entry, partial trace and final 384-term trace used below.
    assert 384*3*18**12 < 2**63
    k = np.kron(np.eye(64, dtype=np.int64), np.diag(d))
    output = [graph_trace(k)]
    for _ in range(12):
        k = comm(h, k)
        output.append(graph_trace(k))
    return output


@lru_cache(None)
def unitary(t):
    *_, eigenvalues, eigenvectors = system()
    return (eigenvectors*np.exp(-1j*float(t)*eigenvalues))@eigenvectors.conj().T


def single_reduction(rho, site):
    shape = [2]*12
    order = [site]+[i for i in range(6) if i != site]
    return rho.reshape(shape).transpose(order+[i+6 for i in order]).reshape(2, 32, 2, 32).trace(axis1=1, axis2=3)


def local(n, site, pauli):
    out = np.ones((1, 1), dtype=complex)
    for i in range(n):
        out = np.kron(out, pauli if i == site else np.eye(2))
    return out


def contrast_coefficients():
    s = system()[4]
    return {k: Q(-((-1)**(k//2))*int(np.trace(s@derivatives()[k])),
                 144*math.factorial(k)) for k in range(2, 13, 2)}


class Audit(unittest.TestCase):
    def close(self, a, b, tol=2e-10):
        self.assertLess(np.linalg.norm(a-b), tol)

    def test_01_axis_free_sources_and_initial_stationary_graph(self):
        trees, _, f, d, s, singlet, triplet, _, _ = system()
        for source in (singlet, triplet):
            self.assertAlmostEqual(float(np.trace(source)), 1)
            self.assertGreaterEqual(np.linalg.eigvalsh(source).min(), -1e-13)
            for site in range(6):
                self.close(single_reduction(source, site), np.eye(2)/2)
            for pauli in (np.array([[0, 1], [1, 0]]),
                          np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])):
                collective = sum(local(6, site, pauli) for site in range(6))
                self.close(comm(collective, source), np.zeros((64, 64)), 1e-13)
        self.close(singlet/4+3*triplet/4, np.eye(64)/64, 1e-13)
        self.assertAlmostEqual(trnorm(singlet-triplet)/2, 1)
        graph = np.eye(6)/6
        self.close(comm(f, graph), np.zeros((6, 6)), 1e-13)
        for edge in itertools.combinations(range(6), 2):
            n = np.diag([edge in tree for tree in trees])
            self.close(comm(n, graph), np.zeros((6, 6)), 1e-13)
        self.assertEqual(Q(int(sum(d)), 6), Q(8, 3))
        OBS['sources'] = dict(
            singlet='(I-S_05)/32', triplet='(I+S_05)/96',
            graph='I_6/6', initial_mean_distance='8/3',
            every_single_data_marginal='I_2/2', collective_SU2_invariant=True,
            weighted_source_mixture_is_maximally_mixed=True,
            graph_commutes_with_F_and_every_edge_projector=True,
            both_sources_have_same_graph_and_no_external_axis=True,
            singlet_and_triplet_preparation_not_generated_by_this_round=True)

    def test_02_complete_operator_jet_and_coupling_polynomial(self):
        _, h, f, d, s, _, _, _, _ = system()
        hf = np.kron(np.eye(64, dtype=np.int64), f)
        hd = h-hf
        pieces = {0: np.kron(np.eye(64, dtype=np.int64), np.diag(d))}
        low_grade_count = 0
        for order in range(1, 7):
            next_pieces = {}
            for grade, a in pieces.items():
                for added, term in ((0, hd), (1, hf)):
                    value = comm(term, a)
                    next_pieces[grade+added] = next_pieces.get(grade+added, 0)+value
            pieces = next_pieces
            if order <= 5:
                for value in pieces.values():
                    self.assertFalse(np.any(graph_trace(value)))
                    low_grade_count += 1
        coefficients = {grade: graph_trace(a) for grade, a in pieces.items()}
        self.assertEqual([grade for grade, a in coefficients.items() if np.any(a)], [2, 3])
        self.assertTrue(np.array_equal(sum(coefficients.values()), derivatives()[6]))
        trace2, trace3 = (int(np.trace(s@coefficients[k])) for k in (2, 3))
        self.assertEqual((trace2, trace3), (3840, -15360))
        self.assertTrue(np.any(4*coefficients[2]+coefficients[3]))
        self.assertTrue(np.any(derivatives()[7]))
        # Independent block recursion, with the graph indices outside the data blocks.
        blocks = hd.reshape(64, 6, 64, 6).transpose(1, 3, 0, 2)
        a = np.zeros((6, 6, 64, 64), dtype=np.int64)
        for g in range(6):
            a[g, g] = d[g]*np.eye(64, dtype=np.int64)
        for order in range(1, 7):
            b = np.zeros_like(a)
            for g, k in itertools.product(range(6), repeat=2):
                b[g, k] = blocks[g, g]@a[g, k]-a[g, k]@blocks[k, k]
                for q in range(6):
                    b[g, k] += f[g, q]*a[q, k]-a[g, q]*f[q, k]
            a = b
            self.assertTrue(np.array_equal(sum(a[g, g] for g in range(6)), derivatives()[order]))
        OBS['operator_certificate'] = dict(
            low_orders_zero=list(range(1, 6)), exact_low_order_coupling_grades_checked=low_grade_count,
            full_data_dimension=64, full_joint_dimension=384,
            sixth_nonzero_grades=['kappa^2*J^4', 'kappa^3*J^3'],
            coefficient_matrix_entry_maxima=[int(np.abs(coefficients[k]).max()) for k in (2, 3)],
            trace_S05_coefficients=[trace2, trace3],
            singlet_sixth_derivative='20*kappa^2*J^3*(J-4*kappa)',
            triplet_sixth_derivative='-20*kappa^2*J^3*(J-4*kappa)/3',
            sixth_J1_kappa1=[-60, 20], seventh_operator_nonzero=True,
            J4kappa_cancels_these_source_coefficients_not_entire_A6=True,
            independent_graph_block_recursion_equal=True,
            integer_arithmetic_bound=str(384*3*18**12), signed_int64_limit=str(2**63))

    def test_03_rational_finite_window(self):
        coefficients = contrast_coefficients()
        expected = {2: Q(0), 4: Q(0), 6: Q(-1, 9), 8: Q(37, 180),
                    10: Q(-561, 2800), 12: Q(551539, 4082400)}
        self.assertEqual(coefficients, expected)
        lo, hi = Q(1, 32), Q(1, 16)
        x = 18*hi
        tail_per_t6 = Q(18)**14*hi**8/Q(math.factorial(14))/(1-x*x/Q(240))
        upper = Q(-1, 9)+Q(37, 180)*hi**2+Q(551539, 4082400)*hi**6+tail_per_t6
        lower = Q(-1, 9)-Q(561, 2800)*hi**4-tail_per_t6
        self.assertLessEqual(upper, Q(-1, 10))
        self.assertGreaterEqual(lower, Q(-1, 8))
        gap = lo**6/10
        self.assertEqual(gap, Q(1, 10737418240))
        OBS['strict_window'] = dict(
            contrast='mean_D_singlet - mean_D_triplet',
            Taylor_coefficients={str(k): str(c) for k, c in coefficients.items()},
            scalar_source_contrast_even_in_time=True,
            full_operator_or_arbitrary_reference_not_claimed_even=True,
            remainder='(18*t)^14/[14!*(1-(18*t)^2/240)]',
            closed_window=[str(lo), str(hi)],
            whole_window_contrast_bounds=['-t^6/8', '-t^6/10'],
            uniform_positive_gap=str(gap),
            scaled_upper_bound=str(upper), scaled_lower_bound=str(lower),
            no_time_peak_search=True)

    def test_04_full_evolution_changes_graph_population_and_distance(self):
        _, _, _, d, _, singlet, triplet, _, _ = system()
        t = Q(1, 16)
        u = unitary(t)
        rows, outputs = [], []
        for label, source in (('singlet', singlet), ('triplet', triplet)):
            initial = np.kron(source, np.eye(6)/6)
            final = u@initial@u.conj().T
            graph = data_trace(final)
            distance = float(np.trace(np.diag(d)@graph).real)
            self.assertGreater(np.linalg.norm(graph.diagonal()-np.ones(6)/6), 1e-11)
            rows.append(dict(source=label, graph_population=[short(x.real) for x in graph.diagonal()],
                             mean_distance=short(distance),
                             graph_change_trace_distance=short(trnorm(graph-np.eye(6)/6)/2)))
            outputs.append(distance)
        contrast = outputs[0]-outputs[1]
        polynomial = sum(float(c*t**k) for k, c in contrast_coefficients().items())
        x = 18*t
        tail = x**14/Q(math.factorial(14))/(1-x*x/Q(240))
        self.assertLess(abs(contrast-polynomial), float(tail)+3e-14)
        self.assertLess(contrast, -float(t**6/10))
        self.assertGreater(contrast, -float(t**6/8))
        self.assertLess(outputs[0], 8/3)
        self.assertGreater(outputs[1], 8/3)
        self.assertAlmostEqual(outputs[0]/4+3*outputs[1]/4, 8/3, places=13)
        OBS['actual_backreaction'] = dict(
            time=str(t), J=1, kappa=1, output=rows,
            actual_distance_contrast=short(contrast), rational_Taylor_value=short(polynomial),
            rational_remainder_bound=str(tail),
            initially_maximally_mixed_graph_is_not_invariant_for_all_data=True,
            graph_population_changes_not_only_graph_phase=True,
            weighted_maximally_mixed_data_mixture_remains_stationary=True,
            no_graph_phase_preparation_or_external_axis_rotation=True)

    def test_05_arbitrary_reference_effect_and_uniform_remainder(self):
        _, _, _, d, _, _, _, _, _ = system()
        t = Q(1, 16)
        u = unitary(t)
        observable = np.kron(np.eye(64), np.diag(d))
        effect = graph_trace(u.conj().T@observable@u)/6
        approximate = np.eye(64)*(8/3)-float(t**6)*derivatives()[6]/4320
        x = 18*t
        bound = x**7/Q(2*math.factorial(7))/(1-x/8)
        operator_error = float(np.linalg.norm(effect-approximate, 2))
        self.assertLess(operator_error, float(bound))
        self.assertGreaterEqual(np.linalg.eigvalsh(effect).min(), 2-1e-12)
        self.assertLessEqual(np.linalg.eigvalsh(effect).max(), 3+1e-12)
        raw = np.array([[complex(((3*i+2*r)%7)-3, ((i+3*r)%5)-2)
                         for r in range(3)] for i in range(64)])
        psi = raw/np.linalg.norm(raw)
        input_columns = np.einsum('dr,gh->dgrh', psi, np.eye(6)).reshape(384, 18)
        out = (u@input_columns).reshape(64, 6, 3, 6)
        actual = np.einsum('dgrh,g,dgsh->rs', out, d, out.conj())/6
        by_effect = psi.T@effect.T@psi.conj()
        approximate_reference = psi.T@approximate.T@psi.conj()
        self.close(actual, by_effect)
        error = trnorm(actual-approximate_reference)
        self.assertLess(error, float(bound))
        OBS['unknown_reference'] = dict(
            time=str(t), reference_dimension=3,
            effective_data_observable='M_D(t)=Tr_G[U^dagger*(I_D tensor D)*U]/6',
            sixth_approximation='(8/3)*I_D - t^6*A6/4320',
            general_remainder_starts_at_order=7,
            uniform_operator_remainder_bound=str(bound),
            operator_error=short(operator_error),
            reference_distance_moment_error=short(error),
            direct_full_U_and_effect_reference_moments_equal=True,
            reference_initially_may_be_entangled_with_unknown_data=True,
            graph_initially_independent_and_maximally_mixed=True)

    def test_06_actual_readout_handoff_and_finite_resource_certificate(self):
        gap = Q(1, 10737418240)
        read_error = gap/4
        paths, weight, m, h = 4, 10, 3, Q(25)
        per_path = read_error/weight
        t = per_path/(48*m*h*h)
        gamma = per_path*t/(8*m)
        x = 2*h*t
        remainder = x*x/(2*(1-x/3))
        delta = (remainder+gamma)/t
        bias = (1+delta)**m-1
        self.assertLess(weight*bias, read_error/2)
        self.assertGreater(sum(Q(8)**k/Q(math.factorial(k)) for k in range(25)), 1600)
        copies = math.ceil(64/(per_path**2*t**(2*m)))
        self.assertGreaterEqual(Q(copies), 64/(per_path**2*t**(2*m)))
        # Two distance estimates, each error <= gap/4, still have a fixed sign.
        self.assertEqual(gap-2*read_error, gap/2)
        # Exact register permutation: active data become internal storage while
        # new probes become active. This is an additional local-SWAP permission.
        original = np.arange(2*3*2).reshape(2, 3, 2)
        fresh = np.array([1, 1j])/math.sqrt(2)
        combined = np.einsum('dgr,p->dgpr', original, fresh)
        handed = combined.transpose(2, 1, 0, 3)
        target = np.einsum('p,sgr->pgsr', fresh, original)
        self.assertTrue(np.array_equal(handed, target))
        OBS['readout_connection'] = dict(
            inherited_protocol_round=472,
            natural_wait_window=['1/32', '1/16'], uniform_distance_gap=str(gap),
            per_source_distance_estimation_error=str(read_error),
            probability_of_either_estimation_failure_at_most='1/100',
            known_department_path_count=paths, sum_path_lengths=weight,
            probe_wait=str(t), complete_instrument_diamond_budget=str(gamma),
            copies_per_path_per_source=str(copies), total_source_trials=str(2*paths*copies),
            retained_distinguishing_gap=str(gap/2),
            exact_local_storage_relabeling_verified=True,
            old_data_and_all_graph_correlations_move_to_internal_storage=True,
            local_swap_preparation_readout_clock_storage_are_explicit_inputs=True,
            actual_instruments_fixed_without_history_feedback=True,
            storage_duration_and_error_in_complete_instrument_budget=True,
            no_cloning_of_one_unknown_graph_no_free_reset=True,
            no_actual_huge_sampling_experiment=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=473, baseline_round=472, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(original_H_no_intermediate_control_during_natural_wait=True,
            initial_graph_maximally_mixed_no_graph_phase_engineered=True,
            internal_SU2_invariant_data_relation_changes_actual_distance_population=True,
            finite_six_tree_operator_certificate_not_all_tree_low_order_theorem=True,
            arbitrary_reference_effect_certificate=True,
            source_preparations_and_followup_readout_remain_inputs=True,
            independent_source_copies_and_internal_storage_counted=True,
            no_universal_attraction_claim=True,
            strong_marginal_stability_not_imposed_as_cognitive_axiom=True,
            four_leaf_department_has_only_two_independent_shape_probabilities=True,
            displacement_or_three_dimensional_space_derived=False,
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
