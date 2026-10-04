"""Round 510: finite historical labels and their complete collective instrument.

Local quantum-walk memory is an existing method. New claims here are scoped to
the unchanged moving graph drift, formation budgets, and the complete readout.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import coherent_graph_mean_obstruction as old
from distributed_role_reader import exponential, operators, opnorm

HERE = Path(__file__).resolve().parent
TARGET = HERE/'historical_membership_interface_results.json'
OBS = {}


@lru_cache(maxsize=1)
def model():
    family, _, _, _, h, _, _, _ = old.model(2)
    n, graph_count, d = 6, len(family), len(h)
    modes = 2**n
    signs = np.array([[1-2*((s >> v) & 1) for v in range(n)] for s in range(modes)])
    hadamard = np.array([[(-1)**((s & z).bit_count()) for s in range(modes)] for z in range(modes)])
    return dict(h=h, trees=family, n=n, G=graph_count, d=d, modes=modes,
                signs=signs, hadamard=hadamard,
                q=np.array([np.repeat([(z >> v) & 1 for v in range(n)], graph_count)
                            for z in range(modes)], dtype=np.int64),
                counts=np.array([z.bit_count() for z in range(modes)]))


@lru_cache(maxsize=8)
def history(time, coupling=Q(1)):
    a = model()
    unitaries = np.array([exponential(float(time)*(a['h']+float(coupling)*
                           np.diag(np.repeat(sign, a['G'])))) for sign in a['signs']])
    # Initial Z-memory vacuum and final Z-memory basis each supply 1/sqrt(2^N).
    return np.einsum('zs,sab->zab', a['hadamard'], unitaries)/a['modes']


def apply_historical(state, coupling=1):
    a = model()
    out = np.einsum('ab,mbc->mac', a['h'], state)
    for v in range(a['n']):
        rows = slice(v*a['G'], (v+1)*a['G'])
        flipped = np.arange(a['modes']) ^ (1 << v)
        out[:, rows, :] += coupling*state[flipped, rows, :]
    return out


def effect(blocks, mask=None, weights=None):
    if mask is None:
        mask = np.ones(len(blocks), dtype=bool)
    if weights is None:
        weights = np.ones(len(blocks))
    return sum(weights[j]*blocks[j].conj().T @ blocks[j] for j in np.flatnonzero(mask))


def tail_bound(x, cutoff):
    assert x < cutoff+1
    return x**cutoff/Q(math.factorial(cutoff))/(1-x/Q(cutoff+1))


def rank_two_trace(a, b):
    columns = np.column_stack((a.ravel(), b.ravel()))/np.sqrt(a.shape[-1])
    _, triangular = np.linalg.qr(columns)
    reduced = triangular @ np.diag([1, -1]) @ triangular.conj().T
    return float(np.abs(np.linalg.eigvalsh(reduced)).sum())


class Audit(unittest.TestCase):
    def test_01_complete_history_isometry_and_independent_generator(self):
        a = model()
        duration = Q(1, 8)
        exact = history(duration)
        initial = np.zeros_like(exact)
        initial[0] = np.eye(a['d'])
        term, taylor = initial.copy(), initial.copy()
        for order in range(1, 25):
            term = (-1j*float(duration)/order)*apply_historical(term)
            taylor += term
        difference = opnorm((taylor-exact).reshape(a['modes']*a['d'], a['d']))
        completeness = opnorm(effect(exact)-np.eye(a['d']))
        self.assertLess(difference, 2e-13)
        self.assertLess(completeness, 2e-13)
        self.assertTrue(np.array_equal(a['hadamard'] @ a['hadamard'].T,
                                       a['modes']*np.eye(a['modes'], dtype=np.int64)))
        OBS['history_implementation'] = dict(DG_dimension=a['d'], memory_qubits=a['n'],
            total_dimension=a['modes']*a['d'], duration=str(duration),
            independent_tensor_action_error=old.short(difference),
            unknown_input_isometry_error=old.short(completeness),
            original_h_not_switched_off=True)

    def test_02_record_count_and_operator_tail_bounds(self):
        a = model()
        rows = []
        for duration in (Q(1, 8), Q(1, 2), Q(2)):
            blocks = history(duration)
            number_effect = effect(blocks, weights=a['counts'])
            largest_mean = float(np.linalg.eigvalsh(number_effect)[-1])
            bound = min(duration, duration**2)
            self.assertLessEqual(largest_mean, float(bound)+1e-12)
            tails = []
            for cutoff in (2, 3, 4):
                tail_effect = effect(blocks, mask=a['counts'] >= cutoff)
                largest_tail = float(np.linalg.eigvalsh(tail_effect)[-1])
                series = tail_bound(duration, cutoff)
                probability_bound = min(Q(1), bound/cutoff, series**2)
                self.assertLessEqual(largest_tail, float(probability_bound)+1e-12)
                tails.append(dict(at_least_marks=cutoff, actual_operator_probability=old.short(largest_tail),
                                  rational_probability_upper=str(probability_bound)))
            rows.append(dict(duration=str(duration), maximum_unknown_input_mean_marks=old.short(largest_mean),
                mean_mark_upper=str(bound), tails=tails))
        OBS['formation_budget'] = dict(cases=rows,
            mean_bound='min(abs(g) T, g^2 T^2)',
            tail_bound='min(1, mean_upper/K, [sum_{j>=K} (abs(g)T)^j/j!]^2)',
            high_probability_K_marks_requires='abs(g) T >= (1-eta) K',
            single_excitation_only=True, size_uniform_operator_proof=True)

    def test_03_nonroot_labels_have_a_full_graph_reference_certificate(self):
        a = model()
        duration = Q(1, 1024)
        rows = []
        for root in (0, 2):
            initial = np.zeros((a['modes'], a['d'], a['G']), dtype=np.int64)
            initial[0, root*a['G']:(root+1)*a['G']] = np.eye(a['G'], dtype=np.int64)
            first = apply_historical(initial)
            second = apply_historical(first)
            nonroot = np.array([bool(z & ~(1 << root)) for z in range(a['modes'])])
            self.assertFalse(np.any(first[nonroot]))
            expected = np.zeros_like(initial)
            for v in range(a['n']):
                if v != root:
                    incident = np.diag([tuple(sorted((root, v))) in tree for tree in a['trees']]).astype(int)
                    expected[1 << v, v*a['G']:(v+1)*a['G']] = incident
            projected = np.where(nonroot[:, None, None], second, 0)
            self.assertTrue(np.array_equal(projected, expected))
            degree = 3 if root < 2 else 1
            gram = expected.reshape(-1, a['G']).T @ expected.reshape(-1, a['G'])
            self.assertTrue(np.array_equal(gram, degree*np.eye(a['G'], dtype=np.int64)))
            remainder = Q(11**3, 6)*duration**3
            self.assertLessEqual(remainder, duration**2/4)
            lower = duration**4/16
            columns = history(duration)[:, :, root*a['G']:(root+1)*a['G']]
            spectrum = np.linalg.eigvalsh(effect(columns, mask=nonroot))
            self.assertGreater(float(spectrum[0]), float(lower))
            rows.append(dict(root=root, root_degree=degree, second_order_gram_exact=degree,
                duration=str(duration), universal_nonroot_probability_lower=str(lower),
                observed_min_over_all_graph_inputs=old.short(spectrum[0]),
                observed_max_over_all_graph_inputs=old.short(spectrum[-1])))
        OBS['actual_nonroot_source'] = dict(cases=rows, arbitrary_unknown_graph_reference=True,
            time_window_not_claimed_uniform_in_N=True)

    def test_04_history_selected_collective_reader_on_all_memory_inputs(self):
        a = model()
        tau = Q(1, 64)
        maximum_commutator = maximum_dilation_error = completeness = 0.0
        actual, target = [], []
        for bits in range(a['modes']):
            q = np.diag(a['q'][bits])
            self.assertTrue(np.array_equal(q @ q, q))
            maximum_commutator = max(maximum_commutator, opnorm(a['h'] @ q-q @ a['h']))
            u0, u1, ks, ideal = operators(a['h'], q, tau)
            maximum_dilation_error = max(maximum_dilation_error,
                                         opnorm(u1-(np.eye(a['d'])-2*q) @ u0)/np.sqrt(2))
            completeness = max(completeness, opnorm(sum(k.conj().T @ k for k in ks)-np.eye(a['d'])))
            actual.append(ks)
            target.append(ideal)
        self.assertLessEqual(maximum_commutator, 3+1e-13)
        analytic_bound = 3*np.pi*float(tau)/np.sqrt(2)
        self.assertLessEqual(2*maximum_dilation_error, analytic_bound)
        self.assertLess(completeness, 1e-13)
        # Compose the actual product-source history isometry with both instruments.
        source = history(Q(1, 2))
        combined_difference = 0.0
        for result in (0, 1):
            aa = np.array([actual[z][result] @ source[z] for z in range(a['modes'])])
            bb = np.array([target[z][result] @ source[z] for z in range(a['modes'])])
            combined_difference += rank_two_trace(aa.reshape(-1, a['d']), bb.reshape(-1, a['d']))
        self.assertLess(combined_difference, analytic_bound)
        OBS['complete_member_reader'] = dict(read_duration=str(tau),
            all_quantum_memory_sectors_checked=a['modes'],
            maximum_commutator=old.short(maximum_commutator),
            completeness_error=old.short(completeness),
            direct_dilation_diamond_upper=old.short(2*maximum_dilation_error),
            size_uniform_diamond_upper=old.short(analytic_bound),
            history_then_read_complete_choi_trace_error=old.short(combined_difference),
            arbitrary_entangled_memory_input_covered_by_block_operator_proof=True,
            extra_local_occupation_memory_cat_coupling_declared=True)

    def test_05_classical_correlations_suffice_for_this_read_statistics(self):
        a = model()
        root, duration, tau = 2, Q(1, 2), Q(1, 64)
        initial = np.zeros(a['d'])
        initial[root*a['G']:(root+1)*a['G']] = 1/np.sqrt(a['G'])
        psi = history(duration) @ initial
        probabilities = np.abs(psi)**2
        direct = float(np.sum(a['q']*probabilities))
        position = probabilities.sum(axis=0).reshape(a['n'], a['G']).sum(axis=1)
        memory_population = probabilities.sum(axis=1)
        member_population = np.array([sum(memory_population[z]*((z >> v) & 1)
                                          for z in range(a['modes'])) for v in range(a['n'])])
        independent = float(position @ member_population)
        self.assertGreater(abs(direct-independent), 1e-3)
        coherent_read = dephased_read = 0.0
        for z in range(a['modes']):
            _, _, ks, _ = operators(a['h'], np.diag(a['q'][z]), tau)
            coherent_read += float(np.linalg.norm(ks[1] @ psi[z])**2)
            if memory_population[z] > 1e-25:
                conditional = psi[z]/np.sqrt(memory_population[z])
                dephased_read += float(memory_population[z]*np.linalg.norm(ks[1] @ conditional)**2)
        self.assertAlmostEqual(coherent_read, dephased_read, places=13)
        amplitude = np.sqrt(memory_population)
        coherence_distance = float(np.abs(np.linalg.eigvalsh(
            np.outer(amplitude, amplitude)-np.diag(memory_population))).sum())
        self.assertGreater(coherence_distance, 0.1)
        reduced_DG = psi.T @ psi.conj()
        free = exponential(float(duration)*a['h']) @ initial
        disturbance = float(np.abs(np.linalg.eigvalsh(reduced_DG-np.outer(free, free.conj()))).sum())
        self.assertGreater(disturbance, 1e-3)
        OBS['interpretation_boundaries'] = dict(history_time=str(duration),
            genuine_joint_membership_probability=old.short(direct),
            independent_DG_and_memory_replacement_probability=old.short(independent),
            actual_read_probability=old.short(coherent_read),
            dephased_but_correlated_read_probability=old.short(dephased_read),
            full_state_coherence_trace_difference=old.short(coherence_distance),
            DG_disturbance_relative_to_free_wait=old.short(disturbance),
            coherent_membership_necessary_for_these_statistics=False,
            passive_undisturbed_DG_history_claimed=False)

    def test_06_reversibility_is_not_a_permanent_visited_set(self):
        # Boundary example: a fixed occupied site, no transport. Existing local
        # memory mechanism rotates twice and erases its Z mark, without data loss.
        x = np.array([[0, 1], [1, 0]])
        first = exponential(np.pi/2*x) @ np.array([1, 0])
        second = exponential(np.pi*x) @ np.array([1, 0])
        self.assertAlmostEqual(abs(first[1])**2, 1)
        self.assertLess(abs(second[1])**2, 1e-28)
        certificates = []
        for marks, eta in ((10, Q(1, 10)), (100, Q(1, 100))):
            needed = (1-eta)*marks
            certificates.append(dict(required_marks=marks, success_at_least=str(1-eta),
                necessary_single_excitation_abs_g_times_T=str(needed)))
        OBS['resource_and_semantic_scope'] = dict(permanent_visit_bit_claimed=False,
            isolation_required_to_freeze_member_populations=True,
            autonomous_isolation_or_clock_derived=False, necessary_budgets=certificates,
            preassigned_membership_table_used=False,
            selected_source_initial_state_and_empty_memory_are_inputs=True,
            collective_cat_source_still_requires_508_contract_or_independent_input=True,
            no_spatial_dimension_selection=True)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(round=510, scientific_baseline_round=509,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'coherent_graph_mean_obstruction.py', 'distributed_role_reader.py',
            'role_covariant_graph_obstruction.py', 'branching_tree_distance_audit.py',
            'research_note_508.md', 'research_note_509.md')},
        scope=dict(original_h_continues_during_label_formation=True,
            member_labels_generated_by_actual_finite_interaction=True,
            size_uniform_formation_mean_and_tail_bounds=True,
            complete_unknown_memory_collective_instrument=True,
            no_preassigned_membership_table=True,
            DG_free_evolution_preserved_during_writing=False,
            permanent_or_current_topology_membership_claimed=False,
            quantum_member_coherence_required_for_read_statistics=False,
            pure_original_swap_only_implementation=False,
            autonomous_control_sources_derived=False,
            multiple_complete_macro_subjects_generated=False,
            dimension_three_generated=False, full_GR_goal_completed=False,
            phase_closure_triggered=False), observations=OBS)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    data = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(data, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(data, ensure_ascii=False, indent=2))
