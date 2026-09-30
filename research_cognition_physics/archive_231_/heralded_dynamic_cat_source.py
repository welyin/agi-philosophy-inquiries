"""Round 508: product-source cat filtering while the unknown graph keeps moving.

Declared inputs: graph-controlled pair-parity detectors, product |+> meters,
finite scheduled couplings, isolation and final detector records. Success is
heralded, exponentially rare, and is not a deterministic quantum channel.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import io
import itertools
import json
from pathlib import Path
import platform
import unittest

import numpy as np

import coherent_graph_mean_obstruction as old
from distributed_role_reader import exponential, opnorm

HERE = Path(__file__).resolve().parent
TARGET = HERE/'heralded_dynamic_cat_source_results.json'
OBS = {}


@lru_cache(maxsize=1)
def finite_model():
    trees, _, _, _, h, _, _, _ = old.model(2)
    n, d, graph_count = 6, len(h), len(trees)
    pairs = list(itertools.combinations(range(n), 2))
    bound = 10
    layer_time = Q(1, 4*bound)
    tau = layer_time/len(pairs)
    theta = np.pi/3
    free = exponential(float(tau)*h)
    pair_ops = []
    for edge in pairs:
        present = np.array([edge in tree for tree in trees], dtype=np.int64)
        projector = np.diag(np.tile(present, n))
        plus = exponential(float(tau)*h+theta*projector)
        minus = exponential(float(tau)*h-theta*projector)
        no_click = (plus+minus)/2
        click = 1j*(plus-minus)/2
        pair_ops.append((projector, no_click, click))
    layers, ideals = [], []
    for bits in range(2**n):
        actual = np.eye(d, dtype=complex)
        ideal = np.eye(d)
        for edge, (projector, no_click, _) in zip(pairs, pair_ops):
            differs = ((bits >> edge[0]) ^ (bits >> edge[1])) & 1
            actual = (no_click if differs else free) @ actual
            if differs:
                ideal = (np.eye(d)-projector/2) @ ideal
        layers.append(actual)
        ideals.append(ideal)
    return dict(trees=trees, h=h, n=n, d=d, graph_count=graph_count,
                pairs=pairs, bound=bound, T=layer_time, tau=tau,
                free=free, pair_ops=pair_ops, layers=layers, ideals=ideals)


def finite_budget(i, error):
    n = 2*i+2
    p0 = Q(1, 2**(n-1))
    q = Q(3, 4)
    layers = 0
    while 4*(1-p0)*q**(2*layers) > error**2*p0:
        layers += 1
    pairs = n*(n-1)//2
    bound = 4*i+2
    return dict(I=i, N=n, target_conditional_trace_error=str(error),
                success_floor=str(p0), layers=layers, pair_count=pairs,
                detector_qubits=layers*pairs, record_bits=layers*pairs,
                cat_meter_qubits=n, layer_time=str(Q(1, 4*bound)),
                total_coupling_time=str(Q(layers, 4*bound)),
                pair_time=str(Q(1, 4*bound*pairs)),
                coupling_over_pi=str(Q(4*bound*pairs, 3)),
                exact_squared_error_bound=str(4*(1-p0)*q**(2*layers)/p0))


class Audit(unittest.TestCase):
    def test_01_common_dark_space_on_all_small_trees(self):
        rows = []
        for i in (2, 3, 4):
            trees = old.previous.ensemble(i)[0]
            n = 2*i+2
            patterns = np.arange(2**n, dtype=np.int64)
            minimum = n
            for tree in trees:
                cuts = np.zeros(2**n, dtype=np.int64)
                for u, v in tree:
                    cuts += ((patterns >> u) ^ (patterns >> v)) & 1
                self.assertTrue(np.array_equal(np.flatnonzero(cuts == 0), [0, 2**n-1]))
                minimum = min(minimum, int(cuts[1:-1].min()))
            self.assertEqual(minimum, 1)
            rows.append(dict(I=i, N=n, graph_count=len(trees),
                             tested_graph_pattern_pairs=len(trees)*2**n,
                             common_kernel_patterns=2, minimum_positive_cut_count=minimum))
        OBS['common_kernel'] = rows

    def test_02_finite_detector_blocks_are_complete(self):
        a = finite_model()
        d, h = a['d'], a['h']
        eye = np.eye(d)
        y = np.array([[0, -1j], [1j, 0]])
        complete_error = block_error = duhamel_error = 0.0
        for projector, no_click, click in a['pair_ops']:
            full = exponential(float(a['tau'])*np.kron(h, np.eye(2))
                               +np.pi/3*np.kron(projector, y))
            block_error = max(block_error, opnorm(full[0::2, 0::2]-no_click),
                              opnorm(full[1::2, 0::2]-click))
            complete_error = max(complete_error,
                opnorm(no_click.conj().T @ no_click+click.conj().T @ click-eye))
            duhamel_error = max(duhamel_error, opnorm(no_click-(eye-projector/2)))
        self.assertLess(block_error, 2e-13)
        self.assertLess(complete_error, 2e-13)
        self.assertLess(duhamel_error, float(a['tau'])*a['bound'])
        self.assertLessEqual(opnorm(h), a['bound']+1e-12)
        OBS['finite_detector'] = dict(DG_dimension=d, joint_detector_dimension=2*d,
            pair_count=len(a['pairs']), pair_time=str(a['tau']),
            maximum_block_error=old.short(block_error),
            maximum_completeness_error=old.short(complete_error),
            maximum_drift_difference=old.short(duhamel_error),
            drift_difference_bound=str(a['tau']*a['bound']))

    def test_03_dark_free_evolution_and_layer_contraction(self):
        a = finite_model()
        free = exponential(float(a['T'])*a['h'])
        dark_error = max(opnorm(a['layers'][j]-free) for j in (0, 2**a['n']-1))
        comparison = max(opnorm(k-b) for k, b in zip(a['layers'], a['ideals']))
        complement_norm = max(opnorm(k) for k in a['layers'][1:-1])
        benchmark_norm = max(opnorm(k) for k in a['ideals'][1:-1])
        self.assertLess(dark_error, 2e-12)
        self.assertLess(comparison, 1/4)
        self.assertLess(complement_norm, 3/4)
        self.assertAlmostEqual(benchmark_norm, 1/2)
        self.assertLess(max(opnorm(a['layers'][j]-a['layers'][-1-j])
                            for j in range(2**a['n'])), 1e-13)
        OBS['layer'] = dict(layer_time=str(a['T']), h_norm_bound=a['bound'],
            dark_free_evolution_error=old.short(dark_error),
            actual_vs_benchmark_norm=old.short(comparison),
            telescoping_upper='1/4', actual_complement_norm=old.short(complement_norm),
            ideal_complement_norm=old.short(benchmark_norm), certified_contraction='3/4',
            global_bit_flip_evenness_preserved=True)

    def test_04_success_preserves_arbitrary_old_reference(self):
        a = finite_model()
        n, d = a['n'], a['d']
        p0 = Q(1, 2**(n-1))
        rows = []
        for layers in (1, 4, 12, 24):
            powers = [np.linalg.matrix_power(k, layers) for k in a['layers']]
            residual_gram = sum(k.conj().T @ k for k in powers[1:-1])/2**n
            residual_squared = max(0, float(np.linalg.eigvalsh(residual_gram)[-1]))
            success_gram = sum(k.conj().T @ k for k in powers)/2**n
            spectrum = np.linalg.eigvalsh(success_gram)
            r_squared = (1-p0)*Q(3, 4)**(2*layers)
            self.assertLessEqual(residual_squared, float(r_squared)+2e-14)
            self.assertGreaterEqual(spectrum[0], float(p0)-2e-13)
            self.assertLessEqual(spectrum[-1], float(p0+r_squared)+2e-13)
            self.assertLess(opnorm(success_gram-float(p0)*np.eye(d)-residual_gram), 3e-13)
            # Maximally entangled old DG/R is a numerical witness, not a diamond proof.
            average_residual = float(np.trace(residual_gram).real)/d
            conditional_reference_distance = 2*np.sqrt(average_residual/(float(p0)+average_residual))
            uniform_bound = 2*np.sqrt(float(r_squared/p0))
            self.assertLessEqual(conditional_reference_distance, uniform_bound+1e-13)
            rows.append(dict(layers=layers, success_probability_min=old.short(spectrum[0]),
                success_probability_max=old.short(spectrum[-1]),
                residual_operator_norm_squared=old.short(residual_squared),
                residual_squared_bound=str(r_squared),
                maximally_entangled_reference_trace_error=old.short(conditional_reference_distance),
                uniform_conditional_trace_error_upper=old.short(min(2, uniform_bound)),
                unnormalized_success_diamond_upper=old.short(
                    2*np.sqrt(float(p0*r_squared))+float(r_squared))))
        OBS['success_and_reference'] = dict(success_floor=str(p0), cases=rows,
            arbitrary_reference_guarantee_from_operator_proof=True,
            normalized_conditional_map_claimed_linear=False)

    def test_05_failure_records_are_complete_and_can_disturb_graph(self):
        a = finite_model()
        d = a['d']
        maximum = 0.0
        # First-click classes form a complete instrument; all later outcomes can
        # be retained and then grouped. No failed branch is silently discarded.
        for bits in range(2**a['n']):
            preceding = np.eye(d, dtype=complex)
            failure_effect = np.zeros((d, d), dtype=complex)
            for edge, (_, no_click, click) in zip(a['pairs'], a['pair_ops']):
                differs = ((bits >> edge[0]) ^ (bits >> edge[1])) & 1
                if differs:
                    branch = click @ preceding
                    failure_effect += branch.conj().T @ branch
                preceding = (no_click if differs else a['free']) @ preceding
            maximum = max(maximum, opnorm(failure_effect+preceding.conj().T @ preceding-np.eye(d)))
        self.assertLess(maximum, 2e-12)
        j = next(j for j, (p, _, _) in enumerate(a['pair_ops']) if 0 < np.trace(p) < d)
        projector, _, click = a['pair_ops'][j]
        present = int(np.flatnonzero(np.diag(projector)[:a['graph_count']] == 1)[0])
        absent = int(np.flatnonzero(np.diag(projector)[:a['graph_count']] == 0)[0])
        psi = np.zeros(d)
        psi[present] = psi[absent] = 1/np.sqrt(2)
        changed = click @ psi
        chance = float(np.vdot(changed, changed).real)
        changed /= np.sqrt(chance)
        target = a['free'] @ psi
        after = float(np.vdot(changed, projector @ changed).real)
        ordinary = float(np.vdot(target, projector @ target).real)
        self.assertGreater(after-ordinary, 0.49)
        OBS['failure_accounting'] = dict(first_click_partition_completeness_error=old.short(maximum),
            witness_edge=a['pairs'][j], witness_is_single_selected_pair_check=True,
            odd_meter_sector_click_probability=old.short(chance),
            product_plus_single_check_click_probability=old.short(chance/2),
            edge_probability_after_click=old.short(after),
            edge_probability_after_free_wait=old.short(ordinary),
            free_retry_preserving_old_graph_justified=False,
            final_deferred_detector_reads_allowed=True)

    def test_06_finite_exact_budgets_and_composition_scope(self):
        rows = [finite_budget(i, eps) for i, eps in
                ((2, Q(1, 10)), (3, Q(1, 100)), (10, Q(1, 100)))]
        for row in rows:
            self.assertLessEqual(Q(row['exact_squared_error_bound']),
                                 Q(row['target_conditional_trace_error'])**2)
            self.assertEqual(row['record_bits'], row['layers']*row['pair_count'])
        OBS['resource_budget'] = dict(cases=rows, read_507_unnormalized_composition=
            '2 sqrt(p0) r + r^2 + p0 pi tau_read',
            read_507_conditional_joint_trace_bound='2 r/sqrt(p0) + pi tau_read',
            finite_508_process_exactly_preserves_506_code=False,
            product_source_phase_calibration_is_input=True,
            pair_addressing_schedule_and_blank_detectors_are_inputs=True,
            all_interaction_idle_time_must_be_in_drift_budget=True,
            dimension_selection_claimed=False)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output, verbosity=0).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(output.getvalue())
    return dict(round=508, scientific_baseline_round=507,
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__,
        dependency_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'coherent_graph_mean_obstruction.py', 'distributed_role_reader.py',
            'role_covariant_graph_obstruction.py', 'branching_tree_distance_audit.py', 'research_note_507.md')},
        scope=dict(original_h_continues_during_preparation=True,
            common_kernel_independent_of_unknown_graph=True,
            product_source_heralded_cat_preparation=True,
            arbitrary_old_reference_uniform_error=True,
            finite_detector_unitary_and_all_failure_records_explicit=True,
            deterministic_preparation_claimed=False, free_nondisturbing_retry_claimed=False,
            pure_original_swap_only_implementation=False,
            exact_506_code_preservation_by_entire_preparation=False,
            autonomous_timing_or_pair_routing_derived=False,
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
