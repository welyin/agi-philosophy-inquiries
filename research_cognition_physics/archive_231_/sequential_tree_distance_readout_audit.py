"""Round 472: finite sequential data instruments read initial tree distances.

Independent scientific baseline 470. The graph-flip Hamiltonian stays on.
Fresh local data preparation, intermediate readout/replacement, timing and an
independent repeated graph-state source are additional operation contracts.
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

TARGET = Path(__file__).with_name('sequential_tree_distance_readout_audit_results.json')
OBS = {}
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def short(x):
    return float(f'{float(x):.12g}')


def tensor(values):
    out = np.ones((1, 1), dtype=complex)
    for value in values:
        out = np.kron(out, value)
    return out


def trace_norm(a):
    assert np.linalg.norm(a-a.conj().T) < 2e-9
    return float(np.abs(np.linalg.eigvalsh((a+a.conj().T)/2)).sum())


def reference_marginal(a, g=6, r=3):
    return a.reshape(g, r, g, r).trace(axis1=0, axis2=2)


def apply_map(superop, rho, g=6, r=3):
    return np.einsum('ghij,irjs->grhs', superop,
                     rho.reshape(g, r, g, r), optimize=True).reshape(g*r, g*r)


def jordan(n):
    g = len(n)
    identity = np.eye(g)
    return (np.einsum('gi,jh->ghij', n, identity)
            + np.einsum('gi,jh->ghij', identity, n))/2


def remainder(x):
    return math.expm1(x)-x


def all_paths(n, a, b):
    middle = [v for v in range(n) if v not in (a, b)]
    for count in range(n-1):
        for inside in itertools.permutations(middle, count):
            yield (a,)+inside+(b,)


def path_edges(path):
    return tuple(old.edge(a, b) for a, b in zip(path[:-1], path[1:]))


@lru_cache(None)
def system():
    trees, _, h, f, distances = old.six_vertex_sector()
    h = np.rint(h.real).astype(np.int64)
    f = np.rint(f.real).astype(np.int64)
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    raw = np.array([[1, 1j, 2], [2, -1j, 1], [1+1j, 2, 0],
                    [-1, 1, 2j], [2, 1j, 1-1j], [1, -2, 1j]], dtype=complex)
    vector = raw.reshape(-1)/math.sqrt(int(np.vdot(raw, raw).real))
    rho = np.outer(vector, vector.conj())
    return trees, h, f, distances, eigenvalues, eigenvectors, rho


@lru_cache(None)
def unitary(t):
    _, _, _, _, eigenvalues, eigenvectors, _ = system()
    return (eigenvectors*np.exp(-1j*float(t)*eigenvalues))@eigenvectors.conj().T


def probe_columns(n, u, v):
    """Pure preparation columns, with their mixture weights already included."""
    rest = [a for a in range(n) if a not in (u, v)]
    rank = 2**len(rest)
    out = np.zeros((2**n, rank), dtype=complex)
    for beta in range(rank):
        base = sum(((beta>>(len(rest)-1-k)) & 1) << (n-1-a)
                   for k, a in enumerate(rest))
        for ub, vb in itertools.product((0, 1), repeat=2):
            index = base | (ub << (n-1-u)) | (vb << (n-1-v))
            out[index, beta] = (1j**vb)/(2*math.sqrt(rank))
    return out


@lru_cache(None)
def instrument(edge, t):
    """Exact original-U Kraus instrument, with w=-z for J=+1."""
    u, v = edge
    columns = probe_columns(6, u, v)
    # beta, final data, output graph, input graph.
    kraus = np.einsum('dgfi,fb->bdgi', unitary(t).reshape(64, 6, 64, 6),
                      columns, optimize=True)
    z = np.array([1-2*((x>>(5-u)) & 1) for x in range(64)])
    pieces = []
    for w in (-1, 1):
        selected = kraus[:, -z == w]
        pieces.append(np.einsum('bdgi,bdhj->ghij', selected,
                                selected.conj(), optimize=True))
    return pieces[0], pieces[1], kraus


def graph_projection(edge):
    return np.diag([int(edge in tree) for tree in system()[0]])


def signed_reference_for_path(path, t):
    state = system()[-1].copy()
    for edge in path_edges(path):
        minus, plus, _ = instrument(edge, t)
        state = apply_map((plus-minus)/float(t), state)
    return reference_marginal(state)


class Audit(unittest.TestCase):
    def close(self, a, b, tol=2e-10):
        self.assertLess(np.linalg.norm(a-b), tol)

    def test_01_exact_local_coefficients_and_complete_first_map(self):
        n, u, v = 6, 0, 1
        eta_raw = tensor([I2+X, I2+Y]+[I2]*4)
        za = tensor([Z]+[I2]*5)
        self.assertEqual(np.trace(za@eta_raw), 0)
        rows = []
        for edge in itertools.combinations(range(n), 2):
            swap = old.core.swap(n, *edge)
            left = np.trace(za@swap@eta_raw)
            right = np.trace(za@eta_raw@swap)
            expected = -32j if edge == (u, v) else 0
            self.assertEqual(left, expected)
            self.assertEqual(right, expected.conjugate() if expected else 0)
            self.assertEqual(np.trace(swap@eta_raw), 32)
            rows.append(dict(edge=list(edge), left_trace_over_64=[short(left.real/64),
                                                                         short(left.imag/64)]))
        _, h, _, _, _, _, _ = system()
        prep = np.kron(eta_raw, np.eye(6))
        za_big = np.kron(za, np.eye(6))
        left = (za_big@h@prep).reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)
        right = (za_big@prep@h).reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)
        projection = graph_projection((u, v))
        self.assertTrue(np.array_equal(left, -32j*projection))
        self.assertTrue(np.array_equal(right, 32j*projection))
        # Weighted record w=-z reverses the sign of the measured Z derivative.
        for i, j in itertools.product(range(6), repeat=2):
            matrix_unit = np.zeros((6, 6)); matrix_unit[i, j] = 1
            actual = 1j*(left@matrix_unit-matrix_unit@right)/64
            self.close(actual, (projection@matrix_unit+matrix_unit@projection)/2, 1e-14)
        OBS['exact_first_map'] = dict(
            all_15_data_edge_coefficients=rows,
            graph_matrix_units_checked=36, full_H_dimension=384,
            instrument_record='w=-sgn(J)*z, data preparation +X_u,+Y_v,I/2 elsewhere',
            signed_derivative='abs(J)*{n_uv,rho}/2',
            graph_F_first_signed_term_zero=True,
            every_SWAP_preparation_mean='1/2',
            mean_data_H='J*(N-1)*I_G/2')

    def test_02_finite_time_actual_CP_instrument_and_Jordan_scope(self):
        edge, t = (0, 1), Q(1, 1024)
        minus, plus, kraus = instrument(edge, t)
        completeness = np.einsum('bdgi,bdgj->ij', kraus.conj(), kraus, optimize=True)
        self.close(completeness, np.eye(6))
        z = np.array([1-2*((x>>5) & 1) for x in range(64)])
        smallest = []
        for w in (-1, 1):
            selected = kraus[:, -z == w].reshape(-1, 36)
            choi = selected.T@selected.conj()
            eigen_min = float(np.linalg.eigvalsh(choi).min())
            self.assertGreaterEqual(eigen_min, -2e-12)
            smallest.append(short(eigen_min))
        n = graph_projection(edge)
        source = next(i for i in range(6) if n[i, i])
        other = next(i for i in range(6) if not n[i, i])
        unit = np.zeros((6, 6)); unit[source, other] = 1
        self.assertTrue(np.array_equal((n@unit+unit@n)/2, unit/2))
        self.assertFalse(np.any(n@unit@n))
        first_error = np.linalg.norm((plus-minus)/float(t)-jordan(n))
        OBS['actual_instrument'] = dict(
            time=str(t), Kraus_count_per_probe=int(np.prod(kraus.shape[:2])),
            Choi_minimum_eigenvalues=smallest,
            trace_preservation_residual=short(np.linalg.norm(completeness-np.eye(6))),
            normalized_signed_superoperator_Frobenius_error=short(first_error),
            finite_time_instrument_from_full_original_U=True,
            no_graph_projector_measurement_or_postselection=True,
            off_diagonal_Jordan_matrix_unit_factor='1/2',
            corresponding_Lueders_factor='0',
            Jordan_sequence_not_a_Lueders_projection_channel=True)

    def test_03_path_records_unknown_reference_and_full_mean_distance(self):
        trees, _, _, distances, _, _, rho = system()
        paths = [p for p in all_paths(6, 0, 5)
                 if any(set(path_edges(p)) <= tree for tree in trees)]
        self.assertEqual(len(paths), 4)
        t, h = Q(1, 65536), 25
        delta = remainder(2*h*float(t))/float(t)
        total_readout = np.zeros((3, 3), dtype=complex)
        total_bound = 0.0
        rows = []
        for path in paths:
            length = len(path)-1
            indicator = np.diag([int(set(path_edges(path)) <= tree) for tree in trees])
            target = reference_marginal(np.kron(indicator, np.eye(3))@rho)
            actual = signed_reference_for_path(path, t)
            error = trace_norm(actual-target)
            bound = (1+delta)**length-1
            self.assertLess(error, bound)
            self.assertGreaterEqual(np.linalg.eigvalsh(target).min(), -1e-12)
            total_readout += length*actual
            total_bound += length*bound
            rows.append(dict(path=list(path), length=length,
                true_probability=short(np.trace(target).real),
                normalized_actual_record_mean=short(np.trace(actual).real),
                reference_trace_norm_error=short(error), proven_bound=short(bound)))
        true_distance = reference_marginal(np.kron(np.diag(distances), np.eye(3))@rho)
        error = trace_norm(total_readout-true_distance)
        self.assertLess(error, total_bound)
        self.assertGreater(error, 1e-9)
        graph = rho.reshape(6, 3, 6, 3).trace(axis1=1, axis2=3)
        self.assertGreater(np.linalg.norm(graph-np.diag(np.diag(graph))), .1)
        OBS['sequential_distance_readout'] = dict(
            time=str(t), full_H_J=1, full_H_kappa=1,
            finite_N_H_bound=h, arbitrary_reference_sample_dimension=3,
            path_records=rows,
            mean_distance_target=short(np.trace(true_distance).real),
            mean_distance_actual_estimator=short(np.trace(total_readout).real),
            complete_reference_distance_moment_error=short(error),
            weighted_proven_bound=short(total_bound),
            no_use_of_six_tree_shortcut_3_minus_Q=True,
            graph_coherence_retained=True)

    def test_04_tree_path_partition_and_adjacency_marginal_boundary(self):
        # Larger fixed-degree department, without simulating its huge Hilbert space.
        n = 8
        trees = {old.prufer_tree(n, p)
                 for p in set(itertools.permutations((1, 1, 2, 2, 3, 3)))}
        self.assertEqual(len(trees), 90)
        paths = list(all_paths(n, 0, 7))
        checked = 0
        for tree in trees:
            events = [p for p in paths if set(path_edges(p)) <= tree]
            self.assertEqual(len(events), 1)
            self.assertEqual(len(events[0])-1, len(old.path_between(tree, n, 0, 7))-1)
            checked += len(paths)
        six_trees, _, _, distances, _, _, _ = system()
        complement = {}
        for i, tree in enumerate(six_trees):
            switched = frozenset(old.edge(2 if a == 1 else 1 if a == 2 else a,
                                         2 if b == 1 else 1 if b == 2 else b)
                                  for a, b in tree)
            complement[i] = six_trees.index(switched)
        close_index = next(i for i, d in enumerate(distances) if d == 2)
        far_index = next(i for i, d in enumerate(distances) if d == 3)
        close = np.zeros(6); far = np.zeros(6)
        for out, index in ((close, close_index), (far, far_index)):
            out[index] += .5; out[complement[index]] += .5
        for edge in itertools.combinations(range(6), 2):
            values = np.array([edge in tree for tree in six_trees], dtype=int)
            self.assertEqual(float(close@values), float(far@values))
        self.assertEqual((float(close@distances), float(far@distances)), (2., 3.))
        OBS['path_polynomial_scope'] = dict(
            larger_N=n, fixed_degree_trees=90,
            paths_per_tree=len(paths), path_event_checks=checked,
            exactly_one_simple_path_event_for_every_tree=True,
            all_adjacency_marginals_equal_but_mean_distance_different=[2, 3],
            close_distribution=close.tolist(), far_distribution=far.tolist(),
            statistical_mean_distance_not_a_single_graph_eigenvalue=True,
            no_spatial_dimension_inferred=True)

    def test_05_nonselective_and_retained_record_backaction(self):
        _, _, f, _, _, _, rho = system()
        t, v = Q(1, 128), Q(15, 2)
        edges = ((0, 1), (1, 2), (2, 5))
        branches = [rho.copy()]
        for edge in edges:
            minus, plus, _ = instrument(edge, t)
            branches = [apply_map(piece, state)
                        for state in branches for piece in (minus, plus)]
        m = len(edges)
        eigenvalues, eigenvectors = np.linalg.eigh(f)
        free_u = (eigenvectors*np.exp(-1j*float(m*t)*eigenvalues))@eigenvectors.conj().T
        free_big = np.kron(free_u, np.eye(3))
        baseline = free_big@rho@free_big.conj().T
        marginal = sum(branches)
        marginal_error = trace_norm(marginal-baseline)/2
        record_error = sum(trace_norm(branch-baseline/(2**m)) for branch in branches)/2
        marginal_bound = min(1, m*remainder(float(2*v*t))/2)
        record_bound = min(1, float(m*v*t))
        self.assertGreater(marginal_error, 1e-8)
        self.assertLess(marginal_error, marginal_bound)
        self.assertLess(record_error, record_bound)
        self.close(reference_marginal(marginal), reference_marginal(rho))
        self.assertTrue(all(np.linalg.eigvalsh((branch+branch.conj().T)/2).min()>-1e-12
                            for branch in branches))
        OBS['backaction'] = dict(
            time_per_probe=str(t), probes=m, classical_record_branches=2**m,
            centered_interaction_bound=str(v),
            graph_reference_trace_distance=short(marginal_error),
            nonselective_bound=short(marginal_bound),
            record_graph_reference_trace_distance=short(record_error),
            retained_record_bound=short(record_bound),
            comparison='free exp(-i*kappa*F*m*t) graph evolution',
            original_reference_marginal_preserved=True,
            graph_reference_state_not_preserved=True,
            conditional_state_not_claimed_undisturbed=True)

    def test_06_finite_statistics_instrument_errors_and_resource_ledger(self):
        n, h, j, m = 6, Q(25), Q(1), 5
        paths = list(all_paths(n, 0, n-1))
        count = len(paths)
        weight = sum(len(path)-1 for path in paths)
        self.assertEqual((count, weight), (65, 261))
        distance_error = Q(1, 10)
        per_path_error = distance_error/weight
        t = per_path_error*j/(48*m*h*h)
        gamma = per_path_error*j*t/(8*m)
        x = 2*h*t
        self.assertLess(x, 1)
        # exp(x)-1-x <= x^2/[2(1-x/3)], by successive Taylor-term ratios.
        r_bound = x*x/(2*(1-x/3))
        delta = (r_bound+gamma)/(j*t)
        bias = (1+delta)**m-1
        self.assertLess(bias, per_path_error/2)
        self.assertLess(weight*bias, distance_error/2)
        # exp(10)>13000=2*65/(1/100), certified by a finite positive sum.
        self.assertGreater(sum(Q(10)**k/Q(math.factorial(k)) for k in range(30)), 13000)
        copies = math.ceil(80/(per_path_error**2*(j*t)**(2*m)))
        self.assertGreaterEqual(Q(copies), 80/(per_path_error**2*(j*t)**(2*m)))
        control_error, control_duration = gamma/2, gamma/(4*h)
        self.assertEqual(control_error+2*h*control_duration, gamma)
        # An explicit readout bit-flip error rescales the signed map, not the graph.
        flip_probability = Q(1, 1000)
        minus, plus, _ = instrument((0, 1), Q(1, 1024))
        noisy_minus = float(1-flip_probability)*minus+float(flip_probability)*plus
        noisy_plus = float(1-flip_probability)*plus+float(flip_probability)*minus
        self.close(noisy_plus-noisy_minus, float(1-2*flip_probability)*(plus-minus))
        OBS['finite_resource_certificate'] = dict(
            N=n, C=3, J=str(j), kappa=1, H_bound=str(h),
            catalog_paths=count, sum_path_lengths=weight, maximum_probes_per_trial=m,
            target_distance_error=str(distance_error), failure_probability='1/100',
            per_path_total_error=str(per_path_error), time_per_probe=str(t),
            complete_instrument_diamond_error_budget=str(gamma),
            actual_stage_instruments_fixed_independently_of_previous_random_records=True,
            intrinsic_control_error_budget=str(control_error),
            total_control_duration_budget=str(control_duration),
            per_path_bias_bound=str(bias),
            independent_copies_per_path=str(copies),
            total_independent_graph_inputs=str(count*copies),
            maximum_probe_stages=str(m*count*copies),
            maximum_local_data_replacements=str(n*m*count*copies),
            exp10_greater_than_13000_exactly_certified=True,
            finite_source_required_not_cloned_or_reused_without_proof=True,
            fresh_probe_storage_per_trial_at_most_n_times_m=n*m,
            repeated_readout_of_one_unknown_graph_not_independent_samples=True,
            fixed_nonzero_hardware_error_cannot_be_removed_by_t_to_zero=True,
            actual_readout_bit_flip_signed_factor=str(1-2*flip_probability),
            statistics_and_timing_numbers_not_experimental_runs=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=472, baseline_round=470, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(original_graph_F_always_on_during_waits=True,
            actual_intermediate_data_instruments_are_new_permissions=True,
            ideal_instantaneous_control_not_autonomously_realized=True,
            complete_instrument_error_includes_control_duration=True,
            actual_instrument_no_feedback_from_previous_records=True,
            fixed_step_channels_trivial_on_prior_records_and_R=True,
            noninteracting_retired_apparatus_required=True,
            arbitrary_initial_graph_coherence_and_reference=True,
            finite_N_tree_mean_distance_accessible_under_stated_permissions=True,
            signed_Jordan_sequence_only_yields_projector_reference_marginal_after_trace_G=True,
            no_graph_measurement_no_postselection_no_free_cloning=True,
            unknown_original_data_not_preserved_by_a_bare_reset=True,
            original_single_wait_observability_problem_solved=False,
            physical_spatial_dimension_derived=False,
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
