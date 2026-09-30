"""Round 474: three relational data components change actual distance means.

Scientific baseline 472. Source and internal directional-reference preparation,
fixed initial graph, readout instruments, clocks, and repeated sources remain
inputs. A rank-three response is not a derivation of spatial dimension.
"""
import argparse
from fractions import Fraction as FQ
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
import sequential_tree_distance_readout_audit as readout

TARGET = Path(__file__).with_name('relational_direction_distance_response_results.json')
OBS = {}
I = np.eye(2, dtype=complex)
X, Y, Z = old.core.PAULI
P = (X, Y, Z)
LEAVES = (0, 3, 4, 5)


def short(x):
    return float(f'{float(x):.11g}')


def tensor(values):
    out = np.ones((1, 1), dtype=complex)
    for value in values:
        out = np.kron(out, value)
    return out


def comm(a, b):
    return a@b-b@a


def data_basis_raw():
    # eta(r) = (eta_raw[0]+sum r_mu eta_raw[mu+1])/64.
    return [tensor([I, source, I, I+X, I+Y, I+Z])
            for source in (I, X, Y, Z)]


def data_state(r):
    raw = data_basis_raw()
    return (raw[0]+sum(float(x)*a for x, a in zip(r, raw[1:])))/64


@lru_cache(None)
def system():
    trees, _, h, f, _ = old.six_vertex_sector()
    h = np.rint(h.real).astype(np.int64)
    f = np.rint(f.real).astype(np.int64)
    population = np.array([[int(old.edge(1, j) in tree) for tree in trees]
                           for j in LEAVES], dtype=np.int64)
    q = population[0]-population[1:]
    eigenvalues, eigenvectors = np.linalg.eigh(h)
    return trees, h, f, population, q, eigenvalues, eigenvectors


@lru_cache(None)
def unitary(t):
    *_, values, vectors = system()
    return (vectors*np.exp(-1j*float(t)*values))@vectors.conj().T


@lru_cache(None)
def effects(t):
    q = system()[4]
    initial = np.kron(np.eye(64), np.ones((6, 1))/math.sqrt(6))
    v = unitary(t)@initial
    return [v.conj().T@(np.tile(diag, 64)[:, None]*v) for diag in q]


def response(t):
    raw = data_basis_raw()
    e = effects(t)
    values = np.array([[np.trace(a@b).real/64 for b in raw] for a in e])
    return values[:, 0], values[:, 1:]


def graph_state(r, t):
    initial = np.kron(np.eye(64), np.ones((6, 1))/math.sqrt(6))
    v = unitary(t)@initial
    out = v@data_state(r)@v.conj().T
    return out.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)


@lru_cache(None)
def collective_generators():
    return [sum(tensor([pauli if k == j else I for k in range(6)])
                for j in range(6))/2 for pauli in P]


def casimir(a):
    return sum(comm(j, comm(j, a)) for j in collective_generators())


def restricted_twirl(a):
    # Valid for the four-active-qubit family, not the full six-qubit algebra.
    out = a.copy()
    for ell in range(1, 5):
        out = out-casimir(out)/(ell*(ell+1))
    return out


class Audit(unittest.TestCase):
    def close(self, a, b, tol=3e-10):
        self.assertLess(float(np.linalg.norm(np.asarray(a)-np.asarray(b))), tol)

    def test_01_exact_tree_covariance_and_distance_response(self):
        trees, h, f, population, q, *_ = system()
        self.assertTrue(np.array_equal(f@np.ones(6, dtype=int), np.full(6, 4)))
        centered = 2*population-1
        gram = centered@centered.T
        self.assertTrue(np.array_equal(gram, 8*np.eye(4, dtype=int)-2))
        for j, k in itertools.product(range(4), repeat=2):
            double = comm(np.diag(population[j]), comm(f, np.diag(population[k])))
            self.assertEqual(int(double.sum()), -2*int(gram[j, k]))
        raws = data_basis_raw()
        second = []
        for diag in q:
            a = np.diag(np.tile(diag, 64))
            first = comm(h, a)
            nested = comm(h, first)
            # Every value below is a Gaussian integer before division by 384.
            values = []
            for raw in raws:
                prep = np.kron(raw, np.ones((6, 6), dtype=int))
                self.assertEqual(np.trace(a@prep), 0)
                self.assertEqual(np.trace(first@prep), 0)
                values.append(np.trace(-nested@prep))
            second.append(values)
        self.assertTrue(np.array_equal(np.asarray(second),
            np.column_stack([np.zeros(3), -512*np.eye(3)])))
        # Internal 1<->2 simultaneously complements each two-leaf assignment.
        renamed = [frozenset(old.edge(2 if a == 1 else 1 if a == 2 else a,
                                      2 if b == 1 else 1 if b == 2 else b)
                             for a, b in tree) for tree in trees]
        graph_perm = [trees.index(tree) for tree in renamed]
        data_perm = np.argmax(old.core.swap(6, 1, 2).real, axis=0)
        full_perm = np.array([6*d+g for d in data_perm for g in graph_perm])
        self.assertTrue(np.array_equal(h[np.ix_(full_perm, full_perm)], h))
        self.assertTrue(np.array_equal(q[:, graph_perm], -q))
        eta0 = raws[0]
        self.assertTrue(np.array_equal(eta0[np.ix_(data_perm, data_perm)], eta0))
        OBS['exact_response'] = dict(
            graph_count=6, total_dimension=384, covariance_diagonal='1/4',
            covariance_off_diagonal='-1/12', graph_F_uniform_eigenvalue=4,
            second_derivative_raw_numerators=np.asarray(second).real.astype(int).tolist(),
            raw_denominator=384,
            Q_definition=['D_13-D_10', 'D_14-D_10', 'D_15-D_10'],
            second_derivative='-(4/3)*kappa*J*r',
            all_time_zero_offset_has_exact_permutation_symmetry=True,
            second_order_rank_equals_affine_span_of_reference_Bloch_vectors=True,
            second_order_rank_not_claimed_as_all_time_rank=True)

    def test_02_uniform_rank_window_and_complete_unitary_inverse(self):
        _, h, _, _, _, *_ = system()
        self.assertLessEqual(float(np.linalg.norm(h, 2)), 9+1e-12)
        tmax = FQ(1, 10000)
        # e^x-1-x-x^2/2 <= x^3/[6(1-x/4)], 0<=x<4.
        x = 18*tmax
        tail = x**3/(6*(1-x/4))
        relative_remainder_bound = 3*tail/(tmax*tmax)
        self.assertLess(relative_remainder_bound, FQ(1, 3))
        bias, a = response(tmax)
        smin = float(np.linalg.svd(a, compute_uv=False)[-1])
        self.close(bias, 0, 2e-15)
        self.assertGreaterEqual(smin, float(tmax*tmax/3))
        r = np.array([.2, -.3, .4])
        direct = np.array([np.trace(e@data_state(r)).real for e in effects(tmax)])
        recovered = np.linalg.solve(a, direct)
        self.close(recovered, r, 2e-7)
        self.close(direct, bias+a@r, 2e-16)
        for t in (FQ(1, 7), FQ(1, 2)):
            self.close(response(t)[0], 0, 2e-14)
        OBS['uniform_injective_window'] = dict(
            window='0 < t <= 1/10000; J=kappa=1', Hamiltonian_norm_upper_bound=9,
            entry_tail='(18t)^3 / (6*(1-18t/4))',
            matrix_remainder_over_t_squared_at_endpoint=str(relative_remainder_bound),
            singular_value_lower_bound='t^2/3',
            endpoint_response_matrix=[[short(v) for v in row] for row in a],
            endpoint_smallest_singular_value=short(smin),
            input_r=r.tolist(), recovered_r=[short(v) for v in recovered],
            exact_affine_response_not_a_truncated_inverse=True,
            all_t_window_certified_by_rational_bound_not_numerical_sampling=True)

    def test_03_collective_twirl_keeps_all_graph_outputs(self):
        r = (.2, -.3, .4)
        eta = data_state(r)
        twirled = restricted_twirl(eta)
        self.close(twirled, twirled.conj().T)
        self.assertAlmostEqual(float(np.trace(twirled).real), 1.)
        eigmin = float(np.linalg.eigvalsh(twirled).min())
        self.assertGreaterEqual(eigmin, -2e-12)
        residual = max(float(np.linalg.norm(comm(j, twirled)))
                       for j in collective_generators())
        self.assertLess(residual, 2e-12)
        for k in range(6):
            self.close(old.core.partial(twirled, [2]*6, [k]), I/2)
        h = system()[1]
        for j in collective_generators():
            self.close(comm(h, np.kron(j, np.eye(6))), 0)
        t = FQ(1, 5)
        initial = np.kron(np.eye(64), np.ones((6, 1))/math.sqrt(6))
        v = unitary(t)@initial
        out = v@eta@v.conj().T
        out_twirl = v@twirled@v.conj().T
        graph = out.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)
        graph_twirl = out_twirl.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)
        self.close(graph, graph_twirl, 3e-13)
        # Orthogonal source states acquire overlapping supports. Hence no
        # channel can perfectly recover every unknown input after this twirl.
        rho_plus = restricted_twirl(data_state((0, 0, 1)))
        rho_minus = restricted_twirl(data_state((0, 0, -1)))
        overlap = float(np.trace(rho_plus@rho_minus).real)
        self.assertAlmostEqual(overlap, 1/96)
        remaining_distance = old.core.distance(rho_plus, rho_minus)
        self.assertLess(remaining_distance, 1.)
        OBS['internal_relational_preparation'] = dict(
            twirl_is_Haar_projection_on_four_active_qubit_family=True,
            Casimir_factors=[2, 6, 12, 20],
            not_a_projection_formula_for_arbitrary_six_qubit_operators=True,
            all_single_qubit_marginals_maximally_mixed=True,
            twirled_state_smallest_eigenvalue=short(eigmin),
            largest_collective_commutator_norm=short(residual),
            entire_graph_density_matrix_difference_norm=short(np.linalg.norm(graph-graph_twirl)),
            all_graph_measurements_preserved_analytically_for_all_t=True,
            absolute_direction_removed_but_relational_reference_resource_input=True,
            orthogonal_source_twirl_overlap=short(overlap),
            orthogonal_source_twirl_trace_distance=short(remaining_distance),
            rank_three_classical_response_not_reversible_unknown_quantum_encoding=True)

    def test_04_actual_single_edge_instruments_read_natural_output(self):
        t, probe_t = FQ(1, 5), FQ(1, 65536)
        rho = graph_state((.2, -.3, .4), t)
        true_n = system()[3]@np.diag(rho).real
        estimated = []
        for leaf in LEAVES:
            minus, plus, kraus = readout.instrument((1, leaf), probe_t)
            record_mean = np.trace(readout.apply_map(plus-minus, rho, g=6, r=1)).real
            estimated.append(record_mean/float(probe_t))
            completeness = np.einsum('bdgi,bdgj->ij', kraus.conj(), kraus, optimize=True)
            self.close(completeness, np.eye(6))
        estimated = np.asarray(estimated)
        bound = (math.expm1(50*float(probe_t))-50*float(probe_t))/float(probe_t)
        self.assertLess(float(np.max(np.abs(estimated-true_n))), bound)
        actual_q = true_n[0]-true_n[1:]
        read_q = estimated[0]-estimated[1:]
        self.assertLess(float(np.linalg.norm(read_q-actual_q)), 2*math.sqrt(3)*bound)
        OBS['actual_distance_readout'] = dict(
            inherited_instrument_round=472, natural_wait=str(t), probe_wait=str(probe_t),
            exact_adjacency_means=[short(v) for v in true_n],
            U_generated_record_estimates=[short(v) for v in estimated],
            largest_adjacency_bias=short(np.max(np.abs(estimated-true_n))),
            certified_single_edge_bias_bound=short(bound),
            true_Q=[short(v) for v in actual_q], measured_Q=[short(v) for v in read_q],
            complete_four_CP_instruments_checked=True,
            direct_graph_observable_measurement_not_assumed=True,
            old_data_correlations_transferred_to_isolated_storage_before_probes=True,
            handoff_and_all_controls_are_additional_explicit_permissions=True)

    def test_05_finite_statistical_and_handoff_resource_certificate(self):
        t, epsilon, h = FQ(1, 10000), FQ(1, 10), FQ(25)
        a = epsilon*t*t/12
        probe = a/(16*h*h)
        gamma = a*probe/8
        x = 2*h*probe
        remainder = x*x/(2*(1-x/3))
        delta = (remainder+gamma)/probe
        self.assertLess(delta, a/2)
        self.assertGreater(sum(FQ(7)**k/math.factorial(k) for k in range(25)), 800)
        copies = math.ceil(56/(a*a*probe*probe))
        self.assertGreaterEqual(FQ(copies), 56/(a*a*probe*probe))
        # ||r_hat-r|| <= 3/t^2 * 2sqrt(3)*a = sqrt(3)/2 epsilon.
        self.assertEqual((6*a/(t*t))**2*3, 3*epsilon*epsilon/4)
        # Finite computation of A: third-order Heisenberg polynomial already
        # meets this example's calibration tolerance; general epsilon permits
        # increasing the finite degree. This does not generate a physical clock.
        calibration = epsilon*t*t/100
        natural_x = 18*t
        calibration_tail = 3*natural_x**4/(24*(1-natural_x/5))
        self.assertLess(calibration_tail, calibration)
        robust_factor = (FQ(7, 24)+FQ(1, 100))/(FQ(1, 3)-FQ(1, 100))
        self.assertEqual(robust_factor, FQ(181, 194))
        self.assertLess(robust_factor, 1)
        gamma_ctrl, control_duration = gamma/2, gamma/(4*h)
        self.assertEqual(gamma_ctrl+2*h*control_duration, gamma)
        OBS['finite_readout_resources'] = dict(
            natural_wait=str(t), desired_r_Euclidean_error=str(epsilon),
            per_edge_total_error=str(a), probe_wait=str(probe),
            complete_instrument_error_budget=str(gamma),
            total_control_duration_budget=str(control_duration),
            control_only_error_budget=str(gamma_ctrl),
            exact_bias_upper_bound=str(delta), copies_per_edge=str(copies),
            total_independent_source_trials=str(4*copies),
            joint_failure_probability_at_most='1/100',
            exact_A_reconstructed_r_error_bound='sqrt(3)*epsilon/2 < epsilon',
            finite_A_calibration_error_budget=str(calibration),
            third_order_A_polynomial_remainder_bound=str(calibration_tail),
            finite_A_robust_r_error_bound='181*epsilon/194 < epsilon',
            finite_model_computation_not_autonomous_physical_calibration=True,
            whole_handoff_plus_probe_error_contract_on_arbitrary_old_data_graph_reference=True,
            initial_data_not_reset_or_cloned=True,
            no_huge_sampling_experiment_actually_run=True,
            clocks_fresh_probes_internal_storage_and_repeated_source_remain_inputs=True)

    def test_06_geometric_nonclosure_and_reference_rank_boundary(self):
        trees, _, _, population, _, *_ = system()
        source_means = [sum(FQ(2-int(n), 6) for n in row) for row in population]
        self.assertEqual(source_means, [FQ(3, 2)]*4)
        leaf_means = [sum(FQ(len(old.path_between(tree, 6, a, b))-1, 6)
                          for tree in trees) for a, b in itertools.combinations(LEAVES, 2)]
        self.assertEqual(leaf_means, [FQ(8, 3)]*6)
        circumradius_squared = FQ(3, 8)*FQ(8, 3)**2
        gap = circumradius_squared-FQ(3, 2)**2
        self.assertEqual(gap, FQ(5, 12))
        references = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        rank = int(np.linalg.matrix_rank(references[1:]-references[0]))
        self.assertEqual(rank, 3)
        collinear = np.array([[0, 0, 0], [0, 0, 1], [0, 0, -.5], [0, 0, .3]])
        self.assertEqual(int(np.linalg.matrix_rank(collinear[1:]-collinear[0])), 1)
        OBS['geometry_and_scope'] = dict(
            initial_source_to_four_leaves=[str(v) for v in source_means],
            initial_all_leaf_pair_means=[str(v) for v in leaf_means],
            squared_circumradius=str(circumradius_squared),
            required_squared_height=str(-gap),
            these_uncalibrated_means_not_lengths_in_any_flat_Euclidean_space=True,
            curved_space_and_independently_justified_calibration_not_excluded=True,
            selected_reference_affine_rank=rank,
            physical_anchors_not_proved_rigid=True,
            rank_three_response_is_only_chosen_preparation_family=True,
            all_location_directions_and_endpoint_identity_not_characterized=True,
            composition_of_displacements_and_cross_basepoint_atlas_not_proved=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=474, baseline_round=472, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(three_internal_relational_components_change_actual_graph_distance_means=True,
            same_fixed_H_during_natural_wait=True,
            internal_reference_preparation_and_readout_hardware_still_inputs=True,
            finite_six_tree_theorem_not_arbitrary_network_theorem=True,
            complete_unknown_quantum_input_not_claimed_recoverable_after_twirl=True,
            spatial_dimension_three_derived=False, full_GR_goal_completed=False,
            phase_closure_triggered=False))


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
