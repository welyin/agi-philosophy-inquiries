"""Round 478: complete current distance-mean quotient at fixed leaf statistics.

Scientific baseline 475. The selected six-tree department and readout family
are inputs. The quotient is a statistical three-cube, not a derived physical
space, a single-shot coordinate reading, or a full quantum-state encoding.
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
import rigid_leaf_reference_audit as previous
import sequential_tree_distance_readout_audit as reader

TARGET = Path(__file__).with_name('stable_reference_coordinate_quotient_results.json')
OBS = {}
LEAVES = (0, 3, 4, 5)
PAIR_ORDER = tuple(itertools.combinations(range(6), 2))


def short(x):
    return float(f'{float(x):.11g}')


@lru_cache(None)
def system():
    trees, h, f = previous.system()
    distances = np.array([previous.distance(a, b) for a, b in PAIR_ORDER], dtype=np.int64)
    plus, minus = [], []
    for b in LEAVES[1:]:
        first = frozenset((0, b))
        second = frozenset(LEAVES)-first
        for target, labels in ((plus, first), (minus, second)):
            target.append(next(i for i, tree in enumerate(trees)
                if frozenset(j for j in LEAVES if tuple(sorted((1, j))) in tree) == labels))
    return h, f, distances, tuple(plus), tuple(minus)


def population(x):
    result = [Q(0)]*6
    for b, (plus, minus) in enumerate(zip(system()[3], system()[4])):
        result[plus] = Q(1, 6)+x[b]/2
        result[minus] = Q(1, 6)-x[b]/2
    return result


def x_to_q(x):
    return [sum(x)-value for value in x]


def q_to_x(q):
    return [sum(q)/2-value for value in q]


def all_means_from_q(q):
    leaf_to_1 = {0: Q(3, 2)-sum(q)/4}
    leaf_to_1.update({b: leaf_to_1[0]+q[k] for k, b in enumerate(LEAVES[1:])})
    result = []
    for a, b in PAIR_ORDER:
        if a in LEAVES and b in LEAVES:
            value = Q(8, 3)
        elif (a, b) == (1, 2):
            value = Q(1)
        else:
            internal = a if a in (1, 2) else b
            leaf = b if a in (1, 2) else a
            value = leaf_to_1[leaf] if internal == 1 else 3-leaf_to_1[leaf]
        result.append(value)
    return result


def rank_fraction(matrix):
    a = [[Q(int(x)) for x in row] for row in matrix]
    row = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(row, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        divisor = a[row][col]
        a[row] = [x/divisor for x in a[row]]
        for j in range(len(a)):
            if j != row:
                scale = a[j][col]
                a[j] = [x-scale*y for x, y in zip(a[j], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def exact_means(p):
    return [sum(Q(int(d))*v for d, v in zip(row, p)) for row in system()[2]]


def central_graph_numerator(a):
    return sum(character*a[np.ix_(gm, gm)]
               for _, _, gm, character in previous.permutations())


class Audit(unittest.TestCase):
    def close(self, a, b, tol=3e-11):
        self.assertLess(float(np.linalg.norm(np.asarray(a)-np.asarray(b))), tol)

    def test_01_complete_mean_map_and_exact_inverse(self):
        _, _, distances, plus, minus = system()
        self.assertEqual(plus, (5, 4, 2))
        self.assertEqual(minus, (0, 1, 3))
        self.assertEqual(rank_fraction(distances), 6)
        matrix = np.ones((3, 3), dtype=np.int64)-np.eye(3, dtype=np.int64)
        self.assertTrue(np.array_equal(matrix@(np.ones((3, 3), dtype=np.int64)
            -2*np.eye(3, dtype=np.int64)), 2*np.eye(3, dtype=np.int64)))
        points = list(itertools.product((Q(-1, 3), Q(0), Q(1, 3)), repeat=3))
        points += [(Q(1, 6), -Q(1, 12), Q(1, 4)), (Q(1, 17), -Q(1, 19), Q(1, 23))]
        for x in points:
            p = population(x)
            q = x_to_q(x)
            self.assertEqual(q_to_x(q), list(x))
            self.assertEqual(exact_means(p), all_means_from_q(q))
            self.assertEqual(sum(p), 1)
            self.assertTrue(all(v >= 0 for v in p))
        OBS['complete_current_mean_classification'] = dict(
            plus_graph_indices=list(plus), minus_graph_indices=list(minus),
            all_15_distance_matrix_exact_rank=6,
            independent_fixed_leaf_class_constraints=2,
            quotient_dimension=3, Q_from_x_matrix=matrix.tolist(), determinant=2,
            inverse='x_b=(sum(Q)-2 Q_b)/2',
            population='p_b_plus/minus=1/6 +/- (sum(Q)-2 Q_b)/4',
            exact_rational_points_checked=len(points),
            all_15_means_determined_by_Q_on_fixed_leaf_distribution=True,
            current_mean_equivalence_iff_graph_populations_equal=True)

    def test_02_necessary_sufficient_domain_and_boundary(self):
        vertices = []
        for signs in itertools.product((-1, 1), repeat=3):
            x = [Q(s, 3) for s in signs]
            p = population(x)
            self.assertEqual(sorted(p), [Q(0)]*3+[Q(1, 3)]*3)
            q = x_to_q(x)
            self.assertTrue(all(abs(sum(q)-2*v) == Q(2, 3) for v in q))
            vertices.append([str(v) for v in q])
            # The same vertex populations also admit a coherent pure graph
            # state; quotient coordinates do not determine original purity.
            psi = np.sqrt(np.asarray(p, dtype=float))
            pure = np.outer(psi, psi)
            self.close(np.diag(pure), np.asarray(p, dtype=float))
            self.close(np.trace(pure@pure), 1.)
        ranks = []
        for saturated in range(4):
            x = [Q(1, 3) if k < saturated else Q(0) for k in range(3)]
            rank = sum(v > 0 for v in population(x))
            self.assertEqual(rank, 6-saturated)
            ranks.append(dict(saturated_coordinate_count=saturated,
                              positive_graph_population_count=rank,
                              section_joint_state_rank=64*rank))
        invalid_x = [Q(1, 3)+Q(1, 100), Q(0), Q(0)]
        invalid_q = x_to_q(invalid_x)
        self.assertLess(min(population(invalid_x)), 0)
        self.assertTrue(any(abs(sum(invalid_q)-2*v) > Q(2, 3) for v in invalid_q))
        OBS['exact_domain'] = dict(
            x_domain='[-1/3,1/3]^3',
            Q_inequalities='abs(sum(Q)-2 Q_b) <= 2/3 for b=3,4,5',
            Q_vertices=vertices, vertices=8, edges=12, faces=6,
            representative_ranks=ranks,
            canonical_diagonal_vertex_representatives_are_mixed=True,
            pure_graph_states_with_same_vertex_populations_checked=True,
            current_coordinates_do_not_determine_original_quantum_purity=True,
            necessity_from_positive_populations_sufficiency_from_explicit_state=True,
            topology_is_the_finite_current_mean_readout_topology=True)

    def test_03_every_coordinate_has_protected_positive_representative(self):
        _, _, _, plus, minus = system()
        raw_basis = [np.eye(6, dtype=np.int64)]
        for p, m in zip(plus, minus):
            a = np.zeros((6, 6), dtype=np.int64)
            a[p, p], a[m, m] = 1, -1
            raw_basis.append(a)
        for a in raw_basis:
            self.assertTrue(np.array_equal(central_graph_numerator(a), np.zeros((6, 6))))
        h = system()[0]
        eigenvalues, eigenvectors = np.linalg.eigh(h)
        t = .37
        u = (eigenvectors*np.exp(-1j*t*eigenvalues))@eigenvectors.conj().T
        rows = []
        for x in ((Q(1, 6), -Q(1, 12), Q(1, 4)), (Q(1, 3),)*3):
            p = np.array(list(map(float, population(x))))
            rho = np.kron(np.eye(64)/64, np.diag(p))
            result = u@rho@u.conj().T
            graph = result.reshape(64, 6, 64, 6).trace(axis1=0, axis2=2)
            pnew = np.diag(graph).real
            for a, b in zip(plus, minus):
                self.close(pnew[a]+pnew[b], 1/3)
            rows.append(dict(initial_x=[str(v) for v in x],
                evolved_x=[short(pnew[a]-pnew[b]) for a, b in zip(plus, minus)],
                largest_leaf_class_error=short(max(abs(pnew[a]+pnew[b]-1/3)
                                                   for a, b in zip(plus, minus)))))
        OBS['protected_section'] = dict(
            section='rho_x=I_D/64 tensor diag(p(x))',
            four_exact_graph_E22_basis_numerators_zero=True,
            data_identity_invariant_under_leaf_permutations=True,
            full_joint_E22_zero_for_every_legal_x=True,
            arbitrary_real_kappa_J_and_all_times_follow_from_475_commutation=True,
            protected_state_subset_surjects_onto_entire_current_mean_cube=True,
            full_U_examples=rows,
            section_is_not_a_lossless_channel_resetting_arbitrary_unknown_input=True,
            dynamic_reachability_from_475_eta_not_proved=True)

    def test_04_readout_metric_and_quantum_information_boundary(self):
        points = list(itertools.product((Q(-1, 3), Q(0), Q(1, 3)), repeat=3))
        for x, y in itertools.product(points, repeat=2):
            px, py = population(x), population(y)
            dx, dy = exact_means(px), exact_means(py)
            operational = max(abs(a-b) for a, b in zip(dx, dy))
            tv = sum(abs(a-b) for a, b in zip(px, py))/2
            l1 = sum(abs(a-b) for a, b in zip(x, y))/2
            self.assertEqual(operational, tv)
            self.assertEqual(tv, l1)
        diagonal = np.eye(6)/6
        coherent = np.ones((6, 6))/6
        self.close(np.diag(diagonal), np.diag(coherent))
        self.close(central_graph_numerator(diagonal), 0)
        self.close(central_graph_numerator(coherent), 0)
        quantum_distance = float(np.abs(np.linalg.eigvalsh(diagonal-coherent)).sum()/2)
        self.assertAlmostEqual(quantum_distance, 5/6)
        # Ordinary coordinate numbers do not encode conditional R statistics,
        # even when each R branch has a protected diagonal representative.
        p_plus = np.array(population((Q(1, 3), Q(0), Q(0))), float)
        p_minus = np.array(population((-Q(1, 3), Q(0), Q(0))), float)
        rho_gr = (np.kron(np.diag(p_plus), np.diag([1., 0.]))
                  +np.kron(np.diag(p_minus), np.diag([0., 1.])))/2
        independent = np.kron(np.eye(6)/6, np.eye(2)/2)
        d10 = system()[2][PAIR_ORDER.index((0, 1))]
        moment = lambda a: np.einsum('g,grgs->rs', d10, a.reshape(6, 2, 6, 2))
        conditional_difference = moment(rho_gr)-moment(independent)
        self.close(conditional_difference, np.diag([-1/12, 1/12]))
        OBS['metric_and_quantum_scope'] = dict(
            exact_rational_pairs_checked=len(points)**2,
            operational_metric='max_15 abs(delta mean distance)=||delta x||_1/2',
            graph_population_total_variation_equals_operational_metric=True,
            canonical_diagonal_section_trace_distance_equals_operational_metric=True,
            general_joint_quantum_trace_distance_not_determined_by_current_means=True,
            same_center_means_graph_trace_distance_witness=short(quantum_distance),
            current_mean_identity_not_full_future_process_equivalence=True,
            scalar_coordinates_do_not_determine_R_conditioned_statistics=True,
            protected_classical_reference_D10_moment_difference=["-1/12", "1/12"],
            operational_distinguishability_norm_not_physical_spatial_length=True)

    def test_05_affine_global_symmetry_scope(self):
        vertices = list(itertools.product((-1, 1), repeat=3))
        maps = set()
        for permutation in itertools.permutations(range(3)):
            for signs in itertools.product((-1, 1), repeat=3):
                matrix = np.zeros((3, 3), dtype=int)
                for row, col in enumerate(permutation):
                    matrix[row, col] = signs[row]
                images = [tuple(matrix@np.array(v)) for v in vertices]
                self.assertEqual(set(images), set(vertices))
                maps.add(tuple(matrix.ravel()))
        self.assertEqual(len(maps), 48)
        # The enumeration verifies realizations; completeness uses the proof
        # that any affine automorphism fixes the center and permutes facets.
        OBS['affine_symmetry_only'] = dict(
            affine_global_automorphism_count=48,
            all_signed_permutations_verified=True,
            completeness_requires_center_and_opposite_facet_proof=True,
            no_nontrivial_continuous_SO3_global_affine_action=True,
            local_or_nonlinear_continuous_SO3_actions_not_excluded=True,
            affine_automorphisms_not_claimed_as_implemented_controls=True)

    def test_06_actual_four_edge_readout_and_finite_error_transfer(self):
        x = np.array([1/6, -1/12, 1/4])
        p = np.array(population(x), dtype=float)
        psi = np.sqrt(p)*np.exp(1j*np.arange(6)/3)
        rho = np.outer(psi, psi.conj())
        tau = Q(1, 65536)
        adjacency, estimated = [], []
        trees = previous.system()[0]
        for leaf in LEAVES:
            n = np.array([int(tuple(sorted((1, leaf))) in tree) for tree in trees])
            adjacency.append(float(n@p))
            minus, plus, _ = reader.instrument((1, leaf), tau)
            signed = reader.apply_map((plus-minus)/float(tau), rho, g=6, r=1)
            estimated.append(float(np.trace(signed).real))
        estimated, adjacency = np.asarray(estimated), np.asarray(adjacency)
        raw_x = np.array([(estimated[0]+estimated[k]
                          -sum(estimated[j] for j in range(1, 4) if j != k))/2
                          for k in range(1, 4)])
        clipped = np.clip(raw_x, -1/3, 1/3)
        edge_error = float(np.max(np.abs(estimated-adjacency)))
        self.assertLessEqual(float(np.max(np.abs(clipped-x))), 2*edge_error+1e-13)
        pop_error = float(np.abs(np.asarray(population(clipped), float)-p).sum()/2)
        self.assertLessEqual(pop_error, 3*edge_error+1e-13)
        epsilon, h = Q(1, 100), Q(25)
        a = epsilon/3
        v = a/(16*h*h)
        gamma = a*v/8
        y = 2*h*v
        remainder = y*y/(2*(1-y/3))
        bias = (remainder+gamma)/v
        self.assertLess(bias, a/2)
        self.assertGreater(sum(Q(7)**k/math.factorial(k) for k in range(25)), 800)
        copies = math.ceil(56/(a*a*v*v))
        self.assertEqual(3*a, epsilon)
        OBS['actual_readout_and_error_transfer'] = dict(
            diagnostic_probe_time=str(tau), coherent_graph_input=True,
            true_x=[short(v) for v in x], reconstructed_clipped_x=[short(v) for v in clipped],
            largest_actual_edge_bias=short(edge_error), actual_population_TV_error=short(pop_error),
            edge_error_a_implies_x_infinity_error_at_most='2a',
            component_clipping_never_increases_coordinate_error=True,
            all_15_mean_errors_and_population_TV_at_most='3a',
            target_total_mean_and_population_error=str(epsilon), per_edge_error=str(a),
            probe_wait=str(v), complete_instrument_diamond_budget=str(gamma),
            control_only_budget=str(gamma/2), total_control_duration_budget=str(gamma/(4*h)),
            copies_per_edge=str(copies), total_independent_trials=str(4*copies),
            joint_failure_probability_at_most='1/100',
            exact_fixed_leaf_distribution_is_a_promise_not_inferred_from_four_edges=True,
            old_data_reference_storage_and_first_handoff_in_full_instrument_contract=True,
            independent_sources_clocks_and_controls_remain_inputs=True,
            no_cloning_single_shot_graph_readout_or_autonomous_preparation_claim=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=478, baseline_round=475, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(complete_current_15_distance_mean_quotient_classified=True,
            selected_six_tree_department_and_fixed_leaf_statistics_are_inputs=True,
            every_legal_coordinate_has_an_E22_zero_positive_representative=True,
            representative_existence_not_dynamic_controllability=True,
            current_readout_topological_dimension_three_not_universe_dimension=True,
            no_requirement_that_position_predict_all_future_internal_information=True,
            displacement_composition_or_cross_basepoint_atlas_not_derived=True,
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
