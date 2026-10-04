"""Round 463: exact participation supports from fixed finite local descriptions.

Baseline 461, independent of round 462. The product-description hypothesis is
explicit. Positive semidefinite factorization bounds exact nonzero supports;
thresholds and globally correlated records are separate, legal countermodels
to extending that hypothesis beyond its stated domain.
"""
import argparse
from fractions import Fraction
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

import exchange_relation_audit as old

TARGET = Path(__file__).with_name('local_descriptor_support_audit_results.json')
OBS = {}


def short(value):
    return float(f'{float(value):.12g}')


def exact_rank(matrix):
    a = [[Fraction(int(x)) for x in row] for row in matrix]
    rank = 0
    for column in range(len(a[0])):
        pivot = next((row for row in range(rank, len(a)) if a[row][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][column]
        a[rank] = [x/scale for x in a[rank]]
        for row in range(len(a)):
            if row != rank and a[row][column]:
                scale = a[row][column]
                a[row] = [x-scale*y for x, y in zip(a[row], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def pair_probabilities(states, effect):
    return np.array([[np.trace(effect@np.kron(a, b)).real for b in states] for a in states])


def induced_partner_effect(effect, rho):
    d = rho.shape[0]
    return np.einsum('iajb,ba->ij', effect.reshape(d, d, d, d), rho)


def state_support_cover(states):
    summed = np.zeros_like(states[0])
    selected, rank = [], 0
    for i, rho in enumerate(states):
        candidate = summed+rho
        next_rank = int(np.linalg.matrix_rank(candidate, tol=1e-11))
        if next_rank > rank:
            summed = candidate
            selected.append(i)
            rank = next_rank
    return selected, rank


def support_graph(matrix, threshold=1e-10):
    adjacency = matrix > threshold
    np.fill_diagonal(adjacency, False)
    return adjacency


def connected(adjacency):
    visited, stack = {0}, [0]
    while stack:
        v = stack.pop()
        for other in np.flatnonzero(adjacency[v]):
            if int(other) not in visited:
                visited.add(int(other))
                stack.append(int(other))
    return len(visited) == len(adjacency)


def equality_effect(d):
    return np.diag([int(a == b) for a in range(d) for b in range(d)]).astype(complex)


class Audit(unittest.TestCase):
    def close(self, left, right, tolerance=6e-12):
        self.assertLess(float(np.linalg.norm(left-right)), tolerance)

    def test_01_quantum_bilinear_factorization_can_have_rank_d_squared(self):
        directions = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])/math.sqrt(3)
        states = [(np.eye(2)+sum(r[mu]*old.PAULI[mu] for mu in range(3)))/2 for r in directions]
        effect = (np.eye(4)+old.swap(2, 0, 1))/2
        basis = [np.eye(2)/math.sqrt(2)]+[p/math.sqrt(2) for p in old.PAULI]
        coordinates = np.array([[np.trace(b@rho).real for b in basis] for rho in states])
        kernel = np.array([[np.trace(effect@np.kron(a, b)).real for b in basis] for a in basis])
        matrix = pair_probabilities(states, effect)
        self.close(matrix, coordinates@kernel@coordinates.T)
        certificate = np.eye(4, dtype=int)+2*np.ones((4, 4), dtype=int)
        self.close(3*matrix, certificate)
        self.assertEqual(exact_rank(certificate), 4)
        for j, rho in enumerate(states):
            f = induced_partner_effect(effect, rho)
            self.assertGreater(np.linalg.eigvalsh(f).min(), -1e-12)
            self.close(matrix[:, j], np.array([np.trace(a@f).real for a in states]))
        OBS['bilinear_and_positive_factorizations'] = dict(local_dimension=2,
            nonorthogonal_continuous_quantum_states_used=True,
            linear_descriptor_span_dimension=4, exact_probability_matrix_rank=4,
            three_times_matrix=certificate.tolist(),
            diagonal_defined_by_same_product_formula=True,
            unknown_state_cloning_operation_not_assumed=True,
            rank_bound='rank M <= dim span_R{rho_i} <= d^2',
            positive_factorization='M_ij = Tr(rho_i F_j), F_j >= 0')

    def test_02_positive_state_support_cover_works_without_orthogonal_labels(self):
        vectors = [np.array([1, 0, 0]), np.array([1, 1, 0])/math.sqrt(2),
                   np.array([1, 1j, 0])/math.sqrt(2)]+[np.array([0, 0, 1])]*3
        states = [np.outer(v, v.conj()) for v in vectors]
        pa, pb = np.diag([1, 1, 0]), np.diag([0, 0, 1])
        effect = np.kron(pa, pa)+np.kron(pb, pb)
        matrix = pair_probabilities(states, effect)
        expected = np.zeros((6, 6))
        expected[:3, :3] = 1
        expected[3:, 3:] = 1
        self.close(matrix, expected)
        selected, dimension = state_support_cover(states)
        self.assertEqual(selected, [0, 1, 3])
        self.assertEqual(dimension, 3)
        covered = np.any(matrix[selected] > .5, axis=0)
        self.assertTrue(covered.all())
        self.assertLessEqual(len(selected), dimension)
        self.assertEqual(int(support_graph(matrix).sum(axis=1).max()), 2)
        self.assertLessEqual(len(states), dimension*(2+1))
        OBS['positive_cover'] = dict(local_dimension=3, local_support_span_dimension=3,
            vertices=6, maximum_degree=2, selected_support_spanning_rows=selected,
            every_nonzero_column_covered=True,
            selected_rows_not_assumed_linear_basis_or_orthogonal_states=True,
            general_bound='N_nonisolated <= s*(Delta+1) <= d*(Delta+1)',
            rank_only_weaker_bound='N_nonisolated <= r*(Delta+1)')

    def test_03_exact_bound_attainment_and_dimension_growth_boundary(self):
        d, group_size = 3, 4
        effect = equality_effect(d)
        self.close(effect@effect, effect)
        states = [np.diag([int(k == a) for k in range(d)])
                  for a in range(d) for _ in range(group_size)]
        matrix = pair_probabilities(states, effect)
        integer = np.kron(np.eye(d, dtype=int), np.ones((group_size, group_size), dtype=int))
        self.close(matrix, integer, 1e-15)
        adjacency = support_graph(matrix)
        degree = int(adjacency.sum(axis=1).max())
        self.assertEqual(degree, 3)
        self.assertEqual(len(states), d*(degree+1))
        self.assertEqual(exact_rank(integer), d)
        self.assertFalse(connected(adjacency))
        # Dimension growth can encode a prescribed connected cycle; the
        # cycle is now input in the effect table, not generated geometry.
        n = 7
        cycle = np.array([[int((i-j) % n in (1, n-1)) for j in range(n)] for i in range(n)])
        cycle_effect = np.diag(cycle.ravel())
        labels = [np.diag([int(i == k) for i in range(n)]) for k in range(n)]
        self.close(pair_probabilities(labels, cycle_effect), cycle, 1e-15)
        self.assertTrue(connected(cycle.astype(bool)))
        self.assertEqual(exact_rank(cycle), 7)
        # Isolated nodes must be excluded from the counted N.
        active = np.array([1, 1]+[0]*18)
        graph_with_isolates = np.outer(active, active)
        np.fill_diagonal(graph_with_isolates, 0)
        self.assertEqual(np.count_nonzero(graph_with_isolates.sum(axis=1)), 2)
        isolate_states = [np.diag([a, 1-a]) for a in active]
        isolate_effect = np.diag([1, 0, 0, 0])
        self.close(support_graph(pair_probabilities(isolate_states, isolate_effect)).astype(int),
                   graph_with_isolates, 1e-15)
        OBS['sharpness_and_boundary'] = dict(
            attaining_example=dict(d=3, vertices=12, degree=3, components=3,
                                   graph='three disjoint K4', exact_bound_saturated=True),
            connected_sharpness_claimed=False,
            dimension_growth_example=dict(d=7, vertices=7, degree=2,
                                          graph='C7', graph_is_input_in_effect=True),
            arbitrary_isolated_nodes_are_not_bounded=True,
            no_programming_theorem_used=False)

    def test_04_small_full_rank_noise_destroys_exact_zeros(self):
        d, epsilon = 3, Fraction(1, 10)
        # Exact diagonal-state numerators: rho_a = diag(28,1,1)/30.
        numerators = [[28 if a == k else 1 for k in range(d)] for a in range(d)]
        probabilities = [[sum(Fraction(x*y, 900) for x, y in zip(a, b))
                          for b in numerators] for a in numerators]
        self.assertEqual(probabilities[0][0], Fraction(131, 150))
        self.assertEqual(probabilities[0][1], Fraction(19, 300))
        lower = epsilon**2*Fraction(d, d*d)
        self.assertEqual(lower, Fraction(1, 300))
        self.assertGreater(min(x for row in probabilities for x in row), lower)
        effect = equality_effect(d)
        states = [np.diag(np.array(row)/30) for row in numerators for _ in range(4)]
        matrix = pair_probabilities(states, effect)
        self.assertEqual(int(support_graph(matrix).sum(axis=1).min()), 11)
        OBS['full_rank_noise'] = dict(epsilon='1/10', d=3,
            same_label_probability='131/150', different_label_probability='19/300',
            general_positive_floor='epsilon^2 Tr(E)/d^2', example_floor='1/300',
            formerly_three_K4_now_exact_K12=True,
            exact_zero_sparsity_not_robust_under_arbitrarily_small_depolarization=True,
            threshold_graph_not_claimed_destroyed=True)

    def test_05_fixed_qubit_effect_and_fixed_threshold_allow_arbitrary_cycles(self):
        effect = (np.eye(4)+old.swap(2, 0, 1))/2
        threshold = .9
        rows = []
        for n in (12, 24):
            c1, c2 = math.cos(2*math.pi/n), math.cos(4*math.pi/n)
            eta = (c1+c2)/2
            radius_squared = (3/5)/eta
            self.assertGreater(eta, 3/5)
            self.assertLess(radius_squared, 1)
            radius = math.sqrt(radius_squared)
            states = [(np.eye(2)+radius*(math.cos(2*math.pi*k/n)*old.PAULI[0]
                                       +math.sin(2*math.pi*k/n)*old.PAULI[1]))/2 for k in range(n)]
            matrix = pair_probabilities(states, effect)
            direct = np.array([[(3+radius_squared*math.cos(2*math.pi*(i-j)/n))/4
                               for j in range(n)] for i in range(n)])
            self.close(matrix, direct, 2e-11)
            adjacency = support_graph(matrix, threshold)
            cycle = np.array([[(i-j) % n in (1, n-1) for j in range(n)] for i in range(n)])
            self.assertTrue(np.array_equal(adjacency, cycle))
            self.assertTrue(connected(adjacency))
            self.assertEqual(int(adjacency.sum(axis=1).max()), 2)
            self.assertEqual(int(support_graph(matrix).sum(axis=1).min()), n-1)
            margin = radius_squared*(c1-c2)/8
            self.assertAlmostEqual(direct[0, 1]-threshold, margin)
            self.assertAlmostEqual(threshold-direct[0, 2], margin)
            self.assertLessEqual(margin, 3*math.pi**2/(4*n*n))
            self.assertEqual(int(np.linalg.matrix_rank(matrix, tol=1e-9)), 3)
            self.assertEqual(exact_rank(adjacency.astype(int)), n-2)
            rows.append(dict(vertices=n, local_dimension=2, threshold='9/10',
                radius_squared=short(radius_squared), probability_margin=short(margin),
                exact_support_graph=f'K{n}', threshold_support_graph=f'C{n}',
                probability_matrix_rank=3, threshold_adjacency_rank=n-2))
        OBS['threshold_counterfamily'] = dict(fixed_effect='(I+SWAP)/2', fixed_dimension=2,
            fixed_threshold='9/10', analytic_range='all integers N >= 12',
            examples=rows, threshold_not_changed_with_N=True,
            local_preparation_radius_and_angles_depend_on_N=True,
            gap_asymptotic='delta_N ~ 9 pi^2/(20 N^2)',
            gap_upper='delta_N <= 3 pi^2/(4 N^2)',
            nonzero_support_theorem_not_violated=True,
            state_dependent_Hamiltonian_weights_not_implemented=True,
            spatial_dimension_not_inferred=True)

    def test_06_correlated_node_records_can_encode_sparse_mean_support(self):
        n = 8
        edges = [(i, (i+1) % n) for i in range(n)]
        configurations = np.zeros((n, n), dtype=int)
        for row, (i, j) in enumerate(edges):
            configurations[row, i] = configurations[row, j] = 1
        counts = configurations.T@configurations
        actual = counts/n
        expected = np.array([[2 if i == j else int((i-j) % n in (1, n-1))
                              for j in range(n)] for i in range(n)])/n
        self.close(actual, expected, 1e-15)
        adjacency = support_graph(actual)
        self.assertTrue(connected(adjacency))
        self.assertEqual(int(adjacency.sum(axis=1).max()), 2)
        self.assertGreater(n, 2*(2+1))
        local_one = configurations.sum(axis=0)/n
        self.close(local_one, np.ones(n)/4)
        product_prediction = np.outer(local_one, local_one)
        self.assertEqual(int(support_graph(product_prediction).sum(axis=1).min()), n-1)
        for config in configurations:
            branch = np.outer(config, config)
            np.fill_diagonal(branch, 0)
            self.assertEqual(int(branch.sum()/2), 1)
        effect = np.diag([0, 0, 0, 1])
        self.close(effect@effect, effect)
        OBS['correlated_support_boundary'] = dict(vertices=n, local_dimension=2,
            fixed_effect='|11><11|', global_state='uniform classical mixture of one active endpoint-pair per graph edge',
            off_diagonal_edge_probability='1/8', nonedge_probability='0',
            all_local_states_identical='diag(3/4,1/4)',
            product_description_prediction_for_every_distinct_pair='1/16',
            actual_mean_support='C8', independent_prediction_support='K8',
            each_configuration_active_graph='one K2 and isolated nodes',
            graph_is_stored_in_global_correlations_not_generated=True,
            arbitrary_graph_extension_analytic=True,
            mixture_support_not_same_as_simultaneous_active_configuration=True,
            support_graph_not_identified_with_causal_propagation=True)


def run():
    OBS.clear()
    output = io.StringIO()
    result = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=463, baseline_round=461, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(exact_bound='N_nonisolated <= d*(Delta+1)',
            nonorthogonal_continuous_descriptions_allowed=True,
            local_product_description_contract_is_an_additional_assumption=True,
            arbitrary_correlated_global_states_excluded_from_bound=True,
            finite_threshold_graphs_excluded_from_exact_zero_bound=True,
            fixed_E_measurement_hardware_or_natural_Hamiltonian_derived=False,
            independent_from_round_462=True,
            all_cognitive_networks_have_finite_size_claimed=False,
            spatial_dimension_derived=False, full_GR_goal_completed=False,
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
