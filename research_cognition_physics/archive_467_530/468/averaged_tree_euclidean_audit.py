"""Round 468: all-label Euclidean means of fixed-degree trees.

The classification is analytic. Exact integer/Fraction audits check its tree
identity, signed projection equality, and path-order boundary. A cross-degree
prepared stationary family under the old Hamiltonian is a necessary scope
counterexample, not an autonomous preparation or selection of three dimensions.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import branching_tree_distance_audit as old

TARGET = Path(__file__).with_name('averaged_tree_euclidean_audit_results.json')
OBS = {}


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for rest in compositions(total-first, length-1):
                yield (first,)+rest


def weighted_distances(n, edges):
    adj = [[] for _ in range(n)]
    for u, v, weight in edges:
        assert weight > 0
        adj[u].append((v, Q(weight)))
        adj[v].append((u, Q(weight)))
    assert len(edges) == n-1
    d = np.full((n, n), Q(0), dtype=object)
    for source in range(n):
        known = {source: Q(0)}
        todo = [source]
        while todo:
            u = todo.pop()
            for v, w in adj[u]:
                if v not in known:
                    known[v] = known[u]+w
                    todo.append(v)
        assert len(known) == n
        for v, value in known.items():
            d[source, v] = value
    return d, np.array([2-len(a) for a in adj], dtype=np.int64)


def ordinary_distances(tree, n):
    return np.array([[len(old.path_between(tree, n, a, b))-1
                      for b in range(n)] for a in range(n)], dtype=np.int64)


def all_paths(n):
    return {frozenset(old.edge(a, b) for a, b in zip(order, order[1:]))
            for order in itertools.permutations(range(n))}


def path_metric(order, weights):
    return weighted_distances(len(order), [(a, b, w) for a, b, w in
                               zip(order, order[1:], weights)])[0]


def average(ds, weights):
    return sum((w*d for w, d in zip(weights, ds)), np.full(ds[0].shape, Q(0), dtype=object))


def binary_six_trees():
    result = []
    for left in itertools.combinations(range(4), 2):
        tree = {old.edge(4, 5)}
        tree.update(old.edge(i, 4 if i in left else 5) for i in range(4))
        result.append((left, frozenset(tree)))
    return result


class Audit(unittest.TestCase):
    def test_01_weighted_tree_and_mean_degree_identities(self):
        checked, grouped, mixed = 0, 0, 0
        for n in range(2, 6):
            groups = defaultdict(list)
            for word in itertools.product(range(n), repeat=n-2):
                tree = old.prufer_tree(n, word)
                edges = [(a, b, Q(k+2, n+1)) for k, (a, b) in enumerate(sorted(tree))]
                d, tau = weighted_distances(n, edges)
                length = sum((w for _, _, w in edges), Q(0))
                self.assertEqual(int(tau.sum()), 2)
                self.assertTrue(np.array_equal(d@tau, np.full(n, length, dtype=object)))
                groups[tuple(tau)].append((d, length))
                checked += 1
            for tup, data in groups.items():
                tau = np.array(tup, dtype=np.int64)
                first, last = data[0], data[-1]
                m = (first[0]+last[0])/2
                c = (first[1]+last[1])/2
                self.assertTrue(np.array_equal(m@tau, np.full(n, c, dtype=object)))
                for i in range(n):
                    b = tau.copy(); b[i] -= 1
                    self.assertEqual(int(b.sum()), 1)
                    self.assertEqual(b@m@b, 0)
                grouped += 1
                mixed += len(data) > 1
        OBS['degree_identity'] = dict(weighted_trees_checked=checked,
            degree_vectors_checked=grouped, nontrivial_mixture_classes=mixed,
            exact_formula='D*tau=L*1, tau_i=2-deg_i, sum(tau)=2',
            mean_formula='M*tau=E[L]*1',
            all_signed_integer_quadratic_forms_zero=True,
            weighted_identity_and_classification_are_analytically_proved=True)

    def test_02_integer_projection_equality_certificate(self):
        vectors, forms, surviving = 0, 0, 0
        for n in range(2, 8):
            positions = np.array([i*(i+1)//2 for i in range(n)], dtype=np.int64)
            line = abs(positions[:, None]-positions[None, :])
            gaps = np.diff(positions)
            for counts in compositions(n-2, n):
                tau = 1-np.array(counts, dtype=np.int64)
                energies = []
                for i in range(n):
                    b = tau.copy(); b[i] -= 1
                    prefix = np.cumsum(b)[:-1]
                    direct = int(b@line@b)
                    decomposed = int(2*np.sum(gaps*prefix*(1-prefix)))
                    self.assertEqual(direct, decomposed)
                    self.assertLessEqual(direct, 0)
                    energies.append(direct)
                    forms += 1
                expected = np.zeros(n, dtype=np.int64); expected[0] = expected[-1] = 1
                equality = all(e == 0 for e in energies)
                self.assertEqual(equality, bool(np.array_equal(tau, expected)))
                surviving += equality
                vectors += 1
        OBS['integer_projection'] = dict(tree_degree_vectors_checked=vectors,
            exact_signed_forms_checked=forms,
            equality_degree_vector_count=surviving,
            one_endpoints_only_vector_per_N=True,
            identity='b^T D_line b=2 sum_k gap_k B_k(1-B_k)',
            integration_over_directions_used_in_proof_not_numerical_quadrature=True,
            finite_certificates_not_substituted_for_general_proof=True)

    def test_03_path_order_is_necessary_and_sufficient(self):
        order = (0, 3, 1, 4, 2, 5)
        lengths = [[Q(i+1, 3) for i in range(5)],
                   [Q(7-i, 5) for i in range(5)],
                   [Q((i+1)**2+1, 7) for i in range(5)]]
        ds = [path_metric(order, w) for w in lengths]
        ws = [Q(1, 6), Q(1, 3), Q(1, 2)]
        m = average(ds, ws)
        positions = m[0, :]
        line = abs(positions[:, None]-positions[None, :])
        self.assertTrue(np.array_equal(m, line))
        changed = [path_metric((0, 1, 2, 3), [Q(1)]*3),
                   path_metric((0, 2, 1, 3), [Q(1)]*3)]
        other = average(changed, [Q(1, 3), Q(2, 3)])
        gap = other[1, 2]-abs(other[0, 1]-other[0, 2])
        self.assertEqual(gap, Q(2, 3))
        self.assertTrue(all(other[0, i]+other[i, 3] == other[0, 3] for i in range(4)))
        near = []
        for eps in (Q(1, 10), Q(1, 100), Q(1, 1000)):
            star, _ = weighted_distances(4, [(0, 1, Q(1)), (0, 2, Q(1)), (0, 3, eps)])
            points = np.array([[0., 0.], [-1., 0.], [1., 0.], [0., float(eps)]])
            euclidean = np.linalg.norm(points[:, None, :]-points[None, :, :], axis=-1)
            off = ~np.eye(4, dtype=bool)
            additive = float(np.max(abs(star.astype(float)-euclidean)))
            relative = float(np.max(abs(star.astype(float)[off]-euclidean[off])/euclidean[off]))
            # Exact inequalities certify 0<1+eps-sqrt(1+eps^2)<eps.
            self.assertGreater((1+eps)**2, 1+eps**2)
            self.assertGreater(1+eps**2, 1)
            self.assertLess(additive, float(eps))
            self.assertLess(relative, float(eps))
            near.append(dict(epsilon=str(eps), max_additive_error=float(format(additive, '.12g')),
                max_pairwise_relative_error=float(format(relative, '.12g'))))
        OBS['path_order_boundary'] = dict(
            same_order_variable_positive_lengths_give_exact_line=True,
            changed_order_with_fixed_endpoints_Jensen_gap=str(gap),
            endpoints_triangle_equalities_hold=True,
            unit_edge_case_requires_one_path_with_probability_one=True,
            arbitrary_continuous_weight_distributions_covered_by_proof=True,
            positive_edge_star_near_line_examples=near,
            exact_obstruction_not_uniform_over_degenerating_edge_lengths=True,
            short_leaf_and_hub_coalesce_only_in_limit=True)

    def test_04_cross_degree_uniform_paths_and_simplex(self):
        cases = []
        for n in range(3, 8):
            trees = all_paths(n)
            total = np.zeros((n, n), dtype=np.int64)
            degrees = set()
            endpoints = np.zeros(n, dtype=np.int64)
            for tree in trees:
                changes = old.tree_flips(tree, n)
                self.assertEqual(sum(changes.values()), n-3)
                self.assertTrue(all(t in trees for t in changes))
                self.assertTrue(all(a == 1 for a in changes.values()))
                total += ordinary_distances(tree, n)
                deg = tuple(map(len, old.adjacency(tree, n)))
                degrees.add(deg)
                endpoints += np.array(deg) == 1
            count = math.factorial(n)//2
            off = np.ones((n, n), dtype=np.int64)-np.eye(n, dtype=np.int64)
            self.assertEqual(len(trees), count)
            self.assertTrue(np.array_equal(3*total, (n+1)*count*off))
            self.assertTrue(np.array_equal(n*endpoints, np.full(n, 2*count)))
            self.assertEqual(len(degrees), math.comb(n, 2))
            projector_n = n*np.eye(n, dtype=np.int64)-np.ones((n, n), dtype=np.int64)
            self.assertTrue(np.array_equal(projector_n@projector_n, n*projector_n))
            self.assertEqual(int(np.trace(projector_n)), n*(n-1))
            cases.append(dict(N=n, paths=count, distinct_degree_vectors=len(degrees),
                F_row_sum=n-3, mean_pair_distance=str(Q(n+1, 3)),
                required_Euclidean_dimension=n-1))
        OBS['cross_degree_stationary_family'] = dict(cases=cases,
            uniform_path_state_eigenvalue_of_F='N-3 for N>=3',
            simplex_rank_certificate='(N I-11^T)^2=N(N I-11^T), trace=N(N-1)',
            chosen_initial_distribution_is_input=True,
            dimension_three_selected=False,
            state_compatible_with_existing_degree_preserving_H=True)

    def test_05_original_full_H_preserves_symmetric_data_and_reference(self):
        n = 4
        graphs = sorted(old.mask(t, n) for t in all_paths(n))
        h, f = old.previous.full_model(n, graphs)
        self.assertTrue(np.array_equal(f@np.ones(len(graphs), dtype=np.int64),
                                       np.ones(len(graphs), dtype=np.int64)))
        dicke = np.zeros((2**n, n+1), dtype=np.int64)
        for state in range(2**n):
            dicke[state, state.bit_count()] = 1
        embed = np.kron(dicke, np.ones((len(graphs), 1), dtype=np.int64))
        eigenvalue = 3*n-5  # J=kappa=mu=1, the old full_model convention.
        residual = h@embed-eigenvalue*embed
        self.assertEqual(int(np.max(abs(residual))), 0)
        norms = len(graphs)*np.array([math.comb(n, k) for k in range(n+1)])
        self.assertTrue(np.array_equal(embed.T@embed, np.diag(norms)))
        normalized = embed/np.sqrt(norms)[None, :]
        rng = np.random.default_rng(468)
        logical_reference = rng.normal(size=(n+1, 3))+1j*rng.normal(size=(n+1, 3))
        logical_reference /= np.linalg.norm(logical_reference)
        joint = normalized@logical_reference
        self.assertLess(np.linalg.norm(h@joint-eigenvalue*joint), 2e-14)
        OBS['unchanged_full_H'] = dict(N=n, full_dimension=len(h),
            graph_states=len(graphs), symmetric_data_dimension=n+1,
            exact_integer_intertwining_residual=0,
            parameters=dict(mu=1, J=1, kappa=1), eigenvalue=eigenvalue,
            general_eigenvalue='kappa*(N-3)+(mu+J)*(N-1)',
            arbitrary_symmetric_data_and_internal_reference_preserved=True,
            all_unknown_raw_data_claimed=False,
            SU2_invariant_preparation_option='P_sym/(N+1)',
            preparation_protocol_derived=False)

    def test_06_six_point_squared_distance_certificate(self):
        trees = binary_six_trees()
        ds = [ordinary_distances(t, 6).astype(object) for _, t in trees]
        x = np.array([1, 1, 1, 1, -2, -2], dtype=np.int64)
        self.assertEqual(int(x.sum()), 0)
        values = []
        groups = [(j, 5-j) for j in range(3)]
        for counts in compositions(5, 6):
            ws = [Q(k, 5) for k in counts]
            m = average(ds, ws)
            q = [ws[a]+ws[b] for a, b in groups]
            u = [sum((w for w, (left, _) in zip(ws, trees) if i in left), Q(0))
                 for i in range(4)]
            self.assertEqual(sum(q), 1)
            self.assertEqual(sum(u), 2)
            self.assertTrue(all(0 <= a <= 1 for a in u))
            value = x@(m*m)@x
            expression = 28+4*sum(a*a for a in q)-8*sum(a*a for a in u)
            self.assertEqual(value, expression)
            self.assertGreaterEqual(value, Q(40, 3))
            values.append(value)
        uniform = average(ds, [Q(1, 6)]*6)
        self.assertEqual(x@(uniform*uniform)@x, Q(64, 3))
        OBS['all_six_vertices_diagnostic'] = dict(
            degree_vector=[1, 1, 1, 1, 3, 3], trees=6,
            rational_populations_checked=len(values),
            exact_expression='x^T(M o M)x=28+4 sum q_j^2-8 sum u_i^2',
            universal_conservative_lower_bound='40/3',
            lower_bound_claimed_sharp=False,
            uniform_population_value='64/3',
            grid_minimum_only=str(min(values)),
            nonpositive_squared_distance_form_required_by_Euclidean_lengths=True,
            retained_hubs_expose_leaf_only_interface_omission=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=468, baseline_round=467, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(all_labels_are_retained=True,
            theorem_requires_fixed_degree_of_each_label=True,
            edge_lengths_strictly_positive_almost_surely=True,
            expected_total_tree_length_finite=True,
            ordinary_Euclidean_lengths_not_squared_lengths_classified=True,
            any_finite_Euclidean_dimension_allowed=True,
            arbitrary_joint_graph_data_reference_state_extension=True,
            round465_allows_cross_degree_superposition=True,
            cross_degree_example_does_not_disprove_fixed_degree_theorem=True,
            uniform_in_N_approximation_obstruction_proved=False,
            quotiented_or_different_measurement_geometries_excluded=False,
            stationary_prepared_mean_not_autonomously_selected=True,
            actual_internal_distance_reader_derived=False,
            dimension_three_selected=False, spatial_dimension_derived=False,
            full_GR_goal_completed=False, phase_closure_triggered=False,
            historical_files_changed=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--dry-run', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert result == json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as output:
            output.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
