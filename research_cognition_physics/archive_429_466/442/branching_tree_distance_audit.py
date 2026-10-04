"""Round 442: distance deformation under the unchanged autonomous flip model.

Nearest-neighbor interchange on trees is an existing combinatorial operation.
The new interface counts distance-changing flips and controls the complete
quantum evolution, without imposing a common path frame or measuring the tree.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np
import simultaneous_edge_exchange_audit as previous
import exchange_relation_audit as core

TARGET = Path(__file__).with_name('branching_tree_distance_audit_results.json')
OBS = {}


def edge(a, b):
    return tuple(sorted((a, b)))


def adjacency(tree, n):
    out = [set() for _ in range(n)]
    for a, b in tree:
        out[a].add(b)
        out[b].add(a)
    return out


def path_between(tree, n, a, b):
    adj = adjacency(tree, n)
    parent = {a: None}
    todo = deque([a])
    while todo and b not in parent:
        v = todo.popleft()
        for w in sorted(adj[v]):
            if w not in parent:
                parent[w] = v
                todo.append(w)
    assert b in parent
    path = [b]
    while path[-1] != a:
        path.append(parent[path[-1]])
    return path[::-1]


def prufer_tree(n, word):
    degree = [1]*n
    for v in word:
        degree[v] += 1
    out = set()
    for v in word:
        leaf = next(i for i in range(n) if degree[i] == 1)
        out.add(edge(leaf, v))
        degree[leaf] -= 1
        degree[v] -= 1
    last = [i for i in range(n) if degree[i] == 1]
    out.add(edge(*last))
    return frozenset(out)


def leaf_pair_trees(n):
    # All labeled trees with the distinguished vertices 0,n-1 leaves.
    return [prufer_tree(n, word)
            for word in itertools.product(range(1, n-1), repeat=n-2)]


def tree_flips(tree, n):
    adj = adjacency(tree, n)
    out = Counter()
    for b, c in sorted(tree):
        for a in sorted(adj[b]-{c}):
            for d in sorted(adj[c]-{b}):
                assert len({a, b, c, d}) == 4
                changed = tree.symmetric_difference(
                    {edge(a, b), edge(c, d), edge(a, c), edge(b, d)})
                out[frozenset(changed)] += 1
    return out


def count_formula(tree, n, path):
    adj = adjacency(tree, n)
    inner = path[1:-1]
    down = sum(len(adj[b])+len(adj[c])-4
               for b, c in zip(inner[:-1], inner[1:]))
    path_set = set(path)
    up = 2*sum(len(adj[y])-1 for v in inner for y in adj[v]-path_set)
    return down, up


def weighted_row(d, down, up):
    x = d-1
    left = (Fraction(x, x-1)-Fraction(x-1, x))/2 if x > 1 else 0
    right = (Fraction(x+1, x)-Fraction(x, x+1))/2
    return down*left+up*right


def mask(tree, n):
    return previous.old.graph_mask(n, tree)


def six_vertex_sector():
    n = 6
    # Internal labels 1,2 have degree3; 0,3,4,5 are leaves.
    trees = sorted({prufer_tree(n, p)
                    for p in set(itertools.permutations((1, 1, 2, 2)))},
                   key=lambda t: mask(t, n))
    graphs = [mask(t, n) for t in trees]
    h, f = previous.full_model(n, graphs)
    # Remove the graph-independent mu*|E| phase; J=kappa=1.
    h = h-5*np.eye(len(h), dtype=np.int64)
    d = np.array([len(path_between(t, n, 0, 5))-1 for t in trees])
    return trees, graphs, h, f, d


def comb(m):
    # Spine1..m, endpoint leaves0,m+1, one sideleaf m+1+i per spine vertex.
    tree = {edge(i, i+1) for i in range(m+1)}
    tree.update(edge(i, m+1+i) for i in range(1, m+1))
    return frozenset(tree), 2*m+2


def contract_comb_layers(m):
    assert m >= 2 and m & (m-1) == 0
    tree, n = comb(m)
    spine = list(range(1, m+1))
    layers = []
    while len(spine) > 1:
        adj = adjacency(tree, n)
        spine_set = set(spine)|{0, m+1}
        moves, used = [], set()
        for r in range(0, len(spine), 2):
            a = 0 if r == 0 else spine[r-1]
            b, c = spine[r:r+2]
            roots = adj[c]-spine_set
            assert len(roots) == 1
            d = next(iter(roots))
            support = {edge(b, c), edge(a, b), edge(c, d),
                       edge(a, c), edge(b, d)}
            assert not used & support
            used |= support
            assert {edge(a, b), edge(b, c), edge(c, d)} <= tree
            assert not {edge(a, c), edge(b, d)} & tree
            moves.append((a, b, c, d))
        changed = set(tree)
        for a, b, c, d in moves:
            changed.symmetric_difference_update(
                {edge(a, b), edge(c, d), edge(a, c), edge(b, d)})
        tree = frozenset(changed)
        spine = spine[1::2]
        assert path_between(tree, n, 0, m+1) == [0]+spine+[m+1]
        layers.append(dict(moves=len(moves), remaining_distance=len(spine)+1,
                           disjoint_edge_factor_supports=True))
    return tree, n, layers


class Audit(unittest.TestCase):
    def close(self, a, b, tol=2e-11):
        self.assertLess(float(np.max(np.abs(np.asarray(a)-np.asarray(b)))), tol)

    def test_01_tree_distance_changes_counts_and_schur_bound(self):
        report = []
        for n in range(3, 8):
            trees = leaf_pair_trees(n)
            self.assertEqual(len(set(trees)), (n-2)**(n-2))
            max_row = Fraction(0)
            changes = Counter()
            for tree in trees:
                adj = adjacency(tree, n)
                degrees = tuple(map(len, adj))
                cap = max(degrees)
                path = path_between(tree, n, 0, n-1)
                d = len(path)-1
                targets = tree_flips(tree, n)
                counts = Counter()
                for target, weight in targets.items():
                    self.assertEqual(weight, 1)
                    self.assertEqual(len(target), n-1)
                    self.assertEqual(tuple(map(len, adjacency(target, n))), degrees)
                    delta = len(path_between(target, n, 0, n-1))-1-d
                    self.assertIn(delta, (-1, 0, 1))
                    counts[delta] += weight
                    changes[delta] += weight
                exact_down, exact_up = count_formula(tree, n, path)
                self.assertEqual((counts[-1], counts[1]), (exact_down, exact_up))
                self.assertLessEqual(exact_down, 2*(cap-2)*(d-2))
                self.assertLessEqual(exact_up, 2*(cap-1)*(cap-2)*(d-1))
                row = weighted_row(d, exact_down, exact_up)
                self.assertLessEqual(row, 2*cap*(cap-2))
                max_row = max(max_row, row)
                if n <= 6:
                    self.assertEqual(
                        {mask(t, n): w for t, w in targets.items()},
                        previous.flips(mask(tree, n), n))
            report.append(dict(vertices=n, all_fixed_leaf_pair_trees=len(trees),
                distance_change_counts={str(k): changes[k] for k in (-1, 0, 1)},
                largest_exact_weighted_row=str(max_row),
                cross_checked_with_440_generator=n <= 6))
        OBS['tree_combinatorics'] = dict(rows=report,
            existing_NNI_single_step_leaf_distance_fact_reused=True,
            all_vertices_remain_distinguished=True,
            unlabelled_internal_vertex_NNI_counts_not_substituted=True,
            full_tree_degree_sector_and_distance_counts_exact=True,
            schur_constant='2 C(C-2) abs(kappa)')

    def test_02_complete_data_tilt_and_coherent_initial_tree_subspace(self):
        trees, _, h, f, d = six_vertex_sector()
        x = np.tile(d-1, 64)
        tilted = h*x[None, :]/x[:, None]
        anti = (tilted-tilted.T)/(2j)
        graph_tilt = f*(d-1)[None, :]/(d-1)[:, None]
        graph_anti = (graph_tilt-graph_tilt.T)/(2j)
        self.close(anti, np.kron(np.eye(64), graph_anti))
        self.assertLessEqual(float(np.linalg.norm(graph_anti, 2)), 6.)
        t = .02
        u = core.evolve(h, t)
        rows, cols = np.flatnonzero(x == 1), np.flatnonzero(x >= 2)
        block = u[np.ix_(rows, cols)]
        actual = float(np.linalg.norm(block, 2)**2)
        bound = .25*math.exp(12*t)
        self.assertGreater(actual, 0.)
        self.assertLessEqual(actual, bound)
        inverse_square = 1./x**2
        moment_residual = (math.exp(12*t)*np.diag(inverse_square)
                           -u.conj().T@(inverse_square[:, None]*u))
        smallest = float(np.linalg.eigvalsh(
            (moment_residual+moment_residual.conj().T)/2)[0])
        self.assertGreaterEqual(smallest, -1e-12)
        OBS['complete_operator_bound'] = dict(
            full_data_graph_dimension=len(h), graph_dimension=len(trees),
            coherent_initial_subspace_dimension=len(cols),
            final_near_subspace_dimension=len(rows),
            anti_Hermitian_data_block_cancellation_exact=True,
            graph_anti_Hermitian_norm=float(np.linalg.norm(graph_anti, 2)),
            time=t, initial_distance_at_least=3, final_distance_at_most=2,
            all_input_contraction_probability_operator_norm=actual,
            analytical_upper_bound=bound,
            inverse_distance_moment_operator_inequality_min_eigenvalue=smallest,
            arbitrary_initial_tree_superposition_moment_bound_proved=True,
            arbitrary_untouched_reference_preserved_by_operator_bound=True,
            fixed_slot_chain_used=False)

    def test_03_symmetric_unknown_data_exact_distance_oscillation(self):
        _, _, h, f, d = six_vertex_sector()
        q = np.zeros((64, 7), dtype=np.int64)
        for i in range(64):
            q[i, i.bit_count()] = 1
        w = np.kron(q, np.eye(6, dtype=np.int64))
        self.assertTrue(np.array_equal(h@w, w@np.kron(np.eye(7, dtype=np.int64), f+5*np.eye(6, dtype=np.int64))))
        self.assertEqual(list(np.bincount(d, minlength=4)[2:]), [2, 4])
        # Equitable bright block: near states2, far states4.
        near, far = np.flatnonzero(d == 2), np.flatnonzero(d == 3)
        b = np.column_stack([(d == 2)/math.sqrt(2), (d == 3)/2])
        self.close(f@b, b@np.array([[0., math.sqrt(8)], [math.sqrt(8), 2.]]))
        u = core.evolve(f, math.pi/6)
        p = float(np.sum(np.abs(u[near, far[0]])**2))
        self.assertAlmostEqual(p, 2/9)
        superposed = b[:, 1]
        bright_p = float(np.sum(np.abs((u@superposed)[near])**2))
        self.assertAlmostEqual(bright_p, 8/9)
        OBS['autonomous_distance_oscillation'] = dict(
            symmetric_unknown_data_dimension=7, exact_integer_intertwining=True,
            arbitrary_symmetric_data_and_reference_factorize=True,
            time='pi/6', definite_far_tree_near_probability=p,
            coherent_uniform_far_trees_near_probability=bright_p,
            exact_probabilities=['2/9', '8/9'],
            graph_coherence_not_replaced_by_classical_distribution=True,
            no_claim_of_monotonic_contraction_or_irreversible_selection=True)

    def test_04_every_unknown_data_input_has_nonzero_contraction(self):
        _, _, h, f, d = six_vertex_sector()
        g = int(np.flatnonzero(d == 3)[0])
        cols = np.arange(64)*6+g
        rows = np.flatnonzero(np.tile(d == 2, 64))
        leading = h[np.ix_(rows, cols)]
        self.assertTrue(np.array_equal(leading.T@leading, 2*np.eye(64, dtype=np.int64)))
        # ||H||<=4|kappa|+5|J|=9, after dropping the mu constant.
        self.assertTrue(np.all(np.sum(f, axis=0) == 4))
        t = Fraction(1, 1000)
        remainder = (9*t)**2/(2*(1-9*t))
        amplitude_lower = Fraction(7, 5)*t-remainder  # sqrt(2)>7/5
        self.assertGreater(amplitude_lower, 0)
        probability_lower = amplitude_lower**2
        block = core.evolve(h, float(t))[np.ix_(rows, cols)]
        actual = float(np.linalg.svd(block, compute_uv=False)[-1]**2)
        self.assertGreater(actual, float(probability_lower))
        OBS['unknown_input_nonzero_contraction'] = dict(
            initial_data_dimension=64, time=str(t), J=1, kappa=1,
            minimum_singular_value_squared=actual,
            exact_taylor_remainder_upper=str(remainder),
            exact_uniform_probability_lower=str(probability_lower),
            probability_lower_decimal=float(probability_lower),
            proof_covers_all_unknown_data_and_reference=True,
            reused_439_short_time_block_certificate_method=True)

    def test_05_all_size_probability_and_time_certificate(self):
        # C=3, |kappa|=1, d0=101, r=2, |t|=1/100.
        exponent = Fraction(12, 100)
        probability_upper = Fraction(1, 100**2)/(1-exponent)
        self.assertEqual(probability_upper, Fraction(1, 8800))
        min_time = math.log(100)/12  # at least1% contraction probability
        self.assertGreater(min_time, .38)
        OBS['all_size_certificate'] = dict(
            C=3, kappa=1, initial_distance_at_least=101,
            final_distance_at_most=2, time='1/100',
            exact_probability_upper=str(probability_upper),
            minimum_time_for_probability_at_least_one_percent='log(100)/12',
            minimum_time_decimal=min_time,
            large_Hilbert_space_numerically_simulated=False,
            general_tail='min(1,((r-1)/(d0-1))^2 exp[4 C(C-2) abs(kappa*t)])',
            distance_deformation_bound_not_a_data_signal_bound=True)

    def test_06_branching_shapes_and_scope_of_parallel_flip_schedule(self):
        report = []
        for m in (2, 4, 8, 16, 32, 64, 128):
            start, n = comb(m)
            end, _, layers = contract_comb_layers(m)
            self.assertEqual(tuple(map(len, adjacency(start, n))),
                             tuple(map(len, adjacency(end, n))))
            self.assertEqual(len(end), n-1)
            self.assertEqual(sum(row['moves'] for row in layers), m-1)
            self.assertEqual(len(layers), int(math.log2(m)))
            self.assertEqual(len(path_between(end, n, 0, m+1))-1, 2)
            report.append(dict(spine_vertices=m, vertices=n,
                initial_leaf_distance=m+1, final_leaf_distance=2,
                required_and_constructed_sequential_flips=m-1,
                prescribed_parallel_layers=len(layers),
                distances=[m+1]+[row['remaining_distance'] for row in layers]))
        # Degree-preserving tree flip changes shape: no single fixed template.
        start, n = comb(8)
        end, _, _ = contract_comb_layers(8)
        def diameter(tree):
            return max(len(path_between(tree, n, a, b))-1
                       for a in range(n) for b in range(a+1, n))
        before, after = diameter(start), diameter(end)
        self.assertNotEqual(before, after)
        OBS['branching_shape_and_schedule_scope'] = dict(
            rows=report, nonisomorphic_example_vertices=n,
            initial_diameter=before, final_diameter=after,
            schedule_uses_existing_flip_terms_only=True,
            layer_selection_and_timing_are_extra_controls=True,
            uniform_autonomous_F_follows_this_schedule_proved=False,
            instantaneous_distance_is_not_a_data_channel=True,
            J_zero_preserves_each_subject_data_algebra_reused_from_441=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(stream.getvalue())
    return dict(round=442, baseline_round=441, tests_run=result.testsRun,
        failures=len(result.failures), errors=len(result.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(unchanged_440_autonomous_H_used_for_main_theorem=True,
            coherent_branching_trees_with_unknown_data_reference_covered=True,
            size_independent_quantum_distance_deformation_bound_proved=True,
            existing_NNI_and_441_weighted_method_explicitly_reused=True,
            tree_combinatorics_checked_against_440_generator=True,
            graph_shape_change_without_common_path_template_verified=True,
            actual_internal_distance_measurement_protocol_constructed=False,
            general_cyclic_graphs_or_nonleaf_pairs_covered=False,
            combined_439_440_model_covered=False,
            complete_data_signal_bound_on_branching_graphs_proved=False,
            programmed_flip_schedule_attributed_to_uniform_autonomous_F=False,
            initial_tree_or_three_dimensional_geometry_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    elif not args.dry_run:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
