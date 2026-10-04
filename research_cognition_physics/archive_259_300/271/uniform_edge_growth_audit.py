"""Round 271: fair uniform mutual proposals can grow cuts too slowly.

Exact moment identities use a supplied cubic seed, not the degree-two boundary
of the earlier single-edge seed. Logical proposal rounds are not physical time.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import random
import unittest
import numpy as np
from paired_interface_growth_audit import port_state, refine, check_ports
from ancestral_quotient_audit import graph_data
from contact_route_closure_audit import adjacency
from vertex_split_cycle_audit import edges

A = 4/3
B = 10/9


def prism(length):
    edge_list = [(rail*length+i, rail*length+(i+1)%length)
                 for rail in (0, 1) for i in range(length)]
    edge_list += [(i, length+i) for i in range(length)]
    return adjacency(2*length, edge_list)


def cube():
    return adjacency(8, [(v, v^(1 << bit)) for v in range(8)
                         for bit in range(3) if v < (v^(1 << bit))])


def graph_cases():
    return [('K4', adjacency(4, itertools.combinations(range(4), 2))),
            ('triangular_prism', prism(3)),
            ('K33', adjacency(6, itertools.product(range(3), range(3, 6)))),
            ('cube', cube())]


def mutual_round(state, rng):
    choices = {u: rng.choice(tuple(sorted(row))) for u, row in sorted(state.items())}
    partners = {u: state[u][p][0] for u, p in choices.items()}
    matching = [(u, v) for u, v in sorted(partners.items()) if u < v and partners[v] == u]
    for u, v in matching:
        refine(state, u, v)
    return matching


def cut_data(state, selected_roots):
    selected = {u for u in state if u[0] in selected_roots}
    boundary = sum((u in selected) != (v in selected)
                   for u, row in state.items() for v, _ in row.values() if u < v)
    return len(state), len(selected), boundary


def exact_weighted_variance(graph, weights):
    es = edges(graph)
    square_sum = sum(w*w for w in weights)
    adjacent_products = sum(weights[i]*weights[j] for i, e in enumerate(es)
                            for j, f in enumerate(es) if i < j and set(e) & set(f))
    return Fraction(8*square_sum-2*adjacent_products, 81)


@lru_cache(maxsize=None)
def exact_cases():
    result = []
    for name, graph in graph_cases():
        es = edges(graph)
        proposals = list(itertools.product(*graph))
        indicators = np.array([[int(choice[u] == v and choice[v] == u) for u, v in es]
                               for choice in proposals], dtype=int)
        for selected in ({0}, set(range(len(graph)//2))):
            volume_w = np.array([int(u in selected)+int(v in selected) for u, v in es])
            cut_w = np.array([int((u in selected) != (v in selected)) for u, v in es])
            delta_volume = indicators@volume_w
            delta_cut = indicators@cut_w
            mean_v = Fraction(int(delta_volume.sum()), len(proposals))
            mean_c = Fraction(int(delta_cut.sum()), len(proposals))
            variance = Fraction(int((delta_volume**2).sum()), len(proposals))-mean_v**2
            expected_variance = exact_weighted_variance(graph, list(map(int, volume_w)))
            assert mean_v == Fraction(len(selected), 3)
            assert mean_c == Fraction(int(cut_w.sum()), 9)
            assert variance == expected_variance
            result.append({'graph': name, 'vertices': len(graph), 'selected_vertices': len(selected),
                           'cut': int(cut_w.sum()), 'proposal_configurations': len(proposals),
                           'mean_volume_increment': str(mean_v), 'mean_cut_increment': str(mean_c),
                           'variance_volume_increment': str(variance)})
    return result


@lru_cache(maxsize=None)
def sampled_histories():
    repeats, rounds, length = 128, 8, 64
    initial = prism(length)
    selected = {rail*length+i for rail in (0, 1) for i in range(length//2)}
    measurements = {t: [] for t in (0, 4, 8)}
    operations = 0
    for value in range(repeats):
        rng = random.Random(271000+value)
        state = port_state(initial)
        measurements[0].append(cut_data(state, selected))
        for t in range(1, rounds+1):
            matching = mutual_round(state, rng)
            operations += len(matching)
            if t in measurements:
                check_ports(state)
                assert all(len(row) == 3 for row in state.values())
                measurements[t].append(cut_data(state, selected))
    rows = []
    for t, entries in measurements.items():
        values = np.array(entries, dtype=float)
        prediction = np.array([128*A**t, 64*A**t, 4*B**t])
        rows.append({'round': t, 'sample_means_N_S_cut': values.mean(axis=0).tolist(),
                     'predicted_means_N_S_cut': prediction.tolist(),
                     'sample_standard_errors_N_S_cut': (values.std(axis=0, ddof=1)/np.sqrt(repeats)).tolist()})
    return {'seed': '128-vertex circular ladder, 64 vertices on each ancestral side',
            'repeats': repeats, 'simulated_rounds': rounds, 'paired_rewrites': operations,
            'rows': rows}


def failure_bound(t, side_size=64, initial_cut=4, c=0.1):
    ratio = B/A**(2/3)
    return max(0.0, 1-12/side_size-initial_cut/(c*(side_size/2)**(2/3))*ratio**t)


def report():
    return {'round': 271, 'exact_proposal_enumeration': exact_cases(),
            'uniform_edge_selection_probability': str(Fraction(1, 9)),
            'volume_mean_multiplier': str(Fraction(4, 3)),
            'cut_mean_multiplier': str(Fraction(10, 9)),
            'exponent_relating_the_two_means_not_a_spatial_dimension': math.log(B)/math.log(A),
            'cut_to_three_dimensional_benchmark_mean_rate_ratio': B/A**(2/3),
            'monte_carlo_moment_checks': sampled_histories(),
            'proved_failure_probability_lower_bounds': [
                {'logical_round': t, 'initial_vertices_per_side': 64, 'initial_cut': 4,
                 'benchmark_c': 0.1, 'failure_probability_at_least': failure_bound(t)}
                for t in (20, 40, 60, 80, 120)],
            'large_round_bounds_are_analytic_not_large_graph_simulations': True,
            'scope': 'With an explicitly cubic seed and independent uniform neighbor proposals, every persistent edge is eventually selected almost surely. Nevertheless the mean cut grows by 10/9 while each side grows by 4/3; a second-moment and Markov bound gives finite-time failure probabilities for a uniform regular-3D cut benchmark. This is not a dimension estimate or an impossibility theorem for all local growth laws.'}


class Audit(unittest.TestCase):
    def test_01_exact_conditional_moments(self):
        self.assertEqual(len(exact_cases()), 8)
        self.assertEqual(sum(row['proposal_configurations'] for row in exact_cases()), 16200)

    def test_02_independence_for_disjoint_edges_and_exclusion_for_incident_edges(self):
        graph = cube()
        es = edges(graph)
        total = 3**len(graph)
        counts = np.zeros((len(es), len(es)), dtype=int)
        for choices in itertools.product(*graph):
            vector = np.array([int(choices[u] == v and choices[v] == u) for u, v in es])
            counts += np.outer(vector, vector)
        for i, e in enumerate(es):
            for j, f in enumerate(es):
                expected = Fraction(1, 9) if i == j else (Fraction(0) if set(e)&set(f) else Fraction(1, 81))
                self.assertEqual(Fraction(int(counts[i, j]), total), expected)

    def test_03_variance_bound_for_all_small_ancestral_subsets(self):
        for _, graph in graph_cases():
            es = edges(graph)
            for mask in range(1, 1 << len(graph)):
                selected = {i for i in range(len(graph)) if mask >> i & 1}
                weights = [int(u in selected)+int(v in selected) for u, v in es]
                variance = exact_weighted_variance(graph, weights)
                self.assertGreaterEqual(variance, 0)
                self.assertLessEqual(variance, Fraction(2*len(selected), 3))

    def test_04_every_matching_preserves_cubic_degree_and_cut_increment(self):
        graph = cube()
        selected = set(range(4))
        for value in range(32):
            state = port_state(graph)
            before = cut_data(state, selected)
            matching = mutual_round(state, random.Random(value))
            after = cut_data(state, selected)
            self.assertTrue(all(len(row) == 3 for row in state.values()))
            self.assertEqual(after[0]-before[0], 2*len(matching))
            self.assertEqual(after[1]-before[1], sum(int(u[0] in selected)+int(v[0] in selected) for u, v in matching))
            self.assertEqual(after[2]-before[2], sum((u[0] in selected)!=(v[0] in selected) for u, v in matching))

    def test_05_sampled_moments_match_predictions_with_reported_uncertainty(self):
        for row in sampled_histories()['rows']:
            for actual, expected, sem in zip(row['sample_means_N_S_cut'],
                                             row['predicted_means_N_S_cut'],
                                             row['sample_standard_errors_N_S_cut']):
                self.assertLessEqual(abs(actual-expected), 6*sem+1e-10)

    def test_06_second_moment_geometric_sum(self):
        for t in (1, 2, 8, 40):
            direct = sum(Fraction(2, 3)*Fraction(4, 3)**j/Fraction(4, 3)**(2*j+2) for j in range(t))
            closed = Fraction(3, 2)*(1-Fraction(3, 4)**t)
            self.assertEqual(direct, closed)

    def test_07_failure_bound_uses_rates_not_a_dimension_fit(self):
        self.assertLess(B/A**(2/3), 1)
        self.assertGreater(failure_bound(60), 0.79)
        self.assertLess(failure_bound(60), 0.80)
        self.assertGreater(failure_bound(120), failure_bound(60))
        self.assertLessEqual(failure_bound(120), 1-12/64)

    def test_08_persistent_edge_survival_tail(self):
        # Track an old edge through reattachment, rather than treating renamed
        # endpoints as disappearance or resetting its waiting time.
        repeats, rounds = 128, 8
        survived = np.zeros(rounds, dtype=int)
        seed_graph = graph_cases()[0][1]
        for value in range(repeats):
            rng = random.Random(2719000+value)
            state = port_state(seed_graph)
            old_u, old_v = (0,), (1,)
            for t in range(rounds):
                matching = mutual_round(state, rng)
                if tuple(sorted((old_u, old_v))) in matching:
                    break
                left = [old_u] if old_u in state else [w for w in state if w[:-1] == old_u]
                right = {old_v} if old_v in state else {w for w in state if w[:-1] == old_v}
                inherited = [(u, v) for u in left for v, _ in state[u].values() if v in right]
                self.assertEqual(len(inherited), 1)
                old_u, old_v = inherited[0]
                survived[t] += 1
        for t, count in enumerate(survived, start=1):
            p = (8/9)**t
            self.assertLessEqual(abs(int(count)-repeats*p), 6*np.sqrt(repeats*p*(1-p)))

    def test_09_channel_and_storage_ledger(self):
        for row in sampled_histories()['rows']:
            n, side, cut = row['predicted_means_N_S_cut']
            expected_operations = (n-128)/2
            expected_edges = 192+3*expected_operations
            self.assertAlmostEqual(expected_edges, 1.5*n)
            self.assertAlmostEqual(side, n/2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['runtime'] = {'python': platform.python_version(), 'numpy': np.__version__}
    result['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    payload = json.dumps(result, ensure_ascii=False, indent=2)+'\n'
    target = Path(__file__).with_name('uniform_edge_growth_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
