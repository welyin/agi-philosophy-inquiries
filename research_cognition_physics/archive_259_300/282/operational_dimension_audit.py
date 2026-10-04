"""Round 282: degree-three networks and dimension measured through travel distance.

Supplied decorated lattices are countermodels to degree = dimension, not a
growth law. Every internal routing hop is included in the measured distance.
"""
from collections import deque
import itertools
import math
import unittest
import numpy as np
from growing_stream_audit import main
from fisher_hanoi_scaling_audit import hanoi


def decorated_torus(dimension, side):
    assert dimension >= 2 and side >= 3
    cells = list(itertools.product(range(side), repeat=dimension))
    labels = {x: j for j, x in enumerate(cells)}
    ports = 2*dimension
    rows = []
    for x in cells:
        index = labels[x]
        for p in range(ports):
            y = list(x)
            y[p//2] = (y[p//2]+(1 if p % 2 == 0 else -1)) % side
            neighbor = labels[tuple(y)]*ports+(p ^ 1)
            rows.append((index*ports+(p-1) % ports,
                         index*ports+(p+1) % ports, neighbor))
    return rows, cells


def distances(graph, source):
    result = np.full(len(graph), -1, dtype=int)
    result[source] = 0
    queue = deque([source])
    while queue:
        v = queue.popleft()
        for w in graph[v]:
            if result[w] < 0:
                result[w] = result[v]+1
                queue.append(w)
    return result


def coarse_distances(cells, origin, side):
    a = np.asarray(cells)
    delta = abs(a-np.asarray(origin))
    return np.minimum(delta, side-delta).sum(axis=1)


def lattice_ball(dimension, radius):
    if radius < 0:
        return 0
    return sum(2**j*math.comb(dimension, j)*math.comb(radius, j)
               for j in range(min(dimension, radius)+1))


def pulse_arrival(graph, source, target):
    # One classical probe experiment: uniform edge ticks, reliable relay,
    # no other traffic; elapsed time counts all micro-node forwarding.
    active, visited = {source}, {source}
    elapsed = 0
    while target not in visited:
        active = {w for v in active for w in graph[v]}-visited
        if not active:
            raise ValueError('Unreachable target')
        visited |= active
        elapsed += 1
    return elapsed


def report():
    cases = []
    for d, side in ((2, 9), (3, 7), (4, 5)):
        graph, cells = decorated_torus(d, side)
        fine = distances(graph, 0)
        coarse = np.repeat(coarse_distances(cells, cells[0], side), 2*d)
        cases.append({'supplied_macro_dimension': d, 'periodic_side': side,
                      'micro_nodes': len(graph), 'degree': 3,
                      'micro_edges': 3*len(graph)//2,
                      'cell_size': 2*d, 'maximum_cell_internal_distance': d,
                      'actual_radius_counts': {str(r): int(sum(fine <= r)) for r in (2, 4, 8)},
                      'source_eccentricity': int(max(fine)),
                      'distance_lower_bound_passed': bool(np.all(fine >= coarse)),
                      'distance_upper_bound_passed': bool(np.all(fine <= (d+1)*coarse+d)),
                      'infinite_volume_exponent_proved_by_bounds': d,
                      'coarse_lattice_ball_radius_100': lattice_ball(d, 100)})
    hausdorff = math.log(3)/math.log(2)
    alpha = hausdorff/3
    return {'round': 282,
            'scope': 'Operational dimension countermodels under supplied topology, clock ticks and counting measure; no selection or derivation of physical three-dimensional space.',
            'degree_three_examples': cases,
            'old_model_dimension': {'volume_exponent': hausdorff,
                                    'spectral_dimension_from_round262': 2*math.log(3)/math.log(5)},
            'nonlinear_metric_relabeling': {'alpha_chosen_to_force_three': alpha,
                                           'resulting_volume_exponent': hausdorff/alpha,
                                           'two_unit_edges_additive_delay': 2.,
                                           'relabeled_two_step_distance': 2**alpha,
                                           'has_unit_edge_travel_time_realization_on_original_graph': False},
            'conclusion': 'Degree bounds do not select operational dimension; bounded local grouping preserves volume exponent, while nonlinear relabeling needs a new measurement mechanism.'}


class Audit(unittest.TestCase):
    def test_01_every_micro_node_has_three_neighbors(self):
        for d, side in ((2, 5), (3, 4), (4, 3)):
            graph, _ = decorated_torus(d, side)
            self.assertTrue(all(len(set(row)) == 3 for row in graph))
            self.assertTrue(all(v != w and v in graph[w] for v, row in enumerate(graph) for w in row))
            self.assertTrue(np.all(distances(graph, 0) >= 0))

    def test_02_quotient_has_exactly_the_supplied_lattice_edges(self):
        for d, side in ((2, 5), (3, 4), (4, 3)):
            graph, cells = decorated_torus(d, side)
            index = {x: j for j, x in enumerate(cells)}
            for cell, x in enumerate(cells):
                actual = {w//(2*d) for p in range(2*d) for w in graph[cell*2*d+p]
                          if w//(2*d) != cell}
                expected = set()
                for axis in range(d):
                    for sign in (-1, 1):
                        y = list(x)
                        y[axis] = (y[axis]+sign) % side
                        expected.add(index[tuple(y)])
                self.assertEqual(actual, expected)

    def test_03_all_internal_routing_is_in_distance_bounds(self):
        for d, side in ((2, 5), (3, 4), (4, 3)):
            graph, cells = decorated_torus(d, side)
            for source in (0, 2*d-1, len(graph)//2):
                actual = distances(graph, source)
                coarse = np.repeat(coarse_distances(cells, cells[source//(2*d)], side), 2*d)
                self.assertTrue(np.all(coarse <= actual))
                self.assertTrue(np.all(actual <= (d+1)*coarse+d))

    def test_04_micro_ball_counts_obey_quotient_sandwich(self):
        for d, side in ((2, 9), (3, 7), (4, 5)):
            graph, cells = decorated_torus(d, side)
            actual = distances(graph, 0)
            coarse = coarse_distances(cells, cells[0], side)
            for radius in range(int(max(actual))+1):
                lower_radius = (radius-d)//(d+1)
                lower = 2*d*sum(coarse <= lower_radius)
                upper = 2*d*sum(coarse <= radius)
                self.assertLessEqual(lower, sum(actual <= radius))
                self.assertLessEqual(sum(actual <= radius), upper)

    def test_05_lattice_ball_polynomial_is_independently_counted(self):
        for d in (2, 3, 4):
            for radius in range(1, 5):
                direct = sum(sum(abs(x) for x in point) <= radius
                             for point in itertools.product(range(-radius, radius+1), repeat=d))
                self.assertEqual(lattice_ball(d, radius), direct)
            # Leading coefficient 2^d/d! via exact dth differences.
            sequence = [lattice_ball(d, r) for r in range(d+1)]
            for _ in range(d):
                sequence = [b-a for a, b in zip(sequence, sequence[1:])]
            self.assertEqual(sequence, [2**d])

    def test_06_round_trip_probe_counts_all_hops(self):
        graph, _ = decorated_torus(3, 4)
        source = 0
        computed = distances(graph, source)
        for target in (1, 6, 37, len(graph)-1):
            outward = pulse_arrival(graph, source, target)
            inward = pulse_arrival(graph, target, source)
            self.assertEqual(outward+inward, 2*computed[target])

    def test_07_snowflake_is_a_metric_but_changes_path_additivity(self):
        graph, _ = hanoi(3)
        matrix = np.stack([distances(graph, v) for v in range(len(graph))])
        alpha = math.log(3)/math.log(2)/3
        transformed = matrix.astype(float)**alpha
        for middle in range(len(graph)):
            self.assertTrue(np.all(transformed <= transformed[:, middle, None]+transformed[None, middle, :]+1e-12))
        self.assertLess(2**alpha, 2.)
        self.assertTrue(np.all(transformed[matrix == 1] == 1.))
        self.assertTrue(np.all(transformed[matrix == 2] < matrix[matrix == 2]))

    def test_08_transformed_balls_are_only_rethresholded_old_balls(self):
        graph, _ = hanoi(4)
        old = distances(graph, 0)
        alpha = math.log(3)/math.log(2)/3
        for radius in (1.2, 2.3, 3.4, 4.1):
            lhs = set(np.flatnonzero(old.astype(float)**alpha <= radius))
            rhs = set(np.flatnonzero(old <= radius**(1/alpha)))
            self.assertEqual(lhs, rhs)
        self.assertAlmostEqual((math.log(3)/math.log(2))/alpha, 3.)


if __name__ == '__main__':
    main(__name__, 'operational_dimension_audit', report)
