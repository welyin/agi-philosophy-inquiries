"""Round 266: fixed ancestral interfaces obstruct a uniform 3D bulk cut bound.

The Z^3 comparison is an explicit geometric benchmark, not a derived universe.
Graph volume dimension, cut growth and Hilbert-space dimension are distinct.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import random
import unittest
import numpy as np
from ancestral_quotient_audit import uniform, graph_data, block_diameter
from vertex_split_cycle_audit import (
    EDGE, TRIANGLE, replace_vertex, binary_assignments, edges)


def edge_boundary(graph, selected):
    selected = set(selected)
    return {(u, v) for u, v in edges(graph) if (u in selected) != (v in selected)}


@lru_cache(maxsize=None)
def block_rows():
    result = []
    for extra in range(7):
        graph, words = graph_data(uniform(extra+2))
        block = [i for i, w in enumerate(words) if w[:2] == (0, 1)]
        boundary = len(edge_boundary(graph, block))
        result.append({'extra_depth': extra, 'graph_vertices': len(graph),
                       'block_vertices': len(block), 'outgoing_edges': boundary,
                       'boundary_over_volume_power_two_thirds': boundary/len(block)**(2/3)})
    return result


@lru_cache(maxsize=None)
def mixed_rewrites():
    checks = 0
    sizes = []
    for value in range(16):
        rng = random.Random(value)
        graph = tuple(tuple(j for j in range(4) if j != i) for i in range(4))
        selected = {0}
        for step in range(40):
            v = rng.randrange(len(graph))
            gadget = EDGE if step % 2 else TRIANGLE
            if gadget == EDGE:
                assignment = rng.choice(list(binary_assignments(graph, v)))
            else:
                outputs = rng.sample(range(3), len(graph[v]))
                assignment = dict(zip(graph[v], outputs))
            old_n = len(graph)
            was_inside = v in selected
            graph = replace_vertex(graph, v, gadget, assignment)
            if was_inside:
                selected.update(range(old_n, len(graph)))
            assert len(edge_boundary(graph, selected)) == 3
            checks += 1
        sizes.append(len(selected))
    return {'runs': 16, 'mixed_internal_and_external_rewrites': checks,
            'final_block_sizes': sizes}


def lattice_boundary(selected):
    selected = set(selected)
    result = set()
    for v in selected:
        for axis in range(3):
            for direction in (-1, 1):
                other = list(v)
                other[axis] += direction
                other = tuple(other)
                if other not in selected:
                    result.add((v, other))
    return result


def projections(selected):
    return [len({tuple(v[j] for j in range(3) if j != axis) for v in selected})
            for axis in range(3)]


@lru_cache(maxsize=None)
def lattice_audit():
    points = tuple(itertools.product(range(2), repeat=3))
    count = 0
    for mask in range(1, 1 << len(points)):
        selected = {v for i, v in enumerate(points) if mask >> i & 1}
        p = projections(selected)
        boundary = len(lattice_boundary(selected))
        assert len(selected)**2 <= p[0]*p[1]*p[2]
        assert boundary >= 2*sum(p)
        assert boundary**3 >= 216*len(selected)**2
        count += 1
    return {'nonempty_subsets_of_two_cube': count,
            'boundary_measured_in_infinite_lattice': True}


def matched_volume():
    result = []
    for extra in (3, 6):
        side = 3**(extra//3)
        cube = set(itertools.product(range(side), repeat=3))
        result.append({'vertices_each': len(cube), 'fisher_outgoing_edges': block_rows()[extra]['outgoing_edges'],
                       'lattice_cube_side': side, 'lattice_outgoing_edges': len(lattice_boundary(cube))})
    return result


def report():
    return {'round': 266, 'ancestral_blocks': block_rows(),
            'mixed_gadget_boundary_invariance': mixed_rewrites(),
            'lattice_benchmark': lattice_audit(), 'equal_volume_comparison': matched_volume(),
            'scope': 'Rewrites that only reattach old boundary edges preserve every ancestral cut size. Unbounded descendant blocks therefore fail any uniform positive |boundary| >= c |block|^(2/3) condition. This is a regular 3D bulk benchmark, not a proof against every volume-growth-three graph or every possible emergent spacetime.'}


class Audit(unittest.TestCase):
    def test_01_large_blocks_keep_three_external_edges(self):
        for row in block_rows():
            self.assertEqual(row['outgoing_edges'], 3)
            self.assertEqual(row['block_vertices'], 3**row['extra_depth'])
            self.assertEqual(row['graph_vertices'], 9*row['block_vertices'])

    def test_02_not_an_artifact_of_almost_whole_finite_graph(self):
        for row in block_rows():
            self.assertLess(row['block_vertices'], row['graph_vertices']/2)

    def test_03_adding_binary_gadgets_does_not_amplify_interface(self):
        self.assertEqual(mixed_rewrites()['mixed_internal_and_external_rewrites'], 640)

    def test_04_descendant_blocks_are_connected(self):
        for extra in range(1, 5):
            graph, words = graph_data(uniform(extra+2))
            block = np.array([i for i, w in enumerate(words) if w[:2] == (0, 1)])
            self.assertEqual(block_diameter(graph, block), 2**extra-1)

    def test_05_no_scale_independent_positive_three_dimensional_cut_constant(self):
        rows = block_rows()
        for a, b in zip(rows, rows[1:]):
            self.assertAlmostEqual(b['boundary_over_volume_power_two_thirds']/
                                   a['boundary_over_volume_power_two_thirds'], 3**(-2/3))

    def test_06_lattice_cubes_and_exhaustive_subsets(self):
        self.assertEqual(lattice_audit()['nonempty_subsets_of_two_cube'], 255)
        for side in range(1, 9):
            cube = set(itertools.product(range(side), repeat=3))
            self.assertEqual(len(lattice_boundary(cube)), 6*side**2)

    def test_07_irregular_lattice_sets_obey_same_bound(self):
        rng = random.Random(266)
        points = tuple(itertools.product(range(3), repeat=3))
        for count in range(128):
            selected = set(rng.sample(points, rng.randint(1, len(points))))
            p = projections(selected)
            boundary = len(lattice_boundary(selected))
            self.assertLessEqual(len(selected)**2, p[0]*p[1]*p[2])
            self.assertGreaterEqual(boundary, 2*sum(p))
            self.assertGreaterEqual(boundary**3, 216*len(selected)**2)

    def test_08_equal_volume_does_not_equal_boundary_capacity(self):
        rows = matched_volume()
        self.assertEqual([(x['vertices_each'], x['lattice_outgoing_edges']) for x in rows],
                         [(27, 54), (729, 486)])
        self.assertTrue(all(x['fisher_outgoing_edges'] == 3 for x in rows))


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
    target = Path(__file__).with_name('ancestral_boundary_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
