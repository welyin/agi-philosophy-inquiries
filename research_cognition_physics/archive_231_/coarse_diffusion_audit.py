"""Round 265: bounded-fiber gap comparison does not imply exact Markov closure.

The added dynamics is unit-conductance classical diffusion with generator -L.
Rayleigh/Dirichlet comparison and strong lumpability are established methods.
"""
import argparse
from functools import lru_cache
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
from ancestral_quotient_audit import (
    uniform, partition, graph_data, block_diameter, random_snapshot)
from hanoi_transport_spectrum_audit import gap_formula


def laplacian(graph):
    matrix = np.zeros((len(graph), len(graph)))
    for i, row in enumerate(graph):
        matrix[i, i] = len(row)
        matrix[i, list(row)] = -1
    return matrix


def aggregation(state, depth):
    graph, words, ancestors, labels, blocks = partition(state, depth)
    c = np.zeros((len(blocks), len(graph)))
    c[labels, np.arange(len(graph))] = 1
    masses = c.sum(axis=1)
    u = c.T/masses[None, :]
    coarse_l = laplacian(graph_data(uniform(depth))[0])
    return laplacian(graph), coarse_l, c, u, masses, graph, words, labels, blocks


def heat(matrix, t):
    values, vectors = np.linalg.eigh(matrix)
    return (vectors*np.exp(-t*values)[None, :])@vectors.T


@lru_cache(maxsize=None)
def gap_case(depth, lag, value):
    state = random_snapshot(depth, lag, value)
    fine, coarse, _, _, masses, graph, _, _, blocks = aggregation(state, depth)
    size = int(max(masses))
    diameter = max(block_diameter(graph, block) for block in blocks)
    fine_gap = float(np.linalg.eigvalsh(fine)[1])
    base_gap = float(np.linalg.eigvalsh(coarse)[1])
    factor = 1 if diameter == 0 else 12*size*diameter
    return {'depth': depth, 'allowed_lag': lag, 'seed': value,
            'fine_vertices': len(graph), 'max_fiber_size': size,
            'max_intrinsic_fiber_diameter': diameter,
            'fine_gap': fine_gap, 'coarse_gap': base_gap,
            'ratio': fine_gap/base_gap, 'proved_ratio_lower_bound': 1/factor}


def energy_parts(state, depth, f):
    fine, coarse, _, _, masses, graph, _, labels, blocks = aggregation(state, depth)
    means = np.array([float(np.mean(f[block])) for block in blocks])
    within = sum(float(np.sum((f[block]-means[i])**2)) for i, block in enumerate(blocks))
    between = float(np.sum(masses*(means-np.mean(f))**2))
    internal = external = 0.
    for i, row in enumerate(graph):
        for j in row:
            if i < j:
                contribution = float((f[i]-f[j])**2)
                if labels[i] == labels[j]:
                    internal += contribution
                else:
                    external += contribution
    return {'within': within, 'between': between,
            'variance': float(np.sum((f-np.mean(f))**2)),
            'internal': internal, 'external': external,
            'coarse_energy_of_means': float(means@coarse@means),
            'K': int(max(masses)),
            'D': max(block_diameter(graph, block) for block in blocks),
            'coarse_gap': float(np.linalg.eigvalsh(coarse)[1])}


@lru_cache(maxsize=None)
def memory_example():
    fine, coarse, c, u, masses, _, words, _, _ = aggregation(uniform(2), 1)
    q = -fine
    a = c@q@u
    inside, exit_vertex = words.index((0, 0)), words.index((0, 1))
    one, two = np.eye(len(words))[:, inside], np.eye(len(words))[:, exit_vertex]
    coarse_t = c@heat(fine, .5)@u
    coarse_2t = c@heat(fine, 1.)@u
    memory_curvature = c@q@q@u-a@a
    return {'same_initial_coarse_state': c@one,
            'derivative_from_interior': c@q@one,
            'derivative_from_exit': c@q@two,
            'instantaneous_uniform_lift_generator': a,
            'second_derivative_memory_term': memory_curvature,
            'coarse_transition_at_half': coarse_t,
            'markov_approximation_at_half': heat(coarse/3, .5),
            'semigroup_defect_max': float(np.max(np.abs(coarse_2t-coarse_t@coarse_t))),
            'uniform_initialization_error_max': float(np.max(np.abs(coarse_t-heat(coarse/3, .5))))}


def report():
    memory = {k: v.tolist() if isinstance(v, np.ndarray) else v
              for k, v in memory_example().items()}
    return {'round': 265,
            'gap_comparisons': [gap_case(n, b, s)
                                for n, b, s in itertools.product((1, 2, 3, 4), (1, 2), (17, 29))],
            'gap_comparison_bound': 'lambda_Q/(12 K D) <= lambda_G <= lambda_Q, D>=1',
            'relaxation_exponent_for_uniformly_bounded_lag': math.log(5)/math.log(2),
            'memory_counterexample': memory,
            'scope': 'Uniform bounds on connected fiber size and intrinsic diameter preserve the spectral-gap scale for unit-conductance classical diffusion. Exact block-mass dynamics is nevertheless generally non-Markovian. No full density-of-states or heat-kernel theorem, quantum Hamiltonian, physical clock or gravity is inferred.'}


class Audit(unittest.TestCase):
    def test_01_gaps_obey_two_sided_comparison(self):
        for n, b, value in itertools.product((1, 2, 3, 4), (1, 2), (17, 29)):
            row = gap_case(n, b, value)
            self.assertGreaterEqual(row['ratio']+1e-10, row['proved_ratio_lower_bound'])
            self.assertLessEqual(row['ratio'], 1+1e-10)

    def test_02_uniform_refinement_recovers_existing_gap_formula(self):
        for n, b in itertools.product((1, 2, 3), (0, 1, 2)):
            matrix = laplacian(graph_data(uniform(n+b))[0])
            self.assertAlmostEqual(float(np.linalg.eigvalsh(matrix)[1]), gap_formula(n+b), delta=1e-11)

    def test_03_variance_decomposition_and_internal_bound(self):
        rng = np.random.default_rng(265)
        for n, b in itertools.product((1, 2, 3), (1, 2)):
            state = random_snapshot(n, b, 19)
            for _ in range(8):
                row = energy_parts(state, n, rng.normal(size=len(state)))
                self.assertAlmostEqual(row['variance'], row['within']+row['between'], places=10)
                self.assertLessEqual(row['within'], row['K']*row['D']*row['internal']/2+1e-10)

    def test_04_boundary_energy_and_final_poincare_bound(self):
        rng = np.random.default_rng(1265)
        for n, b in itertools.product((1, 2, 3), (1, 2)):
            state = random_snapshot(n, b, 29)
            for _ in range(8):
                row = energy_parts(state, n, rng.normal(size=len(state)))
                self.assertLessEqual(row['coarse_energy_of_means'],
                                     3*row['external']+9*row['D']*row['internal']+1e-10)
                bound = 12*row['K']*row['D']/row['coarse_gap']*(row['internal']+row['external'])
                self.assertLessEqual(row['variance'], bound+1e-10)

    def test_05_lifted_energy_equals_coarse_energy(self):
        for n in (1, 2, 3):
            state = random_snapshot(n, 2, 31)
            fine, coarse, c, _, _, *_ = aggregation(state, n)
            np.testing.assert_allclose(c@fine@c.T, coarse, atol=1e-13)

    def test_06_mass_weighted_generator_is_only_initial_closure(self):
        for n in (1, 2, 3):
            fine, coarse, c, u, masses, *_ = aggregation(random_snapshot(n, 2, 41), n)
            np.testing.assert_allclose(c@u, np.eye(len(masses)), atol=1e-14)
            a = c@(-fine)@u
            np.testing.assert_allclose(a, -coarse@np.diag(1/masses), atol=1e-14)
            np.testing.assert_allclose(a@masses, 0, atol=1e-14)

    def test_07_same_coarse_initial_state_different_derivative(self):
        result = memory_example()
        np.testing.assert_allclose(result['same_initial_coarse_state'], [1, 0, 0])
        np.testing.assert_allclose(result['derivative_from_interior'], [0, 0, 0])
        np.testing.assert_allclose(result['derivative_from_exit'], [-1, 1, 0])

    def test_08_memory_curvature_is_nonzero(self):
        fine, _, c, u, *_ = aggregation(uniform(2), 1)
        q = -fine
        expected = c@q@(np.eye(len(fine))-u@c)@q@u
        np.testing.assert_allclose(memory_example()['second_derivative_memory_term'], expected, atol=1e-14)
        self.assertGreater(float(np.max(np.abs(expected))), .1)

    def test_09_projected_heat_is_stochastic_but_not_semigroup(self):
        result = memory_example()
        self.assertGreaterEqual(float(result['coarse_transition_at_half'].min()), -1e-14)
        np.testing.assert_allclose(result['coarse_transition_at_half'].sum(axis=0), 1, atol=1e-14)
        self.assertGreater(result['semigroup_defect_max'], 1e-3)
        self.assertGreater(result['uniform_initialization_error_max'], 1e-3)

    def test_10_singleton_blocks_have_exact_closure(self):
        fine, _, c, u, *_ = aggregation(uniform(2), 2)
        np.testing.assert_allclose(c@heat(fine, .8)@u,
                                   (c@heat(fine, .4)@u)@(c@heat(fine, .4)@u), atol=1e-13)
        np.testing.assert_allclose(c@fine@fine@u-(c@fine@u)@(c@fine@u), 0, atol=1e-13)


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
    target = Path(__file__).with_name('coarse_diffusion_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
