"""Round 281: source accounting, passive elimination and relative-mode screening.

Finite connected undirected positive-conductance networks are supplied, not
generated. Static Schur complements do not replace the full memory dynamics.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main


def laplacian(n, edges):
    out = np.zeros((n, n))
    for a, b, weight in edges:
        out[a, a] += weight
        out[b, b] += weight
        out[a, b] -= weight
        out[b, a] -= weight
    return out


def cycle(n):
    assert n >= 3
    return laplacian(n, [(i, (i+1) % n, 1.) for i in range(n)])


def supplied_graph(n):
    edges = [(i, i+1, 0.5+(i % 4)/3) for i in range(n-1)]
    edges += [(i, i+3, 0.7) for i in range(n-3)]
    return laplacian(n, edges)


def reduce_network(full, retained):
    retained = list(retained)
    internal = [i for i in range(len(full)) if i not in retained]
    a = full[np.ix_(retained, retained)]
    c = full[np.ix_(retained, internal)]
    d = full[np.ix_(internal, internal)]
    transfer = -np.linalg.solve(d, c.T).T
    reduced = a+transfer @ c.T
    return reduced, transfer, internal


def two_layers(base, coupling):
    identity = np.eye(len(base))
    return np.block([[base+coupling*identity, -coupling*identity],
                     [-coupling*identity, base+coupling*identity]])


def screened_ring_profile(n, coupling):
    gap = 2*coupling
    root = np.sqrt(gap*(gap+4))
    ratio = 2/(gap+2+root)
    j = np.arange(n)
    return (ratio**j+ratio**(n-j))/((1-ratio**n)*root)


def zero_mean_solve(matrix, source):
    assert abs(np.sum(source)) < 1e-10
    # A rank-one lift fixes the otherwise undetermined constant potential.
    return np.linalg.solve(matrix+np.ones_like(matrix)/len(matrix), source)


def connected(n, pairs):
    seen = {0}
    while True:
        expanded = seen | {v for a, b in pairs for u, v in ((a, b), (b, a)) if u in seen}
        if expanded == seen:
            return len(seen) == n
        seen = expanded


def report():
    cases = []
    coupling = 0.25
    for n in (8, 16, 32, 64):
        lam = 4*np.sin(np.pi/n)**2
        reduced_lam = lam*(lam+2*coupling)/(lam+coupling)
        profile = screened_ring_profile(n, coupling)
        cases.append({'cycle_sites_per_layer': n, 'joint_nodes': 2*n,
                      'cycle_first_positive_eigenvalue': float(lam),
                      'eliminated_layer_first_positive_eigenvalue': float(reduced_lam),
                      'clamped_layer_smallest_eigenvalue': coupling,
                      'relative_mode_smallest_eigenvalue': 2*coupling,
                      'neutral_cross_layer_point_source_energy': float(2*profile[0]),
                      'normalized_complete_graph_positive_gap': 1.})
    base = supplied_graph(9)
    retained = [0, 4, 8]
    reduced, transfer, internal = reduce_network(base, retained)
    source = np.zeros(9)
    source[2], source[7] = 1., -1.
    effective = source[retained]+transfer @ source[internal]
    return {'round': 281,
            'scope': 'Static positive-conductance networks with supplied topology and scalar quadratic energy; no physical mass, field, dimension or natural dynamics derived.',
            'coupling': coupling, 'two_layer_cases': cases,
            'source_transfer_example': {'retained': retained, 'internal': internal,
                                       'original_total': float(sum(source)),
                                       'effective_source': effective.tolist(),
                                       'effective_total': float(sum(effective)),
                                       'maximum_reduced_row_sum_error': float(np.max(abs(reduced.sum(axis=1))))},
            'infinite_chain_relative_channel': {'gap': 0.5, 'geometric_decay_ratio': 0.5,
                                               'one_layer_center_potential': 2/3,
                                               'full_quadratic_energy': 4/3},
            'conclusions': ['Passive elimination preserves constant shifts and total source.',
                            'Pinning an internal layer is a different boundary condition.',
                            'A closed network can screen a neutral relative channel while its total channel stays ungapped.',
                            'A single constant zero mode does not imply a vanishing sequence of positive eigenvalues.']}


class Audit(unittest.TestCase):
    def test_01_all_small_connected_graphs_preserve_laplacian(self):
        count = 0
        for n in range(2, 6):
            possible = list(itertools.combinations(range(n), 2))
            for mask in range(1 << len(possible)):
                edges = [p for j, p in enumerate(possible) if mask >> j & 1]
                if not connected(n, edges):
                    continue
                full = laplacian(n, [(a, b, 1.) for a, b in edges])
                reduced, transfer, _ = reduce_network(full, range((n+1)//2))
                np.testing.assert_allclose(reduced.sum(axis=1), 0, atol=1e-12)
                np.testing.assert_allclose(transfer.sum(axis=0), 1, atol=1e-12)
                self.assertGreaterEqual(float(transfer.min()), -1e-12)
                self.assertGreaterEqual(float(np.linalg.eigvalsh(reduced)[0]), -1e-12)
                offdiag = reduced-np.diag(np.diag(reduced))
                self.assertLessEqual(float(offdiag.max()), 1e-12)
                count += 1
        self.assertEqual(count, 771)

    def test_02_internal_source_total_is_transferred(self):
        for n in (5, 9, 16):
            full = supplied_graph(n)
            retained = [0, n//2, n-1]
            _, transfer, internal = reduce_network(full, retained)
            for node in range(n):
                source = np.eye(n)[node]
                effective = source[retained]+transfer @ source[internal]
                self.assertAlmostEqual(sum(effective), sum(source))
                self.assertGreaterEqual(effective.min(), -1e-12)

    def test_03_source_functional_and_shift_invariance(self):
        rng = np.random.default_rng(281)
        full = supplied_graph(11)
        retained = [0, 5, 10]
        reduced, transfer, internal = reduce_network(full, retained)
        c = full[np.ix_(retained, internal)]
        d = full[np.ix_(internal, internal)]
        source = rng.normal(size=11)
        source -= source.mean()
        effective = source[retained]+transfer @ source[internal]
        b = rng.normal(size=3)
        x = np.zeros(11)
        x[retained] = b
        x[internal] = np.linalg.solve(d, source[internal]-c.T @ b)
        f = x @ full @ x/2-source @ x
        expected = b @ reduced @ b/2-effective @ b-source[internal] @ np.linalg.solve(d, source[internal])/2
        self.assertAlmostEqual(f, expected)
        shifted = x+1.73
        self.assertAlmostEqual(shifted @ full @ shifted/2-source @ shifted, f)

    def test_04_reduced_solution_restores_full_kirchhoff_equation(self):
        rng = np.random.default_rng(28104)
        for n in (5, 12, 21):
            full = supplied_graph(n)
            retained = [0, n//2, n-1]
            reduced, transfer, internal = reduce_network(full, retained)
            source = rng.normal(size=n)
            source -= source.mean()
            b = zero_mean_solve(reduced, source[retained]+transfer @ source[internal])
            x = np.zeros(n)
            x[retained] = b
            x[internal] = np.linalg.solve(full[np.ix_(internal, internal)], source[internal]-full[np.ix_(internal, retained)] @ b)
            np.testing.assert_allclose(full @ x, source, atol=1e-11)

    def test_05_total_and_relative_modes_decouple(self):
        for n in (5, 9, 16):
            base = cycle(n)
            identity = np.eye(n)
            transform = np.block([[identity, identity], [identity, -identity]])/np.sqrt(2)
            for coupling in (0.1, 0.25, 1.):
                full = two_layers(base, coupling)
                expected = np.block([[base, np.zeros_like(base)], [np.zeros_like(base), base+2*coupling*identity]])
                np.testing.assert_allclose(transform.T @ full @ transform, expected, atol=1e-13)

    def test_06_passive_layer_elimination_has_no_positive_constant(self):
        coupling = 0.25
        for n in (8, 16, 32):
            base = cycle(n)
            full = two_layers(base, coupling)
            reduced, _, _ = reduce_network(full, range(n))
            lam, vectors = np.linalg.eigh(base)
            expected = lam*(lam+2*coupling)/(lam+coupling)
            np.testing.assert_allclose(vectors.T @ reduced @ vectors, np.diag(expected), atol=1e-13)
            for value, eigen in zip(expected[1:], lam[1:]):
                self.assertGreaterEqual(value/eigen, 1-1e-12)
                self.assertLessEqual(value/eigen, 2+1e-12)

    def test_07_pinning_is_not_elimination_and_has_recorded_return(self):
        n, coupling = 17, 0.25
        base = cycle(n)
        source = np.eye(n)[0]
        pinned = base+coupling*np.eye(n)
        u = np.linalg.solve(pinned, source)
        full = two_layers(base, coupling)
        currents = full @ np.concatenate([u, np.zeros(n)])
        np.testing.assert_allclose(currents[:n], source, atol=1e-13)
        np.testing.assert_allclose(currents[n:], -coupling*u, atol=1e-13)
        self.assertAlmostEqual(sum(currents[n:]), -1.)
        self.assertAlmostEqual(np.linalg.eigvalsh(pinned)[0], coupling)

    def test_08_neutral_relative_source_has_exact_screened_profile(self):
        coupling = 0.25
        for n in (9, 17, 33):
            base = cycle(n)
            full = two_layers(base, coupling)
            b = np.eye(n)[0]
            source = np.concatenate([b, -b])
            profile = screened_ring_profile(n, coupling)
            predicted = np.concatenate([profile, -profile])
            direct = zero_mean_solve(full, source)
            np.testing.assert_allclose(direct, predicted, atol=1e-12)
            self.assertAlmostEqual(predicted @ full @ predicted, 2*profile[0])
        # Infinite chain: exact recurrence with gamma=0.5 and ratio=0.5.
        for j in range(-20, 21):
            at = lambda k: (2/3)*0.5**abs(k)
            self.assertAlmostEqual(2.5*at(j)-at(j-1)-at(j+1), float(j == 0))

    def test_09_joint_total_conserved_while_uniform_difference_decays(self):
        n, coupling = 7, 0.25
        full = two_layers(cycle(n), coupling)
        values, vectors = np.linalg.eigh(full)
        initial = np.concatenate([np.ones(n), np.zeros(n)])
        for time in (0., 0.5, 2., 8.):
            state = vectors @ (np.exp(-time*values)*(vectors.T @ initial))
            self.assertAlmostEqual(sum(state), n)
            np.testing.assert_allclose(state[:n]-state[n:], np.exp(-2*coupling*time), atol=1e-13)

    def test_10_constant_zero_does_not_force_small_positive_gap(self):
        previous = float('inf')
        for n in (8, 16, 32, 64):
            complete = np.eye(n)-np.ones((n, n))/n
            np.testing.assert_allclose(complete @ np.ones(n), 0, atol=1e-13)
            np.testing.assert_allclose(np.linalg.eigvalsh(complete)[1:], 1., atol=1e-13)
            gap = np.linalg.eigvalsh(cycle(n))[1]
            self.assertAlmostEqual(gap, 4*np.sin(np.pi/n)**2)
            self.assertLess(gap, previous)
            previous = gap


if __name__ == '__main__':
    main(__name__, 'conserved_screening_audit', report)
