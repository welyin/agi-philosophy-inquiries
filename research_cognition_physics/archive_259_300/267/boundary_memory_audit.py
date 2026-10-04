"""Round 267: exact boundary memory, static Kron response and low-frequency mass.

All dynamics here is classical unit-conductance diffusion on fixed graphs.
The three retained corner values are not the block probabilities of round 265.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np
from ancestral_quotient_audit import uniform, graph_data
from coarse_diffusion_audit import laplacian, heat


@lru_cache(maxsize=None)
def case(depth):
    graph, words = graph_data(uniform(depth))
    full = laplacian(graph)
    boundary = np.array([words.index((i,)*depth) for i in range(3)])
    inside = np.array([i for i in range(len(graph)) if i not in boundary], dtype=int)
    bb = full[np.ix_(boundary, boundary)]
    bi = full[np.ix_(boundary, inside)]
    dd = full[np.ix_(inside, inside)]
    if len(inside):
        values, vectors = np.linalg.eigh(dd)
        extension = -np.linalg.solve(dd, bi.T)
    else:
        values, vectors = np.empty(0), np.empty((0, 0))
        extension = np.empty((0, 3))
    kron = bb+bi@extension
    mass = np.eye(3)+extension.T@extension
    return dict(depth=depth, full=full, boundary=boundary, inside=inside, bb=bb,
                bi=bi, dd=dd, values=values, vectors=vectors,
                extension=extension, kron=kron, mass=mass)


def internal_heat(data, time):
    return (data['vectors']*np.exp(-time*data['values']))@data['vectors'].T


def response(data, frequency):
    if not len(data['inside']):
        return frequency*np.eye(3)+data['bb']
    return (frequency*np.eye(3)+data['bb']-
            data['bi']@np.linalg.solve(frequency*np.eye(len(data['inside']))+
                                      data['dd'], data['bi'].T))


def kernel(data, time):
    return data['bi']@internal_heat(data, time)@data['bi'].T


@lru_cache(maxsize=None)
def convolution_audit():
    nodes, weights = np.polynomial.legendre.leggauss(48)
    rows = []
    for depth in (2, 3):
        data = case(depth)
        # A genuinely occupied interior tests the otherwise omitted forcing.
        initial = np.zeros(len(data['full']))
        initial[data['boundary'][0]] = 0.4
        initial[data['inside'][0]] = 0.6
        values, vectors = np.linalg.eigh(data['full'])
        coefficients = vectors.T@initial
        def state(t):
            return vectors@(np.exp(-t*values)*coefficients)
        for time in (0.2, 1.5):
            now = state(time)
            convolution = np.zeros(3)
            for node, weight in zip(nodes, weights):
                past = time*(node+1)/2
                convolution += time*weight/2*(
                    kernel(data, time-past)@state(past)[data['boundary']])
            forcing = -data['bi']@internal_heat(data, time)@initial[data['inside']]
            predicted = -data['bb']@now[data['boundary']]+forcing+convolution
            actual = (-data['full']@now)[data['boundary']]
            rows.append({'depth': depth, 'time': time,
                         'full_vs_memory_derivative_error': float(np.max(np.abs(actual-predicted))),
                         'error_if_initial_interior_forcing_omitted':
                         float(np.max(np.abs(actual-(predicted-forcing))))})
    return rows


def frequency_rows(data):
    gap = float(data['values'][0])
    norm_t_squared = float(np.linalg.norm(data['extension'], 2)**2)
    result = []
    for ratio in (0.001, 0.01, 0.1):
        s = gap*ratio
        exact = response(data, s)
        low = data['kron']+s*data['mass']
        remainder = exact-low
        bound = s*s*norm_t_squared/(gap+s)
        result.append({'s_over_grounded_gap': ratio, 's': s,
                       'identity_capacity_response_error':
                       float(np.linalg.norm(exact-(data['kron']+s*np.eye(3)), 2)),
                       'effective_capacity_response_error': float(np.linalg.norm(remainder, 2)),
                       'proved_remainder_bound': bound})
    return result


def report():
    rows = []
    for depth in range(2, 6):
        data = case(depth)
        resistance = (5/3)**depth-1
        rows.append({'depth': depth, 'vertices': len(data['full']),
                     'grounded_interior_gap': float(data['values'][0]),
                     'kron_edge_conductance': float(-data['kron'][0, 1]),
                     'conductance_from_previous_resistance': 2/(3*resistance),
                     'total_effective_capacity': float(data['mass'].sum()),
                     'frequency_checks': frequency_rows(data)})
    return {'round': 267, 'corner_reductions': rows,
            'nine_vertex_effective_capacity': case(2)['mass'].tolist(),
            'time_domain_memory_checks': convolution_audit(),
            'nine_vertex_stationary_boundary_values': [1/9]*3,
            'naive_three_vertex_stationary_values': [1/3]*3,
            'scope': 'Exact linear elimination gives a memory kernel and an initial-interior forcing. Static Kron reduction preserves boundary Dirichlet response, not the full heat process. A dense effective capacity is the first low-frequency correction with a proved remainder; no quantum dynamics or three-dimensional space is derived.'}


class Audit(unittest.TestCase):
    def test_01_grounded_and_kron_operators(self):
        for depth in range(2, 6):
            data = case(depth)
            self.assertGreater(data['values'][0], 0)
            self.assertLessEqual(data['values'][0], 6/len(data['inside'])+1e-13)
            np.testing.assert_allclose(data['kron'], data['kron'].T, atol=1e-13)
            np.testing.assert_allclose(data['kron'].sum(axis=1), 0, atol=2e-13)
            self.assertGreater(np.linalg.eigvalsh(data['kron'])[1], 0)

    def test_02_dirichlet_minimization(self):
        rng = np.random.default_rng(267)
        for depth in (2, 3, 4):
            data = case(depth)
            for _ in range(10):
                b = rng.normal(size=3)
                harmonic = data['extension']@b
                x = np.zeros(len(data['full']))
                x[data['boundary']] = b
                x[data['inside']] = harmonic
                np.testing.assert_allclose(data['dd']@harmonic+data['bi'].T@b, 0, atol=1e-12)
                self.assertAlmostEqual(float(x@data['full']@x), float(b@data['kron']@b), places=10)
                other = x.copy()
                other[data['inside']] += rng.normal(size=len(harmonic))
                self.assertGreaterEqual(float(other@data['full']@other), float(x@data['full']@x))

    def test_03_effective_resistance_matches_frozen_scaling(self):
        for depth in range(2, 6):
            data = case(depth)
            expected = (5/3)**depth-1
            np.testing.assert_allclose(data['kron'], 2/(3*expected)*(3*np.eye(3)-np.ones((3, 3))),
                                       rtol=1e-10, atol=1e-13)
            removed = data['boundary'][1]
            keep = np.array([i for i in range(len(data['full'])) if i != removed])
            source = np.zeros(len(keep))
            first = int(np.flatnonzero(keep == data['boundary'][0])[0])
            source[first] = 1
            voltage = np.linalg.solve(data['full'][np.ix_(keep, keep)], source)
            self.assertAlmostEqual(float(voltage[first]), expected, places=9)

    def test_04_full_resolvent_matches_dynamic_schur_response(self):
        for depth in (2, 3, 4):
            data = case(depth)
            for s in (0.003, 0.07, 1.0):
                inverse = np.linalg.inv(s*np.eye(len(data['full']))+data['full'])
                np.testing.assert_allclose(inverse[np.ix_(data['boundary'], data['boundary'])],
                                           np.linalg.inv(response(data, s)), rtol=2e-10, atol=1e-10)

    def test_05_time_convolution_and_initial_interior_forcing(self):
        for row in convolution_audit():
            self.assertLess(row['full_vs_memory_derivative_error'], 2e-12)
            self.assertGreater(row['error_if_initial_interior_forcing_omitted'], 0.01)

    def test_06_exact_remainder_identity_sign_and_bound(self):
        for depth in range(2, 6):
            data = case(depth)
            for row in frequency_rows(data):
                s = row['s']
                remainder = response(data, s)-data['kron']-s*data['mass']
                expected = -s*s*data['extension'].T@np.linalg.solve(
                    s*np.eye(len(data['inside']))+data['dd'], data['extension'])
                np.testing.assert_allclose(remainder, expected, atol=3e-13, rtol=1e-8)
                self.assertLessEqual(np.linalg.eigvalsh(remainder)[-1], 1e-13)
                self.assertLessEqual(np.linalg.norm(remainder, 2),
                                     row['proved_remainder_bound']+3e-13)

    def test_07_effective_capacity_counts_interior_storage(self):
        for depth in range(2, 6):
            data = case(depth)
            np.testing.assert_allclose(data['extension'].sum(axis=1), 1, atol=1e-12)
            np.testing.assert_allclose(data['mass'].sum(axis=1), len(data['full'])/3, atol=1e-10)
            self.assertGreaterEqual(np.linalg.eigvalsh(data['mass'])[0], 1-1e-12)
            self.assertAlmostEqual(float(data['mass'].sum()), len(data['full']), places=8)

    def test_08_naive_kron_heat_loses_internal_storage(self):
        for depth in (2, 3):
            data = case(depth)
            gap = np.linalg.eigvalsh(data['full'])[1]
            long_time = 40/gap
            full = heat(data['full'], long_time)[data['boundary'], data['boundary'][0]]
            naive = heat(data['kron'], long_time)[:, 0]
            np.testing.assert_allclose(full, 1/len(data['full']), atol=1e-10)
            np.testing.assert_allclose(naive, 1/3, atol=1e-10)
            self.assertGreater(float(np.max(np.abs(full-naive))), 0.2)

    def test_09_low_frequency_capacity_improves_response(self):
        for depth in range(2, 6):
            for row in frequency_rows(case(depth)):
                self.assertLess(row['effective_capacity_response_error'],
                                0.11*row['identity_capacity_response_error'])

    def test_10_empty_interior_has_no_artificial_memory(self):
        data = case(1)
        np.testing.assert_array_equal(data['kron'], data['full'])
        np.testing.assert_array_equal(data['mass'], np.eye(3))
        for time in (0, 0.5, 3):
            np.testing.assert_array_equal(kernel(data, time), np.zeros((3, 3)))
        np.testing.assert_allclose(response(data, 0.7), 0.7*np.eye(3)+data['full'])


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
    target = Path(__file__).with_name('boundary_memory_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
