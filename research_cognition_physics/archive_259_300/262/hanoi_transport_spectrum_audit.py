"""Round 262: audit known Hanoi spectra, resistance and transport scaling.

Unit conductances and the classical generator -L are added dynamical inputs.
Sources: Grigorchuk-Sunic arXiv:0711.0068 Theorem 1.1;
Alekseyev-Berger arXiv:1304.3780 Section 4. No quantum Hamiltonian is derived.
"""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
from fisher_hanoi_scaling_audit import hanoi


@lru_cache(maxsize=None)
def laplacian(n):
    graph, _ = hanoi(n)
    matrix = np.zeros((len(graph), len(graph)))
    for v, row in enumerate(graph):
        matrix[v, v] = len(row)
        matrix[v, list(row)] = -1
    return matrix


def inverse_minus(x):
    return 2*x/(5+math.sqrt(25-4*x))


def preimages(values):
    return tuple(y for x in values for y in (inverse_minus(x), 5-inverse_minus(x)))


@lru_cache(maxsize=None)
def spectral_atoms(n):
    """(eigenvalue, multiplicity), translating published adjacency spectrum."""
    atoms = [(0., 1)]
    for seed, offset, constant in ((3., 1, 3), (5., 2, -1)):
        values = (seed,)
        for depth in range(n-offset+1):
            multiplicity = (3**(n-depth-1)+constant)//2
            atoms.extend((value, multiplicity) for value in values)
            values = preimages(values)
    return tuple(sorted(atoms))


def gap_formula(n):
    gap = 3.
    for _ in range(n-1):
        gap = inverse_minus(gap)
    return gap


@lru_cache(maxsize=None)
def numerical_spectrum(n):
    return np.linalg.eigvalsh(laplacian(n))


@lru_cache(maxsize=None)
def spectrum_metrics(n):
    values = numerical_spectrum(n)
    predicted = np.repeat([x for x, _ in spectral_atoms(n)],
                          [m for _, m in spectral_atoms(n)])
    return {'n': n, 'vertices': 3**n, 'gap_direct': float(values[1]),
            'gap_formula': gap_formula(n),
            'full_spectrum_max_abs_error': float(np.max(np.abs(values-predicted)))}


@lru_cache(maxsize=None)
def transport_metrics(n):
    graph, words = hanoi(n)
    source, target = 0, words.index((1,)*n)
    keep = [i for i in range(len(graph)) if i != target]
    grounded = laplacian(n)[np.ix_(keep, keep)]
    b = np.zeros(len(keep))
    b[keep.index(source)] = 1
    voltage = np.linalg.solve(grounded, b)
    resistance = float(voltage[keep.index(source)])
    # For the discrete walk choosing uniformly among actual neighbors,
    # (D-A) h = degree, with h[target] = 0.
    degree = np.array([len(graph[i]) for i in keep], dtype=float)
    hitting = float(np.linalg.solve(grounded, degree)[keep.index(source)])
    predicted_r = (5/3)**n-1
    m = 3*(3**n-1)//2
    return {'n': n, 'corner_resistance_direct': resistance,
            'corner_resistance_formula': predicted_r,
            'discrete_corner_hitting_direct': hitting,
            'discrete_corner_hitting_formula': m*predicted_r}


def dos_rows(n=12):
    atoms = spectral_atoms(n)
    threshold = 2.
    rows = []
    for k in range(7):
        count = sum(m for x, m in atoms if 0 < x <= threshold)
        rows.append({'k': k, 'threshold': threshold, 'positive_eigenvalues_below': count,
                     'finite_fraction': count/3**n, 'limiting_fraction': 3.**(-k-1)})
        threshold = inverse_minus(threshold)
    return rows


def report():
    return {'round': 262, 'direct_spectra': [spectrum_metrics(n) for n in range(1, 7)],
            'direct_transport': [transport_metrics(n) for n in range(1, 6)],
            'dos_level': 12, 'dos_thresholds': dos_rows(),
            'gap_ratio_at_level_20': gap_formula(19)/gap_formula(20),
            'relaxation_exponent_log5_over_log2': math.log(5)/math.log(2),
            'dos_spectral_exponent_2log3_over_log5': 2*math.log(3)/math.log(5),
            'scope': 'Exact published Hanoi graph spectral and resistance formulae are checked on the round-261 graph family. Unit conductances, diffusion generator -L and the discrete uniform-neighbor walk are explicit distinct clock conventions. These exponents do not derive physical time, light speed, a quantum Hamiltonian or mixed-growth universality.'}


class Audit(unittest.TestCase):
    def test_01_loop_convention_matches_schreier_operator(self):
        for n in range(1, 7):
            graph, _ = hanoi(n)
            degree = np.array(list(map(len, graph)))
            adjacency = np.diag(degree)-laplacian(n)
            regularized = adjacency+np.diag(3-degree)
            np.testing.assert_array_equal(regularized.sum(axis=1), 3*np.ones(3**n))
            np.testing.assert_array_equal(3*np.eye(3**n)-regularized, laplacian(n))

    def test_02_complete_spectrum_with_multiplicities(self):
        for n in range(1, 7):
            self.assertLess(spectrum_metrics(n)['full_spectrum_max_abs_error'], 1e-10)

    def test_03_gap_decimation_against_direct_eigensolve(self):
        for n in range(1, 7):
            row = spectrum_metrics(n)
            self.assertAlmostEqual(row['gap_direct'], row['gap_formula'], delta=1e-11)
        for n in range(2, 30):
            gap = gap_formula(n)
            self.assertAlmostEqual(gap*(5-gap)/gap_formula(n-1), 1., places=13)

    def test_04_resistance_against_grounded_network(self):
        for n in range(1, 6):
            row = transport_metrics(n)
            self.assertAlmostEqual(row['corner_resistance_direct']/row['corner_resistance_formula'],
                                   1., delta=1e-11)

    def test_05_hitting_equation_and_commute_identity(self):
        for n in range(1, 6):
            row = transport_metrics(n)
            self.assertAlmostEqual(row['discrete_corner_hitting_direct']/row['discrete_corner_hitting_formula'],
                                   1., delta=1e-11)

    def test_06_spectrum_dimension_trace_and_zero_mode(self):
        for n in range(1, 13):
            atoms = spectral_atoms(n)
            self.assertEqual(sum(m for _, m in atoms), 3**n)
            self.assertEqual([(x, m) for x, m in atoms if x == 0], [(0., 1)])
            self.assertAlmostEqual(sum(x*m for x, m in atoms)/(3*(3**n-1)), 1., places=12)

    def test_07_dos_thresholds_against_dense_spectrum(self):
        for n in range(2, 7):
            threshold = 2.
            for k in range(n-1):
                actual = int(np.count_nonzero((numerical_spectrum(n) > 1e-10) &
                                             (numerical_spectrum(n) <= threshold)))
                self.assertEqual(actual, 3**(n-k-1)-1)
                threshold = inverse_minus(threshold)

    def test_08_large_dos_counts_and_limiting_measure(self):
        for row in dos_rows():
            self.assertEqual(row['positive_eigenvalues_below'], 3**(11-row['k'])-1)
            self.assertAlmostEqual(row['limiting_fraction']-row['finite_fraction'], 3.**-12)
        # Sum the two seed families through depth 11; missing tail is (2/3)^12.
        mass = sum(2*2**i/(6*3**i) for i in range(12))
        self.assertAlmostEqual(mass+(2/3)**12, 1.)

    def test_09_relaxation_scaling(self):
        self.assertAlmostEqual(gap_formula(19)/gap_formula(20), 5., places=12)
        self.assertAlmostEqual(math.log(3)/math.log(2)/(math.log(5)/math.log(2)),
                               math.log(3)/math.log(5), places=14)

    def test_10_geometry_alone_does_not_fix_clock(self):
        matrix = laplacian(3)
        for rate in (.1, 7.):
            scaled = np.linalg.eigvalsh(rate*matrix)
            self.assertAlmostEqual(float(scaled[1])/gap_formula(3), rate, places=11)
            np.testing.assert_array_equal(matrix != 0, rate*matrix != 0)


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
    target = Path(__file__).with_name('hanoi_transport_spectrum_audit_results.json')
    if args.write_results:
        if target.exists() and target.read_text(encoding='utf-8') != payload:
            raise RuntimeError('Preserve the existing scientific result.')
        target.write_text(payload, encoding='utf-8')
    print(payload)
