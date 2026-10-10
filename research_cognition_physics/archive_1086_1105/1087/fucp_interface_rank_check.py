"""Read-only finite calibration of the specific tent-interface rank obstruction."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import numpy as np


def determinant(matrix):
    a = [row[:] for row in matrix]
    value = F(1)
    for k in range(len(a)):
        pivot = next(i for i in range(k, len(a)) if a[i][k])
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            value = -value
        diagonal = a[k][k]
        value *= diagonal
        for i in range(k + 1, len(a)):
            ratio = a[i][k] / diagonal
            for j in range(k + 1, len(a)):
                a[i][j] -= ratio * a[k][j]
            a[i][k] = F(0)
    return value


def calculate():
    exact = []
    for n in (2, 3, 5, 9):
        length = F(1, 2)
        h = length / (n - 1)
        matrix = [[1 - h * abs(i - j) for j in range(n)] for i in range(n)]
        actual = determinant(matrix)
        expected = (2 * h) ** (n - 1) * (1 - length / 2)
        assert actual == expected > 0
        threshold = h * (2 - length) / (4 * n * (1 + 2 * h))
        assert threshold == F(3, 16 * n * n)
        exact.append({'N': n, 'determinant': str(actual), 'error_threshold': str(threshold)})
    rng = np.random.default_rng(1087)
    numerical = []
    for n in (2, 4, 8, 16, 32):
        h = 0.5 / (n - 1)
        t = 1 - h * np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        d = np.eye(n) - np.eye(n, k=-1)
        target = np.eye(n) * (2 * h)
        target[0, 0] = 1
        target[0, 1:] = target[1:, 0] = -h
        assert np.max(np.abs(d @ t @ d.T - target)) < 2e-14
        floor = 3 / (16 * n)
        eigmin = float(np.linalg.eigvalsh(t)[0])
        assert eigmin >= floor - 2e-14
        eta = 3 / (32 * n * n)
        noise = rng.uniform(-eta, eta, size=(n, n))
        singular_min = float(np.linalg.svd(t + noise, compute_uv=False)[-1])
        assert singular_min >= floor - n * eta - 2e-14 > 0
        numerical.append({'N': n, 'lambda_min': eigmin, 'analytic_lower_bound': floor,
                          'eta': eta, 'perturbed_sigma_min': singular_min})
    return {'round': 1087, 'passed': True, 'scope': 'specific exact tent interface on one fixed finite complete object',
            'exact_cases': exact, 'numerical_cases': numerical,
            'fixed_positive_error_requires_infinite_dimension': False,
            'FUCP_plus_CO_inconsistent_proved': False,
            'joint_quantum_clock_signal_realization_supplied': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = calculate()
    if args.check:
        saved = json.loads(Path(__file__).with_name('fucp_interface_rank_results.json').read_text('utf8'))
        assert saved.keys() == data.keys()
        for key in data.keys() - {'numerical_cases'}:
            assert saved[key] == data[key], key
        assert len(saved['numerical_cases']) == len(data['numerical_cases'])
        for old, new in zip(saved['numerical_cases'], data['numerical_cases']):
            assert old.keys() == new.keys()
            for key in new:
                assert abs(old[key] - new[key]) < 1e-12, (key, old, new)
    print(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
