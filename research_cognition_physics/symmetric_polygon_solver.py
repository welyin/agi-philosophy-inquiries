"""Optional candidate generator using the two reflection symmetries of the problem."""

import math
import sys
import unittest
from pathlib import Path

import numpy as np

from polygon_martingale_certificate import source_weights
from outcome_updates import WINDOW_ATTENUATION as ALPHA


def orbit_data(size):
    if size < 4 or size % 4:
        raise ValueError("Use N divisible by four.")
    representatives = np.append(np.arange(size // 4 + 1), size)
    stabilizers = np.ones(len(representatives))
    stabilizers[[0, -2]] = 2
    stabilizers[-1] = 4
    return representatives, stabilizers


def reflected_indices(indices, size, element):
    operations = (indices, -indices, size // 2 - indices, indices + size // 2)
    return np.where(indices == size, size, operations[element] % size)


def expand_symmetric_rows(rows, size):
    representatives, stabilizers = orbit_data(size)
    matrix = np.zeros((size + 1, size))
    for element in range(4):
        sources = reflected_indices(representatives, size, element)
        targets = reflected_indices(np.arange(size), size, element)
        np.add.at(matrix, (sources[:, None], targets[None, :]), rows / stabilizers[:, None])
    return matrix


def solve_symmetric_polygon(size, contrast, radius):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / ".research_runtime"))
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    representatives, stabilizers = orbit_data(size)
    weights = source_weights(size, contrast, radius)
    if radius * math.cos(math.pi / size) < ALPHA or weights.min() <= 0:
        raise ValueError("Source geometry must be positive and contain the circle.")
    count = len(representatives)
    variables = np.arange(count * size).reshape(count, size)
    ri, ci = np.indices((count, size))
    angles = np.arange(size) * 2 * math.pi / size
    beta = np.sinc(1 / size)
    unit = np.stack((np.cos(angles), np.sin(angles)), axis=1)
    points = np.vstack((radius * unit, np.zeros(2)))
    targets = beta * unit
    moment_x = np.tile(targets[:, 0], (count, 1))
    moment_y = np.tile(targets[:, 1], (count, 1))
    moment_x[[-2, -1]] = 0
    moment_y[[0, -1]] = 0
    rows = [ri.ravel(), (count + size + ri).ravel(), (2 * count + size + ri).ravel()]
    columns = [variables.ravel()] * 3
    values = [np.ones(count * size), moment_x.ravel(), moment_y.ravel()]
    for element in range(4):
        rows.append((count + reflected_indices(ci, size, element)).ravel())
        columns.append(variables.ravel())
        values.append(np.repeat(1 / stabilizers, size))
    rows.append(np.concatenate((np.arange(count), count + np.arange(size))))
    columns.append(np.full(count + size, count * size))
    values.append(np.concatenate((weights[representatives], np.full(size, 1 / size))))
    matrix = coo_matrix((np.concatenate(values), (np.concatenate(rows), np.concatenate(columns))),
                        shape=(3 * count + size, count * size + 1)).tocsr()
    right = np.concatenate((weights[representatives], np.full(size, 1 / size),
                            weights[representatives] * points[representatives, 0],
                            weights[representatives] * points[representatives, 1]))
    objective = np.zeros(count * size + 1)
    objective[-1] = -1
    result = linprog(objective, A_eq=matrix, b_eq=right, bounds=[(0, None)] * (count * size) + [(0, 1)],
                     method="highs-ipm", options={"primal_feasibility_tolerance": 1e-9,
                                                 "dual_feasibility_tolerance": 1e-9,
                                                 "ipm_optimality_tolerance": 1e-10})
    report = {"size": size, "contrast": contrast, "source_radius": radius,
              "target_arc_mean_radius": float(beta), "source_center_mass": float(weights[-1]),
              "solver": "reflection-reduced HiGHS IPM; candidate only",
              "variables": count * size + 1, "status": int(result.status), "message": result.message}
    if not result.success:
        return report, None
    reserve = float(result.x[-1])
    coupling = expand_symmetric_rows(result.x[:-1].reshape(count, size), size)
    coupling += reserve * weights[:, None] / size
    residual = max(np.max(np.abs(coupling.sum(axis=0) - 1 / size)),
                   np.max(np.abs(coupling.sum(axis=1) - weights)),
                   np.max(np.abs(coupling @ targets - weights[:, None] * points)))
    report.update({"uniform_reserve": reserve, "minimum_entry": float(coupling.min()),
                   "maximum_equality_residual": float(residual)})
    return report, coupling


class SymmetricPolygonTests(unittest.TestCase):
    def test_expansion_preserves_representative_row_masses_and_reflections(self):
        size = 32
        reps, _ = orbit_data(size)
        rows = np.random.default_rng(37).uniform(size=(len(reps), size))
        matrix = expand_symmetric_rows(rows, size)
        np.testing.assert_allclose(matrix.sum(axis=1)[reps], rows.sum(axis=1), atol=1e-14)
        for element in range(4):
            ri = reflected_indices(np.arange(size + 1), size, element)
            ci = reflected_indices(np.arange(size), size, element)
            np.testing.assert_allclose(matrix, matrix[ri[:, None], ci[None, :]], atol=1e-15)

    def test_axis_and_center_stabilizers_enforce_required_moment_symmetry(self):
        size = 32
        reps, _ = orbit_data(size)
        rows = np.random.default_rng(38).uniform(size=(len(reps), size))
        matrix = expand_symmetric_rows(rows, size)
        angles = np.arange(size) * 2 * math.pi / size
        means = matrix @ np.stack((np.cos(angles), np.sin(angles)), axis=1)
        self.assertAlmostEqual(means[0, 1], 0, places=13)
        self.assertAlmostEqual(means[size // 4, 0], 0, places=13)
        np.testing.assert_allclose(means[-1], 0, atol=1e-14)
