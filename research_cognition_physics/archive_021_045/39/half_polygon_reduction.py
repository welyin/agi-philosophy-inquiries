"""Equality-face reduction of the polygon martingale problem to one half-circle."""

import argparse
import json
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path

import numpy as np

from outcome_updates import WINDOW_ATTENUATION as ALPHA
from polygon_martingale_certificate import source_weights, quantize_source


def pole_capacity_certificate(size, contrast, radius):
    from polygon_martingale_certificate import certified_geometry
    from certified_intervals import Interval as I, SCALE
    geo = certified_geometry(size, Fraction(contrast), Fraction(radius))
    capacity = I.rational(2, size) - geo["weights"][-1] - 2 * geo["weights"][size // 4]
    return {"size": size, "contrast": str(contrast), "radius": str(radius),
            "remaining_two_pole_capacity_interval": capacity.floats(),
            "capacity_integer_interval": [str(capacity.lo), str(capacity.hi)],
            "interval_denominator": str(SCALE), "equality_face_excluded": capacity.hi < 0}


def half_geometry(size, contrast, radius):
    if size < 32 or size % 4:
        raise ValueError("Use N divisible by four, N>=32.")
    if radius * math.cos(math.pi / size) < ALPHA:
        raise ValueError("The polygon must contain the source circle.")
    weights = source_weights(size, contrast, radius)
    quarter = size // 4
    source_indices = np.arange(-quarter + 1, quarter) % size
    target_indices = np.arange(-quarter, quarter + 1) % size
    source_angles = np.arange(-quarter + 1, quarter) * 2 * math.pi / size
    target_angles = np.arange(-quarter, quarter + 1) * 2 * math.pi / size
    source_points = radius * np.stack((np.cos(source_angles), np.sin(source_angles)), axis=1)
    beta = np.sinc(1 / size)
    targets = beta * np.stack((np.cos(target_angles), np.sin(target_angles)), axis=1)
    targets[[0, -1], 0] = 0
    source_mass = weights[source_indices]
    total = source_mass.sum()
    target_mass = np.full(len(target_indices), 1 / size)
    # Symmetric sharing of the two poles after servicing all zero-x source rows.
    pole_mass = (total - (size // 2 - 1) / size) / 2
    target_mass[[0, -1]] = pole_mass
    epsilon = 1 - (source_mass @ source_points[:, 0]) / (target_mass @ targets[:, 0])
    scaled_sources = source_points / (1 - epsilon)
    return {"size": size, "contrast": contrast, "radius": radius, "weights": weights,
            "source_indices": source_indices, "target_indices": target_indices,
            "source_mass": source_mass, "target_mass": target_mass,
            "sources": scaled_sources, "targets": targets,
            "epsilon": epsilon, "beta": beta, "total_mass": total,
            "axis_source_mass": weights[-1] + 2 * weights[quarter],
            "pole_mass_each_in_half": pole_mass,
            "axis_radius_ratio": radius / ((1 - epsilon) * beta)}


def expand_half_coupling(geometry, half_matrix):
    size = geometry["size"]
    weights, eps = geometry["weights"], geometry["epsilon"]
    source_indices = geometry["source_indices"]
    target_indices = geometry["target_indices"]
    core = np.zeros((size + 1, size))
    core[np.ix_(source_indices, target_indices)] = half_matrix
    core[np.ix_((size // 2 - source_indices) % size,
                (size // 2 - target_indices) % size)] = half_matrix
    upper, lower = size // 4, 3 * size // 4
    t = geometry["axis_radius_ratio"]
    for i, sign in ((upper, 1), (lower, -1)):
        core[i, upper] = weights[i] * (1 + sign * t) / 2
        core[i, lower] = weights[i] * (1 - sign * t) / 2
    core[-1, [upper, lower]] = weights[-1] / 2
    return (1 - eps) * core + eps * weights[:, None] / size


def solve_half_polygon(size, contrast, radius, return_half=False):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / ".research_runtime"))
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    geo = half_geometry(size, contrast, radius)
    wi, tj = geo["source_mass"], geo["target_mass"]
    sources, targets = geo["sources"], geo["targets"]
    if tj.min() <= 0 or geo["axis_radius_ratio"] > 1 or geo["epsilon"] < 0:
        raise ValueError("The elementary equality-face conditions fail.")
    count, columns_count = len(wi), len(tj)
    variables = np.arange(count * columns_count)
    ri, ci = np.indices((count, columns_count))
    mean = tj @ targets / geo["total_mass"]
    rows = np.concatenate((ri.ravel(), count + ci.ravel(), count + columns_count + ri.ravel(),
                           2 * count + columns_count + ri.ravel(), np.arange(count), count + np.arange(columns_count),
                           count + columns_count + np.arange(count), 2 * count + columns_count + np.arange(count)))
    cols = np.concatenate((variables, variables, variables, variables,
                           np.full(3 * count + columns_count, count * columns_count)))
    values = np.concatenate((np.ones(count * columns_count), np.ones(count * columns_count),
                             np.tile(targets[:, 0], count), np.tile(targets[:, 1], count),
                             wi, tj, wi * mean[0], wi * mean[1]))
    matrix = coo_matrix((values, (rows, cols)), shape=(3 * count + columns_count, count * columns_count + 1)).tocsr()
    rhs = np.concatenate((wi, tj, wi * sources[:, 0], wi * sources[:, 1]))
    objective = np.zeros(count * columns_count + 1)
    objective[-1] = -1
    result = linprog(objective, A_eq=matrix, b_eq=rhs,
                     bounds=[(0, None)] * (count * columns_count) + [(0, 1)], method="highs-ipm",
                     options={"primal_feasibility_tolerance": 1e-9, "dual_feasibility_tolerance": 1e-9,
                              "ipm_optimality_tolerance": 1e-10})
    report = {"size": size, "contrast": contrast, "source_radius": radius,
              "status": int(result.status), "message": result.message,
              "full_uniform_reserve_cap": float(geo["epsilon"]),
              "half_pole_mass_each": float(geo["pole_mass_each_in_half"]),
              "axis_radius_ratio": float(geo["axis_radius_ratio"]),
              "variables": count * columns_count + 1}
    if not result.success:
        return report, None
    delta = result.x[-1]
    half = result.x[:-1].reshape(count, columns_count) + delta * wi[:, None] * tj / geo["total_mass"]
    full = expand_half_coupling(geo, half)
    angles = np.arange(size) * 2 * math.pi / size
    unit = np.stack((np.cos(angles), np.sin(angles)), axis=1)
    full_targets = geo["beta"] * unit
    full_sources = np.vstack((radius * unit, np.zeros(2)))
    residual = max(np.max(np.abs(full.sum(axis=1) - geo["weights"])),
                   np.max(np.abs(full.sum(axis=0) - 1 / size)),
                   np.max(np.abs(full @ full_targets - geo["weights"][:, None] * full_sources)))
    report.update({"half_uniform_reserve": float(delta), "minimum_entry": float(full.min()),
                   "maximum_equality_residual": float(residual)})
    return report, half if return_half else full


class HalfPolygonReductionTests(unittest.TestCase):
    def test_two_marginals_have_the_same_mass_and_vector_mean(self):
        for n, e, r in ((128, .1425, .98992), (256, .144, .9897), (512, .1445, .989635), (2048, .14473, .98961702)):
            g = half_geometry(n, e, r)
            self.assertAlmostEqual(g["source_mass"].sum(), g["target_mass"].sum(), places=14)
            np.testing.assert_allclose(g["source_mass"] @ g["sources"],
                                       g["target_mass"] @ g["targets"], atol=2e-15)
            self.assertGreater(g["pole_mass_each_in_half"], 0)
            self.assertLess(g["axis_radius_ratio"], 1)

    def test_half_pole_budget_matches_explicit_axis_mass_formula(self):
        g = half_geometry(512, .1445, .989635)
        expected = 1 / (2 * 512) - g["weights"][128] / 2 - g["weights"][-1] / 4
        self.assertAlmostEqual(g["pole_mass_each_in_half"], expected, places=15)

    def test_exact_axis_rows_have_correct_mass_and_mean_after_expansion(self):
        g = half_geometry(512, .1445, .989635)
        half = np.zeros((len(g["source_mass"]), len(g["target_mass"])))
        full = expand_half_coupling(g, half)
        angles = np.arange(512) * 2 * math.pi / 512
        targets = g["beta"] * np.stack((np.cos(angles), np.sin(angles)), axis=1)
        rows = [128, 384, 512]
        np.testing.assert_allclose(full[rows].sum(axis=1), g["weights"][rows], atol=2e-17)
        expected = g["weights"][rows, None] * np.array([[0, .989635], [0, -.989635], [0, 0]])
        np.testing.assert_allclose(full[rows] @ targets, expected, atol=2e-17)

    def test_zero_crossing_support_is_preserved_before_uniform_reserve(self):
        g = half_geometry(128, .1425, .98992)
        half = np.ones((len(g["source_mass"]), len(g["target_mass"])))
        full = expand_half_coupling(g, half)
        core = (full - g["epsilon"] * g["weights"][:, None] / 128) / (1 - g["epsilon"])
        cosine = np.cos(np.arange(128) * 2 * math.pi / 128)
        crossing = cosine[:, None] * cosine < -1e-12
        self.assertLess(np.max(np.abs(core[:-1][crossing])), 1e-17)

    def test_negative_pole_capacity_can_occur_despite_positive_projection_budget(self):
        # Moving the auxiliary source radius too far outward consumes pole capacity.
        g = half_geometry(512, .1445, .9905)
        self.assertGreater(g["epsilon"], 0)
        self.assertLess(g["pole_mass_each_in_half"], 0)
        strict = pole_capacity_certificate(512, "0.1445", "0.9905")
        self.assertTrue(strict["equality_face_excluded"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="store_true")
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(HalfPolygonReductionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    summary = {"pole_capacity_certificates": [pole_capacity_certificate(n, e, r) for n, e, r in
                                              ((512, "0.1445", "0.989635"),
                                               (512, "0.1445", "0.9905"),
                                               (2048, "0.14473", "0.98961702"))],
               "automated_checks": {"run": checks.testsRun, "errors": 0, "failures": 0},
               "scope": "The exact equality-face equivalence is proved in note 39. The negative pole-capacity example rules out saturation of the projection reserve bound for that auxiliary radius, not all continuous couplings."}
    if args.write_results:
        Path(__file__).with_name("half_polygon_reduction_results.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    if args.probe:
        from compressed_polygon_certificate import compress_candidate, expand_candidate
        from polygon_martingale_certificate import verify_polygon_data
        rows = []
        for n, e, r in ((128, .1425, .98992), (256, .144, .9897), (512, .1445, .989635)):
            report, matrix = solve_half_polygon(n, e, r)
            if matrix is not None:
                # This is a diagnostic reconstruction; main high-strength certificate is archived separately.
                report["reconstructed_full_certificate"] = verify_polygon_data(expand_candidate(compress_candidate(report, matrix)))
            rows.append(report)
            print(json.dumps(report), flush=True)
        data = {"rows": rows, "automated_checks": {"run": checks.testsRun, "errors": 0, "failures": 0},
                "scope": "Finite probes and independently checked full reconstructions support the equality-face reduction. They do not prove a feasible sequence reaching eta_B or a continuous endpoint coupling."}
        if args.write_results:
            Path(__file__).with_name("half_polygon_probe_results.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
