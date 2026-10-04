"""Optional y-reflection reduction for the critical half-circle linear program."""

import argparse
import json
import sys
import unittest
from pathlib import Path

import numpy as np

from half_polygon_reduction import half_geometry


def expand_half_rows(representatives):
    count, columns = representatives.shape
    if columns != 2 * count + 1:
        raise ValueError("Expected N/4 representative rows and N/2+1 columns.")
    result = np.empty((2 * count - 1, columns))
    result[count - 1] = (representatives[0] + representatives[0, ::-1]) / 2
    for i in range(1, count):
        result[count - 1 + i] = representatives[i]
        result[count - 1 - i] = representatives[i, ::-1]
    return result


def solve_symmetric_half(size, contrast, radius):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / ".research_runtime"))
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix

    geo = half_geometry(size, contrast, radius)
    if geo["target_mass"].min() <= 0 or geo["axis_radius_ratio"] > 1 or geo["epsilon"] < 0:
        raise ValueError("Elementary half-face conditions fail.")
    count, columns = size // 4, size // 2 + 1
    wi = geo["source_mass"][count - 1:]
    sources, targets = geo["sources"][count - 1:], geo["targets"]
    tj = geo["target_mass"]
    mean = tj @ targets / geo["total_mass"]
    i, j = np.indices((count, columns))
    variables = np.arange(count * columns).reshape(count, columns)
    # Row mass: count. Unique nonnegative-y columns: count+1.
    # Row x moments: count. Row y moments except the fixed central row: count-1.
    column_constraint = count + np.abs(j - count)
    x_constraint = 2 * count + 1 + i
    y_constraint = 3 * count + i[1:]
    factors = np.ones((count, columns))
    factors[0] = .5
    factors[:, count] *= 2
    rows = np.concatenate((i.ravel(), column_constraint.ravel(), x_constraint.ravel(), y_constraint.ravel(),
                           np.arange(4 * count)))
    cols = np.concatenate((variables.ravel(), variables.ravel(), variables.ravel(), variables[1:].ravel(),
                           np.full(4 * count, count * columns)))
    values = np.concatenate((np.ones(count * columns), factors.ravel(),
                             np.tile(targets[:, 0], count), np.tile(targets[:, 1], count - 1),
                             wi, tj[count:], wi * mean[0], wi[1:] * mean[1]))
    matrix = coo_matrix((values, (rows, cols)), shape=(4 * count, count * columns + 1)).tocsr()
    rhs = np.concatenate((wi, tj[count:], wi * sources[:, 0], wi[1:] * sources[1:, 1]))
    objective = np.zeros(count * columns + 1)
    objective[-1] = -1
    result = linprog(objective, A_eq=matrix, b_eq=rhs,
                     bounds=[(0, None)] * (count * columns) + [(0, 1)], method="highs-ipm",
                     options={"primal_feasibility_tolerance": 1e-9, "dual_feasibility_tolerance": 1e-9,
                              "ipm_optimality_tolerance": 1e-10, "time_limit": 600})
    report = {"size": size, "contrast": contrast, "source_radius": radius,
              "status": int(result.status), "message": result.message,
              "full_uniform_reserve_cap": float(geo["epsilon"]),
              "half_pole_mass_each": float(geo["pole_mass_each_in_half"]),
              "axis_radius_ratio": float(geo["axis_radius_ratio"]),
              "variables": count * columns + 1, "equality_constraints": 4 * count,
              "solver_reduction": "Reflection in y=0; central source row averaged with its reversal."}
    if not result.success:
        return report, None
    delta = result.x[-1]
    representatives = result.x[:-1].reshape(count, columns) + delta * wi[:, None] * tj / geo["total_mass"]
    half = expand_half_rows(representatives)
    residual = max(float(np.max(np.abs(half.sum(axis=1) - geo["source_mass"]))),
                   float(np.max(np.abs(half.sum(axis=0) - tj))),
                   float(np.max(np.abs(half @ targets - geo["source_mass"][:, None] * geo["sources"]))))
    report.update(half_uniform_reserve=float(delta), minimum_half_entry=float(half.min()),
                  maximum_half_equality_residual=residual)
    return report, half


class SymmetricHalfSolverTests(unittest.TestCase):
    def test_expansion_matches_direct_pairing_and_fixed_row_average(self):
        rows = np.arange(8 * 17, dtype=float).reshape(8, 17)
        actual = expand_half_rows(rows)
        np.testing.assert_array_equal(actual, actual[::-1, ::-1])
        np.testing.assert_array_equal(actual[7], (rows[0] + rows[0, ::-1]) / 2)
        np.testing.assert_array_equal(actual[8:], rows[1:])

    def test_unique_column_weights_and_row_moments_match_expansion(self):
        count = 8
        rows = (np.arange(count * (2 * count + 1), dtype=float).reshape(count, 2 * count + 1) + 1)**2
        actual = expand_half_rows(rows)
        for j in range(count, 2 * count + 1):
            reflected = 2 * count - j
            expected = (rows[0, j] + rows[0, reflected]) / 2 + np.sum(rows[1:, j] + rows[1:, reflected])
            self.assertEqual(actual[:, j].sum(), expected)
        angle = np.arange(-count, count + 1) * np.pi / (2 * count)
        targets = np.stack((np.cos(angle), np.sin(angle)), axis=1)
        self.assertAlmostEqual(float(actual[count - 1] @ targets[:, 1]), 0, places=9)
        np.testing.assert_allclose(actual[count - 1:] @ targets[:, 0], rows @ targets[:, 0], atol=1e-10)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve-case", choices=("512", "2048"))
    parser.add_argument("--candidate-output", type=Path)
    parser.add_argument("--probe", action="store_true")
    args = parser.parse_args()
    if args.solve_case:
        if args.candidate_output is None or args.candidate_output.exists():
            parser.error("Specify a new candidate output path.")
        from half_face_certificate import CASES, pack_candidate, verify_half_candidate
        report, half = solve_symmetric_half(*CASES[args.solve_case])
        print(json.dumps(report), flush=True)
        if half is None:
            raise SystemExit(1)
        candidate = pack_candidate(report, half)
        certificate = verify_half_candidate(candidate)
        print(json.dumps(certificate), flush=True)
        if not certificate["verified"]:
            raise SystemExit(1)
        args.candidate_output.write_text(json.dumps(candidate, separators=(",", ":")) + "\n", encoding="utf-8")
        return
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SymmetricHalfSolverTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    if args.probe:
        from half_polygon_reduction import solve_half_polygon
        from half_face_certificate import pack_candidate, verify_half_candidate
        for case in ((128, .1425, .98992), (256, .144, .9897)):
            old, _ = solve_half_polygon(*case, return_half=True)
            new, half = solve_symmetric_half(*case)
            if old["status"] != 0 or half is None:
                print(json.dumps({"original": old, "reduced": new}), flush=True)
                raise SystemExit(1)
            certificate = verify_half_candidate(pack_candidate(new, half))
            difference = abs(old["half_uniform_reserve"] - new["half_uniform_reserve"])
            print(json.dumps({"original": old, "reduced": new, "reserve_difference": difference,
                              "certificate": certificate}), flush=True)
            if difference > 1e-8 or not certificate["verified"]:
                raise SystemExit(1)


if __name__ == "__main__":
    main()
