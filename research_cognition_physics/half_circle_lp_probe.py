"""Exploratory finite martingale programs. No discretization is an exact model."""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

from compensated_angle_kernel import preconvolution_density
from convex_order_interface import absolute_threshold
from outcome_updates import WINDOW_ATTENUATION as ALPHA


def solve_circle(size, contrast):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / ".research_runtime"))
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    angles = (np.arange(size) + 0.5) * 2 * math.pi / size
    vectors = np.stack((np.cos(angles), np.sin(angles)), axis=-1)
    source = preconvolution_density(angles, contrast)
    source /= source.sum()
    target = np.full(size, 1 / size)
    incoming, outgoing = np.indices((size, size))
    variables = np.arange(size * size)
    rows = np.concatenate((incoming.ravel(), size + outgoing.ravel(),
                           2 * size + incoming.ravel(), 3 * size + incoming.ravel(),
                           2 * size + np.arange(size), 3 * size + np.arange(size)))
    columns = np.concatenate((variables, variables, variables, variables,
                              np.full(size, size * size), np.full(size, size * size)))
    values = np.concatenate((np.ones(size * size), np.ones(size * size),
                             np.tile(vectors[:, 0], size), np.tile(vectors[:, 1], size),
                             -source * vectors[:, 0], -source * vectors[:, 1]))
    matrix = coo_matrix((values, (rows, columns)), shape=(4 * size, size * size + 1)).tocsr()
    objective = np.zeros(size * size + 1)
    objective[-1] = -1
    right = np.concatenate((source, target, np.zeros(2 * size)))
    result = linprog(objective, A_eq=matrix, b_eq=right,
                     bounds=[(0, None)] * (size * size) + [(0, 1)], method="highs",
                     options={"dual_feasibility_tolerance": 1e-9, "primal_feasibility_tolerance": 1e-9})
    record = {"size": size, "contrast": contrast, "status": int(result.status), "message": result.message}
    if result.success:
        bound = float(result.x[-1])
        homogeneous_bound = float(np.mean(np.abs(vectors[:, 0])) / (source @ np.abs(vectors[:, 0])))
        record.update({"maximum_radius_discrete": bound, "actual_alpha": ALPHA,
                       "discrete_absolute_projection_bound": homogeneous_bound,
                       "maximum_equality_residual": float(np.max(np.abs(matrix @ result.x - right))),
                       "positive_coupling_entries": int(np.sum(result.x[:-1] > 1e-10))})
        dual = result.eqlin.marginals
        planes = np.stack((dual[2 * size:3 * size], dual[3 * size:], dual[:size]), axis=-1)
        record["dual_planes"] = planes.tolist()
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[32, 64, 128])
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    eta = sum(absolute_threshold()) / 2
    rows = []
    for size in args.sizes:
        result = solve_circle(size, eta)
        rows.append(result)
        print(json.dumps({key: value for key, value in result.items() if key != "dual_planes"}), flush=True)
    if args.write_results:
        Path(__file__).with_name("half_circle_lp_probe_results.json").write_text(
            json.dumps({"scope": "Finite grid exploration only, not a continuum feasibility proof.", "rows": rows}, indent=2) + "\n",
            encoding="utf-8")


if __name__ == "__main__":
    main()
