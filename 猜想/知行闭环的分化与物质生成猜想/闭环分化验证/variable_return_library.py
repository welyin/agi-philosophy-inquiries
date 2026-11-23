"""Exact optimal shape and finite preparation libraries for delayed readout."""

import argparse
import json
from math import atan2, cos, hypot, isclose, pi, sin, sqrt
from pathlib import Path

from late_branch_memory import optimal_qubit


def library_scores(angle):
    # Readout may depend on the grid point, the question, and the returned
    # bit. No measurement depends on the relative frequency of those cases.
    a = sin(angle) / sqrt(2)
    c = cos(angle)
    without_return = (1 + a) / 2
    with_return = (3 + hypot(a, c)) / 4
    return without_return, with_return


def run():
    samples = 601
    geometry = [optimal_qubit(i / (samples - 1)) for i in range(samples)]
    diagonal = [c * c - 2 * a * a for i, (a, c, _) in enumerate(geometry)
                if 2 * (samples - 1) < 3 * i and i < samples - 1]
    assert all(later > earlier for earlier, later in zip(diagonal, diagonal[1:]))
    data = {}
    for intervals in (2, 4, 8, 16):
        grid = [j * pi / (2 * intervals) for j in range(intervals + 1)]
        readers = [library_scores(angle) for angle in grid]
        bound = sin(pi / (8 * intervals))
        # The score's first derivative vanishes at the interior optimum.
        # In phi=asin(a), the loss is at most (1-cos(delta_phi))/4;
        # |d phi / d angle| <= 1/sqrt(2). Boundary optima are on the grid.
        quadratic_bound = pi ** 2 / (256 * intervals ** 2)
        worst = 0.
        for i, (a, c, optimum) in enumerate(geometry):
            p = i / (samples - 1)
            angle = atan2(sqrt(2) * a, c)
            index = min(range(len(grid)), key=lambda j: abs(grid[j] - angle))
            without_return, with_return = readers[index]
            loss = optimum - ((1 - p) * without_return + p * with_return)
            # The corresponding pure states have Bloch angle |angle-grid[j]|
            # for each input word, including the flipped sign of the parity axis.
            distance = hypot(sqrt(2) * a - sin(grid[index]), c - cos(grid[index])) / 2
            assert isclose(distance, sin(abs(angle - grid[index]) / 2), abs_tol=1e-12)
            assert -1e-12 <= loss <= distance + 1e-12 and distance <= bound + 1e-12
            assert loss <= quadratic_bound + 1e-12
            worst = max(worst, loss)
        data[str(intervals)] = {"observed_max_loss": round(worst, 9),
                                "uniform_loss_bound": round(bound, 9),
                                "quadratic_uniform_loss_bound": round(quadratic_bound, 9)}
    return {"probability_samples": samples, "strictly_ordered_interior_pairs": len(diagonal) - 1,
            "libraries_by_angular_intervals": data}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    options = parser.parse_args()
    result = run()
    target = Path(__file__).with_name("variable_return_library_results.json")
    if options.check:
        assert json.loads(target.read_text(encoding="utf-8")) == result
    else:
        with target.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
