"""Independent small certificates for 1051; no author function import, no new group."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_checks_results.json"


def null_vector(matrix):
    a = [[F(x) for x in row] for row in matrix]
    rows, cols = len(a), len(a[0])
    pivots, r = [], 0
    for c in range(cols):
        found = next((i for i in range(r, rows) if a[i][c]), None)
        if found is None:
            continue
        a[r], a[found] = a[found], a[r]
        scale = a[r][c]
        a[r] = [x / scale for x in a[r]]
        for i in range(rows):
            if i != r:
                scale = a[i][c]
                a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = next(c for c in range(cols) if c not in pivots)
    v = [F(0)] * cols
    v[free] = F(1)
    for i, c in enumerate(pivots):
        v[c] = -a[i][free]
    denom = math.lcm(*(x.denominator for x in v))
    z = [int(x * denom) for x in v]
    divisor = math.gcd(*z)
    return [x // divisor for x in z]


def run():
    receipt_path = HERE / "research_round_1051_checks.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    for name, digest in receipt["owned_sha256"].items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
    rows = []
    for A in [[[1, 1, 1]], [[2, 3, 5]], [[1, 2, -3, 1], [0, 1, 1, 2]]]:
        v = null_vector(A)
        assert all(sum(a * b for a, b in zip(row, v)) == 0 for row in A)
        j = next(i for i, x in enumerate(v) if x)
        M = len(v)
        pj = min(F(1, 2), F(1, 2**(j + 2)))
        pM = min(F(1, 2), F(1, 2**(M + 1)))
        assert pj * (1 - pj) >= pM * (1 - pM)
        # |0> and |j> entangled with a two-dimensional passive reference.
        psi = np.array([math.sqrt(1 - float(pj)), 0, 0, math.sqrt(float(pj))])
        rho = np.outer(psi, psi)
        theta = np.pi / v[j]
        unitary = np.diag([1, 1, np.exp(1j * theta * v[j]), np.exp(1j * theta * v[j])])
        opposite = unitary @ rho @ unitary.conj().T
        distance = float(np.linalg.svd(rho - opposite, compute_uv=False).sum() / 2)
        # Exact finite roots-of-unity mean agrees with the relevant Haar block.
        order = 2 * max(abs(x) for x in v) + 1
        averaged = np.zeros((4, 4), complex)
        for k in range(order):
            phase = np.exp(2j * np.pi * k * v[j] / order)
            u = np.diag([1, 1, phase, phase])
            averaged += u @ rho @ u.conj().T / order
        loss = float(np.linalg.svd(rho - averaged, compute_uv=False).sum() / 2)
        assert abs(distance**2 - float(4 * pj * (1 - pj))) < 1e-12
        assert abs(loss**2 - float(pj * (1 - pj))) < 1e-12
        assert np.linalg.norm(averaged - np.diag(np.diag(rho))) < 1e-12
        rows.append({"integer_character_matrix": A, "primitive_kernel_circle": v,
                     "nontrivial_coordinate_one_based": j + 1,
                     "budget": "1/2", "p_j": str(pj),
                     "orbit_distance_squared_exact": str(4 * pj * (1 - pj)),
                     "average_error_squared_exact": str(pj * (1 - pj)),
                     "dimension_threshold_squared_exact": str(pM * (1 - pM)),
                     "numerical_average_loss": loss})
    # Independent exact source divergence coefficients; no truncated vector routine.
    delta = F(1, 100)
    costs = [2 * delta**2 * (1 - delta**2) * n / (1 - F(1, 2)**n)
             for n in [4, 16, 64]]
    assert costs[0] < costs[1] < costs[2]
    assert all(x >= 2 * delta**2 * (1 - delta**2) * n
               for x, n in zip(costs, [4, 16, 64]))
    return {"round": 1051, "reviewer": "next_structural_selection",
            "all_passed": True, "additional_science_groups": 0,
            "no_author_science_functions_imported": True,
            "receipt_sha256": hashlib.sha256(receipt_path.read_bytes()).hexdigest(),
            "owned_sha256": receipt["owned_sha256"],
            "kernel_circle_and_reference_certificates": rows,
            "source_costs_exact_for_lengths_4_16_64": [str(x) for x in costs],
            "source_distance_squared_exact": str(delta**2 * (1 - delta**2)),
            "general_quotient_dimension_and_infinite_source_proved_analytically": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        old = json.loads(OUT.read_text(encoding="utf-8"))
        assert result == old
    print(json.dumps({"round": 1051, "independent_checks_passed": True,
                      "additional_science_groups": 0,
                      "mode": "exclusive_write" if args.write else "read_only"}))


if __name__ == "__main__":
    main()
