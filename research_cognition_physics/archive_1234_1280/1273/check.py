"""Round 1273. Check old carrier algebra and new additive-response diagnostic.

Default is read-only and compares a saved result. --write is for first creation.
No new pressure material, geometry evolution or gravitational source is simulated.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "1272" / "protocol_check.py"
spec = importlib.util.spec_from_file_location("r1272_reused", SOURCE)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def rank_exact(rows):
    a = [[F(v) for v in row] for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        divisor = a[rank][col]
        a[rank] = [v / divisor for v in a[rank]]
        for i in range(len(a)):
            if i != rank:
                factor = a[i][col]
                a[i] = [x - factor * y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def compute():
    energy = old.ENERGY
    selected = energy == 2
    assert selected.sum() == 5760
    gate_checks = []
    for gate in old.GATES:
        perm = gate["permutation"]
        assert np.array_equal(energy[perm], energy)
        # A nonlinear f is only an implementation check; the proof covers any f.
        fenergy = energy * energy + 3 * energy + 7
        assert np.array_equal(fenergy[perm], fenergy)
        assert np.all(fenergy[selected] == 17)
        gate_checks.append({"name": gate["name"], "energy_invariant": True})

    incidence = [[1, 0, -1, 0], [0, 1, 0, -1], [1, -1, 0, 0]]
    assert rank_exact(incidence) == 3
    assert all(sum(row) == 0 for row in incidence)
    # Four basis coefficients, checked on every old active basis entry.
    observed = []
    for j in (0, 1, 2):
        gate = old.GATES[j]
        differences = (old.COORD[:4, gate["permutation"]] - old.COORD[:4]).T
        distinct = sorted(set(map(tuple, differences[selected].tolist())))
        nonzero = [v for v in distinct if any(v)]
        assert len(nonzero) == 2
        assert tuple(incidence[j]) in nonzero
        assert tuple(-v for v in incidence[j]) in nonzero
        observed.append({"gate": gate["name"], "nonzero_rows": [list(v) for v in nonzero]})

    weights = ((1, 1, 1, 1), (1, 0, 0, 0), (1, -1, 0, 0), (2, 3, 5, 7))
    diagnostics = []
    max_error = 0.0
    for k in (0, 1):
        for ell in (0, 1, 2):
            final, snapshots = old.run(old.ready(k, ell))
            before, after = snapshots[2], snapshots[3]
            assert np.sum(np.abs(final[~selected]) ** 2) == 0
            for weights_i in weights:
                c = np.asarray(weights_i) @ old.COORD[:4]
                delta = float(np.sum(c[:, None] * (np.abs(after) ** 2 - np.abs(before) ** 2)))
                exact = F(3, 40) * (weights_i[0] - weights_i[1]) if k else F(0)
                error = abs(delta - float(exact))
                max_error = max(max_error, error)
                assert error < 5e-14
                diagnostics.append({"K": k, "L": ell, "weights": list(weights_i),
                                    "exact_delta_C": str(exact)})
    # These dimensionless choices illustrate the theorem's budget arithmetic,
    # not a certified new physical interaction W.
    d1 = d2 = F(1, 10**9)
    eta1, eta2 = F(1, 10000), F(1, 20000)
    response_bound = d1 + d2 + eta1 + eta2
    assert response_bound == F(75001, 500000000)
    return {
        "round": 1273,
        "status": "PASS",
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "carrier_dimension": old.DIM,
        "fixed_energy_sector_dimension": int(selected.sum()),
        "gate_checks": gate_checks,
        "exchange_incidence_rank": rank_exact(incidence),
        "linear_response_kernel_dimension": 1,
        "permutation_response_rows": observed,
        "old_protocol_diagnostics": diagnostics,
        "diagnostic_numerical_error_below": 5e-14,
        "finite_bound_arithmetic_only": str(response_bound),
        "new_coupled_geometry_prediction_claimed": False,
        "analytic_theorem_not_replaced_by_sampling": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = compute()
    target = HERE / "results.json"
    if args.write:
        assert not target.exists(), "Do not overwrite accepted results"
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    else:
        assert result == json.loads(target.read_text(encoding="utf-8"))
    print(json.dumps({"round": 1273, "status": "PASS", "gates": 9,
                      "additive_diagnostics": len(result["old_protocol_diagnostics"]),
                      "source_unchanged": True, "mode": "write" if args.write else "read_only"}))
