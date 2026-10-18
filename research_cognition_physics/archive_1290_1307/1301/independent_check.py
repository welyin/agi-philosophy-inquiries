"""Independent exact flow and matrix diagnostics for round 1301.

No author code/results are imported. Unbounded integer-charge exclusion is
certified by exact polynomial coefficients, not a charge-range search.
Default and --check are read-only; --write exclusively creates this result.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from functools import reduce
import hashlib
from itertools import product
import json
from math import gcd
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"


def positive_compositions(total):
    if total == 0:
        yield ()
    else:
        for first in range(1, total + 1):
            for rest in positive_compositions(total - first):
                yield (first,) + rest


def all_flow_shapes(total):
    """Positive compositions, with at most one zero inserted in each gap."""
    answer = set()
    for positive in positive_compositions(total):
        for gaps in product((False, True), repeat=len(positive) - 1):
            out = []
            for i, value in enumerate(positive):
                out.append(value)
                if i < len(gaps) and gaps[i]:
                    out.append(0)
            flow = tuple(out)
            assert flow[0] > 0 and flow[-1] > 0
            assert not any(a == b == 0 for a, b in zip(flow, flow[1:]))
            assert len(flow) <= 2 * total - 1
            answer.add(flow)
    assert len(answer) == 3 ** (total - 1)
    return sorted(answer, key=lambda x: (len(x), x))


def multiplicities(flow):
    return [a + b for a, b in zip((0,) + tuple(flow), tuple(flow) + (0,))]


def cubic_polynomial(flow):
    # A_1 = sum (-1)^i*x_i = 0 cancels the k^2 and common terms.
    assert sum((-1) ** i * x for i, x in enumerate(flow)) == 0
    m1 = sum((-1) ** i * i * x for i, x in enumerate(flow))
    m2 = sum((-1) ** i * i * i * x for i, x in enumerate(flow))
    return 6 * m1, 3 * (m1 + m2)


def exact_lower_bound():
    output = {}
    expected = {
        2: {(1, 1): (-6, -6)},
        4: {
            (2, 2): (-12, -12),
            (1, 2, 1): (0, 6),
            (1, 0, 1, 2): (-24, -54),
            (1, 1, 1, 1): (-12, -24),
            (2, 1, 0, 1): (-24, -42),
            (1, 1, 0, 1, 1): (0, 18),
            (1, 0, 1, 1, 0, 1): (-36, -108),
        },
    }
    for total in (2, 4):
        rows = []
        certified = {}
        for flow in all_flow_shapes(total):
            if sum((-1) ** i * x for i, x in enumerate(flow)):
                continue
            slope, intercept = cubic_polynomial(flow)
            if slope:
                root = Fraction(-intercept, slope)
                assert root < 1
                no_root_reason = "unique root " + str(root) + " is below k >= 1"
            else:
                assert intercept != 0
                no_root_reason = "nonzero constant"
            # Diagnostic cross-check of the exact coefficient identity;
            # these sample k do not carry the all-k exclusion.
            for k in (1, 2, 17, 1000003):
                charges = [(-1) ** (i + 1) * (k + i)
                           for i in range(len(flow) + 1)]
                n = multiplicities(flow)
                assert sum(q * count for q, count in zip(charges, n)) == 0
                assert sum(q ** 3 * count for q, count in zip(charges, n)) == slope*k+intercept
            certified[flow] = (slope, intercept)
            rows.append({"flow": list(flow), "cubic_slope": slope,
                         "cubic_intercept": intercept, "exclusion": no_root_reason})
        assert certified == expected[total]
        output[str(2 * total)] = {
            "all_compositions_with_optional_gap_zero": 3 ** (total - 1),
            "linear_anomaly_zero_shapes": len(rows),
            "all_integer_k_ge_1_excluded": True, "certificate": rows}
    return output


def path_matching(charges, flow):
    n = multiplicities(flow)
    groups, offset = [], 0
    for count in n:
        groups.append(list(range(offset, offset + count)))
        offset += count
    remaining = [group.copy() for group in groups]
    expanded = [q for q, count in zip(charges, n) for _ in range(count)]
    pairs = []
    for edge, count in enumerate(flow):
        for _ in range(count):
            left, right = remaining[edge].pop(0), remaining[edge+1].pop(0)
            pairs.append((left, right))
    assert not any(remaining)
    return expanded, pairs


def graph_components(values, h):
    unseen = set(values)
    components = []
    while unseen:
        seen, frontier = set(), [min(unseen)]
        while frontier:
            q = frontier.pop()
            if q in seen:
                continue
            seen.add(q)
            frontier += [r for r in unseen if r != q and abs(q+r) == h and r not in seen]
        unseen -= seen
        components.append(sorted(seen))
    return sorted(components)


def inspect_matrix(name, charges, pairs, h, connected):
    size = len(charges)
    assert size % 2 == 0
    assert 0 not in charges
    assert not set(charges).intersection(-q for q in charges)
    assert reduce(gcd, [abs(q) for q in charges] + [h]) == 1
    assert sum(charges) == 0
    assert sum(q ** 3 for q in charges) == 0
    used = [i for pair in pairs for i in pair]
    assert sorted(used) == list(range(size))
    base = np.zeros((size, size))
    perturbation = np.zeros_like(base)
    for a, b in pairs:
        assert abs(charges[a]+charges[b]) == h
        base[a, b] = base[b, a] = 1.
    assert np.array_equal(base @ base, np.eye(size))
    for a in range(size):
        for b in range(a+1, size):
            if abs(charges[a]+charges[b]) == h:
                # Independent positive rational coefficients; <= 7/8.
                perturbation[a, b] = perturbation[b, a] = (1+(3*a+5*b) % 7)/8
    epsilon = Fraction(1, 32*size)
    actual = base + float(epsilon)*perturbation
    error_bound = float(epsilon)*np.linalg.norm(perturbation, "fro")
    assert error_bound < 1/32
    singular = np.linalg.svd(actual, compute_uv=False)
    assert min(singular) >= 31/32 - 1e-12
    assert np.allclose(actual, actual.T, atol=0, rtol=0)
    for a in range(size):
        for b in range(size):
            allowed = a != b and abs(charges[a]+charges[b]) == h
            assert bool(actual[a, b] != 0) == allowed
    components = graph_components(charges, h)
    assert (len(components) == 1) == connected
    return {
        "name": name, "Weyl_number": size, "h": h, "charges_with_multiplicity": charges,
        "linear_anomaly": sum(charges), "cubic_anomaly": sum(q**3 for q in charges),
        "complete_chirality": True, "primitive_combined_integer_charge_lattice": True,
        "Higgs_reference_stabilizer_order": h,
        "actual_role_components": components, "all_allowed_blocks_nonzero": True,
        "epsilon": str(epsilon), "uniform_analytic_singular_value_lower_bound": "31/32",
        "computed_smallest_singular_value": round(float(min(singular)), 12),
        "computed_frobenius_perturbation_bound": round(float(error_bound), 12),
        "mass_rank": int(np.linalg.matrix_rank(actual)),
        "finite_matrix_not_actual_instrument": True}


TEMPLATES = {
    12: ([2, -3, 4, -5, 6, -7], [2, 2, 0, 1, 1]),
    16: ([-1, 2, -3, 4, -5], [1, 3, 3, 1]),
    20: ([-1, 2, -3, 4], [3, 5, 2]),
}


def templates_and_deletions():
    output = {}
    for size, (charges, flow) in TEMPLATES.items():
        expanded, pairs = path_matching(charges, flow)
        output[str(size)] = inspect_matrix("connected_template_"+str(size), expanded, pairs, 1, True)
        output[str(size)]["distinct_path"] = charges
        output[str(size)]["matching_flow"] = flow
        assert len(expanded) == size
    charges = [2, -1, 10, -9, 4, -5, 8, -9]
    output["delete_connectedness"] = inspect_matrix(
        "eight_disconnected_h1", charges, [(0,1),(2,3),(4,5),(6,7)], 1, False)
    expanded, pairs = path_matching([9,-5,1,3,-7], [1,0,1,2])
    output["delete_trivial_stabilizer"] = inspect_matrix("eight_connected_h4", expanded, pairs, 4, True)
    assert len(expanded) == 8
    return output


def decomposition(p):
    assert p >= 3
    if p % 3 == 0:
        counts = {12: p//3, 16: 0, 20: 0}
    elif p % 3 == 1:
        counts = {12: (p-4)//3, 16: 1, 20: 0}
    else:
        counts = {12: (p-5)//3, 16: 0, 20: 1}
    assert all(n >= 0 for n in counts.values())
    assert sum(size*n for size,n in counts.items()) == 4*p
    return counts


def all_number_construction():
    common_branch = {-1,2,-3,4,-5,6,-7}
    supports = [set(charges) for charges, _ in TEMPLATES.values()]
    assert all(s <= common_branch for s in supports)
    assert not common_branch.intersection(-q for q in common_branch)
    assert all(a & b for a in supports for b in supports)
    # Only finite diagnostics of the symbolic mod-3 construction.
    rows = []
    for p in (3,4,5,6,7,8,11,16):
        charges, pairs = [], []
        counts = decomposition(p)
        for size, count in counts.items():
            base_charges, flow = TEMPLATES[size]
            expanded, base_pairs = path_matching(base_charges, flow)
            for _ in range(count):
                offset = len(charges)
                charges += expanded
                pairs += [(a+offset,b+offset) for a,b in base_pairs]
        row = inspect_matrix("4p_"+str(p), charges, pairs, 1, True)
        rows.append({"p": p, "template_copies": {str(size): count for size, count in counts.items()},
                     "Weyl_number": row["Weyl_number"],
                     "mass_rank": row["mass_rank"],
                     "smallest_singular_value": row["computed_smallest_singular_value"],
                     "actual_role_connected": len(row["actual_role_components"]) == 1})
    return {"symbolic_scope": "all N=4p with integer p>=3",
            "three_residue_formulas": ["p=3a", "p=3a+4", "p=3a+5"],
            "supports_pairwise_overlap": True, "no_opposite_charges_in_common_branch": True,
            "finite_diagnostics": rows}


def payload():
    admission = HERE / "../_shared/round1301_admission.md"
    return {
        "status": "PASS_INDEPENDENT_EXACT_FLOW_AND_FINITE_MATRIX_CHECKS",
        "lower_bound_all_k": exact_lower_bound(),
        "templates_and_deletions": templates_and_deletions(),
        "full_number_set_construction": all_number_construction(),
        "new_cognitive_axioms": 0,
        "actual_instrument_or_full_quantum_global_theory_certified": False,
        "admission_sha256": hashlib.sha256(admission.read_bytes()).hexdigest(),
        "independent_check_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group()
    flags.add_argument("--write", action="store_true")
    flags.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = payload()
    if args.write:
        with OUT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("CREATED", result["status"])
    else:
        if json.loads(OUT.read_text(encoding="utf-8")) != result:
            raise SystemExit("MISMATCH: independent saved result differs")
        print("READ_ONLY_MATCH", result["status"])


if __name__ == "__main__":
    main()

