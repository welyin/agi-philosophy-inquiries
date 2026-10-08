"""1041 author-local independent finite checks; explicitly bounded freeze scope.

Default read-only comparison. --write exclusively creates first receipt.
Does not count a new science group or certify external independent review.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
PHYSICS = HERE.parents[1]
RECEIPT = HERE/"research_round_1041_checks.json"
OWN = ("../research_note_1041.md", "common_relay_time.py",
       "common_relay_time_results.json", "selection_audit.md",
       "input_dependency_update.md", "NEXT.md", "review.md",
       "verify_round1041.py")


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def det(matrix):
    a = [r[:] for r in matrix]; value = F(1); n = len(a)
    for i in range(n):
        j = next((j for j in range(i, n) if a[j][i]), None)
        if j is None: return F(0)
        if j != i:
            a[i], a[j] = a[j], a[i]; value *= -1
        p = a[i][i]; value *= p
        for j in range(i+1, n):
            ratio = a[j][i]/p
            a[j] = [x-ratio*y for x, y in zip(a[j], a[i])]
    return value


def direct_ring(eps, flips):
    # Direct assignment enumeration; not author science's matrix multiplication.
    weight = F(0)
    for values in itertools.product((0, 1), repeat=3):
        term = F(1)
        for j in range(3):
            wanted = values[j] ^ flips[j]
            term *= 1-eps[j] if values[(j+1) % 3] == wanted else eps[j]
        weight += term
    return weight


def exact_checks(science):
    axes = [[F(1), F(0), F(0), F(0)],
            [F(-3, 5), F(4, 5), F(0), F(0)],
            [F(-3, 5), F(-4, 5), F(0), F(0)]]
    metrics = [[[F(i == j)-F(17, 16)*u[i]*u[j] for j in range(4)]
               for i in range(4)] for u in axes]
    assert all(det(g) == F(-1, 16) for g in metrics)
    pair_checks = 0
    for i, j in itertools.combinations(range(3), 2):
        n = [a+b for a, b in zip(axes[i], axes[j])]
        for u in (axes[i], axes[j]):
            a = dot(n, u); b2 = dot(n, n)-a*a
            b = F(4, 5) if i == 0 else F(24, 25)
            assert b*b == b2
            v = [4*x-(y-a*x)/b for x, y in zip(u, n)]
            assert dot(u, v) == 4 and dot(v, v) == 17
            projection = dot(n, v)
            assert projection == 4*a-b
            assert projection**2/(17*dot(n, n)) > F(1, 25)
            pair_checks += 1
    legs = [[F(6, 5), F(0), F(0), F(0)], axes[1], axes[2]]
    assert all(sum(v[k] for v in legs) == 0 for k in range(4))
    for u, v, g in zip(axes, legs, metrics):
        assert dot(u, v) > 0 and sum(v[i]*g[i][j]*v[j] for i in range(4) for j in range(4)) < 0
    cert = science["finite_certificate_reduction"]
    points = [[F(x) for x in p] for p in cert["points"]]
    weights = [F(x) for x in cert["weights"]]
    assert 1 < len(points) <= 5 and sum(weights) == 1 and all(w > 0 for w in weights)
    assert all(sum(w*p[k] for w, p in zip(weights, points)) == 0 for k in range(4))
    for p in points:
        assert any(dot(u, p) > 0 and dot(u, p)**2 > 16*(dot(p, p)-dot(u, p)**2) for u in axes)
    cases = 0
    for eps in ((F(0),)*3, (F(1, 100), F(1, 20), F(1, 10)), (F(1, 10),)*3,
                (F(1, 2), F(1, 10), F(1, 20))):
        prod = (1-2*eps[0])*(1-2*eps[1])*(1-2*eps[2])
        for flips in itertools.product((0, 1), repeat=3):
            assert direct_ring(eps, flips) == 1+(-1)**sum(flips)*prod
            cases += 1
    functions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    counts = {0: 0, 1: 0, 2: 0}
    for local in itertools.product(functions, repeat=3):
        z = sum(all(values[(j+1) % 3] == local[j][values[j]] for j in range(3))
                for values in itertools.product((0, 1), repeat=3))
        counts[z] += 1
    assert counts == {0: 4, 1: 56, 2: 4}
    assert [F(r["max_space_per_clock"]) for r in science["distinct_compatible_cones"]] == [F(1), F(1, 2), F(1, 4)]
    return {"independent_metric_determinants": 3, "exact_full_cone_extremizers": pair_checks,
            "reduced_zero_certificate_verified_segments": len(points),
            "direct_assignment_heterogeneous_noise_cases": cases,
            "all_deterministic_triples": 64,
            "deterministic_normalization_histogram": {str(k): v for k, v in counts.items()}}


def links(allow_missing_receipt=False):
    count = 0
    for name in OWN:
        p = (HERE/name).resolve()
        if p.suffix != ".md": continue
        for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", p.read_text(encoding="utf8")):
            target = target.strip().strip("<>")
            if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*://", target) or target.startswith("#"): continue
            q = (p.parent/unquote(target.split("#", 1)[0])).resolve()
            if allow_missing_receipt and q == RECEIPT: continue
            assert q.exists(), (name, target)
            count += 1
    return count


def build_receipt():
    result = json.loads((HERE/"common_relay_time_results.json").read_text(encoding="utf8"))
    run = subprocess.run([sys.executable, "-B", "-X", "utf8", str(HERE/"common_relay_time.py")],
                         check=True, capture_output=True, text=True, encoding="utf8")
    output = json.loads(run.stdout.strip())
    assert output["all_checks_passed"] and result["round"] == 1041
    assert result["new_calibration_groups"] == 1 and result["new_adopted_cognitive_axioms"] == 0
    assert not result["goal_completed"]
    history = result["historical_sha256"]
    assert len(history) == 17
    for path, expected in history.items():
        assert hashlib.sha256((PHYSICS/path).read_bytes()).hexdigest() == expected, path
    links(allow_missing_receipt=True)
    return {"round": 1041, "date": "2026-10-08", "all_checks_passed": True,
            "kind": "author_second_algorithm_and_root_admission_pre_review",
            "external_independent_review_certified_here": False,
            "additional_science_groups": 0, "same_round_science_groups": 1,
            "science_read_only_compare_passed": True, "independent_checks": exact_checks(result),
            "historical_sha256": history,
            "owned_sha256": {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in OWN},
            "freeze_scope": "explicit_author_assets_only_excludes_root_and_independent_reviewer_files"}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--write", action="store_true")
    args = parser.parse_args(); result = build_receipt()
    if args.write:
        with RECEIPT.open("x", encoding="utf8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2); f.write("\n")
    else:
        assert result == json.loads(RECEIPT.read_text(encoding="utf8"))
    print(json.dumps({"round": 1041, "all_checks_passed": True,
                      "mode": "exclusive_write" if args.write else "read_only_compare",
                      "owned_files": len(OWN), "historical_files": len(result["historical_sha256"]),
                      "local_links": links(), "checks": result["independent_checks"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
