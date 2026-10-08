"""Author-local exact rederivation and narrowly scoped receipt for round 1040.

Default mode compares without writing. --write exclusively creates the receipt.
The exact checks below do not call the scientific constitutive implementation.
They are author checks, not an external independent-agent review.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
PHYSICS = HERE.parents[1]
RECEIPT = HERE / "research_round_1040_checks.json"
OWN = (
    "../research_note_1040.md",
    "nonlinear_common_cone.py",
    "nonlinear_common_cone_results.json",
    "selection_audit.md",
    "input_dependency_update.md",
    "NEXT.md",
    "review.md",
    "verify_round1040.py",
)


def matvec(a, x):
    return [sum((v*w for v, w in zip(row, x)), F(0)) for row in a]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def exact_null_screen():
    """Build the Ampere rows directly, separately from the science matrix code."""
    backgrounds = [(F(1), F(0)), (F(0), F(2)),
                   (F(3, 5), F(4, 5)), (F(-2, 3), F(5, 7))]
    hessians = [(F(0), F(0), F(0)), (F(1), F(0), F(0)),
                (F(0), F(1), F(0)), (F(0), F(0), F(1)),
                (F(2, 7), F(-3, 8), F(4, 9))]
    cases = 0
    for e, b in backgrounds:
        J = [[-b, e], [e, b]]
        r2 = e*e+b*b
        assert r2 > 0
        for z in ([F(1), F(0)], [F(0), F(1)]):
            assert matvec(J, matvec(J, z)) == [r2*x for x in z]
        for ss, sp, pp in hessians:
            for u, w in [(F(1), F(0)), (F(0), F(1)), (F(2, 5), F(-3, 7))]:
                LS, LP = F(7, 3), F(-2, 9)
                dS, dP = e*w-b*u, b*w+e*u
                dLS, dLP = ss*dS+sp*dP, sp*dS+pp*dP
                Dy = LS*u-LP*w
                Dz = LS*w+LP*u+e*dLS+b*dLP
                Hy = -LS*w-LP*u
                Hz = LS*u-LP*w+b*dLS-e*dLP
                residual = [Hz-Dy, Hy+Dz]
                expected = matvec(J, matvec([[ss, sp], [sp, pp]], matvec(J, [u, w])))
                expected[0] *= -1
                assert residual == expected
                if ss == sp == pp == 0:
                    assert residual == [0, 0]
                cases += 1
    return cases


def exact_pure_bi():
    """Exact two transverse solutions of k x mu^-1(k x e)+omega² eps e=0."""
    count = 0
    for R in (F(5, 4), F(3, 2), F(2)):
        eps, mu = [1/R, 1/R, R], [1/R, 1/R, 1/R**3]
        assert all(x > 0 for x in eps+mu)
        for k in ([F(1), F(0), F(0)], [F(0), F(0), F(1)],
                  [F(3, 5), F(0), F(4, 5)]):
            w2 = k[2]**2+k[0]**2/R**2
            modes = [[F(0), F(1), F(0)], [k[2], F(0), -k[0]/R**2]]
            assert dot(modes[0], modes[1]) == 0 and dot(modes[1], modes[1]) > 0
            for e in modes:
                D = [x*y for x, y in zip(eps, e)]
                assert dot(k, D) == 0
                ke = cross(k, e)
                magnetic = cross(k, [x*y for x, y in zip(mu, ke)])
                assert [x+w2*y for x, y in zip(magnetic, D)] == [0, 0, 0]
                count += 1
        T, B2 = F(1), R*R-1
        rho, pressure, speed = T*(R-1), T*(1-1/R), 1/R
        assert speed == 1/(1+rho/T)
        assert rho >= pressure > 0
        assert 1-speed**2 == B2/(T+B2)
    gap, tolerance = F(4, 5)-F(1, 2), F(1, 100)
    assert gap == F(3, 10) and gap-2*tolerance == F(7, 25)
    lower_T = (1-tolerance)/tolerance
    assert lower_T == 99 and 1/(lower_T+1) == tolerance
    return count


def local_links():
    count = 0
    for name in OWN:
        path = (HERE/name).resolve()
        if path.suffix != ".md":
            continue
        content = path.read_text(encoding="utf8")
        for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", content):
            target = target.strip().strip("<>")
            if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*://", target) or target.startswith("#"):
                continue
            target = unquote(target.split("#", 1)[0])
            assert (path.parent/target).resolve().exists(), (name, target)
            count += 1
    return count


def build_receipt():
    run = subprocess.run([sys.executable, "-B", "-X", "utf8",
                          str(HERE/"nonlinear_common_cone.py")],
                         capture_output=True, text=True, encoding="utf8", check=True)
    science_run = json.loads(run.stdout.strip())
    result = json.loads((HERE/"nonlinear_common_cone_results.json").read_text(encoding="utf8"))
    assert result["round"] == 1040 and result["scientific_baseline"] == 1038
    assert result["new_calibration_groups"] == 1 and not result["goal_completed"]
    assert result["new_adopted_cognitive_axioms"] == 0
    history = result["historical_sha256"]
    assert len(history) == 18
    for path, expected in history.items():
        assert hashlib.sha256((PHYSICS/path).read_bytes()).hexdigest() == expected, path
    # The receipt link is deliberately not required before its exclusive first write.
    links = local_links() if RECEIPT.exists() else None
    return {
        "round": 1040, "date": "2026-10-08", "all_checks_passed": True,
        "review_kind": "author_local_rederivation_plus_root_admission_pre_review",
        "external_independent_final_review_certified_here": False,
        "new_calibration_groups": 0,
        "same_round_science_calibration_groups": 1,
        "exact_direct_Maxwell_null_screen_checks": exact_null_screen(),
        "exact_BI_physical_polarization_checks": exact_pure_bi(),
        "exact_difference": "3/10", "conditional_error_margin": "7/25",
        "science_read_only_compare_passed": science_run["all_checks_passed"],
        "max_full_symbol_residual": science_run["max_full_symbol_residual"],
        "max_same_action_stress_residual": science_run["max_stress_residual"],
        "historical_sha256": history,
        "owned_sha256": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in OWN},
        "freeze_scope": "explicit_author_assets_only_no_navigation_no_parallel_no_root_integration",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = build_receipt()
    if args.write:
        with RECEIPT.open("x", encoding="utf8") as f:
            json.dump(receipt, f, ensure_ascii=False, indent=2)
            f.write("\n")
    else:
        assert receipt == json.loads(RECEIPT.read_text(encoding="utf8"))
    links = local_links()
    print(json.dumps({"round": 1040, "all_checks_passed": True,
                      "mode": "exclusive_write" if args.write else "read_only_compare",
                      "owned_files": len(OWN), "historical_files": len(receipt["historical_sha256"]),
                      "local_links": links,
                      "exact_null_screen_checks": receipt["exact_direct_Maxwell_null_screen_checks"],
                      "exact_BI_polarization_checks": receipt["exact_BI_physical_polarization_checks"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
