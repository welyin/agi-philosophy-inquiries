"""1264: one exact-arithmetic author group for the original record/source budget.

No full-field truncation and no second reconstruction of the frozen 1253 material.
Default is read-only. --write-results creates results.json exclusively.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / "results.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ceil_root(x: F, degree: int) -> int:
    """Exact smallest integer n with n**degree >= x; no float root."""
    lo, hi = 0, 1
    while hi**degree * x.denominator < x.numerator:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**degree * x.denominator >= x.numerator:
            hi = mid
        else:
            lo = mid
    assert hi**degree * x.denominator >= x.numerator
    assert hi == 0 or (hi - 1)**degree * x.denominator < x.numerator
    return hi


def calculate() -> dict:
    assert ROOT.name == "AGI的哲学思考"
    assert (ROOT / "research_cognition_physics" / "archive_1234_" / "1253").is_dir()

    # Exact trigonometric orthogonality on [-pi/2, pi/2].
    # cos^4(theta) = 3/8 + (1/2) cos(2 theta) + (1/8) cos(4 theta).
    coeff = {0: F(3, 8), 2: F(1, 2), 4: F(1, 8)}
    norm_over_pi = coeff[0]**2 + sum(c*c / 2 for k, c in coeff.items() if k)
    first_over_pi = sum(F(k*k, 2)*c*c for k, c in coeff.items() if k)
    second_over_pi = sum(F(k**4, 2)*c*c for k, c in coeff.items() if k)
    d1 = first_over_pi / norm_over_pi
    d2 = second_over_pi / norm_over_pi
    lap = 3*d2 + 6*d1*d1
    normalization = F(128, 35)*norm_over_pi
    assert norm_over_pi == F(35, 128)
    assert normalization == 1
    assert d1 == F(16, 7) and d2 == F(512, 35)
    assert lap == F(18432, 245) < 81
    # pi < 4 then ||Delta psi_L|| < 36/L^2.
    assert 9*4 == 36

    r = F(1, 10000)
    assert F(3, 4) + r/2 < F(5, 6)  # all total ages in the 1253 prior window
    assert F(16) < 19 < F(81, 4)    # original 26 < m < 32, width exactly 5
    C2 = F(91062363423308015531, 3)
    C4 = F(610242449441917720898626279056227, 3)
    L = max(10000, ceil_root(F(625, 52)*C2, 2),
            ceil_root(F(125, 281216)*C4, 4))
    R = 2*L
    assert L == 19100606975 and R == 38201213950
    source2 = F(5, 208)*C2/R**2
    source4 = C4/F(281216*R**4)
    source_D = source2 + source4
    assert source2 <= F(1, 2000) and source4 <= F(1, 2000)
    assert source_D <= F(1, 1000)

    ideal_contrast = F(9, 16)*r
    packet_each = F(18, 13*L**2)
    packet_both = 2*packet_each
    source_both = 4*(r*source_D + 33*r*r)
    contrast_H = ideal_contrast - packet_both
    contrast_Q = contrast_H - source_both
    universal_H = ideal_contrast - F(36, 13*10000**2)
    universal_Q = universal_H - 4*(r/F(1000) + 33*r*r)
    assert universal_H == F(73089, 1300000000)
    assert universal_Q == F(70853, 1300000000)
    assert contrast_Q >= universal_Q > F(1, 20000)
    assert contrast_H >= universal_H > F(1, 20000)

    kinetic2_upper = F(192, 7*L**2)
    initial_H2_upper = (32*(1+r))**2 + kinetic2_upper
    assert initial_H2_upper < 33**2
    luders_work = F(5, 2)*(1+r)**2
    assert luders_work == F(100020001, 40000000)

    # These files are inherited inputs, not a repeated material experiment.
    old = HERE.parent / "1253"
    old_result = json.loads((old / "results.json").read_text(encoding="utf-8-sig"))
    assert old_result["status"] == "PASS"
    sources = json.loads((HERE / "source_results.json").read_text(encoding="utf-8-sig"))
    assert sources["status"] == "PASS_AUTHOR_WITH_SCOPE"
    assert F(sources["C2_upper"]) == C2 and F(sources["C4_upper"]) == C4
    assert int(sources["finite_L"]) == L and int(sources["finite_R"]) == R
    assert F(sources["source_difference_upper"]) == source_D
    assert F(sources["contrast_lower"]) == contrast_Q
    assert sources["source_bindings"]["source_proof"] == sha(HERE / "source_proof.md")
    assert sources["source_bindings"]["source_code"] == sha(HERE / "source_check.py")

    paths = [
        HERE / "check.py", HERE / "record_proof.md",
        HERE / "source_proof.md", HERE / "source_check.py", HERE / "source_results.json",
        old / "proof.md", old / "check.py", old / "results.json",
    ]
    bindings = {p.relative_to(ROOT).as_posix(): sha(p) for p in paths}
    exact = {
        "r": r, "one_dimensional_norm_over_pi": norm_over_pi,
        "packet_normalization": normalization,
        "derivative_first_ratio": d1, "derivative_second_ratio": d2,
        "three_dimensional_laplacian_coefficient": lap,
        "C2_profile_upper": C2, "C4_profile_upper": C4,
        "source_D_second_order_upper": source2,
        "source_D_fourth_order_upper": source4,
        "source_D_upper": source_D,
        "ideal_record_contrast_lower": ideal_contrast,
        "each_packet_probability_error_upper": packet_each,
        "two_packet_probability_error_upper": packet_both,
        "two_source_probability_error_upper": source_both,
        "actual_H_record_contrast_lower": contrast_H,
        "linear_Q_record_contrast_lower": contrast_Q,
        "universal_H_contrast_lower": universal_H,
        "universal_Q_contrast_lower": universal_Q,
        "declared_strict_contrast_lower": F(1, 20000),
        "common_free_retiming_one_task_error_lower": contrast_H/2,
        "initial_material_H_squared_upper": initial_H2_upper,
        "nonselective_native_luders_material_energy_change_upper": luders_work,
    }
    return {
        "status": "PASS_AUTHOR_WITH_SCOPE", "round": 1264, "author_groups": 1,
        "group": "finite_original_record_same_source_joint_budget",
        "finite_L": str(L), "finite_R": str(R),
        "rational_authoritative": {k: str(v) for k, v in exact.items()},
        "display_only": {
            "source_D": float(source_D),
            "actual_H_contrast_lower": float(contrast_H),
            "linear_Q_contrast_lower": float(contrast_Q),
            "universal_joint_contrast_lower": float(universal_Q),
            "luders_material_energy_bound": float(luders_work),
        },
        "scope": [
            "two known normal preparations, no free path-label system",
            "original 1253 M/J/Mr and arbitrary reference retained",
            "same static KG parent and source comparison, not generated geometry",
            "native compact support is not local KG Cauchy support",
            "common coordinate endpoint and native effect remain declared contracts",
            "Luders bound concerns total material energy, not full instrument sources",
        ],
        "bindings_sha256": bindings,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_results:
        with RESULT.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        mode = "EXCLUSIVE_FIRST_WRITE"
    else:
        saved = json.loads(RESULT.read_text(encoding="utf-8-sig"))
        if saved != result:
            raise AssertionError("saved 1264 record results or bound inputs do not match")
        mode = "READ_ONLY_MATCH"
    print(json.dumps({
        "status": result["status"], "mode": mode, "author_groups": 1,
        "finite_L": result["finite_L"],
        "joint_contrast_lower": result["rational_authoritative"]["universal_Q_contrast_lower"],
        "declared_contrast": "1/20000",
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()

