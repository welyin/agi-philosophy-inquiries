"""Non-numbered finite jet/WKB check; no SLE or absolute source certificate.

Default is read-only exact replay. --write creates the result exclusively.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents)
            if (p / "research_cognition_physics").is_dir()
            and (p / "scripts").is_dir())
STAGE = ROOT / "research_cognition_physics" / "archive_1234_"
SCRIPT = Path(__file__).resolve()
MEMO = SCRIPT.with_suffix(".md")
RESULT = SCRIPT.with_suffix(".json")
GATE = STAGE / "_admission/quantum_identity_gate_after1270/original_field_reference_gate_after1270.md"


def ratio(x):
    x = F(x)
    return {"numerator": str(x.numerator), "denominator": str(x.denominator)}


def ceiling(x):
    x = F(x)
    return -((-x.numerator) // x.denominator)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def derivative(poly):
    """Formal d/deta: q'=U1, Uj'=U(j+1); exact Laurent polynomials."""
    out = {}
    for powers, coef in poly.items():
        n = powers[0]
        if n:
            key = list(powers)
            key[0] -= 1
            key[1] += 1
            key = tuple(key)
            out[key] = out.get(key, F(0)) + coef * n
        for j in range(1, 4):
            if powers[j]:
                key = list(powers)
                key[j] -= 1
                key[j + 1] += 1
                key = tuple(key)
                out[key] = out.get(key, F(0)) + coef * powers[j]
        assert powers[4] == 0, "Only two derivatives of r0 are requested."
    return {key: val for key, val in out.items() if val}


def output():
    spec = importlib.util.spec_from_file_location("old1250_wkb_screen", STAGE / "1250/check.py")
    old = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(old)
    saved = json.loads((STAGE / "1250/results.json").read_text(encoding="utf-8"))
    branches = []
    uniform = [0] * 7
    for index, br in enumerate(saved["group_same_history_null"]["branches"]):
        z = [old.I(*v) for v in br["hold"]["all_time_tube"]]
        mc = old.I(*br["frozen_t2"]["M_c"])

        def conformal_rhs(q):
            return [q[6] * v for v in old.bg_rhs(q, mc)]

        cs = old.coefficients(z, conformal_rhs, 6)
        q = [old.Jet(v) for v in cs]
        potential = q[6] * q[6] * (1 - 2 * q[1] * q[1] - old.bg_rhs(q, mc)[1])
        derivative_bounds = []
        for order, coefficient in enumerate(potential.c):
            val = coefficient * math.factorial(order)
            lo, hi = F(val.lo, old.S), F(val.hi, old.S)
            bound = max(abs(lo), abs(hi))
            upper = ceiling(bound)
            derivative_bounds.append({
                "order": order, "lower": ratio(lo), "upper": ratio(hi),
                "absolute_bound": ratio(bound), "integer_upper": upper,
            })
            uniform[order] = max(uniform[order], upper)
        assert F(potential.c[0].lo, old.S) > 6
        assert F(potential.c[0].hi, old.S) < 8
        branches.append({"saved_branch_index": index, "derivatives": derivative_bounds})
    assert uniform == [8, 14, 25, 88, 596, 5373, 55860]

    r0 = {(-2, 2, 0, 0, 0): F(5, 16), (-1, 0, 1, 0, 0): F(-1, 4)}
    dr0 = {
        (-1, 0, 0, 1, 0): F(-1, 4),
        (-2, 1, 1, 0, 0): F(7, 8),
        (-3, 3, 0, 0, 0): F(-5, 8),
    }
    d2r0 = {
        (-1, 0, 0, 0, 1): F(-1, 4),
        (-2, 0, 2, 0, 0): F(7, 8),
        (-2, 1, 0, 1, 0): F(9, 8),
        (-3, 2, 1, 0, 0): F(-29, 8),
        (-4, 4, 0, 0, 0): F(15, 8),
    }
    assert derivative(r0) == dr0
    assert derivative(dr0) == d2r0

    k = 600
    b1, b2, b3, b4 = uniform[1:5]
    c0 = F(b2, 4) + F(5 * b1**2, 16 * k**2)
    c1 = F(b3, 4) + F(7 * b1 * b2, 8 * k**2) + F(5 * b1**3, 8 * k**4)
    c2 = (F(b4, 4) + F(7 * b2**2 + 9 * b1 * b3, 8 * k**2)
          + F(29 * b1**2 * b2, 8 * k**4) + F(15 * b1**4, 8 * k**6))
    assert c0 < 7 and c1 < 23 and c2 < 150
    assert F(k**2) - F(7, k**2) > F(k**2, 2)
    assert F(5, 16) * (8 * 14 * 23) + F(25 * 14, 4) == F(1785, 2)
    assert F(5, 16) * (4 * 23**2 + 42 * 14**2) == F(12935, 4)
    corrected = 75 + F(1785, 2 * k**2) + F(12935, 4 * k**4)
    assert corrected < 76
    bindings = {
        "memo": MEMO, "code": SCRIPT, "old_1250_code": STAGE / "1250/check.py",
        "old_1250_results": STAGE / "1250/results.json",
        "old_1258_proof": STAGE / "1258/proof.md", "published_reference_gate": GATE,
    }
    return {
        "status": "PARTIAL_WKB_BOUND_NO_ABSOLUTE_SOURCE_CERTIFICATE",
        "source_bindings": {key: {"path": str(path.relative_to(ROOT)).replace("\\", "/"),
                                  "sha256": sha(path)} for key, path in bindings.items()},
        "scope": {
            "new_science_groups": 0, "background_ODE_rerun": False,
            "field_mode_integration": False, "absolute_reference_source_certified": False,
            "SLE_tail_certified": False, "torus_subtraction_matching_certified": False,
            "thermal_ratio_certified": False, "old_fixed_scale_pass_claimed": False,
            "zero_mode_changed": False,
        },
        "background_interval": ["2", "9/4"],
        "small_torus_period": "1/100",
        "uniform_U_lower_strict": 6,
        "uniform_U_upper_strict": 8,
        "uniform_derivative_integer_upper": uniform,
        "branches": branches,
        "wkb_checks": {
            "nonzero_k_lower_conservative": k,
            "r0_derivative_identity_exact": True,
            "r0_second_derivative_identity_exact": True,
            "k_squared_r0_upper": ratio(c0),
            "k_squared_r0_prime_upper": ratio(c1),
            "k_squared_r0_second_upper": ratio(c2),
            "r0_integer_constants": [7, 23, 150],
            "Q_positive_uniform": True,
            "k_fourth_r1_upper": ratio(corrected),
            "r1_integer_constant_strict": 76,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = output()
    if args.write:
        with RESULT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
        print("FIRST_WRITE_PARTIAL_WKB_BOUND")
    else:
        saved = json.loads(RESULT.read_text(encoding="utf-8"))
        assert saved == result, "Saved result differs."
        print("READ_ONLY_MATCH_PARTIAL_WKB_BOUND")
    print("absolute_reference_source_certified=False; SLE_tail_certified=False")
    print("U: 6<U<8; derivative bounds:", result["uniform_derivative_integer_upper"][1:])
    print("For all k>=600: |r1|<76/k^4; no field-mode or background integration.")


if __name__ == "__main__":
    main()
