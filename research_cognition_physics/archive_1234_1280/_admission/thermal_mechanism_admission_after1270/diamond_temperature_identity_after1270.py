"""Exact, bounded admission check; no field/instrument simulation.

Default execution compares the recomputed result with the saved adjacent JSON.
--write is only for first author generation and refuses to overwrite.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sys


def derive():
    # Dimensionless coordinates H*eta, H*r and H*ell; O uses strict bounds.
    eta_low, eta_high, r_high = -F(9, 8), -F(7, 8), F(1, 8)
    u_low, v_high = eta_low-r_high, eta_high+r_high
    ell = (F(1), F(2))
    beta_E = tuple(1-1/(2*x) for x in ell)
    gap = abs(beta_E[1]-beta_E[0])
    minimax = gap/2
    midpoint = sum(beta_E)/2
    assert all(u_low > -2*x and v_high < 0 for x in ell)
    assert minimax == F(1, 8)
    assert max(abs(midpoint-x) for x in beta_E) == minimax
    # The continuum bound is proved by the triangle inequality in the note.
    return {
        "status": "PASS_EXACT_REGION_IDENTITY_DIAGNOSTIC",
        "H_ell": list(map(str, ell)),
        "strict_lower_H_u": str(u_low),
        "strict_upper_H_v": str(v_high),
        "beta_E_at_T_zero": list(map(str, beta_E)),
        "sum_error_lower_bound": str(gap),
        "maximum_error_lower_bound": str(minimax),
        "minimax_logratio": str(midpoint),
        "scope": "Same local state and complete instrument, only region label changes",
        "physical_probabilities_computed": False,
        "actual_instrument_realized": False,
        "new_science_rounds": 0,
        "new_axioms": 0,
    }


if __name__ == "__main__":
    result = derive()
    dst = Path(__file__).with_suffix(".json")
    if sys.argv[1:] == ["--write"]:
        with dst.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print("WROTE_EXACT_REGION_IDENTITY_DIAGNOSTIC")
    elif not sys.argv[1:]:
        assert json.loads(dst.read_text(encoding="utf-8-sig")) == result
        print("READ_ONLY_MATCH_EXACT_REGION_IDENTITY_DIAGNOSTIC")
    else:
        raise SystemExit("Expected no arguments or --write")
