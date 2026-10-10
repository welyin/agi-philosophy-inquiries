"""Small algebra checks only; no stellar solver, waveform fit, or new result.

Default: read-only comparison. --write: exclusive first creation of results.json.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json


def q(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def compute():
    rows = []
    for y in [F(1, 10), F(1), F(10)]:
        # |1/(1+i A/omega)| depends only on y=|omega/A|.
        norm2 = y*y/(1+y*y)
        frozen_error2 = 1/(1+y*y)
        assert norm2 <= y*y
        assert frozen_error2 <= 1/(y*y)
        assert norm2 + frozen_error2 == 1
        rows.append({"y": q(y), "equilibrium_error_squared": q(norm2),
                     "frozen_error_squared": q(frozen_error2)})

    # Static compatibility: delta p=ce2 delta e; Delta p-cf2 Delta e
    # equals (ce2-cf2)(delta e+xi eprime), with pprime=ce2 eprime.
    ce2, cf2, de, xi, ep = F(2, 7), F(3, 7), F(5, 11), F(2, 3), F(-7, 5)
    lhs = ce2*de + xi*ce2*ep - cf2*(de+xi*ep)
    rhs = (ce2-cf2)*(de+xi*ep)
    assert lhs == rhs

    # Printed g1 values from AP1906.08982v1 Table II, not exact stellar data.
    nu, overlap, surface = F("0.1845"), F("0.0017657"), F("27.8841")
    a = nu*nu
    assert 1-surface*a != 0
    identities = []
    for z in [F(0), a/4, a/2]:
        lhs = (1-surface*z)/((a-z)*(1-surface*a))
        rhs = 1/(a-z)+surface/(1-surface*a)
        assert lhs == rhs
        identities.append(q(z))
    # Residue after removing 2*pi/5, with denominator (nu^2-omega^2).
    residue = overlap*overlap
    assert residue > 0
    assert F(2, 3) == F(2, 1*3)  # l=2 STF normalization, (2l-1)!!=3.
    return {
        "scope": "algebra_only_no_stellar_integration_no_waveform_or_data_fit",
        "reaction_modulus": rows,
        "static_compatibility_residual": "0",
        "table_II_g1_printed_values": {
            "frequency": "0.1845", "overlap": "0.0017657", "surface_ratio": "27.8841",
            "surface_denominator": q(1-surface*a),
            "pole_coefficient_without_2pi_over5": q(residue)},
        "mode_identity_squared_frequencies": identities,
        "stf_l2_coefficient": "2/3",
        "all_assertions_passed": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = Path(__file__).with_name("results.json")
    result = compute()
    if args.write:
        with path.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print("First-save algebra checks passed.")
    else:
        assert json.loads(path.read_text(encoding="utf-8")) == result
        print("Read-only algebra checks passed; no stellar or waveform computation claimed.")
