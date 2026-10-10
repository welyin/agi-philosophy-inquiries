"""Count-0 conventions only: no solar solution, likelihood or new empirical fit."""
from fractions import Fraction as F
from decimal import Decimal
import json
from pathlib import Path
import sys


TABLE = [
    ("pp", "5.960e10", "13.099", "0.91864"),
    ("Be7", "4.854e9", "12.552", "0.071693"),
    ("pep", "1.425e8", "11.920", "0.0019987"),
    ("N13", "2.795e8", "12.658", "0.0041630"),
    ("O15", "2.067e8", "12.368", "0.0030082"),
    ("F17", "5.350e6", "12.365", "0.000077841"),
    ("B8", "5.025e6", "6.6305", "0.000039205"),
    ("hep", "7.950e3", "3.7355", "0.000000034944"),
]


def rounded_interval(text):
    value = Decimal(text)
    half_unit = F(10) ** value.as_tuple().exponent / 2
    return F(value)-half_unit, F(value)+half_unit


def calculate():
    photon = F("8.4984e11")
    plo, phi = rounded_interval("8.4984e11")
    rows = []
    for name, flux, alpha, beta in TABLE:
        flo, fhi = rounded_interval(flux)
        alo, ahi = rounded_interval(alpha)
        blo, bhi = rounded_interval(beta)
        lower, upper = alo*flo/phi, ahi*fhi/plo
        overlap = max(lower, blo) <= min(upper, bhi)
        assert overlap, name
        derived = F(flux)*F(alpha)/photon
        rows.append({"reaction": name, "printed_beta": beta,
                     "beta_from_printed_inputs": float(derived),
                     "rounding_intervals_overlap": overlap})
    # Basis coefficients of arbitrary (r_He4,r_He3,r_N14,Sdot,Lother).
    # Lgamma = Lnuc - Lnu - Sdot - Lother; Lnu cancels identically.
    Q4, Q3, Q14 = F("26.73097"), F("6.936"), F("11.710")
    c = (F(2), F(1), F(1), F(0), F(0))
    Q = (Q4, Q3, Q14, F(0), F(0))
    lhs = tuple(Q4*x/2 for x in c)
    gamma = (Q4, Q3, Q14, F(-1), F(-1))
    corrections = (F(0), Q4/2-Q3, Q4/2-Q14, F(1), F(1))
    rhs = tuple(x+y for x, y in zip(gamma, corrections))
    assert lhs == rhs
    assert corrections[1] > 0 and corrections[2] > 0
    lower, upper = F("1.038")-F("0.060"), F("1.038")+F("0.069")
    assert lower < 1 < upper
    return {
        "scope": "count0: rounded table/linear energy ledger/published interval only",
        "source_version": "2311.16226v2 Table 1 and Eq 20; no-LC CNO-Rfixed",
        "table_checks": rows,
        "sum_printed_beta_at_reference_f_equal_one": float(sum(F(r[3]) for r in TABLE)),
        "reference_sum_is_not_imposed_as_exact_unity": True,
        "energy_ledger_basis_order": ["r_He4", "r_He3", "r_N14", "Sdot", "Lother"],
        "linear_identity_coefficients": [str(x) for x in lhs],
        "correction_coefficients": [str(x) for x in corrections],
        "published_1sigma_interval_exact": [str(lower), str(upper)],
        "published_interval_contains_unity": True,
        "no_new_pvalue_or_stellar_prediction": True,
        "all_assertions_passed": True,
    }


if __name__ == "__main__":
    result = calculate()
    target = Path(__file__).with_name("results.json")
    if sys.argv[1:] == ["--save-exclusive"]:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    elif sys.argv[1:]:
        raise SystemExit("Use no arguments for read-only checks, or --save-exclusive once.")
    else:
        assert result == json.loads(target.read_text(encoding="utf-8"))
    print("PASS: solar energy ledger conventions and published no-LC interval; count0.")
