"""Read-only reproduction of published summaries, not raw-data reanalysis."""
from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


def calculate() -> dict:
    # Page--Geilker PRL 47, p.981. AB -> -1, BA -> +1.
    labels = [-1, -1, -1, 1, 1, -1, 1, -1, 1, 1]
    values = [
        "-61.3", "-63.9", "-36.0", "69.2", "36.1",
        "-48.8", "46.4", "-45.2", "51.3", "59.6",
    ]
    y = list(map(F, values))
    n = len(y)
    xbar = F(sum(labels), n)
    ybar = sum(y) / n
    xx = sum((F(x) - xbar) ** 2 for x in labels)
    yy = sum((v - ybar) ** 2 for v in y)
    xy = sum((F(x) - xbar) * (v - ybar) for x, v in zip(labels, y))
    r2 = xy * xy / (xx * yy)
    assert n == 10 and sum(labels) == 0
    assert r2 == F(744769, 777349)

    # Exact normalization of 35/32 * (1-r^2)^3 on [-1,1].
    normalization = F(35, 32) * 2 * (1 - F(3, 3) + F(3, 5) - F(1, 7))
    assert normalization == 1
    with localcontext() as ctx:
        ctx.prec = 65
        d = lambda q: Decimal(q.numerator) / Decimal(q.denominator)
        r = d(r2).sqrt()
        antiderivative = Decimal(35) / 32 * (
            r - r**3 + Decimal(3) * r**5 / 5 - r**7 / 7
        )
        one_tail = Decimal("0.5") - antiderivative
        assert Decimal("4.29e-7") < one_tail < Decimal("4.30e-7")
        numeric = {
            "pearson_r": format(r, ".35g"),
            "gaussian_null_one_tail": format(one_tail, ".35g"),
            "gaussian_null_two_tail": format(2 * one_tail, ".35g"),
        }

    # Conditional exchangeability diagnostic, not the paper's test.
    scores = [
        sum((1 if j in selected else -1) * v for j, v in enumerate(y))
        for selected in combinations(range(n), 5)
    ]
    observed_score = sum(x * v for x, v in zip(labels, y))
    one_count = sum(score >= observed_score for score in scores)
    two_count = sum(abs(score) >= abs(observed_score) for score in scores)
    assert len(scores) == 252 and one_count == 1 and two_count == 2
    return {
        "status": "published_summary_reproduction_only",
        "input_source": "Page--Geilker, Phys. Rev. Lett. 47 (1981), p.981",
        "labels_AB_minus_BA_plus": labels,
        "equilibrium_change_cm_strings": values,
        "exact": {
            "pearson_r_squared": str(r2),
            "sample_mean_cm": str(ybar),
            "signed_sample_mean_cm": str(observed_score / n),
            "null_density_integral": str(normalization),
        },
        "decimal_65_digit_calculation": numeric,
        "conditional_permutation_diagnostic": {
            "assumption": "fixed five/five labels exchangeable with the ten readings",
            "assignments": len(scores),
            "one_tail_count": one_count,
            "two_tail_count": two_count,
            "one_tail_probability": str(F(one_count, len(scores))),
            "two_tail_probability": str(F(two_count, len(scores))),
            "not_original_paper_test": True,
        },
        "limits": [
            "No refit of raw torsion time series or independent experimental calibration.",
            "Continuous-null p values require an additional spherical/Gaussian noise contract.",
            "Exchangeability was not experimentally certified here.",
            "Neither p value is a probability that all classical gravity is false.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--save", action="store_true", help="exclusive first save only")
    args = parser.parse_args()
    target = Path(__file__).with_name("results.json")
    result = calculate()
    if args.save:
        with target.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("Saved new results.json; no prior result overwritten.")
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert result == saved, "saved result differs from default calculation"
        print("PASS: published ten summaries, exact r^2, null integral and 252 assignments.")
        print("Read-only; no new experiment, raw-data fit or gravity-wide exclusion.")


if __name__ == "__main__":
    main()
