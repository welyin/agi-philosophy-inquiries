"""Bound a top-induced scalar quartic using the finite logarithmic running interval."""

import argparse
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
K = 8 * math.pi**2
COEFFICIENTS = (4.1, -19 / 6, -7.0)


def read(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8-sig"))


def envelope(uv_inverse_y2, barrier, gauge, mz, cells=1024):
    length = math.log(1.22e19 / mz)
    width = length / cells
    w_right = uv_inverse_y2
    least_gap = float("inf")
    least_fixed_value = float("inf")
    for index in reversed(range(cells)):
        limits = []
        for beta, (central, error) in zip(COEFFICIENTS, gauge):
            corners = [1 / ((central + sign * 2 * error)**-2 - beta * t / K)
                       for sign in (-1, 1)
                       for t in (index * width, (index + 1) * width)]
            limits.append((min(corners) - 1e-12, max(corners) + 1e-12))
        (g1_lo, g1_hi), (g2_lo, g2_hi), (g3_lo, g3_hi) = limits
        damping = .85 * g1_lo + 2.25 * g2_lo + 8 * g3_lo
        assert uv_inverse_y2 * (.85 * g1_hi + 2.25 * g2_hi + 8 * g3_hi) < 4.5
        exponent = -damping * width / K
        w_left = (math.exp(exponent) * w_right
                  - 4.5 / damping * math.expm1(exponent) + 1e-12)
        w_cell = max(w_left, w_right)
        linear = 16 * g3_hi - 4.5 * g2_lo - .1 * g1_lo
        constant = 3 / 8 * (2 * g2_hi**2 + (g2_hi + .6 * g1_hi)**2)
        assert linear > 0
        coefficient = 3 + linear * w_cell
        gap_at_zero = 6 - constant * w_cell**2
        least_gap = min(least_gap, gap_at_zero - coefficient * barrier - 24 * barrier**2)
        cell_root = 2 * gap_at_zero / (math.sqrt(coefficient**2 + 96 * gap_at_zero) + coefficient)
        least_fixed_value = min(least_fixed_value, cell_root)
        w_right = w_left

    result = dict(inverse_y2_uv_cap=uv_inverse_y2, cells=cells, q_barrier=barrier,
                  y2_mz_lower=1 / w_right, q_beta_normalized_negative_gap=least_gap,
                  smallest_positive_comparison_root=least_fixed_value,
                  quartic_mz_lower=barrier / w_right)
    assert least_gap > 0
    if uv_inverse_y2:
        shortest = math.log(1e3 / mz)
        integral = K / 4.5 * math.log1p(4.5 / K * shortest / uv_inverse_y2)
        decrease = least_gap * integral / (2 * K)
        assert decrease > barrier
        result.update(integral_y2_lower=integral, q_forced_decrease_if_below=decrease)
    return result


def certificate():
    low = read("../1369/inputs.json")
    quartic = read("../1396/inputs.json")
    gauge = [low[name] for name in ("g1_GUT", "g2", "g3")]
    mz = low["mu_ref_GeV"]
    upper = (quartic["lambda_paper"][0] + 2 * quartic["lambda_paper"][1]) / quartic["lambda_paper_factor"]
    pole = envelope(0.0, .2625, gauge, mz)
    finite = envelope(1 / (4 * math.pi), .195, gauge, mz)
    for branch in (pole, finite):
        branch["downward_additive_matching_needed"] = branch["quartic_mz_lower"] - upper
    return dict(positive_quartic_upper_from_mz_data=upper, pole=pole, finite_matching=finite,
                asymptotic_ratio_roots=[(-3 + sign * math.sqrt(585)) / 48 for sign in (-1, 1)])


def compare(actual, saved):
    if isinstance(actual, dict):
        return actual.keys() == saved.keys() and all(compare(v, saved[k]) for k, v in actual.items())
    if isinstance(actual, list):
        return len(actual) == len(saved) and all(compare(x, y) for x, y in zip(actual, saved))
    return math.isclose(actual, saved, rel_tol=2e-10, abs_tol=2e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write", action="store_true")
    args = parser.parse_args()
    answer = certificate()
    filename = HERE / "independent_results.json"
    if args.write:
        with filename.open("x", encoding="utf-8") as handle:
            json.dump(answer, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
    else:
        assert compare(answer, read(filename.name)), "Independent certificate changed"
        print("PASS", answer["pole"]["quartic_mz_lower"], answer["finite_matching"]["quartic_mz_lower"])
