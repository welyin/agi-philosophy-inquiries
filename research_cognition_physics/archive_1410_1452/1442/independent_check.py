"""Rational mass-interval bounds from all three complementary minors."""
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent


def bounds(values):
    return values[0] - 2 * values[1], values[0] + 2 * values[-1]


def sqrt_bounds(x):
    grid = 10**14
    k = math.isqrt(x.numerator * grid**2 // x.denominator)
    return F(k, grid), F(k + 1, grid)


def verify():
    data = json.loads((HERE.parent / "1386/inputs.json").read_text(encoding="utf8"), parse_float=F)
    e, mu, tau = map(bounds, data["yukawa"]["charged_leptons"])
    f = mu[0]**2 / (e[1]*tau[1])
    a_lo, a_hi = bounds(data["DayaBay_amplitude"])
    w_lo = (1 - sqrt_bounds(1-a_lo)[1]) / 2
    w_hi = (1 - sqrt_bounds(1-a_hi)[0]) / 2
    light_norm = (e[1]**2 + mu[1]**2) / tau[0]**2
    n = 1024
    results = {}
    for order in ("NO", "IO"):
        osc = data[order]
        solar = bounds(osc["solar_gap"])
        atm = bounds(osc["atmospheric_gap"])
        s12 = bounds(osc["solar_sine_squared"])
        s23 = bounds(osc["atmospheric_sine_squared"])
        e_w = [(1-s12[1])*(1-w_hi), s12[0]*(1-w_hi), w_lo]
        tau3_lo, tau3_hi = (1-s23[1])*(1-w_hi), (1-s23[0])*(1-w_lo)
        a = solar[0] / (solar[0] + atm[1]), solar[1] / (solar[1] + atm[0])
        min_f2, min_r2 = None, None
        for j in range(n):
            lo, hi = F(j, n), F(j+1, n)
            if order == "NO":
                middle_lo, middle_hi = (lo+(1-lo)*a[0], hi+(1-hi)*a[1])
                lower, upper = (lo, middle_lo, 1), (hi, middle_hi, 1)
                tau_upper = tau3_hi + (1-tau3_hi)*middle_hi
            else:
                middle_lo, middle_hi = (1-(1-lo)*a[1], 1-(1-hi)*a[0])
                lower, upper = (middle_lo, 1, lo), (middle_hi, 1, hi)
                tau_upper = 1 - tau3_lo*(1-hi)
            d = min(1, tau_upper + light_norm)
            minor_sum = sum(e_w[i]*lower[(i+1)%3]*lower[(i+2)%3] for i in range(3))
            f2 = f**2 * minor_sum**3 / (d**3 * math.prod(upper))
            r2 = (mu[0]/tau[1])**2 * minor_sum / d**2
            min_f2 = f2 if min_f2 is None else min(f2, min_f2)
            min_r2 = r2 if min_r2 is None else min(r2, min_r2)
        results[order] = {"F_lower_rounded_down": float(sqrt_bounds(min_f2)[0]),
                          "middle_over_heavy_lower_rounded_down": float(sqrt_bounds(min_r2)[0])}
        assert min_f2 > 1 and min_r2 > F(75, 1000000)**2
    return {"cells_per_order": n, "squared_mass_ratio_domain": [0,1], "bounds": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write-new", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    res = verify()
    file = HERE / "independent_results.json"
    if args.write_new:
        with file.open("x", encoding="utf8") as stream:
            json.dump(res, stream, indent=2)
            stream.write("\n")
    else:
        assert res == json.loads(file.read_text(encoding="utf8"))
    print(res)
