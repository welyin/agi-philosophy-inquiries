"""Fixed fourth moment and the observable single-Higgs portal rate increment."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from recoil import window


HERE = Path(__file__).resolve().parent
DEN = 8 * math.pi**2


def read(filename):
    with (HERE / filename).open(encoding="utf-8-sig") as handle:
        return json.load(handle)


def loop_kernel(ratio, order):
    if not 0 <= ratio < 4:
        raise ValueError("use a closed scalar-pair threshold")
    x, weights = np.polynomial.legendre.leggauss(order)
    weights = weights / 2
    mid = (x + 1) / 2
    u = mid * (1 - mid)
    result = float(weights @ (u / (1 - ratio*u)))
    # Exact beta integrals for the Taylor coefficients, with geometric tail.
    coefficient = 1 / 6
    total = coefficient
    for k in range(1, 50):
        coefficient *= ratio * (k+1)**2 / ((2*k+2)*(2*k+3))
        total += coefficient
    tail = coefficient * (ratio/4) / (1 - ratio/4)
    assert abs(result - total) <= tail + 3e-14
    return result


def response(share, fourth, v, mh, order):
    assert all(s > 0 for s in share) and math.isclose(sum(share), 1)
    masses = [(s*fourth)**0.25 for s in share]
    squared_over_v2 = [(mass/v)**2 for mass in masses]
    terms = [weight*loop_kernel((mh/mass)**2, order)/DEN
             for weight, mass in zip(squared_over_v2, masses)]
    zero = sum(squared_over_v2)/(6*DEN)
    exact = sum(terms)
    return {"masses_GeV": masses, "portal_coefficients": squared_over_v2,
            "suppression_at_mass": exact, "suppression_at_zero": zero}


def bounds(v, mh, deficit, order):
    # deficit is -J0/v^4, a is the pole mass squared in units of v^2.
    a = (mh/v)**2
    y = math.sqrt(deficit + DEN*a)
    assert y > a/2
    return {"finite_mass": y*loop_kernel(a/y, order)/DEN, "zero": y/(6*DEN)}


def calculate():
    new = read("inputs.json")
    old = read(new["weak_and_pole_inputs"])
    trace = read(new["spectrum_inputs"])["SM_trace_over_v4"]
    v, dv = old["matching"]["v_GeV"]
    mh, dm = old["atlas"]["mass_GeV"]
    order = new["quadrature_order"]
    deficit = -trace["center"]
    fourth = deficit*v**4 + DEN*v**2*mh**2
    assert math.isclose(fourth**0.25, 541.1715234452621, rel_tol=1e-10)
    cases = {str(n): response([1/n]*n, fourth, v, mh, order)
             for n in new["equal_species"]}
    cases["unequal_check"] = response([0.1, 0.3, 0.6], fourth, v, mh, order)
    minimum = bounds(v, mh, deficit, order)
    upper_trace = trace["marginal_box"][1]
    lower_trace = trace["marginal_box"][0]
    low = bounds(v+dv, mh-dm, -upper_trace, order)
    high = bounds(v-dv, mh+dm, -lower_trace, order)
    assert math.isclose(minimum["finite_mass"], cases["1"]["suppression_at_mass"])
    assert all(c["suppression_at_mass"] >= minimum["finite_mass"]-1e-13
               for c in cases.values())
    for case in cases.values():
        factor = 1-case["suppression_at_mass"]
        # A channel ratio removes the common factor for arbitrary positive widths.
        assert math.isclose((factor*.02)/(factor*.58), .02/.58)
    electroweak = read("../1369/inputs.json")
    for case in cases.values():
        case["real_pair_recoil"] = window(case["masses_GeV"], v, mh,
                                          electroweak["g1_GUT"][0], electroweak["g2"][0])
    single = cases["1"]["real_pair_recoil"]
    for case in cases.values():
        assert case["real_pair_recoil"]["rate_fb"] >= single["rate_fb"]-1e-18
    gap = cases["4"]["suppression_at_mass"] - minimum["finite_mass"]
    return {"fourth_GeV4": fourth, "fourth_root_GeV": fourth**0.25,
            "minimum": minimum, "box_lower": low, "box_upper": high,
            "fixed_same_moment": cases,
            "equal_model_gap": gap, "half_gap_between_models": gap/2,
            "half_gap_vs_zero": minimum["finite_mass"]/2,
            "future_relative_precision_from_simulation_only":
                new["expected_machine_precision"]["relative"]}


def compare(a, b):
    if isinstance(a, dict):
        assert set(a) == set(b)
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    else:
        assert math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-15)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    opts = parser.add_mutually_exclusive_group()
    opts.add_argument("--write-new", action="store_true")
    opts.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_new:
        with (HERE / "results.json").open("x", encoding="utf-8") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    elif args.check:
        compare(read("results.json"), result)
        print("1374: portal integral, sharp mass-distribution bound and held-out rate verified")
    else:
        print(json.dumps(result, indent=2, allow_nan=False))
