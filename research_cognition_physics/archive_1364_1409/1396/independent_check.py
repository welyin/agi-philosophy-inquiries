"""Independent upward flow with a linear fractional map for the scalar quartic."""
import argparse
import json
from fractions import Fraction
from math import atan2, exp, pi, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "../1369/inputs.json").read_text(encoding="utf-8"))
K = 16*pi*pi


def derivative(state):
    h, w, g, t, b, tau = state[:6]
    trace = 3*(t*t+b*b) + tau*tau
    gauge_yukawa = [41/6*h**3, -19/6*w**3, -7*g**3,
                    t*(trace+1.5*(t*t-b*b)-17/12*h*h-9/4*w*w-8*g*g),
                    b*(trace+1.5*(b*b-t*t)-5/12*h*h-9/4*w*w-8*g*g),
                    tau*(trace+1.5*tau*tau-15/4*h*h-9/4*w*w)]
    linear = 4*trace - 9*w*w - 3*h*h
    constant = 3/8*(2*w**4+(h*h+w*w)**2)-6*(t**4+b**4)-2*tau**4
    a, b12, c, d = state[6:]
    return np.array(gauge_yukawa+[-24*c, -24*d,
                                 constant*a+linear*c, constant*b12+linear*d])/K


def step(state, length):
    k1 = derivative(state)
    k2 = derivative(state+length*k1/2)
    k3 = derivative(state+length*k2/2)
    k4 = derivative(state+length*k3)
    return state+length*(k1+2*k2+2*k3+k4)/6


def crossing(aligned, shared, ds):
    base = np.array([sqrt(3/5)*DATA["g1_GUT"][0], DATA["g2"][0], DATA["g3"][0],
                     DATA["yt"][0], DATA["yb"][0] if aligned else 0,
                     DATA["ytau"][0] if aligned else 0, 1., 0., 0., 1.])
    def residual(state):
        return state[3]**2 + (state[4]**2 if shared else 0) - state[1]**2
    s = 0
    while s < 34:
        advance = step(base, ds)
        if residual(advance) <= 0:
            break
        s += ds
        base = advance
    else:
        raise ValueError("No crossing inside declared domain")
    lo, hi = 0., ds
    for _ in range(48):
        mid = (lo+hi)/2
        if residual(step(base, mid)) > 0:
            lo = mid
        else:
            hi = mid
    event = step(base, (lo+hi)/2)
    h, w, _, t, b, _ = event[:6]
    a, bd, c, d = event[6:]
    assert a*d-bd*c > 0 and abs(residual(event)) < 1e-11
    predicted = -c/d
    observed = .55853/4
    local_sensitivity = (a*d-bd*c)/d**2
    return {"scale_GeV": DATA["mu_ref_GeV"]*exp(s+(lo+hi)/2),
            "lambda0_for_zero_boundary": float(predicted),
            "boundary_lambda_for_observed": float((c+d*observed)/(a+bd*observed)),
            "d_lambda0_d_boundary_at_zero": float(local_sensitivity),
            "hypercharge_mix_sin2": float(h*h/(w*w+h*h)),
            "needed_extra_gX": 1/sqrt(1/h**2-1/(3*w**2)),
            "t_to_w": float(t/w), "b_to_w": float(b/w),
            "angle_t_b_if_shared_radians": atan2(b, t) if shared else None}


def calculate():
    # Fundamental and conjugate pairs must both contain the actual light doublet.
    basic = np.array([Fraction(1,6), Fraction(-1,3)])
    assert tuple(-basic + Fraction(1,3)) == (Fraction(1,6), Fraction(2,3))
    output = {}
    for name, aligned, shared in [("top_only", False, False),
                                  ("three_yukawa_top_boundary", True, False),
                                  ("three_yukawa_shared_doublet", True, True)]:
        small = crossing(aligned, shared, .025)
        large = crossing(aligned, shared, .05)
        assert abs(small["scale_GeV"]/large["scale_GeV"]-1) < 5e-8
        assert abs(small["lambda0_for_zero_boundary"]-large["lambda0_for_zero_boundary"]) < 2e-9
        output[name] = small
    return output


def check(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            check(actual[key], saved[key])
    elif actual is None:
        assert saved is None
    else:
        assert np.isclose(actual, saved, rtol=2e-9, atol=3e-11)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.check:
        check(result, json.loads((HERE/"independent_results.json").read_text(encoding="utf-8")))
        print("1396 independent hypercharge and shared-doublet boundary checks passed")
    else:
        print(json.dumps(result, indent=2, allow_nan=False))
