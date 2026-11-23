"""Recompute joint beta-kinematic and light-exchange certificates in eV."""
import argparse
import json
from math import acos, degrees, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def offsets(order, solar, atm):
    if order == "NO":
        return np.array([0.0, solar, atm])
    return np.array([atm-solar, atm, 0.0])


def prediction(order, mass, params):
    s, t, solar, atm = params
    p = np.array([(1-s)*(1-t), s*(1-t), t])
    masses = np.sqrt(mass**2 + offsets(order, solar, atm))
    lengths = p * masses
    beta = sqrt(p @ masses**2)
    low = max(0.0, 2*max(lengths) - sum(lengths))
    return beta, [float(low), float(sum(lengths))], lengths


def zero_window(params):
    def boundaries(x):
        _, _, (a, b, c) = prediction("NO", x, params)
        return a+c-b, a-b-c

    def root(which):
        lo, hi = 0.0, 0.05
        assert boundaries(lo)[which] < 0 < boundaries(hi)[which]
        for _ in range(60):
            mid = (lo+hi)/2
            if boundaries(mid)[which] < 0:
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2

    return [root(0), root(1)]


def fixed_osc_distance_check(order, params, x):
    s, t, solar, atm = params
    weights = np.array([(1-s)*(1-t), s*(1-t), t])
    row = np.sqrt(weights)
    second = np.cross(row, [0., 0., 1.])
    second /= np.linalg.norm(second)
    u = np.array([row, second, np.cross(row, second)]) @ np.diag(np.exp(1j*np.array([0, .64, 1.07])))
    m0 = np.sqrt(offsets(order, solar, atm))
    target = u.conj() @ np.diag(np.sqrt(x**2 + m0**2)) @ u.conj().T
    base = u.conj() @ np.diag(m0) @ u.conj().T
    assert np.linalg.matrix_rank(base, tol=1e-12) == 2
    assert np.allclose(np.diff(np.linalg.svd(base, compute_uv=False)**2),
                       np.diff(np.linalg.svd(target, compute_uv=False)**2), atol=1e-13)
    assert np.isclose(np.linalg.norm(target-base, 2), x, atol=1e-12)
    assert np.isclose(np.linalg.svd(target, compute_uv=False)[-1], x, atol=1e-12)


def evaluate():
    inputs = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    data = {}
    for order in ["NO", "IO"]:
        raw = inputs[order]
        s, t, delta, atm = [raw[k] for k in ["s12", "s13", "solar", "atm_abs"]]
        center = [s[0], t[0], delta[0], atm[0]]
        beta0, bb0, terms = prediction(order, 0., center)
        katrin = inputs["katrin"]["m_beta_upper_eV"]
        entry = {"m_beta_rank2": beta0, "m_bb_rank2_window": bb0,
                 "permitted_lightest_under_KATRIN_central_osc": sqrt(katrin**2-beta0**2)}
        for x in [0.001, 0.025, 0.112]:
            fixed_osc_distance_check(order, center, x)
        if order == "NO":
            win = zero_window(center)
            # Lower root increases with s12 and solar and decreases with s13 and atm.
            low_lo = zero_window([s[1], t[2], delta[1], atm[2]])[0]
            low_hi = zero_window([s[2], t[1], delta[2], atm[1]])[0]
            high_lo = zero_window([s[1], t[1], delta[1], atm[1]])[1]
            high_hi = zero_window([s[2], t[2], delta[2], atm[2]])[1]
            entry.update({"cancel_window_m0":win, "cancel_window_m_beta":[
                prediction(order, x, center)[0] for x in win],
                "lower_cancel_edge_over_box":[low_lo,low_hi],
                "upper_cancel_edge_over_box":[high_lo,high_hi]})
            # Explicit noncollinear phase witness halfway through the cancellation interval.
            _, _, (a, b, c) = prediction(order, sum(win)/2, center)
            phi = acos((c*c-a*a-b*b)/(2*a*b))
            third = -(a+b*np.exp(1j*phi))/c
            assert abs(abs(third)-1) < 1e-12
            assert abs(a+b*np.exp(1j*phi)+c*third) < 1e-13
        else:
            gamma = sqrt(atm[1]/(atm[1]-delta[2]))
            coeff = (1-t[2])*(1-s[2]*(1+gamma)) - t[2]
            assert coeff > 0
            a, b = terms[:2]
            limit = inputs["kamland_zen"]["m_bb_upper_eV_endpoints_by_nme"][0]
            cosmax = (limit**2 - a*a-b*b)/(2*a*b)
            entry.update({"all_x_box_floor_certificate": coeff * sqrt(atm[1]-delta[2]),
                          "phase_for_28meV_rank2_IO": {
                              "cos_max": cosmax, "minimum_absolute_phase_deg": degrees(acos(cosmax))}})
        data[order] = entry
    return data


def compare(x, y):
    if isinstance(x, dict):
        assert x.keys() == y.keys()
        for k in x:
            compare(x[k], y[k])
    else:
        np.testing.assert_allclose(x, y, rtol=1e-10, atol=1e-13)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    path = HERE/"results.json"
    result = evaluate()
    if args.write_new:
        with path.open("x", encoding="utf-8") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
    print("1370 matched-oscillation distance, phase bands and cancellation checks passed.")
