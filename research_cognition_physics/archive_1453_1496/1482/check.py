"""Replay the two-ordering, inclusive lepton-current certificate."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def null_direction(order, a, b, c, delta):
    """a,b,c denote sin²(theta_12,23,13)."""
    phase = np.exp(1j * delta)
    if order == "IO":
        return np.array([math.sqrt(c) / phase, math.sqrt(b * (1-c)),
                         math.sqrt((1-b) * (1-c))])
    return np.array([math.sqrt((1-a)*(1-c)),
                     -math.sqrt(a*(1-b)) - math.sqrt((1-a)*b*c)*phase,
                     math.sqrt(a*b) - math.sqrt((1-a)*(1-b)*c)*phase])


def rate(x, i, j, s):
    return 1 + 2*s*x[3-i-j] + s*s*(1-x[i])*(1-x[j])


def audit(n, s):
    x = np.abs(n)**2
    f = np.array([[0, n[2], -n[1]], [-n[2], 0, n[0]], [n[1], -n[0], 0]])
    assert np.allclose(f @ n, 0)
    assert np.isclose(abs((f.conj().T @ f)[0, 1])**2, x[0]*x[1])
    for i, j in [(2, 1), (2, 0), (1, 0)]:
        current = -s * np.outer(f[j].conj(), f[i])
        current[i, j] += 1
        assert np.isclose(np.linalg.norm(current)**2, rate(x, i, j, s))
    return math.sqrt(rate(x, 2, 1, s) / rate(x, 2, 0, s))


def certificate():
    cfg = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    osc = json.loads((HERE / cfg["oscillations"]).read_text(encoding="utf-8"))
    amplitude, err = osc["DayaBay_amplitude"]
    widths = cfg["oscillation_width_count"]
    c_lo, c_hi = [(1-math.sqrt(1-A))/2 for A in (amplitude-widths*err, amplitude+widths*err)]
    b_tau = cfg["Belle_II"]
    tau_upper = b_tau["ratio_central"] + b_tau["ratio_width"]
    result = {"theta13_sine_squared_box": [c_lo, c_hi]}
    for order in ["NO", "IO"]:
        p = osc[order]
        a, ea = p["solar_sine_squared"]
        a_lo, a_hi = a-widths*ea, a+widths*ea
        b, eb = p["atmospheric_sine_squared"][:2]
        b_lo = b - widths*eb
        b_hi = b + widths*p["atmospheric_sine_squared"][-1]
        if order == "NO":
            e_lo, e_hi = (1-a_hi)*(1-c_hi), (1-a_lo)*(1-c_lo)
            positive = math.sqrt(a_lo*(1-b_hi)) - math.sqrt((1-a_lo)*b_hi*c_hi)
            m_lo, m_hi = max(0, positive)**2, 1-e_lo
        else:
            e_lo, e_hi = c_lo, c_hi
            m_lo, m_hi = b_lo*(1-c_hi), b_hi*(1-c_lo)
        q = math.sqrt(24*math.pi * cfg["MEG_II"]["BR_90CL_upper"] *
                      cfg["inverse_alpha_upper"] / (cfg["matched_dipole_rate_lower_factor"]*e_lo*m_lo))
        s_cap = q/(1-q)
        gap = max(abs(e_hi-m_lo), abs(m_hi-e_lo))
        cap = gap*(s_cap+s_cap*s_cap/2)
        ratio_interval = [1, 1+cap] if order == "NO" else [1-cap, 1]
        s_witness = cfg["illustrative_contact_strength"]
        n = null_direction(order, a, b, (c_hi+c_lo)/2, math.pi/3)
        x = np.abs(n)**2
        ratio = audit(n, s_witness)
        for delta in np.linspace(0, 2*math.pi, 17):
            for theta in [(a_lo,b_hi,c_hi), (a_hi,b_lo,c_lo)]:
                for strength in (0., s_cap, .3, 10.):
                    audit(null_direction(order, *theta, delta), strength)
        assert ratio_interval[0] <= ratio <= ratio_interval[1]
        assert ratio_interval[0] > tau_upper
        result[order] = {
            "x_e_outer": [e_lo, e_hi], "x_mu_outer": [m_lo, m_hi],
            "contact_upper": s_cap, "tau_ratio_outer": ratio_interval,
            "minimum_relative_downshift_to_Belle_one_width": 1-tau_upper/ratio_interval[0],
            "witness_tau_ratio": ratio,
            "witness_h_dipole_over_inclusive_muon_contact_rate":
                s_witness**2*x[0]*x[1]/(137.036*24*math.pi*rate(x,1,0,s_witness))}
    return result


def check_saved(actual, saved):
    if isinstance(actual, dict):
        return actual.keys() == saved.keys() and all(check_saved(v, saved[k]) for k, v in actual.items())
    return np.allclose(actual, saved, rtol=1e-11, atol=1e-18)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--write", action="store_true")
    options.add_argument("--check", action="store_true")
    result = certificate()
    target = HERE / "results.json"
    if parser.parse_args().write:
        with target.open("x", encoding="utf-8") as output:
            json.dump(result, output, indent=2, allow_nan=False)
            output.write("\n")
    else:
        assert check_saved(result, json.loads(target.read_text(encoding="utf-8")))
    print(json.dumps(result, indent=2))
