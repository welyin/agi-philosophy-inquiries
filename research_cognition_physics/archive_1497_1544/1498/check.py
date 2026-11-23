"""Uniform four-meson flavor-ray certificates and explicit blind directions."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
PAIRS = ((0, 1), (0, 2), (1, 2))


def read(name):
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def interval(parts, n):
    center, *widths = parts
    return center - n * sum(widths), center + n * sum(widths)


def outer_moduli(inputs):
    tree = read(inputs["charged_current_inputs"])
    n = inputs["half_widths_of_sum"]
    product = interval(tree["K_product"], n)
    form = interval(tree["K_form"], n)
    us = [product[0] / form[1], product[1] / form[0]]
    cb = interval(read(inputs["belle_bgl"])["BelleII_Vcb_BGL"], n)
    ratio = interval(tree["ub_over_cb"], n)
    ub = [ratio[0] * cb[0], ratio[1] * cb[1]]
    c13min = math.sqrt(1 - ub[1]**2)
    s12max, s23max = us[1] / c13min, cb[1] / c13min
    c12min, c23min = math.sqrt(1-s12max**2), math.sqrt(1-s23max**2)
    lower = np.array([
        [math.sqrt(1 - us[1]**2 - ub[1]**2), us[0], ub[0]],
        [us[0] * c23min - s23max * ub[1],
         c12min * c23min - s12max * s23max * ub[1], cb[0]],
    ])
    upper = np.array([[1, us[1], ub[1]], [s12max + s23max*ub[1], 1, cb[1]]])
    assert np.min(lower) > 0 and np.all(upper >= lower)
    return lower, upper


def ckm(s12, s13, s23, phase):
    c12, c13, c23 = [math.sqrt(1-s*s) for s in (s12, s13, s23)]
    z = np.exp(1j * phase)
    return np.array([
        [c12*c13, s12*c13, s13*z.conjugate()],
        [-s12*c23-c12*s23*s13*z, c12*c23-s12*s23*s13*z, s23*c13],
        [s12*s23-c12*c23*s13*z, -c12*s23-s12*c23*s13*z, c23*c13],
    ])


def dominant_norm(eps):
    return math.sqrt((1 + math.sqrt(1 - 8*eps**2)) / 2)


def certified_floor(gap):
    lo, hi = 0.0, 0.01
    assert gap(lo) > 0 and gap(hi) < 0
    for _ in range(60):
        mid = (lo + hi) / 2
        if gap(mid) > 0:
            lo = mid
        else:
            hi = mid
    # Return a conservative decimal strictly inside the excluded interval.
    value = (math.floor(lo * 10**12) - 1) / 10**12
    assert gap(value) > 0
    return value


def unrestricted_gap(eps, lower, upper):
    t = dominant_norm(eps)
    floors = []
    for j in range(3):
        rest = [i for i in range(3) if i != j]
        bound = np.maximum(0., lower[:, j]*t - eps/t * upper[:, rest].sum(axis=1))
        floors.append(float(np.prod(bound)))
    return min(floors) - eps


def linear_terms(v, j):
    rest = [i for i in range(3) if i != j]
    z = v[0, j] * v[1, j].conjugate()
    phase = z.conjugate() / abs(z)
    first = phase * v[0, rest] * v[1, j].conjugate()
    second = (phase * v[0, j] * v[1, rest].conjugate()).conjugate()
    return rest, z, first + second


def fixed_matrix_gap(eps, v):
    t = dominant_norm(eps)
    floors = []
    for j in range(3):
        rest, z, slopes = linear_terms(v, j)
        quadratic = np.prod(np.abs(v[:2, rest]).sum(axis=1))
        floors.append(t*t*abs(z) - eps * sum(abs(slopes)) - (eps/t)**2 * quadratic)
    return min(floors) - eps


def channels(v, d):
    u = v @ d
    return [float(abs(d[i]*d[j])) for i, j in PAIRS] + [float(abs(u[0]*u[1]))]


def witness(v):
    rest, _, slopes = linear_terms(v, 2)
    def direction(q):
        d = np.zeros(3, complex)
        d[2] = math.sqrt(1 - 2*q*q)
        d[rest] = -q * slopes.conjugate() / np.abs(slopes)
        return d
    lo, hi = 0., 0.001
    for _ in range(60):
        mid = (lo + hi) / 2
        terms = channels(v, direction(mid))
        if terms[-1] > max(terms[:-1]):
            lo = mid
        else:
            hi = mid
    d = direction(hi)
    assert np.isclose(np.vdot(d, d), 1)
    return {"down_line_re_im": [[z.real, z.imag] for z in d],
            "channel_moduli": channels(v, d), "upper_bound": max(channels(v, d))}


def certificate():
    inputs = read("inputs.json")
    lower, upper = outer_moduli(inputs)
    tree_floor = certified_floor(lambda e: unrestricted_gap(e, lower, upper))
    angles = read(inputs["standard_model_reference"])["quark_sines_and_phase"]
    v = ckm(*[pair[0] for pair in angles])
    assert np.allclose(v @ v.conjugate().T, np.eye(3))
    central_floor = certified_floor(lambda e: fixed_matrix_gap(e, v))
    point = witness(v)
    assert central_floor < point["upper_bound"]
    blind = {}
    labels = ["K", "Bd", "Bs", "D"]
    for channel, (i, j) in enumerate(PAIRS):
        d = np.zeros(3, complex)
        d[i], d[j] = v[0, j], -v[0, i]
        d /= np.linalg.norm(d)
        amplitudes = channels(v, d)
        assert amplitudes[channel] > 0 and max(amplitudes[:channel]+amplitudes[channel+1:]) < 1e-14
        blind[labels[channel]] = amplitudes
    blind["D"] = channels(v, np.array([0.,0.,1.]))
    return {"tree_only_moduli_outer_low": lower.tolist(), "tree_only_moduli_outer_high": upper.tolist(),
            "tree_all_directions_strict_lower": tree_floor,
            "tree_max_abs_C_over_F0_lower": tree_floor**2 / 2,
            "matched_reference_strict_lower": central_floor, "matched_reference_witness": point,
            "central_width": point["upper_bound"]-central_floor, "one_channel_only": blind}


def equal(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(equal(x, b[k]) for k, x in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(equal(x, y) for x, y in zip(a, b))
    return math.isclose(a, b, abs_tol=1e-15, rel_tol=1e-10)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = p.parse_args()
    cert = certificate()
    path = ROOT / "results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as out:
            json.dump(cert, out, indent=2)
            out.write("\n")
    elif not equal(cert, read("results.json")):
        raise ValueError("Four-channel certificate differs")
    print(json.dumps(cert, indent=2))
