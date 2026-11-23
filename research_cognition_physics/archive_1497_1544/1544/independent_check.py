"""Cross-check the three leptonic probes of the two neutral eigenstates."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"
VEV = (math.sqrt(2) * 1.1663787e-5) ** -.5


def pole_angle(a):
    return (1 - a / (1 + math.sqrt(1 - a*a))) / 4


def angle_boxes():
    # E158 used the leptonic effective angle at BOTH ends of its transport.
    shift = .2381 - .23149
    shift_margin = 2 * (.0006 + .00015)
    sigma_e = math.hypot(.0010, .0008)
    lo_s, hi_s = [pole_angle(.15056 + sign*2*.00239) for sign in (1, -1)]
    lo_e = 2*(.2397 - 2*sigma_e - (shift + shift_margin) - hi_s)
    hi_e = 2*(.2397 + 2*sigma_e - (shift - shift_margin) - lo_s)
    return {"sld_two_width": [lo_s, hi_s], "e158_central_transfer": shift,
            "transfer_margin_work_box": shift_margin, "E_box": [lo_e, hi_e],
            "E_center": 2*(.2397 - shift - pole_angle(.15056))}


def probe(gb, gr, gl, sqrt_w, q2=0.):
    # Start in the original (WL3, WR3, B) basis; project out the photon.
    gg = np.array([gl, gr, gb])
    photon = 1 / gg
    photon /= np.linalg.norm(photon)
    broken = np.array([[gl, -gr, 0.], [0., gr, -gb]])
    mass2 = broken.T @ np.diag([VEV**2/4, sqrt_w**2]) @ broken
    eig, rot = np.linalg.eigh(mass2)
    assert np.linalg.norm(np.cross(photon, rot[:, 0])) < 1e-10
    # Columns carry T3L, T3R and (B-L)/2 for eL, eR and nuL.
    currents = np.array([[-.5, 0, -.5], [0., -.5, -.5], [.5, 0., -.5]]) * gg
    massive = currents @ rot[:, 1:]
    coeff = VEV**2 * (massive / (eig[1:] + q2)) @ massive.T
    qe = coeff[1, 1] - coeff[0, 0]
    gv = (coeff[2, 0] + coeff[2, 1]) / 2
    ga = (coeff[2, 0] - coeff[2, 1]) / 2
    left, right = massive[:2, 0]
    S = (1 - (left + right)/(left - right))/4
    gy = (1/gr**2 + 1/gb**2)**-.5
    s = gy**2 / (gy**2 + gl**2)
    r, z, d = gr**2/(gr**2 + gb**2), eig[1]/(sqrt_w**2*(gr**2+gb**2)), VEV**2/(4*sqrt_w**2)
    E, L = qe/2 + .5 - 2*S, gv + .5 - 2*S
    if q2 == 0.:
        common = d*(1 - 2*r*(1-S))/(1-z+r*r*d)
        np.testing.assert_allclose([qe, ga, E, L, (1-z)*E - r*(E-L)],
                                   [-1+4*s, -.5, -r*common, (1-z-r)*common, 0], atol=3e-12)
    return {"E": E, "L": L, "S": S, "gV": gv, "gA": ga, "Qe": qe,
            "r": r, "z": z, "d": d, "e_diagonal_sum": coeff[0, 0] + coeff[1, 1],
            "light_neutral_mass2_GeV2": eig[1]}


def run():
    box = angle_boxes()
    gauge = json.loads((HERE.parent/"1369"/"inputs.json").read_text(encoding="utf-8"))
    gy, gl, gb = math.sqrt(.6)*gauge["g1_GUT"][0], gauge["g2"][0], math.sqrt(1.5)*gauge["g3"][0]
    gr = (gy**-2 - gb**-2)**-.5
    finite_shift = {}
    examples = {}
    for scale in (270, 450, 3000):
        zero = probe(gb, gr, gl, scale)
        q2max = 2 * .000511 * 48
        finite = probe(gb, gr, gl, scale, q2max)
        bound = q2max / zero["light_neutral_mass2_GeV2"] * zero["e_diagonal_sum"]
        assert abs(zero["Qe"] - finite["Qe"]) <= bound
        examples[str(scale)] = zero
        finite_shift[str(scale)] = {"delta_Qe_bound": bound, "q2max_GeV2": q2max}
    assert examples["270"]["E"] < box["E_box"][0] < examples["450"]["E"] < 0
    lo_r, hi_r = .055, .065
    lo_s = box["sld_two_width"][0]
    m = -box["E_box"][0]
    alpha = lo_r*(1-2*hi_r*(1-lo_s))
    dmax = m/(alpha-m*hi_r**2)
    # With matched S fixed, a limiting blind family approaches r=0.
    return {"angle_boxes":box, "conditional_r_box": [lo_r, hi_r],
            "conditional_d_upper":dmax, "conditional_sqrtW_lower_GeV":VEV/(2*math.sqrt(dmax)),
            "tree_examples_GeV": examples, "finite_momentum_propagator_bounds": finite_shift}


def compare(now, saved):
    if isinstance(now, dict):
        assert now.keys() == saved.keys()
        for key, value in now.items():
            compare(value, saved[key])
    else:
        np.testing.assert_allclose(now, saved, atol=2e-11, rtol=2e-9)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf-8") as fp:
            json.dump(result, fp, indent=2, allow_nan=False)
            fp.write("\n")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf-8")))
        print("1544 independent: 3 neutral currents, Moller transport and rate windows verified")
