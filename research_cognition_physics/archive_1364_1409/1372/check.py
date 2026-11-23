"""Reproduce the leading single-flat-direction spectrum; input and masses are in GeV."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
LOOP = 8 * math.pi**2


def read(name):
    with (HERE / name).open(encoding="utf-8-sig") as handle:
        return json.load(handle)


def matching_couplings():
    wide = read("../1369/inputs.json")
    up = read("../1368/inputs.json")["MZ"]
    light = read("../1367/inputs.json")
    values = {key: wide[key] for key in ["g1_GUT", "g2", "yb", "yt", "ytau"]}
    values.update(yu=up["up"], yc=up["charm"])
    point = light["points"]["MZ"]
    for i, label in enumerate(light["order"]):
        values[label] = [point[field][i] * light["units"][i]
                         for field in ["central", "half_width"]]
    return values


def known_trace(values, corner):
    # For corner=+1 maximise the signed supertrace over the whole marginal box.
    def value(key):
        center, width = values[key]
        orientation = 1 if key.startswith("g") else -1
        return center + corner * orientation * width

    g2 = value("g2")
    gy = math.sqrt(3 / 5) * value("g1_GUT")
    vectors = (6 * g2**4 + 3 * (g2**2 + gy**2)**2) / 16
    quarks = 3 * sum(value(k)**4 for k in ("yu", "yc", "yt", "yd", "ys", "yb"))
    leptons = sum(value(k)**4 for k in ("ye", "ymu", "ytau"))
    return vectors - quarks - leptons


def check_block_trace():
    # Genuinely complex rectangular Dirac and symmetric Majorana blocks.
    d = np.array([[1+2j, 2], [-1j, 3+1j], [2-2j, -1]], dtype=complex)
    r = np.array([[3+2j, -1+3j], [-1+3j, 1j]])
    mass = np.block([[np.zeros((3, 3)), d], [d.T, r]])
    direct = float(sum(np.linalg.svd(mass, compute_uv=False)**4))
    h = d.conj().T @ d
    k = r.conj().T @ r
    expansion = float(np.trace(2*h@h + 4*h@k + k@k).real)
    assert math.isclose(expansion, direct, rel_tol=1e-13)
    return expansion


def calculate():
    inputs = read("inputs.json")
    vals = matching_couplings()
    low, center, high = [known_trace(vals, c) for c in (-1, 0, 1)]
    assert low < center < high < 0
    v, dv = inputs["matching"]["v_GeV"]
    m_h, dm_h = inputs["atlas"]["mass_GeV"]
    neutrino = read(inputs["neutrino"]["source"])["NO"]
    light = [math.sqrt(neutrino[key][0]) * 1e-9 for key in ("solar", "atm_abs")]
    setup = inputs["competition"]
    heavy = setup["heavy_over_v"] * v
    f = math.hypot(v, setup["singlet_vev_over_v"] * v)
    light_moment = 2 * sum(m**4 for m in light)
    neutral = light_moment + 4 * heavy**4
    for mass in light:
        block = np.array([[0, math.sqrt(mass * heavy)],
                          [math.sqrt(mass * heavy), heavy - mass]])
        assert np.allclose(np.linalg.eigvalsh(block), [-mass, heavy],
                           rtol=1e-12, atol=1e-24)

    witnesses = []
    for coupling in setup["lambda"]:
        transverse_mass = 2 * math.sqrt(coupling) * v
        curvature = (center*v**4 + transverse_mass**4 - neutral) / (LOOP*f**2)
        witnesses.append({"lambda": coupling, "transverse_mass_GeV": transverse_mass,
                          "curvature_GeV2": curvature,
                          "radial_mass_if_positive_GeV":
                              math.sqrt(curvature) if curvature > 0 else None})

    def roots(vv, deficit, mh):
        return {"stability_boundary_GeV": vv*deficit**0.25,
                "Higgs_matching_boundary_GeV":
                    (deficit*vv**4 + LOOP*vv**2*mh**2)**0.25}

    old = inputs["historical_check"]
    vector_bound = 6*old["upper_W_GeV"]**4 + 3*old["upper_Z_GeV"]**4
    historic_upper = vector_bound - 12*old["lower_top_GeV"]**4
    assert historic_upper < 0
    return {
        "SM_trace_over_v4": {"center": center, "marginal_box": [low, high]},
        "pre_Higgs_direct_sign": {"top_threshold_GeV": (vector_bound/12)**0.25,
                                  "upper_trace_proxy_GeV4": historic_upper},
        "bosonic_fourth_root_with_f_equal_v_no_neutral": {
            "center": roots(v, -center, m_h),
            "box_low": roots(v-dv, -high, m_h-dm_h),
            "box_high": roots(v+dv, -low, m_h+dm_h)},
        "test_of_complex_block_trace": check_block_trace(),
        "doublet_singlet": {"f_GeV": f, "kappa_V": v/f, "heavy_each_GeV": heavy,
                            "light_GeV": light, "neutral_fourth_moment_over_v4": neutral/v**4,
                            "boson_fourth_root_for_ATLAS_GeV":
                                (neutral-center*v**4+LOOP*f**2*m_h**2)**0.25,
                            "witnesses": witnesses}}


def compare(expected, actual):
    if isinstance(expected, dict):
        assert expected.keys() == actual.keys()
        for key in expected:
            compare(expected[key], actual[key])
    elif isinstance(expected, list):
        assert len(expected) == len(actual)
        for left, right in zip(expected, actual):
            compare(left, right)
    elif isinstance(expected, (int, float)):
        assert math.isclose(expected, actual, rel_tol=5e-12, abs_tol=1e-27)
    else:
        assert expected == actual


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    opts = parser.add_mutually_exclusive_group()
    opts.add_argument("--write-new", action="store_true")
    opts.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_new:
        with (HERE / "results.json").open("x", encoding="utf-8") as output:
            json.dump(result, output, indent=2, allow_nan=False)
            output.write("\n")
    elif args.check:
        compare(read("results.json"), result)
        print("1372: spectral bounds, full Majorana trace and scalar witnesses verified")
    else:
        print(json.dumps(result, indent=2, allow_nan=False))
