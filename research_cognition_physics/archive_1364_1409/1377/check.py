"""Higgs kinetic saturation, sequential colored families and a width-independent readout."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def inherited_spectrum():
    path = HERE.parent / "1373/check.py"
    spec = importlib.util.spec_from_file_location("round1373_scalar", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _, _, fermions = mod.common_inputs(mod.read("inputs.json"))
    return mod.loop_functions, [mass for mass, colors, _ in fermions if colors == 3]


def evaluate(data):
    forms, masses = inherited_spectrum()
    mh = data["on_shell_h_GeV"]
    amplitude = sum(forms(mh, m)[1] for m in masses)
    assert len(masses) == 6 and amplitude.real > 0
    enhancement = lambda shift: abs(amplitude + shift)**2 / abs(amplitude)**2
    min_pure = enhancement(8/3)
    scan = {str(m): enhancement(2*forms(mh, m)[1])
            for m in data["sequential_heavy_tests_GeV"]}
    assert all(v > min_pure for v in scan.values())
    obs = data["atlas_four_leptons"]
    norm = obs["central_normalization"]
    bb_cap = (norm["inclusive_bbh_pb"] * norm["br_zz"]
              / norm["stage_gg_plus_bb_times_br_pb"])
    min_mixed = (1-bb_cap)*min_pure + bb_cap
    g, w = obs["gg_plus_bb"], obs["vbf"]
    lims = {}
    for widths in (1, 2):
        r = (g["mu"] + widths * math.hypot(*g["width_plus"])) / (
            w["mu"] - widths * math.hypot(*w["width_minus"]))
        upper_pure = (r - bb_cap) / (1 - bb_cap)
        shift = math.sqrt(upper_pure * abs(amplitude)**2 - amplitude.imag**2) - amplitude.real
        lims[str(widths)] = {"upper_ratio": r, "upper_sum_nonnegative_heavy_kappa": 3*shift/4}
        assert r < min_mixed
    # A second Higgs doublet carries the extra family mass.
    y = data["synthetic_septet_fraction"]
    z, x = math.sqrt(15)*y, math.sqrt(1-16*y*y)
    vev = np.array([x, z, y/4])
    light = np.array([x, 0, 4*y])
    weights = np.array([1, 1, 16])
    k_v = (weights * vev) @ light
    k_2v = weights @ (light**2)
    deficit = np.sum(weights * (light-k_v*vev)**2)
    np.testing.assert_allclose([np.linalg.norm(light), weights@(vev**2), k_v], 1)
    np.testing.assert_allclose([k_2v - 1, deficit], 240*y*y)
    assert deficit >= z*z/(1-z*z)
    # Test exact mass-matrix alignment through cancellations and off-diagonal Yukawa entries.
    y1 = np.array([[1, 2j], [0.3, -0.7]])
    y2 = np.array([[-1, -1j], [0.1j, 0.2]])
    h_sat = np.r_[vev, math.sqrt(1-vev@vev)]
    w_ext = np.r_[weights, 0.]
    vev_ext = np.r_[vev, 0.]
    np.testing.assert_allclose([w_ext @ h_sat**2, (w_ext*vev_ext)@h_sat], 1)
    np.testing.assert_allclose(y1*h_sat[0] + y2*h_sat[1], y1*x + y2*z)
    return {"SM_quark_loop_real_imag": [amplitude.real, amplitude.imag],
            "LO_pure_gg_min_one_extra_family": min_pure, "LO_mass_scan": scan,
            "central_upper_bb_fraction": bb_cap, "min_gg_plus_bb_over_vbf_ratio": min_mixed,
            "atlas_central_ratio": g["mu"]/w["mu"], "marginal_rectangles": lims,
            "synthetic_shared_h_without_saturation": {"vev_fraction_second_doublet": z,
                "k_known_f": light[0]/x, "k_extra_f": light[1]/z,
                "k_W": k_v, "k_2W": k_2v}}


def compare(a, b):
    if isinstance(a, dict):
        assert set(a) == set(b)
        for key in a:
            compare(a[key], b[key])
    else:
        np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-11)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = evaluate(json.loads((HERE/"inputs.json").read_text(encoding="utf-8")))
    target = HERE / "results.json"
    if args.write:
        with target.open("x", encoding="utf-8") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(result, json.loads(target.read_text(encoding="utf-8")))
        print("1377: matrix alignment, extra-family envelope and the ATLAS ratio verified.")
