"""A finite SU(5) family, its anomalous currents and the split-carrier alternative."""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def indices(weight):
    c, s, y = weight["color"], weight["weak"], weight["y6"]
    return np.array([s * int(c == 3), c * int(s == 2), c * s * y**2])


def ratio(vec):
    k, w, b = map(int, vec)
    return Fraction(b + 18 * w, 18 * k)


def spectral_radius_mod_scalar(a, b):
    eig_difference = np.sort(a) - np.sort(b)
    return float(np.ptp(eig_difference)) / 2


def compute():
    data = read(HERE / "inputs.json")
    ten = sum(indices(s) for s in data["SM_weights_10"].values())
    five = sum(indices(s) for s in data["SM_weights_5bar"].values())
    q_only = -indices(data["SM_weights_10"]["Q"])
    np.testing.assert_equal(ten, [3, 3, 90])
    np.testing.assert_equal(five, [1, 1, 30])
    prior = read(HERE / data["old_period"])
    np.testing.assert_equal(q_only, prior["texture_K3_K2_K1"])
    partial = data["partial_gamma_witness"]
    families, pairs = partial["Nfamily"], partial["n_Higgs_pairs"]
    assert 2 * families + pairs + 5 == 15
    d2 = Fraction(partial["Ddiag_squared_over_g2"])
    f2 = Fraction(partial["f44_squared_over_g2"])
    matter = [4 * d2 - Fraction(24, 5)] * families
    bar_h = matter + [Fraction(24, 5) * (f2 - 1)]
    assert all(x == 0 for x in matter + bar_h)
    weighted_gamma = -5 * Fraction(24, 5)
    assert weighted_gamma == -partial["adjoint_dim"] * int(five[0])

    q = np.array(data["ordinary_family_spectra"]["Q"], dtype=float)
    u = np.array(data["ordinary_family_spectra"]["uc"], dtype=float)
    mid = (q + u) / 2
    metric = spectral_radius_mod_scalar
    relaxed = max(metric(q, mid), metric(u, mid))
    integer = min(max(metric(q, x), metric(u, x)) for x in itertools.combinations_with_replacement(range(4), 3))
    assert relaxed == metric(q, u) / 2 == 0.25 and integer == 0.5

    # Four parent tens; first Q in a charged carrier, first u^c in a different neutral one.
    qmap = np.eye(4)[:, [0, 1, 2]]
    umap = np.eye(4)[:, [1, 2, 3]]
    current = umap.conj().T @ qmap
    parent = np.diag([-1., 0., 0., 0.])
    np.testing.assert_equal(umap.T @ parent @ umap @ current, current @ qmap.T @ parent @ qmap)
    principal_values = np.linalg.svd(current, compute_uv=False)
    np.testing.assert_equal(principal_values, [1, 1, 0])

    calib = read(HERE / data["old_QCD"])
    widths = calib["window_widths"]
    c, dc = calib["photon_lattice"]
    cmin, cmax = c - widths * dc, c + widths * dc
    t, dt1, dt2 = calib["chisus_4throot_MeV"]
    lo, hi = (t - widths * (dt1 + dt2)) / 1000, (t + widths * (dt1 + dt2)) / 1000
    pref = 1 / (2 * math.pi * calib["alpha_inverse"])
    lower_norm, upper_norm = pref / hi**2, pref / lo**2
    qratio, gratio = float(ratio(q_only)), float(ratio(ten))
    assert qratio < cmin < cmax < gratio
    qsl = [(cmin - qratio) * lower_norm, (cmax - qratio) * upper_norm]
    gsl = [(gratio - cmax) * lower_norm, (gratio - cmin) * upper_norm]
    np.testing.assert_allclose(qsl, prior["unsigned_slope_GeV_minus2"]["target_5_over_3"])
    gap_ind = gsl[0] - qsl[1]
    gap_shared = lower_norm * (gratio + qratio - 2 * cmax)
    return {
        "complete_ten": ten.tolist(), "complete_five": five.tolist(), "Q1": q_only.tolist(),
        "r_complete": str(ratio(ten)), "r_Q1": str(ratio(q_only)), "r_defect": str(ratio(ten) - ratio(q_only)),
        "Q1_defects_w_minus_k_b_minus_30k": [int(q_only[1] - q_only[0]), int(q_only[2] - 30 * q_only[0])],
        "partial_finiteness_H1_axion": {"gamma_bar5_matter_and_barH_times_16pi2_over_g2": list(map(float, matter + bar_h)),
                                       "PQ_weighted_dim_trace_gamma_times_16pi2_over_g2": float(weighted_gamma),
                                       "k_w_b": five.tolist()},
        "flavour_minimax_opnorm_mod_scalar": {"continuous": relaxed, "unit_period": float(integer)},
        "split_tens_current_singular_values": principal_values.tolist(),
        "qcd_slope_GeV_minus2": {"Q1": qsl, "standard_SU5_anomalous_branch": gsl,
                                 "separate_calibrations_gap": gap_ind, "shared_calibration_min_gap": gap_shared,
                                 "separate_calibrations_per_hypothesis_additive_halfwidth_threshold": gap_ind / 2,
                                 "separate_calibrations_per_hypothesis_relative_halfwidth_threshold": gap_ind / (gsl[0] + qsl[1]),
                                 "shared_calibration_per_hypothesis_relative_halfwidth_threshold": gratio + qratio - 2 * cmax,
                                 "shared_calibration_slope_ratio": [(gratio - cmax)/(cmax - qratio),
                                                                    (gratio - cmin)/(cmin - qratio)]},
    }


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-13), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    mode = cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = cli.parse_args()
    path = HERE / "results.json"
    if args.write:
        with path.open("x", encoding="utf8") as out:
            json.dump(compute(), out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(compute(), read(path))
        print("1512 finite trace, compact charges, split carriers and QCD response verified.")
