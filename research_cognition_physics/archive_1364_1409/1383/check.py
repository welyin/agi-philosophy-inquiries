"""Historical CKM row and independent pion/tau charged-current measurements."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def interval(values, width):
    radius = width * sum(values[1:])
    return values[0]-radius, values[0]+radius


def upper_ub(old, width):
    rate_hi = interval(old["F_times_cb"], width)[1]
    ratio_hi = interval(old["ub_over_cb"], width)[1]
    return rate_hi * ratio_hi / old["F_lower"]


def kaon_weights(kaons):
    errors = np.array(kaons["width"])
    covariance = np.outer(errors, errors) * np.array(kaons["correlation"])
    raw = np.linalg.solve(covariance, np.ones(len(errors)))
    weights = raw / raw.sum()
    assert np.min(weights) > 0 and np.linalg.eigvalsh(covariance).min() > 0
    return {"weights": weights.tolist(),
            "mu_weight": sum(weights[i] for i in (1,4)),
            "rounded_table_mean": float(weights @ kaons["product"]),
            "rounded_table_error": float(np.sqrt(weights @ covariance @ weights))}


def pion_interval(pion, width):
    obs = interval(pion["exp_in_units_1e_4"], width)
    sm = interval(pion["same_paper_SM"], width)
    return obs[0]/sm[1], obs[1]/sm[0]


def minimal_shifts(smax):
    return {"mu_amplitude_increase": max(0.0, 1/math.sqrt(smax)-1),
            "common_quark_amplitude_reduction": max(0.0, 1-math.sqrt(smax))}


def evaluate():
    cfg = read(ROOT / "inputs.json")
    old = read(ROOT / cfg["inherited"])
    kaons, pion = cfg["KLOE"], cfg["PIENU"]
    tau = cfg["HFLAV_tau"]["pure_leptonic_eq_6_effective_g_tau_over_g_mu"]
    u = old["ud_beta"][0]
    f = old["K_form"][0]
    ub_central = old["F_times_cb"][0] * old["ub_over_cb"][0] / old["F_example"]
    s0 = u*u + (kaons["official_average"][0]/f)**2 + ub_central**2
    root_s = math.sqrt(s0)
    central_e_sum = u*u + (kaons["product"][0]/f)**2 + ub_central**2
    out = {"kloe": kaon_weights(kaons), "row_nominal": s0,
           "nominal_minimal_changes": minimal_shifts(s0),
           "pure_mu_contact_tau_ratio": root_s, "pure_quark_vertex_tau_ratio": 1.0,
           "pion_ratio_central": pion_interval(pion, 0)[0], "bounds": {},
           "KL_e_nominal_complete_electron_row": central_e_sum}
    for width in (1,2):
        dlow, dhigh = interval(old["ud_beta"], width)
        flow, fhigh = interval(old["K_form"], width)
        bhigh = upper_ub(old, width)
        pion_box = pion_interval(pion, width)
        # An unresolved e/mu rate behaves between the pure-tag endpoints.
        electron_b_high = bhigh * max(1, math.sqrt(pion_box[1]))
        thresholds = {}
        channels = {"official_mixed": kaons["official_average"],
                    "KL_e": (kaons["product"][0], kaons["width"][0]),
                    "KL_mu": (kaons["product"][1], kaons["width"][1])}
        for label, values in channels.items():
            klow, khigh = interval(values, width)
            smin = dlow**2 + (klow/fhigh)**2
            smax = dhigh**2 + (khigh/flow)**2 + bhigh**2
            thresholds[label] = {"row": [smin, smax], "minimal_changes": minimal_shifts(smax)}
        k_high = interval(channels["KL_e"], width)[1]
        s_e_upper = dhigh**2 + (k_high/flow)**2 + electron_b_high**2
        t_low = interval(tau, width)[0]
        thresholds["tau_vs_complete_electronic_row"] = {
            "pion_Ae_over_Amu": pion_box, "nominal_b_cap": bhigh,
            "electronic_b_cap": electron_b_high, "R_tau_upper_if_q_one": math.sqrt(s_e_upper),
            "R_tau_lower_data": t_low, "q_upper_if_Atau_at_most_one": math.sqrt(s_e_upper)/t_low,
            "q_one_intervals_disjoint": bool(t_low > math.sqrt(s_e_upper))}
        out["bounds"][str(width)] = thresholds
    # Two different low-energy theories giving exactly the same nominal semileptonic readouts.
    true_v = np.array([u, kaons["official_average"][0]/f, ub_central]) / root_s
    assert abs(true_v @ true_v - 1) < 1e-14
    for quark, mu in ((root_s, 1.0), (1, 1/root_s)):
        assert abs((quark/mu)**2 * np.dot(true_v, true_v) - s0) < 1e-14
    # Existence with two accessible deficits equal and tau's current unit-normalized.
    # Uses the mean kaon product, not an exact refit of the five observed central values.
    common_a = 1/tau[0]**2
    q = root_s/tau[0]
    effective = q / math.sqrt(common_a)
    assert abs(effective - root_s) < 1e-14
    out["example_with_tau_central_and_pion_universality"] = {
        "Ae_equals_Amu": common_a, "Atau": 1.0, "mu_contact": 0.0,
        "common_quark_amplitude_q": q, "unit_true_quark_row": true_v.tolist(),
        "effective_tau_ratio": 1/math.sqrt(common_a)}
    return out


def same(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, float):
        return math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-14)
    return a == b


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    outcome = evaluate()
    path = ROOT / "results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(outcome, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
    else:
        assert same(outcome, read(path))
    print("1383 OK")
