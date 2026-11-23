"""One right-electron target, a tagged neutrino column, and its charged partners."""
import argparse
import json
import math
from itertools import product
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "inputs.json").read_text(encoding="utf-8"))


def box(pair, widths):
    return pair[0] - widths * pair[1], pair[0] + widths * pair[1]


def inputs_for_order(solar):
    case = DATA["residual"]
    p_mu = math.cos(math.pi * case["irrep_j"] / case["d"]) ** 2
    corners = product(box(case["fast"], case["widths"]), box(solar, case["widths"]))
    p_e = [(1 + math.sqrt(1 - fast)) * (1 - slow) / 2 for fast, slow in corners]
    return p_mu, [min(p_e), max(p_e)]


def bound(p_mu, electron, eta, response_mu, response_tau):
    # eta = accepted fraction / fraction assumed in the published muon limit.
    ratio_limit = DATA["mu3e_ratio_90"] / (eta * response_mu)
    tau = DATA["ordinary_tau_e"]
    tau_e = box([tau["value"], tau["combined"]], DATA["residual"]["widths"])
    p_min, p_max = electron
    c_max = math.sqrt(ratio_limit / (p_min * p_mu))
    ratio = [tau_e[0] * (1 - p_max - p_mu) / p_mu,
             tau_e[1] * (1 - p_min - p_mu) / p_mu]
    def diagonal_difference(p):
        return math.sqrt(ratio_limit / (p * p_mu)) * abs(p + 2 * p_mu - 1)
    mu_tau = max(diagonal_difference(p) for p in electron)
    mu_e = (p_max - p_mu) * math.sqrt(ratio_limit / (p_max * p_mu))
    tau_upper = DATA["mu3e_ratio_90"] * ratio[1] * response_tau / (eta * response_mu)
    return {
        "common_splitting_abs_max": c_max,
        "mu_minus_tau_diagonal_abs_max": mu_tau,
        "mu_minus_e_diagonal_abs_max": mu_e,
        "mu_beam_muon_line_shift_abs_max": c_max * p_mu,
        "mu_beam_inelastic_variance_max": c_max**2 * p_mu * (1 - p_mu),
        "tauB_over_mu_partial_ratio_LO": ratio,
        "tauB_upper": tau_upper,
        "direct_Belle_limit_over_upper": DATA["tau3e_Belle_90"] / tau_upper,
        "contact_scale_TeV": 1 / math.sqrt(2 * math.sqrt(2) * DATA["fermi_GeV_minus_2"] * c_max) / 1000,
    }


def witness(p_mu, p_e):
    # Both phases are unmeasured by the rates used here.
    line = np.sqrt([p_e, p_mu, 1 - p_mu - p_e]) * np.exp(1j * np.array([0, .8, -.4]))
    P = np.outer(line, line.conj())
    M = np.diag([0., 1., 0.])
    identity = np.eye(3)
    assert np.max(np.abs(P @ P - P)) < 1e-14
    assert np.isclose(np.linalg.norm(M @ P - P @ M, ord=2), math.sqrt(p_mu * (1 - p_mu)))
    charm = DATA["charm"]
    effective = (charm["gV"][0] - charm["gA"][0]) / 2
    c = .9 * math.sqrt(DATA["mu3e_ratio_90"] / (p_mu * p_e))
    variance = c**2 * p_mu * (1 - p_mu)
    # 'q' below denotes the combined common SM+new right coefficient.
    common_right = [sign * math.sqrt(effective**2 - variance) - c * p_mu
                    for sign in (1, -1)]
    for q in common_right:
        Q = q * identity + c * P
        assert np.isclose((Q @ Q)[1, 1].real, effective**2, atol=1e-14)
        assert np.isclose(abs(Q[0, 2] / Q[0, 1])**2, (1 - p_mu - p_e) / p_mu)
    return {
        "p_e": p_e, "p_mu": p_mu, "p_tau": 1 - p_mu - p_e, "c": c,
        "two_real_common_right_amplitudes": common_right,
        "common_only_witness": {"c": 0, "right": effective, "both_LFV_widths": 0},
        "incompatible_tags_commutator_norm": math.sqrt(p_mu * (1 - p_mu)),
        "limiting_mu_ratio": abs(c * P[0, 1])**2,
    }


def calculate():
    sectors = {}
    stress = DATA["stress"]
    f = DATA["mu3e_ratio_90"]
    for key in ("solar_no", "solar_io_prior"):
        mu, electron = inputs_for_order(DATA["residual"][key])
        nominal = bound(mu, electron, 1., 1., 1.)
        # Even arbitrary common charge drops out of the LFV amplitudes.
        c_trial = stress["minimum_scalar_tag_split_example"]
        sectors[key] = {
            "p_mu": mu, "p_e": electron, "nominal": nominal,
            "stress": bound(mu, electron, stress["minimum_mu_efficiency_ratio"],
                            stress["minimum_mu_rate_transfer"], stress["maximum_tau_rate_transfer"]),
            "mu_eff_times_rate_needed_for_1pct_eigenvalue_split":
                f / (c_trial**2 * mu * electron[0]),
        }
    charm = DATA["charm"]
    gv, ga = box(charm["gV"], 2), box(charm["gA"], 2)
    right = [(gv[0] - ga[1]) / 2, (gv[1] - ga[0]) / 2]
    mu, _ = inputs_for_order(DATA["residual"]["solar_no"])
    central_e = (1 + math.sqrt(1 - DATA["residual"]["fast"][0])) / 2 * (1 - DATA["residual"]["solar_no"][0])
    assert sectors["solar_no"]["nominal"]["tauB_upper"] < 1.39e-13
    return {
        "order": sectors, "charm_right_magnitude_box_for_published_branch": right,
        "standalone_c_bound_with_free_common_and_mapped_right_branch":
            right[1] / math.sqrt(mu * (1 - mu)),
        "witness": witness(mu, central_e),
    }


def stable(item):
    if isinstance(item, dict):
        return {k: stable(v) for k, v in item.items()}
    if isinstance(item, list):
        return [stable(v) for v in item]
    return float(f"{item:.11g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    opts = parser.parse_args()
    result = stable(calculate())
    out = ROOT / "results.json"
    if opts.check:
        if json.loads(out.read_text(encoding="utf-8")) != result:
            raise SystemExit("1617 expected results differ")
    else:
        with out.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print("1617 projector, three-electron transfer, inclusive right-electron response verified")
