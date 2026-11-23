"""Round 1609: integrate the Dirac-triplet Born rate, then test event-level acceptances."""
import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


@lru_cache(maxsize=None)
def quadrature(n):
    abscissas, weights = np.polynomial.legendre.leggauss(n)
    return (1 + abscissas) / 2, weights / 2


def pdf_plane(inputs, scale):
    by_scale = inputs["xf_by_Q_GeV"]
    grid_q = np.array([float(q) for q in by_scale])
    index = np.searchsorted(grid_q, scale)
    assert 0 < index < len(grid_q)
    fraction = math.log(scale / grid_q[index - 1]) / math.log(grid_q[index] / grid_q[index - 1])
    planes = list(by_scale.values())
    xf = (1 - fraction) * np.array(planes[index - 1]) + fraction * np.array(planes[index])
    alpha = np.interp(np.log(scale), np.log(inputs["alpha_Q_GeV"]), inputs["alpha_values"])
    return np.log(inputs["x_nodes"]), xf, alpha


def born(inputs, mass, scale_factor=1.0, n=128):
    log_nodes, xf, alpha = pdf_plane(inputs, mass * scale_factor)
    coord, weight = quadrature(n)
    length = 2 * math.log(inputs["collider_GeV"] / (2 * mass))
    # t=log(shat/(4M²)), the second integration variable is -log(x1).
    t = length * coord**2
    remaining = (length - t)[:, None]
    log_x1 = -remaining * coord[None, :]
    log_x2 = -remaining - log_x1
    assert min(log_x1.min(), log_x2.min()) >= log_nodes[0]
    shat = 4 * mass**2 * np.exp(t)
    r = mass**2 / shat
    beta = np.sqrt(1 - 4 * r)
    log_ratio = 2 * np.arctanh(beta)
    common = math.pi * alpha**2 / shat
    sigma_qq = (8 / 27) * common * beta * (1 + 2 * r)
    sigma_gg = common / 3 * ((1 + 4*r + r*r)*log_ratio - (7 + 31*r)*beta/4)
    assert np.all(sigma_qq > 0) and np.all(sigma_gg > 0)

    def density(flavour, log_x):
        col = inputs["pdg_ids"].index(flavour)
        return np.interp(log_x, log_nodes, xf[:, col])

    luminosity_gg = density(21, log_x1) * density(21, log_x2)
    luminosity_qq = np.zeros_like(log_x1)
    for flavour in range(1, 6):
        for sign in [-1, 1]:
            luminosity_qq += density(sign * flavour, log_x1) * density(-sign * flavour, log_x2)
    measure = ((2 * length * coord * weight)[:, None] * remaining
               * weight[None, :] * inputs["gev_inverse_squared_to_fb"])
    gg = float(np.sum(measure * luminosity_gg * sigma_gg[:, None]))
    qq = float(np.sum(measure * luminosity_qq * sigma_qq[:, None]))
    return {"gg_fb":gg, "qqbar_fb":qq, "total_fb":gg + qq, "alpha_s":float(alpha)}


def solve(inputs):
    lumi = inputs["luminosity_fb_inverse"]
    limit = inputs["table9"]["S95obs"] / lumi
    cms = inputs["cms_table1_1100_fb"]
    central = cms["central"]
    lum_low = lumi * (1 - inputs["quoted_luminosity_fraction"])
    eta = limit / central
    output = {
        "visible_limit_from_displayed_S95_fb": limit,
        "zero_count_background_free_bound_events": -math.log(0.05),
        "zero_count_visible_bound_fb": -math.log(0.05) / lumi,
        "light_pair_yield_before_selection": central * lumi,
        "light_efficiency_using_central": eta,
        "light_efficiency_using_quote_lower_edge": limit / (central - cms["error"]),
        "light_efficiency_lower_cross_section_and_luminosity": 3 / (lum_low * (central - cms["error"])),
        "response_examples_only": {
            "expected_if_double_track_efficiency_0.0005": central*lumi*.0005,
            "expected_if_double_track_efficiency_0.002": central*lumi*.002,
            "opposite_fixed_sign_independent_legs_min_p_if_other_eff_1": math.sqrt(eta),
            "same_min_p_if_other_eff_0.10": math.sqrt(eta/.10),
        },
        "born_grid": {}, "double_track_efficiency_2100": {},
    }
    for mass in (1100, 2000, 2100):
        grid = {}
        for factor in (.5, 1., 2.):
            grid[str(factor)] = born(inputs, mass, factor, n=160)
        coarse = born(inputs, mass, n=80)["total_fb"]
        central_rate = grid["1.0"]["total_fb"]
        grid["relative_80_to_160_central"] = abs(coarse / central_rate - 1)
        assert grid["relative_80_to_160_central"] < .002
        output["born_grid"][str(mass)] = grid
    heavy = output["born_grid"]["2100"]
    for scale in ("0.5", "1.0", "2.0"):
        output["double_track_efficiency_2100"][scale] = limit / heavy[scale]["total_fb"]
    output["minimal_sigma2_fb_at_given_efficiency"] = {str(e): limit/e for e in (.05, .1, .25, .5, 1.)}
    output["born_1100_to_CMS_central_ratio"] = output["born_grid"]["1100"]["1.0"]["total_fb"]/central
    output["born_2000_to_Fuks_Shao_LO_ratio"] = output["born_grid"]["2000"]["1.0"]["total_fb"]/.08980
    output["positive_two_species_working_example_events"] = lumi*(central*.0009 + heavy["1.0"]["total_fb"]*.10)
    assert 0 < eta < .002
    return output


def compare(current, stored):
    if isinstance(current, dict):
        assert current.keys() == stored.keys()
        for key, value in current.items():
            compare(value, stored[key])
    else:
        assert math.isclose(current, stored, rel_tol=2e-9, abs_tol=1e-13), (current, stored)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = solve(json.loads((HERE / "inputs.json").read_text(encoding="utf-8-sig")))
    target = HERE / "results.json"
    if args.check:
        compare(result, json.loads(target.read_text(encoding="utf-8-sig")))
        print("1609 OK: CT14lo, three masses, scale variation and independent ATLAS event counts")
    else:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
