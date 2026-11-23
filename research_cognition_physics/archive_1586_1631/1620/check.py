"""Distinguish a zero up source from a partly induced mass using two decay extractions."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def band(mean, error, n):
    return [mean-n*error, mean+n*error]


def q_of(s_ratio, up_over_down):
    return 0.5 * math.sqrt((s_ratio*s_ratio-1)*(1+up_over_down)/(1-up_over_down))


def r_of(s_ratio, q):
    k = (s_ratio*s_ratio-1)/(4*q*q)
    return (1-k)/(1+k)


def slope_test(s_high, ms_high, b_high, q_low):
    largest_r = ms_high*b_high
    max_q = q_of(s_high, largest_r)
    needed = r_of(s_high, q_low)/ms_high
    return {"Q_upper_zero_direct_source": max_q,
            "min_induced_coefficient_for_data_GeV_inv": needed,
            "min_opposite_sea_slope_GeV_inv": needed-b_high,
            "strong_isospin_signal_lower_relative_to_reported": (q_low/max_q)**4,
            "eta_amplitude_scale_upper_if_rescaled": (max_q/q_low)**2}


def cofactor_point(microscopic, strength):
    u, d, s = microscopic
    return microscopic + strength*np.array([d*s, u*s, u*d])


def calculate():
    cfg = read(HERE/"inputs.json")
    n = cfg["error_half_widths"]
    geo = read(HERE/cfg["geometry_inputs"])
    ss, stat, sys = geo["ms_over_average"]
    s_band = band(ss, math.hypot(stat, sys), n)
    ms, ms_err = read(HERE/cfg["strange"])["lattice"]["ms"]
    ms_band = band(ms/1000, ms_err/1000, n)
    b, e1, e2 = read(HERE/cfg["slope"])["sea_slope_GeV_inv"]
    b_band = band(b, e1+e2, n)
    channels = {"eta": geo["eta_q"], "eta_prime": cfg["eta_prime_Q"]}
    zero = {}
    for name, (central, error) in channels.items():
        q_low = central - n*error
        zero[name] = slope_test(s_band[1], ms_band[1], b_band[1], q_low)
        assert zero[name]["Q_upper_zero_direct_source"] < q_low
        assert abs(q_of(s_band[1], r_of(s_band[1], q_low))-q_low) < 1e-12

    # Freeze the earlier HKS fit. No direct up Yukawa nor either eta decay is used here.
    result = read(HERE/cfg["geometry"])["per_scale"]["MZ"]["0"]["ordered_positive_root_ranges"]
    y_up = sum(result[0])/2
    y = read(HERE/cfg["yukawa"])["yukawa"]
    yd, ys = [pair[0] for pair in y["down"][:2]]
    m0 = np.array([y_up, yd, ys]) * (ms/1000/ys)
    u0, d0, s0 = m0
    raw_s = 2*s0/(u0+d0)
    kappa = (raw_s-ss)/(ss*s0-raw_s*u0*d0/s0)
    u_eff, d_eff, s_eff = cofactor_point(m0, kappa)
    actual_s = 2*s_eff/(u_eff+d_eff)
    actual_r = u_eff/d_eff
    actual_q = q_of(actual_s, actual_r)
    assert math.isclose(actual_s, ss, rel_tol=1e-13)
    assert b_band[0] <= kappa <= b_band[1]
    for central, error in channels.values():
        assert abs(actual_q-central) < n*error
    assert np.max(m0) == s0 and u0 > 0
    assert abs((u0/d0+kappa*s0)/(1+kappa*s0*u0/d0)-actual_r) < 1e-13
    kappa_interval_from_S = [
        (raw_s-v)/(v*s0-raw_s*u0*d0/s0) for v in reversed(s_band)]
    # Independent comparison with the old short-distance up holdout and the isoscalar response.
    old = read(HERE/cfg["geometry"])["per_scale"]["MZ"]["2"]
    assert old["held_out_yu"][0] > old["ordered_positive_root_ranges"][0][1]
    assert raw_s > s_band[1]
    return {"isoscalar_band": s_band, "strange_GeV_band": ms_band, "sea_slope_GeV_inv_band": b_band,
            "zero_short_distance": zero,
            "finite_direct_HKS_example": {
                "short_masses_GeV": m0.tolist(), "induced_kappa_GeV_inv": kappa,
                "at_fixed_seed_kappa_interval_from_S_only": kappa_interval_from_S,
                "raw_S": raw_s, "raw_Q": q_of(raw_s, u0/d0), "effective_S": actual_s,
                "effective_r": actual_r, "effective_Q": actual_q,
                "fraction_of_effective_up_from_cofactor": (u_eff-u0)/u_eff,
                "effective_masses_GeV": [u_eff, d_eff, s_eff],
                "needed_interpretation_shift_in_S_at_central_inputs": raw_s/actual_s-1,
                "short_up_yukawa_high": old["ordered_positive_root_ranges"][0][1],
                "short_up_heldout_low": old["held_out_yu"][0] },
            "source_roles": ["eta_width uses absolute dispersive normalization",
                "eta_prime uses charged-neutral Dalitz shapes with separate floated normalizations",
                "sea slope is total, opposite-signed non-topological terms can hide a larger induced term"]}


def same(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(same(i, j) for i, j in zip(a,b))
    if isinstance(a, str):
        return a == b
    return math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-13)


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    target = HERE/"results.json"
    if args.write:
        target.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    else:
        assert same(result, read(target))
    print("1620: partial and zero-source spectra, widths and strange-sea slope verified")


if __name__ == "__main__":
    main()
