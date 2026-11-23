"""Recompute the odd colored fermion portal, inclusive parton widths and rate thresholds."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def read(relative):
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def phase_space(r):
    if r == 0:
        return 1.0
    u = r * r
    return 1 - 8*u + 8*u**3 - u**4 - 12*u*u*math.log(u)


def quadrature(r):
    nodes, weights = np.polynomial.legendre.leggauss(64)
    t = (nodes + 1)/2
    xmin, xmax = (1-r)**2, (1+r)**2
    x = xmin*(1-t*t)
    integrand = x*(1+r*r-x)*np.sqrt((xmin-x)*(xmax-x))
    return 12*np.dot(weights/2, integrand*2*xmin*t)


def compute():
    inp = read("inputs.json")
    old = read(inp["frozen_texture"])
    frozen = read(inp["frozen_mass"])["conditional_masses_TeV_descending"][-1]*1000
    eps = old["scalar_magnitude_over_cutoff"]
    a = (eps*np.exp(1j*np.pi*old["scalar_angle_over_pi"]))**3
    matrix = np.array([[0, a], [-a, old["colored_M_T_over_cutoff"]]])
    _, norm, right_dagger = np.linalg.svd(matrix)
    mass = norm[-1]*old["cutoff_TeV"]*1000
    right_weight = abs(right_dagger[-1, 1])**2
    assert abs(mass-frozen) < 1e-7
    assert abs(right_weight-(1-4/math.sqrt(17))/2) < 1e-14
    # Charges are (Z12, aux Z3, superfield R, partner parity exponent, 6Y).
    q = (0, -1, 1, 0, 1)
    tbar = (0, 0, 1, 1, 2)
    j = (0, 1, 1, 1, 0)
    sums = [2*a-b-c for a, b, c in zip(q,tbar,j)]
    assert sums[0] % 12 == sums[1] % 3 == sums[2] == sums[3] % 2 == sums[4] == 0
    widths = []
    denominator = 6144*math.pi**3*inp["cutoff_GeV"]**4
    for r in inp["j_over_d"]:
        phase = phase_space(r)
        assert abs(phase-quadrature(r)) < 1e-13
        for coeff in inp["component_coefficients"]:
            gamma = inp["color_channels"]*right_weight*coeff**2*mass**5*phase/denominator
            length = inp["hbar_c_GeV_mm"]/gamma
            widths.append(dict(mass_ratio=r, coefficient=coeff, phase=phase,
                               width_GeV=gamma, proper_length_mm=length))
    dataset = inp["cms_pair_production_13TeV"]
    masses, rates = dataset["mass_GeV"], dataset["sigma_fb"]
    weight = (mass-masses[0])/(masses[1]-masses[0])
    sigma = rates[0]*(rates[1]/rates[0])**weight
    dv = inp["atlas_dv"]
    long_lived = inp["atlas_hscp"]
    dv_required = {key: value/sigma for key,value in dv["visible_95_fb"].items()}
    hscp_required = {key: value/(long_lived["lumi_inverse_fb"]*sigma)
                     for key,value in long_lived["long_lifetime_observed_s95"].items()}
    max_trans_boost = math.sqrt((inp["sqrt_s_GeV"]/(2*mass))**2-1)
    min_frac = min(dv_required.values())
    blind_length = dv["min_primary_transverse_separation_mm"]/(max_trans_boost*math.log(2/min_frac))
    half = next(row["proper_length_mm"] for row in widths
                if row["mass_ratio"]==0.5 and row["coefficient"]==1)
    return dict(light_mass_GeV=mass, light_Tbar_probability=right_weight,
                kahler_charge_sums=sums, width_examples=widths,
                mass_half_coefficient_for_1mm=math.sqrt(half),
                mass_half_coefficient_for_1m=math.sqrt(half/1000),
                interpolated_production_fb=sigma,
                DV_exclusion_required_event_Aeff=dv_required,
                long_lifetime_exclusion_required_event_Aeff=hscp_required,
                max_13TeV_transverse_boost=max_trans_boost,
                genuine_D_DV_rate_below_limit_for_ctau_mm_less_than=blind_length)


def round_result(obj):
    if isinstance(obj, dict):
        return {key: round_result(value) for key,value in obj.items()}
    if isinstance(obj, list):
        return [round_result(value) for value in obj]
    if isinstance(obj, (float, np.floating)):
        return float(f"{obj:.11g}")
    return obj


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = round_result(compute())
    destination = ROOT / "results.json"
    if args.check:
        assert read("results.json") == result
        print("1694: frozen mass, charges, phase space and event thresholds verified")
    else:
        destination.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
