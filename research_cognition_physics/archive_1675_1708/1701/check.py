"""Space-like electromagnetic response and pair cuts from the previously matched spectra."""
from argparse import ArgumentParser
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


CFG = read("inputs.json")
ALPHA = 1/CFG["inverse_alpha0"]
X, WEIGHTS = np.polynomial.legendre.leggauss(64)
X = (1 + X)/2
WEIGHTS /= 2


def response(mass, endpoints, scalar=False):
    a, b = endpoints
    u = X*(1-X)
    # This ratio cancels the subtraction constant without losing small differences.
    log_ratio = np.log1p((b-a)*u/(mass*mass+a*u))
    integrand = (1-2*X)**2/8 if scalar else u
    return 2*ALPHA/math.pi * np.dot(WEIGHTS, integrand*log_ratio)


def gauge_log_mass(m, soft):
    return 2*math.log(m)/3 + math.log(math.hypot(m, soft))/3


def bisect_partner(low, equal, soft):
    total = 2*gauge_log_mass(equal, soft)
    hi = equal
    while gauge_log_mass(low,soft) + gauge_log_mass(hi,soft) < total:
        hi *= 2
    lo = equal
    for _ in range(90):
        middle = math.sqrt(lo*hi)
        if gauge_log_mass(low,soft) + gauge_log_mass(middle,soft) < total:
            lo = middle
        else:
            hi = middle
    return (lo+hi)/2


def run():
    frozen = read(CFG["charged_gauge_higgs"])["shared_witness"]
    zero = read(CFG["other_zero_boundaries"])["actual_zero_boundary_solutions"]
    weak_rows = [("30TeV", frozen["soft_GeV"], frozen["weak_fermions_GeV"][0])]
    weak_rows += [(name, p["soft_GeV"], p["gauge"]["equal_weak_fermions_GeV"])
                  for name, p in zero.items()]
    masses_color = frozen["colored_fermions_GeV"]
    floor = read(CFG["mz_reference"])["mu_ref_GeV"]
    output = {}
    for name, ms, m in weak_rows:
        windows = {}
        for experiment in ("opal", "l3_mid", "l3"):
            interval = CFG[experiment+"_transfer_GeV2"]
            assert floor**2 > math.sqrt(math.prod(interval))/4
            weak = 2*response(m, interval)
            others = 4*response(math.hypot(m,ms), interval, scalar=True)
            ordinary = 2*response(ms,interval)+17*response(ms,interval,scalar=True)
            color = sum((response(d,interval)+2*response(math.hypot(d,ms),interval,scalar=True))/3
                        for d in masses_color)
            second = bisect_partner(floor, m, ms)
            split = response(floor,interval)+response(second,interval)
            assert split >= weak
            # Convexity in the conserved log threshold bounds the fermion endpoints.
            # For either scalar eigenmass >= ms, its positive response is at most the value at ms.
            min_sum = weak + color + ordinary
            max_sum = split + 4*response(ms,interval,scalar=True) + color + ordinary
            windows[experiment] = dict(weak_fermions=weak, weak_scalars=others, color=color,
                                       ordinary_superpartners=ordinary, full_benchmark=weak+others+color+ordinary,
                                       full_interval_for_every_gauge_compatible_split=[min_sum,max_sum],
                                       split_at_minimum_mass_GeV=[floor,second], split_weak_fermions=split)
        threshold = 2*m
        velocity = math.sqrt(max(0, 1-(threshold/CFG["timelike_energy_GeV"])**2))
        one_cut = velocity*(3-velocity*velocity)/2
        output[name] = dict(equal_weak_mass=m, pair_cut_GeV=threshold, alpha_windows=windows,
                            equal_mass_300GeV_R_EM=2*one_cut,
                            any_split_300GeV_R_EM_lower_bound=one_cut)
    opal_error = math.hypot(*CFG["opal_difference"][1:])
    weak = [row["alpha_windows"]["opal"]["full_benchmark"] for row in output.values()]
    high = [row["alpha_windows"]["l3"]["full_benchmark"] for row in output.values()]
    mid = [row["alpha_windows"]["l3_mid"]["full_benchmark"]/ALPHA for row in output.values()]
    bounds = [row["alpha_windows"]["l3"]["full_interval_for_every_gauge_compatible_split"]
              for row in output.values()]
    certified_gap = bounds[0][0] - max(pair[1] for pair in bounds[1:])
    assert certified_gap > 0.00002
    assert high[0] > high[1] > high[2]
    # Low-energy expansion of the Dirac and complex scalar responses.
    for scalar,denom in [(False,15),(True,120)]:
        derivative = response(1000,(0,.0001),scalar)
        assert abs(derivative/(ALPHA/math.pi*.0001/(1000**2*denom))-1)<1e-8
    return dict(witnesses=output, OPAL_quadrature_sigma=opal_error,
                largest_OPAL_pair_separation_over_sigma=(max(weak)-min(weak))/opal_error,
                L3_mid_inverse_alpha_benchmark_shifts=mid,
                largest_L3_mid_pair_separation_over_sigma=(max(mid)-min(mid))/CFG["l3_mid_inverse_alpha_difference"][1],
                L3_finite_transfer_30TeV_minimum_gap_from_other_two_full_families=certified_gap,
                L3_two_heavier_benchmark_separation=high[1]-high[2],
                L3_two_lighter_benchmark_separation=high[0]-high[1])


if __name__ == "__main__":
    cli=ArgumentParser()
    cli.add_argument("--check",action="store_true")
    cfg=cli.parse_args()
    def stable(x):
        if isinstance(x,dict):
            return {key:stable(value) for key,value in x.items()}
        if isinstance(x,list):
            return [stable(value) for value in x]
        return float(f"{x:.10g}")
    result=stable(run())
    path=HERE/"results.json"
    if cfg.check:
        assert result == json.loads(path.read_bytes())
        print("PASS: three finite transfer windows, weighted-mass splitting floor and three pair thresholds")
    else:
        with path.open("x",encoding="utf8") as fp:
            fp.write(json.dumps(result,indent=2)+"\n")
