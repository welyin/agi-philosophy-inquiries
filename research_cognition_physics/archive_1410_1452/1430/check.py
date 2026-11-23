"""Certified scalar signal envelopes in the NA62 phase-space normalization."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def integrals(k, pion, mu, electron, cells):
    splus, sminus = (mu + electron)**2, (mu - electron)**2
    smax = (k - pion)**2
    grid = np.linspace(splus, smax, cells + 1)
    steps = np.diff(grid)
    recoil = np.sqrt(np.maximum(0., ((k + pion)**2 - grid) * (smax - grid)))
    opening = np.sqrt(np.maximum(0., (1. - splus/grid) * (1. - sminus/grid)))
    low = recoil[1:] * opening[:-1]
    high = recoil[:-1] * opening[1:]
    # Monotone factors give cell bounds. Padding is much larger than summation roundoff.
    integrals = []
    for power in (0, 1):
        lower = np.dot(steps, low * (grid[:-1] - sminus)**power)
        upper = np.dot(steps, high * (grid[1:] - sminus)**power)
        integrals.append([float(lower)*(1-1e-10), float(upper)*(1+1e-10)])
    return integrals


def evaluate():
    inp = read(HERE / "inputs.json")
    calib = read(HERE / inp["angles_masses"])
    atlas = read(HERE / inp["atlas_et"])
    k, pion = inp["MK_GeV"], inp["Mpi_GeV"]
    mu, electron = inp["mu_pole_GeV"], inp["e_pole_GeV"]
    h, v = inp["mh_GeV"], inp["v_GeV"]
    smax, sminus = (k - pion)**2, (mu - electron)**2
    number = atlas["width_multipliers"]
    lower = lambda item: item[0] - number*item[1]
    upper = lambda item: item[0] + number*item[-1]
    amp = lower(calib["DayaBay_amplitude"])
    x = (1 - math.sqrt(1-amp))/2
    y = upper(calib["NO"]["atmospheric_sine_squared"])
    q = upper(calib["quark_sines_and_phase"][0])
    dy, sy = calib["yukawa"]["down"][:2]
    ey, my, ty = calib["yukawa"]["charged_leptons"]
    sd_ratio = (lower(sy)+upper(dy))/(lower(sy)-upper(dy))
    mass_factor = (upper(ey)**2+upper(my)**2)/(lower(ey)**2+lower(ty)**2)**2 * 2/v**2
    flavor = q*q*(1-q*q) * y/(x*(1-x)*(1-y)**2)
    transfer = 1 - atlas["leptonic_relative_rate_factor_half_width"]
    atlas_rate = atlas["atlas"]["MC_marginal95"]["et"]
    numerator = (8*math.pi*atlas_rate/(h*transfer))**2
    propagator = (k*k-pion*pion)*inp["shape_abs_upper"]*sd_ratio/(2*(h*h-smax))
    phase_norm = inp["tauK_s"] / (inp["hbar_GeV_s"] * 256*math.pi**3*k**3)
    coefficient = numerator * mass_factor * flavor * propagator**2 * phase_norm
    bounds = integrals(k,pion,mu,electron,inp["cells"])
    physical = coefficient*bounds[1][1]
    flat = coefficient*bounds[0][1]*(smax-sminus)
    na62 = inp["NA62"]
    scenarios = {}
    for label, width in inp["hwidth_GeV"].items():
        equivalent = flat * width**2
        scenarios[label] = {"physical_B_upper":physical*width**2, "accepted_flat_B_upper":equivalent,
                            "mean_events_upper_nominal_SES":equivalent/na62["ses"],
                            "mean_events_upper_SES_minus_one_width":equivalent/(na62["ses"]-na62["ses_error"]),
                            "observed_limit_divided_by_count_envelope":na62["flat_upper90"]/equivalent}
    # A fully specified real phase, central masses and a piecewise linear shape.
    # Only the physical rate is calculated here. The user's detector response is unknown.
    center_x = (1 - math.sqrt(1-calib["DayaBay_amplitude"][0]))/2
    center_y = calib["NO"]["atmospheric_sine_squared"][0]
    center_q = calib["quark_sines_and_phase"][0][0]
    me, mm, mt = [v*i[0]/math.sqrt(2) for i in (ey,my,ty)]
    sd_center = (sy[0]+dy[0])/(sy[0]-dy[0])
    grid = np.linspace((mu+electron)**2,smax,inp["cells"]+1)
    mid = (grid[1:]+grid[:-1])/2
    shape = np.interp(mid,*np.transpose(inp["ETMC_table_IV_s_f0"]))
    omega = np.sqrt(((k+pion)**2-mid)*(smax-mid)*(1-(mu+electron)**2/mid)*(1-sminus/mid))
    real_trace = (mm**2+me**2)*(mid-mu**2-electron**2)-4*mm*me*mu*electron
    assert np.all(real_trace > 0)
    projection = center_q**2*(1-center_q**2)*center_y/(center_x*(1-center_x)*(1-center_y)**2)
    propagator_sq = ((k*k-pion*pion)*sd_center*shape/(2*(h*h-mid)))**2
    central = (8*math.pi*atlas_rate/h)**2*phase_norm*projection/(mt**2+me**2)**2
    example = float(central*np.dot(np.diff(grid), omega*real_trace*propagator_sq))
    assert 0 < physical < flat and example < physical
    assert scenarios["matched_offshell"]["mean_events_upper_SES_minus_one_width"] < .05
    assert scenarios["direct"]["accepted_flat_B_upper"] > na62["flat_upper90"]
    assert smax < inp["ETMC_table_IV_s_f0"][-1][0]
    return {"smax_GeV2":smax,"box_corner":{"x_min":x,"y_max":y,"q_max":q,"sd_scalar_ratio_max":sd_ratio},
            "monotone_integral_bounds_I0_I1":bounds,"physical_B_upper_per_GeV2_width":physical,
            "accepted_flat_B_upper_per_GeV2_width":flat,"scenarios":scenarios,
            "width_GeV_guaranteeing_below_flat_upper":math.sqrt(na62["flat_upper90"]/flat),
            "illustration_central_inputs_real_source_saturating_et":{label:example*width**2 for label,width in inp["hwidth_GeV"].items()},
            "swapped_row_h_only_B":0.0}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = evaluate()
    target = HERE / "results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    elif args.check:
        def compare(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(compare(a[i],b[i]) for i in a)
            if isinstance(a,list):
                return len(a)==len(b) and all(compare(x,y) for x,y in zip(a,b))
            return math.isclose(a,b,rel_tol=2e-11,abs_tol=1e-28)
        assert compare(result,read(target)), "Stored results differ."
        print("1430 verified")
    else:
        print(json.dumps(result,indent=2))
