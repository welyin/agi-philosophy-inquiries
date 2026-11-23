"""Single three-pass charged boundary and fixed-gain neutral continuation."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


CFG = read("inputs.json")


def module(key):
    return runpy.run_path(str(HERE / CFG[key]))


KAON = {k: tuple(r["central_ratio_and_marginal_error"]) for k, r in
        read(CFG["kaon_inputs"])["experimental_comparisons"].items()}
HAD = module("hadron")
REFERENCE = module("flavour")
DIRECT = REFERENCE["DIRECT"]
LEPTONS = REFERENCE["LEPTONS"]


def one_width(value, central, lower, upper=None):
    return (value-central) / ((lower if upper is None else upper) if value>central else lower)


def moduli(offsets):
    return {k: DIRECT[k][0]+offsets[k[1:]]*math.hypot(*DIRECT[k][1:]) for k in ("Vus","Vub","Vcb")}


def axial_ratio(offsets):
    v = moduli(offsets)
    return v["Vus"]/math.sqrt(1-v["Vus"]**2-v["Vub"]**2)


def strange_offset(width, vector):
    v = DIRECT["Vus"] if vector=="original" else CFG["alternative_vector_Vus"]
    val = v[0] + width*math.hypot(*v[1:])
    return (val-DIRECT["Vus"][0])/math.hypot(*DIRECT["Vus"][1:])


def axial_width(ratio_key, ub_offset=None, vector="original"):
    center, error = KAON[ratio_key]
    lo, hi = 0., 4.
    for _ in range(60):
        mid = (hi+lo)/2
        offsets = dict(us=strange_offset(mid, vector), ub=mid if ub_offset is None else ub_offset, cb=0.)
        if axial_ratio(offsets) < center-mid*error:
            lo = mid
        else:
            hi = mid
    return hi


def shoot(weak, up_shift, step):
    """Two simultaneous boundary residuals, with the three charged-lepton masses fixed."""
    flavour = module("flavour")
    up = HAD["UP"][0]+up_shift*HAD["UP"][1]
    flavour["FLOW"]["low"]["yukawa"]["up"][0][0] = up
    def run(point):
        row = flavour["case"]({**weak, "strange":point[0], "bottom":point[1]}, True, 63., step=step)
        return np.log([row["x"]/row["Vud_10PeV"]**CFG["k"], row["tau_MZ"]/LEPTONS[2,0]]), row
    point = np.array([-1.2,.8])
    r, row = run(point)
    jac = np.column_stack([(run(point+dx)[0]-r)/.05 for dx in .05*np.eye(2)])
    for _ in range(10):
        if max(abs(r)) < 2e-11:
            break
        change = np.linalg.solve(jac,-r)
        point += change
        new_r,row = run(point)
        jac += np.outer(new_r-r-jac@change,change)/(change@change)
        r = new_r
    assert max(abs(r)) < 2e-11
    dsb = np.array(row["MZ_down_yukawas"])
    dsb_width = (dsb-flavour["DOWN"][:,0])/flavour["DOWN"][:,1]
    neutral = module("neutral")
    original = neutral["load"]
    def boundary(path):
        charged = original(path)
        charged["FLOW"]["low"]["yukawa"]["up"][0][0] = up
        return charged
    neutral["endpoint"].__globals__["load"] = boundary
    end = neutral["endpoint"]("explicit_open",step,explicit_offsets={**weak,"strange":point[0],"bottom":point[1]})
    assert abs(end["row"]["gamma_MZ_deg"]-row["gamma_MZ_deg"]) < 1e-7
    assert abs(end["row"]["x"]-row["x"]) < 1e-8
    obs = neutral["observables"](end["C"],end["b"],neutral["FROZEN"][0])
    a = read(CFG["reactor_and_solar"])
    t = read(CFG["atmospheric_angle"])["atmospheric_NOvA_upper_octant"]
    a_boxes = [a["DayaBay_sin2_2theta13"],a["JUNO_sin2_theta12_NO"],[t["central"],t["low"],t["high"]]]
    gamma = read(CFG["direct_gamma"])["time_integrated_gamma_deg"]
    pulls = dict(**weak, u=up_shift, **dict(zip(("d","s","b"),dsb_width)),
                 **{k:one_width(axial_ratio(weak),*v) for k,v in KAON.items()},
                 gamma=one_width(row["gamma_MZ_deg"],gamma["central"],gamma["central"]-gamma["68"][0],gamma["68"][1]-gamma["central"]),
                 **{k:one_width(v,*widths) for k,v,widths in zip(("angle13","angle12","angle23"), obs["angle_readings"],a_boxes)},
                 dm31=one_width(obs["dm31_eV2"],*neutral["FREQ"]["juno_NO_dm31"]),
                 dm32=one_width(obs["dm32_eV2"],*neutral["FREQ"]["accelerator_NO_dm32"]))
    raw = HAD["observables"](up,dsb)
    out = dict(weak_offsets=weak, up_MZ=up, moduli_MZ=moduli(weak), axial_ratio=axial_ratio(weak),
               down_MZ=dsb, high_down=row["high_down_yukawas"], high_e_mu_tau=[*row["high_e_mu"],row["tau_10PeV"]],
               x=row["x"], c=row["Vud_10PeV"], tau_MZ=row["tau_MZ"], gamma=row["gamma_MZ_deg"],
               positive_soft_tan_beta=row["positive_zero_stop_tan_beta"], return_residual=row["return_relative"],
               boundary_log_residuals=r, observables=obs, solar_eV2=neutral["SOL"], signed_widths=pulls,
               largest_width=max(abs(v) for v in pulls.values()), hadrons_at_central_matching=HAD["hadrons"](raw))
    return out


def hadron_budgets(raw, at_width):
    out = {}
    for n, width_name in ((1,"one"),(2,"two"),(at_width,"at_witness")):
        def expand(x, interval):
            return (x+n*(interval[0]-x),x+n*(interval[1]-x))
        B = expand(HAD["HM"][0],HAD["BS"])
        Q = expand(HAD["ETA"][0],HAD["QB"])
        BMW = expand(HAD["LIN"][0],HAD["SB"])
        MILC = expand(HAD["CFG"]["MILC_S"][0],HAD["SS"])
        r = expand(HAD["CFG"]["MILC_light_ratio"][0],HAD["RR"])
        for label,Slight,R in (("BMW",BMW,None),("MILC",MILC,r),("both",(max(BMW[0],MILC[0]),min(BMW[1],MILC[1])),r)):
            test = dict(Sbox=Slight,rbox=R,Qbox=Q,bsbox=B)
            if Slight[0]>Slight[1]:
                out[label+"_"+width_name] = None
                continue
            budget = HAD["threshold"](raw,**test)
            witness = HAD["feasible"](raw,budget+1e-9,**test)
            ratios = np.array([witness["d_over_s"],witness["u_over_d"],witness["b_over_s"]])
            factors = raw/ratios/HAD["MATCH"]
            assert max(abs(factors-1)) <= budget+1.1e-9
            out[label+"_"+width_name] = dict(min_budget=budget,example=witness,
                                               K_over_K_reference=factors)
    return out


def kaon_projection(width, key, vector):
    """Projection onto Vus of the axial, direct strange and beauty boxes."""
    r, dr = KAON[key]
    bmax = moduli(dict(us=0.,ub=width,cb=0.))["Vub"]
    bmin = moduli(dict(us=0.,ub=-width,cb=0.))["Vub"]
    def inverted(a,b):
        return a*math.sqrt((1-b*b)/(1+a*a))
    s_input = DIRECT["Vus"] if vector=="original" else CFG["alternative_vector_Vus"]
    s,e = s_input[0], math.hypot(*s_input[1:])
    lo = max(s-width*e,inverted(r-width*dr,bmax))
    hi = min(s+width*e,inverted(r+width*dr,bmin))
    return [lo,hi] if lo<=hi else None


def main():
    rows = {}
    for key, (axis,vector,ub,cb,up) in CFG["candidates"].items():
        width = axial_width(axis,ub,vector)
        weak = dict(us=strange_offset(width,vector),ub=ub,cb=cb)
        row = shoot(weak,up,CFG["step"])
        touched = row["signed_widths"].copy()
        touched["us_original_reference"] = touched.pop("us")
        touched["us_alternative_reference"] = one_width(row["moduli_MZ"]["Vus"],CFG["alternative_vector_Vus"][0],
                                                         math.hypot(*CFG["alternative_vector_Vus"][1:]))
        used = touched.copy()
        if axis=="RBC_axial":
            used.pop("RM123_axial")
        if vector=="original":
            used.pop("us_alternative_reference")
        # The alternative also passes the original Vus requirement; both can be imposed.
        row["signed_widths"] = touched
        row["largest_width_in_named_box"] = max(abs(v) for v in used.values())
        row.pop("largest_width")
        assert row["largest_width_in_named_box"] < width+1e-8
        assert row["observables"]["SK_T2K_FC68_90"] == [False,True]
        assert abs(used["gamma"])<1
        raw = HAD["observables"](row["up_MZ"],row["down_MZ"])
        row["hadron_relative_matching"] = hadron_budgets(raw,width)
        assert row["hadron_relative_matching"]["both_at_witness"]["min_budget"] < CFG["matching_ceiling"]
        rows[key] = row
    core = rows["old_ub_cb_and_central_up"]
    assert all(core["observables"]["angles_in_box"]) and all(core["observables"]["dm31_and_32_in_box_at_central_solar"])
    # Different integrator resolution on both charged and neutral boundaries.
    half = shoot(core["weak_offsets"],0,CFG["step"]/2)
    refined = {"max_mass_width_shift":max(abs(half["signed_widths"][k]-core["signed_widths"][k]) for k in ("d","s","b")),
               "dm31_shift":abs(half["observables"]["dm31_eV2"]-core["observables"]["dm31_eV2"]),
               "gamma_shift_deg":abs(half["gamma"]-core["gamma"])}
    assert refined["max_mass_width_shift"]<1e-5 and refined["dm31_shift"]<2e-9 and refined["gamma_shift_deg"]<1e-5
    # A second solver brackets the bottom and strange masses successively.
    old = module("three_pass")
    old["WEAK"]["1736_original_calibrations"] = core["weak_offsets"]
    other = old["solve"](CFG["k"],"1736_original_calibrations",CFG["step"])
    nested = max(abs(np.array(other["down_strange_bottom_marginal_widths"]) -
                     np.array([core["signed_widths"][k] for k in ("d","s","b")])))
    assert nested<1e-6
    bounds = {axis+"_"+v:axial_width(axis,None,v) for axis in ("RM123_axial","RBC_axial") for v in ("original","alternative")}
    boxes = {axis+"_"+v:kaon_projection(2.,axis,v) for axis in ("RM123_axial","RBC_axial") for v in ("original","alternative")}
    assert boxes["RM123_axial_alternative"] is None and boxes["RM123_axial_original"] is not None
    return HAD["tidy"]({"lower_bounds_in_exact_unitary_matching":bounds, "two_width_Vus_projection":boxes,
                        "matching_reference":HAD["MATCH"],"witnesses":rows,
                        "half_step_check":refined,"nested_bracket_mass_width_difference":nested})


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args = cli.parse_args()
    result = main()
    path = HERE / "results.json"
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1736 verified joint kaon, charged flavour, neutrino and hadron window bounds")
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")
