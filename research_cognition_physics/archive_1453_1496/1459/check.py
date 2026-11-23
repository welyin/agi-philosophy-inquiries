"""A lifetime and beta-endpoint test using the same electroweak triplet."""
import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def coefficient(mass, wmass, gf):
    x = (wmass/mass)**2
    if 4*x >= 1:
        raise ValueError("two on-shell W bosons required")
    return gf**2*mass**4*math.sqrt(1-4*x)*(1-4*x+12*x*x)/(32*math.pi**2)


def calculate():
    cfg = read("inputs.json")
    gaps = read(cfg["gaps"])["juno_207_days"]
    signed = read(cfg["oscillation"])["orders"]
    gf = read(cfg["fermi"])["charged_current_check"]["GF_GeVminus2"]
    solar_min = gaps["solar_21"][0]-2.05*gaps["solar_21"][1]
    rows = {}
    for order in ("NO","IO"):
        a = gaps["normal_atmospheric_31" if order=="NO" else "inverted_abs_31"]
        # Minimum sum of squared masses at zero lightest mass.
        trace_min = (a[0]-2*a[1])*(1 if order=="NO" else 2)+solar_min
        envelope = signed[order]["G_interval_eV2"]
        entries = []
        for m in cfg["mass_GeV"]:
            k = coefficient(m,cfg["w_mass_GeV"],gf)
            signed_products = [k*g*1e-18 for g in envelope]
            total_floor = 2*math.sqrt(k*trace_min*1e-18)
            entries.append({"mass_GeV":m, "signed_width_product_GeV2":signed_products,
                            "maximum_proper_length_mm":cfg["hbarc_GeV_mm"]/total_floor,
                            "minimum_total_width_GeV":total_floor,
                            "electron_weighted_product_ceiling_GeV2":k*cfg["katrin_mass_eV"]**2*1e-18})
        rows[order] = {"sum_mass_squared_floor_eV2":trace_min,"benchmarks":entries}

    # Inject independently chosen extra decay, compare total rate inferred from fractions to true rate.
    m = 500.
    vev = 7e-5
    mbeta2 = .00014
    g = .00235
    trace = g+3*mbeta2
    k = coefficient(m,cfg["w_mass_GeV"],gf)
    leptons = m*trace*1e-18/(16*math.pi*vev**2)
    ww = k*trace*1e-18/leptons
    total = (leptons+ww)/(1-cfg["unseen_branch_for_witness"])
    re = mbeta2/trace
    inferred = math.sqrt(k*g*1e-18/((ww/total)*(leptons/total)*(1-3*re)))
    assert math.isclose(total,inferred,rel_tol=1e-13)
    assert math.isclose(ww*leptons*re/(k*1e-18),mbeta2,rel_tol=1e-13)
    return {"orders":rows,"common_source_witness":{"M_GeV":m,"unseen_branch":cfg["unseen_branch_for_witness"],
            "total_width_GeV":total,"width_from_fractions_GeV":inferred,
            "independent_endpoint_mass_squared_eV2":mbeta2}}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    opt=p.parse_args()
    res=calculate()
    if opt.write_new:
        with (ROOT/"results.json").open("x",encoding="utf-8") as out:
            json.dump(res,out,indent=2)
            out.write("\n")
    else:
        assert res==read("results.json")
    for order, result in res["orders"].items():
        print(order,[(x["mass_GeV"],x["maximum_proper_length_mm"]) for x in result["benchmarks"]])
