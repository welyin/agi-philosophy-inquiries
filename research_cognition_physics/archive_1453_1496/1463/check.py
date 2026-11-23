"""Cross-family tree rank test using scalar chirality and independent t and B data."""
import argparse
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def allowed_f(rate_over_norm, scalar_weight, pmax):
    # Minimize over an arbitrary complex axial background inside the specified disk.
    scalar_only = math.sqrt(rate_over_norm)/scalar_weight
    if pmax == "free" or scalar_only <= pmax:
        return scalar_only
    w2 = scalar_weight**2
    return (pmax + math.sqrt((1+w2)*rate_over_norm-w2*pmax**2))/(1+w2)


def evaluate():
    cfg=read("inputs.json")
    old=read(cfg["bottom_benchmark"])
    top=read(cfg["top_benchmark"])["benchmark"]
    phys=old["const_GeV"]
    mt,mh,tot=top["m_top_GeV"],top["m_h_GeV"],top["top_width_GeV"]
    v=1/math.sqrt(math.sqrt(2)*phys["GF"])
    mb,mm=phys["M_Bd"],phys["m_mu"]
    beta=math.sqrt(1-4*(mm/mb)**2)
    lifetime=old["tau_ps"]["d"][0]-2*old["tau_ps"]["d"][1]
    norm=lifetime/cfg["hbar_GeV_ps"]*mb*beta/(8*math.pi)*(cfg["fBd_floor_GeV"]*mm/(2*v**2))**2
    scalar_weight=beta*mb**2/mh**2*cfg["mass_asymmetry_floor"]
    kt=mt**3/(32*math.pi*v**2)*(1-(mh/mt)**2)**2
    up=math.sqrt(cfg["cms_hu_branch_limit"]*tot/kt)
    row=old["LHCb_run123"]["d_in_1e10"]
    rates={"cms_published_95":old["CMS_d_95_CL_upper_in_1e10"]*1e-10,
           "lhcb_2_marginal_widths":(row[0]+2*(row[2]+row[4]))*1e-10}
    results={}
    for key,rate in rates.items():
        bounds={}
        for pmax in cfg["other_pseudoscalar_limits"]:
            f=allowed_f(rate/norm,scalar_weight,pmax)
            axial_floor=0 if pmax=="free" else max(0,f-pmax)
            assert math.isclose(norm*(axial_floor**2+(scalar_weight*f)**2),rate,rel_tol=1e-10)
            bounds[str(pmax)]={"F_db_max":f,"max_abs_Rud_times_Rtb":up*f}
        results[key]=bounds
    ratios=read("../1386/inputs.json")["yukawa"]["down"]
    q=(ratios[0][0]+2*ratios[0][1])/(ratios[2][0]-2*ratios[2][1])
    assert (1-q)/(1+q)>cfg["mass_asymmetry_floor"]
    return {"tree_v_GeV":v,"width_t_to_hu_per_unit_F_GeV":kt,
            "ATLAS_single_analysis_times_CMS_Bd_free_p":math.sqrt(cfg["atlas_hu_multilepton_only"]/cfg["cms_hu_branch_limit"])*results["cms_published_95"]["free"]["max_abs_Rud_times_Rtb"],
            "top_Fut_upper_at_nominal_response":up,
            "Bd_branch_norm_floor":norm,"scalar_weight_floor":scalar_weight,
            "max_product_rescales_with_top_rate_budget":[math.sqrt(a) for a in cfg["higgs_rate_transfer"]],
            "lower_factor_for_Bd_branch_times_top_partialwidth_GeV":kt*norm*scalar_weight**2,
            "bounds":results}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    q=p.add_mutually_exclusive_group(required=True)
    q.add_argument("--check",action="store_true")
    q.add_argument("--write-new",action="store_true")
    opt=p.parse_args()
    res=evaluate()
    if opt.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as f:
            json.dump(res,f,indent=2)
            f.write("\n")
    else:
        assert res==read("results.json")
    print(json.dumps(res["bounds"]))
