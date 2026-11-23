"""Rare B calibration and independent charged-current bound for weak-singlet down quarks."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def absolute_test(row, kmin, kmax, c, C, weight, hfactor):
    n_low, n_high = row["abs_N"]
    p = row["abs_P_max"]
    # P bounds apply to the pseudoscalar amplitude even when a scalar is present.
    l_min = max(0., n_low-p)/(1-kmin)
    charged_min = max(0., l_min*(1/c-kmin)-p)
    l_max = c*(C+p)/(1-c*kmin)
    n_max = p+(1-kmin)*l_max
    free_f_min = max(0., (kmin*n_low-p)/(1-kmin))
    free_f_max = (kmax*n_high+p)/(1-kmax)
    f_from_mu = (p+kmax*c*C)/(1-kmax*c)
    e0_lever = (1-kmin)/(1-c*kmin)
    p_lever = 1+c*e0_lever
    p_needed = (n_low-e0_lever*c*C)/p_lever
    return {"N_measured_min":n_low, "P_measured_max":p,
            "min_abs_L_for_both":l_min, "charged_norm_required":charged_min,
            "max_N_predicted_from_mu_and_charged":n_max,
            "minimum_extra_N_amplitude":max(0.,n_low-n_max),
            "minimum_extra_L_if_only_L_changes":max(0.,l_min*(1-c*kmin)-c*(C+p)),
            "mu_transfer_floor_for_separation":weight*(p/p_needed)**2 if p_needed>0 else None,
            "F_bounds_if_free_L_and_both_rates":[free_f_min,free_f_max],
            "F_upper_mu_and_charged_only":f_from_mu,
            "tree_hbs_over_reference_hbb_upper":hfactor*f_from_mu**2}


def calculate():
    cfg = read("inputs.json")
    inputs = read(cfg["ratio_inputs"])
    coeff = inputs["sm_coefficient_ratio"]
    kmin,kmax = coeff["range"]
    tree = cfg["charged_current_caps"]
    C = math.hypot(tree["abs_Kcb"],tree["abs_Kub"])
    c = cfg["top_coefficient_cap"]
    absdata = read(cfg["absolute_calibration"])
    published = read(cfg["absolute_inputs"])
    ref = published["HPQCD_tableI_eq15"]
    c_ref = absdata["nu_reference_amplitude"]/ref["lambda_ref"]
    assert 0<kmin<kmax<1 and c_ref < c < 1
    w = {"CMS":inputs["CMS"]["accepted_light_to_heavy_weight_range"][0],
         "LHCb_R3":published["run3_working_light_response"]}
    massbox = cfg["higgs_tree_work_box"]
    mh = read(cfg["higgs_inputs"])["kinematic_GeV"]["H"]
    hfactor = (1+massbox["m_s_over_m_b_upper"]**2)/(1-4*(massbox["m_b_upper_GeV"]/mh)**2)**1.5
    absolute = {}
    for channel,weight in w.items():
        for widths in (1,2):
            label = f"{channel}_{widths}"
            absolute[label] = absolute_test(absdata["working_boxes"][label],kmin,kmax,c,C,weight,hfactor)
    rel = read(cfg["relative_calibration"])["comparisons"]["CMS"]
    disk = {}
    for widths,data in rel.items():
        up=data["accepted_mu_rate"][1]
        low=data["Belle_Rnu"][0]
        radius=math.sqrt(up/w["CMS"])
        max_amp=1+kmax*(radius-1)
        assert radius>1
        disk[widths]={"upper_Rnu_free_clock":max_amp**2, "observed_Rnu_floor":low,
                      "critical_time_weight":up/(1+(math.sqrt(low)-1)/kmax)**2}
    assert absolute["LHCb_R3_2"]["charged_norm_required"]>C
    assert absolute["CMS_2"]["minimum_extra_N_amplitude"]==0
    assert disk["2"]["upper_Rnu_free_clock"]>=disk["2"]["observed_Rnu_floor"]
    return {"charged_norm_cap":C,"ordinary_top_factor_reference":c_ref,
            "higgs_kinematic_ratio_upper":hfactor,"absolute":absolute,"relative_CMS_comparison":disk}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    val=calculate()
    if args.check:
        if val!=read("results.json"):
            raise SystemExit("1632: saved results differ")
        print("1632 PASS: axial sign, free normalization, charged closure and Higgs width bounds")
    else:
        with (HERE/"results.json").open("x",encoding="utf-8") as out:
            out.write(json.dumps(val,indent=2)+"\n")
