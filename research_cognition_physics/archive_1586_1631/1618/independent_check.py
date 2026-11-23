"""Compare neutral-meson scalar decay with the inclusive Higgs width for the same rotations."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def calculate():
    cfg = read(HERE/"inputs.json")
    weak = read(HERE/cfg["weak"])
    yd = read(HERE/cfg["yukawa"])["yukawa"]["down"]
    yuk_up = read(HERE/cfg["yukawa"])["yukawa"]["up"]
    h = cfg["kinematics_GeV"]
    widths = cfg["marginal_widths"]
    higgs_lo = h["h_center"] - widths*h["h_sigma"]
    higgs_hi = h["h_center"] + widths*h["h_sigma"]
    v = cfg["v_GeV"]
    # CMS central and three upward marginal errors, in units of 1e-9.
    cms_upper = (3.83 + widths*(.38+.19+.14))*1e-9
    min_rate_response = .8
    y_s_max = .11
    f_bs_lo = .2307-widths*.0013
    m_bs_lo = 5.36
    mu = .105658
    tau_lo_natural = 1.45e-12 / 6.582119569e-25
    meson_factor = (tau_lo_natural*f_bs_lo**2*m_bs_lo**5*mu**2
                    *(1-4*mu**2/m_bs_lo**2)**1.5
                    /(32*math.pi*higgs_hi**4*v**4))
    minus_limit = math.sqrt((1+y_s_max)*cms_upper/(min_rate_response*meson_factor))

    # VBF WW supplies a total width limit under the quoted numerator/acceptance budget.
    vbf = cfg["ATLAS_VBF_WW_pb"]
    higgs_width = (cfg["h_reference_width_GeV"] * cfg["transfer"]["VBF_WW_numerator_and_acceptance_max"]
                   *(vbf["SM"]+widths*vbf["SM_error"])/(vbf["value"]-widths*vbf["lower_error"]))
    bs_phase = (1-(h["b_plus_s_upper"]/higgs_lo)**2)**2
    tree_fraction = cfg["transfer"]["tree_FCNC_width_min"]
    coupling_norm = math.sqrt(8*math.pi*higgs_width/(3*higgs_lo*tree_fraction*bs_phase))
    # As in the author test, the quoted running Yukawas are used with the declared rate budgets.
    mass_gap_over_v = (yd[2][0]-widths*yd[2][1]-yd[1][0]-widths*yd[1][1])/math.sqrt(2)
    plus_limit = math.sqrt(2)*coupling_norm/mass_gap_over_v
    pair_left = (minus_limit+plus_limit)/2
    def upper(key):
        row = weak[key]
        return row["central"]+widths*sum(row["errors"])
    down_velocity = pair_left*math.hypot(1,upper("td_over_ts"))

    # Independent evaluation of the up-sector bound using the inclusive top width and tree Wb.
    mt_min, mt_max = h["top"]
    w_min, w_max = h["W"]
    sb_max = h["b_plus_s_upper"]
    sin_delta = upper("cb")*math.hypot(1,upper("ub_over_cb"))
    wb_shape = ((1-(sb_max/mt_min)**2)**2 + (w_min/mt_max)**2 - 2*(w_max/mt_min)**4)
    wb_phase = math.sqrt((1-((w_max+sb_max)/mt_min)**2)*(1-(w_max/mt_min)**2))
    wb_floor = (cfg["GF_lower_GeV_inverse2"]*mt_min**3/(8*math.pi*math.sqrt(2))
                *(1-sin_delta**2)*wb_shape*wb_phase*cfg["transfer"]["tree_Wb_width_min"])
    yt = yuk_up[2][0]-widths*yuk_up[2][1]
    yc = yuk_up[1][0]+widths*yuk_up[1][1]
    ratio_h, ratio_c = higgs_hi/mt_min, h["charm_upper"]/mt_min
    kinematic = math.sqrt((1-(ratio_h+ratio_c)**2)*(1-ratio_h**2))*(1-ratio_h**2-2*ratio_c)
    min_top_coefficient = (mt_min/(32*math.pi)*kinematic*(yt**2-yc**2)**2
                           /(2*(yt**2+yc**2)*(1+upper("ub_over_cb")**2)))
    up_velocity = math.sqrt((cfg["CMS_top_same_sign"]["top_width_upper_GeV"]-wb_floor)
                            /(tree_fraction*min_top_coefficient))
    required = 2*(math.pi-math.asin(sin_delta))
    for left,right in [(2+1j,-.7+.3j),(1j,1j)]:
        b,s = 3.,.12
        ys = b*left-s*right
        yb_conj = -s*left+b*right
        assert abs((ys-yb_conj)/(b+s)-(left-right)) < 1e-12
        assert abs((ys+yb_conj)/(b-s)-(left+right)) < 1e-12
    assert up_velocity+down_velocity < required
    return {"CMS_branching_upper_working_box":cms_upper,
            "scalar_branching_coefficient_lower":meson_factor,
            "minimal_Bs_rate_response_assumed":min_rate_response,
            "max_chiral_difference_from_Bs":minus_limit,"max_chiral_sum_from_Higgs":plus_limit,
            "max_down_velocity":down_velocity,"max_up_velocity_without_FCNC_search":up_velocity,
            "combined_upper":up_velocity+down_velocity,"same_sign_lower":required}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    result = calculate()
    path = HERE/"independent_results.json"
    if parser.parse_args().check:
        old = read(path)
        assert old.keys()==result.keys()
        assert all(math.isclose(old[k],v,rel_tol=1e-12,abs_tol=1e-20) for k,v in result.items())
        print("1618 PASS: scalar Bs rate, chiral blind direction, VBF and top inclusive widths")
    else:
        with path.open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(result,indent=2)+"\n")
