"""Real K rare rates, their matched complex amplitude, and a single-doublet right row."""
import argparse
import json
from math import pi, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def solve():
    p = read("inputs.json")
    rare = read(p["historical_rare"])["B_rare"]
    charged = read(p["Fuu_and_charged_product"])
    fu = charged["weak_charges_Fu_upper"]
    target = charged["charged_windows"][1]["abs_RdRs_lower"]
    lam = p["lambda_reference"]
    h = pi*p["sW2"]*p["inverse_alpha_MZ"]/lam**5
    bcentral = rare["central"]*1e11
    bmin, bmax = (rare["central"]-2*rare["minus"])*1e11, (rare["central"]+2*rare["plus"])*1e11
    kplus, kp_err = p["kplus_e11"]
    klong, kl_err = p["kL_e11"]
    electrom = 1 + p["delta_EM"]
    kp0, kpl, kph = electrom*kplus, electrom*(kplus-2*kp_err), electrom*(kplus+2*kp_err)
    rmin, rmax = sqrt(bmin/kph), sqrt(bmax/kpl)
    kl_max = (klong+2*kl_err)*rmax**2
    assert kl_max < p["KOTO_90CL_e11"]
    sm, sme = p["charged_BV_e11"]
    neutral, _ = p["long_BV_e11"]
    yi = sqrt(neutral/klong)
    smamp = complex(-sqrt(sm/kp0-yi**2),yi)
    assert abs(kp0*abs(smamp)**2-sm)<1e-12
    mirror = []
    for sx, sy in [(1,1),(-1,1),(1,-1),(-1,-1)]:
        delta = (smamp - complex(sx*smamp.real,sy*smamp.imag))/h
        tot = smamp - h*delta
        assert abs(kp0*abs(tot)**2-sm) < 1.e-12 and abs(klong*tot.imag**2-neutral)<1.e-12
        mirror.append([delta.real,delta.imag])
    sm_bound = sqrt((sm+2*sme)/kpl)
    s12,s23,s13 = p["left_angle_conditional_upper_s12_s23_s13"]
    tbound = (s12*s23+s13)*(s23+s12*s13)
    cbound = s12+s23*s13
    xtop = p["Xtop"][0]+2*p["Xtop"][1]
    charm = p["Pcharm"][0]+2*p["Pcharm"][1]
    loose_tree = tbound*xtop/lam**5 + cbound*charm/lam
    budgets = {}
    for tag, cap in [("mixing_matched",sm_bound),("broad_left_angle_contract",loose_tree)]:
        fbound = (rmax+cap)/h
        budgets[tag] = {"background_amplitude_cap":cap,
                        "right_Z_sd_modulus_cap":fbound,"right_W_product_cap":fu*fbound,
                        "extra_amplitude_needed_to_reach_charged_floor":h*target/fu-rmax-cap,
                        "equal_scale_normalization_min_fraction_to_keep_gap":(fu*fbound/target)**2}
        assert fu*fbound<target
    # Two independent doublet directions. A common u vector overlaps with each
    # orthogonal d,s vector, while their dot product is zero.
    a,b=.05,.03
    return {"normalizations":{"h":h,"pure_right_B_coefficient":kp0*1e-11*h*h},
            "two_width_measured_annulus_total_amplitude":[rmin,rmax],
            "neutral_limit_on_abs_imag_amplitude":sqrt(p["KOTO_90CL_e11"]/(klong-2*kl_err)),
            "max_neutral_over_whole_charged_range_e11":kl_max,
            "threshold_to_cut_central_charged_circle_e11":klong*bcentral/kp0,
            "standard_signed_center":[smamp.real,smamp.imag],
            "four_right_Fsd_for_same_nominal_SM_rates":mirror,
            "branches":budgets,"target_charged_product":target,
            "two_doublets_example":{"right_Z_ds":0,"right_W_u_d_and_s":a*b,
                                     "product":(a*b)**2,"up_diagonal":2*a*a}}


def formatted(data):
    return json.loads(json.dumps(data),parse_float=lambda x:float(f"{float(x):.11g}"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=formatted(solve())
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2)
            fp.write("\n")
    else:
        assert result==read("results.json")
        print("1461: rare charged and neutral kaons, signs and product budgets verified")
