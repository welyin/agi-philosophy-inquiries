"""Match all five saved flavour solutions to low-energy quark mass ratios."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
HAD = runpy.run_path(str(HERE / CFG["mass_boxes"]))
CKM = runpy.run_path(str(HERE / CFG["physical_ckm"]))["ckm"]


def fractional_range(value, center, box, n):
    return [value/(center+n*(box[1]-center)),value/(center+n*(box[0]-center))]


def expand(center, box, n):
    return [center+n*(v-center) for v in box]


def light_projection(width, both):
    """Sharp range of D=d/s for the low-energy data with u free."""
    h = HAD
    s0,s1 = expand(h["CFG"]["MILC_S"][0],h["SS"],width)
    r0,r1 = expand(h["CFG"]["MILC_light_ratio"][0],h["RR"],width)
    z0,z1 = (1-r1)/(1+r1),(1-r0)/(1+r0)
    if both:
        q0,q1 = expand(h["ETA"][0],h["QB"],width)
        b0,b1 = expand(h["LIN"][0],h["SB"],width)
        s0 = max(s0,b0,math.sqrt(1+4*q0*q0*z0))
        s1 = min(s1,b1,math.sqrt(1+4*q1*q1*z1))
        if s0>s1:
            return None
        assert s1*s1+1 < 4*q0*q0
        return [(1+max(z0,(s1*s1-1)/(4*q1*q1)))/s1,
                (1+min(z1,(s0*s0-1)/(4*q0*q0)))/s0]
    return [(1+z0)/s1,(1+z1)/s0]


def source_thresholds():
    ref = CFG["SMDR_MZ"]
    eft = np.array(ref["eft_b_s_d_u_MeV"])
    sm = np.array(ref["sm_b_s_d_u_MeV"])
    rnd = np.array(ref["rounding_radii_MeV"])
    fac = sm/eft
    lower, upper = (sm-rnd)/(eft+rnd), (sm+rnd)/(eft-rnd)
    bs_rounding = [lower[0]/upper[1],upper[0]/lower[1]]
    yt = math.sqrt(2)*ref["table_mass_top"]/ref["mass_table_vF"]
    t = yt**2*ref["running_v"]**2/2
    w = (ref["g2"]*ref["running_v"]/2)**2
    a = w/(t-w)
    h = 0.75*yt**2*(5/6-math.log(t/ref["MZ"]**2)+a*a*math.log(t/w)-a)/(16*math.pi**2)
    src = CFG["SMDR_200"]
    low = np.array(src["low_ds_u_MeV"])
    y = np.array(src["y_d_s_u"])
    ds,ud = y[0]/y[1]*low[1]/low[0], y[2]/y[0]*low[0]/low[2]
    # d ln(yu/yd)/dt = -gY^2/(16 pi^2) in the light, diagonal limit.
    by = 41/6
    gy = src["hypercharge_running"]
    inv_gz2 = 1/gy**2-2*by*math.log(src["MZ"]/src["Q"])/(16*math.pi**2)
    ud_at_Z = ud*(math.sqrt(1/inv_gz2)/gy)**(-1/by)
    assert abs(ds-1)<2e-6
    return {"ds_at_200_over_low":ds,"ud_at_200_over_low":ud,"ud_at_Z_light_approx":ud_at_Z,
            "HZ_bs_full_over_eft":fac[0]/fac[1],"HZ_bs_print_rounding_range":bs_rounding,
            "HZ_ds_full_over_eft_printed":fac[2]/fac[1],"HZ_ud_full_vs_2GeV_printed":sm[3]/sm[2]*low[0]/low[2],
            "HZ_ud_print_rounding_range":[(sm[3]-rnd[3])/(sm[2]+rnd[2])*low[0]/low[2],
                                           (sm[3]+rnd[3])/(sm[2]-rnd[2])*low[0]/low[2]],
            "top_minus_light_one_loop_in_decoupling":h}


def thresholds(witness, t):
    weak = witness["moduli_MZ"]
    ratio = CKM(witness["weak_offsets"],witness["gamma"])["ratio"]**2
    pb = 1-weak["Vub"]**2-weak["Vcb"]**2
    ps = (1-pb)/(1+ratio)
    pd = ratio*ps
    h = t["top_minus_light_one_loop_in_decoupling"]
    # High / low convention. The published CKM-diagonal mass ratio already
    # includes the top correction on b. Add only the orientation difference.
    return np.array([1+(ps-pd)*h, t["ud_at_Z_light_approx"]*(1+pd*h),
                     t["HZ_bs_full_over_eft"]+(ps-pb+1)*h]), [pd,ps,pb]


def qcdqed_tests(ds, ud, n):
    data = CFG["direct_QCD_QED"]
    xs, es = data["QCDSF_s_over_d"]
    xd, ed = data["QCDSF_u_over_d_Dashen"]
    f, _ = data["QCDSF_ratio_Dashen_to_MS2"]
    # Conservative lowest centre AND envelope for a downward <1% conversion.
    dashen_pull = (ud/f-xd)/ed
    def pull(row):
        x,*errors=row
        return (ud-x)/math.hypot(*errors)
    p = {"QCDSF_s_over_d":(1/ds-xs)/es,"QCDSF_u_over_d_most_favorable":dashen_pull,
         "ETM_u_over_d":pull(data["ETM_u_over_d"]),"BMW_u_over_d":pull(data["BMW_u_over_d"])}
    return {"s_over_d_pred":1/ds,"u_over_d_pred":ud,"signed_widths":p,
            "in_original_n":{key:abs(v)<=n for key,v in p.items()}}


def analyze(w, theory):
    n = w["largest_width_in_named_box"]
    raw = HAD["observables"](w["up_MZ"],w["down_MZ"])
    K, prob = thresholds(w,theory)
    result = {"width":n, "Vt_probabilities":prob, "perturbative_K_ds_ud_bs":K,
              "one_loop_Kbs":1+(prob[1]-prob[2])*theory["top_minus_light_one_loop_in_decoupling"]}
    for label, both in (("MILC_only",False),("all_light",True)):
        target = light_projection(n,both)
        needed = raw[0]/np.array(target[::-1])
        result[label+"_Kds_range"] = needed
        result[label+"_minimum_residual_ds"] = needed[0]/K[0]-1
    result["heavy_Kbs_range"] = fractional_range(raw[2],HAD["HM"][0],HAD["BS"],n)
    bs_low,bs_high = result["heavy_Kbs_range"]
    result["heavy_contains_perturbative_central"] = bool(bs_low<=K[2]<=bs_high)
    result["direct_QCD_QED"] = qcdqed_tests(raw[0]/K[0],raw[1]/K[1],n)
    S = expand(HAD["CFG"]["MILC_S"][0],HAD["SS"],n)
    bmw = expand(HAD["LIN"][0],HAD["SB"],n)
    S = (max(S[0],bmw[0]),min(S[1],bmw[1]))
    constraints = dict(Sbox=S,rbox=expand(HAD["CFG"]["MILC_light_ratio"][0],HAD["RR"],n),
                       Qbox=expand(HAD["ETA"][0],HAD["QB"],n),bsbox=expand(HAD["HM"][0],HAD["BS"],n))
    # Change coordinates in R1734's exact physical reachability test, whose
    # reference is internal to the old function. Only the raw/reference ratio enters.
    relative = raw*HAD["MATCH"]/K
    eps = HAD["threshold"](relative,**constraints)
    common = HAD["feasible"](relative,eps+1e-9,**constraints)
    factors = raw/(K*np.array([common["d_over_s"],common["u_over_d"],common["b_over_s"]]))
    assert max(abs(factors-1))<=eps+1.1e-9
    result.update(min_joint_finite_residual=eps, joint_witness=common,finite_multipliers_at_witness=factors)
    assert result["all_light_minimum_residual_ds"]>.011
    assert result["direct_QCD_QED"]["in_original_n"]["QCDSF_s_over_d"]
    assert not result["direct_QCD_QED"]["in_original_n"]["QCDSF_u_over_d_most_favorable"]
    return result


def main():
    th = source_thresholds()
    points = json.loads((HERE/CFG["points"]).read_bytes())["witnesses"]
    sols = {name:analyze(row,th) for name,row in points.items()}
    assert len(sols)==5
    assert sum(row["heavy_contains_perturbative_central"] for row in sols.values())==4
    return HAD["tidy"]({"thresholds":th,"same_1736_witnesses":sols})


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=main()
    path=HERE/"results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1737 checked low-energy ratio projections and weak threshold transfers at five fixed points")
    else:
        path.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n",encoding="utf8")
