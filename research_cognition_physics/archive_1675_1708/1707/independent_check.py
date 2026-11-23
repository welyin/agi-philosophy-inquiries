"""Invert the Bd CP phase and Bd/Bs clock ratio using exact CKM triangle lengths."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"
DATA = {"golden_S":[.717,.013,.008], "xi_4sea":[1.216,.016],
        "dm_d_inverse_ps":[.5050,.0021,.0010], "dm_s_inverse_ps":[17.7683,.0051,.0032],
        "mass_d_MeV":[5279.74,.30,.10], "mass_s_MeV":[5366.85,.19,.13]}
URLS = ["https://arxiv.org/abs/2309.09728","https://arxiv.org/abs/1907.01025",
        "https://arxiv.org/abs/1604.03475","https://arxiv.org/abs/2104.04421",
        "https://doi.org/10.1103/PhysRevLett.122.191804"]


def shifted(value, width):
    return value[0]+width*math.hypot(*value[1:])


def ratio_at(t, xi=None):
    # Simultaneously interpolates between the extrema of the independently widened input ranges.
    row = {k: shifted(v, -t if k in ("dm_s_inverse_ps", "mass_d_MeV") else t)
           for k,v in DATA.items() if k != "golden_S"}
    if xi is not None:
        row["xi_4sea"] = shifted(xi, t)
    ratio = row["xi_4sea"]*math.sqrt(row["dm_d_inverse_ps"]*row["mass_s_MeV"]/
                                   (row["dm_s_inverse_ps"]*row["mass_d_MeV"]))
    return ratio, row


def inverse_roots(us, cb, ratio, golden):
    s,b = us**2,cb**2
    h = ratio**2/(1+ratio**2)
    u = Polynomial([0.,1.])
    q = u+b
    M = b*(u+s-h*q)
    T = h*q*(1-q)
    W = u*(1-u-s)
    C = M+T-W
    beta1 = math.asin(golden)/2
    rows = []
    for beta in (beta1, math.pi/2-beta1):
        equation = C*C - 4*math.cos(beta)**2*M*T
        brackets = ([(1e-5,2e-5),(.0003,.00035),(.941,.943),(.949,.951)] if beta==beta1
                    else [(.00009,.00011),(.00021,.00025),(.944,.945),(.947,.948)])
        # These disjoint sign changes exhaust the degree-four polynomial at the central inputs.
        assert equation.degree()==4 and all(equation(a)*equation(b)<0 for a,b in brackets)
        for root in equation.roots():
            if abs(root.imag)>1e-9 or not (0<root.real<min(1-s,1-b)):
                continue
            v = float(root.real)
            if min(M(v),T(v),W(v),C(v))<=0:
                continue
            m,t,w = [math.sqrt(poly(v)) for poly in (M,T,W)]
            if abs((m*m+t*t-w*w)/(2*m*t)-math.cos(beta))>1e-7:
                continue
            gamma = math.acos((m*m+w*w-t*t)/(2*m*w))
            rows.append({"beta_deg":math.degrees(beta),"Vub":math.sqrt(v),
                       "Vub_over_Vcb":math.sqrt(v)/cb,"gamma_deg":math.degrees(gamma)})
    assert len(rows)==4
    return rows


def from_tree(us, ub, cb, gamma_deg):
    gamma=math.radians(gamma_deg)
    c13=math.sqrt(1-ub**2)
    s12,s23=us/c13,cb/c13
    c12,c23=math.sqrt(1-s12**2),math.sqrt(1-s23**2)
    # Multiplication by a real positive number does not affect beta.
    phase=gamma+math.asin(c12/s12*s23/c23*ub*math.sin(gamma))
    p=cmath.exp(1j*phase)
    cd=-s12*c23-c12*s23*ub*p
    td=s12*s23-c12*c23*ub*p
    ts=-c12*s23-s12*c23*ub*p
    beta=cmath.phase(-cd/td)
    actual_gamma=cmath.phase(-c12*c13*ub*p/(cd*cb))
    assert abs(actual_gamma-gamma)<2e-13
    return math.sin(2*beta),abs(td/ts)


def check_point(label,row,ub,us,cb,interval,sgate,ggate):
    gamma=row["gamma_MZ_deg"]
    pred_s, pred_r=from_tree(us,ub,cb,gamma)
    want=min(max(pred_r,interval[0]),interval[1])
    need=(want/pred_r)**2
    eps=abs(need-1)/(need+1)
    item={"label":label,"Vub":ub,"gamma_deg":gamma,"golden_S":pred_s,"Vtd_over_Vts":pred_r,
          "gamma_inside_one_width":ggate[0]<=gamma<=ggate[1],
          "S_inside_one_width":sgate[0]<=pred_s<=sgate[1],
          "frequency_lattice_inside_one_width":interval[0]<=pred_r<=interval[1],
          "minimum_per_flavour_relative_clock_budget":eps}
    if item["frequency_lattice_inside_one_width"]:
        lo,hi=-1.,1.
        for _ in range(80):
            mid=(lo+hi)/2
            if ratio_at(mid)[0]<pred_r:lo=mid
            else:hi=mid
        r, data=ratio_at((lo+hi)/2)
        assert abs(r-pred_r)<1e-13
        item["clock_witness_within_one_width"]=data
    return item


def compute():
    r2=json.loads((HERE/"../1702/results.json").read_bytes())
    cfg=json.loads((HERE/"../1702/inputs.json").read_bytes())["direct_semileptonic_moduli"]
    us,ub,cb=(cfg[k][0] for k in ("Vus","Vub","Vcb"))
    center=ratio_at(0)[0]
    one=[ratio_at(t)[0] for t in (-1,1)]
    two=[ratio_at(t)[0] for t in (-2,2)]
    sg=[shifted(DATA["golden_S"],t) for t in (-1,1)]
    gg=[61.,65.5]
    central=inverse_roots(us,cb,center,DATA["golden_S"][0])
    # Three input widths and their correlations across determinations are not combined into a likelihood.
    points = []
    for k,row in r2["charged_endpoints"]["direct"].items():
        points.append(check_point("central_k"+k,row,ub,us,cb,one,sg,gg))
    for k,row in r2["down_one_width_lower"].items():
        points.append(check_point("lower_down_"+k,row,ub,us,cb,one,sg,gg))
    for k,row in r2["independent_Vub_sensitivity"].items():
        offset=int(k.rsplit("_",1)[1])
        points.append(check_point(k,row,shifted(cfg["Vub"],offset),us,cb,one,sg,gg))
    passing=[p["label"] for p in points if p["gamma_inside_one_width"] and
             p["S_inside_one_width"] and p["frequency_lattice_inside_one_width"]]
    assert passing==["k2_Vub_widths_-1"]
    alternative=next(p for p in points if p["label"]=="k1_Vub_widths_-1")
    assert 0<alternative["minimum_per_flavour_relative_clock_budget"]<.01
    uslo=shifted(cfg["Vus"],-2)
    cblo,cbhi=[shifted(cfg["Vcb"],k) for k in (-2,2)]
    ubhi=shifted(cfg["Vub"],2)
    s23max=cbhi/math.sqrt(1-ubhi**2)
    cdlo=uslo*math.sqrt(1-s23max**2)-s23max*ubhi
    rbhi=ubhi/(cdlo*cblo)
    min_big_sin=math.cos(math.asin(shifted(DATA["golden_S"],2))/2)
    assert min_big_sin>rbhi
    fnal_milc=[1.206,.018,.006]
    fnal_one=[ratio_at(t, fnal_milc)[0] for t in (-1,1)]
    fnal_points=[p["label"] for p in points if p["gamma_inside_one_width"] and p["S_inside_one_width"]
                 and fnal_one[0]<=p["Vtd_over_Vts"]<=fnal_one[1]]
    assert fnal_points==["k1_Vub_widths_-1","k2_Vub_widths_-1"]
    return {"measurements_in_order_as_source_urls":DATA,"source_urls":URLS,
            "SM_clock_ratio_Vtd_over_Vts":center,"one_width_clock_ratio_range":one,
            "two_width_clock_ratio_range":two,"one_width_golden_S_range":sg,
            "alternative_FNAL_MILC_2plus1_with_omitted_charm_budget":{"xi":fnal_milc,
                "source":"https://arxiv.org/abs/1602.03560", "one_width_ratio":fnal_one,
                "joint_one_width_pass_without_extra_clock_residual":fnal_points},
            "one_width_direct_gamma_deg_range":gg,"all_central_physical_triangle_inversions":central,
            "two_width_mirror_gate":{"Rb_upper":rbhi,"large_beta_sine_lower":min_big_sin,
                 "direct_Vub_upper":ubhi,"direct_Vcd_lower":cdlo},
            "frozen_direct_tree_points":points,"one_width_joint_pass":passing,
            "one_width_gamma_S_clock_with_one_percent_each_flavour":[p["label"] for p in points if
             p["gamma_inside_one_width"] and p["S_inside_one_width"] and
             p["minimum_per_flavour_relative_clock_budget"]<=.01]}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=compute()
    if args.check:
        assert data==json.loads(OUT.read_bytes())
        print("1707 independent two-clock inversion and flavour-budget checks passed")
    else:
        OUT.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")

