"""Higgs mass-projector and soft hypercharge budget on the 1518 scale."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent

def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))

def required(gap_lambda, gauge, top, hmax, extra, ratio_y):
    bino,weak,_=gauge
    d=(weak**2+.6*bino**2)/8
    minimum_top=2*(top/hmax)**2-1
    minimum_abs=math.sqrt(max(0.,(gap_lambda-extra)/d))
    assert minimum_top>0 and minimum_abs<1
    a=max(minimum_top,minimum_abs)
    tan=math.sqrt((1+a)/(1-a))
    return dict(D=d,top_min_signed_response=minimum_top,
                min_abs_from_quartic=minimum_abs,min_joint_signed_response=a,
                min_up_fraction=(1+a)/2,min_tan_beta=tan,
                min_exchange_odd_operator_over_gap=a/2,
                min_common_incomplete_mass_error_over_gap=a/18,
                conditional_charged_Higgs_R_over_L_lower=ratio_y*tan**2)

def fixed_spectrum(asymmetry, gap, soft):
    # Only first two families of Q, E, D move. Products of masses stay invariant.
    r_minus_inv=math.sqrt(asymmetry*gap/(3*soft**2))
    r=(r_minus_inv+math.sqrt(r_minus_inv**2+4))/2
    shifted=np.array([r,1/r,1.])*soft
    baseline=np.ones(3)*soft
    masses=np.array([shifted,baseline,shifted,shifted,baseline])
    # Q, U, E, D, L: complex scalar coefficients in GUT-normalized order.
    beta=np.array([[.1,1.5,1],[.8,0,.5],[.6,0,0],[.2,0,.5],[.3,.5,0]])/3
    log_shifts=np.log(masses/soft).sum(axis=1)@beta/(2*math.pi)
    trace=(np.array([1,-2,1,1,-1])@np.square(masses)).sum()
    assert np.linalg.norm(log_shifts)<1e-14
    np.testing.assert_allclose(trace/gap,asymmetry)
    angle=math.asin(math.sqrt((1+asymmetry)/2))
    real_light=np.array([math.sin(angle),math.cos(angle)])
    heavy=np.array([math.cos(angle),-math.sin(angle)])
    mass_matrix=gap*np.outer(heavy,heavy)
    np.testing.assert_allclose((mass_matrix[1,1]-mass_matrix[0,0])/gap,asymmetry)
    np.testing.assert_allclose(np.linalg.eigvalsh(mass_matrix)/gap,[0,1],atol=2e-16)
    relabel=np.array([[0,1],[1,0]])
    distance=np.linalg.norm(mass_matrix-relabel@mass_matrix@relabel,2)/(2*gap)
    np.testing.assert_allclose(distance,asymmetry/2)
    return dict(moved_first_two_scalar_masses_GeV=shifted[:2].tolist(),
                soft_hypercharge_matter_sum_GeV2=trace,
                delta_inverse_gauge=log_shifts.tolist(),
                symmetric_phase_mass_matrix_GeV2=mass_matrix.tolist())

def calculate():
    cfg=read("inputs.json")
    previous=read(cfg["stop_data"])
    transport=read(cfg["top_transport"])["windows"]
    soft=cfg["comparison_soft_scale_GeV"]
    assert read(cfg["gauge_and_spectrum"])["witness"]["soft_GeV"]==soft
    cases={}
    for name,win,wide in [("central","1518_central",False),("outer_two_width","1518_slow_growth",True)]:
        w=transport[win]
        bound=w[cfg["endpoint_test"]]
        gap=previous["bounded_cases"]["1"]
        missing=gap["necessary_other_threshold_whole_input_box" if wide else "missing_central"]
        top0,terr=cfg["low_matched_yt_DR"]
        yt,ye=cfg["low_yt_MS_table"]
        yb,be=cfg["low_yb_MS_table"]
        n=2 if wide else 0
        arg=dict(gap_lambda=missing,gauge=w["gauge_at_low"],top=(top0-n*terr)/(1+cfg["Delta_t"]),
                 hmax=bound["quadrature_bracket"][1] if wide else bound["maximum_initial_top_for_endpoint"],
                 ratio_y=(yb-n*be)/(yt+n*ye))
        cases[name]={str(b):required(**arg,extra=b) for b in cfg["positive_extra_quartic_budgets"]}
    a=cases["central"]["0"]["min_joint_signed_response"]
    abox=cases["outer_two_width"]["0"]["min_joint_signed_response"]
    massgap=cfg["illustrative_gap_GeV2"]
    # For the hypercharge-independent aligned-family sum rule, subtract copies
    # of the explicitly gauge-universal field equations. The Yukawa traces are labelled.
    # Each row: Hu,Hd,U3-U1,D3-D1,E3-E1; columns: top,bottom,tau,(3/5)g1^2 S.
    beta_rows=np.array([[6,0,0,1],[0,6,2,-1],[4,0,0,0],[0,4,0,0],[0,0,4,0]])
    residual=np.array([-1,1,1.5,-1.5,-.5])@beta_rows
    np.testing.assert_allclose(residual,[0,0,0,-2])
    return dict(finite_budget_cases=cases,conditional_low_spectrum_example=fixed_spectrum(a,massgap,soft),
        necessary_abs_soft_trace_GeV2_for_registered_gap=abox*massgap,
        per_matrix_min_error_GeV2_if_all_share_tolerance=abox*massgap/18,
        first_RH_up_mass_lower_GeV_if_top_only_universal_and_U3_at_soft=math.sqrt(soft**2+2*abox*massgap/3),
        remaining_coefficients_before_adding_S_over_11=residual.tolist())

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            compare(a[key],b[key])
    else:
        np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-9)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    args=p.add_mutually_exclusive_group(required=True)
    args.add_argument("--write-new",action="store_true")
    args.add_argument("--check",action="store_true")
    options=p.parse_args()
    output=calculate()
    if options.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as h:
            json.dump(output,h,indent=2)
            h.write("\n")
    else:
        compare(output,read("results.json"))
        print("1524 OK: spectral response, soft trace, two data contracts and compensated gauge matching")
