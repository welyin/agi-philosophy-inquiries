"""Two neutral weak charges, one right-handed doublet, and held-out deuteron asymmetries."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def vertices(p,c):
    # Rows store a*x+b*y <= h. x,y are the first light state's RH weak fractions.
    forms=[(1,0,1),(-1,0,0),(0,1,1),(0,-1,0),
           (2,-1,p[1]),(-2,1,-p[0]),(188,-211,c[1]),(-188,211,-c[0])]
    points=[]
    for f,g in itertools.combinations(forms,2):
        matrix=np.array([f[:2],g[:2]])
        if abs(np.linalg.det(matrix))<1e-10:
            continue
        point=np.linalg.solve(matrix,[f[2],g[2]])
        if all(np.dot(row[:2],point)<=row[2]+2e-14 for row in forms):
            points.append(np.clip(point,0,1))
    assert points
    return points


def sharp_upper(p,c):
    x=(211*p-c)/234
    y=(188*p-2*c)/234
    assert 0<=x<=1 and 0<=y<=1
    return x,y,math.sqrt(x*y)


def charged_completion(x,y,common=1500.):
    # An explicit single-H renormalizable left-column alignment for each right mixing.
    data=[]
    for fraction,light_mass in [(x,.002),(y,.005)]:
        heavy=math.sqrt((common**2-light_mass**2*fraction)/(1-fraction))
        a=np.array([math.sqrt(fraction),math.sqrt(1-fraction)])
        rotation=np.array([[a[1],a[0]],[-a[0],a[1]]])
        vector=np.array([light_mass,heavy])*a/common
        left=np.array([[vector[1],-vector[0]],[vector[0],vector[1]]])
        mass=left@np.diag([light_mass,heavy])@rotation
        assert np.allclose(mass[:,1],[0,common],atol=1e-9)
        data.append(heavy)
    return data


def deuteron(sum_xy,data):
    window=[]
    for (b1,b3),asm,sm_err,obs in zip(data["table21_valence_derivatives_ppm"],
                                      data["table21_MSTW_SM_intercept_ppm"],
                                      data["illustrative_remaining_SM_intercept_budget_ppm"],
                                      data["table18_observed_ppm_stat_sys"]):
        s2=data["electron_sin2_for_new_axial_response"]
        eps=data["additional_absolute_response_error_as_fraction_of_sum_of_term_magnitudes"]
        slopes=[]
        for angle in s2:
            k=.5-2*angle
            nominal=-b1/2+b3*k
            e=eps*(b1/2+b3*abs(k))
            slopes.extend([nominal-e,nominal+e])
        offsets=[slope*s for slope in slopes for s in sum_xy]
        observed_err=math.hypot(obs[1],obs[2])
        window.append({"new_shift_ppm":[min(offsets),max(offsets)],
                       "shift_plus_baseline_ppm":[asm+min(offsets)-sm_err,asm+max(offsets)+sm_err],
                       "direct_observed_one_width_ppm":[obs[0]-observed_err,obs[0]+observed_err]})
    return window


def compute():
    cfg=read("inputs.json")
    hist=read(cfg["historical_p_and_Cs"])
    qp,qp_err=hist["proton"]["pv_electron_scattering_only"]
    stdp,stdp_err=hist["proton"]["SM"]
    stdc,stdc_err=cfg["baseline_Cs"]
    results={}
    for model,(qc,exp_err,atom_err) in cfg["atomic_branches"].items():
        rows=[]
        for w in cfg["width_multipliers"]:
            ep=w*(qp_err+stdp_err)
            ec=w*(exp_err+atom_err+stdc_err)
            p=(qp-stdp-ep,qp-stdp+ep)
            c=(qc-stdc-ec,qc-stdc+ec)
            pts=vertices(p,c)
            x,y,upper=sharp_upper(p[1],c[0])
            assert any(np.linalg.norm(pt-[x,y])<1e-11 for pt in pts)
            upper_check=max(math.sqrt(a*b) for a,b in pts)
            lower=min(math.sqrt(a*b) for a,b in pts)
            assert math.isclose(upper,upper_check,rel_tol=1e-11) and lower<1e-10
            combo=[2*z[0]+z[1] for z in pts]
            s_bounds=[min(combo),max(combo)]
            forecast=sharp_upper(ep,-ec)[2]
            rows.append({"width_multiplier":w,"neutral_deviation_rectangle":[p,c],
                         "abs_right_ud_bounds":[lower,upper],"attaining_powers":[x,y],
                         "example_heavy_masses_GeV":charged_completion(x,y),
                         "two_x_plus_y_bounds":s_bounds,
                         "deuteron_holdouts":deuteron(s_bounds,cfg["PVDIS"]),
                         "zero_center_upper_with_same_widths":forecast,
                         "equal_improvement_factor_for_absR_1e_minus3":.001/forecast})
        results[model]=rows
    # For any nucleus, shifts induced by a common weak-angle displacement
    # plus x=2y vanish for delta sin^2=3y/4, leaving a nonzero right W.
    a,b,z=.04,.02,.015
    assert abs(2*a-b-4*z)<1e-14 and abs(a-2*b)<1e-14
    results["free_weak_angle_witness"]={"x":a,"y":b,"delta_sin2":z,"abs_right_ud":math.sqrt(a*b)}
    x,y,r=sharp_upper(.001,.025)
    s=2*x+y
    direct=[]
    for (d1,d3),sm,obs in zip(cfg["PVDIS"]["table21_valence_derivatives_ppm"],
                               cfg["PVDIS"]["table21_MSTW_SM_intercept_ppm"],
                               cfg["PVDIS"]["table18_observed_ppm_stat_sys"]):
        prediction=sm+(-d1/2+d3*(.5-2*.235))*s
        assert abs(prediction-obs[0])<=math.hypot(obs[1],obs[2])
        direct.append(prediction)
    for model in cfg["atomic_branches"]:
        pre,cs=results[model][0]["neutral_deviation_rectangle"]
        assert pre[0]<=.001<=pre[1] and cs[0]<=.025<=cs[1]
    results["single_witness_for_both_DIS_bins"]={"p":.001,"c":.025,"abs_right_ud":r,
                                                "two_nominal_asymmetries_ppm":direct}
    return results


def same(a,b):
    if isinstance(a,dict) and isinstance(b,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,(list,tuple)) and isinstance(b,(list,tuple)):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-11,abs_tol=2e-13)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=compute()
    path=HERE/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            f.write(json.dumps(data,indent=2,allow_nan=False)+"\n")
    else:
        assert same(data,read("results.json"))
    print("PASS neutral-only inverse, vector-doublet completion, and direct DIS holdouts")
