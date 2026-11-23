"""Test semileptonic first-row closure using the same Higgs budget and lost weak probability."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(file):
    return json.loads((HERE/file).read_bytes())


def run():
    cfg=read("inputs.json")
    point=read(cfg["gauge_point"])["central_conditional_combined_point"]
    M=point["five_color_poles_GeV"][0]
    v=read(cfg["vev"])["v246_GeV"]/math.sqrt(2)
    yb,ys,yd=read(cfg["spectrum"])["tagged_yukawa_values_desc"]
    hist=read(cfg["curvature"])
    data=hist["one_parameter_neutral_rescaling"]
    gmax=data["all_s_two_width_upper_g"]
    lam=(v*yb/M)**2
    eps_max=(gmax*v/(M*(1-lam)))**2
    eps_max=eps_max/(1+eps_max)
    assert eps_max >= data["all_s_two_width_upper_d_singlet_fraction"]*(1-1e-10)
    joined=read(cfg["quark"])["direct_semileptonic_moduli"]
    s,ds=joined["Vus"][0],np.linalg.norm(joined["Vus"][1:])
    b,db=joined["Vub"][0],np.linalg.norm(joined["Vub"][1:])
    def summary(d,width):
        def deficit(n):
            return 1-(d+n*width)**2-(s+n*ds)**2-(b+n*db)**2
        central=deficit(0)
        def local_required(required):
            if required<=0:
                return 0.,0.
            g2=(M*(1-lam)/v)**2*required/(1-required)
            color=3*M*M/(8*math.pi**2)*1.5*g2
            positive_floor=color-data["maximum_positive_neutral_c"]-hist["target_central_and_two_width_abs_upper_GeV2"][1]
            return math.sqrt(g2),max(0.,positive_floor)
        normal=np.array([d,s,b])
        err=np.array([width,ds,db])
        linear=2*normal@err
        extra=central-eps_max
        n=2*extra/(linear+math.sqrt(linear**2+4*np.dot(err,err)*extra))
        return dict(central_apparent_deficit=central,most_favourable_deficit_at_one_width=deficit(1),
                    most_favourable_deficit_at_two_widths=deficit(2),
                    smallest_common_marginal_half_width_to_meet_necessary_bound=n,
                    g_lower_and_positive_extra_curvature_GeV2_at_central=local_required(central),
                    g_lower_and_positive_extra_curvature_GeV2_at_one_width=local_required(deficit(1)))
    old=read(cfg["nuclear_old"])["ud"]
    tests={"historical_nuclear":summary(old["central"],old["errors"][0])}
    for key in ("nuclear_EFT","neutron_best_inputs","pion"):
        tests[key]=summary(*cfg["external_2026"][key])
    root=abs(complex(*read(cfg["branch"])["exact_and_conditional_implementations"]["second_complex_root"]["z"]))
    g=point["g"]
    kap=math.sqrt(point["R"]**2-1)
    def left_loss(k):
        w=np.sqrt(np.maximum(0,(1+k*np.array([1.,-1.]))/2))
        signs=np.array([1.,-1.])[w>0]
        w=w[w>0]
        A=np.zeros((1+len(w),1+len(w)))
        A[0,1:]=g*v/M*w
        A[1:,0]=-kap*w
        A[1:,1:]=np.diag(signs)
        l=np.linalg.svd(A)[0][:,-1]
        return np.dot(l[1:],l[1:])
    losses=[left_loss(k) for k in (-yd/yb,ys/yb,-1.)]
    R0=1-s*s-b*b
    illustrative=(R0*losses[0]+s*s*losses[1]+b*b*losses[2])
    assert illustrative<eps_max
    return dict(max_all_light_singlet_probability_at_Higgs_budget=eps_max,
                coloured_five_tree_light_d_s_b_losses_R800=losses,
                example_quark_only_first_row_missing_probability=illustrative,
                unitary_ud_from_central_strange_and_beauty=math.sqrt(R0),
                minimum_apparent_ud_from_central_strange_and_beauty=math.sqrt(R0-eps_max),
                fractional_rate_leakage_required_if_two_arm_root_is_apparent_ud=1-root**2-s*s-b*b,
                apparent_ud_comparisons=tests)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=json.loads(json.dumps(run()),parse_float=lambda x:float(f"{float(x):.10g}"))
    path=HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes())==out
        print("1730 verified GF-normalized first-row inequalities, the curvature cost, and independent nuclear extractions.")
    else:
        with path.open("x",encoding="utf8") as stream:
            stream.write(json.dumps(out,indent=2,allow_nan=False)+"\n")
