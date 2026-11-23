"""Exact chiral triangular matching and independent leptonic tests."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8"))


def weak_ratio(active_sq, shat):
    return ((shat-active_sq/2)**2+shat**2)/((shat-.5)**2+shat**2)


def witness(T,m,H,shat,vev):
    # Choose the off-diagonal, diagonal heavy mass by the exact trace and determinant.
    bare_m=m/T
    bare_M=T*H
    mix=math.sqrt((1-T*T)*(H*H-bare_m*bare_m))
    output={}
    for name,mat in [("higgs",np.array([[bare_m,mix],[0.,bare_M]])),
                     ("soft",np.array([[bare_m,0.],[mix,bare_M]]))]:
        left, masses, right=np.linalg.svd(mat)
        assert np.allclose(masses,[H,m],atol=1e-12)
        l=abs(left[0,-1])**2
        r=abs(right[-1,0])**2
        h=left[:,-1] @ (np.diag([1.,0.])@mat) @ right[-1]/m
        assert abs(h-l)<1.e-11
        expected=(1-T*T)/(1-(m/H)**2)
        assert abs((1-(l if name=="higgs" else r))-expected)<1e-12
        if name=="soft":
            assert np.isclose(abs(left[1,-1])**2,(m/H/T)**2*expected,rtol=1e-7,atol=1e-20)
        output[name]={"charged_current_mu_over_e":math.sqrt(l),
                      "relative_Z_width":weak_ratio(l,shat),"h_mumu_over_m_by_v":l,
                      "left_singlet_probability":1-l,"right_heavy_label_probability":1-r}
    return {"matrix_diagonals_GeV":[bare_m,bare_M],"offdiag_GeV":mix,
            "higgs_offdiag_yukawa_if_upper":math.sqrt(2)*mix/vev,"branches":output}


def evaluate():
    cfg=read("inputs.json")
    tex=read(cfg["texture"])["textures"][cfg["texture_key"]]
    matched=read(cfg["running_matching"])
    lo,hi=matched["T_box_with_transport"]
    center=matched["T_QED"]
    belle=read(cfg["tau"])["Belle_II"]
    n=cfg["widths"]
    b=belle["ratio_central"]
    bel=[b-n*belle["ratio_width"],b+n*belle["ratio_width"]]
    se,se_err=read(cfg["electron_pole"])["SLD_sin2_eff_e"]
    slo,shi=se-n*se_err,se+n*se_err
    e,ee=cfg["lep"]["R_e"]
    u,eu=cfg["lep"]["R_mu"]
    z=[(e-n*ee)/(u+n*eu),(e+n*ee)/(u-n*eu)]
    # Conservative all-angle enclosure for 0<=|l_SM|²<=hi².
    worst_left=max(shi,abs(slo-hi*hi/2),abs(shi-hi*hi/2))
    zmax=(worst_left**2+shi**2)/((shi-.5)**2+shi**2)
    fmin,fmax=cfg["amplitude_or_rate_transfer"]
    assert hi*fmax<bel[0] and zmax*fmax<z[0]
    vev=read(cfg["MZ_vev"])["matching"]["v_GeV"][0]
    mass_mu=read(cfg["MZ_yukawa"])["yukawa"]["charged_leptons"][1][0]*vev/math.sqrt(2)
    point=witness(center,mass_mu,cfg["illustrative_physical_M_GeV"],se,vev)
    legacy=witness(tex["real_minus_witness"]["lepton_ratio_transfer"],cfg["mu_mass_GeV"],
                   cfg["illustrative_physical_M_GeV"],se,cfg["example_vev_GeV"])
    soft=point["branches"]["soft"]
    assert bel[0]<soft["charged_current_mu_over_e"]<bel[1]
    assert z[0]<soft["relative_Z_width"]<z[1]
    for T in [lo,hi]:
        for mass in [50.,cfg["illustrative_physical_M_GeV"],1000.]:
            witness(T,mass_mu,mass,se,vev)
    return {"T_center":center,"T_window":[lo,hi],"Belle_tau_two_widths":bel,
            "LEP_Z_mu_e_two_widths":z,"higgs_only": {
                "bare_W_max":hi,"bare_Z_max_conservative":zmax,
                "minimum_additional_factor_W":bel[0]/hi,
                "minimum_additional_factor_Z":z[0]/zmax,
                "weak_singlet_leak_probability_floor":1-hi*hi},
            "MZ_mass_example_mu_GeV":mass_mu,"T_window_pure_QED":matched["T_box_pureQED"],
            "common_tree_singular_mass_example":point,
            "old_pole_ratio_algebraic_example_only":legacy}


def audit():
    def clean(x):
        if isinstance(x,dict): return {k:clean(v) for k,v in x.items()}
        if isinstance(x,list): return [clean(v) for v in x]
        return float(f"{x:.12g}")
    return clean(evaluate())


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("1488: chiral transfer and both measured weak responses verified.")
    else:
        print(json.dumps(out, indent=2))
