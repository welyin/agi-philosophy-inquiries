"""Strong isospin and heavy-meson holdouts of the charged three-link solution."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
def read(path): return json.loads((HERE/path).read_bytes())
def api(path): return runpy.run_path(str(HERE/path))
CFG=read("inputs.json")
Y=read(CFG["running_inputs"])["yukawa"]
UP=np.array(Y["up"][0])
FROZEN=read(CFG["solutions"])
key=CFG["primary"]
LIN=np.array(CFG["BMW_S"])
ETA=np.array(CFG["eta_Q"])
HM=np.array(CFG["MILC_b_over_s"])

def bounds(c,low,high=None):
    return c-low,c+(low if high is None else high)

def mass_interval(x):
    c,stat,lo,hi,aux,em=x
    return bounds(c,math.hypot(stat,lo,aux,em),math.hypot(stat,hi,aux,em))

def reference():
    low=CFG["running_low_MeV"]
    d,s,b=[r[0] for r in Y["down"]]
    b4,b5=CFG["MILC_conversion_b4_at2_b5_atmb_MeV"]
    return np.array([(d/s)/(low["d"]/low["s"]),(UP[0]/d)/(low["u"]/low["d"]),
                     (b/s)/(low["b5_at_mb"]*b4/b5/low["s"])])

MATCH=reference()

def masses(S,Q):
    delta=(S*S-1)/(4*Q*Q)
    return (1+delta)/S,(1-delta)/(1+delta)

SS=mass_interval(CFG["MILC_S"])
RR=mass_interval(CFG["MILC_light_ratio"])
BS=bounds(*HM)
SB=bounds(LIN[0],np.linalg.norm(LIN[1:]))
QB=bounds(*ETA)

def observables(u,dsb):
    d,s,b=dsb
    return np.array([d/s,u/d,b/s])

def hadrons(raw,reference=MATCH):
    d,r,B=raw/reference
    S=2/(d*(1+r))
    Q=.5*math.sqrt((S*S-1)*(1+r)/(1-r))
    return {"d_over_s":d,"u_over_d":r,"S":S,"Q":Q,"b_over_s":B}

def feasible(raw,budget,Sbox=None,rbox=None,Qbox=None,bsbox=None):
    """Exact one-variable elimination. Compare each lower quadratic to each upper."""
    ranges=[(x/m/(1+budget),x/m/(1-budget)) for x,m in zip(raw,MATCH)]
    (dl,dh),(rl,rh),(bl,bh)=ranges
    if bsbox and (bl>bsbox[1] or bh<bsbox[0]):
        return None
    if rbox:
        rl,rh=max(rl,rbox[0]),min(rh,rbox[1])
    if rh<rl:
        return None
    zl,zh=(1-rh)/(1+rh),(1-rl)/(1+rl)
    slo,shi=Sbox if Sbox else (1.001,60.)
    lower=[np.array([zl,0.,0.]),np.array([-1.,dl,0.])]
    upper=[np.array([zh,0.,0.]),np.array([-1.,dh,0.])]
    if Qbox:
        for Q,lines in ((Qbox[1],lower),(Qbox[0],upper)):
            c=1/(4*Q**2)
            lines.append(np.array([-c,0.,c]))
    cuts=[slo,shi]
    for L in lower:
        for H in upper:
            for x in np.polynomial.polynomial.polyroots(L-H):
                if abs(x.imag)<1e-12 and slo<=x.real<=shi:
                    cuts.append(x.real)
    cuts.sort()
    for S in cuts+[(a+b)/2 for a,b in zip(cuts,cuts[1:])]:
        xx=np.array([1.,S,S*S])
        a=max(l@xx for l in lower)
        b=min(h@xx for h in upper)
        if a<=b+2e-14:
            z=(a+b)/2
            return {"S":S,"u_over_d":(1-z)/(1+z),"Q":math.sqrt((S*S-1)/(4*z)),
                    "d_over_s":(1+z)/S,"b_over_s":max(bsbox[0],bl) if bsbox else bl}
    return None

def threshold(raw,**kw):
    if feasible(raw,0.,**kw):
        return 0.
    lo,hi=0.,.30
    assert feasible(raw,hi,**kw)
    for _ in range(45):
        mid=(lo+hi)/2
        if feasible(raw,mid,**kw):
            hi=mid
        else:
            lo=mid
    return hi

def weak_case(yu,step,include_neutral):
    # Both readers own separate model instances. Set the same physical up input before
    # solving the charged boundary, its neutral lifting and the 10 PeV Schur kernel.
    flavor=api(CFG["model"])
    flavor["PARENT"]["FLOW"]["low"]["yukawa"]["up"][0][0]=float(yu)
    row=flavor["solve"](3,key,step)
    mt=observables(yu,row["down_strange_bottom_MZ"])
    out={"u_MZ":yu,"up_shift_widths":(yu-UP[0])/UP[1],"dsb_MZ":row["down_strange_bottom_MZ"],
         "charged_widths":row["down_strange_bottom_marginal_widths"],"gamma":row["gamma_MZ_deg"],
         "hadron":hadrons(mt),"matching_relative_budget_needed_for_BMW_eta_and_heavy":threshold(mt,Sbox=SB,Qbox=QB,bsbox=BS),
         "matching_relative_budget_needed_for_MILC_S_ud_eta_and_heavy":threshold(mt,Sbox=SS,rbox=RR,Qbox=QB,bsbox=BS)}
    if include_neutral:
        neutral=api(CFG["neutral"])
        get=neutral["load"]
        def modified(path):
            charged=get(path)
            charged["FLOW"]["low"]["yukawa"]["up"][0][0]=float(yu)
            return charged
        neutral["endpoint"].__globals__["load"]=modified
        wd,ws,wb=row["down_strange_bottom_marginal_widths"]
        pt=neutral["endpoint"]("explicit_open",step,explicit_offsets={**flavor["WEAK"][key],"strange":ws,"bottom":wb})
        assert abs(pt["row"]["gamma_MZ_deg"]-row["gamma_MZ_deg"])<1e-6
        obs=neutral["observables"](pt["C"],pt["b"],neutral["FROZEN"][0])
        out["frozen_neutral"]={k:obs[k] for k in ("angles_in_box","angle_readings","dm31_and_32_in_box_at_central_solar","dm31_eV2","delta_deg","SK_T2K_FC68_90")}
        freq=neutral["FREQ"]
        widths=[bounds(*freq[k]) for k in ("juno_NO_dm21","juno_NO_dm31","accelerator_NO_dm32")]
        seeds=[neutral["SOL"],obs["dm31_eV2"],obs["dm32_eV2"]]
        lo=max(a/z for (a,b),z in zip(widths,seeds))
        hi=min(b/z for (a,b),z in zip(widths,seeds))
        out["frozen_neutral"]["one_solar_scale_can_place_all_three_clocks_in_boxes"]=bool(lo<=hi)
        out["if_recalibrating_two_clocks"]=[]
        for rho in neutral["clocks"](pt["C"],pt["b"],pt["n"]):
            refit=neutral["observables"](pt["C"],pt["b"],rho)
            out["if_recalibrating_two_clocks"].append({k:refit[k] for k in
                                        ("rho","angles_in_box","angle_readings","delta_deg","SK_T2K_FC68_90")})
    return out

def main():
    edge=FROZEN["exploratory_k3"][key]
    raw=observables(UP[0],edge["down_strange_bottom_MZ"])
    # In the physical box, D decreases in S,Q: S^2+1 < 4Q^2.
    assert (SB[1]**2+1)<4*QB[0]**2
    dv=[masses(SB[1],QB[1])[0],masses(SB[0],QB[0])[0]]
    lattice_d=[2/(SS[1]*(1+RR[1])),2/(SS[0]*(1+RR[0]))]
    def rate(row):
        u=observables(UP[0],row["down_strange_bottom_MZ"])
        h=hadrons(u)
        bs=[u[2]/BS[1],u[2]/BS[0]]
        return {"k":row["k"],"Q_if_BMW_S_central":math.sqrt((LIN[0]**2-1)/(4*(h["d_over_s"]*LIN[0]-1))),
                "marginal_MZ_over_QCD_bs_required":bs,"min_relative_bs_matching_budget":threshold(u,bsbox=BS),
                "hadron":h}
    byk={k:[rate(r) for r in entries] for k,entries in FROZEN["solutions"].items()}
    flavor=api(CFG["model"])
    byk[key]=[rate(flavor["solve"](k,key,CFG["step"])) for k in flavor["CONFIG"]["path_lengths"]]
    np.testing.assert_allclose(byk[key][3]["hadron"]["d_over_s"],raw[0]/MATCH[0],rtol=1e-9)
    registered=[weak_case(UP[0]+j*UP[1],CFG["step"],True) for j in (-1,0,1)]
    target=CFG["MILC_S"][0]
    # Substituting a precise strange/average lattice input for the original up
    # input is a separate solution. Keep both threshold coefficients fixed.
    lo,hi=UP[0]-5*UP[1],UP[0]
    l,h=weak_case(lo,CFG["step"],False),registered[1]
    assert l["hadron"]["S"]>target>h["hadron"]["S"]
    for _ in range(8):
        u=lo+(hi-lo)*(l["hadron"]["S"]-target)/(l["hadron"]["S"]-h["hadron"]["S"])
        w=weak_case(u,CFG["step"],False)
        if abs(w["hadron"]["S"]-target)<1e-8:
            break
        if w["hadron"]["S"]>target:
            lo,l=u,w
        else:
            hi,h=u,w
    else:
        raise AssertionError("isospin calibration did not converge")
    extra=weak_case(u,CFG["step"],True)
    rep=weak_case(UP[0],CFG["step"]/2,True)
    err=abs(rep["frozen_neutral"]["dm31_eV2"]-registered[1]["frozen_neutral"]["dm31_eV2"])
    assert err<2e-9
    return {"MZ_over_QCD_ratio_ds_ud_bs_reference":MATCH,"BMW_S_interval":SB,"MILC_S_interval":SS,
            "MILC_mu_md_interval":RR,"eta_Q_interval":QB,"BMW_eta_down_over_s_eliminating_up":dv,
            "MILC_down_over_s_eliminating_up":lattice_d,"central_projected_eta_Q":rate(edge)["Q_if_BMW_S_central"],
            "lengths":byk,"primary_weak_up_input_sensitivity":registered,
            "primary_fixed_MILC_S_instead_of_up":extra,"halfstep_delta31_difference_eV2":err,
            "threshold_examples_same_point":{str(e):feasible(raw,e,Sbox=SS,rbox=RR,Qbox=QB,bsbox=BS)
                                   for e in CFG["match_test_half_widths"]}}

def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)): return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.11g}")
    if isinstance(x,np.bool_): return bool(x)
    return x

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    result=tidy(main())
    path=HERE/"results.json"
    if check:
        assert result==json.loads(path.read_bytes())
        print("1734 shared charged and neutral transport, hadronic ratio boxes and recalibration verified")
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
