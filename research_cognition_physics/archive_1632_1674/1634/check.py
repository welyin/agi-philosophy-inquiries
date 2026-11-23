"""Gauge threshold intersections after adding the bottom-current doublet."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

D=Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def meeting(mu,gY,g2,g3,m1,m2,md):
    i1,i2,i3=3/(5*gY**2),1/g2**2,1/g3**2
    s=math.log(m1/mu)+math.log(m2/mu)
    d=math.log(md/mu)
    t=(8*math.pi**2*(i1-i3)-.4*s+2*d)/12.3
    w=((7.8+4/3)*t-4*s/15-4*d/3-8*math.pi**2*(i1-i2))/(4/3)
    return {"U_GeV":mu*math.exp(t),"equivalent_triplets":w/t,"P_GeV":mu*math.exp(t-w/2)}


def run():
    cfg=read(D/"inputs.json")
    old=read(D/cfg["baseline"])
    reg=read((D/cfg["baseline"]).parent/old["reference"])
    fit=read(D/cfg["Z_fit"])
    mu=np.array(fit["central"])
    sig=np.array(fit["error"])
    cov=np.outer(sig,sig)*np.array(fit["correlation"])
    a=np.array([0,0,-1/3,1])
    b=np.array([-1,0,1,0])
    pmin,pmax=runpy.run_path(str(D/cfg["round1633"]))["project_ratio"](mu,cov,a,b,2.)
    m,v=cfg["mb_GeV"],cfg["vev_GeV"]
    def ceiling(cap):
        F=v*cap/math.sqrt(2)*math.sqrt((1-pmin)/pmin)
        return .5*(F+math.hypot(F,2*m))
    cap=math.sqrt(4*math.pi)
    mhi=ceiling(cap)
    atlas=cfg["conditional_ATLAS_lower_GeV"]
    required=(math.sqrt(2)/v)*(atlas-m*m/atlas)*math.sqrt(pmin/(1-pmin))
    ref=reg["reference"]
    g0=ref["g3"]
    change=reg["matching"]["g3_per_alpha"]*old["widths"]*reg["marginal_width"]["alpha_s"]
    strong=[g0*(1-change),g0*(1+change)]
    m1box,m2=old["two_B_masses_GeV"]
    def point(g3,m1,md):
        return meeting(old["matching_scale_GeV"],ref["gY"],ref["g2"],g3,m1,m2,md)
    scenarios={}
    bridge=read(D/cfg["bridge_results"])["matching"]["scenarios"]["ATLAS_single_muon"]
    oldP=read(D/"../1611/results.json")["whole_triplet_geometric_mean_GeV"]
    width=bridge["windows"][-1]
    assert width["factor"]==1.1
    for name,mdlo in (("formal",200),("conditional_ATLAS",atlas)):
        ends=[point(*corner) for corner in itertools.product(strong,m1box,(mdlo,mhi))]
        bounds={key:[min(pt[key] for pt in ends),max(pt[key] for pt in ends)] for key in ends[0]}
        assert 1<bounds["equivalent_triplets"][0]<=bounds["equivalent_triplets"][1]<2
        plo,phi=bounds["P_GeV"]
        bands={o:[lo*(oldP[1]/phi)**2,hi*(oldP[0]/plo)**2] for o,(lo,hi) in width["Z_invariant"].items()}
        scenarios[name]={**bounds,"Z_invariant_1p10":bands,
                         "gauge_part_neutral_mixing_ceiling_y1":v*v/(2*plo**2)}
    # Keep the identical Z point in the two examples. It is the extremum of D_R.
    z=mu-2*(cov@a)/math.sqrt(a@cov@a)
    p=(a@z)/(b@z)
    dl=np.array([-1,1,2/3,0])@z
    def quark(g,m1,md):
        prob_s=dl/(b@z)-2*(m/md)**2*p
        # Both original weak-singlet masses are hard masses in the symmetric phase.
        # Give each ordinary-left mixing equal squared weight in this witness.
        hard_s=np.array([m1,m2])
        ls=np.sqrt(prob_s/2)*np.ones(2)
        rs=m*ls/hard_s
        ld=m*math.sqrt(p)/md
        rb=math.sqrt(1-p-rs@rs)
        lb=math.sqrt(1-ls@ls-ld*ld)
        xl=(m*rs-hard_s*ls)/lb
        xr=(m*ld-md*math.sqrt(p))/rb
        mb0=(m*lb-xl@rs)/rb
        matrix=np.diag([mb0,m1,m2,md])
        matrix[0,1:3]=xl
        matrix[3,0]=xr
        l=np.r_[lb,ls,ld]
        r=np.r_[rb,rs,math.sqrt(p)]
        assert np.linalg.norm(matrix@r-m*l)<1e-11 and np.linalg.norm(matrix.T@l-m*r)<1e-11
        yuk=math.sqrt(2)*abs(xr)/v
        assert yuk<cap
        bb=(1-prob_s-p)**2
        cms=read(D/"../1602/inputs.json")["CMS_ratio_bb_ZZ_normalized"]
        assert cms[0]-2*cms[1] <= bb <= cms[0]+2*cms[2]
        return {**point(g,m1,md),"bb_ZZ_over_SM":bb,"right_doublet_p":p,"charged_exotic_GeV":md,"physical_bottomlike_GeV":np.sort(np.linalg.svd(matrix,compute_uv=False)).tolist(),"doublet_Y":yuk}
    light=read(D/cfg["oscillation"])["box"]
    alo,ahi=light["abs31_eV2"]
    rlo,rhi=light["solar_over_abs31"]
    no=quark(strong[1],m1box[0],1701.)
    io=quark(strong[0],m1box[1],5460.)
    mno=np.sqrt([ahi*rhi,ahi])*1e-9
    mio=np.sqrt([alo,alo*(1+rlo)])*1e-9
    def invariant(mass,gate,factor):
        # Exact 2x2 Majorana pairs; light masses are input Takagi eigenvalues.
        return float(np.prod(mass/(factor*gate+mass)))
    fmin,fmax=1/1.1,1.1
    qno=invariant(mno,no["P_GeV"],fmin)
    low,high=fmin,fmax
    assert invariant(mio,io["P_GeV"],low)>qno>invariant(mio,io["P_GeV"],high)
    for _ in range(64):
        mid=(low+high)/2
        if invariant(mio,io["P_GeV"],mid)>qno: low=mid
        else: high=mid
    qio=invariant(mio,io["P_GeV"],high)
    assert abs(qio/qno-1)<1e-13
    return {"right_probability_range":[pmin,pmax],"MD_ceiling_GeV":{"y1":ceiling(1),"ysqrt4pi":mhi},
            "necessary_Y_at_1700":required,"scenarios":scenarios,
            "equal_q_witness":{"NO":{**no,"gate_to_bare":fmin,"qZ":qno},"IO":{**io,"gate_to_bare":high,"qZ":qio}}}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=run()
    target=D/"results.json"
    if args.check:
        assert result==read(target)
        print("1634 PASS: doublet bound, eight gate corners per window and equal-q mass realizations")
    else:
        with target.open("w",encoding="utf-8") as f:
            f.write(json.dumps(result,indent=2)+"\n")
