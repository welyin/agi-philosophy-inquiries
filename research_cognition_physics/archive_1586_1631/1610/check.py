"""Round 1610: the prompt singlet extension; frozen results are read-only on --check."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = 16*math.pi**2


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def rotation(s):
    c = math.sqrt(1-s*s)
    v = np.array([1,1])/math.sqrt(2)
    return np.block([[np.array([[c]]), (s*v)[None,:]],
                     [(-s*v)[:,None], np.eye(2)-s*s/(1+c)*np.outer(v,v)]])


def vector_factor(M, mf, mv):
    x, y = (mf/M)**2, (mv/M)**2
    return math.sqrt((1-x-y)**2-4*x*y)*((1-x)**2+(1+x)*y-2*y*y)


def higgs_factor(M, mf, mh):
    x, y = (mf/M)**2, (mh/M)**2
    return math.sqrt((1-x-y)**2-4*x*y)*((1+x)*(1+x-y)+4*x)


def widths(masses,s,cfg,ckm):
    v = cfg["v_GeV"]
    sm = cfg["kinematic_GeV"]
    row = rotation(s)[0]
    decays = []
    for i,M in enumerate(masses,1):
        pref = row[i]**2*M**3/(32*math.pi*v*v)
        channels={}
        for quark,weight in zip(("u","c","t"),ckm):
            channels[quark+"W"]=2*pref*weight*vector_factor(M,sm[quark],sm["W"])
        for n,partner in enumerate([sm["b"],*masses[:i-1]]):
            label = "b" if n==0 else "B1"
            channels[label+"Z"]=pref*row[n]**2*vector_factor(M,partner,sm["Z"])
            channels[label+"h"]=pref*row[n]**2*higgs_factor(M,partner,sm["H"])
        total=sum(channels.values())
        probs={key:x/total for key,x in channels.items()}
        retained=min(probs["tW"]/.5,probs["bZ"]/.25,probs["bh"]/.25)**2
        decays.append({"mass_GeV":M,"width_GeV":total,"proper_length_mm":cfg["hbarc_GeV_mm"]/total,
                       "gamma_over_mass":total/M,"branching":probs,
                       "every_nominal_pair_channel_weight_ratio_floor":retained})
    return decays


def gates_and_yukawa(masses,s,cfg):
    v=cfg["v_GeV"]
    mb=cfg["low_yb_matching"]*v/math.sqrt(2)
    block=rotation(s)@np.diag([mb,*masses])
    _, sigmas, vdag=np.linalg.svd(block[1:,:],full_matrices=True)
    y=math.sqrt(2)/v*block[0,:]@vdag.T
    assert 0 < sigmas[1] < sigmas[0]
    return sorted(zip(sigmas.tolist(),(y[:2]**2).tolist())),abs(float(y[2]))


def shoot(yt,log_end,heavy,yb,overlap,ref,steps):
    gates=[(math.log(m/200),y) for m,y in heavy]
    g0=np.array([ref["gY"],ref["g2"],ref["g3"]])
    beta=np.array([41/6,-19/6,-7.])
    addition=np.array([4/9,0.,2/3])
    def rhs(t,state):
        u,d,p,l=state
        a,b,c=1/(1/g0**2-2*(beta*t+addition*sum(max(0,t-tg) for tg,_ in gates))/LOOP)
        return np.array([u*(9*u+(6-3*p)*d-17*a/6-9*b/2-16*c),
                         d*(9*d+(6-3*p)*u-5*a/6-9*b/2-16*c),
                         -3*p*(1-p)*(u+d),
                         24*l*l+(12*(u+d)-9*b-3*a)*l+3*(2*b*b+(a+b)**2)/8-6*(u*u+d*d)])/LOOP
    state=np.array([yt**2,yb**2,overlap,ref["lambda"]])
    start=0.
    for end,jump in [*gates,(log_end,0)]:
        count=max(2,math.ceil((end-start)*steps/log_end))
        h=(end-start)/count
        for j in range(count):
            t=start+j*h
            a=rhs(t,state); b=rhs(t+h/2,state+h*a/2)
            c=rhs(t+h/2,state+h*b/2); d=rhs(t+h,state+h*c)
            state+=h*(a+2*b+2*c+d)/6
        state[1]+=jump
        start=end
    slope=rhs(log_end,state)
    delta=1.e-4
    curvature=(rhs(log_end+delta,state+delta*slope)[3]-rhs(log_end-delta,state-delta*slope)[3])/(2*delta)
    return np.array([state[3],slope[3]*LOOP]),float(curvature)


def critical(heavy,yb,overlap,ref,steps):
    z=np.array([.91,math.log(1.e16/200)])
    for _ in range(10):
        def fun(v):
            return shoot(v[0],v[1],heavy,yb,overlap,ref,steps)[0]
        residual=fun(z)
        if max(abs(residual))<5.e-11:
            break
        plus=np.array([[2.e-6,0],[0,1.e-3]])
        jac=np.column_stack([(fun(z+dv)-fun(z-dv))/(2*max(dv)) for dv in plus])
        z-=np.linalg.solve(jac,residual)
    residual,curvature=shoot(*z,heavy,yb,overlap,ref,steps)
    assert max(abs(residual))<5.e-10 and curvature>0
    return {"yt_200_GeV":float(z[0]),"high_scale_GeV":200*math.exp(z[1]),
            "positive_local_curvature":curvature}


def evaluate():
    cfg=read("inputs.json")
    old=read(cfg["baseline"])["two_weak_singlet_control"]["gauge_basis_five_masses_GeV"]
    assert cfg["heavy_examples_GeV"][0]==old[-2:]
    weak=read(cfg["angles_masses"])
    _,s13,s23,_=[a[0] for a in weak["quark_sines_and_phase"]]
    ckm=[s13**2,(1-s13**2)*s23**2,(1-s13**2)*(1-s23**2)]
    assert abs(sum(ckm)-1)<2e-15
    smin,smax=cfg["mixing_s_window"]
    required=cfg["hbarc_GeV_mm"]*cfg["beta_gamma_cap"]*math.log(1/(1-cfg["probability_decay_required"]))/cfg["lab_radius_mm"]
    examples={}
    ref=read(cfg["rg"])["reference"]
    pure=critical([],cfg["low_yb_matching"],ckm[2],ref,176)
    top_only=critical([],0,ckm[2],ref,176)
    assert abs(top_only["yt_200_GeV"]-.89939)<1.e-5
    for masses in cfg["heavy_examples_GeV"]:
        dec=widths(masses,smin,cfg,ckm)
        # Uniform analytic bound using only the tW width for both poles.
        floor=min(m**3*ckm[2]*vector_factor(m,cfg["kinematic_GeV"]["t"],cfg["kinematic_GeV"]["W"])
                  /(32*math.pi*cfg["v_GeV"]**2) for m in masses)
        s_prompt=math.sqrt(required/floor)
        assert s_prompt<smin and smax**2<=cfg["light_relative_probability_budget"]
        gates,yl=gates_and_yukawa(masses,smax,cfg)
        r=critical(gates,yl,ckm[2],ref,176)
        cross=critical(gates,yl,ckm[2],ref,88)
        assert abs(r["yt_200_GeV"]-cross["yt_200_GeV"])<3e-8
        assert abs(r["high_scale_GeV"]/cross["high_scale_GeV"]-1)<2e-5
        aligned=critical([(m,0.) for m in masses],cfg["low_yb_matching"],ckm[2],ref,176)
        examples[str(masses[0])]={
            "poles_at_s_min":dec,"s_sufficient_prompt_all_poles":s_prompt,
            "required_width_GeV":required,"gamma_over_mass_max_at_smax":
                max(p["gamma_over_mass"] for p in widths(masses,smax,cfg,ckm)),
            "max_heavy_hard_threshold_relative_change":max(abs(mh/mp-1) for (mh,_),mp in zip(gates,masses)),
            "high_gate_yukawa_sq":sum(y for _,y in gates),
            "ratio_light_W_shift_at_smax":smax*smax/(1+math.sqrt(1-smax*smax)),
            "rg_mixed_at_smax":r,"rg_aligned_same_gates":aligned,
            "mu_root_shift_from_no_extra":r["yt_200_GeV"]-pure["yt_200_GeV"]}
    return {"SM_top_only_cross_check":top_only,"SM_same_top_bottom_contract":pure,"examples":examples}


def equal(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=2e-7,abs_tol=2e-15)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=evaluate()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as fp:
            json.dump(out,fp,ensure_ascii=False,indent=2,allow_nan=False)
            fp.write("\n")
    elif not equal(out,read("results.json")):
        raise SystemExit("1610 reproduction failed")
    print("1610 finite rates, prompt criterion and matched two-Yukawa running verified")
