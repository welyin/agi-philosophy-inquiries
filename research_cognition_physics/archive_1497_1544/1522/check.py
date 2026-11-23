"""Same scale Higgs quartic and gauge threshold for complete scalar tens."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = 16*math.pi**2
GROUP = np.array([4.1,-19/6,-7])

def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))

def central_state(rge, scale, n):
    x=np.array([rge.MASS[key][0]**2 for key in rge.KEYS]+
               [rge.HIGGS["lambda_paper"][0]/4])
    step=math.log(scale/rge.MZ)/n
    for i in range(n):
        t=i*step
        a=rge.beta(t,x)
        b=rge.beta(t+step/2,x+a*step/2)
        c=rge.beta(t+step/2,x+b*step/2)
        d=rge.beta(t+step,x+c*step)
        x+=step*(a+2*b+2*c+d)/6
    return x

def enclose(cfg, scale, n=512):
    # Natural interval Euler enclosure, with a verified speed bound across each entire step.
    dat=read(cfg["gauge_yukawa"])
    src=read(cfg["quartic"])
    half=cfg["symmetric_low_energy_input_widths"]
    def borders(keys):
        p=np.array([dat[key] for key in keys])
        return (p[:,0]-half*p[:,1])**2,(p[:,0]+half*p[:,1])**2
    gl,gh=borders(["g1_GUT","g2","g3"])
    yl,yh=borders(["yt","yb","ytau"])
    l,err=src["lambda_paper"]
    low=np.r_[yl,(l-half*err)/4]
    high=np.r_[yh,(l+half*err)/4]
    dt=math.log(scale/dat["mu_ref_GeV"])/n
    speed=np.array([.30,.00025,.00003,.09])
    gain=np.array([[4.5,1.5,1],[1.5,4.5,1],[3,3,2.5]])
    loss=np.array([[17/12,9/4,8],[5/12,9/4,8],[15/4,9/4,0]])
    for i in range(n):
        left=low-speed*dt-1e-13
        right=high+speed*dt+1e-13
        assert left.min()>0
        gg=np.array([g/(1-2*GROUP*g*t/LOOP) for g in (gl,gh) for t in (i*dt,(i+1)*dt)])
        g0=gg.min(axis=0)-1e-13
        g1=gg.max(axis=0)+1e-13
        g0[0]*=.6
        g1[0]*=.6
        fl=gain@left[:3]-loss@g1
        fh=gain@right[:3]-loss@g0
        possible=np.array([2*y*f for y in (left[:3],right[:3]) for f in (fl,fh)])
        ll,hh=left[-1],right[-1]
        trace=np.array([3.,3.,1.])
        def gauge4(g):
            return 3/8*(2*g[1]**2+(g[0]+g[1])**2)
        scalar_min=(24*ll**2+4*ll*(trace@left[:3])-(9*g1[1]+3*g1[0])*hh
                    -2*(trace@right[:3]**2)+gauge4(g0))
        scalar_max=(24*hh**2+4*hh*(trace@right[:3])-(9*g0[1]+3*g0[0])*ll
                    -2*(trace@left[:3]**2)+gauge4(g1))
        dlo=np.r_[possible.min(axis=0),scalar_min]/LOOP
        dhi=np.r_[possible.max(axis=0),scalar_max]/LOOP
        assert np.all(np.maximum(abs(dlo),abs(dhi))+1e-12<speed)
        low+=dt*dlo-1e-12
        high+=dt*dhi+1e-12
    return low,high

def stop_response(coef, ratio, x2):
    return coef*(2*math.log(ratio)+x2*(1-x2/12))

def hgg_and_masses(mt, mass, x2):
    z=(mt/mass)**2
    m2=mass**2*np.array([1+z-math.sqrt(x2*z),1+z+math.sqrt(x2*z)])
    assert m2.min()>0
    amplitude=z*(2+2*z-x2)/(4*((1+z)**2-z*x2))
    return dict(physical_tree_stop_masses_GeV=np.sqrt(m2).tolist(),
                relative_to_heavy_top_amplitude=amplitude)

def calculate():
    cfg=read("inputs.json")
    scale=read(cfg["gauge_witness"])["witness"]["soft_GeV"]
    assert read(cfg["same_angle"])["tan_beta"]==1
    spec=importlib.util.spec_from_file_location("rge1494",HERE/cfg["previous_RGE"])
    rge=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rge)
    x=central_state(rge,scale,256)
    np.testing.assert_allclose(x,central_state(rge,scale,512),atol=1e-13,rtol=1e-12)
    low,high=enclose(cfg,scale)
    assert np.all((low<x)&(x<high))
    yt=math.sqrt(x[0])
    coeff=3*x[0]**2/LOOP
    coeff_upper=3*high[0]**2/LOOP
    mt=yt*cfg["vev_GF_GeV"]/math.sqrt(2)
    out={}
    for s in cfg["complete_ten_mass_ratio_range"]:
        ten=np.array([scale/math.sqrt(s),scale/math.sqrt(s),scale*s])
        # Per complex scalar, indices of Q, u^c, e^c divided by three.
        beta=np.array([[.1,1.5,1],[.8,0,.5],[.6,0,0]])/3
        np.testing.assert_allclose(beta.sum(axis=0),.5)
        delta=np.log(ten/scale).sum()*beta.sum(axis=0)/(2*math.pi)
        assert max(abs(delta))<1e-14
        cap=stop_response(coeff,s,6)
        out[str(s)]=dict(tens_GeV=ten.tolist(),delta_inverse_gauge=delta.tolist(),
            stop_positive_max=cap,missing_central=x[3]-cap,
            necessary_other_threshold_whole_input_box=low[3]-stop_response(coeff_upper,s,6),
            no_mixing=hgg_and_masses(mt,ten[-1],0),
            maximal_quartic=hgg_and_masses(mt,ten[-1],6))
    z=(mt/scale)**2
    x2_blind=2+2*z
    null=hgg_and_masses(mt,scale,x2_blind)
    same=12-x2_blind
    np.testing.assert_allclose(stop_response(coeff,1,same),stop_response(coeff,1,x2_blind))
    comparison=hgg_and_masses(mt,scale,same)
    assert abs(null["relative_to_heavy_top_amplitude"])<1e-15
    return dict(scale_GeV=scale,yt=yt,quartic=x[3],gauge_fixed_coupling=read(cfg["gauge_witness"])["witness"]["gauge_unified"],
        coefficient=coeff,two_width_ODE_enclosure=dict(yt=[math.sqrt(low[0]),math.sqrt(high[0])],quartic=[low[3],high[3]]),
        published_table_1TeV_yt=cfg["TeV_table_yt"],bounded_cases=out,
        first_s_formal_central=math.exp((x[3]/coeff-3)/2),
        first_s_formal_necessary_at_input_box=math.exp((low[3]/coeff_upper-3)/2),
        equal_quartic_distinct_gluons=dict(x_squared=[x2_blind,same],delta_lambda=stop_response(coeff,1,same),
            amplitude_factors=[null["relative_to_heavy_top_amplitude"],comparison["relative_to_heavy_top_amplitude"]]))

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            compare(a[key],b[key])
    else:
        np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-13)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--write-new",action="store_true")
    g.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=calculate()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as h:
            json.dump(result,h,indent=2)
            h.write("\n")
    else:
        compare(result,read("results.json"))
        print("1522 OK: all three shared thresholds, quartic and gluon response verified")
