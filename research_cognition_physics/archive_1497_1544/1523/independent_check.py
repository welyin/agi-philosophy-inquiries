"""Independent comparison integral, delayed rescue cost, and gauge thresholds."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent
K=8*math.pi**2
C=np.array([13/15,3.,16/3])
B=np.array([33/5,1.,-3.])
CAP=4*math.pi
YS=.8464
TAB=np.array([.468064,.63746,1.0489])
YT2=2*YS**2


def segment(t,g2,b):
    u=np.asarray(t)[...,None]
    w=1-u*g2*b/K
    assert np.all(w>0),"gauge coupling pole"
    val=np.zeros_like(w)
    for i,bb in enumerate(b):
        val[...,i] = -C[i]*g2[i]*u[...,0]/K if bb==0 else C[i]*np.log(w[...,i])/bb
    return np.exp(val.sum(axis=-1))


def envelope(t,g,mu=1e3,n=0):
    cut=max(0.,math.log(mu/1e3))
    left=np.minimum(t,cut)
    right=np.maximum(0.,t-cut)
    gc=g*g/(1-g*g*B*min(cut,float(np.max(t)))/K)
    return segment(left,g*g,B)*segment(right,gc,B+n)


def integral(t,g,mu=1e3,n=0):
    # The exact kernel is positive decreasing before the gauge pole.
    grid=np.linspace(0.,t,16385)
    e=envelope(grid,g,mu,n)
    assert np.all(np.diff(e)<=1e-14)
    step=t/(len(grid)-1)
    lo=6*step*float(np.sum(e[1:]))/K
    hi=6*step*float(np.sum(e[:-1]))/K
    return lo,hi,float(e[-1])


def requirement(t, g, mu=1e3,n=0):
    lo,hi,e=integral(t,g,mu,n)
    return np.array([lo+e/CAP,hi+e/CAP])*YT2


def full_third_family(t,g,mu,n):
    # Third family only, aligned top/bottom. Gauge running across the step is exact.
    s=np.array([YT2,2*.01372**2,2*.010119**2])
    peak=s[0]
    x0=max(0.,math.log(mu/1e3))
    def deriv(time,state):
        gauge=g*g/(1-g*g*(B*time+n*max(0.,time-x0))/K)
        top,bot,tau=state
        a1,a2,a3=gauge
        return state*np.array([6*top+bot-C@gauge,
                               6*bot+top+tau-7*a1/15-3*a2-16*a3/3,
                               4*tau+3*bot-9*a1/5-3*a2])/K
    for start,end in ((0.,min(x0,t)),(min(x0,t),t)):
        if end<=start:
            continue
        dt=(end-start)/1024
        for i in range(1024):
            x=start+i*dt
            a=deriv(x,s)
            b=deriv(x+dt/2,s+dt*a/2)
            c=deriv(x+dt/2,s+dt*b/2)
            d=deriv(x+dt,s+dt*c)
            s+=dt*(a+2*b+2*c+d)/6
            peak=max(peak,s[0])
    return float(peak)


def audit():
    read=lambda p:json.loads((BASE/p).read_text(encoding="utf8"))
    old=read("../1518/results.json")
    data=read("../1369/inputs.json")
    t=math.log(old["witness"]["even_doublets_GeV"][1]/1e3)
    g0=np.array([data[k][0] for k in ["g1_GUT","g2","g3"]])
    inv=1/g0**2-np.array([41/10,-19/6,-7])*math.log(1e3/data["mu_ref_GeV"])/K
    inv-=np.array([0.,2,3])/(48*math.pi**2)
    gs=inv**-.5
    xmass=old["central"]["MX_GeV"]
    alpha_inv=old["witness"]["gauge_inverse_from_each"][0]
    out={}
    for label,g in (("published_DR",TAB),("same_1518_gauge",gs)):
        lo,hi,e=integral(t,g)
        D=(lo+hi)/2
        f0=YT2*(D+e/CAP)
        delays={}
        for mu in [1e3,1e6,1e8,1e9,1e10]:
            tcut=math.log(mu/1e3)
            a,_,ec=integral(tcut,g)
            if YT2*(a+ec/CAP)>=1:
                delays[str(int(mu))]={"already_crosses_cap":True}
            else:
                # Conservative version with the lower Darboux sums.
                delays[str(int(mu))]={"integrated_negative_beta_at_least":
                    max(0.,math.log(YT2*(lo-a+e/CAP)/(1-YT2*a)))}
        scenarios={}
        for mu in [1e3,1e4,1e5]:
            for n in range(1,10):
                budget=requirement(t,g,mu,n)
                if budget[1] < 1:
                    assert requirement(t,g,mu,n-1)[0]>1
                    top_peak=full_third_family(t,g,mu,n)
                    assert top_peak<CAP
                    scenarios[str(int(mu))]={"min_pairs_in_window":n,
                        "necessary_budget_without_last_pair":requirement(t,g,mu,n-1).tolist(),
                        "budget_with_all_pairs":budget.tolist(),"aligned_peak_yt_squared":top_peak,
                        "inverse_alpha5_if_other_masses_fixed":alpha_inv-n*math.log(xmass/mu)/(2*math.pi)
                            if label=="same_1518_gauge" else None}
                    break
        out[label]={"gauge_at_1TeV":g.tolist(),"unmodified_endpoint_budget":requirement(t,g).tolist(),
                    "central_delta_t_needed_at_least":math.sqrt(f0)-1,
                    "delayed_negative_correction":delays,
                    "countermodels_with_complete_multiplets":scenarios}
    return out

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
