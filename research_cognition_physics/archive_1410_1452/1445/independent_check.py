"""A real tridiagonal reconstruction and a tagged neutral-B time asymmetry."""
import argparse
import itertools
import json
from math import asin, cos, degrees, hypot, pi, radians, sin, sqrt
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def state(ell,b,a,c):
    x=abs(ell)
    m=np.array([[c*x**4,b*x**3,0.],[b*x**3,b*x*x,a*x*x],[0.,a*x*x,1.]])
    ev,vec=np.linalg.eigh(m)
    order=np.argsort(abs(ev))
    return ev[order],vec[:,order]


def quantities(pars,orientation=1):
    scale,b,a,cu,cd=pars
    uv,uo=state(scale*sin(pi/8),b,a,cu)
    dv,do=state(scale*cos(pi/8),b,a,cd)
    v=uo.T@np.diag([orientation*1j,-1,1])@do
    calibrated=np.array([uv[0]/uv[1],dv[0]/dv[1],dv[1]/dv[2],abs(v[0,1]),abs(v[1,2])])
    gamma=np.angle(-v[0,0]*v[0,2].conjugate()/(v[1,0]*v[1,2].conjugate()),deg=True)
    beta=np.angle(-v[1,0]*v[1,2].conjugate()/(v[2,0]*v[2,2].conjugate()),deg=True)
    J=float((v[0,0]*v[1,1]*v[1,0].conjugate()*v[0,1].conjugate()).imag)
    return calibrated,float(gamma%180),float(beta),J


def solve(target):
    b=np.sign(target[2])*.43
    scale=.225
    cu=b*(1+target[0]/(scale*sin(pi/8))**2)
    cd=b*(1+target[1]/(scale*cos(pi/8))**2)
    pars=np.array([scale,b,.85,cu,cd])
    def residual(p):
        return (quantities(p)[0]-target)/abs(target)
    for _ in range(28):
        r=residual(pars)
        if max(abs(r))<1e-11:
            return pars
        basis=np.eye(5)*1e-6
        jac=np.column_stack([(residual(pars+h)-residual(pars-h))/(2e-6) for h in basis])
        step=np.linalg.solve(jac,-r)
        for n in range(18):
            trial=pars+step/2**n
            if np.linalg.norm(residual(trial))<np.linalg.norm(r):
                pars=trial
                break
        else:
            raise ValueError("Local Newton iteration failed")
    raise ValueError("No local convergence")


def minimum_phase(phase,lower,upper):
    intervals=[(lower,upper),(180-upper,180-lower)]
    return min(max(lo+360*j-phase, phase-hi-360*j, 0.)
               for lo,hi in intervals for j in (-1,0,1))


def run():
    data=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    m=data["mass_ratios"]
    base=np.array([m["mu_mc"][0],m["md_ms"][0],m["ms_mb"][0],
                   data["K_product"][0]/data["K_form"][0],
                   data["BelleII_Vcb"][0]*np.divide(*data["run_cb_ratio"])])
    err_s=hypot(.013,.008)
    cmax=.008+2*hypot(.012,.003)
    low=degrees(asin(.717-2*err_s))
    high=degrees(asin((.717+2*err_s)/sqrt(1-cmax*cmax)))
    rows=[]
    for middle,ua,da in itertools.product((-1,1),repeat=3):
        target=base*np.array([ua,da,middle,1,1])
        pars=solve(target)
        primary,gamma,beta,j=quantities(pars)
        other,mirror,beta2,j2=quantities(pars,-1)
        assert np.allclose(primary,target,atol=1e-12)
        assert np.allclose(primary,other,atol=1e-14)
        assert np.isclose(gamma+mirror,180) and np.isclose(beta,-beta2) and np.isclose(j,-j2)
        rows.append({"signs":[middle,ua,da],"gamma":gamma,"mirror_gamma":mirror,
                     "beta":beta,"S_if_no_extra_phase":sin(radians(2*beta)),
                     "phase_offset_floor_deg":minimum_phase(2*beta,low,high),
                     "mirror_phase_offset_floor_deg":minimum_phase(2*beta2,low,high)})
    return {"two_width_allowed_low_phase_branch_degrees":[low,high],
            "central_eight_independent_branches":rows,
            "minimal_mirror_extra_phase_across_eight_deg":min(r["mirror_phase_offset_floor_deg"] for r in rows)}


def same(x,y):
    if isinstance(x,dict):
        return x.keys()==y.keys() and all(same(x[k],y[k]) for k in x)
    if isinstance(x,list):
        return len(x)==len(y) and all(same(a,b) for a,b in zip(x,y))
    return np.isclose(x,y,rtol=1e-9,atol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=run()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    else:
        saved=json.loads((HERE/"independent_results.json").read_text(encoding="utf-8"))
        assert same(out,saved)
        print("1445 eight real-tridiagonal roots, signed gamma branches and neutral-B time test verified")
