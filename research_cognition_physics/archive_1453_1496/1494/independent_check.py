"""Run the measured third-family data independently and resolve all heavy charged modes."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


BASE=Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8-sig"))


def run(scale,steps):
    inp=read("../1369/inputs.json")
    lam=read("../1396/inputs.json")["lambda_paper"][0]/4
    initial_g=np.array([math.sqrt(.6)*inp["g1_GUT"][0],inp["g2"][0],inp["g3"][0]])
    initial_b=np.array([41/6,-19/6,-7.])
    def gauge(t):
        return 1/(1/initial_g**2-initial_b*t/(8*math.pi**2))
    z=np.array([math.log(inp[k][0]) for k in ["yt","yb","ytau"]]+[lam])
    def derivative(t,z):
        yy,w,c=gauge(t)
        u,d,l=np.exp(2*z[:3])
        lam=z[3]
        ydot=np.array([4.5*u+1.5*d+l-17*yy/12-9*w/4-8*c,
                       4.5*d+1.5*u+l-5*yy/12-9*w/4-8*c,
                       2.5*l+3*d+3*u-15*yy/4-9*w/4])
        quartic=(24*lam**2+lam*(12*u+12*d+4*l-9*w-3*yy)-
                 6*u**2-6*d**2-2*l**2+3*(2*w**2+(w+yy)**2)/8)
        return np.append(ydot,quartic)/(16*math.pi**2)
    last=math.log(scale/inp["mu_ref_GeV"])
    step=last/steps
    for i in range(steps):
        t=i*step
        a=derivative(t,z)
        b=derivative(t+step/2,z+a*step/2)
        c=derivative(t+step/2,z+b*step/2)
        d=derivative(t+step,z+c*step)
        z+=step*(a+2*b+2*c+d)/6
    return np.r_[np.exp(z[:3]),z[3],gauge(last)[:2]]


def audit():
    scale=3e6
    data=run(scale,512)
    np.testing.assert_allclose(data,run(scale,1024),rtol=1.e-10)
    yt,yb,yl,lam,gy2,g22=data
    pphys=yt**2/(yt**2+yb**2)
    D=(gy2+g22)/8
    branch={}
    for s in [-1,1]:
        p=(1+s*math.sqrt(lam/D))/2
        ymax=max(yt/math.sqrt(p),yb/math.sqrt(1-p))
        rem=1-(yt**2+yb**2)/ymax**2
        # Label active up, spectator up, active conjugate down, spectator conjugate down.
        light=np.array([yt/ymax,0.,yb/ymax,math.sqrt(rem)])
        signs=np.diag([1.,1.,-1.,-1.])
        assert abs(np.linalg.norm(light)-1)<1.e-12
        assert abs(D*(light@signs@light)**2-lam)<1e-12
        # A real Householder basis; the other three columns are physical heavy charged modes.
        e=np.eye(4)
        delta=e[:,0]-light
        modes=e-2*np.outer(delta,delta)/(delta@delta)
        np.testing.assert_allclose(modes[:,0],light,atol=2.e-15)
        mass2=scale**2*(e-np.outer(light,light))
        np.testing.assert_allclose(mass2@modes[:,1:]/scale**2,modes[:,1:],atol=2.e-15)
        gt,gb=ymax*modes[0,1:],ymax*modes[2,1:]
        heavy=float(gt@gt+gb@gb)
        assert abs(heavy-(2*ymax**2-yt**2-yb**2))<1e-13
        # Conservative tree-level kinematic lower bound at M>=3 PeV, mt+mb<=180 GeV.
        correction=(1-(180/scale)**2)**2
        branch[str(s)]={"p_up_total":p,"tan_beta":math.sqrt(p/(1-p)),
                        "Yukawa_equal_min":ymax,"inactive_doublet_weight_min":rem,
                        "heavy_bR_tL_norm_squared":float(gb@gb),"heavy_bL_tR_norm_squared":float(gt@gt),
                        "sum_tb_partial_width_over_mass_floor":3/(16*math.pi)*correction*heavy}
    single_D=D*(2*pphys-1)**2
    return {"scale_GeV":scale,"third_family_yukawa":data[:3].tolist(),
            "required_quartic":lam,"quartic_unit_D":D,"single_pair_p_from_equal_Y":pphys,
            "single_pair_quartic_minus_running":single_D-lam,"branches":branch}


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
