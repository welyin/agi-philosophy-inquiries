"""Exact below-threshold one-loop kernels and a complete charged-family cancellation."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent


def fermion(pole,heavy):
    tau=(pole/(2*heavy))**2
    if not 0<tau<1:
        raise ValueError("below pair threshold required")
    return 2*(tau+(tau-1)*math.asin(math.sqrt(tau))**2)/tau**2


def direct_integral(pole,heavy):
    roots,weights=np.polynomial.legendre.leggauss(36)
    x=(roots[:,None]+1)/2
    y=(roots[None,:]+1)*(1-x)/2
    xy=4*x*y
    tau=(pole/(2*heavy))**2
    return float(np.sum((weights[:,None]*weights[None,:])*(1-x)*(1-xy)/(1-tau*xy)))


def family(h,u,d,l,s):
    A,B,C=(fermion(h,i) for i in (u,d,l))
    for base,m in zip((A,B,C),(u,d,l)):
        assert abs(base-direct_integral(h,m))<1e-11
    x=(-s*(B-A)+math.hypot(s*(B-A),2*math.sqrt(A*B)))/2
    ku,kd=x/A,-x/B
    assert math.isclose(ku*kd-s*(ku+kd),-1.)
    cbma=math.sqrt(1-s*s)
    tanb=cbma/(ku-s)
    assert math.isclose(s-cbma*tanb,kd)
    beta=math.atan(tanb)
    glu=ku*A+kd*B
    photon=4*ku*A/3+kd*B/3
    total=photon+kd*C
    assert abs(glu)<1e-12 and abs(photon-x)<1e-12
    lo,hi=sorted([A,B])
    assert lo-1e-13<=abs(x)<=hi+1e-13
    assert ku**2+kd**2>=1+(lo/hi)**2-1e-12
    return {"masses_u_d_lepton_GeV":[u,d,l],"electroweak_projection":s,"tan_beta":tanb,
            "kappa_up":ku,"kappa_down":kd,"yukawa_up_down_lepton":
                [math.sqrt(2)*m/(246.22*vev) for m,vev in [(u,math.sin(beta)),(d,math.cos(beta)),(l,math.cos(beta))]],
            "extra_on_shell_gg":float(glu),"extra_quarks_on_shell_photon":photon,
            "extra_full_on_shell_photon":total,
            "soft_delta_c_gghh":-(ku*ku+kd*kd)/3}


def audit():
    same=family(125.11,180.,130.,130.,.6)
    different=family(125.11,180.,130.,240.,.6)
    assert abs(same["extra_full_on_shell_photon"])<1e-12
    assert abs(different["extra_full_on_shell_photon"])>1e-2
    for s in np.linspace(-.99,.99,7):
        for hi,lo in [(1.999,4/3+.0001),(1.7,1.33+0.01)]:
            for A,B in [(hi,lo),(lo,hi)]:
                for sign in (-1,1):
                    x=(-s*(B-A)+sign*math.hypot(s*(B-A),2*math.sqrt(B*A)))/2
                    assert min(A,B)<=abs(x)<=max(A,B)
                    assert 1/A**2*x**2+1/B**2*x**2>=13/9
    return {"gamma_quarks_lower_per_sequential_family_in_fermion_kernel_units":4/3,
            "Q_lower_per_family_when_physical_ggh_invisible":13/9,
            "abs_soft_delta_c_gghh_floor_per_family":13/27,
            "complete_charged_family_example":same,"lepton_mass_changed_only":different}


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
