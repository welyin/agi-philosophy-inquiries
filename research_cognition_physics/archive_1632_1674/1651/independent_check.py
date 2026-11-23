"""Pole, static response and coherent h/H/A exchange from one CP-even 2HDM tree Hessian."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent
V=246.22
MH=125.1


def potential(M,c):
    s=math.sqrt(1-c*c)
    delta=(M**2-MH**2)/V**2
    z1=(s*s*MH**2+c*c*M**2)/V**2
    z6=-s*c*delta
    z5=-c*c*delta
    Y2=M**2-V**2
    hessian=np.array([[z1*V**2,z6*V**2],[z6*V**2,M**2+z5*V**2]])
    m2,Q=np.linalg.eigh(hessian)
    np.testing.assert_allclose(m2,[MH**2,M**2],rtol=3e-14)
    if Q[0,0]<0: Q[:,0]*=-1
    if Q[0,1]<0: Q[:,1]*=-1
    np.testing.assert_allclose(Q,[[s,c],[c,-s]],atol=1e-12)
    slope_static=-hessian[0,1]/hessian[1,1]
    slope_pole=-hessian[0,1]/(hessian[1,1]-MH**2)
    assert abs(slope_pole-c/s)<1e-13
    # Exact complete-square sufficient condition. Negative coefficient need not mean instability.
    square_certificate=2+2*(z5-z6*z6/z1)
    t=.3
    quartic_ray=z1/2-2*abs(z6)*t+(2+2*z5)*t*t+1.5*t**4
    return {"s":s,"c":c,"Z1":z1,"Z6":z6,"Z5_eq_Z4":z5,"Y2_GeV2":Y2,
            "global_square_certificate":square_certificate,"quartic_on_aligned_t_03":quartic_ray,
            "static_response":slope_static,"pole_response":slope_pole,
            "static_relative_to_pole":slope_static/slope_pole}


def wilson(m_h,alpha):
    logratio=2*math.log(m_h/163)
    loop=alpha/(4*math.pi)
    return np.array([-.064*(1+loop*(4.5-3*logratio)),
                     -.107*loop*(.125-logratio/12),
                     -.11*(-1.5*loop),.18*(1-loop)])


def exchange(M,point,c_ref,common):
    r,d,gR,gL,alpha=common
    c,s=point["c"],point["s"]
    Ph,PH=wilson(MH,alpha),wilson(M,alpha)
    results={}
    for name,edge in [("real",r),("phase",d)]:
        data={}
        for calibration,ref in [("fixed_heavy_yukawa",c_ref),("fixed_light_vertex",c)]:
            # The prefactor recovers the round-1648 +2 physical light FC vertex at the reference point.
            rr=2*edge*gR/ref
            ll=2*edge.conjugate()*gL/ref
            def amp(hadron,m,rh,lh):
                return -(hadron*np.array([.5*(rh*rh+lh*lh)]*2+[rh*lh]*2)).conjugate()/m**2
            h=amp(Ph,MH,c*rr,c*ll)
            H=amp(PH,M,-s*rr,-s*ll)
            A=amp(PH,M,1j*rr,-1j*ll)
            both=h+H+A
            directed=-(np.r_[((rr**2+ll**2)/2)*c*c*(Ph[:2]/MH**2-PH[:2]/M**2),
                               (rr*ll)*(c*c*Ph[2:]/MH**2+(2-c*c)*PH[2:]/M**2)]).conj()
            np.testing.assert_allclose(both,directed,rtol=2e-12,atol=1e-27)
            zs=sum(both)
            axis=zs/abs(zs)
            proj=min(sum(np.array(margins)*both/axis).real
                     for margins in itertools.product((.8,1.2),repeat=4))
            yb,err=read("../1386/inputs.json")["yukawa"]["down"][2]
            low=(.95*(1-2*err/yb))**2*proj
            def pair(x):
                return [x.real,x.imag]
            data[calibration]={"h":pair(sum(h)),"H":pair(sum(H)),"A":pair(sum(A)),
                               "total":pair(zs),"absolute_total_over_h":abs(zs/sum(h)),
                               "LR_enhancement":1+(2-c*c)*MH**2/(c*c*M**2),
                               "shared_response_lower_abs_total_GeV":low,
                               "physical_h_path_slope":2*c/ref}
        results[name]=data
    return results


def read(path):
    return json.loads((HERE/path).read_bytes())


def clean(obj):
    if isinstance(obj,dict): return {k:clean(v) for k,v in obj.items()}
    if isinstance(obj,list): return [clean(v) for v in obj]
    return float(f"{obj:.11g}")


def solve():
    # Reuse the already independently checked 163-GeV chiral and phase dictionary, never a 1651 result.
    common=runpy.run_path(str(HERE.parent/"1649"/"independent_check.py"))["couplings"]()
    c0=.05
    p=potential(600,c0)
    sc=-p["Z6"]*V**2/(1000**2-MH**2)
    c1=math.sqrt(2)*sc/math.sqrt(1+math.sqrt(1-4*sc*sc))
    points={"600_fixed":p,"1000_fixed_Z6":potential(1000,c1)}
    values={}
    for name,M in [("600_fixed",600),("1000_fixed_Z6",1000)]:
        cfg=points[name]
        assert cfg["Y2_GeV2"]>0 and cfg["global_square_certificate"]>0
        values[name]=exchange(M,cfg,c0,common)
        assert all(values[name][tag]["fixed_heavy_yukawa"]["LR_enhancement"]>30 for tag in ("phase","real"))
    unstable=potential(1000,c0)
    assert unstable["quartic_on_aligned_t_03"]<0
    return clean({"global_points":points,"unchecked_local_only_alternative":unstable,
                  "coherent_amplitudes_GeV":values})


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=solve()
    path=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1651 independent PASS: two complete-square vacua, static/pole transfer, coherent h/H/A kaon matching")
    else:
        with path.open("x",encoding="utf8") as f:
            f.write(json.dumps(result,indent=2)+"\n")

