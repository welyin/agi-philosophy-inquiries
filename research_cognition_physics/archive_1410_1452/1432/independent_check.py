"""Ray geometry, fixed-Dirac-phase ellipse, and two distinct kinds of charged-port repairs."""
import argparse
import itertools
import json
from math import acos, asin, atan, cos, degrees, pi, sin, sqrt, tan
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def get(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def box(row,mult):
    return [row[0]-mult*row[1],row[0]+mult*row[-1]]


def pmns(a,b,c,d):
    s,t,r=map(sin,(a,b,c))
    C,B,A=map(cos,(a,b,c))
    p=np.exp(1j*d)
    return np.array([[C*B,s*B,t/p],[-s*A-C*t*r*p,C*A-s*t*r*p,r*B],
                     [s*r-C*t*A*p,-C*r-s*t*A*p,A*B]])


def repairs(te,h,z):
    # Actual orthonormal mass columns at delta=0 and the cancelling first 12 phase.
    U=pmns(.55,h,z,0.)
    row=cos(te)*U[0]-sin(te)*U[1]
    t=row[2].real
    assert t>0
    tau=atan(t/abs(U[2,2]))
    corrected=cos(tau)*row-sin(tau)*U[2]
    projected=row.copy()
    projected[2]=0.
    projected/=np.linalg.norm(projected)
    assert abs(corrected[2])<1e-14
    assert np.isclose(abs(np.vdot(row,projected)),sqrt(1-t*t))
    return [t,degrees(asin(t)),degrees(tau),degrees(atan(tan(h)/sin(z))-te)]


def calc():
    cfg=get("inputs.json")
    low=get(cfg["charged_low"])
    reactor=get(cfg["reactors"])
    m1,m2=low["yukawa"]["charged_leptons"][:2]
    k=low["NO"]["atmospheric_sine_squared"]
    db=reactor["daya_bay"]
    rt=[db["sin_sq_2theta13"],db["one_sigma"]]
    js=reactor["juno_primary"]
    solar=asin(sqrt(js["sin_sq_theta12"]))
    out={}
    for mult in (0,2):
        aa,bb,kk,dd=[box(r,mult) for r in (m1,m2,k,rt)]
        ebox=[atan(sqrt(aa[0]/bb[1])),atan(sqrt(aa[1]/bb[0]))]
        hbox=[.5*asin(sqrt(i)) for i in dd]
        zbox=[asin(sqrt(i)) for i in kk]
        assert sin(zbox[0])*tan(hbox[0])>tan(ebox[1])
        allrows=np.array([repairs(*xs) for xs in itertools.product(ebox,hbox,zbox)])
        out[str(mult)]={"min_intrinsic_third_amplitude":float(allrows[:,0].min()),
                   "min_extra_arbitrary_ray_angle_deg":float(allrows[:,1].min()),
                   "min_extra_left_13_angle_deg":float(allrows[:,2].min()),
                   "min_extra_pure_12_angle_deg":float(allrows[:,3].min())}
    te=atan(sqrt(m1[0]/m2[0]))
    h=.5*asin(sqrt(rt[0]))
    z=asin(sqrt(k[0]))
    ellipse={}
    for d in (0.,.8,pi/2):
        U=pmns(solar,h,z,d)
        e,m=U[:2]
        base=cos(te)**2*abs(e)**2+sin(te)**2*abs(m)**2
        coeff=sin(2*te)*e.conj()*m
        C=np.column_stack((coeff.real,-coeff.imag))[1:,:]
        J=(e[1]*m[2]*e[2].conj()*m[1].conj()).imag
        assert np.isclose(abs(np.linalg.det(C)),(sin(2*te)**2)*abs(J),atol=1e-14)
        for phi in np.linspace(0,2*pi,19):
            row=cos(te)*e+np.exp(1j*phi)*sin(te)*m
            p=abs(row)**2
            assert np.allclose(p[1:]-base[1:],C@np.array([cos(phi),sin(phi)]))
            if abs(J)>1e-5:
                assert np.isclose(np.linalg.norm(np.linalg.solve(C,p[1:]-base[1:])),1.)
        phi=pi-d
        row=cos(te)*e+np.exp(1j*phi)*sin(te)*m
        assert np.isclose(abs(row[2]),out["0"]["min_intrinsic_third_amplitude"])
        ellipse[f"{d:.10g}"]={"J":float(J),"area_enclosed_in_p2_p3":pi*abs(np.linalg.det(C)),
                             "solar_deg_at_minimum_t":degrees(atan(abs(row[1]/row[0])))}
    # Construct upper and lower boundaries at each t by rephasing in the projected plane.
    for psi in np.linspace(0.,pi,13):
        A=cos(te)*cos(h)-sin(te)*sin(z)*sin(h)*np.exp(1j*psi)
        for sign in (1,-1):
            phi=np.angle(A)+(0 if sign==1 else pi)
            U=pmns(solar,h,z,psi-phi)
            row=cos(te)*U[0]+sin(te)*np.exp(1j*phi)*U[1]
            radius=asin(sin(te)*cos(z)/sqrt(1-abs(row[2])**2))
            assert np.isclose(atan(abs(row[1]/row[0])),solar+sign*radius)
    return {"repairs":out,"fixed_delta_ellipses":ellipse}


def matches(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(matches(a[k],b[k]) for k in a)
    return np.allclose(a,b,rtol=1e-10,atol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group()
    flags.add_argument("--write",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calc()
    if args.write:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert matches(result,get("independent_results.json"))
        for k,val in result["repairs"].items():
            author=get("results.json")[k]
            assert np.isclose(author["intrinsic_s13"][0],val["min_intrinsic_third_amplitude"])
            assert np.isclose(author["intrinsic_zero13_extra_charged_angle_deg"],val["min_extra_pure_12_angle_deg"])
        print("1432 joint extrema, fixed-phase ellipse and sharp rotation costs verified")
