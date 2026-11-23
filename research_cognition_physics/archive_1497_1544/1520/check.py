"""Absolute joint constraint from two inclusive proton searches and a single lost direction."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
REF=1e15


def load(path):
    return json.loads((HERE/path).read_text(encoding="utf8"))


def rows(cfg,widths,weakest):
    old=load(cfg["previous_nucleons"])
    theta12,theta13=load(cfg["angles"])["high_scale"]["angles_theta12_theta13_theta23_delta"][:2]
    x=theta12[0]+widths*theta12[1]
    y=theta13[0]+widths*theta13[1]
    ud=math.cos(x)*math.cos(y)
    us=math.sin(x)*math.cos(y)
    def ff(entry,sign):
        return abs(entry[0])+sign*widths*sum(entry[1:])
    pion=ff(cfg["pion_W0_GeV2"],-1)
    us_d=ff(old["W0_GeV2"]["us_d"],-1)
    ud_s=ff(old["W0_GeV2"]["ud_s"],1 if weakest else -1)
    assert min(pion,us_d,ud_s)>0
    def scale(meson,tau):
        mp=old["mass_GeV"]["p"]
        f=mp*(1-(meson/mp)**2)**2/(32*math.pi)
        return math.sqrt(tau*cfg["year_seconds"]/cfg["hbar_GeV_s"]*f)/REF**2
    spi=scale(cfg["pion_mass_GeV"],cfg["pion_lifetime_years"])
    sk=scale(old["mass_GeV"]["K+"],old["proton"]["partial_lifetime_lower_yr"])
    return (spi*ud*pion,sk*us*ud_s,sk*ud*us_d)


def optimum(A,B,C,z):
    def residual(u,d,s,t):
        return z*z*(u*d+s*t)**2+(u*s-t*d)**2
    def evaluate(theta):
        d,s=math.cos(theta),math.sin(theta)
        return residual(A,d,s,0),residual(B,d,s,C)
    lo,hi=0.,math.atan2(C,B)
    if evaluate(lo)[0]>=evaluate(lo)[1]:
        theta=lo
    elif evaluate(hi)[0]<=evaluate(hi)[1]:
        theta=hi
    else:
        for _ in range(75):
            mid=(lo+hi)/2
            pi_rate,k_rate=evaluate(mid)
            if pi_rate<k_rate:
                lo=mid
            else:
                hi=mid
        theta=(lo+hi)/2
    val=max(evaluate(theta))
    if z==0:
        exact=(A*C)**2/((A+B)**2+C*C)
        np.testing.assert_allclose(val,exact,rtol=2e-14)
    return {"contact_scale_GeV":REF*val**.25,"lost_direction_angle_degrees":math.degrees(theta),
            "relative_to_experimental_limits_at_ref":list(evaluate(theta))}


def fixed_mass_example(cfg,abc):
    # Old fit fixes both right heavy directions. An independent real left rotation
    # drops out of the right mass eigenvectors and need not be refitted here.
    base=(HERE/cfg["parent"]).resolve()
    spec=importlib.util.spec_from_file_location("mass1516",base/"check.py")
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    core=mod.parent_core(json.loads((base/"inputs.json").read_text(encoding="utf8")))
    s,d,e,hda,hea,_=core
    phase=12.
    v,x,_=mod.mixing_core(core,phase)
    z=np.exp(1j*math.radians(phase))
    m=hda[1:]/abs(1+2*x*z)
    hd=np.r_[hda[0],m*(1+2*x*z)]
    he=np.r_[hda[0],m*(1-3*x*z)]
    hd,he=hd/np.linalg.norm(hd),he/np.linalg.norm(he)
    yd=np.column_stack([np.zeros(3),np.diag(s)])
    P=np.eye(4)-np.outer(hd,hd.conj())
    vr=(P@yd.conj().T@v)/d
    np.testing.assert_allclose(vr.conj().T@vr,np.eye(3),atol=3e-10)
    w=vr.conj().T@he
    weights=abs(w)
    overlap=abs(hd.conj()@he)
    np.testing.assert_allclose(weights@weights,1-overlap**2,atol=1e-9)
    # Use the actual accepted real up frame for the flavor witness. The
    # generic row bound above uses a separate high-scale CKM calibration.
    saved=json.loads((base/"results.json").read_text(encoding="utf8"))
    frame=saved["solutions"]["direct"]["12"]["one_passing_frame"]["O"]
    ckm=np.array(frame)@v
    ud,us=abs(ckm[0,0]),abs(ckm[0,1])
    a12,a13=load(cfg["angles"])["high_scale"]["angles_theta12_theta13_theta23_delta"][:2]
    udref=math.cos(a12[0])*math.cos(a13[0])
    usref=math.sin(a12[0])*math.cos(a13[0])
    A,B,C=abc[0]*ud/udref,abc[1]*us/usref,abc[2]*ud/udref
    pi=A*A*(1-weights[0]**2)
    k_floor=B*B+C*C-(B*weights[0]+C*weights[1])**2
    k_top=B*B+C*C-(B*weights[0]-C*weights[1])**2
    return {"matched_abs_Vud_Vus":[float(ud),float(us)],
            "loss_amplitudes_d_s_b":weights.tolist(),"smallest_gauge_singular_value":float(overlap),
            "pion_fraction_of_fully_visible_rate":float(1-weights[0]**2),
            "kaon_rate_at_ref_arbitrary_relative_phase":[float(k_floor),float(k_top)],
            "contact_scale_from_pion_GeV":REF*pi**.25,
            "absolute_bound_with_both_GeV":REF*max(pi,k_floor)**.25}


def calculation():
    cfg=load("inputs.json")
    center=rows(cfg,0,True)
    result={"reference_GeV":REF}
    for n in (0,1,2):
        abc=rows(cfg,n,True)
        result[f"worst_case_widths_{n}"]=optimum(*abc,0)
    result["central_if_z_known"]={str(z):optimum(*center,z) for z in (.05922517833,.5,1.)}
    result["previous_flavor_example_at_central_matching"]=fixed_mass_example(cfg,center)
    # Explicit geometry for a second lost direction, and attainable first-defect witnesses.
    E=np.eye(5)
    X=E[:,:3].T@E[:,[3,4,2]]
    np.testing.assert_equal(X,np.diag([0,0,1]))
    assert np.linalg.norm(np.array([[center[0],0,0],[center[1],center[2],0]])@X)==0
    return result


def check(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            check(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=5e-9,atol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-new",action="store_true")
    modes.add_argument("--check",action="store_true")
    opt=parser.parse_args()
    result=calculation()
    if opt.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        check(result,load("results.json"))
    print("Verified the two-source minimax, one-partner projection and matched weak-flavor witness.")
