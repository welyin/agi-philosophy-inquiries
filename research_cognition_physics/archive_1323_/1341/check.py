"""Reproduce round 1341: inverse rays, historical calibration, rational box."""
from __future__ import annotations
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np

def unit(v):
    return v / np.linalg.norm(v)

def spectrum(matrix):
    return np.linalg.eigh((matrix+matrix.conj().T)/2)

def adjacent_rays(A,B):
    vals,vecs=spectrum(A+B)
    if vals.min()<=1e-10:
        raise ValueError("The supports must span the three dimensional flavor space.")
    root=(vecs*np.sqrt(vals))@vecs.conj().T
    inv=np.linalg.inv(root)
    lam,E=spectrum(inv@A@inv)
    V=root@E
    norms=np.linalg.norm(V,axis=0)
    rays=V/norms
    return rays,lam,norms**2

def overlap(V):
    G=np.abs(V.conj().T@V)
    return G[0,1],G[1,2]

def ckm(delta):
    s13=.00431
    c13=math.sqrt(1-s13*s13)
    s12=.2257/c13
    s23=.0416/c13
    c12=math.sqrt(1-s12*s12)
    c23=math.sqrt(1-s23*s23)
    p=np.exp(1j*delta)
    return np.array([[c12*c13,s12*c13,s13/p],
            [-s12*c23-c12*s23*s13*p,c12*c23-s12*s23*s13*p,s23*c13],
            [s12*s23-c12*c23*s13*p,-c12*s23-s12*c23*s13*p,c23*c13]])

def beta_observable(delta):
    V=ckm(delta)
    beta=np.angle(-V[1,0]*V[1,2].conjugate()/(V[2,0]*V[2,2].conjugate()))
    return np.sin(2*beta)

def bisect(f,left,right):
    lo=f(left)
    assert lo*f(right)<0
    for _ in range(60):
        mid=(left+right)*.5
        if lo*f(mid)>0: left=mid
        else: right=mid
    return (left+right)*.5

def h(z,eps):
    return (1-eps)/math.sqrt((1+z*z)*(1+(eps/z)**2))

def rational_box():
    q=Fraction
    # Loose rational enclosures of standard CKM sines and cosines.
    assert 1-q("0.00461")**2>q("0.9999")**2
    assert q(".2236")>q(".223") and q(".2278")/q(".9999")<q(".229")
    assert q(".0410")>q(".0409") and q(".0422")/q(".9999")<q(".0423")
    assert 1-q(".0423")**2>q(".99")**2
    assert q(".223")*q(".99")-q(".0423")*q(".00461")>q(".22")
    assert q(".229")+q(".0423")*q(".00461")<q(".23")
    assert (q(".223")*q(".0409")-q(".00461"))/q(".23")>q(".019")
    assert (q(".229")*q(".0423")+q(".00461"))/q(".22")<q(".07")
    assert q(".00461")/q(".2236")<q(".0207")
    eu=(q(".703")-q(".00085"))/(q("168.7")-q(".00177"))
    ed=(q(".04")-q(".00414"))/(q("2.98")-q(".00171"))
    assert eu<q(".0042") and ed>q(".012")
    def u_lower_squared(e):
        return (1-e)**2/((1+q(".07")**2)*(1+(e/q(".019"))**2))
    hi_d_sq=1/(1+(q(".012")/q(".0207"))**2)
    assert u_lower_squared(q(".0042")) > q(".96")**2 and hi_d_sq < q(".87")**2
    assert u_lower_squared(q(".009")) > q(".87")**2
    assert q(".009")*q("168.7")>q("1.51")
    return dict(h_up_minimum_guarantee=.96,h_down_maximum_guarantee=.87,
                held_out_charm_prediction_lower_GeV=1.51,observed_charm_upper_GeV=.703,
                cp_phase_scope="all real; rational inequalities, not sampled confidence estimates")

def compute():
    rng=np.random.default_rng(1341)
    maximum=0
    accepted=0
    for _ in range(32):
        data=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        U=np.linalg.qr(data)[0]
        r=unit(rng.normal(size=3)+1j*rng.normal(size=3))
        W=np.column_stack([r,U@r,U@U@r])
        a,b,c,d=rng.uniform(.1,3,4)
        A=(W*np.array([a,b,0]))@W.conj().T
        B=(W*np.array([0,c,d]))@W.conj().T
        F,val,w=adjacent_rays(A,B)
        delta=np.abs(np.array(overlap(F))-abs(np.vdot(r,U@r))).max()
        maximum=max(maximum,float(delta))
        assert delta<2e-10 and val[0]>-1e-10 and abs(val[2]-1)<1e-10
        assert np.linalg.norm((F*(val*w))@F.conj().T-A)<1e-11
        assert np.linalg.norm((F*((1-val)*w))@F.conj().T-B)<1e-11
        accepted+=1
    mu,mc,mt=.00127,.619,171.7
    md,ms,mb=.0029,.055,2.89
    A=np.diag([0,mc-mu,mt-mu])
    roots=[bisect(lambda x:beta_observable(x)-.687,l,r) for l,r in ((.52,.79),(1.57,1.84))]
    calibrations=[]
    cmp=0
    for phase in roots+list(np.linspace(-np.pi,np.pi,65)):
        V=ckm(phase);B=(V*np.array([0,ms-md,mb-md]))@V.conj().T
        rays,ev,w=adjacent_rays(A,B); gd,gu=overlap(rays)
        zu=abs(V[2,0]/V[1,0]); zd=abs(V[0,2]/V[0,1])
        eu=(mc-mu)/(mt-mu);ed=(ms-md)/(mb-md)
        assert abs(gd-h(zd,ed))<2e-11 and abs(gu-h(zu,eu))<2e-11
        assert gu>.96 and gd<.87
        if phase in roots:
            prediction=mu+(mt-mu)*bisect(lambda e:h(zu,e)-gd,0,1)
            calibrations.append(dict(ckm_phase_radians=phase,sin2beta=float(beta_observable(phase)),
                    overlap_from_down=float(gd),overlap_from_up=float(gu),
                    prediction_for_charm_GeV=prediction,held_out_charm_GeV=mc))
        cmp+=1
    return dict(round=1341,status="PASS",positive_orbit_cases=accepted,maximum_positive_error=maximum,
                full_matrix_comparisons=cmp,central_examples=calibrations,entire_box=rational_box())

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-10)
    else: assert a==b

if __name__=="__main__":
    p=argparse.ArgumentParser()
    modes=p.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-new",action="store_true");modes.add_argument("--check",action="store_true")
    args=p.parse_args();result=compute();path=Path(__file__).with_name("results.json")
    if args.write_new:
        with path.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,json.loads(path.read_text(encoding="utf8")))
    print(json.dumps(result,ensure_ascii=False))

