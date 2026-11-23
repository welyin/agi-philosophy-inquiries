"""Single-spurion determinant response and CKM-calibrated simultaneous flavour bound."""
import argparse
import json
from math import asin, cos, sin, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def ckm(v):
    (s,t,r), phase = np.sin(v[:3]), np.exp(1j*v[3])
    c,b,a=np.cos(v[:3])
    return np.array([[c*b,s*b,t/phase],[-s*a-c*r*t*phase,c*a-s*r*t*phase,r*b],
                     [s*r-c*a*t*phase,-c*r-s*a*t*phase,a*b]])


def sine_box(data, low_scale, width):
    vals,errs=np.array(data)[:3].T
    lo,hi=vals-width*errs, vals+width*errs
    return (lo,hi) if low_scale else (np.sin(lo),np.sin(hi))


def gap_box(data,low_scale,width):
    lo,hi=sine_box(data,low_scale,width)
    s,t,r=lo
    S,T,R=hi
    c,b,a=np.sqrt(1-hi**2)
    C,B,A=np.sqrt(1-lo**2)
    floors=np.array([[c*b,s*b,t],[s*a-C*R*T,c*a-S*R*T,r*b],
                     [s*r-C*A*T,c*r-S*A*T,a*b]])
    assert np.min(floors)>0
    edge=(floors**2)@(floors**2).T
    lap=np.diag(edge.sum(axis=1))-edge
    gap=np.linalg.eigvalsh(lap)[1]
    return gap/(1+sqrt(1-gap))


def core(a):
    u=np.exp(1j*a)
    return np.array([[1,1,1],[1,u,u*u],[1,-1,1]])


def ratio_and_static_example(m):
    cen,err=np.array(m).T
    def ratio(y):
        return np.prod(y[:3])/np.prod(y[3:])
    R=ratio(cen)
    d=(np.pi/2,asin(1/R)) if R>1 else (asin(R),np.pi/2)
    val=abs(np.linalg.det(core(d[0]))/np.linalg.det(core(d[1])))
    assert np.isclose(val,R)
    step=2*err*np.array([1,1,1,-1,-1,-1])
    return {"central":R,"two_width":[ratio(cen-step),ratio(cen+step)],
            "unit_modulus_core_phases_rad":list(d)}


def explicit_witness(V):
    # Symmetric up tensor and a separately rotated down right frame.
    U=ckm([.51,.26,.39,1.12])
    T=np.diag([3.,2.,0.])
    B=np.diag([1.,0.,0.])
    matpairs=[]
    norms=0.
    for left,right,masses,rcharge in ((U,U.conj(),[.02,.25,.8],T),
                         (U@V,ckm([.45,.16,.3,.6]),[.01,.12,.65],B)):
        ml=np.diag(masses)
        l=left.conj().T@T@left
        r=right.conj().T@rcharge@right
        g=l@ml+ml@r
        norm=0.
        for i in range(3):
            for j in range(i+1,3):
                mi,mj=masses[i],masses[j]
                inv=(mj*g[i,j]-mi*g[j,i].conjugate())/(mj*mj-mi*mi)
                assert abs(inv-l[i,j])<1e-12
                norm+=2*abs(inv)**2
        norms+=norm
        matpairs.append(np.trace(np.linalg.solve(ml,g)).real)
    assert np.isclose(matpairs[1],np.trace(T)+np.trace(B))
    return norms


def results():
    cfg=read("inputs.json")
    low=read(cfg["low"])
    high=read(cfg["high"])["high_scale"]
    mass=read(cfg["volume"])["scales"]
    stats={}
    for scale, angles, m in (
        ("MZ",low["quark_sines_and_phase"],low["yukawa"]["down"]+low["yukawa"]["charged_leptons"]),
        ("GUT_SM",high["angles_theta12_theta13_theta23_delta"],mass["1e16"])):
        cent=np.array(angles)[:,0].copy()
        is_low=scale=="MZ"
        if is_low:
            cent[:3]=np.arcsin(cent[:3])
        V=ckm(cent)
        assert np.allclose(V@V.conj().T,np.eye(3))
        orth_gap=1-np.linalg.svd(abs(V)**2,compute_uv=False)[1]
        lower=gap_box(angles,is_low,cfg["width"])
        assert 0<lower<=orth_gap
        test=explicit_witness(V)
        assert test >= orth_gap*(9+4-25/3)
        stats[scale]={"volume":ratio_and_static_example(m),"projection_gap_central":orth_gap,
                      "rigorous_two_width_gap_floor":lower,
                      "total_offdiagonal_charge_norm_squared_floor":(2/3)*lower,
                      "at_least_one_reconstructed_pair_magnitude_floor":sqrt(lower/18)}
    # Holomorphic triangular support from charges t=(2,1,0), f=(-2,-1,0).
    power=np.array([2,1,0])[:,None]+np.array([-2,-1,0])[None,:]
    def mono(eps):
        return np.where(power>=0,eps**power,0.)
    assert np.isclose(np.linalg.det(mono(.21)),1)
    e1,e2=.12,.23
    exps=np.array([3,2,0])[:,None]+np.array([1,0,0])[None,:]
    cores=(core(.7), core(1.1))
    for c in cores:
        yy=(e2**exps)*c
        deriv=exps*yy
        assert np.allclose(np.trace(np.linalg.solve(yy,deriv)),6)
        h=1e-6
        central=((e2*np.exp(h))**exps*c-(e2*np.exp(-h))**exps*c)/(2*h)
        assert np.allclose(central,deriv)
    def ratio_at(e):
        return abs(np.linalg.det(e**exps*cores[0])/np.linalg.det(e**exps*cores[1]))
    assert np.isclose(ratio_at(e1),ratio_at(e2))
    # Non-holomorphic neutral insertion acts on d only, yet obeys the same charges.
    compare=((1+e2**2)/(1+e1**2))**3
    return {"same_scale":stats,"allowed_neutral_dressing_double_ratio_eps2_over_eps1":compare,
            "radial_dressing_trace_offset_at_eps2":6*e2**2/(1+e2**2)}


def matches(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(matches(a[k],b[k]) for k in a)
    return np.allclose(a,b,rtol=1e-10,atol=1e-12)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    modes=cli.add_mutually_exclusive_group()
    modes.add_argument("--write-new",action="store_true")
    modes.add_argument("--check",action="store_true")
    args=cli.parse_args()
    output=results()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as f:
            json.dump(output,f,indent=2)
            f.write("\n")
    else:
        assert matches(output,read("results.json"))
        print("1431 single flavon, cross-sector derivative and CKM projections verified")
