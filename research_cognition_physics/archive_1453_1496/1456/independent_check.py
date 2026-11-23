"""Takagi, finite-cell bounds, and a controlled triplet-flavor transfer test."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def pmns(x,y,z,delta):
    a,b,c=np.sqrt([x,y,z])
    t,u,v=np.sqrt([1-x,1-y,1-z])
    phase=np.exp(1j*delta)
    return np.array([[t*u,a*u,b/phase],
                     [-a*v-t*c*b*phase,t*v-a*c*b*phase,c*u],
                     [a*c-t*v*b*phase,-t*c-a*v*b*phase,v*u]])


def box(order):
    # Historical disappearance boxes; the absolute IO atmospheric upper
    # side comes from the lower side of the signed negative gap.
    small=lambda p:(1-math.sqrt(1-p))/2
    return [(.3036-2.05*.0064,.3036+2.05*.0064),
            (small(.0851-2*.0024),small(.0851+2*.0024)),
            ((.402,.624) if order=="NO" else (.401,.621)),
            (7.388e-5-2.05*7.8e-7,7.388e-5+2.05*7.8e-7),
            ((.002509-2*.000025,.002509+2*.000027) if order=="NO"
             else (.002482-2*.000025,.002482+2*.000026))]


def cell(order,box):
    x,y,z,s,a=box
    high_bg=s[1]*math.sqrt(x[1]*(1-x[1])*(1-z[0]))
    scale=a[0]-s[1]*x[1] if order=="NO" else a[0]+s[0]*x[0]
    p=math.sqrt(1-y[1])*(math.sqrt(z[0]*y[0])*scale-high_bg)
    solar=(s[0]*(1-3*x[1]*(1-y[0])),s[1]*(1-3*x[0]*(1-y[1])))
    air=(a[0]*(1-3*y[1]),a[1]*(1-3*y[0]))
    g=(solar[0]+air[0],solar[1]+air[1]) if order=="NO" else (solar[0]-air[1],solar[1]-air[0])
    return p,g


def cell_bounds(order):
    axes=[np.linspace(low,high,5) for low,high in box(order)]
    subaxes=[list(zip(axis[:-1],axis[1:])) for axis in axes]
    offset=5e-5
    best_exact=best_transfer=math.inf
    for sub in itertools.product(*subaxes):
        p,g=cell(order,sub)
        gleast=min(abs(v) for v in g)
        assert g[0]*g[1]>0 and gleast>4*offset
        gmost=max(abs(v) for v in g)
        best_exact=min(best_exact,p/gmost)
        best_transfer=min(best_transfer,(p-offset)/(gmost+4*offset))
    norm=math.sqrt(1.5e-13*192*math.pi*137.035999)*1.1663787e-5
    return {"certified_ratio": best_exact,
            "ratio_with_5e_minus5_eV2_transfer": best_transfer,
            "equal_1TeV_masses_width_GeV_with_transfer":norm*1000**3/(9*8*math.pi*best_transfer)}


def positive_root(b):
    w,v=np.linalg.eigh(b)
    assert min(w)>0
    return (v*np.sqrt(w))@v.T.conj()


def exact_complex_check(order):
    x,y,z,s,a=[(lo+hi)/2 for lo,hi in box(order)]
    a=a if order=="NO" else -a
    for delta in np.linspace(0,2*math.pi,13):
        u=pmns(x,y,z,delta)@np.diag(np.exp(1j*np.array([.41,-.13,.32])))
        a1=.032**2+max(0,-a)
        m=(u.conj()*np.sqrt([a1,a1+s,a1+a]))@u.T.conj()
        B=m.T.conj()@m
        P=math.sqrt(1-y)*(s*math.sqrt(x*(1-x)*(1-z))+
                          math.sqrt(y*z)*(a-s*x)*np.exp(-1j*delta))
        G=s*(1-3*x*(1-y))+a*(1-3*y)
        assert abs(B[0,1]-P)<2e-18 and abs(np.trace(B)-3*B[0,0]-G)<1e-17
        width=np.abs(m)**2/(4*math.pi)
        i,j=np.diag_indices(3)
        width[i,j]*=.5
        W=np.triu(width).sum()-1.5*(2*width[0,0]+width[0,1]+width[0,2])
        assert abs(8*math.pi*W-G)<1e-17


def erase_example():
    # CP-conserving input is one allowed branch, which permits positive real roots.
    x,y,z,s,a=[(lo+hi)/2 for lo,hi in box("NO")]
    u=pmns(x,y,z,math.pi).real
    a1=.03**2
    B=(u*[a1,a1+s,a1+a])@u.T
    observed=positive_root(B)
    changed=B.copy()
    changed[0,1]=changed[1,0]=0
    triplet=positive_root(changed)
    E=changed-B
    transfer=np.linalg.norm(triplet-observed,2)/np.linalg.norm(observed,2)
    assert transfer < .1
    assert np.isclose(np.ptp(np.linalg.eigvalsh(E))/2,abs(B[0,1]))
    assert abs((triplet.T@triplet)[0,1]) < 1e-17
    G=np.trace(changed)-3*changed[0,0]
    return {"m_lightest_eV":.03,"relative_mass_norm_change":transfer,
            "traceless_quotient_change_eV2":abs(B[0,1]),"unchanged_G_eV2":G}


def calculate():
    ans={}
    for order in ("NO","IO"):
        exact_complex_check(order)
        ans[order]=cell_bounds(order)
    ans["controlled_second_source"]=erase_example()
    return ans


def same(a,b):
    if isinstance(a,dict) and isinstance(b,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    return math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-15)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=calculate()
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf-8") as output:
            output.write(json.dumps(data,indent=2,allow_nan=False)+"\n")
    else:
        assert same(data,json.loads(target.read_text(encoding="utf-8")))
    print("PASS signed widths, coherent phases, finite-box and second-source tests")
