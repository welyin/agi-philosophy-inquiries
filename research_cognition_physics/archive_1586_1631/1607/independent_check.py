"""Arbitrary right-handed vectorlike doublets and an isosinglet control."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as Poly

HERE = Path(__file__).resolve().parent
OUT = HERE / "independent_results.json"
WEAK = HERE.parent / "1386" / "inputs.json"


def ckm(x):
    u,v,w = x[:3]
    a,b,c = np.sqrt(1-x[:3]**2)
    z = np.exp(1j*x[3])
    return np.array([[a*b,u*b,v/z],[-u*c-a*w*v*z,a*c-u*w*v*z,w*b],
                     [u*w-a*c*v*z,-a*w-u*c*v*z,c*b]])


def pair_floor(table, angle):
    ym,yp = [table[:,0]+sign*2*table[:,1] for sign in (-1,1)]
    sm,sp = [angle[:3,0]+sign*2*angle[:3,1] for sign in (-1,1)]
    cm = np.sqrt(1-sp**2)
    dmin = np.array([cm[0]*cm[1], sm[0]*cm[2]-sp[1]*sp[2], sm[0]*sm[2]-sp[1]])
    smax = np.array([sp[0],1.,sp[2]+sp[0]*sp[1]])
    bmax = np.array([sp[1],sp[2],1.])
    tangent = np.hypot(yp[0]*smax/ym[1],yp[0]*bmax/ym[2])/dmin
    angle_bound = np.arctan(tangent)
    return {l: math.cos(angle_bound[i]+angle_bound[j])
            for l,i,j in [("uc",0,1),("ut",0,2),("ct",1,2)]}


def pole_bound(q):
    x=Poly([0,1])
    c=1-2*x/3
    # Dropping mu_max <= 1 gives a rigorous bound for all q, even
    # when the optimum cannot be realized as an isometric compression.
    down_r=c*c*(q*q/(2*(1+q*q)))
    num=3*(2*(.5-2*x/3)**2+3*(.5-x/3)**2+down_r)
    den=(.5-x)**2+x*x
    ends=[.228,.234]
    stat=(num.deriv()*den-num*den.deriv()).roots()
    window=ends+[r.real for r in stat if abs(r.imag)<1e-11 and ends[0]<r.real<ends[1]]
    z=min(window,key=lambda a:num(a)/den(a))
    rho=num(z)/den(z)
    r=(1-q)/(1+q)
    mu_max=c(z)*(1+r)/(1+r*r)
    return {"tree_width_ratio_min":rho,"minimizing_sin2":z,"relaxed_optimum_mu_max":mu_max,
            "electroweak_pole_min_at_097":.97*rho, "escape_transfer_at_most":20.817/rho}


def charges_witness(q):
    # Singular strengths on the three ordinary gauge-singlet inputs.
    z=.234
    c=1-2*z/3
    r=(1-q)/(1+q)
    hi=c*(1+r)/(1+r*r)
    ev=np.array([hi,c,r*hi])
    assert 0 < min(ev) < max(ev) < 1
    t=np.diag(np.sqrt(ev))
    tags=np.array([[1,1],[0,0],[1,-1]])/np.sqrt(2)
    a,b=(t@tags).T
    rec=np.vdot(a,b)/np.linalg.norm(a)/np.linalg.norm(b)
    np.testing.assert_allclose(q,rec)
    full=np.vstack((t,np.diag(np.sqrt(1-ev))))
    np.testing.assert_allclose(full.T@full,np.eye(3))
    d_r=(ev-c)/2
    exact=3*(2*(.5-2*z/3)**2+3*(.5-z/3)**2+d_r@d_r)/((.5-z)**2+z*z)
    np.testing.assert_allclose(exact,pole_bound(q)["tree_width_ratio_min"])
    return {"compressed_singlet_strengths":ev.tolist(),"test_overlap":rec,
            "combined_down_right_strength":d_r@d_r,"full_Gram_relaxation_ratio":exact}


def singlet_control(v,down):
    # Choose a smooth, unequal-success completion (two hidden dimensions).
    b=v.conj().T/down[:,None]
    b/=np.linalg.norm(b,axis=0)
    gram=b.conj().T@b
    seed=np.sqrt([.8,.9,1.1])
    h=seed[:,None]*gram*seed
    success=seed**2/np.linalg.eigvalsh(h)[-1]
    t=b*np.sqrt(success)
    eig,vec=np.linalg.eigh(np.eye(3)-t.conj().T@t)
    extra=np.sqrt(eig[1:])[:,None]*vec[:,1:].conj().T
    tags=np.vstack((t,extra))
    u,_,_=np.linalg.svd(tags,full_matrices=True)
    ur=np.column_stack((tags,u[:,3:])).conj().T
    ul=np.eye(5,dtype=complex)
    ul[:3,:3]=v
    # SM Higgs vev derivative covers exactly the doublet rows.
    vev=246.
    m=np.r_[vev*down/np.sqrt(2),1100.,2100.]
    y=ul@np.diag(m)@ur.conj().T
    mask=np.diag([1.,1.,1.,0.,0.])
    np.testing.assert_allclose(ur.conj().T@ur,np.eye(5),atol=2e-14)
    np.testing.assert_allclose(ur.conj().T[:3,:3],t,atol=2e-13)
    np.testing.assert_allclose(ul.conj().T@y@ur,np.diag(m),atol=2e-11)
    np.testing.assert_allclose(ul.conj().T@mask@y@ur,np.diag(np.r_[m[:3],0,0]),atol=2e-11)
    np.testing.assert_allclose(ul.conj().T@mask@ul,mask,atol=2e-13)
    np.testing.assert_allclose(ur.conj().T@np.eye(5)@ur,np.eye(5),atol=2e-13)
    return {"light_fractions":success.tolist(),"gram_hidden":eig.tolist(),
            "gauge_basis_five_masses_GeV":m.tolist(),"singular_mass_error_GeV":float(max(abs(np.linalg.svd(y,compute_uv=False)-m[::-1]))),
            "three_overlap_moduli":[float(abs(gram[i,j])) for i,j in [(0,1),(0,2),(1,2)]],
            "cyclic_imag":(gram[0,1]*gram[1,2]*gram[2,0]).imag}


def calculate():
    data=json.loads(WEAK.read_text(encoding="utf8"))
    y=np.array(data["yukawa"]["down"])
    x=np.array(data["quark_sines_and_phase"])
    pairs=pair_floor(y,x)
    joint=min(pairs.values())
    assert joint > .7814
    test=pole_bound(joint)
    assert test["electroweak_pole_min_at_097"]>20.817
    low,high=.5,.9
    for _ in range(55):
        mid=(low+high)/2
        if pole_bound(mid)["electroweak_pole_min_at_097"]>20.817:
            high=mid
        else:
            low=mid
    return {"MZ_pair_angle_floor":pairs,"any_doublet_number_all_pairs":test,
            "overlap_needed_at_097":high,"propagation_angle_margin_degrees":math.degrees(math.acos(high)-math.acos(joint)),
            "q_090":pole_bound(.9),"gram_only_sharp_control_q090":charges_witness(.9),
            "two_weak_singlet_control":singlet_control(ckm(x[:,0]),y[:,0])}


def compare(a,b):
    if isinstance(b,dict):
        assert a.keys()==b.keys()
        for k in b:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=1e-10,atol=2e-11)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group()
    g.add_argument('--write',action='store_true')
    g.add_argument('--check',action='store_true')
    args=p.parse_args()
    current=calculate()
    if args.write:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(current,f,indent=2,allow_nan=False)
            f.write('\n')
    else:
        compare(current,json.loads(OUT.read_text(encoding='utf8')))
        print("1607 independent: all weak-doublet counts and a two-singlet control checked")
