"""Full flavour matrix evolution and bidirectional charged-frame shooting."""
import argparse
import json
import math
from pathlib import Path
import runpy
from functools import lru_cache

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


cfg=read("inputs.json")
physics=runpy.run_path(str(HERE/"../1699/check.py"))
texture=runpy.run_path(str(HERE/cfg["texture_source"]))
low=read(cfg["low_masses_and_old_angles"])
TAN=read(cfg["soft_projection"])["ordinary_finite_matching"]["minimum_root_high"]
MZ=physics["weak"]["mu_ref_GeV"]
MS=physics["soft"]
END=cfg["upper_scale_GeV"]
LOOP=16*np.pi**2
I=np.eye(3)
fields=physics["spec"]


def central(key):
    return np.array(low["yukawa"][key])[:,0]


def angles(kind, vub_shift=0.):
    if kind=="old":
        return np.arcsin(np.array(low["quark_sines_and_phase"][:3])[:,0])
    vals=cfg["direct_semileptonic_moduli"]
    us,ub,cb=(vals[k][0] for k in ("Vus","Vub","Vcb"))
    ub+=vub_shift*np.linalg.norm(vals["Vub"][1:])
    cb13=math.sqrt(1-ub*ub)
    return np.arcsin([us/cb13,ub,cb/cb13])


def low_ckm(gamma, theta):
    a,b,c=theta
    phase=np.radians(gamma)+np.arcsin(np.cos(a)*np.tan(c)*np.sin(b)*np.sin(np.radians(gamma))/np.sin(a))
    c1,c2,c3=np.cos(theta)
    s1,s2,s3=np.sin(theta)
    p=np.exp(1j*phase)
    return np.array([[c1*c2,s1*c2,s2*p.conjugate()],
        [-s1*c3-c1*s3*s2*p,c1*c3-s1*s3*s2*p,s3*c2],
        [s1*s3-c1*c3*s2*p,-c1*s3-s1*c3*s2*p,c3*c2]])


@lru_cache(maxsize=4)
def segments(step):
    g=np.array([physics["weak"][k][0] for k in ("g1_GUT","g2","g3")])
    poles=fields["neutral_poles_GeV"]
    edges=sorted(set([MZ,MS,END]+poles+[m for m,_ in physics["thresholds"] if MZ<m<END]))
    groups=[]
    for a,b in zip(edges,edges[1:]):
        if a==MS:
            g=1/np.sqrt(1/g**2-np.array([0.,2.,3.])/(48*np.pi**2))
        susy=a>=MS
        bs=np.array([33/5,1,-3.] if susy else [41/10,-19/6,-7.])
        for m,bb in physics["thresholds"]:
            if m<=a:
                bs+=bb
        steps=max(1,math.ceil(math.log(b/a)/step))
        dt=math.log(b/a)/LOOP/steps
        series=1/(1/g[None,:]**2-2*np.arange(2*steps+1)[:,None]*dt*bs/2)
        groups.append((a,susy,dt,series))
        g=np.sqrt(series[-1])
    return groups


def beta(y,g,susy):
    # Rows are the left doublet. The three charged matrices have one padding column.
    h=y@y.conj().transpose(0,2,1)
    u,d,e,n=h
    a,b,c,f=h.trace(axis1=1,axis2=2).real
    g1,g2,g3=g
    if susy:
        terms=np.stack([3*u+d,3*d+u,3*e+n,3*n+e])
        tr=np.array([3*a+f-13*g1/15-3*g2-16*g3/3,
                     3*b+c-7*g1/15-3*g2-16*g3/3,
                     3*b+c-9*g1/5-3*g2,3*a+f-3*g1/5-3*g2])
    else:
        terms=1.5*np.stack([u-d,d-u,e-n,n-e])
        tr=3*a+3*b+c+f - np.array([17*g1/20+9*g2/4+8*g3,
               g1/4+9*g2/4+8*g3,9*g1/4+9*g2/4,9*g1/20+9*g2/4])
    return (terms+tr[:,None,None]*I)@y


def shoot(gamma,theta,frame,step=.42,d_z=0.):
    y=np.zeros((4,3,4),complex)
    y[0,:,:3]=np.diag(central("up"))
    down=central("down")
    down[0]+=d_z*low["yukawa"]["down"][0][1]
    y[1,:,:3]=low_ckm(gamma,theta)@np.diag(down)
    y[2,:,:3]=frame@np.diag(central("charged_leptons"))
    at_soft=None
    quartic=physics["quartic"]["lambda_paper"][0]/physics["quartic"]["lambda_paper_factor"]
    def dl(xx,g,la):
        h=xx@xx.conj().transpose(0,2,1)
        W=np.array([3,3,1,1])
        tr=W@h.trace(axis1=1,axis2=2).real
        tr2=W@(h@h).trace(axis1=1,axis2=2).real
        gy,g2=.6*g[0],g[1]
        return 24*la**2+(4*tr-9*g2-3*gy)*la-2*tr2+3/8*(2*g2**2+(g2+gy)**2)
    for a,susy,dt,gg in segments(step):
        if a==MS:
            at_soft=y.copy()
            cs=1/math.hypot(1,TAN)
            y/=np.array([TAN*cs,cs,cs,TAN*cs])[:,None,None]
        for j,m in enumerate(fields["neutral_poles_GeV"]):
            if m==a:
                y[3,:,j]=physics["Y0"][:,j]
        for k in range(0,len(gg)-1,2):
            g=gg[k+1]
            f1=beta(y,gg[k],susy)
            f2=beta(y+dt*f1/2,g,susy)
            f3=beta(y+dt*f2/2,g,susy)
            f4=beta(y+dt*f3,gg[k+2],susy)
            if not susy:
                d1=dl(y,gg[k],quartic)
                d2=dl(y+dt*f1/2,g,quartic+dt*d1/2)
                d3=dl(y+dt*f2/2,g,quartic+dt*d2/2)
                d4=dl(y+dt*f3,gg[k+2],quartic+dt*d3)
                quartic+=dt*(d1+2*d2+2*d3+d4)/6
            y=y+dt*(f1+2*f2+2*f3+f4)/6
    return y,at_soft,quartic


def basis(matrix):
    u,s,_=np.linalg.svd(matrix[:,:3])
    return u[:,::-1],s[::-1]


def extraction(y):
    (uu,su),(ud,sd),(ue,se)=[basis(a) for a in y[:3]]
    v=uu.conj().T@ud
    cb=math.sqrt(1-abs(v[0,2])**2)
    ts=np.arcsin([abs(v[0,1])/cb,abs(v[0,2]),abs(v[1,2])/cb])
    ga=np.degrees(np.angle(-v[0,0]*v[0,2].conjugate()/(v[1,0]*v[1,2].conjugate())))
    return np.r_[sd,ts],se,ga,ue,su


def texture_frame(y,k):
    q,l,g,cur,up=extraction(y)
    edges,_,roots=texture["parameters"](g,q,l)
    size,scale,_=roots[1]
    t=size*np.exp(1j*k*texture["phi"])
    charged=scale*texture["OLD"]["LOW"]["mtx"](edges,[1,1,3],t).T
    target,sv=basis(charged)
    residual=texture["residual"](g,k,1,q,l)
    return residual,cur,target,charged,size,scale


def mapped_frame(frame,actual,target):
    # Align each nondegenerate mass eigenline, discarding arbitrary eigenvector phases.
    ov=np.diag(target.conj().T@actual)
    target=target*(ov/abs(ov))
    return target@actual.conj().T@frame


def solve(k,kind,seed=60.,vub_shift=0.,step=.42,d_z=0.):
    g=seed
    frame=np.eye(3,dtype=complex)
    th=angles(kind,vub_shift)
    for i in range(18):
        y,soft,lam=shoot(g,th,frame,step,d_z)
        x,cur,target,charged,size,scale=texture_frame(y,k)
        move=mapped_frame(frame,cur,target)
        err=np.linalg.norm(move-frame)
        if abs(x)<5e-11 and err<5e-10:
            break
        high,l,gh,_,_=extraction(y)
        # The static partial derivative preconditions the shooting map. The zero
        # itself is always checked on a fresh fully transported state.
        der=(texture["residual"](gh+.001,k,1,high,l)-texture["residual"](gh-.001,k,1,high,l))/.002
        g-=x/der
        frame=move
    else:
        raise RuntimeError(("full matrix shooting failed",k,kind,i,g,x,err))
    q,l,gh,_,u=extraction(y)
    sv=basis(charged)[1]
    assert max(abs(sv/l-1))<2e-8
    assert max(np.linalg.norm(np.outer(cur[:,j],cur[:,j].conj())-
                              np.outer(target[:,j],target[:,j].conj())) for j in range(3))<1e-8
    base=read(cfg["soft_projection"])["soft_interface"].copy()
    base["yt"]=basis(soft[0])[1][-1]
    w=read(cfg["gamma_holdout"])["LHCb_gamma_2026_deg"]
    edges=texture["inputs_at_gamma"](gh,q)[0]
    down_texture=q[2]*texture["OLD"]["LOW"]["mtx"](edges,[1,1,1],1.)
    quark_test=np.zeros_like(y)
    quark_test[0,:,:3]=np.diag(u)
    quark_test[1,:,:3]=down_texture
    qq,_,gg,_,_=extraction(quark_test)
    np.testing.assert_allclose(qq,q,rtol=2e-10)
    assert abs(gg-gh)<1e-8
    _,_,right=np.linalg.svd(charged)
    return {"gamma_MZ_deg":g,"gamma_10PeV_deg":gh,"source_k":k,"large_modulus":size,
            "GJ_common_scale":scale,"max_charged_spectrum_relative_error":max(abs(sv/l-1)),
            "distance_in_direct_gamma_marginal_widths":abs(g-w[0])/(w[1] if g<w[0] else w[2]),
            "lambda_30TeV_full_flavour":lam,"remaining_zero_stop_quartic":lam-physics["finite_quartic"](base,TAN,0.),
            "frame_step_norm":err,"mass_invariant_residual":x,
            "Y_up_at_10PeV":u.tolist(),"Y_down_at_10PeV":q[:3].tolist(),
            "Y_charged_at_10PeV":l.tolist(),"CKM_angles_at_10PeV_rad":q[3:].tolist(),
            "charged_left_MZ":pack(frame),"charged_left_10PeV":pack(target),
            "charged_right_10PeV":pack(right.conj().T[:,::-1]),
            "down_texture_10PeV":pack(down_texture),"up_texture_diag_10PeV":u.tolist(),
            "charged_texture_10PeV":pack(charged),"neutral_columns_10PeV":pack(y[3]),
            "top_soft":basis(soft[0])[1][-1],"bottom_soft":basis(soft[1])[1][-1],
            "iterations":i+1}


def pack(x):
    return {"real":x.real.tolist(),"imag":x.imag.tolist()}


def main():
    outs={}
    for kind in ("direct","old"):
        outs[kind]={str(k):solve(k,kind) for k in (0,1,2)}
    fine=solve(1,"direct",step=.21)
    shifted={f"direct_k{k}":solve(k,"direct",d_z=-1.) for k in (1,2)}
    r=math.sqrt(7)/5
    def gain(q):
        return q*(1+r)/(q*q+r)
    max_gain=gain(math.sqrt(r))
    spectra={}
    for key,wit in shifted.items():
        s=wit["large_modulus"]
        disc=math.sqrt((1+r)**2-4*s*s*r)
        qlo,qhi=((1+r-disc)/(2*s),(1+r+disc)/(2*s))
        np.testing.assert_allclose([gain(qlo),gain(qhi),gain(r)],[s,s,1.],atol=3e-15)
        spectra[key]={"charged_bare_moduli":[qlo,qhi],"real_return_plus_branch_weights":
                      [(1+qlo)/2,(1+qhi)/2],"gain_shortfall_if_fixed_1_over_sqrt3":s/gain(1/math.sqrt(3))-1}
    shifts={}
    for k in (0,1,2):
        for n in (-1.,1.):
            v=solve(k,"direct",vub_shift=n)
            shifts[f"k{k}_Vub_widths_{n:g}"]={key:v[key] for key in ("gamma_MZ_deg","large_modulus")}
    return {"charged_endpoints":outs,"halfstep_gamma_displacement_deg":
            abs(fine["gamma_MZ_deg"]-outs["direct"]["1"]["gamma_MZ_deg"]),"tan_beta":TAN,
            "down_one_width_lower":shifted,"independent_Vub_sensitivity":shifts,
            "one_pole_with_unit_neutral":{"r":r,"max_gain":max_gain,"fixed_two_path_gain":gain(1/math.sqrt(3)),
                                          "witnesses":spectra}}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=main()
    dest=HERE/"results.json"
    if args.check:
        assert json.loads(dest.read_bytes())==result
        print("1702 verified")
    else:
        dest.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
