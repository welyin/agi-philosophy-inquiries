"""Minimal rank-one correction on a tagged plane. Yukawas are dimensionless."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def pmatrix(t12,t13,t23,delta):
    a,b,c = np.sin([t12,t13,t23])
    A,B,C = np.cos([t12,t13,t23])
    z = np.exp(1j*delta)
    return np.array([[A*B, a*B, b/z],[-a*C-A*c*b*z, A*C-a*c*b*z,c*B],
                     [a*c-A*C*b*z,-A*c-a*C*b*z,C*B]])


def plane(v):
    return np.linalg.svd(v.conj().reshape(1,3),full_matrices=True)[2].conj().T[:,1:]


def increasing_svd(matrix):
    left,values,right = np.linalg.svd(matrix,full_matrices=False)
    return left[:,::-1],values[::-1],right[::-1,:]


def lift(D, charged, b, n, k):
    B, N = plane(b), plane(n)
    ld, ds, rd = increasing_svd(D@B)
    raw = np.diag(charged)@N
    c = np.linalg.svd(raw,compute_uv=False)[::-1]
    a = k*ds
    assert max(a[0],c[0]) <= min(a[1],c[1])*(1+1e-12)
    signs = np.where(c >= a,1.,-1.)
    perm = np.argsort(signs*a)
    old = (signs*a)[perm]
    new = np.sort(signs*c)
    assert np.all(new+1e-12>=old) and new[0]<=old[1]+1e-12
    size = np.array([-np.prod(x-new)/(x-old[1-i]) for i,x in enumerate(old)])
    assert np.min(size)>-1e-14
    col = np.sqrt(np.maximum(size,0))
    target = (ld*signs)[:,perm]@(np.diag(old)+np.outer(col,col))@rd[perm,:]

    ur,sr,vr = np.linalg.svd(raw,full_matrices=True)
    ut,st,vt = np.linalg.svd(target,full_matrices=True)
    assert np.max(abs(sr-st)) < 1e-13
    L = ut@ur.conj().T
    R = N@vr.conj().T@vt@B.conj().T+np.outer(n,b.conj())
    E = L@np.diag(charged)@R
    F = (target-k*D@B)@B.conj().T
    first = np.outer((E-k*D)@b,b.conj())
    error = np.linalg.norm(E-k*D-first-F,2)
    assert np.linalg.norm(R.conj().T@R-np.eye(3))<1e-12
    assert np.linalg.norm(R@b-n)<1e-12 and error<1e-13
    assert np.max(abs(np.linalg.svd(E,compute_uv=False)[::-1]-charged))<1e-13
    fvals=np.linalg.svd(F,compute_uv=False)
    exact=float(np.sum(abs(c-a)))
    assert fvals[1]<1e-13 and abs(fvals[0]-exact)<1e-13
    return {"k":float(k),"projected_down":ds.tolist(),"projected_charged_lepton":c.tolist(),
            "second_source_minimum_and_attained_norm":exact,"relative_to_muon":exact/charged[1],
            "first_source_norm":float(np.linalg.norm(first,2)),"max_reconstruction_error":float(error)}


def box(row,n):
    return row[0]-n*row[1],row[0]+n*row[1]


def lower_bounds(m,angles,osc,atm,widths,theta_degrees):
    (d,s,b,e,mu,tau)=[box(row,widths) for row in m]
    lo,hi=np.array([box(row,widths) for row in angles]).T
    s12,s13,s23=np.sin(lo[:3]); c12,c13,c23=np.cos(hi[:3])
    tden=s12*s23*min(abs(math.sin(lo[3])),abs(math.sin(hi[3])))
    # Bottom fraction of w: use phase-space-free bounds on |V_jb/V_jd| first.
    dlows=np.array([c12*c13, s12*c23-np.sin(hi[1])*np.sin(hi[2]), tden])
    assert min(dlows)>0 and 0<lo[3]<hi[3]<math.pi/2
    bratio=np.array([np.sin(hi[1]),np.sin(hi[2]),1.])/dlows*d[1]/b[0]
    wb2max=bratio**2/(1+bratio**2)
    ratio=float(np.max(s[1]/(b[0]*np.sqrt(1-wb2max))))
    nelo=math.sqrt((1-osc["NO"]["s12"][2])*(1-osc["NO"]["s13"][2]))
    maxrot=math.radians(theta_degrees)
    nelo=max(0.,math.cos(min(math.pi/2,math.acos(nelo)+maxrot)))
    free=max(0.,nelo*mu[0]-ratio*tau[1])
    fixed={}
    feasible={}
    for o in ("NO","IO"):
        n=osc[o]
        ne=math.sqrt((1-n["s12"][2])*(1-n["s13"][2])) if o=="NO" else math.sqrt(n["s13"][1])
        ne=max(0.,math.cos(min(math.pi/2,math.acos(ne)+maxrot)))
        taumax=(math.sqrt(n["s12"][2]*atm[o][2])+math.sqrt((1-n["s12"][1])*(1-atm[o][1])*n["s13"][2])
                if o=="NO" else math.sqrt((1-atm[o][1])*(1-n["s13"][1])))
        taumin_angle=max(0.,math.acos(taumax)-maxrot)
        fixed[o]=max(0.,mu[0]*ne-s[1])+max(0.,tau[0]*math.sin(taumin_angle)-b[1])
        s2lo=b[0]*math.sqrt(1-max(wb2max))
        c2lo=tau[0]*math.sin(taumin_angle)
        feasible[o]=bool(s[1]<c2lo and mu[1]<s2lo)
    return {"quark_ratio_s1_over_s2_ceiling":ratio,"normal_free_k_norm_floor":free,
            "normal_free_k_in_units_of_actual_muon_floor":free/mu[1],"unit_k_norm_floor":fixed,
            "unit_k_interval_overlap_certified_for_whole_box":feasible}


def calculate():
    cfg=read(HERE/"inputs.json")
    def older(key):
        return read(HERE.parent/cfg[key])
    values=older("spectra")["scales"][cfg["high_scale"]]
    angles=older("quark")["high_scale"]["angles_theta12_theta13_theta23_delta"]
    osc,atm=older("oscillations"),older("atmospheric")["theta23_sin2"]
    V=pmatrix(*(row[0] for row in angles))
    y=np.array([row[0] for row in values[:3]])
    ell=np.array([row[0] for row in values[3:]])
    D=V*y
    examples={}
    for o,i in [("NO",0),("IO",2)]:
        n=osc[o]
        U=pmatrix(*(math.asin(math.sqrt(x)) for x in (n["s12"][0],n["s13"][0],atm[o][0])),
                  math.radians(cfg["example_lepton_phase_deg"][o]))
        null=U[:,i]
        examples[o]={}
        for j,label in enumerate(["u","c","t"]):
            b=V[j].conj()/y
            b/=np.linalg.norm(b)
            ds=np.linalg.svd(D@plane(b),compute_uv=False)
            cs=np.linalg.svd(np.diag(ell)@plane(null),compute_uv=False)
            # Match the larger of the two compressed singular values.
            examples[o][label]={"unit_k":lift(D,ell,b,null,1.),
                                 "optimum_k":lift(D,ell,b,null,cs[0]/ds[0])}
            # Exercise signed interlacing when one compressed mass rises and the
            # other falls, and when both decrease.
            slopes=cs/ds
            for ktest in (float(np.mean(slopes)),1.05*max(slopes)):
                lift(D,ell,b,null,ktest)
    bounds={str(angle):lower_bounds(values,angles,osc,atm,cfg["quark_widths"],angle)
            for angle in cfg["neutral_transport_examples_deg"]}
    assert bounds["0"]["normal_free_k_norm_floor"]>0
    assert all(bounds["0"]["unit_k_interval_overlap_certified_for_whole_box"].values())
    return {"full_box":bounds,"central_exact_examples":examples}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    dest=HERE/"results.json"
    if args.check:
        def close(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(close(v,b[k]) for k,v in a.items())
            if isinstance(a,list):
                return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
            return math.isclose(a,b,abs_tol=3e-13,rel_tol=2e-10)
        assert close(result,read(dest))
    else:
        with dest.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    print(result["full_box"]["0"])
