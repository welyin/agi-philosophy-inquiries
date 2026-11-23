"""Oppositely ordered flavor tensors on the actual last two colored mass rows."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())
CFG=read("inputs.json")
PORT=runpy.run_path(str(HERE/CFG["flavor"]))


def color_surface(b,F,a):
    y=b["yW"].copy()
    y[2]=a*F@y[2]
    y[3]=F@y[3]
    for scale,susy,dt,grid in reversed(b["segments"]):
        y,_=b["model"]["advance"](y,np.eye(3,dtype=complex),-dt,grid[::-1],b["pt"]["flow"])
        if scale==b["Mcolor"]:
            break
    else:
        raise RuntimeError("missing physical colored-doublet threshold")
    Y=y[1,:,:3]
    left,s,rightH=np.linalg.svd(Y)
    up=np.linalg.svd(y[0,:,:3])[0][:,-1]
    c=abs(up.conj()@left[:,-1])
    f=np.array([c,math.sqrt((1+c)/2)])
    right=rightH.conj().T[:,::-1]
    X=np.sqrt(f**-2-1)[:,None]*right[:,[0,2]].conj().T
    S=np.eye(3)-2*np.outer(right[:,1],right[:,1].conj())
    assert np.linalg.norm(X@S-X)<5e-15
    return Y,X,f


def window(fun,limit):
    limits=[]
    for sign in (-1,1):
        lo,hi=0.,1.
        assert abs(fun(sign*hi))>limit
        for _ in range(85):
            mid=(lo+hi)/2
            if abs(fun(sign*mid))>limit:hi=mid
            else:lo=mid
        limits.append(sign*lo)
    return limits


def mat(x):
    return [x.real,x.imag]


def evaluate(step):
    b=PORT["OLD"]["mother"](step)
    raw=[]
    for s in CFG["relative_signs"]:
        x,F,_,_=PORT["terms"](b,0,s)
        a,_=PORT["daughter"](b,F)
        raw.append((x.conj().T@x,*color_surface(b,F,a)))
    Q=b["right_quark_to_GJ"]
    # Right quarks and left leptons carry conjugate family indices in this identification.
    def convert(a):return Q@a.conj()@Q.conj().T
    le=np.linalg.svd(b["yW"][2,:,:3])[0][:,1]
    S=convert(np.eye(3)-2*np.outer(le,le.conj()))
    K0,K1=[convert(item[0]) for item in raw]
    np.testing.assert_allclose(K1,S@K0@S,atol=1e-13)
    Y,X,f=raw[1][1:]
    T=S@K0
    np.testing.assert_allclose(S@K1,T.conj().T,atol=1e-13)
    G=X.conj().T@X
    B=X@T@X.conj().T
    B0=(B+B.conj().T)/2
    Bi=(B-B.conj().T)/2
    imag_tr=np.trace(Bi).imag
    adj0=np.trace(B0)*np.eye(2)-B0
    imag_quad=np.trace(adj0@Bi).imag
    base,quad=np.linalg.det(B0).real,np.linalg.det(Bi).real
    eps=read(CFG["edm"])["conditional_edm"]["theta_limit_rad"]
    def z(h,lam):
        return 1.-h*np.trace(B0).real+h*h*(base+lam*lam*quad)+1j*h*lam*(-imag_tr+h*imag_quad)
    mass=b["Mcolor"]
    u=246.22/math.sqrt(2)
    def full(Y,X,T,h):
        vals,vs=np.linalg.eigh(np.eye(2)+X@X.conj().T)
        C=(vs/np.sqrt(vals))@vs.conj().T
        return np.block([[u*Y,h*u*Y@T@X.conj().T],[mass*C@X,mass*C]])
    flat=full(Y,X,T,0.)
    d0=np.linalg.det(flat)
    rel=[]
    for h in (-1.,1./(16*np.pi**2),1.):
        for lam in (-1.,-.375,0.,.375,1.):
            M=(T+T.conj().T)/2+lam*(T-T.conj().T)/2
            direct=np.linalg.det(full(Y,X,M,h))/d0
            rel.append(abs(direct-z(h,lam)))
            np.testing.assert_allclose(direct,z(h,-lam).conjugate(),atol=1e-12,rtol=0)
    assert max(rel)<3e-12
    assert np.linalg.eigvalsh(B0)[-1]<.002
    symmetric_window=window(lambda h:np.angle(z(h,1.)),eps)
    loop=1/(16*np.pi**2)
    J=X.conj().T@np.linalg.solve(np.eye(2)+X@X.conj().T,X)
    H=Y.conj().T@Y
    assert np.linalg.norm(J@H-H@J)<3e-14*np.linalg.norm(H)
    def coeffs(t):
        return np.array([np.trace(J@H).real,2*np.trace(J@H@t).real,np.trace(J@t.conj().T@H@t).real])
    plus,minus=[coeffs(t) for t in (T,T.conj().T)]
    np.testing.assert_allclose(plus[:2],minus[:2],atol=2e-16*plus[0])
    physical=[]
    for t,coef in ((T,plus),(T.conj().T,minus)):
        heavy=np.vstack([X.conj().T,np.eye(2)])@np.diag(f)
        np.testing.assert_allclose(heavy.conj().T@heavy,np.eye(2),atol=1e-13)
        hh=np.column_stack([Y,Y@t@X.conj().T])@heavy
        direct=np.linalg.norm(hh)**2
        np.testing.assert_allclose(direct,np.sum(coef),atol=1e-13*plus[0])
        physical.append(direct)
    # Uniform boundary on |h lambda| for amplitudes individually limited to unity.
    remax=1+abs(np.trace(B0).real)+abs(base)+abs(quad)
    remin=1-abs(np.trace(B0).real)-abs(base)-abs(quad)
    assert remin>.99 and abs(imag_tr)>abs(imag_quad)
    exact_one=window(lambda d:np.angle(z(1.,d)),eps)
    Ms=[full(Y,X,t,1e-5) for t in (T,T.conj().T)]
    masses=[np.linalg.svd(z0,compute_uv=False) for z0 in [flat]+Ms]
    comparison=raw[0]
    changedX=comparison[2]
    z_with_recomputed=np.linalg.det(np.eye(2)-changedX@T@changedX.conj().T)
    first7=[np.angle(z(loop**n,1.)) for n in (1,2,3)]
    assert min(abs(np.array(first7[:2])))>eps and abs(first7[2])<eps
    return dict(scales_weak_color_GeV=[b["Mweak"],mass],color_factors=f,
        endomorphism_Wilson_is_an_extra_input=True,conditional_theta_limit=eps,
        heavy_port_trace=mat(np.trace(B)),heavy_port_det=mat(np.linalg.det(B)),
        ordered_theta_at_eta_one=[np.angle(z(1.,1.)),np.angle(z(1.,-1.))],
        ordered_slope_from_forward_to_signed=-imag_tr,full_five_by_five_max_relative_cert_error=max(rel),
        full_theta_test_eta_1e5=[np.angle(z(1e-5,1)),np.angle(z(1e-5,-1))],
        ordered_eta_connected_interval=symmetric_window,imbalance_connected_interval_at_eta_one=exact_one,
        all_box_necessary_abs_h_lambda=math.tan(eps)*remax/(abs(imag_tr)-abs(imag_quad)),
        all_box_sufficient_abs_h_lambda=math.tan(eps)*remin/(abs(imag_tr)+abs(imag_quad)),
        norm_and_valid_reciprocal_Geigenvalues=np.linalg.eigvalsh(B0),
        possible_loop_count_scalings_not_computed_loops=first7,
        heavy_Higgs_sum_coefficients_in_h=dict(positive=plus,signed=minus),
        heavy_Higgs_summed_values_at_h_one=physical,
        heavy_Higgs_order_rate_difference_normalized_to_no_return=(plus[2]-minus[2])/plus[0],
        fractional_mass_changes_at_eta_1e5=[(m/masses[0]-1)[::-1] for m in masses[1:]],
        independent_lower_pole_rerun_norm_change_G=np.linalg.norm(changedX.conj().T@changedX-G),
        independent_lower_pole_rerun_theta_change=np.angle(z_with_recomputed/z(1.,1.)))


def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [clean(z) for z in x]
    if isinstance(x,(float,np.floating)):return float(f"{x:.11g}")
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    arg=parser.parse_args()
    c,f=[evaluate(s) for s in CFG["steps"]]
    err=abs(c["ordered_slope_from_forward_to_signed"]/f["ordered_slope_from_forward_to_signed"]-1)
    assert err<5e-5
    out=clean(dict(**f,step_halving_relative_slope=err))
    if arg.check:
        assert out==read("results.json")
        print("1752 verified complex-conjugate colored determinants, protected reciprocal return and heavy Higgs budget")
    else:
        (HERE/"results.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n",encoding="utf8")
