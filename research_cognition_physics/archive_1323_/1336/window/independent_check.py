"""Independent checks for the two-port common-isometry constraint, round 1336."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

def shape_angle_bound(rho_u, rho_d):
    denom = math.sqrt(1.0-rho_u**2)-rho_d
    if denom <= 0:
        return None
    eta = (rho_d+rho_u)/denom
    if eta+rho_u >= 1:
        return None
    return eta/(1-eta-rho_u)

def finite_bound(rho_u, rho_d, k_u, k_d):
    shape = shape_angle_bound(rho_u*math.hypot(1,k_u),rho_d*math.hypot(1,k_d))
    assert shape is not None and 0 < k_u < 1 and 0 < k_d < 1
    return shape + (k_u*k_d)**2/(4*(1-k_u**2)*(1-k_d**2))

def root_pos(a):
    vals,vecs=np.linalg.eigh(a)
    assert vals.min()>-1e-12
    return (vecs*np.sqrt(np.maximum(vals,0)))@vecs.conj().T

def random_frame(rng,n):
    m=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    q,_=np.linalg.qr(m)
    return q

def analyze_sector(frame,port,k):
    n=3
    g=np.kron(port.conj()[None,:],np.eye(n))
    A=g@frame
    _,values,vh=np.linalg.svd(A)
    v=vh.conj().T
    orthogonal=np.linalg.qr(frame,mode='complete')[0][:,n:]
    mass=np.column_stack([orthogonal,k*g.conj().T])
    u,sing,_=np.linalg.svd(mass)
    order=np.argsort(sing)[:n]
    masses=sing[order]
    actual_light=u[:,order]
    ratio_leading=values[1]/values[0]
    ratio_full=masses[1]/masses[2]
    assert ratio_leading<=ratio_full*math.hypot(1,k)+1e-10

    vals=values[::-1]**2
    discriminant=np.sqrt((1+k*k)**2-4*k*k*vals)
    lambda_expected=2*k*k*vals/(1+k*k+discriminant)
    assert np.allclose(masses**2,lambda_expected,atol=7e-15,rtol=1e-12)
    leakage=np.linalg.norm(orthogonal.conj().T@actual_light,ord=2)
    assert leakage<=k**2/(2*(1-k*k))+2e-12
    return v,actual_light,ratio_leading,ratio_full

def compute():
    rng=np.random.default_rng(61733241)
    n=3
    residual=0.
    max_full=0.
    used=0
    for size in (.03,.24,.85):
        for rho in (.003,.014,.037):
            for fraction in (0.,.01,.1,.7):
                # Actual hierarchy and common coordinate rotation vary independently.
                top=random_frame(rng,n)@np.diag([size,size*rho,size*rho*.37])@random_frame(rng,n).conj().T
                bottom=random_frame(rng,n)@root_pos(np.eye(n)-top.conj().T@top)
                W=np.vstack([top,bottom])
                assert np.linalg.norm(W.conj().T@W-np.eye(3),ord=2)<3e-13
                sine=fraction*size*rho
                phase=rng.uniform(0,2*np.pi)
                ports=[np.array([1.,0.]), np.array([math.sqrt(1-sine**2),np.exp(1j*phase)*sine])]
                for ku,kd in ((.08,.16),(.25,1/3),(1/3,1/3)):
                    vu,uu,ru,fu=analyze_sector(W,ports[0],ku)
                    vd,ud,rd,fd=analyze_sector(W,ports[1],kd)
                    angle=np.linalg.norm(vu[:,1:].conj().T@vd[:,0])
                    weak_norm=np.linalg.norm(uu[:,:2].conj().T@ud[:,2])
                    bound=shape_angle_bound(ru,rd)
                    assert bound is not None and angle <= bound+3e-11
                    assert weak_norm<=finite_bound(fu,fd,ku,kd)+3e-11
                    gap=ku**2*kd**2/(4*(1-ku**2)*(1-kd**2))
                    assert weak_norm<=angle+gap+3e-11
                    residual=max(residual,float(angle-bound),float(weak_norm-angle-gap))
                    max_full=max(max_full,float(weak_norm))
                    used+=1

    ku=kd=1/3
    physical_upper=[.703/168.7,.071/2.80]
    promoted=[physical_upper[0]*math.hypot(1,ku),physical_upper[1]*math.hypot(1,kd)]
    lower_ckm=.0416-.0006
    bound=finite_bound(*physical_upper,ku,kd)
    assert bound<.0371<lower_ckm
    # Keep both sector spectra and readouts fixed. A third relation channel holds the residue.
    # A rotation on the *domain* changes the common left-handed heavy directions.
    first=.55*np.diag([.00015,.0035,1.])
    second=.55*np.diag([.0009,.018,1.])
    comp=[]
    for th in (0.,.08,.19):
        c,s=math.cos(th),math.sin(th)
        rot=np.array([[1,0,0],[0,c,-s],[0,s,c]])
        b=second@rot.T
        residue=root_pos(np.eye(n)-first.T@first-b.T@b)
        W=np.vstack([first,b,residue])
        assert np.allclose(W.T@W,np.eye(n),atol=2e-14,rtol=0)
        v1,u1,a1,f1=analyze_sector(W,np.array([1,0,0]),.1)
        v2,u2,a2,f2=analyze_sector(W,np.array([0,1,0]),.1)
        angle=float(np.linalg.norm(v1[:,1:].T@v2[:,0]))
        assert abs(angle-abs(s))<2e-14
        comp.append(dict(rotation=th,leading_mixing=angle,finite_left_current=float(np.linalg.norm(u1[:,:2].T@u2[:,2])),rho_full_u=f1,rho_full_d=f2))
    assert abs(comp[-1]["finite_left_current"]-comp[0]["finite_left_current"])>.18
    assert np.ptp([r["rho_full_u"] for r in comp])<1e-12
    assert np.ptp([r["rho_full_d"] for r in comp])<1e-12
    return dict(round=1336,status="PASS",seed=61733241,isometry_port_mass_comparisons=used,
                max_inequality_violation=residual,largest_tested_two_level_current=max_full,
                historical_data_source="Xing, Zhang and Zhou arXiv:0712.1419 Table II and (50)",
                historical_interval_comparison_scope="direct common-scale mass/left-current dictionary; component input ranges, no new-physics data fit",
                historical_upper_physical_ratios=physical_upper,
                finite_to_leading_upper_ratios=promoted,leading_angle_upper_from_mass_box=shape_angle_bound(*promoted),
                finite_current_upper=bound,direct_historical_Vcb_lower=lower_ckm,
                third_carrier_dimension_examples=comp,
                competitor_scope="arbitrary isometry allowed with extra relation port and matching residual; fixed per-sector singular spectra, free relative orientation, no autonomous formation or generation of real parameters")

def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for key in a: compare(a[key],b[key])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,abs_tol=1e-10,rel_tol=1e-8),(a,b)
    else:
        assert a==b,(a,b)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    val=compute()
    target=Path(__file__).with_name("independent_results.json")
    if args.write_new:
        with target.open("x",encoding="utf8",newline="\n") as f:
            json.dump(val,f,indent=2,ensure_ascii=False)
            f.write("\n")
    else:
        compare(val,json.loads(target.read_text(encoding="utf8")))
    print(json.dumps({key:val[key] for key in ("status","isometry_port_mass_comparisons","finite_current_upper","direct_historical_Vcb_lower")}))

