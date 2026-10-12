"""Local-density matching qualification. No UV completion or full-field simulation."""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss

BASE = Path(__file__).resolve().parent
OUT = BASE / "after1148_local_density_results.json"
MATERIAL = BASE.parent / "1147" / "material_bound.md"

def j0(x):
    return np.sinc(np.asarray(x)/np.pi)

def j1(x):
    x=np.asarray(x,dtype=float)
    y=np.empty_like(x)
    small=np.abs(x)<1e-3
    z=x[small]
    y[small]=z/3-z**3/30+z**5/840-z**7/45360
    z=x[~small]
    y[~small]=np.sin(z)/z**2-np.cos(z)/z
    return y

def ground(m):
    lo,hi=2.0,3.0
    for _ in range(100):
        q=(lo+hi)/2
        en=math.hypot(q,m)
        f=float(j0([q])[0]-q/(en+m)*j1([q])[0])
        if f>0: lo=q
        else: hi=q
    q=(lo+hi)/2
    en=math.hypot(q,m)
    return q,en,q/(en+m)

def radial(n, m=1.001):
    q,en,b=ground(m)
    x,w=leggauss(n)
    r=(x+1)/2
    w=w/2
    F=j0(q*r)
    G=b*j1(q*r)
    N2=1/(4*np.pi*np.sum(w*r*r*(F*F+G*G)))
    rho=N2*(F*F-G*G)
    scalar=4*np.pi*np.sum(w*r*r*rho)
    qs=4*np.pi*np.sum(w*r*r*rho*np.sinh(r)/r)
    qrho=qs/scalar
    rb=r/2
    wb=w/2
    bump=np.exp(-1/(1-4*rb*rb))
    qb=float(np.sum(wb*rb*rb*bump*np.sinh(rb)/rb)/np.sum(wb*rb*rb*bump))
    a=.5
    qball=3*(a*math.cosh(a)-math.sinh(a))/a**3
    qg=qball*qb
    kernel=math.exp(-3)/(12*np.pi)
    ratio=(qrho/qg)**2
    angle=np.pi*(ratio-1)
    # Independent angular integral for the exterior potential of actual rho.
    z,v=leggauss(min(n,128))
    d=3.0
    dist=np.sqrt(d*d+r[:,None]**2+2*d*r[:,None]*z[None,:])
    yker=np.exp(-dist)/(4*np.pi*dist)
    potential=2*np.pi*np.sum(w*r*r*(rho/scalar)*(yker@v))
    predicted=kernel*qrho
    scalar_closed=(en/2+m*(en-1))/(en*(en-1)+m/2)
    row={
        "quadrature_order":n,
        "q":q,"E":en,"b":b,"N_squared":float(N2),
        "scalar_integral":float(scalar),
        "scalar_closed":float(scalar_closed),
        "scalar_formula_residual":float(abs(scalar-scalar_closed)),
        "normalized_density_yukawa_moment":float(qrho),
        "unnormalized_scalar_yukawa_moment":float(qs),
        "ball_yukawa_moment":float(qball),
        "bump_yukawa_moment":qb,
        "old_g_yukawa_moment":float(qg),
        "normalized_local_density_cross_kernel":float(kernel*qrho**2),
        "old_g_cross_kernel":float(kernel*qg**2),
        "matching_constant_a_over_kappa":float(qg/qrho),
        "charge_only_matching_cross_ratio":float(ratio),
        "asymptotic_original_gate_phase_mismatch":float(angle),
        "asymptotic_code_isometry_mismatch_after_local_phase_alignment":float(2*abs(np.sin(angle/2))),
        "exterior_potential_direct_angular":float(potential),
        "exterior_potential_factorized":float(predicted),
        "exterior_potential_residual":float(abs(potential-predicted)),
        "density_boundary_derivative":float(4*N2*j0([q])[0]**2*(1-en))
    }
    assert row["scalar_formula_residual"]<1e-12
    assert row["exterior_potential_residual"]<1e-12
    assert row["density_boundary_derivative"]<0
    assert np.min(rho)>0
    return row

def spin_group():
    q,en,b=ground(1.001)
    N2=radial(128)["N_squared"]
    pauli=[
        np.array([[0,1],[1,0]],complex),
        np.array([[0,-1j],[1j,0]],complex),
        np.array([[1,0],[0,-1]],complex)]
    beta=np.diag([1,1,-1,-1])
    rows=[]
    for r,theta,phi in [(.1,.2,.3),(.5,1.2,2.1),(.9,2.8,-.4)]:
        nhat=[math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta)]
        sn=sum(x*s for x,s in zip(nhat,pauli))
        F=float(j0([q*r])[0]); G=float(b*j1([q*r])[0])
        U=math.sqrt(N2)*np.vstack([F*np.eye(2),1j*G*sn])
        scalar=U.conj().T@beta@U
        target=N2*(F*F-G*G)*np.eye(2)
        err=float(np.linalg.norm(scalar-target,2))
        assert err<1e-14
        rows.append({"r":r,"theta":theta,"phi":phi,"matrix_scalar_residual":err})
    return {"pass":True,"kind":"actual_ground_doublet_pointwise_scalar_density","rows":rows}

def compute():
    rows=[radial(n) for n in [64,128,256,512]]
    keys=["normalized_density_yukawa_moment","old_g_yukawa_moment","matching_constant_a_over_kappa"]
    spread={k:max(v[k] for v in rows)-min(v[k] for v in rows) for k in keys}
    assert max(spread.values())<1e-12
    return {
        "qualification_only":True,"formal_round":None,"formal_check_increment":0,
        "material_source_sha256":hashlib.sha256(MATERIAL.read_bytes()).hexdigest(),
        "all_passed":True,
        "scope":"spin compression and finite static cross coefficient; not local parent existence",
        "groups":[spin_group(),{"pass":True,"kind":"radial_density_and_static_exterior_matching",
                              "rows":rows,"quadrature_spread":spread,
                              "interval_certified":False}],
        "heat_trace_evaluated":False,
        "full_sea_or_field_truncation_used":False,
        "self_energy_or_clock_shift_computed":False,
        "missing_for_clock_shift":"same finite-effective/renormalized prescription and local matching conditions",
        "parent_full_state_long_time_error_certified":False
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write-results",action="store_true")
    a=ap.parse_args()
    data=compute()
    output=json.dumps(data,ensure_ascii=False,indent=2)+"\n"
    if a.write_results:
        if OUT.exists() and OUT.read_text(encoding="utf-8")!=output:
            raise SystemExit("refuse overwrite different qualification result")
        OUT.write_text(output,encoding="utf-8")
        print("saved",str(OUT))
    else:
        assert json.loads(OUT.read_text(encoding="utf-8"))==data
        print("PASS: two qualification groups; formal increment 0")
if __name__=="__main__":
    main()
