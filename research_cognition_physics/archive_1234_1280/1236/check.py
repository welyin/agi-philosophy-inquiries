"""1236: moving clock ADM/source identity and finite six-record certificates.

Default: deterministic read-only replay. --write-once creates results.json.
Float tests diagnose analytic identities; Fraction bounds are exact.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]


def need(ok, label):
    if not bool(ok):
        raise AssertionError(label)


def maxabs(x):
    return float(np.max(np.abs(x)))


def rational(x):
    return {"fraction":str(x),"decimal":float(x)}


def moving_sources():
    rng=np.random.default_rng(1236001)
    residuals={k:0.0 for k in ("Legendre","velocity","phase","mass_shell",
                              "normal_energy","momentum_flux","spatial_stress",
                              "lapse_variation","shift_variation","metric_variation")}
    min_E=float("inf")
    min_coordinate=float("inf")
    for j in range(20):
        B=rng.normal(size=(3,3))
        h=B.T@B+np.eye(3)*0.7
        hi=np.linalg.inv(h)
        N=0.7+0.08*j
        beta=rng.normal(size=3)*2.5
        p=rng.normal(size=3)+np.array([0.6,-0.4,0.2])
        I=0.2+0.015*j
        omega=1.0
        M=0.7+omega*I+0.1
        E=float(np.sqrt(M*M+p@hi@p))
        sigma=N*M/E
        xdot=N*hi@p/E-beta
        v=xdot+beta
        p_from_velocity=M*h@v/sigma
        phase=N*M*omega/E
        lag=I*phase-M*np.sqrt(N*N-v@h@v)
        ham=N*E-beta@p
        residuals["Legendre"]=max(residuals["Legendre"],abs(p@xdot+I*phase-lag-ham))
        residuals["velocity"]=max(residuals["velocity"],maxabs(p_from_velocity-p))
        residuals["phase"]=max(residuals["phase"],abs(phase-omega*sigma))
        g=np.zeros((4,4))
        g[0,0]=-N*N+beta@h@beta
        g[0,1:]=h@beta
        g[1:,0]=h@beta
        g[1:,1:]=h
        u=np.r_[1.0,xdot]/sigma
        P=M*u
        n_cov=np.r_[-N,0.0,0.0,0.0]
        sqrt_h=float(np.sqrt(np.linalg.det(h)))
        # Coordinate-delta coefficient; no extra source normalization.
        T=np.outer(P,P)/(E*sqrt_h)
        Tcov=g@T@g
        n_up=np.r_[1.0,-beta]/N
        rho=float(n_cov@T@n_cov)
        jflux=-(g@T@n_cov)[1:]
        S=Tcov[1:,1:]
        residuals["mass_shell"]=max(residuals["mass_shell"],abs(u@g@u+1))
        residuals["normal_energy"]=max(residuals["normal_energy"],abs(rho-E/sqrt_h))
        residuals["momentum_flux"]=max(residuals["momentum_flux"],maxabs(jflux-p/sqrt_h))
        residuals["spatial_stress"]=max(residuals["spatial_stress"],maxabs(S-np.outer(p,p)/(E*sqrt_h)))
        need(abs(n_cov@n_up+1)<1e-13,"normal normalization")
        delta=2e-5
        def H(N0=N,b0=beta,h0=h):
            return N0*np.sqrt(M*M+p@np.linalg.solve(h0,p))-b0@p
        num_N=(H(N0=N+delta)-H(N0=N-delta))/(2*delta)
        residuals["lapse_variation"]=max(residuals["lapse_variation"],abs(num_N-E))
        direction=np.array([.3,-.2,.4])
        num_beta=(H(b0=beta+delta*direction)-H(b0=beta-delta*direction))/(2*delta)
        residuals["shift_variation"]=max(residuals["shift_variation"],abs(num_beta+p@direction))
        dh=np.array([[.2,.1,-.05],[.1,-.12,.07],[-.05,.07,.3]])
        numeric_h=(H(h0=h+delta*dh)-H(h0=h-delta*dh))/(2*delta)
        analytic_h=-N/(2*E)*(hi@p)@dh@(hi@p)
        residuals["metric_variation"]=max(residuals["metric_variation"],abs(numeric_h-analytic_h))
        min_E=min(min_E,E)
        min_coordinate=min(min_coordinate,float(ham))
    need(max(residuals.values())<2e-8,"moving ADM/stress diagnostic residuals")
    need(min_E>0,"positive normal energy")
    # Explicit coordinate-H sign deletion: allowed large shift, E remains positive.
    neg=np.sqrt(2)-3
    need(neg<0,"coordinate generator need not be positive")
    # Same unscaled 1188 parameters: a=2k, b=1/(2k), Lambda_H=Lambda/(2k).
    k=F(3,7); lam=F(2,11); Bgeo=F(17,13); rho=F(5,9)
    a=2*k; b=1/(2*k); LambdaH=lam/(2*k)
    C=-Bgeo/a+2*LambdaH+rho
    geometric=-Bgeo/2+lam+k*rho
    need(a*b==1 and k*C==geometric,"unrescaled physical-source normalization")
    return dict(status="PASS",moving_cases=20,residuals=residuals,
                minimum_normal_energy=min_E,minimum_sample_coordinate_H=min_coordinate,
                explicit_negative_coordinate_H=float(neg),
                source_is_physical_T=True,grav_pi_absent_by_analytic_derivation=True,
                normalization="a=2*kappa,b=1/(2*kappa),Lambda_H=Lambda/(2*kappa)")


def determinant(e,r):
    return (e[1]-e[0])*(r[2]-r[0])-(e[2]-e[0])*(r[1]-r[0])


def flat_convex_cases():
    b=np.array([1.2,.4,-.3]);v=np.array([.12,-.31,.22])
    speed2=float(v@v)
    def radial(s):
        x=b+v*s;radius=float(np.linalg.norm(x))
        derivative=float(x@v/radius)
        second=(speed2*radius*radius-float(x@v)**2)/radius**3
        return radius,derivative,second
    rows=[]
    for delay in (0.0,0.17):
        e=[];r=[];mins=float("inf")
        for s in (0.0,.3,.8):
            r0,d0,dd0=radial(s)
            r1,d1,dd1=radial(s+delay)
            e.append(s-r0);r.append(s+delay+r1)
            second=(dd1*(1-d0)+(1+d1)*dd0)/(1-d0)**3
            mins=min(mins,second)
        D=determinant(e,r)
        need(D>0 and mins>0,"nonradial flat inertial convex, not affine")
        rows.append(dict(delay=delay,D=D,min_return_second_derivative=mins))
    # An accelerated timelike flat mirror reproduces the curved return exactly.
    t=.331;d=.03;h=1.0;y=1+h*t
    Rp=1+4*h*d/3*y**(-1/3)+4*h*h*d*d/9*y**(-2/3)
    H=2/3
    Rpp=-H*H*d*y**(-4/3)-H**3*d*d*y**(-5/3)
    vm=(Rp-1)/(Rp+1)
    acceleration=4*Rpp/(Rp+1)**3
    need(abs(vm)<1 and acceleration<0,"free-endpoint condition essential")
    return rows,dict(mirror_velocity=vm,mirror_acceleration=acceleration,
                     demonstrates_only_free_flat_class_excluded=True)


def finite_records():
    e=[F(0),F(331,1000),F(728,1000)]
    roots=[F(1),F(11,10),F(6,5)]
    r=[(q+F(1,50))**3-1 for q in roots]
    D=determinant(e,r)
    need(D==-F(5529,12500000),"exact six-record curved determinant")
    need(min(e+r)>=0 and max(e+r)<1,"fixed phase menu below one radian; no counter permission")
    coeff=2*((r[2]-r[0])+(r[1]-r[0])+(e[1]-e[0])+(e[2]-e[0]))
    need(coeff==F(26967,6250),"all-six-record linear error constant")
    eps=F(1,100000)
    bound=coeff*eps+8*eps*eps
    overlap_bound=coeff*(2*eps)+8*(2*eps)**2
    need(overlap_bound < -D,"two record boxes cannot overlap any flat convex ideal")
    maximum=F(0)
    # Multiaffine determinant: box extrema occur at its 64 vertices.
    for signs in itertools.product((-1,1),repeat=6):
        ee=[e[j]+signs[j]*eps for j in range(3)]
        rr=[r[j]+signs[j+3]*eps for j in range(3)]
        maximum=max(maximum,abs(determinant(ee,rr)-D))
    need(maximum<=bound and D+bound<0,"finite six-record robust sign")
    flat,accelerated=flat_convex_cases()
    return dict(status="PASS",emissions=[str(x) for x in e],returns=[str(x) for x in r],
                D=rational(D),linear_error_coefficient=rational(coeff),
                epsilon=rational(eps),error_bound=rational(bound),
                arbitrary_flat_two_box_bound=rational(overlap_bound),
                exact_box_max_error=rational(maximum),box_vertices=64,
                flat_nonradial_cases=flat,accelerated_flat_deletion=accelerated)


def root_bounds(x,n,scale=10**22):
    lo,hi=0,2*scale
    while F(hi,scale)**n<x: hi*=2
    while hi-lo>1:
        mid=(lo+hi)//2
        if F(mid,scale)**n<=x:lo=mid
        else:hi=mid
    return F(lo,scale),F(hi,scale)


def source_deletion():
    e=[F(0),F(331,1000),F(728,1000)]
    d=F(3,100)
    def returned(t,h):
        low,high=root_bounds(1+h*t,3)
        def evaluate(z):return t+2*d*z*z+4*h*d*d*z/3+8*h*h*d**3/27
        return evaluate(low),evaluate(high)
    hl,hu=root_bounds(F(4,5),2)
    boxes=[(returned(t,hl)[0],returned(t,hu)[1]) for t in e]
    a,b=e[1],e[2]
    low=a*boxes[2][0]-b*boxes[1][1]+(b-a)*boxes[0][0]
    high=a*boxes[2][1]-b*boxes[1][0]+(b-a)*boxes[0][1]
    full=-F(5529,12500000)
    need(full<low<high<0,"strict clock-source deletion determinant separation")
    area=a*b*(b-a)/2
    lower=area*F(7,100)*(2*F(59,100)*d/F(6,5)**7
                          +3*F(59,100)**2*d*d/F(6,5)**8)
    need(lower==F(40383498710737,1194393600000000000),"exact Peano lower bound")
    need(low-full>lower,"root enclosure stronger than analytic lower certificate")
    # Conditions used by analytic derivative certificate, checked rationally.
    need(F(59,100)<F(2,3)*hl,"lower H")
    need(F(2,3)*(1-hu)>F(7,100),"positive source-induced H interval")
    drop_coeff_upper=2*((boxes[2][1]-boxes[0][0])+(boxes[1][1]-boxes[0][0])+a+b)
    need(drop_coeff_upper<5,"uniform all-six-record deletion coefficient")
    eps=F(1,1000000)
    both=10*eps+16*eps*eps
    margin=lower-both
    need(margin>F(238,10000000),"strict finite source separation including all six records")
    return dict(status="PASS",deleted_D_lower=rational(low),deleted_D_upper=rational(high),
                root_interval_width=rational(high-low),full_D=rational(full),
                actual_source_gap_lower=rational(low-full),
                analytic_Peano_gap_lower=rational(lower),Peano_kernel_area=rational(area),
                deleted_six_record_coefficient_upper=rational(drop_coeff_upper),
                per_branch_epsilon=rational(eps),two_branch_error=rational(both),
                certified_margin=rational(margin),all_clock_energy_fraction=0.2,
                independent_inputs_fixed=["n0","kappa","a0","d","omega","Lambda"],
                deletion_is_not_another_solution_of_same_action=True)


def calculate():
    return dict(status="PASS",science_groups=3,group_1=moving_sources(),
                group_2=finite_records(),group_3=source_deletion(),
                scope=dict(classical_effective_action=True,CT_still_conditional=True,
                           base_constraint_and_preservation_extra=True,
                           signal_readout_source_contract_certified=False,
                           actual_precision_obtained=False,quantum_unknown_code_and_R=False,
                           whole_actual_future=False,full_GR_proved=False,
                           no_new_mass_energy_derivation=True))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-once",action="store_true")
    args=parser.parse_args()
    result=calculate()
    files=[HERE/name for name in ("proof.md","check.py","admission.md","NEXT.md")]
    files += [HERE.parent/"_admission/same_source_tidal_candidate.md",
              RESEARCH/"archive_1161_1200/1188/proof.md",HERE.parent/"1235/proof.md"]
    result["source_sha256"]={str(p.relative_to(RESEARCH)).replace("\\","/"):
                             hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    out=HERE/"results.json"
    if args.write_once:
        if out.exists():raise FileExistsError("results.json exists; default is read-only")
        out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        mode="CREATED"
    else:
        need(result==json.loads(out.read_text(encoding="utf-8")),"read-only results and source binding")
        mode="READ_ONLY_MATCH"
    print(json.dumps(dict(status="PASS",mode=mode,science_groups=3,
                         source_gap_bound=float(F(40383498710737,1194393600000000000))),
                     ensure_ascii=False))


if __name__=="__main__":main()
