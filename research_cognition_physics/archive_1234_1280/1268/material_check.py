"""1268 material principal-symbol check: one author group, no PDE simulation.

Default readonly exact recomputation. --write-results exclusively creates JSON.
The finite parameter samples are diagnostics, not the actual tiny-epsilon task.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import sqrt, exp
from pathlib import Path

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
RESULT=HERE/"material_results.json"

def fs(x):
    x=F(x)
    return str(x.numerator) if x.denominator==1 else f"{x.numerator}/{x.denominator}"

def sha(path):
    return sha256(path.read_bytes()).hexdigest()

def jread(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def principal_sample(kind,a,c):
    gamma=F(1)
    C0=F(3,10)             # original 3 k^2
    s2=F(1,10)/a**2       # local physical rod gradient squared
    B=3*s2
    n=a**-3
    if kind=="sqrt":
        f=a
        fb=-f/(2*B)
        fbb=3*f/(4*B**2)
    elif kind=="linear":
        f=B/C0
        fb=1/C0
        fbb=F(0)
    else:
        raise ValueError(kind)
    mbase=F(1)
    mbar=mbase+c*f
    A=gamma+2*n*c*fb
    D=A+4*n*c*fbb*s2
    K=n*mbar
    C2=4*n*n*c*c*fb*fb*s2
    assert min(A,D,K,mbar)>0 and C2>0

    # Derive the principal coefficients from the actual quadratic expansion.
    # phi_long=s*u. Then all coefficients including the mixed Cr are rational.
    def L2(ut,ux,xt,xx):
        dB1=2*s2*ux
        dB2=s2*(ux*ux-ut*ut)
        field=gamma*s2*(ut*ut-ux*ux)/2
        dust=K*xt*xt/2-n*c*(fb*dB2+fbb*dB1*dB1/2)
        advection_after_parts=2*n*c*fb*s2*xx*ux
        return field+dust+advection_after_parts
    Ar=2*L2(F(1),F(0),F(0),F(0))
    Dr=-2*L2(F(0),F(1),F(0),F(0))
    Kr=2*L2(F(0),F(0),F(1),F(0))
    Cr=L2(F(0),F(1),F(0),F(1))-L2(F(0),F(1),F(0),F(0))
    assert Ar==A*s2 and Dr==D*s2 and Kr==K
    assert Cr==2*n*c*fb*s2 and Cr*Cr==C2*s2
    for y in (F(-1),F(0),F(1),F(2),F(3,7)):
        determinant=(Dr-Ar*y)*(-K*y)-Cr*Cr
        polynomial=s2*(A*K*y*y-D*K*y-C2)
        assert determinant==polynomial
    # Continuous root signs are exact: leading positive and constant negative.
    poly0=-C2
    poly1=A*K-D*K-C2
    assert poly0<0
    if kind=="sqrt":
        assert A==gamma-n*c*f/B and D==gamma
        assert (D/A>1)==(c>0)
        if c>0:
            assert poly1<0
    else:
        assert D==A and poly1==-C2<0  # full positive root >1

    # Stable root display; analytic bounds use exact C2/(K D).
    Af,Df,Kf,C2f=map(float,(A,D,K,C2))
    disc=sqrt(Df*Df+4*Af*C2f/Kf)
    yplus=(Df+disc)/(2*Af)
    yminus=-2*C2f/(Kf*(disc+Df))
    assert yminus<0<yplus
    growth_squared_upper=C2/(K*D)
    qmax=F(10)
    T=F(1,10)
    rate_time_squared_upper=qmax*qmax*T*T*growth_squared_upper
    assert rate_time_squared_upper<F(1,10000)
    gamma_times_T=float(qmax*T)*sqrt(-yminus)
    assert gamma_times_T<0.01
    assert exp(gamma_times_T)<100/99
    return {
        "kind":kind, "a":fs(a), "c":fs(c), "gamma":fs(gamma),
        "Bphi":fs(B), "n":fs(n), "mbar":fs(mbar),
        "f":fs(f), "f_B":fs(fb), "f_BB":fs(fbb),
        "A":fs(A), "D":fs(D), "K":fs(K), "C_squared":fs(C2),
        "quadratic_action_coefficients":{
            "A_scaled":fs(Ar),"D_scaled":fs(Dr),
            "K":fs(Kr),"C_scaled":fs(Cr)},
        "polynomial_at_zero":fs(poly0),"polynomial_at_one":fs(poly1),
        "positive_root_display":f"{yplus:.14g}",
        "negative_root_display":f"{yminus:.14g}",
        "finite_frozen_example":{
            "qmax":fs(qmax),"T":fs(T),
            "Gamma_T_squared_upper":fs(rate_time_squared_upper),
            "Gamma_T_display":f"{gamma_times_T:.12g}",
            "mode_exponential_growth_bound":"<100/99",
            "scope":"one constant-coefficient eigenmode, not full evolution norm"}
    }

def equal_reduction_samples():
    rows=[]
    C0=F(3,10)
    for a in (F(1),F(2),F(8,3)):
        for c in (F(1,1000),F(-1,1000)):
            n=a**-3
            B=C0/a**2
            f=a
            fb=-f/(2*B)
            d=c/C0
            hmass=a**3*n*c*f
            hvol=a**3*d*B
            rho_mass=n*c*f
            rho_vol=d*B
            # 2 n c f_B S_ij, isotropic S_ij=B delta_ij/3.
            pressure_mass=2*n*c*fb*B/3
            pressure_vol=-d*B/3
            assert hmass==hvol==c*a
            assert rho_mass==rho_vol==c/a**2
            assert pressure_mass==pressure_vol==-c/(3*a**2)
            shear_mass=2*n*c*fb
            shear_vol=2*d
            assert shear_mass==-d and shear_vol==2*d
            assert shear_mass-shear_vol==-3*d
            assert c*f!=0  # density derivative differs from zero volume response
            Gamma=1+2*d
            assert Gamma>0
            rows.append({
                "a":fs(a),"c":fs(c),"Bphi":fs(B),"n":fs(n),
                "common_h":fs(hmass),"common_rho":fs(rho_mass),
                "common_pressure":fs(pressure_mass),
                "mass_density_derivative":fs(c*f),
                "volume_density_derivative":"0",
                "mass_traceless_shear_coefficient":fs(shear_mass),
                "volume_traceless_shear_coefficient":fs(shear_vol),
                "volume_gamma_effective":fs(Gamma)})
    # Exact pulse coefficient in positivity criterion (original k^2=1/10).
    assert 2*F(35,16)/C0==F(175,12)
    return rows

def result():
    sources={
        "1268/admission.md":HERE/"admission.md",
        "1268/material_proof.md":HERE/"material_proof.md",
        "1268/material_check.py":Path(__file__).resolve(),
        "1245/proof.md":STAGE/"1245/proof.md",
        "1247/proof.md":STAGE/"1247/proof.md",
        "1249/canonical.md":STAGE/"1249/canonical.md",
        "1267/geometry_proof.md":STAGE/"1267/geometry_proof.md",
        "_admission/minimal_common_task_after1267/candidate_screen.md":
            STAGE/"_admission/minimal_common_task_after1267/candidate_screen.md",
        "_admission/minimal_common_task_after1267/independent_review.md":
            STAGE/"_admission/minimal_common_task_after1267/independent_review.md",
        "_admission/local_density_reader_scope_review_after1251.md":
            STAGE/"_admission/local_density_reader_scope_review_after1251.md"}
    principal=[]
    for kind in ("sqrt","linear"):
        for a in (F(1),F(2),F(8,3)):
            for c in (F(1,1000),F(-1,1000)):
                principal.append(principal_sample(kind,a,c))
    return {
        "round":1268,"part":"material_principal_and_reduction",
        "status":"PASS_AUTHOR_WITH_SCOPE","author_groups":1,
        "arithmetic":"Fraction identities/signs; stable floating roots display only",
        "principal_samples":principal,
        "equal_reduction_and_distinct_response":equal_reduction_samples(),
        "analytic_certificates":{
            "dynamic_dust_polynomial":"A*K*y^2-D*K*y-C^2",
            "opposite_root_signs":"A,K>0 and C!=0",
            "linear_B_mass_positive_root_above_one":"D=A and polynomial(1)=-C^2",
            "finite_frozen_growth_rate_bound":"Gamma(q)<=|q|*|C|/sqrt(K*D), D>0",
            "homogeneous_reduction":"h_mass=h_volume=c*a",
            "homogeneous_sources":"rho=c/a^2, P=-c/(3*a^2)",
            "traceless_shear_difference":"-3*d deltaS",
            "volume_positivity_pulse_coefficient":"175/12"
        },
        "scope":[
            "explicit homogeneous comoving jet and longitudinal perturbation block",
            "dynamic pressureless dust retained; no artificial restoring pressure",
            "untruncated principal obstruction is not a finite-band task refutation",
            "finite frozen mode growth is not complete time-dependent propagation bound",
            "volume coefficient changes material extension; full theta/P identity not certified",
            "same reduced readout and homogeneous stress do not select spatial extension",
            "finite samples are not original tiny-epsilon preparation/clock/device parameters",
            "no new PDE evolution, QFT, physical causal implementation or GR completion"
        ],
        "source_bindings":{name:sha(path) for name,path in sources.items()}
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--write-results",action="store_true")
    args=p.parse_args()
    current=result()
    if args.write_results:
        with RESULT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(current,f,ensure_ascii=False,indent=2,sort_keys=True)
            f.write("\n")
        print("PASS_AUTHOR_WITH_SCOPE FIRST_WRITE_EXCLUSIVE")
    else:
        saved=jread(RESULT)
        if saved!=current:
            keys=[k for k in set(current)|set(saved) if current.get(k)!=saved.get(k)]
            raise SystemExit("MISMATCH: "+", ".join(sorted(keys)))
        print("PASS_AUTHOR_WITH_SCOPE READ_ONLY_MATCH")
    print("author_groups=1")
    print("dynamic mass-reader negative longitudinal branch: exact polynomial sign")
    print("equal homogeneous h/rho/P, different density and shear response")
    print("finite frozen example growth<100/99; not an actual task certificate")

if __name__=="__main__":
    main()
