"""983: inherited SM anomaly data -> ONE conserved FLRW effective source.

The low-branch comparison is within the declared free-conformal, b=0
semiclassical equations. It is NOT a full SM, quantum-geometry, or UV bound.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse,hashlib,json,math
import numpy as np

HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'sm_curvature_radiation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rational(x):return dict(exact=str(x),decimal=float(x))
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-7,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)
def integrate(f,a,b,n):
    x,w=np.polynomial.legendre.leggauss(n)
    return float((b-a)/2*np.dot(w,f((b-a)/2*x+(a+b)/2)))

def run():
    # P981 table: per generation q_L,u_R,d_R,l_L,e_R; no sterile Weyl field.
    per_generation={'q_L':6,'u_R':3,'d_R':3,'l_L':2,'e_R':1}
    nw=3*sum(per_generation.values());ns=4;nv=8+3+1
    aa=(F(ns)+F(11,2)*nw+62*nv)/360
    cc=F(ns+3*nw+12*nv,120);gstar=F(ns+2*nv)+F(7,8)*2*nw
    assert (ns,nw,nv)==(4,45,12)
    assert aa==F(1991,720) and cc==F(283,120) and gstar==F(427,4)
    # Independent conversion of 1804.02020 (95)-(97) at xi=1/6.
    xi=F(1,6)
    alpha1=2*xi*xi-F(2,3)*xi-F(277,144)
    alpha2=F(571,90);alpha3=-F(293,720)
    bR=alpha1+(alpha2+alpha3)/3
    bW=alpha2/2+2*alpha3;bE=-alpha2/2-alpha3
    assert (bR,bW,bE)==(0,cc,-aa)
    # Old extended matter counts are NOT the chosen parent.
    olda=(F(5)+F(11,2)*48+62*12)/360
    assert olda-aa==F(7,144)

    # Coefficient algebra: trace = -24 A (H^4+H^2 Hdot).
    # Conservation plus this trace fixes rho_A=6 A H^4+C/a^4.
    rho_coeff=F(6);p_h4=F(-6);p_h2hdot=F(-8)
    assert -rho_coeff+3*p_h4==-24 and 3*p_h2hdot==-24
    assert 4*rho_coeff+3*p_h2hdot==0 and 3*(rho_coeff+p_h4)==0
    # On zero-order radiation Hdot=-2 H^2: p_A=10 A H0^4.
    pred_coeff=p_h4-2*p_h2hdot
    assert pred_coeff==10 and pred_coeff/rho_coeff==F(5,3)
    assert -8*rho_coeff+3*(rho_coeff+pred_coeff)==0
    wrong_rad_pressure=(-8*rho_coeff+3*(rho_coeff+rho_coeff/3))/rho_coeff
    missing_pressure=(-8*rho_coeff+3*rho_coeff)/rho_coeff
    assert wrong_rad_pressure==-4 and missing_pressure==-5

    # Universal bound for 0<=z<=1/100, not inferred from sample points.
    zmax=F(1,100);xbound=F(51,50)
    assert xbound-1-zmax*xbound*xbound>0 and xbound<1/(2*zmax)
    assert xbound*xbound*(xbound+1)<3
    # x_low-(1+z) = z^2*x_low^2*(x_low+1), so it is <3 z^2.
    # Integral of 3 z_*^2 a^-7/2 from 1 to a_f is <=z_*^2/4.
    zstar=F(1,10000);af=F(2)
    time_bound=zstar*zstar*(1-af**-6)/4
    A=float(aa)/(16*math.pi**2)
    Hstar=1.;M2=2*A/float(zstar);rho_star=3*M2
    samples=[];max_residual=0.
    for scale in np.linspace(1,2,33):
        y0=scale**-4;z=float(zstar)*scale**-4
        y=y0*(1+z);H=math.sqrt(y);Hdot=-2*y0-4*y0*z
        rho0=rho_star*scale**-4;p0=rho0/3
        rho1=6*A*y0*y0;p1=10*A*y0*y0
        R=6*(Hdot+2*y);E4=24*y*(y+Hdot)
        friedmann=abs(3*M2*y-rho0-rho1)/rho0
        acceleration=abs(-M2*(2*Hdot+3*y)-p0-p1)/rho0
        trace_red=abs(-M2*R-(-rho1+3*p1))/rho0
        conservation=abs(-4*H*rho0-8*H*rho1+3*H*(rho0+rho1+p0+p1))/(H*rho0)
        max_residual=max(max_residual,friedmann,acceleration,trace_red,conservation)
        # Residual against the un-reduced b=0 anomaly equation; still not all-order physics.
        friedmann_full=(3*M2*y-rho0-6*A*y*y)/rho0
        trace_full=(-M2*R+A*E4)/rho0
        assert abs(friedmann_full-(-2*z*z-z**3))<2e-15
        assert abs(trace_full-(-16*z*z-12*z**3))<1e-14
        if scale in (1.,1.5,2.):
            samples.append(dict(scale=scale,H2_reduced=y,R_reduced=R,
                rho_radiation=rho0,rho_anomaly=rho1,pressure_anomaly=p1,
                fractional_density_correction=rho1/rho0,
                unreduced_Friedmann_residual_over_rho0=friedmann_full,
                unreduced_trace_residual_over_rho0=trace_full))
    assert max_residual<2e-14
    zs=float(zstar)
    xlow=lambda z:2/(1+np.sqrt(1-4*z))
    tau_low=[integrate(lambda s:s/np.sqrt(xlow(zs*s**-4)),1,2,n) for n in (48,80)]
    tau_red=.5*(math.sqrt(16+zs)-math.sqrt(1+zs))
    tau_first=1.5+zs/4*(.25-1)
    assert 0<tau_red-tau_low[1]<float(time_bound)
    assert abs(tau_low[1]-tau_low[0])<5e-14
    assert 0<tau_red-tau_first<zs*zs/16
    with localcontext() as ctx:
        ctx.prec=65;z=Decimal(1)/10000
        low=2/(1+(1-4*z).sqrt());high=(1+(1-4*z).sqrt())/(2*z)
        # z*x_high is about 1: the high root is not in the perturbative window.
        assert z*high>Decimal('.98')
        assert 0<low-(1+z)<3*z*z
        branch=dict(z0=str(z),x_low=str(low),x_high=str(high),
            low_minus_first_order=str(low-1-z),z_times_high=str(z*high))
    sources=[RESEARCH_PATH for RESEARCH_PATH in (
        STAGE.parent/'archive_531_553/research_note_553.md',
        STAGE.parent/'archive_585_628/research_note_601.md',
        STAGE.parent/'archive_629_652/research_note_631.md',
        STAGE/'research_note_964.md',STAGE/'research_note_973.md',
        STAGE/'981/drafts/common_parent_contract_v1.md',STAGE/'research_note_982.md')]
    return dict(round=983,all_scientific_checks_passed=True,
        inherited_matter=dict(per_generation_Weyl=per_generation,real_scalars=ns,
            Weyl_components=nw,gauge_vectors=nv,a=rational(aa),c=rational(cc),
            thermal_gstar=rational(gstar),beta_curvature_times_16pi2=dict(bR=str(bR),bW=str(bW),bE=str(bE)),
            old_extended_parent_a_difference=rational(olda-aa)),
        analytic=dict(anomaly_density_coefficient=6,anomaly_pressure_on_radiation_coefficient=10,
            correction_equation_of_state=rational(F(5,3)),
            wrong_radiation_pressure_conservation_coefficient=str(wrong_rad_pressure),
            missing_pressure_conservation_coefficient=str(missing_pressure),
            universal_z_max=str(zmax),low_root_upper=str(xbound),
            low_root_remainder_upper_coefficient=3,
            proper_time_error_upper=rational(time_bound),
            Friedmann_residual_polynomial=['0','0','-2','-1'],
            trace_residual_polynomial=['0','0','-16','-12']),
        benchmark=dict(zstar=str(zstar),Hstar=Hstar,M_squared=M2,A_anomaly=A,
            finite_scale_interval=[1,2],shared_equations_max_residual=max_residual,
            samples=samples,dimensionless_proper_time_low_branch=tau_low,
            dimensionless_proper_time_reduced=tau_red,
            dimensionless_proper_time_first_order=tau_first,
            dimensionless_proper_time_classical=1.5,
            reduced_vs_low_time_difference=tau_red-tau_low[1],branches=branch),
        scope=dict(inherited_coefficients_are_not_new_discoveries=True,
            same_source_pressure_trace_and_geometry=True,
            free_massless_conformal_SM_order_only=True,
            boundary_and_state_constant_specified=True,
            b_zero_scheme_used_for_resummed_comparison=True,
            universal_all_scheme_first_order_radiation_under_EFT_power_counting=True,
            full_interacting_SM_thermal_state=False,
            all_quantum_source_to_classical_geometry_error=False,
            full_physical_truncation_error_certified=False,
            high_curvature_root_adopted=False,cosmic_arrow_generated=False,
            full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sources})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes',)},ensure_ascii=False,indent=2))
