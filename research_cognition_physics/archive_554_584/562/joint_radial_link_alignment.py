"""562: full-source alignment cost when a covariant edge is used as a radial readout.

The all-state lower bound is analytic. Numerical integrals check normal physical
sources of the same H, not frozen matter or a truncated evolution. A lattice
spacing used to compare readouts is an additional interface input.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_gauge_link_reference as link
import joint_finite_time_gauge_probe as inherited

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_radial_link_alignment_results.json'
A=.5; B=1.; K=.8; RHO=.6; RMAX=1.4


@lru_cache(None)
def legendre(n):
    return np.polynomial.legendre.leggauss(n)


def angular(lam,n=320):
    x,w=legendre(n);chi=math.pi*(x+1)/2;z=np.cos(chi)
    base=w*np.sin(chi)**2*np.exp(2*lam*(z-1))
    prob=base/base.sum()
    mean=float(prob@z);sin2=float(prob@(1-z*z))
    deficit=float(prob@(1-z))
    unit_kinetic=lam*lam*sin2
    by_laplacian=float(prob@(3*lam*z-lam*lam*(1-z*z)))
    assert abs(unit_kinetic-by_laplacian)<2e-9*max(1.,unit_kinetic)
    assert abs(unit_kinetic-1.5*lam*mean)<2e-9*max(1.,unit_kinetic)
    return dict(lam=lam,m=mean,one_minus_z2=sin2,deficit=deficit,
        sphere_kinetic=unit_kinetic,kinetic_from_laplacian=by_laplacian,
        electric_energy=B*unit_kinetic/4,
        electric_energy_times_deficit=B*unit_kinetic*deficit/4,
        lambda_times_deficit=lam*deficit)


def radial(n=240):
    x,w=legendre(n);r=1+.4*x
    f=np.exp(-1/(1-x*x))
    weights=.4*w*r**3*f*f;weights/=weights.sum()
    lp=-2*x/(.4*(1-x*x)**2)
    lpp=-2/(.4**2*(1-x*x)**2)-8*x*x/(.4**2*(1-x*x)**3)
    grad=float(weights@(lp*lp))
    lap=float(weights@(-lpp-lp*lp-3*lp/r))
    assert abs(grad-lap)<1e-10
    return dict(mean=float(weights@r),square=float(weights@(r*r)),
        fourth=float(weights@(r**4)),inverse_square=float(weights@(1/(r*r))),
        radial_gradient_energy_without_a=grad,radial_laplacian_energy_without_a=lap)


def invariant_checks():
    rng=np.random.default_rng(562);errors=dict(invariance=0.,decomposition=0.,
        group_gradient=0.,endpoint_gradient=0.,group_laplacian=0.)
    eps=2e-4
    for _ in range(48):
        x,y=rng.normal(size=(2,4));x/=np.linalg.norm(x);y/=np.linalg.norm(y)
        q=rng.normal(size=(3,4));q/=np.linalg.norm(q,axis=1)[:,None]
        U,G,H=map(link.left,q);z=float(x@U@y)
        moved=float((G@x)@(G@U@H.T)@(H@y))
        errors['invariance']=max(errors['invariance'],abs(moved-z))
        r1,r2=rng.uniform(RHO,RMAX,size=2)
        F=K/2*np.sum((r1*x-U@(r2*y))**2)
        decoded=K*((r1-r2)**2/2+r1*r2*(1-z))
        errors['decomposition']=max(errors['decomposition'],abs(F-decoded))
        gz=np.array([x@j@U@y for j in link.J])
        errors['group_gradient']=max(errors['group_gradient'],abs(gz@gz-(1-z*z)/4))
        gx=U@y-z*x;gy=U.T@x-z*y
        errors['endpoint_gradient']=max(errors['endpoint_gradient'],
            abs(gx@gx-(1-z*z)),abs(gy@gy-(1-z*z)))
        lap=0.
        for j in link.J:
            plus=math.cos(eps/2)*np.eye(4)+2*math.sin(eps/2)*j
            minus=plus.T
            lap+=(x@plus@U@y+x@minus@U@y-2*z)/eps**2
        errors['group_laplacian']=max(errors['group_laplacian'],abs(lap+.75*z))
    assert max(v for k,v in errors.items() if k!='group_laplacian')<2e-14
    assert errors['group_laplacian']<8e-8
    return errors


def noncentral_sources():
    # Exact Gauss invariant relative quaternion w = conj(X1/r1) U (X2/r2).
    # Its Haar marginal is S3. These are not all class-function sources.
    xx,ww=legendre(96);chi=math.pi*(xx+1)/2
    zz,vv=legendre(16);phi=2*math.pi*np.arange(24)/24
    z=np.broadcast_to(np.cos(chi)[:,None,None],(96,16,24)).reshape(-1)
    sin=np.sin(chi)[:,None,None]
    yy=np.sqrt(1-zz**2)[None,:,None]
    w1=np.broadcast_to(sin*yy*np.cos(phi)[None,None,:],(96,16,24)).reshape(-1)
    w2=np.broadcast_to(sin*yy*np.sin(phi)[None,None,:],(96,16,24)).reshape(-1)
    w3=np.broadcast_to(sin*zz[None,:,None],(96,16,24)).reshape(-1)
    pts=np.array([z,w1,w2,w3]).T
    weights=np.broadcast_to((ww*np.sin(chi)**2)[:,None,None]*vv[None,:,None]/48,
                            (96,16,24)).reshape(-1)
    rows=[]
    for lam,eta,xi in ((0.,1.,.7),(.4,1.1,.6),(2.,.3,-.8),(6.,1.5,.5)):
        fac=np.exp(lam*(z-1));P=1+1j*eta*w1+xi*w2
        wave=fac*P
        ambient=np.zeros_like(pts,dtype=complex)
        ambient[:,0]=lam*wave;ambient[:,1]=1j*eta*fac;ambient[:,2]=xi*fac
        tangent=ambient-pts*np.sum(pts*ambient,axis=1)[:,None]
        norm=float(weights@(abs(wave)**2))
        m=float(weights@(z*abs(wave)**2))/norm
        sin2=float(weights@((1-z*z)*abs(wave)**2))/norm
        kinetic=float(weights@np.sum(abs(tangent)**2,axis=1))/norm
        lower=9*m*m/(4*sin2)
        assert kinetic+1e-12>=lower
        rows.append(dict(lam=lam,eta=eta,xi=xi,m=m,sphere_kinetic=kinetic,
            lower_before_Jensen=lower,lower_after_Jensen=9*m*m/(4*(1-m*m))))
    return rows


def source_budget(lam,n=320,nr=240):
    ag=angular(lam,n);rr=radial(nr);L,u,c=inherited.constants()
    s2=.5;s4=.75
    W=.5*(L[0,0]*(rr['fourth']-2*u[0]*rr['square']+u[0]**2)
        +2*L[0,1]*(rr['square']-u[0])*(s2-u[1])
        +L[1,1]*(s4-2*u[1]*s2+u[1]**2))
    radial_energy=2*A*rr['radial_gradient_energy_without_a']
    singlet_energy=A
    angular_coefficient=2*A*rr['inverse_square']+B/4
    angular_energy=angular_coefficient*ag['sphere_kinetic']
    edge_radial=K*(rr['square']-rr['mean']**2)
    edge_angular=K*rr['mean']**2*ag['deficit']
    E=radial_energy+singlet_energy+angular_energy+W+edge_radial+edge_angular
    full_B=B+8*A/RMAX**2
    lower=9*full_B/16*ag['m']**2/(1-ag['m']**2)
    assert E>=lower and angular_energy>=lower
    return dict(**ag,radial_moments=rr,radial_energy=radial_energy,
        singlet_energy=singlet_energy,angular_coefficient=angular_coefficient,
        angular_energy=angular_energy,potential_W=float(W),radial_edge_energy=edge_radial,
        angular_edge_energy=edge_angular,full_source_energy=float(E),
        universal_annulus_energy_lower=float(lower),
        deficit_lower_using_full_budget=1-math.sqrt(16*E/(16*E+9*full_B)),
        normalized_angular_energy_deficit=angular_energy*ag['deficit']/angular_coefficient)


def complete_H_stencil():
    L,u,_=inherited.constants();rng=np.random.default_rng(1562)
    lam=1.2;errors=[]
    def wave(q,U):
        r1=np.linalg.norm(q[:4]);r2=np.linalg.norm(q[5:9])
        z=q[:4]@U@q[5:9]/(r1*r2)
        t1=(r1-1)/.4;t2=(r2-1)/.4
        if abs(t1)>=1 or abs(t2)>=1:return 0.
        return math.exp(-1/(1-t1*t1)-1/(1-t2*t2)
            -(q[4]**2+q[9]**2)/2+lam*(z-1))
    for _ in range(12):
        q=rng.normal(size=10)
        for o in (0,5):q[o:o+4]*=rng.uniform(.92,1.08)/np.linalg.norm(q[o:o+4])
        q[[4,9]]*=.5
        quat=rng.normal(size=4);quat/=np.linalg.norm(quat);U=link.left(quat)
        r1=np.linalg.norm(q[:4]);r2=np.linalg.norm(q[5:9]);z=q[:4]@U@q[5:9]/(r1*r2)
        kinetic=0.;W=0.
        for o,r in ((0,r1),(5,r2)):
            t=(r-1)/.4;lp=-2*t/(.4*(1-t*t)**2)
            lpp=-2/(.4**2*(1-t*t)**2)-8*t*t/(.4**2*(1-t*t)**3)
            kinetic+=-A*(lpp+lp*lp+3*lp/r)+A*(1-q[o+4]**2)
            d=np.array([r*r,q[o+4]**2])-u;W+=d@L@d/4
        coef=A*(1/r1**2+1/r2**2)+B/4
        kinetic+=coef*(3*lam*z-lam*lam*(1-z*z))
        F=K/2*np.sum((q[:4]-U@q[5:9])**2)
        target=kinetic+W+F;h=1e-4;psi=wave(q,U);lap=glap=0.
        for j in range(10):
            shift=np.zeros(10);shift[j]=h
            lap+=(wave(q+shift,U)+wave(q-shift,U)-2*psi)/h**2
        for J in link.J:
            plus=math.cos(h/2)*np.eye(4)+2*math.sin(h/2)*J
            glap+=(wave(q,plus@U)+wave(q,plus.T@U)-2*psi)/h**2
        actual=(-A*lap-B*glap)/psi+W+F
        errors.append(abs(actual-target)/max(1.,abs(target)))
    assert max(errors)<2e-6
    return dict(samples=12,step=1e-4,maximum_relative_residual=max(errors),
        all_ten_material_and_three_group_derivatives_included=True)


def scaling():
    full_B=B+8*A/RMAX**2;eta=.1;rows=[]
    for eps in (.5,.25,.125,.0625):
        x=eta*eps**2/RHO**2
        necessary=9*full_B/16*(1-x)**2/(x*(2-x))
        electric_only=9*B/16*(1-x)**2/(x*(2-x))
        rows.append(dict(epsilon=eps,allowed_radial_readout_bias=eta,
            maximum_alignment_deficit=x,necessary_source_energy=necessary,
            electric_only_necessary_energy=electric_only,
            rescaled_energy=necessary*eps**2))
    budget=100.
    d=1-math.sqrt(16*budget/(16*budget+9*full_B))
    if_B_zero=8*A/RMAX**2
    d0=1-math.sqrt(16*budget/(16*budget+9*if_B_zero))
    assert d>0 and d0>0
    return dict(rho=RHO,upper_radius=RMAX,a=A,b=B,effective_B=full_B,
        fixed_total_energy_budget=budget,deficit_lower=d,
        deficit_lower_even_with_b_zero=d0,necessary_energy_rows=rows,
        limiting_rescaled_energy=9*full_B*RHO**2/(32*eta),
        statement='Necessary for the uncorrected covariant readout to approximate radial energy; not all continuum limits.')


def run():
    identity=invariant_checks();noncentral=noncentral_sources()
    stencil=complete_H_stencil()
    rows=[source_budget(lam) for lam in (.5,2.,8.,32.,128.,512.,2048.)]
    other=[source_budget(row['lam'],480,360) for row in rows]
    fields=('m','deficit','full_source_energy','angular_energy')
    residual=max(abs(r[k]-s[k])/max(1.,abs(r[k])) for r,s in zip(rows,other) for k in fields)
    assert residual<1e-10
    assert abs(rows[-1]['lambda_times_deficit']-.75)<1e-4
    assert abs(rows[-1]['electric_energy_times_deficit']-9/32)<2e-4
    assert abs(rows[-1]['normalized_angular_energy_deficit']-9/8)<6e-4
    requirements=scaling()
    checks=['same_gauge_invariant_edge_and_normalized_angular_geometry',
        'general_angular_energy_lower_bound_with_noncentral_complex_physical_sources',
        'smooth_annular_full_matter_source_and_radial_integration_by_parts',
        'complete_H_differential_action_with_both_endpoint_and_link_kinetics',
        'concentrating_Gauss_family_and_independent_quadrature',
        'readout_scale_energy_requirements_and_nonzero_material_floor']
    deps=('joint_gauge_link_reference.py','joint_finite_time_gauge_probe.py',
        'joint_singlet_common_mass_rg_results.json','research_note_548.md','research_note_556.md',
        'research_note_558.md','joint_finite_duration_clock_readout_results.json')
    return dict(round=562,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        angular_identity_checks=identity,noncentral_source_checks=noncentral,
        full_H_stencil=stencil,normal_Gauss_source_family=rows,
        independent_quadrature_relative_residual=residual,scaling_requirements=requirements,
        scope=dict(all_state_alignment_lower_bound_is_analytic=True,
            mixed_correlated_sources_and_passive_references_allowed=True,
            finite_amplitude_annulus_and_fixed_H_budget_are_explicit_inputs=True,
            concentrating_states_are_full_normal_sources_not_frozen_matter=True,
            only_uncorrected_covariant_to_radial_readout_identification_is_limited=True,
            separate_radius_readouts_and_angle_subtraction_not_excluded=True,
            epsilon_is_added_interface_scale_not_derived_space=True,
            no_radial_smooth_continuum_from_alignment_alone=True,
            no_three_dimensions_GR_or_unified_completion=True,
            numerical_integrals_not_machine_interval_certificates=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=562,tests=result['tests_run'],stencil=result['full_H_stencil'],
        quadrature=result['independent_quadrature_relative_residual'],
        last_family=result['normal_Gauss_source_family'][-1],
        scaling=result['scaling_requirements']),ensure_ascii=False))
