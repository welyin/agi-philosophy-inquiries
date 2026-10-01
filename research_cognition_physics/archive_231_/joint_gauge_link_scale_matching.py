"""558: common fast-link energy/observable matching and an instrument limit.

The full-source family is an initial-state construction, not a proof of
Born-Oppenheimer dynamics. All model inputs of 557 remain explicit.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import joint_matter_energy_moment_control as matter

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauge_link_scale_matching_results.json'


def spectral(alpha,b,N=24):
    n=np.arange(N+1,dtype=float); d=n*(n+2)/4
    C=np.diag(np.ones(N)/2,1)+np.diag(np.ones(N)/2,-1)
    K=b*np.diag(d)-alpha*C
    ev,U=np.linalg.eigh(K); psi=U[:,0]
    if psi[0]<0: psi=-psi
    z=float(psi@C@psi)
    z2=float(np.linalg.norm(C@psi)**2+psi[-1]**2/4)
    electric=float(b*np.sum(d*psi*psi))
    # Normalized real eigenvector derivative, in the fixed character basis.
    dc=U[:,1:].T@C@psi
    derivative=U[:,1:]@(dc/(ev[1:]-ev[0]))
    delta=3*b/4
    tail_floor=b*(N+1)*(N+3)/4-alpha
    assert tail_floor>0
    adjusted=K.copy(); adjusted[-1,-1]-=alpha*alpha/(4*tail_floor)
    lower=float(min(0,np.linalg.eigvalsh(adjusted)[0]))
    mu=float(ev[0])
    assert lower<=mu+1e-12
    out=dict(alpha=alpha,b=b,lambda_fast=mu,cosine_mean=z,cosine_second=z2,
        edge_shift=-alpha*z,electric_energy=electric,
        derivative_norm=float(derivative@derivative),
        outside_constant_weight=float(max(0,1-psi[0]**2)),
        last_retained_coefficient=float(psi[-1]),
        truncation_lower_bound_floating=lower,truncation_upper_bound_floating=mu,
        exact_tail_floor=tail_floor)
    if alpha<delta:
        out['second_order_error_bound']=alpha**4/(16*delta*(delta-alpha))*(1/(delta-alpha)+1/(2*b-alpha))
        s=alpha/(2*(delta-alpha)); q=s*s/(1+s*s)
        out['constant_weight_loss_bound']=q
        assert out['outside_constant_weight']<=q+1e-12
        assert abs(z2-.25)<=math.sqrt(q)+1e-12
        assert abs(mu+alpha*alpha/(3*b))<=out['second_order_error_bound']+1e-13
    return out,psi,derivative


def character_check():
    nodes,w=np.polynomial.legendre.leggauss(160)
    theta=np.pi*(nodes+1)/2; haar=w*np.sin(theta)**2
    n=np.arange(9)[:,None]+1
    chars=np.sin(n*theta)/np.sin(theta)
    gram=(chars*haar)@chars.T
    C=(chars*(haar*np.cos(theta)))@chars.T
    wanted=np.diag(np.ones(8)/2,1)+np.diag(np.ones(8)/2,-1)
    err=float(max(np.max(abs(gram-np.eye(9))),np.max(abs(C-wanted))))
    assert err<1e-12
    # The top retained character has genuine C^2 moment 1/2, not 1/4.
    actual=float(haar@(np.cos(theta)**2*chars[-1]**2))
    naive=float(np.sum(wanted[-1]**2))
    assert abs(actual-.5)<1e-12 and naive==.25
    return dict(orthogonality_and_multiplication_residual=err,
                top_character_actual_second=actual,naive_truncated_square=naive)


def independent_reading(alpha,b,nu,N=24):
    row,psi,_=spectral(alpha,b,N)
    x,w=np.polynomial.legendre.leggauss(180)
    theta=np.pi*(x+1)/2; weights=w*np.sin(theta)**2
    n=np.arange(N+1)[:,None]+1
    basis=np.sin(n*theta)/np.sin(theta)
    der=(n*np.cos(n*theta)*np.sin(theta)-np.sin(n*theta)*np.cos(theta))/np.sin(theta)**2
    wave=psi@basis; wave_der=psi@der
    density=weights*wave**2; z=np.cos(theta)
    norm=float(np.sum(density)); cm=float(density@z); c2=float(density@(z*z))
    el=float(b/4*(weights@(wave_der*wave_der)))
    assert abs(norm-1)<1e-12
    residual=max(abs(cm-row['cosine_mean']),abs(c2-row['cosine_second']),abs(el-row['electric_energy']))
    assert residual<1e-10
    # Equal endpoint radii make F0=alpha. This affects the read mean but not
    # the electric injection, which depends only on the U derivative.
    F=alpha*(1-z)
    gh,gw=np.polynomial.hermite.hermgauss(32); gw=gw/math.sqrt(math.pi)
    noise=math.sqrt(2)*nu*gh
    Y=F[:,None]+noise
    mean=float(density@(Y@gw))
    second=float(density@((Y*Y)@gw))
    # Integrate the actual conditional wavefunction derivatives:
    # d(K_y psi)/K_y = psi' + psi*(y-F)*F'/(2 nu^2).
    changed=wave_der[:,None]+wave[:,None]*noise[None,:]*(alpha*np.sin(theta))[:,None]/(2*nu*nu)
    after=float(b/4*(weights@((changed*changed)@gw)))
    injection=after-el
    predicted=b*alpha*alpha*(1-row['cosine_second'])/(16*nu*nu)
    residual=max(residual,abs(injection-predicted))
    variance=second-mean*mean
    target=alpha*alpha*(row['cosine_second']-row['cosine_mean']**2)+nu*nu
    assert abs(variance-target)<1e-12 and abs(injection-predicted)<1e-10
    if 'constant_weight_loss_bound' in row:
        lower=b*alpha*alpha*max(0,.75-math.sqrt(row['constant_weight_loss_bound']))/(16*nu*nu)
        assert injection>=lower-1e-10
    else: lower=None
    row.update(noise_variance=nu*nu,actual_record_mean=mean,
        actual_record_variance=variance,naive_one_band_record_variance=nu*nu,
        instantaneous_electric_injection=injection,injection_per_b=injection/b,
        asymptotic_injection_per_b=3*alpha*alpha/(64*nu*nu),
        explicit_injection_lower_bound=lower,independent_probability_and_energy_residual=residual)
    return row


def full_source(b,quadrature=32,N=20):
    # Smooth compact radial amplitude exp[-1/(1-t^2)] at r=1+.5t;
    # each singlet has a normalized Gaussian wavefunction, omega_s=1.3.
    x,w=np.polynomial.legendre.leggauss(quadrature)
    r=1+.5*x
    amp=np.exp(-1/(1-x*x))
    rw=w*.5*r**3*amp**2; rw/=sum(rw)
    logder=-2*x/(.5*(1-x*x)**2)
    v=1.;kappa=.8;a=.5;omega=1.3;nu=.5
    L,C,u=matter.matter()
    r2=float(rw@(r*r));r4=float(rw@(r**4));s2=1/(2*omega);s4=3*s2*s2
    oneW=v/4*(L[0,0]*(r4-2*u[0]*r2+u[0]**2)
        +2*L[0,1]*(r2-u[0])*(s2-u[1])+L[1,1]*(s4-2*u[1]*s2+u[1]**2))
    phiT=a*(2*float(rw@(logder*logder))+omega)
    base=phiT+2*oneW+kappa*v*r2
    extra_matter=0.;electric=0.;edge=0.;read_cost=0.;alpha2=0.
    for i,r1 in enumerate(r):
        for j,r2local in enumerate(r):
            weight=rw[i]*rw[j];alpha=kappa*v*r1*r2local
            row,_,_=spectral(alpha,b,N)
            # Full X gradients of a normalized class-function ground state.
            extra_matter+=weight*a*((kappa*v)**2*(r1*r1+r2local*r2local)*row['derivative_norm']
                +4*(1/(r1*r1)+1/(r2local*r2local))*row['electric_energy']/b)
            electric+=weight*row['electric_energy']
            edge+=weight*row['edge_shift']
            read_cost+=weight*b*alpha*alpha*(1-row['cosine_second'])/(16*nu*nu)
            alpha2+=weight*alpha*alpha
    E1=base+extra_matter+electric+edge
    return dict(b=b,quadrature=quadrature,character_cutoff=N,
        initial_E1=float(E1),limiting_initial_energy=float(base),
        induced_matter_kinetic=extra_matter,initial_electric_energy=electric,
        edge_mean_shift=edge,instantaneous_electric_injection=read_cost,
        injection_per_b=read_cost/b,asymptotic_injection_per_b=3*alpha2/(64*nu*nu),
        support_radii=[.5,1.5],noise_variance=nu*nu)


def run():
    checks=[];character=character_check()
    checks.append('Haar_character_basis_and_exact_second_moment_boundary')
    spectral_rows=[]
    for alpha in (.025,.05,.1,.2,.4,1.,4.,8.):
        row,_,_=spectral(alpha,1.,24)
        other,_,_=spectral(alpha,1.,36)
        row['independent_cutoff_difference']=abs(row['lambda_fast']-other['lambda_fast'])
        assert row['independent_cutoff_difference']<1e-11
        assert row['truncation_lower_bound_floating']-1e-11<=other['lambda_fast']<=row['truncation_upper_bound_floating']+1e-11
        assert abs(row['lambda_fast']-row['edge_shift']-row['electric_energy'])<1e-11
        h=1e-4
        lp=spectral(alpha+h,1.,24)[0]['lambda_fast']
        lm=spectral(alpha-h,1.,24)[0]['lambda_fast']
        derivative=(lp-lm)/(2*h)
        assert abs(derivative+row['cosine_mean'])<1e-8
        row['independent_Hellmann_Feynman_residual']=abs(derivative+row['cosine_mean'])
        spectral_rows.append(row)
    checks.append('two_stage_Schur_error_bound_and_full_infinite_tail_bracket')
    checks.append('same_coupling_derivative_matches_edge_and_electric_observables')
    for alpha in (.005,.01,.02):
        row,_,_=spectral(alpha,1.,24)
        coefficient=(row['lambda_fast']+alpha*alpha/3)/alpha**4
        assert abs(coefficient-5/54)<3e-5
    # At equal radii, exact F0+lambda >=0. Extending second order to alpha=8
    # would instead be negative, so it is not a global effective potential.
    outside=spectral_rows[-1]
    assert 8+outside['lambda_fast']>0 and 8-64/3<0
    checks.append('fourth_order_coefficient_and_outside_window_counterexample')
    records=[independent_reading(.3,b,.5) for b in (1.,2.,10.,50.,100.)]
    assert records[-1]['actual_record_variance']>.2724
    assert abs(records[-1]['injection_per_b']-records[-1]['asymptotic_injection_per_b'])<1e-6
    checks.append('actual_Gaussian_record_and_independent_post_instrument_group_energy')
    sources=[full_source(b) for b in (4.,10.,40.,100.)]
    other=full_source(100.,40,24)
    residual=abs(sources[-1]['initial_E1']-other['initial_E1'])
    # Compact bump quadrature, not a floating-point interval certificate.
    assert residual<1e-3
    assert all(abs(row['initial_E1']-row['limiting_initial_energy'])<.2 for row in sources)
    assert sources[-1]['instantaneous_electric_injection']>20*sources[0]['instantaneous_electric_injection']
    checks.append('normalized_full_Gauss_source_family_bounded_initial_energy_and_growing_read_cost')
    deps=('joint_matter_energy_moment_control.py','joint_singlet_common_mass_rg_results.json',
          'joint_gauge_link_reference.py','joint_gauge_link_reference_results.json',
          'round558_drafts/fast_link_matching_probe.py','round558_drafts/fast_link_matching_probe_results.json')
    return dict(round=558,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        character_check=character,fast_link_spectrum=spectral_rows,
        matched_records=records,full_source_family=sources,
        full_source_quadrature_difference=residual,
        scope=dict(given_graph_SU2_and_ideal_instrument_inherited=True,
          fast_spectral_surface_not_full_slow_dynamics=True,
          full_source_family_is_initial_state_not_stationary_or_low_band_dynamics=True,
          same_action_determines_mean_observables_but_not_all_instrument_statistics=True,
          instantaneous_fixed_resolution_family_not_low_energy_closed=True,
          slow_or_adiabatic_detectors_not_excluded=True,
          analytic_tail_bounds_have_only_floating_evaluations=True,
          no_spacetime_gravity_or_unified_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=558,tests=result['tests_run'],
         full_source_energy=[r['initial_E1'] for r in result['full_source_family']],
         electric_injection=[r['instantaneous_electric_injection'] for r in result['full_source_family']],
         quadrature_difference=result['full_source_quadrature_difference']),ensure_ascii=False))

