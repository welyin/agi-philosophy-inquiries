"""561: a positive bounded-record time response under the full source H.

Full-H existence follows from the polynomial-domain Duhamel certificate.
Numerical compressed-model records are explicitly secondary diagnostics.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import joint_matter_energy_moment_control as poly
import joint_probe_time_resolution as preceding
import joint_gauge_link_reference as link

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_finite_duration_clock_readout_results.json'
spec=importlib.util.spec_from_file_location('source_probe561',HERE/'round561_drafts/short_time_source_probe.py')
source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)


def value(P,q):
    return sum(c*math.prod(x**n for x,n in zip(q,powers)) for powers,c in P.items())


def differential_check(model):
    P,H,J,F,Gamma=model['P'],model['H'],model['J'],model['F'],model['Gamma']
    HP=H(P);JP=J(P);GP=poly.mul(Gamma,P);rng=np.random.default_rng(561)
    a,b,k,eta=.5,1.,.8,.3;L,u,_=preceding.old.constants()
    theta=.4;norm=model['norm'];residuals=[]
    def wave(q,U):
        x,y=q[:4],q[5:9]
        f=k/2*np.sum((x-U@y)**2)
        return np.exp(-q@q/2)*(1+1j*eta*(f-1.6))*np.exp(-1j*theta*f)/math.sqrt(norm)
    for _ in range(8):
        q=.5*rng.normal(size=10);I=np.eye(4);base=wave(q,I);W=0.
        for o in (0,5):
            d=np.array([q[o:o+4]@q[o:o+4],q[o+4]**2])-u
            W+=.25*d@L@d
        f=k/2*np.sum((q[:4]-q[5:9])**2)
        target=np.exp(-q@q/2-1j*theta*f)*(value(HP,q)+theta*value(JP,q)+theta*theta*value(GP,q))/math.sqrt(norm)
        local_errors=[]
        for h in (4e-4,2e-4):
            lap=0j;glap=0j
            for i in range(10):
                shift=np.zeros(10);shift[i]=h
                lap+=(wave(q+shift,I)+wave(q-shift,I)-2*base)/h**2
            for gen in link.J:
                plus=math.cos(h/2)*I+2*math.sin(h/2)*gen
                minus=math.cos(h/2)*I-2*math.sin(h/2)*gen
                glap+=(wave(q,plus)+wave(q,minus)-2*base)/h**2
            actual=-a*lap-b*glap+(W+f)*base
            local_errors.append(abs(actual-target)/max(1.,abs(target)))
        residuals.append(local_errors)
    worst=max(r[1] for r in residuals)
    assert worst<1e-6
    # Independently check the terminating conjugation commutator.
    comm=poly.scale(poly.add(poly.mul(F,J(P)),poly.scale(J(poly.mul(F,P)),-1)),1j)
    diff=poly.add(comm,poly.scale(poly.mul(Gamma,P),-2))
    coeff=max([abs(c) for c in diff.values()]+[0.])
    assert coeff<1e-12
    return dict(complete_10_coordinate_and_3_group_stencil_residual=worst,
        stencil_steps=[4e-4,2e-4],each_residual=residuals,
        double_commutator_coefficient_residual=float(coeff))


def norm_check(model):
    P=model['P'];H=model['H'];J=model['J'];up=model['upper']
    polys={'Hpsi':H(P),'Jpsi':J(P),'Gamma_psi':poly.mul(model['Gamma'],P)}
    rows={}
    for name,p in polys.items():
        actual=math.sqrt(poly.norm2(p,model['variance'])/model['norm'])
        bound=up(p);assert actual<=bound
        rows[name]=dict(exact_Gaussian_moment_norm_floating=actual,monomial_triangle_upper=bound)
    return rows


def threshold_integral(n):
    x,w=np.polynomial.legendre.leggauss(n)
    theta=.8;D=1.9;sigmaQ=.5;eta=.3;norm=1.1152
    f=24*theta*(x+1);weights=24*theta*w
    density=f/theta**2*np.exp(-f/theta)
    G=np.array([.5*(1+math.erf((z-2*theta)/(sigmaQ*math.sqrt(2)))) for z in f])
    dG=np.exp(-(f-2*theta)**2/(2*sigmaQ*sigmaQ))/(sigmaQ*math.sqrt(2*math.pi))
    probability=float(weights@(density*(1+eta*eta*(f-2*theta)**2)*G))/norm
    flux=2*eta*D/norm*float(weights@(density*f*dG))
    # Conditional full-H density derivative: 2 eta D (F/theta-2)/norm.
    schrodinger=2*eta*D/norm*float(weights@(density*(f/theta-2)*G))
    return dict(initial_threshold_probability=probability,
        slope_from_probability_current=flux,slope_from_full_H_density_derivative=schrodinger)


def bounded_record_check(certificate):
    row=threshold_integral(240);other=threshold_integral(320)
    residual=max(abs(row[k]-other[k]) for k in row)
    assert residual<1e-11
    assert abs(row['slope_from_probability_current']-row['slope_from_full_H_density_derivative'])<1e-11
    assert 0<row['initial_threshold_probability']<1
    assert row['slope_from_probability_current']>certificate['instantaneous_threshold_slope_positive_bound']
    row.update(independent_quadrature_difference=residual,
        analytic_positive_lower_bound=certificate['instantaneous_threshold_slope_positive_bound'],
        integration_cutoff_F=48*.8,
        quadrature_is_diagnostic_positive_lower_bound_is_analytic=True)
    return row


def unitary(H,F,p,tau):
    E,V=np.linalg.eigh(tau*H+p*F)
    return (V*np.exp(-1j*E))@V.T


def compressed_bit(tau,nq=256):
    H,F,psi,c,info=preceding.compression(8)
    hp=H@psi;zeta=.6;y0=1.6
    x,w=np.polynomial.hermite.hermgauss(nq);w/=math.sqrt(math.pi)
    characteristic=0j;derivative=0j;norm_error=0.
    for p,weight in zip(zeta/2+math.sqrt(2)*x,w):
        U=unitary(H,F,p,tau);V=unitary(H,F,p-zeta,tau)
        a=U@psi;b=V@psi;ah=U@hp;bh=V@hp
        characteristic+=weight*np.vdot(a,b)
        derivative+=weight*1j*(np.vdot(ah,b)-np.vdot(a,bh))
        U0=unitary(H,F,p,0.)
        norm_error=max(norm_error,float(np.linalg.norm((U-U0)@psi)))
    factor=np.exp(-zeta*zeta/8-1j*zeta*y0)
    prob=(1+(factor*characteristic).imag)/2;slope=(factor*derivative).imag/2
    assert 0<prob<1 and slope>0
    assert norm_error<=tau*np.linalg.norm(H,2)+1e-12
    return dict(duration=tau,binary_effect='(1+sin(0.6*(Q_compensated-1.6)))/2',
        probability=float(prob),initial_time_derivative=float(slope),
        maximum_sampled_impulse_state_difference=norm_error,
        finite_matrix_Duhamel_operator_bound=float(tau*np.linalg.norm(H,2)))


def run():
    certificate=source.run()
    assert certificate==json.loads(source.TARGET.read_text('utf8'))
    model=source.build();differential=differential_check(model);norms=norm_check(model)
    record=bounded_record_check(certificate)
    tau=certificate['sufficient_chosen_duration'];c=model['c'];M=certificate['sufficient_mass']
    E=certificate['source_E1'];mu1=math.sqrt(2/math.pi);mu3=2*mu1
    source_change=3*E*mu1/tau+2*(E+c['A'])/tau**2+c['B_squared']*mu3/tau**3
    assert abs(1-M*c['B1']/(2*tau*tau)-.5)<1e-14
    budget=dict(duration=tau,mass=M,initial_source_energy=E,
        initial_pointer_energy=certificate['initial_pointer_energy'],
        sufficient_source_energy_change_upper=source_change,
        total_H_lower_bound=-M*c['B0']/(2*tau*tau),remaining_W_coercivity=.5,
        switch_on_mean_work=0.,switch_off_mean_work_equals_source_energy_change=True,
        unmodeled_controller_and_final_reader_resources_remain=True)
    rows=[compressed_bit(t) for t in (tau,.001,.01)]
    changes=[];convergence=[]
    for row in rows:
        estimates={n:compressed_bit(row['duration'],n) for n in (48,64,96,128,192,320)}
        other=estimates[320]
        changes.append(max(abs(row[k]-other[k]) for k in ('probability','initial_time_derivative')))
        estimates[256]=row
        convergence.append(dict(duration=row['duration'],orders={str(n):{
            k:estimates[n][k] for k in ('probability','initial_time_derivative')}
            for n in sorted(estimates)}))
    assert max(changes)<1e-11
    # Independent limiting calculation of this bounded compressed readout.
    H,F,psi,_,_=preceding.compression(8);E,V=np.linalg.eigh(F)
    C=(V*np.exp(1j*.6*E))@V.T
    slope0=.5*(np.exp(-.6**2/8-1j*.6*1.6)*1j*np.vdot(psi,(H@C-C@H)@psi)).imag
    assert abs(rows[0]['initial_time_derivative']-slope0)<1e-5
    checks=['full_source_polynomial_domain_and_inherited_nonstationary_witness',
        'full_differential_H_and_terminating_gauge_phase_conjugation',
        'Gaussian_source_norms_and_conservative_monomial_upper_bounds',
        'bounded_threshold_probability_current_and_independent_Schrodinger_derivative',
        'finite_duration_positive_response_and_complete_declared_energy_budget',
        'actual_compressed_joint_binary_record_and_independent_pointer_quadrature']
    deps=('joint_probe_time_resolution.py','joint_probe_time_resolution_results.json',
        'joint_finite_time_gauge_probe.py','joint_matter_energy_moment_control.py',
        'joint_gauge_link_reference.py','joint_singlet_common_mass_rg_results.json',
        'round561_drafts/short_time_source_probe.py','round561_drafts/short_time_source_probe_results.json')
    return dict(round=561,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        full_source_certificate=certificate,differential_conjugation=differential,
        Gaussian_norm_cross_checks=norms,actual_bounded_threshold=record,energy_budget=budget,
        compressed_binary_record_diagnostics=rows,
        pointer_quadrature_convergence=convergence,
        initial_48_vs_64_quadrature_failed_tolerance=1e-11,
        compressed_impulse_slope=float(slope0),pointer_quadrature_difference=max(changes),
        scope=dict(full_H_finite_duration_positive_slope_proved_analytically=True,
            high_source_norms_and_final_compensated_quadrature_required=True,
            mass_family_and_high_pointer_energy_are_explicit_inputs=True,
            no_optimal_cost_or_autonomous_complete_apparatus_claim=True,
            secondary_Galerkin_sine_bit_is_different_from_full_source_threshold=True,
            local_time_response_not_four_coordinates_or_spacetime=True,
            floating_evaluations_not_machine_interval_certificates=True,
            no_spacetime_gravity_or_unified_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=561,tests=result['tests_run'],threshold=result['actual_bounded_threshold'],
        budget=result['energy_budget'],compressed=result['compressed_binary_record_diagnostics']),ensure_ascii=False))
