"""565: fixed-parameter large-electric-energy limit of matter and its radial reader.

The full-model convergence statement is analytic via closed forms.
Finite form compressions here are diagnostics, not an untruncated rate certificate.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_full_matter_large_b_limit_results.json'
spec=importlib.util.spec_from_file_location('probe565',HERE/'round565_drafts/large_b_matter_probe.py')
probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)


def independent_moments():
    n=3;polys=probe.laguerre_polys(n);worst=0.;rows={}
    # Direct radial-coordinate integral avoids the endpoint issue of x^(-1).
    x,w=np.polynomial.legendre.leggauss(196);r=6*(x+1);weight=6*w*2*r**3*np.exp(-r*r)
    values=np.array([np.polynomial.polynomial.polyval(r*r,p) for p in polys])
    for power in (-1.,0.,.5,1.,2.):
        moment=(values*(weight*r**(2*power)))@values.T
        exact=probe.radial_moment(n,power)
        error=float(np.max(abs(moment-exact)));worst=max(worst,error);rows[str(power)]=error
    assert worst<2e-11
    # Independent Hermite-function quadrature versus extended coordinate powers.
    ns=3;x,w=np.polynomial.hermite.hermgauss(12);w/=math.sqrt(math.pi)
    values=np.array([np.polynomial.hermite.hermval(x,[0]*j+[1])/math.sqrt(2**j*math.factorial(j)) for j in range(ns)])
    q=np.diag(np.sqrt(np.arange(1,ns+4)/2),1)+np.diag(np.sqrt(np.arange(1,ns+4)/2),-1)
    for power in (2,4):
        actual=(values*(w*x**power))@values.T
        target=np.linalg.matrix_power(q,power)[:ns,:ns]
        error=float(np.max(abs(actual-target)));worst=max(worst,error)
        assert error<1e-12
    return dict(radial_power_matrix_residuals=rows,combined_maximum_residual=worst,
        exact_moments_not_squared_truncated_coordinates=True)


def projected_potential():
    L,u,_=probe.inherited.constants();k=.8
    x,w=np.polynomial.legendre.leggauss(160);chi=np.pi*(x+1)/2;z=np.cos(chi)
    angular_weights=w*np.sin(chi)**2
    def W(r,s):
        d=np.array([r*r,s*s])-u;return float(d@L@d/4)
    def averaged(r1,r2):
        full=W(r1,.2)+W(r2,-.3)+k*(r1*r1+r2*r2)/2-k*r1*r2*z
        return float(angular_weights@full)
    r1=(.7,1.4);r2=(.8,1.6)
    mixed=averaged(r1[1],r2[1])-averaged(r1[1],r2[0])-averaged(r1[0],r2[1])+averaged(r1[0],r2[0])
    wrong_retained_radial_gradient=-k*(r1[1]-r1[0])*(r2[1]-r2[0])
    assert abs(mixed)<1e-13 and abs(wrong_retained_radial_gradient)>.1
    m=probe.build();H0=m['H'](0.)
    assert np.linalg.eigvalsh(H0)[0]>0
    assert np.linalg.norm(m['R']@m['penalty']-m['penalty']@m['R'])<1e-13
    for b in (1.,10.,100.):
        projected=m['H'](b)[::m['nang'],::m['nang']]
        assert np.max(abs(projected-m['HZ']))<1e-12
        assert np.max(abs(m['H'](b)-b*m['penalty']-H0))<1e-12
    return dict(Haar_averaged_potential_mixed_difference=mixed,
        wrongly_retained_k_R_mixed_difference=wrong_retained_radial_gradient,
        constant_sector_full_form_equals_sum_of_two_cell_forms=True,
        R_instrument_commutes_with_angular_projection=True,
        source_factorization_not_claimed_for_external_RP_instrument=True)


def common_apparatus_budget(certificate):
    L,u,constants=probe.inherited.constants();ell=constants['ell']
    assert abs(ell-np.linalg.eigvalsh(L)[0])<1e-12
    D0=2*u[0]**2;D1=4/ell;E=certificate['initial_source_energy']
    gaussian_W=.5*(L[0,0]*(6-4*u[0]+u[0]**2)
        +2*L[0,1]*(2-u[0])*(.5-u[1])+L[1,1]*(.75-u[1]+u[1]**2))
    assert abs(E-(2.5+1.6+gaussian_W))<1e-12
    tau=.2;g=1.;sigmaQ=.5;mu2=1.;M=tau*tau/(g*g*D1)
    theta=M*g*g*D1/(2*tau*tau);C=M*g*g*D0/(2*tau*tau)
    pointer_E=mu2/(2*M);whole=E+pointer_E
    mu1=math.sqrt(2/math.pi);mu3=2*mu1;k=.8;Br2=2/ell
    source_change=3*g*E*mu1/(k*tau)+2*g*g*(E/k+u[0])*mu2/(k*tau*tau)+Br2*g**3*mu3/(k*tau**3)
    Rmean=2-9*math.pi/16
    impulse_change=.5*Rmean/sigmaQ**2
    assert abs(theta-.5)<1e-14 and Rmean>0
    # Coercive bound H_tot + C >= (1-theta) H_b yields a b-uniform pulse bound.
    leakage_numerator=4*(whole+C)/(3*(1-theta))
    return dict(tau=tau,g=g,sigmaQ=sigmaQ,M=M,relative_form_bound=math.sqrt(theta),
        retained_W_fraction=1-theta,uniform_lower_bound=-C,
        fixed_initial_pointer_energy=pointer_E,whole_initial_energy=whole,
        same_b_independent_source_energy=E,finite_pulse_source_change_sufficient_upper=source_change,
        pulse_angular_leakage_upper=f'min(1, {leakage_numerator}/b)',
        source_R_mean=Rmean,impulse_matter_energy_increase=impulse_change,
        instantaneous_nonselective_electric_increase=0.,
        fixed_tau_M_g_as_b_changes=True,no_time_or_coupling_rescaling_hidden=True)


def run():
    certificate=probe.run();assert certificate==json.loads(probe.TARGET.read_text('utf8'))
    moments=independent_moments();projection=projected_potential();budget=common_apparatus_budget(certificate)
    checks=['correct_full_radial_singlet_form_moments_by_independent_integrals',
        'Haar_projected_full_potential_and_missing_cross_cell_source_interaction',
        'full_matter_compression_dynamics_energy_and_angular_leakage',
        'finite_duration_continuous_pointer_joint_state_comparison',
        'b_independent_whole_apparatus_coercivity_and_source_budget']
    deps=('round565_drafts/large_b_matter_probe.py','round565_drafts/large_b_matter_probe_results.json',
        'joint_finite_time_gauge_probe.py','joint_radial_prediction_closure.py',
        'joint_radial_prediction_closure_results.json','joint_direct_radial_record.py')
    return dict(round=565,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        independent_form_moments=moments,projected_source=projection,
        finite_form_diagnostics=certificate,common_apparatus_budget=budget,
        scope=dict(full_unbounded_convergence_is_analytic_not_from_Galerkin=True,
            generalized_resolvent_limit_lives_only_on_P0=True,
            real_time_strong_limit_for_fixed_P0_input_compact_time_intervals=True,
            no_full_space_unitary_limit_or_uniform_all_source_rate_claim=True,
            fixed_a_k_W_and_fixed_tau_M_g=True,
            unmeasured_limit_matter_factorizes_but_RP_reader_can_connect_cells=True,
            no_spacetime_dimension_gravity_or_unified_completion=True,
            floating_diagnostics_not_machine_interval_certificates=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=565,tests=result['tests_run'],moments=result['independent_form_moments'],
        projection=result['projected_source'],budget=result['common_apparatus_budget']),ensure_ascii=False))
