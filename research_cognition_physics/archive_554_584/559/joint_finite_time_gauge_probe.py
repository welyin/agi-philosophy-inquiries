"""559: finite-time neutral probe, energy budget and protected-fiber record.

The full unbounded source energy theorem is analytic, not a finite matrix
simulation. Fiber records below retain the source K during the pulse. Switching,
pointer mass, duration and final compensated quadrature remain model inputs.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import joint_matter_energy_moment_control as matter
import joint_gauge_link_reference as link

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_finite_time_gauge_probe_results.json'
spec=importlib.util.spec_from_file_location('probe559',HERE/'round559_drafts/finite_duration_probe.py')
probe=importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)


def constants(kappa=.8,v=1.):
    L,C,u=matter.matter(); ell=float(np.linalg.eigvalsh(L)[0])
    A=2*kappa*v*u[0]; B2=8*kappa*kappa*v/ell
    return L,u,dict(kappa=kappa,v=v,ell=ell,uh=float(u[0]),
        A=float(A),B_squared=B2,B0=float(8*kappa*kappa*v*v*u[0]**2),B1=2*B2)


def potential_checks():
    L,u,c=constants(); rng=np.random.default_rng(559)
    margins=[]; invariance=0.; square_residual=0.
    M,g=1.,1.; tau=10.; coefficient=M*g*g/(2*tau*tau)
    assert coefficient*c['B1']<1
    for amplitude in (.03,.3,1.,3.,30.):
        for _ in range(30):
            x,y=rng.normal(size=(2,4))*amplitude
            s=rng.normal(size=2)*amplitude
            q=rng.normal(size=(3,4)); q/=np.linalg.norm(q,axis=1)[:,None]
            U,G,J=map(link.left,q)
            d1=np.array([x@x,s[0]**2])-u; d2=np.array([y@y,s[1]**2])-u
            W=c['v']/4*(d1@L@d1+d2@L@d2)
            F=c['kappa']*c['v']/2*np.sum((x-U@y)**2)
            Fg=c['kappa']*c['v']/2*np.sum((G@x-G@U@J.T@J@y)**2)
            invariance=max(invariance,abs(F-Fg)/max(1.,F))
            assert F<=c['A']+math.sqrt(c['B_squared']*W)+1e-9
            assert F*F<=c['B0']+c['B1']*W+1e-6
            p=float(rng.normal()*50); z=abs(g*p/tau)
            assert z*F<=c['A']*z+W/2+c['B_squared']*z*z/2+1e-8
            total=W+F+p*p/(2*M)+g*p*F/tau
            complete=W+F-coefficient*F*F+(p+M*g*F/tau)**2/(2*M)
            square_residual=max(square_residual,abs(total-complete)/max(1.,abs(total)))
            lower=(1-coefficient*c['B1'])*W+F-coefficient*c['B0']
            margins.append(float(total-lower))
    assert invariance<1e-12 and square_residual<1e-10 and min(margins)>0
    return dict(constants=c,pointwise_samples=len(margins),relative_invariance_error=invariance,
        completing_square_residual=square_residual,minimum_sampled_lower_bound_margin=min(margins),
        sufficient_duration_threshold=math.sqrt(M*g*g*c['B1']/2),
        note='Samples check coefficients; the all-source form bounds are proved in the report.')


def fiber_dynamics_checks():
    saved=json.loads((HERE/'round559_drafts/finite_duration_probe_results.json').read_text('utf8'))
    current=probe.run(); assert current==saved
    rows=[]; energy_residual=0.
    nodes,w=np.polynomial.hermite.hermgauss(56); w/=math.sqrt(math.pi)
    for b,tau in ((1.,5.),(1.,100.),(10.,5.),(10.,100.)):
        K,F,psi,lam=probe.data(.3,b,28); fbar=float(psi@F@psi)
        avg_cost=0.;off_work=0.
        for p,weight in zip(math.sqrt(2)*nodes,w):
            out,_=probe.propagate(p,tau,K,F,psi)
            dK=float(np.vdot(out,K@out).real-lam)
            dF=float(np.vdot(out,F@out).real-fbar)
            energy_residual=max(energy_residual,abs(dK+p*dF/tau))
            assert abs(dK)<=2*.3*abs(p)/tau+1e-11
            avg_cost+=weight*dK
            off_work-=weight*p*np.vdot(out,F@out).real/tau
        assert abs(avg_cost-off_work)<1e-11
        rows.append(dict(b=b,duration=tau,source_energy_change=float(avg_cost),
            switch_off_mean_work=float(off_work),fixed_fiber_energy_bound=2*.3*math.sqrt(2/math.pi)/tau))
    assert energy_residual<1e-10
    return dict(actual_record_rows=current['rows'],energy_work_rows=rows,
        maximum_conditional_conservation_residual=energy_residual)


def gap_check(b,tau,N=28,nq=80):
    alpha=.3; gamma=3*b/4-alpha; sigma=1.
    K,F,psi,lam=probe.data(alpha,b,N); fbar=float(psi@F@psi)
    nodes,w=np.polynomial.hermite.hermgauss(nq); w/=math.sqrt(math.pi)
    joint2=0.;projection_error=0.;shift_error=0.;window_count=0
    for p,weight in zip(math.sqrt(2)*sigma*nodes,w):
        eps=p/tau; eta=abs(eps)*alpha
        out,_=probe.propagate(p,tau,K,F,psi)
        target=np.exp(-1j*(tau*lam+p*fbar))*psi
        error=float(np.linalg.norm(out-target))
        joint2+=weight*error*error
        if eta<=gamma/4:
            window_count+=1
            E,V=np.linalg.eigh(K+eps*(F-alpha*np.eye(N+1)))
            proj=math.sqrt(max(0.,1-float(psi@V[:,0])**2))
            shift=abs(float(E[0])-lam-eps*(fbar-alpha))
            projection_error=max(projection_error,proj)
            shift_error=max(shift_error,shift)
            assert proj<=4*eta/gamma+1e-7
            assert shift<=eta*eta/(gamma-2*eta)+1e-10
            bound=8*eta/gamma+tau*eta*eta/(gamma-2*eta)
            assert error<=bound+1e-10
    cutoff=tau*gamma/(4*alpha)
    tail=math.erfc(cutoff/(math.sqrt(2)*sigma))
    analytic2=(128*alpha**2*sigma**2+24*alpha**4*sigma**4)/(tau*tau*gamma*gamma)+4*tail
    assert joint2<=analytic2+1e-10
    return dict(b=b,duration=tau,whole_group_gap_lower_bound=gamma,
        numerical_joint_state_norm_error_squared=float(joint2),
        analytic_joint_state_norm_error_squared_bound=analytic2,
        bounded_record_event_error_bound=min(1.,math.sqrt(analytic2)),
        Gaussian_tail_formula_floating=tail,window_quadrature_points=window_count,
        largest_sampled_projection_distance=projection_error,largest_sampled_eigenvalue_remainder=shift_error)


def free_pointer_checks():
    alpha,b,tau,M=.3,10.,100.,1.
    K,F,psi,lam=probe.data(alpha,b,28)
    nodes,w=np.polynomial.hermite.hermgauss(56); w/=math.sqrt(math.pi)
    means=np.zeros(3); seconds=np.full(3,.25); residual=0.
    for p,weight in zip(math.sqrt(2)*nodes,w):
        out,der=probe.propagate(p,tau,K,F,psi)
        phase=np.exp(-1j*tau*p*p/(2*M))
        total=phase*out; dtot=phase*(der-1j*tau*p/M*out)
        corrected=dtot+1j*tau*p/M*total
        residual=max(residual,float(np.linalg.norm(corrected-phase*der)))
        for j,(vec,dvec) in enumerate(((out,der),(total,dtot),(total,corrected))):
            means[j]+=weight*(1j*np.vdot(vec,dvec)).real
            seconds[j]+=weight*np.vdot(dvec,dvec).real
    variances=seconds-means*means
    assert residual<1e-12 and abs(means[0]-means[2])<1e-11
    assert abs(variances[0]-variances[2])<1e-10 and variances[1]>9000
    return dict(mass=M,duration=tau,pointer_initial_and_final_energy=1/(2*M),
        record_means=means.tolist(),record_variances=variances.tolist(),
        order=['free_phase_removed','uncompensated_Q','single_Q_minus_tau_P_over_M'],
        compensation_residual=residual,
        final_single_quadrature_is_input_not_joint_sharp_Q_P_measurement=True)


def source_budget_checks():
    _,_,c=constants(); E1=18.; sigma=1.;g=1.;M=1.
    mu1=math.sqrt(2/math.pi)*sigma;mu2=sigma*sigma;mu3=2*math.sqrt(2/math.pi)*sigma**3
    rows=[]
    for tau in (10.,100.,1000.):
        bound=3*E1*g*mu1/tau+2*(E1+c['A'])*g*g*mu2/tau**2+c['B_squared']*g**3*mu3/tau**3
        assert tau*tau>M*g*g*c['B1']/2
        rows.append(dict(duration=tau,absolute_mean_source_energy_change_bound=bound,
            whole_H_lower_bound=-M*g*g*c['B0']/(2*tau*tau),
            W_coercivity_coefficient=1-M*g*g*c['B1']/(2*tau*tau)))
    assert all(rows[i+1]['absolute_mean_source_energy_change_bound']<rows[i]['absolute_mean_source_energy_change_bound'] for i in (0,1))
    return dict(assumed_full_source_E1_upper_bound=E1,
        pointer_absolute_moments=[mu1,mu2,mu3],rows=rows,
        all_b_positive_at_fixed_E1=True,
        budget_is_conditional_input_not_a_simulated_full_matter_record=True)


def run():
    potential=potential_checks()
    fiber=fiber_dynamics_checks()
    gaps=[gap_check(b,t) for b in (1.,10.) for t in (5.,20.,100.)]
    differences=[]
    for row in gaps:
        other=gap_check(row['b'],row['duration'],N=36,nq=104)
        differences.append(abs(row['numerical_joint_state_norm_error_squared']-other['numerical_joint_state_norm_error_squared']))
    assert max(differences)<1e-10
    pointer=free_pointer_checks();budget=source_budget_checks()
    checks=['unbounded_source_potential_and_Young_bounds',
        'common_total_H_coercivity_and_gauge_invariance',
        'exact_fiber_records_source_energy_and_switching_work',
        'whole_group_gap_bound_joint_state_error_and_Gaussian_tail',
        'finite_pointer_kinetic_and_single_compensated_quadrature',
        'one_full_source_energy_budget_all_b_and_finite_duration']
    deps=('joint_matter_energy_moment_control.py','joint_singlet_common_mass_rg_results.json',
        'joint_gauge_link_reference.py','joint_gauge_link_scale_matching_results.json',
        'round559_drafts/finite_duration_probe.py','round559_drafts/finite_duration_probe_results.json')
    return dict(round=559,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        potential_and_stability=potential,fiber_dynamics=fiber,finite_gap_record_bounds=gaps,
        independent_gap_cutoff_and_quadrature_max_difference=max(differences),
        pointer_free_evolution=pointer,full_source_energy_budget=budget,
        scope=dict(hbar_and_numerical_gain_one=True,full_unbounded_source_energy_theorem_analytic=True,
            full_source_Gauss_preserved_each_pointer_branch=True,
            accurate_record_theorem_only_fixed_matter_fast_ground=True,
            finite_energy_not_sufficient_for_full_dynamic_reference_coordinates=True,
            graph_group_source_preparation_pulse_and_final_quadrature_are_inputs=True,
            controller_energy_bookkept_but_no_autonomous_controller_constructed=True,
            slow_probe_changes_original_instantaneous_instrument=True,
            no_spacetime_gravity_or_unified_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else: assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=559,tests=result['tests_run'],
        duration_threshold=result['potential_and_stability']['sufficient_duration_threshold'],
        long_fiber=result['finite_gap_record_bounds'][-1],
        full_source_budgets=result['full_source_energy_budget']['rows']),ensure_ascii=False))
