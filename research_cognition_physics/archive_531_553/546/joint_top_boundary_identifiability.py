"""546: strict top-boundary identifiability and simultaneous spectral diagnostics.

One-loop fixed gauge trajectory. Synthetic MS targets, not measured pole inputs.
The inverse theorem is analytic; shooting brackets are floating point diagnostics.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_singlet_common_mass_rg as rg
import joint_singlet_spectrum_and_weights as spectral

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_top_boundary_identifiability_results.json'
END=math.log(rg.MU/173.34)


def gauge_combinations(u):
    gy,gw,gc=rg.gauge_squared(-u,rg.XSTAR)
    return 17*gy/12+9*gw/4+8*gc,3*gy+9*gw


def jacobian(u,values):
    q,S,z=values;A,B=gauge_combinations(u)
    return np.array([[2*A-18*q-2*S,-2*q,0.],[-6*S,B/2-6*q-10*S-z,-S],
                     [0.,-2*z,-20*z-2*S]])/rg.LOOP


def rhs(u,state):
    q,S,z=state[:3];A,B=gauge_combinations(u)
    base=[2*q*(A-4.5*q-S),S*(B/2-6*q-5*S-z),-2*z*(5*z+S)]
    if len(state)==6:
        dq,dS,dz=state[3:]
        base += [(2*A-18*q-2*S)*dq-2*q*dS,
                 -6*S*dq+(B/2-6*q-10*S-z)*dS-S*dz,
                 -2*z*dS-(20*z+2*S)*dz]
    return np.array(base)/rg.LOOP


def flow(q0,u,steps=600,variational=False):
    q0=np.asarray(q0,float);S0=(rg.T-3*q0)/rg.R
    state=[q0,S0,np.full_like(q0,rg.T/(2*rg.R))]
    if variational:state += [np.ones_like(q0),np.full_like(q0,-3/rg.R),np.zeros_like(q0)]
    state=np.array(state);dt=u/steps
    for j in range(steps):
        t=j*dt;a=rhs(t,state);b=rhs(t+dt/2,state+dt*a/2)
        c=rhs(t+dt/2,state+dt*b/2);d=rhs(t+dt,state+dt*c)
        state+=dt*(a+2*b+2*c+d)/6
    assert np.all(np.isfinite(state))
    return state


def upper_endpoint(u,steps=2048):
    """Independent Bernoulli quadrature for S=0. Composite Simpson rule."""
    assert steps%2==0
    grid=np.linspace(0,u,steps+1)
    b=np.array([-41/6,19/6,7]);c=np.array([17/12,9/4,8])
    factors=1-b[:,None]*grid[None,:]/(8*math.pi**2*np.array(rg.XSTAR)[:,None])
    assert np.all(factors>0)
    E=np.prod(factors**(-c[:,None]/b[:,None]),axis=0)
    integral=(u/steps)/3*(E[0]+E[-1]+4*np.sum(E[1:-1:2])+2*np.sum(E[2:-1:2]))
    return float(E[-1]/(3/rg.T+9*integral/rg.LOOP))


def shoot(targets,u=END,steps=600,iterations=32):
    targets=np.array(targets,float)
    assert min(targets)>0 and max(targets)<upper_endpoint(u)
    low=np.zeros_like(targets);high=np.full_like(targets,rg.T/3)
    for _ in range(iterations):
        mid=(low+high)/2;got=flow(mid,u,steps)[0]
        lower=got<targets;low=np.where(lower,mid,low);high=np.where(lower,high,mid)
    return low,high,(low+high)/2


def run():
    checks=[];finite_difference_error=0.;equation_error=0.
    transform=np.diag([1,-1,1])
    for u,values in ((0.,[.2,.3,.4]),(END/2,[.4,.2,.3]),(END,[.8,0.,.2]),(0.,[0.,.7,.3])):
        J=jacobian(u,values);metzler=transform@J@transform
        assert min(metzler[i,j] for i in range(3) for j in range(3) if i!=j)>=0
        a=np.array(values);eps=1e-6
        numerical=np.column_stack([(rhs(u,a+eps*np.eye(3)[j])-rhs(u,a-eps*np.eye(3)[j]))/(2*eps) for j in range(3)])
        finite_difference_error=max(finite_difference_error,float(np.max(abs(J-numerical))))
        padded=np.r_[a,.2,.1,.4,.8,.7]
        equation_error=max(equation_error,float(np.max(abs(rhs(u,a)-rg.rhs(u,padded)[:3]))))
    assert finite_difference_error<1e-9 and equation_error<1e-15
    checks.append('independent_variational_Jacobian_cooperative_signs_and_original_equation_agreement')

    sensitivity_rows=[];variation_error=0.;flow_error=0.
    for q0 in (.15,.3,.5):
        value=flow(q0,END,800,True);fine=flow(q0,END,1600,True)
        flow_error=max(flow_error,float(np.max(abs(value-fine))))
        eps=1e-5
        diff=(flow(q0+eps,END,1600)-flow(q0-eps,END,1600))/(2*eps)
        variation_error=max(variation_error,float(np.max(abs(diff-fine[3:]))))
        assert fine[3]>0 and fine[4]<0 and fine[5]>0
        sensitivity_rows.append(dict(q0=q0,endpoint=fine[:3].tolist(),derivative=fine[3:].tolist()))
    assert variation_error<1e-8 and flow_error<1e-10
    checks.append('coupled_sensitivity_flow_step_doubling_and_independent_parameter_differences')

    zero=flow(0.,END,1600);upper=flow(rg.T/3,END,1600)
    closed=upper_endpoint(END,4096);coarse=upper_endpoint(END,2048)
    assert zero[0]==0 and abs(upper[1])<1e-15
    assert abs(upper[0]-closed)<1e-10 and abs(closed-coarse)<1e-10
    # Logistic upper bounds on the chosen finite gauge interval, not bounds on all fields.
    grid=np.linspace(0,END,1001);AB=np.array([gauge_combinations(u) for u in grid])
    qbound=max(rg.T/3,2*max(AB[:,0])/9);sbound=max(rg.T/rg.R,max(AB[:,1])/10)
    sampled=flow(np.linspace(0,rg.T/3,13),END,800)
    assert np.all(sampled[0]>=0) and np.max(sampled[0])<qbound
    assert np.min(sampled[1])>-1e-15 and np.max(sampled[1])<sbound
    assert np.min(sampled[2])>0 and np.max(sampled[2])<rg.T/(2*rg.R)
    checks.append('zero_endpoint_and_independent_Bernoulli_reachable_upper_endpoint')

    targets=np.array([.4,.8,.9,1.]);lo,hi,roots=shoot(targets)
    lo_values=flow(lo,END,1600)[0];hi_values=flow(hi,END,1600)[0]
    refined=flow(roots,END,1600,True)
    assert np.all(lo_values<=targets+1e-11) and np.all(hi_values>=targets-1e-11)
    assert np.max(hi-lo)<2e-10 and np.max(abs(refined[0]-targets))<5e-10
    checks.append('unique_boundary_inverse_numerical_brackets_and_refined_target_residuals')

    rows=[]
    for i,q0 in enumerate(roots):
        coarse=rg.flow(q0,END,800);full=rg.flow(q0,END,1600)
        err=float(np.max(abs(coarse-full)))
        assert err<2e-10 and np.max(abs(full[:3]-refined[:3,i]))<2e-11
        diag=rg.diagnostics(full)
        row=dict(synthetic_target_q=float(targets[i]),q0=float(q0),q0_bracket=[float(lo[i]),float(hi[i])],
            endpoint_q=float(full[0]),target_residual=float(full[0]-targets[i]),
            derivative_of_endpoint_top=float(refined[3,i]),step_doubling_error=err,
            vacuum=diag)
        if diag['strict_double_vev']:row['joint_spectrum']=spectral.row_from_state(full)
        rows.append(row)
    assert not rows[0]['vacuum']['strict_double_vev']
    assert all(r['vacuum']['strict_double_vev'] for r in rows[1:])
    checks.append('same_unique_boundary_joint_spectra_and_noncoexistence_counterexample')

    deps=('research_note_544.md','research_note_545.md','joint_singlet_common_mass_rg.py',
          'joint_singlet_spectrum_and_weights.py','research_round_545_checks.json')
    return dict(round=546,tests_run=len(checks),failures=0,errors=0,checks=checks,
        gauge_path_input_fixed=True,endpoint_log_scale=END,endpoint_scale_GeV=173.34,
        q0_interval=[0,rg.T/3],reachable_upper_top_squared=closed,
        Jacobian_difference_max=finite_difference_error,original_rhs_difference_max=equation_error,
        parameter_difference_max=variation_error,sensitivity_step_doubling_max=flow_error,
        sensitivity_examples=sensitivity_rows,inverse_examples=rows,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(strict_top_monotonicity_proved_for_the_declared_one_loop_family=True,
            other_spectral_quantities_monotonicity_proved=False,
            all_scalar_trajectories_exist_or_stable_proved=False,
            targets_are_synthetic_MS_parameters_not_experimental_fits=True,
            shooting_is_not_validated_interval_arithmetic=True,
            thresholds_pole_extraction_and_unification_completed=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','reachable_upper_top_squared',
        'parameter_difference_max','sensitivity_step_doubling_max')},ensure_ascii=False))
