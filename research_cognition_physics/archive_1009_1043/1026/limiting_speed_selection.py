"""1026: conditional one-loop limiting-speed flow and finite-domain audit.

RK4 refinement differences are cross-checks, not rigorous ODE error bounds.
The hypothetical higher-order envelope is not an established QFT remainder.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "limiting_speed_selection_results.json"
HISTORICAL = [
    "archive_1009_/research_note_1025.md",
    "archive_1009_/1025/input_dependency_update_v0_14.md",
    "archive_1009_/1025/NEXT.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
    *[f"archive_301_341/research_note_{n}.md" for n in (334,335,339,341)],
    "archive_342_369/research_note_352.md",
    "archive_923_934/research_note_924.md",
    "archive_956_989/research_note_957.md",
]
ELLS = np.array([0.,1.,10.,100.,1000.])


def f(a):
    return 1 + 1/a + 4/(3*(1+a)**2)


def b(a):
    return (2*a**3 + 4*a*a + 10*a + 4)/(a*(1+a)**2)


def reduced_rhs(y):
    a,u,z = y
    return np.array([-u*f(a)*(a-1), -b(a)*u*u,
                     4*u*(a-1)/(3*a*(1+a)**2)])


def original_rhs(y):
    """Minus Eq.(18), with variables (c_f,c_b,g), not reduced RHS."""
    cf,cb,g = y
    beta_g = g**3*(3*cb*cf**2+2*cb**2*cf+cb**3+4*cf**3)/(8*math.pi**2*cb*cf**3*(cf+cb)**2)
    beta_b = g*g*(cb*cb-cf*cf)/(8*math.pi**2*cb*cf**3)
    beta_f = g*g*(cf-cb)/(6*math.pi**2*cb*(cf+cb)**2)
    return -np.array([beta_f,beta_b,beta_g])


def reduced_to_original(y):
    a,u,z = y
    cf = math.exp(z)
    return np.array([cf,a*cf,math.sqrt(8*math.pi**2*u*cf**3)])


def original_to_reduced(y):
    cf,cb,g = y
    return np.array([cb/cf,g*g/(8*math.pi**2*cf**3),math.log(cf)])


def trajectory(rhs,y0,u0,density,ells=ELLS):
    if u0 == 0:
        return np.tile(y0,(len(ells),1))
    targets = np.log1p(u0*ells)
    y = np.array(y0,dtype=float)
    answer = [y.copy()]
    position = 0.
    for target in targets[1:]:
        count = max(1,int(math.ceil((target-position)*density)))
        step = (target-position)/count
        def transformed(s,v):
            return math.exp(s)/u0*rhs(s,v)
        for _ in range(count):
            k1 = transformed(position,y)
            k2 = transformed(position+step/2,y+step*k1/2)
            k3 = transformed(position+step/2,y+step*k2/2)
            k4 = transformed(position+step,y+step*k3)
            y += step*(k1+2*k2+2*k3+k4)/6
            position += step
        position = target
        answer.append(y.copy())
    return np.array(answer)


def bounds(a0,u0,ell,eta_a=0.,eta_u=0.):
    lower,upper = min(a0,1),max(a0,1)
    flo,fhi = f(upper)-eta_a*u0,f(lower)+eta_a*u0
    blo,bhi = b(upper)-eta_u*u0,b(lower)+eta_u*u0
    assert min(flo,blo)>0
    return dict(f_min=flo,f_max=fhi,b_min=blo,b_max=bhi,
                u_lower=u0/(1+bhi*u0*ell),u_upper=u0/(1+blo*u0*ell),
                delta_lower=abs(a0-1)*(1+blo*u0*ell)**(-fhi/blo),
                delta_upper=abs(a0-1)*(1+bhi*u0*ell)**(-flo/bhi))


def gauss_common_speed(a0,a1,cf0,order):
    x,w = np.polynomial.legendre.leggauss(order)
    points = (a1-a0)*x/2 + (a1+a0)/2
    integral = (a1-a0)/2*float(w @ (-4/(3*(points+1)**3+4*points)))
    return cf0*math.exp(integral)


def exact_chain_checks():
    cases = []
    for cf,cb,v in [(F(4,5),F(2,5),F(1,100)),(F(6,5),F(9,5),F(1,80)),
                    (F(1),F(1),F(1,50)),(F(2),F(6),F(1,20))]:
        a,u = cb/cf,v/cf**3
        cfprime = 4*v*(cb-cf)/(3*cb*(cf+cb)**2)
        cbprime = -v*(cb*cb-cf*cf)/(cb*cf**3)
        vp = -2*v*v*(3*cb*cf**2+2*cb**2*cf+cb**3+4*cf**3)/(cb*cf**3*(cf+cb)**2)
        ap = cbprime/cf-cb*cfprime/cf**2
        up = vp/cf**3-3*v*cfprime/cf**4
        residuals = [ap+u*f(a)*(a-1),up+b(a)*u*u,
                     cfprime/cf-4*u*(a-1)/(3*a*(1+a)**2)]
        assert all(r==0 for r in residuals)
        cases.append(dict(cf=str(cf),cb=str(cb),v=str(v),residuals=[str(r) for r in residuals]))
    assert f(F(1))==F(7,3) and b(F(1))==5 and f(F(1))/b(F(1))==F(7,15)
    return cases


def check_case(a0,u0,cf0):
    y0 = np.array([a0,u0,math.log(cf0)])
    red = [trajectory(lambda s,y:reduced_rhs(y),y0,u0,k) for k in (128,256,512)]
    orig = trajectory(lambda s,y:original_rhs(y),reduced_to_original(y0),u0,512)
    orig_reduced = np.array([original_to_reduced(y) for y in orig])
    refinement = [float(np.max(np.abs(red[i+1]-red[i]))) for i in (0,1)]
    cross_error = float(np.max(np.abs(red[-1]-orig_reduced)))
    assert cross_error<2e-8 and refinement[-1]<2e-8
    rows = []
    for ell,y,yo in zip(ELLS,red[-1],orig):
        a,u,z = y
        contract = bounds(a0,u0,float(ell))
        delta = abs(a-1)
        assert contract['u_lower']-1e-12<=u<=contract['u_upper']+1e-12
        assert contract['delta_lower']-1e-10<=delta<=contract['delta_upper']+1e-10
        assert min(a0,1)-1e-12<=a<=max(a0,1)+1e-12
        if a0!=1:
            assert (a-1)*(a0-1)>0
        cf=math.exp(z)
        gauss16=gauss_common_speed(a0,a,cf0,16)
        gauss32=gauss_common_speed(a0,a,cf0,32)
        assert abs(cf-gauss32)<2e-9 and abs(gauss16-gauss32)<2e-12
        rows.append(dict(ell=float(ell),a=float(a),u=float(u),cf=cf,cb=a*cf,
                         g=float(yo[2]),absolute_ratio_difference=delta,bounds=contract,
                         speed_integral_error=abs(cf-gauss32),gauss_refinement_error=abs(gauss16-gauss32)))
    if a0==1:
        assert max(abs(row['u']-u0/(1+5*u0*row['ell'])) for row in rows)<2e-10
    if u0==0:
        assert np.array_equal(red[-1],np.tile(y0,(len(ELLS),1)))
    cfs=np.exp(red[-1][:,2]);cbs=red[-1][:,0]*cfs
    direction=np.sign(a0-1)
    assert np.all(direction*np.diff(cfs)>=-1e-12)
    assert np.all(direction*np.diff(cbs)<=1e-12)
    cstar=gauss_common_speed(a0,1.,cf0,32) if u0>0 else None
    if cstar is not None:
        assert min(cf0,a0*cf0)-1e-12<=cstar<=max(cf0,a0*cf0)+1e-12
    return dict(a0=a0,u0=u0,cf0=cf0,rows=rows,original_reduced_max_error=cross_error,
                refinement_max_differences=refinement,
                speed_monotonicity_checked=True,one_loop_mathematical_limit_speed=cstar,
                refinement_difference_is_rigorous_error_bound=False)


def linearization_checks():
    result=[]
    for delta in (-.001,-.0005,.0005,.001):
        u0=.02
        path=trajectory(lambda s,y:reduced_rhs(y),[1+delta,u0,0.],u0,512)
        predicted=delta*(1+5*u0*ELLS)**(-7/15)
        difference=(path[:,0]-1)-predicted
        result.append(dict(delta0=delta,ell=ELLS.tolist(),exact_one_loop_delta=(path[:,0]-1).tolist(),
                           linearized_delta=predicted.tolist(),max_error=float(np.max(np.abs(difference))),
                           max_error_over_delta_squared=float(np.max(np.abs(difference))/delta**2)))
    for full,half in ((result[0],result[1]),(result[3],result[2])):
        ratio=full['max_error']/half['max_error']
        assert 3.95<ratio<4.05
    return result


def remainder_checks():
    result=[]
    for a0 in (.5,2.):
        u0=.02
        def modified(s,y):
            a,u,z=y
            change=reduced_rhs(y)
            change[0]+=u*u*(a-1)*math.sin(s)
            change[1]+=u**3*math.cos(s)
            return change
        paths=[trajectory(modified,[a0,u0,0.],u0,n) for n in (128,256)]
        rows=[]
        for ell,y in zip(ELLS,paths[-1]):
            contract=bounds(a0,u0,float(ell),1,1)
            delta=abs(y[0]-1)
            assert contract['u_lower']-1e-12<=y[1]<=contract['u_upper']+1e-12
            assert contract['delta_lower']-1e-10<=delta<=contract['delta_upper']+1e-10
            rows.append(dict(ell=float(ell),a=float(y[0]),u=float(y[1]),bounds=contract))
        result.append(dict(a0=a0,u0=u0,eta_a=1,eta_u=1,rows=rows,
                           refinement_max_difference=float(np.max(np.abs(paths[1]-paths[0]))),
                           actual_QFT_remainder_certified=False))
    return result


def threshold_checks():
    cases=[]
    for a0 in (.9,1.1,2.):
        u0=.02;epsilon=.01
        c=bounds(a0,u0,0.)
        ratio=abs(a0-1)/epsilon
        necessary=math.expm1(math.log(ratio)*c['b_min']/c['f_max'])/(c['b_min']*u0)
        sufficient=math.expm1(math.log(ratio)*c['b_max']/c['f_min'])/(c['b_max']*u0)
        assert necessary<=sufficient
        low=bounds(a0,u0,necessary)['delta_lower']
        high=bounds(a0,u0,sufficient)['delta_upper']
        assert math.isclose(low,epsilon,rel_tol=1e-12) and math.isclose(high,epsilon,rel_tol=1e-12)
        cases.append(dict(a0=a0,u0=u0,epsilon=epsilon,necessary_ell=necessary,
                          sufficient_ell=sufficient,lower_at_necessary=low,upper_at_sufficient=high,
                          this_is_not_a_certified_physical_RG_range=True))
    return cases


def run():
    exact=exact_chain_checks()
    cases=[check_case(a,u,.8 if index%2==0 else 1.2)
           for index,(a,u) in enumerate((a,u) for a in (.5,.9,1.,1.1,2.) for u in (.005,.02))]
    cases += [check_case(.5,0.,.8),check_case(2.,0.,1.2)]
    linear=linearization_checks()
    remainder=remainder_checks()
    threshold=threshold_checks()
    return dict(round=1026,new_calibration_groups=1,cumulative_test_groups=3803,new_cognitive_axioms=0,
                primary_source='https://arxiv.org/html/1102.0789',primary_equations=[18,19],
                exact_chain_rule_checks=exact,finite_trajectories=cases,near_line_checks=linear,
                hypothetical_remainder_checks=remainder,threshold_checks=threshold,
                near_line_exponent='7/15',finite_domain_scale_variable='ell = log(mu0/mu), not physical time',
                maximum_original_reduced_error=max(c['original_reduced_max_error'] for c in cases),
                maximum_refinement_difference=max(c['refinement_max_differences'][-1] for c in cases),
                maximum_speed_integral_error=max(r['speed_integral_error'] for c in cases for r in c['rows']),
                scheme='two-species positive-speed massless Yukawa one-loop adopted flow',
                one_loop_flow_is_adopted=True,one_loop_beta_functions_rederived_from_loop_integrals=False,
                mathematical_bounds_proved_in_note_not_by_sampling=True,
                numerical_error_rigorously_certified=False,actual_higher_loop_envelope_certified=False,
                physical_common_metric_or_all_Lorentz_operators_derived=False,
                actual_instruments_or_initial_conditions_generated=False,
                infinite_IR_extrapolation_is_required_for_current_effective_goal=False,
                fixed_line_linearization_is_exact_off_line=False,
                exact_finite_scale_equalization_claimed=False,all_scientific_calibrations_passed=True,
                historical_source_sha256={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in HISTORICAL})


def compare(fresh,saved,path='result'):
    if isinstance(fresh,dict):
        assert fresh.keys()==saved.keys(),path
        for key in fresh:compare(fresh[key],saved[key],path+'.'+key)
    elif isinstance(fresh,list):
        assert len(fresh)==len(saved),path
        for i,(a,c) in enumerate(zip(fresh,saved)):compare(a,c,path+f'[{i}]')
    elif isinstance(fresh,float):
        assert math.isfinite(fresh) and math.isclose(fresh,saved,rel_tol=2e-8,abs_tol=3e-12),(path,fresh,saved)
    else:
        assert fresh==saved,(path,fresh,saved)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run()
    if args.write:
        with OUT.open('x',encoding='utf8') as handle:json.dump(result,handle,ensure_ascii=False,indent=2);handle.write('\n')
    else:compare(result,json.loads(OUT.read_text('utf8')))
    print(json.dumps(dict(round=1026,passed=True,trajectories=len(result['finite_trajectories']),
                          max_cross_error=result['maximum_original_reduced_error']),ensure_ascii=False))
