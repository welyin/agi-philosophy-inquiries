"""Round 534: coupled constant backgrounds and a full spectral scale check."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
TARGET=HERE/'spectral_background_feedback_results.json'


def coefficients(n,a,b,c,d,e,f0,f2,f4,cutoff=1):
    return (2*n*f2*cutoff**2-2*f0*c,4*f0*a,
            4*n*f4*cutoff**4-8*f2*cutoff**2*c+2*f0*d,
            -16*a*f2*cutoff**2+8*f0*e,4*f0*b)


def polynomial_data(cs,s):
    A,B,C,D,E=cs
    k=A-B*s; v=C+D*s+E*s*s
    P0=A*D+2*B*C; P1=2*A*E+B*D
    return k,v,P0,P1


def radial_moments(N,wq,wl,ys,M2):
    yv,ye,yu,yd=ys
    return (8*(N*wq+wl),wl*(yv*yv+ye*ye)+N*wq*(yu*yu+yd*yd),
            wl*(yv**4+ye**4)+N*wq*(yu**4+yd**4),
            wl*M2,wl*M2*M2,wl*M2*yv*yv)


def reduced_sphere_action(cs,R,s):
    k,v,_,_=polynomial_data(cs,s)
    # Vol(S4)=384*pi^2/R^2; K=k/(48*pi^2), V=v/(8*pi^2).
    # Euler term is radius/scalar independent and omitted here only.
    return 48*v/R**2-8*k/R


def exact_gaussian_sphere(r):
    """Per finite identity trace, Lambda=1, Phi=0; tail bounds are analytic."""
    stop=max(20,math.ceil(10*r))
    terms=[]; derivatives=[]
    for m in range(2,stop+1):
        multiplicity=F(4,3)*(m**3-m)
        q=float(multiplicity)*math.exp(-(m/r)**2)
        terms.append(q); derivatives.append(q*2*m*m/r**3)
    z=(stop/r)**2
    # Bounds use decreasing x^3 exp(-x^2/r^2) and x^5 exp(-x^2/r^2).
    assert stop >= math.sqrt(2.5)*r
    log_tail_S=math.log(2*r**4/3)+math.log(z+1)-z
    log_tail_derivative=math.log(4*r**3/3)+math.log(z*z+2*z+2)-z
    return dict(radius_times_cutoff=r,terms_through_m=stop,
                exact_action_partial=math.fsum(terms),
                exact_derivative_partial=math.fsum(derivatives),
                positive_first_mode_derivative=64/r**3*math.exp(-4/r**2),
                log_upper_bound_positive_action_tail=log_tail_S,
                log_upper_bound_positive_derivative_tail=log_tail_derivative,
                truncated_action=(2*r**4-2*r*r)/3+11/90,
                truncated_derivative=(8*r**3-4*r)/3)


def rounded(row):
    return {key:round(value,12) if isinstance(value,float) else value for key,value in row.items()}


def run():
    checks=[]; count=0
    for N in (3,5):
        for M2 in (F(0),F(1,2),F(3,2)):
            values=radial_moments(N,F(1),F(2),tuple(map(F,(1,2,3,4))),M2)
            cs=coefficients(*values,F(1),F(1,2),F(1,2))
            A,B,C,D,E=cs
            for s in (F(0),F(1,100),F(1,10)):
                k,v,P0,P1=polynomial_data(cs,s)
                assert k*(D+2*E*s)+2*B*v == P0+P1*s
                count+=1
            edge=A/B
            k,v,P0,P1=polynomial_data(cs,edge)
            assert k==0 and v>0 and P0+P1*edge==2*B*v
    checks.append('exact_coupled_elimination_and_positive_K_endpoint_identity')

    massless_rows=[]
    for N in (3,5):
        for ys in ((1,1,1,1),(1,2,3,4)):
            n,a,b,c,d,e=radial_moments(N,F(1),F(2),tuple(map(F,ys)),F(0))
            cs=coefficients(n,a,b,c,d,e,F(1),F(1,2),F(1,2))
            k,v,P0,P1=polynomial_data(cs,F(0))
            assert P0==32*n*a*(F(1,2)-F(1,4)) and P0>0
            assert P1==8*(n*b-4*a*a) and P1>=0
            R=12*v/k
            assert R==24 and R/F(1,2)==48
            massless_rows.append({'N':N,'Yukawas':list(ys),'P0':str(P0),'P1':str(P1),
                                  'R_over_cutoff_squared':str(R)})
    checks.append('zero_Majorana_branch_has_no_nonzero_constant_radial_solution')

    # Exact polynomial classification controls, including the critical P0=0 case.
    classification=[]
    for C in (F(-1,4),F(0),F(1,4),F(1)):
        cs=(F(2),F(1),C,F(-1),F(1))
        A,B,C,D,E=cs
        _,v_end,P0,P1=polynomial_data(cs,A/B)
        assert v_end>0
        if P0<0:
            s=-P0/P1
            assert P1>0 and 0<s<A/B
            k,v,_,_=polynomial_data(cs,s)
            assert k>0 and k*(D+2*E*s)+2*B*v==0
        else:
            assert min(P0,P0+P1*A/B)>=0
        classification.append({'P0':str(P0),'P1':str(P1),
                               'interior_stationary_point':P0<0})
    checks.append('linear_stationarity_classification_including_zero_P0_control')

    # The exact flat minimum reported in 533 is not a coupled stationary point.
    old=json.loads((HERE/'shared_spectral_budget_results.json').read_text('utf8'))
    s=F(old['flat_radial_minimum_example']['t_squared'])
    cs=coefficients(F(40),F(85),F(1045),F(0),F(0),F(0),F(1),F(1,2),F(1,2))
    k,v,P0,P1=polynomial_data(cs,s)
    A,B,C,D,E=cs
    assert D+2*E*s==0 and P0+P1*s>0
    R=12*v/k
    scalar_residual=6*(D+2*E*s)+B*R
    assert scalar_residual>0
    flat_control={'s':str(s),'R_from_metric_equation':str(R),
                  'scaled_scalar_equation_residual':str(scalar_residual)}
    checks.append('round533_flat_minimum_fails_same_model_scalar_gravity_equations')

    # Smooth cutoff f=exp(-u^32), all four Yukawas nonzero, fixed Majorana.
    p=16; f2=math.gamma(1/p)/(2*p); f4=math.gamma(2/p)/(2*p)
    moments=radial_moments(3,1.,2.,(.1,2.,3.,4.),f2)
    cs=coefficients(*moments,1.,f2,f4)
    A,B,C,D,E=cs
    _,_,P0,P1=polynomial_data(cs,0.)
    s=-P0/P1; k,v,_,_=polynomial_data(cs,s); R=12*v/k
    assert P0<0<P1 and 0<s<A/B and k>0 and v>0
    assert abs(6*(D+2*E*s)+B*R)<1e-9
    # An independent derivative of the sphere action checks both variations.
    eps=1e-25
    dr=reduced_sphere_action(cs,R+eps*1j,s).imag/eps
    ds=reduced_sphere_action(cs,R,s+eps*1j).imag/eps
    assert abs(dr)<1e-10 and abs(ds)<1e-9
    majorana=rounded(dict(p=p,f2=f2,f4=f4,M_squared=f2,s=s,K=k/(48*math.pi**2),
                          P0=P0,P1=P1,R=R,R_over_moment_scale_squared=R/f2,
                          radial_U_second_derivative_in_s=6*P1/k**3))
    checks.append('Majorana_branch_and_independent_sphere_variation_control')

    # Moment-scale curvature diagnostic, invariant under changing width of f.
    scale_rows=[]
    for p in (1,2,8,16):
        f2=math.gamma(1/p)/(2*p); f4=math.gamma(2/p)/(2*p)
        ratio=24*f4/(f2*f2)
        assert ratio>=24
        for width in (.2,1.,5.):
            f2w=width**2*f2;f4w=width**4*f4
            assert abs(24*f4w/f2w/f2w-ratio)<1e-10
        scale_rows.append({'p':p,'R_over_moment_scale_squared':round(ratio,12)})
    checks.append('cutoff_width_invariant_curvature_diagnostic_not_remainder_theorem')

    exact_rows=[]
    for r in (.5,math.sqrt(.5),1.,2.,5.,10.):
        row=exact_gaussian_sphere(r)
        assert row['exact_action_partial']>0 and row['exact_derivative_partial']>0
        assert row['exact_derivative_partial']>=row['positive_first_mode_derivative']>0
        if r>=5:
            assert abs(row['exact_action_partial']/row['truncated_action']-1)<1e-5
        exact_rows.append(rounded(row))
    checks.append('exact_round_sphere_spectrum_with_explicit_positive_tail_bounds')
    fake=exact_gaussian_sphere(math.sqrt(.5))
    assert abs(fake['truncated_derivative'])<1e-14
    assert fake['exact_derivative_partial']>.06 and fake['truncated_action']<0
    assert fake['exact_action_partial']>.002
    checks.append('truncated_symmetric_vacuum_is_not_full_spectral_stationary_point')

    names=['research_note_533.md','shared_spectral_budget.py','shared_spectral_budget_results.json',
           'unified_physics_condition_ledger_533.md','research_note_303.md','research_note_313.md']
    return dict(round=534,tests_run=len(checks),failures=0,errors=0,checks=checks,
        exact_polynomial_cases=count,zero_Majorana_examples=massless_rows,
        algebraic_classification_controls=classification,old_flat_minimum_control=flat_control,
        smooth_cutoff_Majorana_stationary_example=majorana,
        moment_scale_examples=scale_rows,exact_Gaussian_sphere_examples=exact_rows,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        scope=dict(common_four_dimensional_Euclidean_model_supplied=True,
            stationary_constant_radial_truncation_classified=True,
            positive_K_and_strict_cutoff_moment_gap_required=True,
            full_round_sphere_scale_obstruction_for_strictly_decreasing_cutoff=True,
            fixed_cutoff_weights_and_constant_internal_Phi=True,
            no_gauge_flux_fermion_source_volume_constraint_or_extra_action=True,
            full_gravitational_stability_proved=False,
            Lorentzian_or_nonconstant_backgrounds_excluded=False,
            observed_vacuum_or_GR_or_cognition_ruled_out=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
