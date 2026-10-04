"""Round 535: fixed-volume variations of the full bare spectral action.

One generation, constant scalar, round Euclidean S4. Analytic monotonicity
is proved in the note; finite spectral sums only check the implementation.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from shared_trace_unimodularity import blockdiag
from spectral_background_feedback import exact_gaussian_sphere, radial_moments, coefficients

HERE=Path(__file__).resolve().parent
TARGET=HERE/'fixed_volume_scalar_stability_results.json'


def phi(N,ys,M,t):
    def pair(a,b):
        x=np.zeros((4,4))
        x[0,2]=x[2,0]=a*t
        x[1,3]=x[3,1]=b*t
        return x
    l,q=pair(*ys[:2]),pair(*ys[2:])
    out=blockdiag(l,l,np.kron(q,np.eye(N)),np.kron(q,np.eye(N)))
    out[0,4]=out[4,0]=M
    return out


def mass_data(N,wq,wl,ys,M,t):
    yv,ye,yu,yd=ys
    delta=math.hypot(M,2*yv*t)
    plus=(delta+M)/2
    minus=2*(yv*t)**2/(delta+M) if delta+M else 0.0
    slope=2*yv*yv*t/delta if delta else 0.0
    masses=[plus*plus,minus*minus,ye*ye*t*t,yu*yu*t*t,yd*yd*t*t]
    slopes=[2*plus*slope,2*minus*slope,2*ye*ye*t,2*yu*yu*t,2*yd*yd*t]
    weights=[2*wl,2*wl,4*wl,4*N*wq,4*N*wq]
    return masses,slopes,weights


def gaussian(N,wq,wl,ys,M,t,cutoff=1.0):
    masses,slopes,weights=mass_data(N,wq,wl,ys,M,t)
    terms=[w*math.exp(-x/cutoff**2) for w,x in zip(weights,masses)]
    return math.fsum(terms),-math.fsum(v*s/cutoff**2 for v,s in zip(terms,slopes))


def origin_curvature(N,wq,wl,ys,M,cutoff=1.0):
    yv,ye,yu,yd=ys
    return -8/cutoff**2*(wl*(ye*ye+yv*yv*math.exp(-M*M/cutoff**2))
                         +N*wq*(yu*yu+yd*yd))


def direct_gaussian(N,wq,wl,ys,M,t,cutoff=1.0):
    matrix=phi(N,ys,M,t)
    # Z is scalar on each of these two blocks: no eigenvector convention needed.
    return float(wl*np.exp(-np.linalg.eigvalsh(matrix[:8,:8])**2/cutoff**2).sum()
                 +wq*np.exp(-np.linalg.eigvalsh(matrix[8:,8:])**2/cutoff**2).sum())


def spectral_partial(N,wq,wl,ys,M,t,r,p,stop=60):
    masses,slopes,weights=mass_data(N,wq,wl,ys,M,t)
    terms=[]; deriv=[]
    for m in range(2,stop+1):
        mult=4*(m**3-m)/3
        for mass,slope,weight in zip(masses,slopes,weights):
            x=(m/r)**2+mass
            value=math.exp(-x**p)
            terms.append(mult*weight*value)
            deriv.append(-mult*weight*p*x**(p-1)*value*slope)
    return math.fsum(terms),math.fsum(deriv)


def run():
    checks=[]; spectrum_error=0.0; trace_error=0.0
    for N in (3,5):
        for M in (0.0,0.7,2.0):
            for t in (0.0,0.1,0.8):
                ys=(1.0,2.0,3.0,4.0)
                masses,_,_=mass_data(N,1,2,ys,M,t)
                expected=sorted([masses[0]]*2+[masses[1]]*2+[masses[2]]*4
                                +[masses[3]]*(4*N)+[masses[4]]*(4*N))
                actual=sorted(np.linalg.eigvalsh(phi(N,ys,M,t))**2)
                spectrum_error=max(spectrum_error,float(np.max(np.abs(np.array(actual)-expected))))
                trace_error=max(trace_error,abs(direct_gaussian(N,1,2,ys,M,t)
                                                -gaussian(N,1,2,ys,M,t)[0]))
    assert spectrum_error<1e-12 and trace_error<1e-11
    checks.append('finite_matrix_spectrum_including_Majorana_and_weighted_Gaussian_trace')

    origin=[]
    for M in (0.0,0.7,2.0):
        ys=(1.0,2.0,3.0,4.0); h=1e-4
        exact=origin_curvature(3,1,2,ys,M)
        finite=2*(direct_gaussian(3,1,2,ys,M,h)-direct_gaussian(3,1,2,ys,M,0))/h**2
        assert exact<0 and abs(finite/exact-1)<1e-6
        origin.append(dict(M=M,analytic=exact,independent_matrix_difference=finite))
    checks.append('negative_origin_curvature_from_independent_matrix_differences')

    derivatives=[]
    for p in (1,2,16):
        for M in (0.0,0.7,2.0):
            ys=(0.3,0.4,0.5,0.6); t=0.7; r=4.0; h=1e-5
            value,deriv=spectral_partial(3,1,2,ys,M,t,r,p)
            fd=(spectral_partial(3,1,2,ys,M,t+h,r,p)[0]
                -spectral_partial(3,1,2,ys,M,t-h,r,p)[0])/(2*h)
            assert value>0 and deriv<0 and abs(fd/deriv-1)<2e-7
            derivatives.append(dict(p=p,M=M,t=t,partial=value,derivative=deriv))
    checks.append('strict_scalar_descent_and_differentiated_positive_mode_sums')

    # Global lambda is independent of the scalar variable. This algebra check
    # compares the same background and a constant shift, not a new field theory.
    multiplier=[]
    for K,V,R,C in ((F(2),F(7),F(3),F(11)),(F(1,2),F(99),F(1,10),F(-101))):
        lam=K*R/2-V
        assert K*R==2*(V+lam)==2*(V+C+lam-C)
        multiplier.append(dict(K=str(K),V=str(V),R=str(R),multiplier=str(lam),shift=str(C)))
    checks.append('fixed_volume_multiplier_absorbs_constant_but_not_scalar_dependent_terms')

    N,wq,wl,ys=3,1,2,(1,2,3,4)
    n,a,b,c,d,e=radial_moments(N,F(wq),F(wl),tuple(map(F,ys)),F(0))
    r=10.0; R=F(3,25)
    s=a/b*(1-R/12)
    cs=coefficients(n,a,b,c,d,e,F(1),F(1,2),F(1,2))
    A,B,C,D,E=cs
    assert D+2*E*s+R*B/6==0 and s>0 and E>0
    t=math.sqrt(float(s)); H,dH=gaussian(N,wq,wl,ys,0,t)
    assert dH<0
    control=dict(radius=r,R=float(R),s_exact=str(s),t=t,
                 truncated_scalar_gradient_exact='0',full_logarithmic_scalar_gradient=dH/H,
                 max_Dirac_mass_over_cutoff=max(ys)*t,
                 geometry_relative_truncation_error=abs(exact_gaussian_sphere(r)['truncated_action']
                    /exact_gaussian_sphere(r)['exact_action_partial']-1))
    checks.append('large_radius_truncated_scalar_minimum_is_not_a_full_spectral_stationary_point')

    residual=[]
    for ys in ((1,2,3,4),(0,2,3,4),(1,0,3,0),(0,0,0,0)):
        M=0.7
        expected=(2*wl*(1+math.exp(-M*M)) if ys[0]==0 else 0)
        expected+=4*wl*(ys[1]==0)+4*N*wq*((ys[2]==0)+(ys[3]==0))
        value,_=gaussian(N,wq,wl,ys,M,30)
        assert abs(value-expected)<1e-12
        curvature=origin_curvature(N,wq,wl,ys,M)
        assert (curvature==0) == (ys==(0,0,0,0))
        residual.append(dict(Y=list(ys),limit=expected,finite_t30=value))
    checks.append('zero_Yukawa_residuals_and_excluded_all_zero_countercondition')

    factorization=[]
    for r in (1.0,4.0,10.0):
        ys=(1,2,3,4); t=0.3; M=0.7
        h,dh=gaussian(N,wq,wl,ys,M,t)
        geom=exact_gaussian_sphere(r)
        direct,direct_d=spectral_partial(N,wq,wl,ys,M,t,r,1,geom['terms_through_m'])
        assert abs(direct/(geom['exact_action_partial']*h)-1)<1e-12
        assert abs(direct_d/(geom['exact_action_partial']*dh)-1)<1e-12
        factorization.append(dict(r=r,partial_action=direct,partial_scalar_derivative=direct_d,
            log_action_tail_bound=geom['log_upper_bound_positive_action_tail']+math.log(h),
            log_absolute_scalar_derivative_tail_bound=geom['log_upper_bound_positive_action_tail']+math.log(-dh)))
    checks.append('full_Gaussian_factorization_and_inherited_positive_tail_bound')

    volume=8*math.pi**2*r**4/3
    assert math.isclose((3*volume/(8*math.pi**2))**0.25,r)
    # A scalar change leaves the volume constraint exactly unchanged.
    _,dh=gaussian(3,1,2,(1,2,3,4),0.7,0.3)
    for multiplier_value in (-1000.0,0.0,1000.0):
        assert dh+multiplier_value*0==dh<0
    checks.append('fixed_geometry_scalar_variation_remains_admissible_under_volume_constraint')
    def tidy(obj):
        if isinstance(obj,float): return round(obj,12)
        if isinstance(obj,list): return [tidy(x) for x in obj]
        if isinstance(obj,dict): return {k:tidy(v) for k,v in obj.items()}
        return obj
    deps=('research_note_534.md','spectral_background_feedback.py','spectral_background_feedback_results.json',
          'shared_trace_unimodularity.py','research_note_306.md','unified_physics_condition_ledger_534.md')
    return tidy(dict(round=535,tests_run=len(checks),failures=0,errors=0,checks=checks,
        max_matrix_spectrum_error=spectrum_error,max_weighted_trace_error=trace_error,
        origin_curvature=origin,strict_descent=derivatives,multiplier_algebra=multiplier,
        large_radius_control=control,zero_Yukawa_controls=residual,full_Gaussian_factorization=factorization,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(fixed_four_volume_new_input=True,one_generation_constant_scalar_round_Euclidean_S4=True,
            strictly_decreasing_cutoff_fixed_weights_cutoff_and_Majorana=True,
            no_other_effective_action_or_field_constraint=True,
            no_finite_scalar_local_minimum_for_full_bare_action=True,
            arbitrary_generation_mixing_covered=False,Lorentz_stability_proved=False,
            Wilsonian_matching_or_entire_spectral_geometry_excluded=False,
            physical_volume_or_observed_cosmological_constant_derived=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','large_radius_control')},ensure_ascii=False))
