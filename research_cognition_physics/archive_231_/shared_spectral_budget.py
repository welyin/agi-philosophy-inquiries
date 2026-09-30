"""Round 533: joint bare coefficients of one weighted a0+a2+a4 action.

This tests an explicitly supplied Euclidean truncation, not observed vacuum energy.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from shared_trace_unimodularity import structures, rep

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'shared_spectral_budget_results.json'


def reduced_coefficients(weights, squared_eigenvalues, f0, f2, f4, cutoff=F(1)):
    """Return 8*pi^2*V and 48*pi^2*K, using exact rational arithmetic."""
    n = sum(weights)
    t2 = sum(w*x for w, x in zip(weights, squared_eigenvalues))
    t4 = sum(w*x*x for w, x in zip(weights, squared_eigenvalues))
    v = 4*n*f4*cutoff**4 - 4*f2*cutoff**2*t2 + f0*t4
    k = 2*n*f2*cutoff**2 - f0*t2
    variance = t4 - t2*t2/n
    moment_gap = f4 - f2*f2/f0
    # Exact identity; its two final terms are nonnegative under the stated inputs.
    assert v == k*k/(n*f0) + f0*variance + 4*n*cutoff**4*moment_gap
    return n, t2, t4, v, k, variance, moment_gap


def close(a, b):
    assert np.allclose(a, b, rtol=2e-12, atol=1e-10), (a, b)


def run():
    checks = []
    rational_cases = 0
    for weights in ([F(1), F(2)], [F(2, 3), F(7, 4), F(5)]):
        for j in range(8):
            xs = [F((i+j)**2, 7) for i in range(len(weights))]
            for cutoff in (F(1), F(3, 2)):
                n,t2,t4,v,k,var,gap = reduced_coefficients(weights,xs,F(1),F(1,2),F(1,2),cutoff)
                assert var >= 0 and gap > 0 and v >= k*k/n
                rational_cases += 1
    checks.append('exact_square_and_weighted_variance_identity')

    # A decreasing function represented as a positive sum of interval indicators.
    # These nonsmooth examples test moments only, not the smooth heat expansion.
    moment_rows = []
    for atoms in ([(F(1),F(1))], [(F(2),F(1)),(F(3),F(2))],
                  [(F(1,3),F(2)),(F(5,7),F(4)),(F(2),F(5))]):
        f0=sum(c for c,r in atoms)
        f2=sum(c*r*r/2 for c,r in atoms)
        f4=sum(c*r**4/4 for c,r in atoms)
        assert f0*f4 >= f2*f2
        moment_rows.append([str(x) for x in (f0,f2,f4,f4-f2*f2/f0)])
    checks.append('positive_layer_cake_moments_and_single_step_equality')

    # Smooth, even, rapidly decreasing f_p(u)=exp(-u^(2p)).
    smooth_rows=[]
    for p in (1,2,4,8,16,32):
        f2=math.gamma(1/p)/(2*p)
        f4=math.gamma(2/p)/(2*p)
        ratio=f2*f2/f4
        assert 0 < ratio < 1
        smooth_rows.append({'p':p,'f2':round(f2,12),'f4':round(f4,12),
                            'f2_squared_over_f0_f4':round(ratio,12)})
    assert smooth_rows[-1]['f2_squared_over_f0_f4'] > .998
    checks.append('smooth_monotone_cutoffs_approach_moment_equality')

    gauge_rows=[]
    for N in (3,5):
        for wq,wl in ((1.,1.),(1.,2.),(.5,3.)):
            Z,J,_,_=structures(N,wq,wl)
            A=rep(N,0,np.diag([.5,-.5]),np.zeros((N,N)))
            B=A-J@A.conj()@J
            index=float(np.trace(Z@B@B).real)
            n=float(np.trace(Z).real)
            close(index,N*wq+wl)
            close(n,8*index)
            for ng in (1,3):
                f0=2.7
                inverse_gw2=f0*ng*index/(6*math.pi**2)
                close(f0*ng*n,48*math.pi**2*inverse_gw2)
            gauge_rows.append({'N':N,'wq':wq,'wl':wl,'TrZ':n,'full_weak_index':index})
    checks.append('actual_weak_generator_and_shared_gauge_gravity_normalization')

    matrix_rows=[]
    for N in (3,5):
        for wq,wl in ((1.,1.),(1.,2.)):
            for t in (0.,.1,.3):
                Z,_,_,D=structures(N,wq,wl,t)
                n=float(np.trace(Z).real)
                t2=float(np.trace(Z@D@D).real)
                t4=float(np.trace(Z@np.linalg.matrix_power(D,4)).real)
                f0,f2,f4,cutoff=1.,.5,.5,20.
                V=(4*n*f4*cutoff**4-4*f2*cutoff**2*t2+f0*t4)/(8*math.pi**2)
                K=(2*n*f2*cutoff**2-f0*t2)/(48*math.pi**2)
                gw2=48*math.pi**2/(n*f0)
                variance=t4-t2*t2/n
                delta=n*cutoff**4/(2*math.pi**2)*(f4-f2*f2/f0)
                close(V-6*gw2*K*K,f0*variance/(8*math.pi**2)+delta)
                assert K > 0 and V >= 6*gw2*K*K
                matrix_rows.append({'N':N,'wq':wq,'wl':wl,'t':t,
                    'K':round(K,9),'V_over_Meff4':round(V/(4*K*K),9),
                    'joint_lower_bound':round(1.5*gw2,9)})
    checks.append('full_Higgs_Majorana_matrices_obey_joint_vacuum_bound')

    # Exact one-generation Higgs radial coefficients, M=0 for this comparison.
    radial_rows=[]
    for N,wq,wl in ((3,F(1),F(1)),(3,F(1),F(2)),(5,F(2),F(3))):
        for ys in ((1,1,1,1),(1,2,3,4)):
            yv,ye,yu,yd=map(F,ys)
            a=wl*(yv*yv+ye*ye)+N*wq*(yu*yu+yd*yd)
            b=wl*(yv**4+ye**4)+N*wq*(yu**4+yd**4)
            n=8*(N*wq+wl)
            assert a*a <= n*b/4
            # lambda/gw^2; pi^2 and f0 cancel in the same normalization.
            ratio=n*b/(24*a*a)
            assert ratio >= F(1,6)
            assert (ratio == F(1,6)) == (len(set(ys))==1)
            radial_rows.append({'N':N,'wq':str(wq),'wl':str(wl),'Yukawas':list(ys),
                               'lambda_over_gw_squared':str(ratio)})
    checks.append('canonical_radial_Higgs_quartic_lower_bound_and_equality')

    # Nonempty positive-K coefficient domain at a flat-background radial minimum.
    # This is deliberately NOT claimed to solve coupled scalar and metric equations.
    N,wq,wl=3,F(1),F(2)
    a,b,n=F(85),F(1045),F(40)
    f0,f2,f4=F(1),F(1,2),F(1,2)
    s=a/b  # t^2, zero of the derivative in s of the flat radial potential
    weights=[4*wl,4*wl,4*N*wq,4*N*wq]
    xs=[s*j*j for j in (1,2,3,4)]
    _,_,_,v,k,_,_=reduced_coefficients(weights,xs,f0,f2,f4)
    assert -2*a*f2+f0*b*s == 0 and f0*b > 0 and k > 0 and v > 0
    flat_min={'t_squared':str(s),'48_pi_squared_K':str(k),
              '8_pi_squared_V':str(v),'coupled_gravity_solution_claimed':False}
    checks.append('flat_radial_minimum_can_have_positive_effective_EH_coefficient')

    # Positive smooth cutoff f=(1+10u^2) exp(-u^2), NOT monotone.
    # Phi=0 is a coefficient counterexample, not a stable vacuum assertion.
    n,t2,t4,v,k,_,gap=reduced_coefficients([F(32)],[F(0)],F(1),F(11,2),F(21,2))
    assert gap == F(-79,4) and k > 0 and v > 0 and v < k*k/n
    nonmonotone={'f0':'1','f2':'11/2','f4':'21/2','f4_minus_f2_squared_over_f0':str(gap),
                 'ratio_of_V_to_claimed_monotone_bound':str(v/(k*k/n)),
                 'field_stationarity_or_observed_vacuum_claimed':False}
    checks.append('positive_nonmonotone_cutoff_disproves_unconditional_version')

    old=json.loads((HERE/'joint_gauge_matter_constraints_results.json').read_text('utf8'))
    inverse_gw2=old['historical_running']['matched_inverse_couplings'][1]
    bound=1.5/inverse_gw2
    assert .48 < bound < .49
    checks.append('saved_historical_gauge_benchmark_only_no_new_RG_fit')

    names=['research_note_532.md','shared_trace_unimodularity.py',
           'shared_trace_unimodularity_results.json','joint_gauge_matter_constraints_results.json',
           'unified_physics_condition_ledger_532.md']
    return dict(round=533,tests_run=len(checks),failures=0,errors=0,checks=checks,
        exact_rational_cases=rational_cases,layer_cake_examples=moment_rows,
        smooth_cutoff_examples=smooth_rows,weak_generator_examples=gauge_rows,
        full_matrix_examples=matrix_rows,radial_Higgs_examples=radial_rows,
        flat_radial_minimum_example=flat_min,nonmonotone_counterexample=nonmonotone,
        historical_benchmark=dict(inverse_gw_squared=inverse_gw2,
            bare_V_over_Meff4_lower_bound=round(bound,12),current_measurement_or_fit=False),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        scope=dict(supplied_four_dimensional_Euclidean_geometry=True,
            positive_nonincreasing_cutoff_is_extra_assumption=True,
            same_weighted_action_truncated_through_a4=True,
            no_independent_vacuum_or_EH_or_YM_counterterm=True,
            known_moment_inequality_reused=True,
            coupled_metric_scalar_vacuum_solved=False,
            observed_cosmological_constant_excluded=False,
            full_spectral_action_or_cognition_ruled_out=False,
            spacetime_dimension_or_GR_derived=False))


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
