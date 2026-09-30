"""Round 536: constructive joint gauge/Yukawa/quartic boundary conditions.

No RG integration or current experimental fit. The gauge benchmark is read
from the frozen round 531 result. lambda multiplies chi**4/4.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_yukawa_higgs_matching_results.json'


def moments(ys,N,r):
    # Rows are generations; columns are nu,e,u,d; use singular values.
    weights=np.array([r,r,N,N])
    v=np.array(ys,dtype=float)
    return float(np.sum(weights*v*v)),float(np.sum(weights*v**4))


def interval(T,K2,K4,r,ng):
    assert T>0 and r>0 and ng>=1
    remainder=T-K2
    if remainder < 0: return None
    return ((K4+remainder*remainder/(ng*r))/T,
            (K4+remainder*remainder/r)/T)


def complete(T,K2,K4,r,ng,lam):
    bounds=interval(T,K2,K4,r,ng)
    if bounds is None: return None
    lo,hi=bounds
    tol=1e-12*max(1,abs(lo),abs(hi))
    if lam<lo-tol or lam>hi+tol: return None
    S=(T-K2)/r
    Q=(lam*T-K4)/r
    if ng==1: return [S]
    excess=max(0.0,Q-S*S/ng)
    first=S/ng+math.sqrt((ng-1)*excess/ng)
    rest=(S-first)/(ng-1)
    assert rest>=-tol
    return [first]+[max(0.0,rest)]*(ng-1)


def run():
    checks=[]; normalization=[]
    for N in (3,5):
        for ng in (1,3):
            for r in (F(1),F(7,3)):
                g2=F(1,3); T=ng*(N+r)*g2/3
                alpha=1/T  # f0/(2*pi^2)
                raw=[F(j+1,10) for j in range(4*ng)]
                weights=[r,r,F(N),F(N)]*ng
                a=sum(w*y*y for w,y in zip(weights,raw))
                b=sum(w*y**4 for w,y in zip(weights,raw))
                physical_squared=[y*y/(alpha*a) for y in raw]
                K2=sum(w*y for w,y in zip(weights,physical_squared))
                K4=sum(w*y*y for w,y in zip(weights,physical_squared))
                lam=b/(alpha*a*a)
                assert K2==T and K4/K2==lam
                # Overall raw Dirac scale is not an extra physical Yukawa dial.
                assert (25*b)/(alpha*(5*a)**2)==lam
                normalization.append(dict(N=N,ng=ng,r=str(r),T=str(T),lambda_exact=str(lam)))
    checks.append('exact_common_canonical_normalization_and_raw_scale_cancellation')

    assert F(3)*(3+1)*F(1,3)/3==4*F(1,3)
    for rho in (F(0),F(1),F(2)):
        g2=F(1,3); x=4*g2/(3+rho*rho)
        quartic=(3+rho**4)*x*x/(4*g2)
        assert quartic==4*g2*(3+rho**4)/(3+rho*rho)**2
        assert 6*quartic==24*g2*(3+rho**4)/(3+rho*rho)**2
    checks.append('unweighted_literature_limit_and_lambda_h4_over_24_conversion')

    moment_cases=0; worst=0.0
    for ng in (1,2,3,5):
        for r in (F(1,2),F(1),F(7,3)):
            T,K2,K4=F(2),F(3,5),F(1,5)
            lo,hi=interval(T,K2,K4,r,ng)
            for fraction in (F(0),F(1,5),F(1,2),F(1)):
                lam=lo+fraction*(hi-lo)
                z=complete(float(T),float(K2),float(K4),float(r),ng,float(lam))
                assert z is not None and min(z)>=0
                err=max(abs(float(r)*sum(z)+float(K2)-float(T)),
                        abs(float(r)*sum(x*x for x in z)+float(K4)-float(lam*T)))
                worst=max(worst,err); assert err<1e-11
                moment_cases+=1
    checks.append('constructive_necessary_sufficient_neutrino_completion_including_rank_boundaries')

    assert complete(2,3,1,1,3,1) is None
    lo,hi=interval(2,0.5,0.2,2,3)
    assert complete(2,0.5,0.2,2,3,lo-0.01) is None
    assert complete(2,0.5,0.2,2,3,hi+0.01) is None
    assert complete(2,2,0.5,2,3,0.25)==[0.0]*3
    assert complete(2,2,0.5,2,3,0.26) is None
    checks.append('negative_budget_outside_fourth_moment_interval_and_zero_residual_controls')

    rng=np.random.default_rng(536); matrix_err=0.0
    singular=np.array([0.2,0.4,0.8])
    for _ in range(9):
        u,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        v,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
        Y=u@np.diag(singular)@v.conj().T
        A=Y.conj().T@Y
        err=max(abs(np.trace(A)-sum(singular**2)),abs(np.trace(A@A)-sum(singular**4)))
        matrix_err=max(matrix_err,float(err))
    assert matrix_err<1e-12
    checks.append('flavour_unitaries_preserve_the_two_boundary_moments_not_the_full_Majorana_spectrum')

    exact_controls=[]
    for r in (F(1,2),F(1),F(7,3),F(5)):
        T=F(5,3); xmin=T/(3*(1+r))
        minimum=(3*xmin*xmin+(T-3*xmin)**2/(3*r))/T
        assert minimum==T/(3*(1+r))
        maximum=T/min(F(3),r)
        for x in (F(0),xmin,T/6,T/3):
            lo,hi=interval(T,3*x,3*x*x,r,3)
            assert minimum<=lo<=hi<=maximum
            gap=3*(1+r)/(r*T)*(x-xmin)**2
            assert lo-minimum==gap
        exact_controls.append(dict(r=str(r),minimum=str(minimum),maximum=str(maximum)))
    checks.append('sharp_top_plus_three_neutrino_envelope_and_exact_square_remainder')

    old=json.loads((HERE/'joint_gauge_matter_constraints_results.json').read_text('utf8'))['historical_running']
    xY,xw,xc=old['matched_inverse_couplings']
    r=4*xw/xc-3; gw2=1/xw; T=(3+r)*gw2
    assert abs(r-old['positive_weight_ratio_lepton_over_quark'])<2e-9
    assert abs(xY-3*xw+4*xc/3)<2e-9
    f0=2*math.pi**2/T
    reconstructed=[f0*(11+9*r)/(6*math.pi**2),f0*(3+r)/(2*math.pi**2),2*f0/math.pi**2]
    assert max(abs(a-b) for a,b in zip(reconstructed,(xY,xw,xc)))<2e-9
    table=[]
    for yt in (0.0,0.3,0.5,0.7):
        lo,hi=interval(T,3*yt*yt,3*yt**4,r,3)
        table.append(dict(top_Yukawa_input=yt,lambda_min=lo,lambda_max=hi))
    benchmark=dict(source='frozen round 531 historical gauge-only one-loop diagnostic; not a current fit',
        matching_scale_GeV=old['weighted_trace_scale_GeV'],r=r,gw_squared=gw2,T=T,
        all_top_values_lambda_min=T/(3*(1+r)),all_top_values_lambda_max=T/min(3,r),
        maximum_top_Yukawa=math.sqrt(T/3),allowed_intervals=table)
    checks.append('frozen_gauge_benchmark_fixes_same_weight_and_matter_budget_without_refitting')

    yt=0.5; lam=0.2; K2=3*yt*yt; K4=3*yt**4
    z=complete(T,K2,K4,r,3,lam)
    assert z is not None
    ys=np.zeros((3,4)); ys[:,0]=np.sqrt(z); ys[2,2]=yt
    A,B=moments(ys,3,r)
    assert abs(A/T-1)<1e-12 and abs(B/A-lam)<1e-12
    # Choose raw y=Y. Then f0*A=2*pi^2, so canonical conversion is identity.
    assert abs(2*math.pi**2*B/(f0*A*A)-lam)<1e-12
    assert 0.13 < benchmark['all_top_values_lambda_min']
    witness=dict(top_Yukawa=yt,lambda_input=lam,neutrino_Yukawas=list(np.sqrt(z)),
                 reconstructed_second_moment=A,reconstructed_lambda=B/A,
                 diagnostic_lambda_013_excluded_in_top_only_family=True,
                 any_low_energy_or_current_experimental_interpretation=False)
    checks.append('simultaneous_gauge_and_matter_witness_and_strictly_excluded_diagnostic_target')

    def tidy(obj):
        if isinstance(obj,(float,np.floating)):return round(float(obj),12)
        if isinstance(obj,dict):return {k:tidy(v) for k,v in obj.items()}
        if isinstance(obj,list):return [tidy(v) for v in obj]
        return obj
    deps=('research_note_535.md','fixed_volume_scalar_stability_results.json',
          'joint_gauge_matter_constraints_results.json','research_note_532.md','research_note_533.md',
          'unified_physics_condition_ledger_535.md')
    return tidy(dict(round=536,tests_run=len(checks),failures=0,errors=0,checks=checks,
        rational_normalization_cases=normalization,constructive_moment_cases=moment_cases,
        max_completion_residual=worst,max_flavour_matrix_residual=matrix_err,
        exact_top_envelopes=exact_controls,historical_boundary=benchmark,joint_witness=witness,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(ordinary_group_selection_and_weighted_action_branch=True,
            all_ng_right_handed_neutrinos_active_at_matching=True,canonical_fermion_action_supplied=True,
            exact_tree_level_dimension_four_boundary_feasibility=True,
            general_flavour_mixing_only_for_second_and_fourth_Yukawa_moments=True,
            sharp_top_envelope_assumes_other_charged_Yukawas_zero=True,
            RG_evolution_or_thresholds_or_seesaw_observations_checked=False,
            complete_stable_gravitating_vacuum_or_SM_or_cognition_derived=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','historical_boundary','joint_witness')},ensure_ascii=False))
