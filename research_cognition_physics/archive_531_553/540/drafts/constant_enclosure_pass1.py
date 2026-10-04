"""Round 540: analytic all-trajectory lower bound in the declared leading-log branch.

Exact rational downward Euler certification; samples are diagnostics, not the proof.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import protected_pair_rg as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_higgs_flow_obstruction_results.json'
A=F(82,25); C=F(53,200); Z0=F(161,500)
DMAX=F(1,1000000); TIME=F(19,157); STEPS=100; SCALE=10**12


def f(z):return -32*z*z+A*z-C


def exp_upper(x,n=20):
    """Positive rational Taylor sum plus geometric bound on its tail."""
    x=F(x);term=F(1);s=term
    for k in range(1,n+1):term*=x/k;s+=term
    first=term*x/(n+1);ratio=x/(n+2)
    assert 0<=ratio<1
    return s+first/(1-ratio)


def certificate():
    h=TIME/STEPS;z=Z0;values=[z]
    assert 1+h*(-64*Z0+A)>0
    for _ in range(STEPS):
        assert F(17,100)<=z<=Z0 and f(z)<0 and -64*z+A<0
        raw=z+h*f(z);numerator=(raw*SCALE).numerator//(raw*SCALE).denominator
        new=F(numerator,SCALE)
        assert 0<=raw-new<F(1,SCALE) and new<z
        z=new;values.append(z)
    assert z>F(17009,100000)
    return values


def run():
    checks=[];hist,b=old.inputs();rng=np.random.default_rng(540)
    for _ in range(40):
        q,s,lam,D,gy,gw=[F(int(x),100) for x in rng.integers(0,100,size=6)]
        Q=s*s-2*D;B=3*gy+9*gw;C0=F(3,8)*(gy*gy+2*gy*gw+3*gw*gw)
        direct=-24*lam*lam+B*lam-C0-4*(3*q+s)*lam+2*(3*q*q+Q)
        square=-32*lam*lam+B*lam-C0-4*D+6*(q-lam)**2+2*(s-lam)**2
        assert direct==square
    checks.append('exact_completed_square_identity_for_top_and_two_neutrino_columns')

    low=old.gauge_squared(math.log(hist['MZ_GeV']/b['matching_scale_GeV']),hist['matched_inverse_couplings'])
    high=old.gauge_squared(0,hist['matched_inverse_couplings'])
    Bmin=3*low[0]+9*high[1]; Bmax=3*high[0]+9*low[1]
    Cmax=3/8*(high[0]**2+2*high[0]*low[1]+3*low[1]**2)
    # The signs of the inherited one-loop coefficients establish endpoint extrema.
    assert Bmin>3.28 and Bmax<4.31 and Cmax<.26431
    assert Cmax+4*1.7e-6<.265
    T,r=F(str(b['T'])),F(str(b['r'])); w=T/(3+r)
    assert w-2*r*DMAX/T>Z0
    # pi>3.14 is enough for L>157. A finite positive sum certifies e^19>mu*/MZ.
    assert 16*F(314,100)**2>157
    lower_exp19=sum((F(19)**k/F(math.factorial(k)) for k in range(55)),F(0))
    assert lower_exp19>F(str(b['matching_scale_GeV']))/F(str(hist['MZ_GeV']))
    checks.append('conservative_frozen_gauge_and_high_boundary_enclosures')

    det_error=0.
    for _ in range(30):
        Y=rng.normal(size=(3,2))+1j*rng.normal(size=(3,2));G=Y.conj().T@Y
        S=np.trace(G).real;q=.3;gy=.16;gw=.32
        dY=-(1.5*Y@G+(3*q+S-.75*gy-2.25*gw)*Y)
        dG=dY.conj().T@Y+Y.conj().T@dY
        exact=3*gy+9*gw-12*q-7*S
        det_error=max(det_error,float(abs(np.trace(np.linalg.solve(G,dG)).real-exact)))
    assert det_error<1e-11
    growth=exp_upper(F(431,100)*TIME)
    assert growth<F(17,10)
    checks.append('full_matrix_determinant_growth_and_rational_small_breaking_budget')

    vals=certificate();terminal=vals[-1]
    checks.append('one_hundred_exact_downward_rounded_Euler_subsolution_steps')

    k=math.sqrt(128*float(C)-float(A)**2)/64;center=float(A)/64
    exact=center+k*math.tan(math.atan((float(Z0)-center)/k)-32*k*float(TIME))
    assert exact>float(terminal)>.17009
    z=float(Z0);h=float(TIME)/2000
    for _ in range(2000):
        k1=float(f(z));k2=float(f(z+h*k1/2));k3=float(f(z+h*k2/2));k4=float(f(z+h*k3))
        z+=h*(k1+2*k2+2*k3+k4)/6
    assert abs(z-exact)<1e-12
    checks.append('independent_closed_Riccati_expression_and_numerical_crosscheck')

    samples=[]
    for x in (0.,.1,.2,b['gw_squared'],.5,b['T']/3):
        # Avoid a negative last-bit neutrino norm in the endpoint diagnostic.
        if x==b['T']/3:x=np.nextafter(x,0.)
        for dec in (1e3,1e6,1e9):
            row=old.trajectory(x,dec,300)
            assert row['lambda_at_MZ']>float(terminal)
            samples.append(dict(x=x,threshold=dec,lambda_MZ=row['lambda_at_MZ']))
    checks.append('diagnostic_trajectories_cover_varied_top_and_threshold_without_serving_as_proof')

    amplification=exp_upper(A*TIME)
    assert amplification<F(3,2)
    target=F(12604,100000);gap=F(17009,100000)-target
    minimum_budget=gap/F(3,2)
    assert minimum_budget>F(2936,100000)
    checks.append('conditional_negative_correction_budget_with_positive_first_crossing')

    deps=('research_note_539.md','rank_two_breaking_bridge_results.json',
        'protected_pair_rg.py','joint_gauge_matter_constraints_results.json',
        'joint_yukawa_higgs_matching_results.json','unified_physics_condition_ledger_539.md')
    return dict(round=540,tests_run=len(checks),failures=0,errors=0,checks=checks,
        assumptions=dict(high_determinant_at_most=str(DMAX),top_only_charged_Yukawa=True,
            at_most_two_active_neutrino_columns=True,tree_continuous_quartic_thresholds=True,
            one_loop_dimension_four_running=True,frozen_round531_gauge_input=True),
        constants=dict(B_min=float(Bmin),B_max=float(Bmax),C_max=float(Cmax),
            initial_lambda_lower=float(w-2*r*DMAX/T),det_growth_upper=float(growth),
            reduced_time_upper=str(TIME),comparison_A=str(A),comparison_C=str(C)),
        determinant_identity_max_error=det_error,
        exact_certificate=dict(scale=SCALE,steps=STEPS,terminal_fraction=str(terminal),
            all_step_integer_numerators=[int(v*SCALE) for v in vals],strict_reported_floor=.17009,
            comparison_exact_solution=exact,comparison_RK4_error=abs(z-exact)),
        historical_low_benchmark=dict(scale_GeV=173.34,lambda_MS=float(target),
            source='Buttazzo et al 1307.3536 Eq(55), same quartic convention; historical benchmark, not current fit',
            source_url='https://arxiv.org/html/1307.3536',required_gap=float(gap)),
        conditional_revision_budget=dict(negative_initial_shift_plus_negative_jumps_plus_negative_source_integral=True,
            source_integrated_in_reduced_time=True,amplification_upper=float(amplification),
            safe_amplification=1.5,necessary_budget_exact=str(minimum_budget),
            necessary_budget_float=float(minimum_budget),not_a_sufficient_repair=True),
        numerical_samples_not_used_as_proof=samples,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(all_top_and_threshold_choices_in_declared_branch=True,
            covers_round539_examples=True,finite_matching_or_higher_loops_excluded_from_theorem=True,
            full_physical_model_or_cognitive_axioms_disproved=False,
            current_global_experimental_fit_performed=False,
            full_unification_or_GR_derived=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','constants','conditional_revision_budget')},ensure_ascii=False))
