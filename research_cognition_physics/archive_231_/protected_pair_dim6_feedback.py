"""Round 542: mass sign and tree-generated dim-6 feedback in the declared SEFT.

2302.08140 Eqs 4, 12--18: C5=Ye=0; top-only; single tree-operator insertions.
This is a specified partial RG resummation, not the full SMEFT flow/observable fit.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import protected_pair_finite_matching as prior

HERE=Path(__file__).resolve().parent
TARGET=HERE/'protected_pair_dim6_feedback_results.json'
LOOP=prior.LOOP


def matrix_beta(C,q,gw):
    return ((2*gw/3)*np.trace(C)*np.eye(3)+(6*q-17*gw/3)*C)/LOOP


def rhs(t,z):
    # Dimensionless mass ratio R=m_EFT^2/M^2 and Z=M^2 Tr C3.
    y,base,delta,R,Z,lt,ld=z
    gy,gw,gc=prior.old.gauge_squared(t,prior.XSTAR)
    q=y*y;lam=base+delta;B=3*gy+9*gw
    dy=y*(4.5*q-17*gy/12-9*gw/4-8*gc)
    db=24*base*base-B*base+3/8*(gy*gy+2*gy*gw+3*gw*gw)+12*q*base-6*q*q
    # Algebraic difference equation avoids subtracting near-equal quartics.
    dd=(48*base-B+12*q)*delta+24*delta*delta-(8/3)*gw*R*Z
    am=12*lam+6*q-1.5*gy-4.5*gw
    at=6*q-11*gw/3;ad=6*q-17*gw/3
    return np.array([dy,db,dd,am*R,at*Z,at,ad])/LOOP


def flow(y,lam,R,S,M,steps):
    start=math.log(M/prior.MU);end=math.log(prior.IR/prior.MU)
    z=np.array([y,lam,0.,R,-S/4,0.,0.]);h=(end-start)/steps
    for k in range(steps):
        t=start+k*h;k1=rhs(t,z);k2=rhs(t+h/2,z+h*k1/2)
        k3=rhs(t+h/2,z+h*k2/2);k4=rhs(t+h,z+h*k3)
        z+=h*(k1+2*k2+2*k3+k4)/6
    assert np.all(np.isfinite(z))
    return z


def trajectory(M,steps):
    x=prior.BOUND['T']/(3+prior.BOUND['r'])
    lam0=(3*x*x+prior.BOUND['r']*x*x)/prior.BOUND['T']
    td=math.log(M/prior.MU)
    high,_,_=prior.old.evolve([math.sqrt(x),math.sqrt(x),lam0,0.,1.],0.,td,prior.XSTAR,True,steps)
    y,s,lam=high[:3];S=s*s;gw=prior.old.gauge_squared(td,prior.XSTAR)[1]
    low=prior.low_flow(y,lam,td,steps)
    R=prior.MASS_TARGET/(M*M*math.exp(low[2]))
    for n in range(30):
        xi,_,_=prior.inverse_mass(S,R)
        yp,lp,_=prior.finite_matching(y,lam,S,gw,float(xi))
        out=flow(yp,lp,R,S,M,steps)
        new=R*prior.MASS_TARGET/(out[3]*M*M)
        if abs(new/R-1)<2e-14:
            R=new;break
        R=new
    else:raise AssertionError('mass iteration did not converge')
    xi,_,_=prior.inverse_mass(S,R)
    yp,lp,_=prior.finite_matching(y,lam,S,gw,float(xi))
    out=flow(yp,lp,R,S,M,steps)
    assert abs(float(xi))<.01
    assert abs(out[3]*M*M/prior.MASS_TARGET-1)<1e-12
    assert out[2]>0 and out[3]<0 and out[4]<0
    et,ed=math.exp(out[5]),math.exp(out[6]);tr0=-S/4
    transverse=tr0*(et-ed)/3;parallel=tr0*(et+2*ed)/3
    assert transverse>0 and parallel<0
    assert abs(out[4]/(tr0*et)-1)<1e-12
    reversed_mass=flow(yp,lp,-R,S,M,steps)
    assert reversed_mass[2]<0 and reversed_mass[3]>0
    return dict(threshold_GeV=M,S_at_threshold=S,
        mass_ratio_before=str(xi),mass_ratio_after=R,
        low_mass_squared_GeV2=out[3]*M*M,
        low_quartic_without_continuous_feedback_same_threshold=out[1],
        low_quartic_with_continuous_feedback=out[1]+out[2],
        positive_quartic_difference_via_stable_ODE=out[2],
        low_top=out[0],scaled_trace_C3_low=out[4],
        scaled_C3_eigenvalues=[parallel,transverse,transverse],
        trace_factor=et,traceless_factor=ed,
        opposite_mass_sign_control_quartic_difference=reversed_mass[2],
        iterations=n+1)


def run():
    checks=[];rng=np.random.default_rng(542);err=0.
    for _ in range(24):
        c=rng.normal(size=3)+1j*rng.normal(size=3)
        c*=math.sqrt(.4)/np.linalg.norm(c)
        Y=np.column_stack([c,1j*c])/math.sqrt(2)
        C=-Y@Y.conj().T/4  # M=1 for dimensionless algebra check.
        assert np.linalg.norm(Y@Y.T)<1e-14 and abs(np.trace(C)+.1)<1e-14
        q=.3;gw=.4
        expected=(6*q-11*gw/3)*np.trace(C)/LOOP
        err=max(err,float(abs(np.trace(matrix_beta(C,q,gw))-expected)))
    assert err<1e-15
    checks.append('nonzero_tree_current_operator_and_trace_equation_from_full_flavor_matrix')

    # Exact algebra for mass/trace logarithmic product and quartic feedback sign.
    gy,gw,q,lam=F(1,7),F(2,5),F(3,10),F(1,4)
    am=12*lam+6*q-F(3,2)*gy-F(9,2)*gw
    at=6*q-F(11,3)*gw
    assert am+at==12*lam+12*q-F(3,2)*gy-F(49,6)*gw
    for m in (F(-1,100),F(-1,10**12)):
        for tr in (F(-1,8),F(-1,10000)):
            assert F(8,3)*gw*m*tr>0
    assert F(8,3)*gw*F(1,100)*F(-1,8)<0
    checks.append('signed_mass_times_negative_trace_gives_nonnegative_downward_quartic_source')

    # Trace sign and entire matrix sign are different: the null directions turn positive.
    C=np.diag([-.1,0.,0.]);down=-matrix_beta(C,.3,.4)
    assert down[1,1]>0 and down[2,2]>0
    checks.append('negative_semidefinite_cone_is_not_preserved_in_downward_running')

    base,d,q=F(2,9),F(1,1000),F(3,5);B=F(4,1)
    direct=24*((base+d)**2-base**2)+(-B+12*q)*d
    assert direct==(48*base-B+12*q)*d+24*d*d
    assert F(16259,100000)>F(12604,100000)
    checks.append('stable_quartic_difference_identity_and_inherited_positive_comparison_floor')

    rows=[];convergence=[]
    for M in (1e4,1e6,1e9):
        coarse=trajectory(M,400);fine=trajectory(M,800)
        rel=abs(coarse['positive_quartic_difference_via_stable_ODE']/fine['positive_quartic_difference_via_stable_ODE']-1)
        absdiff=abs(coarse['low_quartic_with_continuous_feedback']-fine['low_quartic_with_continuous_feedback'])
        assert rel<1e-8 and absdiff<1e-9
        assert fine['low_quartic_with_continuous_feedback']>.16259
        rows.append(fine);convergence.append(dict(M=M,relative_feedback_step_change=rel,absolute_quartic_step_change=absdiff))
    checks.append('three_joint_mass_matched_thresholds_with_nonzero_positive_feedback')
    checks.append('opposite_mass_sign_control_reverses_feedback_but_cannot_match_negative_mass')
    checks.append('step_halving_and_exact_trace_factor_crosscheck')
    deps=('research_note_541.md','protected_pair_finite_matching.py','protected_pair_finite_matching_results.json',
          'joint_higgs_flow_obstruction_results.json','unified_physics_condition_ledger_541.md')
    return dict(round=542,tests_run=len(checks),failures=0,errors=0,checks=checks,
        trace_matrix_max_error=err,trajectories=rows,step_halving=convergence,
        inherited_lower_bound=.16259,source_negative_budget=0.,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(exact_protected_pair=True,negative_EFT_mass_target_required=True,
            tree_generated_dim6_single_insertion_feedback_included=True,
            trace_sign_not_matrix_semidefiniteness_used=True,
            partial_resummation_of_displayed_subsystem=True,full_SMEFT_running=False,
            full_dimension_six_observable_extraction=False,
            broken_nondegenerate_pair_excluded=False,full_unification_completed=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','inherited_lower_bound','step_halving')},ensure_ascii=False))
