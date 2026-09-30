"""Round 538: the same protected pair and spectral boundary under RG.

One-loop dimension-four MS running plus tree-level simultaneous decoupling.
No finite one-loop thresholds, current data fit, or nonzero kappa running.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'protected_pair_rg_results.json'
LOOP=16*math.pi**2


def inputs():
    a=json.loads((HERE/'joint_gauge_matter_constraints_results.json').read_text('utf8'))['historical_running']
    b=json.loads((HERE/'joint_yukawa_higgs_matching_results.json').read_text('utf8'))['historical_boundary']
    return a,b


def gauge_squared(t,xstar):
    return 1/(np.array(xstar)+np.array([-41/6,19/6,7])*t/(8*math.pi**2))


def rhs(t,state,xstar,active):
    u1,u2,u3=gauge_squared(t,xstar)
    y,s,l,integral,ratio=state
    n=s*s if active else 0.0
    return np.array([y*(4.5*y*y+n-17*u1/12-9*u2/4-8*u3),
        s*(3*y*y+2.5*n-3*u1/4-9*u2/4) if active else 0.0,
        24*l*l-(3*u1+9*u2)*l+3/8*(u1*u1+2*u1*u2+3*u2*u2)
            +4*(3*y*y+n)*l-2*(3*y**4+n*n),
        n,ratio*n])/LOOP


def evolve(state,t0,t1,xstar,active,steps):
    a=np.array(state,float); h=(t1-t0)/steps
    mins=a.copy(); maxs=a.copy()
    for j in range(steps):
        t=t0+j*h
        k1=rhs(t,a,xstar,active)
        k2=rhs(t+h/2,a+h*k1/2,xstar,active)
        k3=rhs(t+h/2,a+h*k2/2,xstar,active)
        k4=rhs(t+h,a+h*k3,xstar,active)
        a+=h*(k1+2*k2+2*k3+k4)/6
        mins=np.minimum(mins,a);maxs=np.maximum(maxs,a)
        assert np.all(np.isfinite(a)) and min(a[:3])>=0
    return a,mins,maxs


def trajectory(x,decoupling,steps=800,omit_neutrino=False):
    old,b=inputs(); mu=b['matching_scale_GeV']; r,T=b['r'],b['T']
    xstar=old['matched_inverse_couplings']; end=math.log(old['MZ_GeV']/mu)
    S=(T-3*x)/r; lam=(3*x*x+r*S*S)/T
    assert S>=0 and old['MZ_GeV']<decoupling<mu
    initial=[math.sqrt(x),math.sqrt(S),lam,0.0,1.0]
    td=math.log(decoupling/mu)
    high,minimum,maximum=evolve(initial,0,td,xstar,not omit_neutrino,steps)
    low,lo_min,lo_max=evolve(high,td,end,xstar,False,steps)
    # The carried s below threshold is a frozen diagnostic, not an EFT coupling.
    Mstar=decoupling*math.exp(-high[3])
    assert Mstar<=mu and math.isclose(math.log(high[4]),high[3],abs_tol=3e-12)
    return dict(top_squared_at_matching=x,neutrino_squared_at_matching=S,lambda_at_matching=lam,
        decoupling_scale_GeV=decoupling,Majorana_at_matching_GeV=Mstar,
        Majorana_ratio_at_threshold=high[4],log_Majorana_running=high[3],
        top_at_threshold=high[0],neutrino_at_threshold=high[1],lambda_at_threshold=high[2],
        top_at_MZ=low[0],lambda_at_MZ=low[2],
        sampled_min_lambda=min(minimum[2],lo_min[2]),
        sampled_max_dimensionless_coupling=max(maximum[0],maximum[1],maximum[2],lo_max[0],lo_max[2]),
        neutrino_omitted_above_threshold=omit_neutrino)


def run():
    old,b=inputs();xstar=old['matched_inverse_couplings']; checks=[]
    rng=np.random.default_rng(538);matrix_error=0.0
    for _ in range(12):
        c=rng.normal(size=3)+1j*rng.normal(size=3)
        c*=0.6/np.linalg.norm(c); s2=float(np.vdot(c,c).real)
        Y=np.column_stack([c,1j*c,np.zeros(3)])/math.sqrt(2)
        M=np.diag([1.3,1.3,0.7]).astype(complex);H=Y.conj().T@Y
        dM=H.T@M+M@H
        expected=np.diag([1.3*s2,1.3*s2,0])
        matrix_error=max(matrix_error,float(np.max(abs(dM-expected))))
        y=0.5;u1,u2,_=gauge_squared(0,xstar)
        a=3*y*y+s2-3*u1/4-9*u2/4
        dY=1.5*Y@H+a*Y
        factor=3*y*y+2.5*s2-3*u1/4-9*u2/4
        assert np.allclose(dY,factor*Y,atol=1e-14)
        inv=np.linalg.inv(M)
        C=Y@inv@Y.T
        dC=dY@inv@Y.T+Y@inv@dY.T-Y@inv@dM@inv@Y.T
        assert np.linalg.norm(C)<1e-15 and np.linalg.norm(dC)<1e-14
    assert matrix_error<1e-14
    checks.append('independent_matrix_beta_reduction_and_zero_Weinberg_tangency')

    # All arithmetic exact: lambda_src=6*lambda in 1204.0328 Eq(8.12).
    y,s,l,u1,u2=map(F,['0.5','0.6','0.2','0.16','0.32'])
    src=6*l
    source_beta=(4*src*src-(3*u1+9*u2)*src+F(9,4)*(u1*u1+2*u1*u2+3*u2*u2)
                 +4*(3*y*y+s*s)*src-12*(3*y**4+s**4))/6
    canonical=(24*l*l-(3*u1+9*u2)*l+F(3,8)*(u1*u1+2*u1*u2+3*u2*u2)
               +4*(3*y*y+s*s)*l-2*(3*y**4+s**4))
    assert source_beta==canonical
    checks.append('exact_quartic_normalization_conversion')

    S=b['gw_squared'];c=np.array([math.sqrt(S),0,0],complex)
    Y=np.column_stack([c,1j*c,np.zeros(3)])/math.sqrt(2)
    individual=[np.outer(Y[:,j],Y[:,j]) for j in (0,1)]
    assert np.linalg.norm(individual[0])>0 and np.linalg.norm(sum(individual))<1e-15
    # Every common loop weight multiplies the same zero column sum.
    for weight in (0.1,1.0,7.0):assert np.linalg.norm(weight*sum(individual))<1e-14
    checks.append('degenerate_pair_simultaneous_tree_threshold_and_nonzero_single_column_control')

    initial=np.array([math.sqrt(S),math.sqrt(S),S,0,1])
    active=rhs(0,initial,xstar,True); absent=rhs(0,initial,xstar,False)
    assert math.isclose((active[2]-absent[2])*LOOP,2*S*S,abs_tol=1e-14)
    assert math.isclose((active[0]-absent[0])*LOOP,initial[0]*S,abs_tol=1e-14)
    loop_effect=dict(Weinberg_coefficient='zero',neutrino_contribution_to_beta_lambda=active[2]-absent[2],
                     neutrino_contribution_to_beta_top=active[0]-absent[0])
    checks.append('zero_Weinberg_does_not_cancel_Hermitian_loop_invariants')

    rows=[]; convergence=[]
    for dec in (1e3,1e6,1e9,1e10):
        coarse=trajectory(S,dec,400); fine=trajectory(S,dec,800)
        error=max(abs(coarse[k]-fine[k]) for k in ('top_at_MZ','lambda_at_MZ','log_Majorana_running'))
        assert error<1e-9
        rows.append(fine);convergence.append(dict(decoupling=dec,max_step_halving_change=error))
    checks.append('same_boundary_leading_log_trajectories_with_self_consistent_running_mass_threshold')

    absent=trajectory(S,1e3,800,True)
    deltas=[dict(decoupling=row['decoupling_scale_GeV'],
                 top_change=row['top_at_MZ']-absent['top_at_MZ'],
                 lambda_change=row['lambda_at_MZ']-absent['lambda_at_MZ']) for row in rows]
    assert all(d['top_change']<0 and d['lambda_change']<0 for d in deltas)
    assert rows[0]['lambda_at_MZ']<rows[-1]['lambda_at_MZ']
    checks.append('explicit_neutrino_omission_changes_joint_low_scale_output_for_tested_thresholds')

    family=[]
    for x in (0.1,0.2,S,0.4,0.5):
        row=trajectory(x,1e6,800)
        family.append(row)
        assert math.isclose(3*x+b['r']*row['neutrino_squared_at_matching'],b['T'],abs_tol=1e-12)
        assert row['sampled_min_lambda']>0 and row['sampled_max_dimensionless_coupling']<2
    checks.append('several_common_boundary_points_evolve_without_independent_Higgs_retuning')

    # Gauge relation must be a boundary condition, not imposed throughout the flow.
    slope=F(-41,6)-3*F(19,6)+F(4,3)*7
    assert slope==-7
    t=-1.0;gs=gauge_squared(t,xstar);inv=1/gs
    delta0=xstar[0]-3*xstar[1]+F(4,3)*xstar[2]
    assert abs((inv[0]-3*inv[1]+4*inv[2]/3)-float(delta0)-float(slope)*t/(8*math.pi**2))<1e-14
    checks.append('inherited_gauge_plane_is_not_reimposed_at_every_scale')

    def tidy(x):
        if isinstance(x,(float,np.floating)):return float(format(float(x),'.13g'))
        if isinstance(x,list):return [tidy(v) for v in x]
        if isinstance(x,dict):return {k:tidy(v) for k,v in x.items()}
        return x
    deps=('research_note_537.md','seesaw_boundary_compatibility_results.json',
          'joint_gauge_matter_constraints_results.json','joint_yukawa_higgs_matching_results.json',
          'unified_physics_condition_ledger_537.md')
    return tidy(dict(round=538,tests_run=len(checks),failures=0,errors=0,checks=checks,
        matrix_beta_max_error=matrix_error,local_loop_effect=loop_effect,
        common_initial_boundary=dict(mu_GeV=b['matching_scale_GeV'],r=b['r'],top=math.sqrt(S),
            nu_singular_value=math.sqrt(S),quartic=S,low_scale_GeV=old['MZ_GeV']),
        threshold_trajectories=rows,step_halving=convergence,wrong_neutrino_omission_control=absent,
        difference_from_omission=deltas,common_boundary_family_samples=family,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(exact_protected_two_column_branch=True,top_only_charged_Yukawa_approximation=True,
            one_loop_dimension_four_MS_running_and_tree_matching=True,
            simultaneous_running_mass_pair_threshold=True,
            fixed_historical_gauge_input_not_current_fit=True,
            finite_one_loop_threshold_corrections_included=False,
            nonzero_Weinberg_RG_normalization_resolved=False,
            nonzero_neutrino_observations_or_low_energy_pole_masses_fitted=False,
            numerical_samples_not_global_parameter_exclusion=True,
            GR_or_spacetime_or_full_unification_derived=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        body=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(body)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','local_loop_effect','threshold_trajectories')},ensure_ascii=False))
