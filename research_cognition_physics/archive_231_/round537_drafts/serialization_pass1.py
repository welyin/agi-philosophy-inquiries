"""Round 537: optimal tree-level seesaw norm at fixed Yukawa singular values.

This is a common-scale coefficient audit, not RG evolution or an oscillation fit.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'seesaw_boundary_compatibility_results.json'


def witness(z):
    a,b,c=np.sqrt(z)
    return np.array([[a/math.sqrt(2),1j*a/math.sqrt(2),0],
                     [0,0,b],[c/math.sqrt(2),-1j*c/math.sqrt(2),0]],complex)


def lower_norm(z,Mmax):
    assert z[0]>=z[1]>=z[2]>=0 and Mmax>0
    return max(z[1],math.sqrt(z[0]*z[2]))/Mmax


def haar(rng):
    q,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    return q


def run():
    checks=[]; exact=[]; max_err=0.0
    for z in ((1.0,0.3,0.1),(1.0,0.01,0.001),(0.8,0.8,0.8),(2.0,0.0,0.0),(0,0,0)):
        Y=witness(z)
        actual=np.linalg.svd(Y,compute_uv=False)**2
        assert np.allclose(actual,z,rtol=1e-12,atol=1e-14)
        for Mmax in (0.7,1.0,2.0):
            C=Y@Y.T/Mmax
            value=float(np.linalg.norm(C,2)); target=lower_norm(z,Mmax)
            max_err=max(max_err,abs(value-target))
            assert abs(value-target)<1e-12
            exact.append(dict(z=list(z),Mmax=Mmax,norm=value,analytic=target))
    checks.append('explicit_complex_flavour_construction_attains_both_singular_value_bounds')

    rng=np.random.default_rng(537); ratios=[]
    for _ in range(80):
        s=np.sort(rng.uniform(0.02,1.2,3))[::-1]
        Y=haar(rng)@np.diag(s)@haar(rng)
        masses=rng.uniform(0.2,2.0,3); Mmax=2.0
        A=Y@np.diag(1/np.sqrt(masses)); sa=np.linalg.svd(A,compute_uv=False)
        assert np.all(sa+1e-12>=s/math.sqrt(Mmax))
        value=float(np.linalg.norm(A@A.T,2))
        bound=max(sa[1]**2,sa[0]*sa[2])
        assert value+1e-12>=bound>=lower_norm(s*s,Mmax)-1e-12
        ratios.append(value/lower_norm(s*s,Mmax))
    checks.append('nondegenerate_Majorana_and_generic_complex_flavour_lower_bound_diagnostics')

    rational=[]
    for T in (F(1),F(2)):
        for r in (F(1,2),F(1),F(7,3)):
            for x in (F(0),T/12,T/6,T/3):
                S=(T-3*x)/r
                rank1=(3*x*x+r*S*S)/T
                assert rank1-T/(3+r)==3*(3+r)/(r*T)*(x-T/(3+r))**2
                rational.append(dict(T=str(T),r=str(r),x=str(x),rank1_lambda=str(rank1)))
    checks.append('zero_Weinberg_rank_one_boundary_and_exact_quartic_floor')

    finite=[]
    for z in ((F(1),F(1,100),F(1,10000)),(F(1,10),F(1,100),F(1,200)),(F(0),F(0),F(0))):
        eps=max(z[1],float(z[0]*z[2])**0.5)
        S=sum(z); Q=sum(a*a for a in z)
        assert float(Q)>=float(S*S)-4*float(S)*eps-1e-15
        r=F(7,3); x=F(1,4); T=3*x+r*S
        lam=(3*x*x+r*Q)/T
        assert float(lam)>=float(T/(3+r))-4*eps-1e-14
        finite.append(dict(z=list(map(str,z)),epsilon=eps,lambda_value=float(lam),
                           necessary_floor=float(T/(3+r))-4*eps))
    checks.append('finite_Weinberg_budget_implies_joint_quartic_lower_bound')

    old=json.loads((HERE/'joint_yukawa_higgs_matching_results.json').read_text('utf8'))
    base=old['historical_boundary']; joint=old['joint_witness']
    Mmax=base['matching_scale_GeV']; z=sorted([x*x for x in joint['neutrino_Yukawas']],reverse=True)
    v_reference=246.0; delta_eV_input=0.1
    eta=2*(delta_eV_input*1e-9)/v_reference**2
    eps=Mmax*eta
    minC=lower_norm(z,Mmax)
    conversion=minC*v_reference**2/2*1e9
    assert minC>eta and joint['lambda_input']<base['gw_squared']-4*eps
    diagnostic=dict(source='frozen 536 boundary; eta is an added common-scale diagnostic input, not a run experimental bound',
        Mmax_GeV=Mmax,eta_per_GeV=eta,eta_times_Mmax=eps,
        reference_v_GeV=v_reference,reference_mass_budget_eV=delta_eV_input,
        previous_witness_z=z,minimum_Weinberg_norm_per_GeV=minC,
        formal_mass_conversion_eV=conversion,
        previous_witness_excluded_for_this_budget=True,
        top_family_lambda_necessary_floor=base['gw_squared']-4*eps,
        rank_one_zero_budget_lambda_min=base['gw_squared'])
    checks.append('previous_joint_witness_fails_added_small_Weinberg_budget_even_with_complex_cancellations')

    T,r=base['T'],base['r']; x=T/(3+r); S=(T-3*x)/r
    Y=witness((S,0.0,0.0)); C=Y@Y.T/Mmax
    lam=(3*x*x+r*S*S)/T
    assert np.linalg.norm(C,2)<1e-24 and abs(lam-base['gw_squared'])<1e-12
    assert abs(3*x+r*np.trace(Y@Y.conj().T).real-T)<1e-12
    # Real diagonal Y has the same moments but cannot make C zero.
    Yreal=np.diag([math.sqrt(S),0,0])
    assert np.linalg.norm(Yreal@Yreal.T/Mmax,2)>eta
    surviving=dict(top_Yukawa=math.sqrt(x),neutrino_singular_values=[math.sqrt(S),0,0],
        lambda_value=lam,analytic_Weinberg_coefficient='zero',
        generic_tree_level_existence=True,observed_nonzero_oscillation_spectrum_reproduced=False)
    checks.append('surviving_same_gauge_boundary_with_rank_one_complex_cancellation_and_real_control')

    # Exact rank, not an ill-conditioned SVD of a physical 10^10 GeV matrix:
    # if D D^T=0 and MR invertible, block elimination gives rank(full)=3.
    D=np.array([[1,1j,0],[2,2j,0],[0,0,0]],complex)
    assert np.array_equal(D@D.T,np.zeros((3,3),complex))
    heavy=2*np.eye(3)
    full=np.block([[np.zeros((3,3)),D],[D.T,heavy]])
    assert np.linalg.matrix_rank(full,tol=1e-12)==3
    checks.append('exact_zero_cancellation_has_three_massless_full_matrix_modes')

    # A small nonzero third column gives rank-two Yukawa, a controlled nonzero C.
    perturb=eps/4
    Y=witness((S-perturb,perturb,0))
    C=Y@Y.T/Mmax
    assert abs(np.trace(Y@Y.conj().T).real-S)<1e-12
    assert 0<np.linalg.norm(C,2)<eta
    lam_pert=(3*x*x+r*((S-perturb)**2+perturb**2))/T
    assert lam_pert>=base['gw_squared']-4*eps
    checks.append('small_nonzero_Weinberg_example_keeps_second_moment_and_stays_near_rank_one_quartic')

    def tidy(obj):
        # Preserve tiny dimensionful coefficients with significant digits.
        if isinstance(obj,(float,np.floating)): return float(format(float(obj),'.13g'))
        if isinstance(obj,list):return [tidy(v) for v in obj]
        if isinstance(obj,dict):return {k:tidy(v) for k,v in obj.items()}
        return obj
    deps=('research_note_536.md','joint_yukawa_higgs_matching.py','joint_yukawa_higgs_matching_results.json',
          'unified_physics_condition_ledger_536.md')
    return tidy(dict(round=537,tests_run=len(checks),failures=0,errors=0,checks=checks,
        exact_norm_constructions=exact,max_attainment_error=max_err,generic_random_cases=len(ratios),
        minimum_generic_ratio_to_bound=min(ratios),rational_rank_one_controls=rational,
        finite_budget_controls=finite,common_scale_diagnostic=diagnostic,surviving_boundary=surviving,
        small_nonzero_perturbation=dict(z=[S-perturb,perturb,0],lambda_value=lam_pert,
            Weinberg_norm_per_GeV=float(np.linalg.norm(C,2))),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(three_RH_neutrinos_positive_bounded_Majorana_masses=True,
            exact_optimal_tree_level_Weinberg_norm_for_fixed_Yukawa_singular_values=True,
            tree_level_matching_at_a_common_scale_supplied=True,
            finite_eta_quartic_floor_necessary_not_sharp=True,
            reference_mass_conversion_not_RG_or_empirical_fit=True,
            radiative_stability_or_PMNS_or_complete_observed_spectrum_verified=False,
            all_spectral_models_or_cognition_excluded=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else: assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','common_scale_diagnostic','surviving_boundary')},ensure_ascii=False))
