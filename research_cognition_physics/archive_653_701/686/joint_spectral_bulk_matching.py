"""686: declared spectral bulk matching and original all-source limit.

Numerics check exact matrices and fixed-source examples. Full H_b/double Haar/S9
convergence is analytic in the note, using the already proved 669/673/677 bounds.
Neither positive Euclidean weights nor convergence is a reflection positivity proof.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_body_weight_source_matching as old
import joint_local_source_lift as body

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_spectral_bulk_matching_results.json'
prior=body.prior
base=body.base


def factors(h,g5,a,layers):
    n=len(h)
    B=2*np.eye(n)+a*g5@h
    phase,lb=np.linalg.slogdet(B)
    assert abs(phase-1)<2e-12
    eps,ha,ev=prior.regulate(h,g5,a,layers)
    absolute=np.abs(ev)
    sbulk=float((layers+1)*lb+layers*np.sum(np.log1p(absolute)))
    logR=float(np.sum(np.log1p(np.exp(-2*layers*np.abs(np.arctanh(ev))))))
    assert 0<=logR<=n*np.log(2)+1e-12
    logdet=float((layers+1)*lb+np.sum(np.logaddexp(layers*np.log1p(ev),layers*np.log1p(-ev))))
    assert abs(logdet-sbulk-logR)<2e-8*max(1,logdet/1e6)
    return eps,sbulk,logR,logdet


def matching_identity_and_limits():
    links,e,phis=prior.fixture()
    _,_,_,H,_=base.kernel(links)
    mat=base.fixed_matrices(e,phis)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    direct=[]
    for layers in (1,3):
        _,S,logR,formula=factors(H,g5,.2,layers)
        K,_,_=body.blocks(g5@H,g5,.2,layers)
        phase,actual=np.linalg.slogdet(K)
        residual=max(abs(actual-S-logR),abs(actual-formula),abs(phase-1))
        assert residual<2e-9
        direct.append(dict(L=layers,original_dimension=len(H),log_residual=logR,
                           direct_determinant_error=float(residual)))
    # Controlled zero modes: not a negative or singular witness for the original gauge model.
    hc=np.diag([-1.,0.,0.,1.]);gc=np.diag([1.,1.,-1.,-1.])
    _,_,zero_log,_=factors(hc,gc,.001,100000)
    assert abs(zero_log-2*np.log(2))<1e-13
    # If aL stays bounded the residual need not approach one, even at a gapped point.
    g1=np.kron(np.eye(16),old.old.spin.G5)
    h1=old.old.frames(0)[3]
    bounded=[]
    target=64*np.log1p(np.exp(-1.))
    for a in (.01,.002,.0004):
        L=round(1/a)
        _,_,r,_=factors(h1,g1,a,L)
        bounded.append(dict(a=a,L=L,aL=a*L,log_residual=r,error=abs(r-target)))
    assert bounded[-1]['error']<.004 and bounded[-1]['error']<bounded[0]['error']
    # Only leading subtraction leaves a BACKGROUND-DEPENDENT finite action in an allowed limit.
    theta=.17
    q=old.old.Q;c=np.cos(q*theta)
    predicted=float(32-np.sum(c*(c+1)))
    truncation=[]
    for a in (.04,.02,.01,.005):
        L=round(1/a**2)
        values=[]
        for t in (theta,0.):
            costheta=np.cos(q*t);d=4+4*a*costheta+a*a
            # Stable S_bulk minus linear part, with the huge n L log 2 removed analytically.
            excess=float(2*np.sum(np.log(d/4))-2*a*np.sum(costheta)+
                4*L*np.sum(np.log((np.sqrt(d)+a)/2)-a*(costheta+1)/2))
            values.append(excess)
        contrast=values[0]-values[1]
        truncation.append(dict(a=a,L=L,La_squared=L*a*a,
                               relative_log_weight_after_linear_subtraction=contrast,
                               analytic_limit=predicted,error=abs(contrast-predicted)))
    assert predicted>0 and truncation[-1]['error']<.06
    assert truncation[-1]['error']<truncation[0]['error']
    return dict(actual_nonflat_K_checks=direct,zero_mode_control_log_residual=zero_log,
        zero_modes_handled_by_measure_zero_not_uniform_claim=True,
        bounded_aL_control=bounded,bounded_aL_limit=target,
        leading_only_matching_original_allowed_rate_counterexample=truncation,
        leading_subtraction_extra_rate_not_required_for_exact_matching=True)


def physical_source_limit():
    links,e,phis=prior.fixture()
    u,v,d,H,gap=base.kernel(links)
    mat=base.fixed_matrices(e,phis)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    target_eps=v@v.conj().T-u@u.conj().T
    exact=prior.soft(target_eps,mat,.37)
    rng=np.random.default_rng(68621)
    z=rng.normal(size=(len(H),4))+1j*rng.normal(size=(len(H),4))
    z/=np.linalg.norm(z,axis=0)
    target=[prior.source_coeff(exact,z[:,:k]) for k in (0,2,4)]
    assert min(abs(x) for x in target)>1e-30
    rows=[]
    for a,L in ((.125,128),(.03125,1024),(.0078125,8192),(.001953125,65536),(.00048828125,524288)):
        eps,S,logR,logdet=factors(H,g5,a,L)
        q=prior.soft(eps,mat,.37)
        R=float(np.exp(logR))
        coefficients=[prior.source_coeff(q,z[:,:k]) for k in (0,2,4)]
        matched=[R*x for x in coefficients]
        errors=[float(abs(x-y)/abs(y)) for x,y in zip(matched,target)]
        # This is an identity check for scalar matching, not an auxiliary positivity test.
        bound_checks=[]
        for raw,actual,want in zip(coefficients,matched,target):
            bound=R*abs(raw-want)+abs(R-1)*abs(want)
            assert abs(actual-want)<=bound+2e-13*abs(want)
            bound_checks.append(float(bound/abs(want)))
        rows.append(dict(a=a,L=L,aL=a*L,original_Wilson_gap=gap,
            log_residual=logR,residual_weight=R,
            matched_coefficients=[base.old.cpair(c) for c in matched],
            relative_target_errors=errors,triangle_relative_bounds=bound_checks))
    assert all(rows[-1]['relative_target_errors'][i]<rows[0]['relative_target_errors'][i] for i in range(3))
    assert max(rows[-1]['relative_target_errors'])<.05, [(r['a'], r['relative_target_errors']) for r in rows]
    return dict(original_nonflat_16_channel_massive_fixture=True,lambda_value=.37,
        source_orders=[0,2,4],exact_target_coefficients=[base.old.cpair(c) for c in target],
        rows=rows,full_Haar_S9_thermal_integral_not_numerically_replaced=True,
        weighted_L1_all_fixed_order_sources_proof_is_analytic=True)


def run():
    deps=('research_note_669.md','research_note_673.md','research_note_677.md',
          'research_note_678.md','research_note_685.md','joint_rational_physical_limit.py',
          'joint_local_source_lift.py','joint_body_weight_source_matching.py',
          'joint_body_weight_source_matching_results.json')
    return dict(date='2026-10-02',round=686,tests_run=2,failures=0,errors=0,
        exact_matching_and_failure_controls=matching_identity_and_limits(),
        same_physical_sources=physical_source_limit(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(declared_nonlocal_spectral_matching=True,
            all_fixed_order_physical_sources_complete_average_limit=True,
            uniform_spectral_gap_required=False,additional_La_squared_rate_required=False,
            uniform_volume_limit_proved=False,normalized_functional_without_nonzero_Z_proved=False,
            generic_gauge_or_geometry_derivatives_proved=False,
            positive_auxiliary_process_or_physical_RP_proved=False,
            original_HF_identity_proved=False,quantum_GR_completed=False,
            old_space_and_full_goal_unchanged=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=686,tests_run=2,all_checks_passed=True,
        final_errors=result['same_physical_sources']['rows'][-1]['relative_target_errors'],
        linear_matching_contrast=result['exact_matching_and_failure_controls']['leading_only_matching_original_allowed_rate_counterexample'][-1])))
