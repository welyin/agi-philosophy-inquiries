"""Round 543: allowed fluctuations versus a jointly normalized singlet candidate.

Explicit extra singlet assumption; same-scale bare a4 potential, not RG/pole fit.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import shared_trace_unimodularity as finite
import protected_pair_finite_matching as matching

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_singlet_origin_and_vacuum_results.json'


def algebra_basis(N):
    zeroq=np.zeros((2,2));zeroc=np.zeros((N,N))
    basis=[finite.rep(N,z,zeroq,zeroc) for z in (1,1j)]
    pauli=(np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))
    basis += [finite.rep(N,0,q,zeroc) for q in [np.eye(2)]+[1j*s for s in pauli]]
    for i in range(N):
        for j in range(N):
            e=np.zeros((N,N),complex);e[i,j]=1
            basis += [finite.rep(N,0,zeroq,z*e) for z in (1,1j)]
    return basis


def neutral_dirac(Y,K,h,s):
    D=np.zeros((12,12),complex)
    D[0:3,3:6]=h*Y.conj().T;D[3:6,0:3]=h*Y
    D[6:9,9:12]=h*Y.T;D[9:12,6:9]=h*Y.conj()
    D[0:3,6:9]=s*K.conj().T;D[6:9,0:3]=s*K
    return D


def couplings(q,S,r,T):
    return (3*q*q+r*S*S)/T,S,T/r


def vacuum(q,S,r,T,R=1.):
    assert q>0 and q>S
    lh,p,ls=couplings(q,S,r,T)
    h2=2*R/q;s2=2*r*R/T*(1-S/q)
    H=2*np.array([[lh*h2,p*math.sqrt(h2*s2)],
                  [p*math.sqrt(h2*s2),ls*s2]])
    return dict(h2=h2,s2=s2,lambda_H=lh,portal=p,lambda_s=ls,
        threshold_drop=p*p/ls,lambda_effective=lh-p*p/ls,
        squared_vev_ratio=s2/h2,Majorana_squared_over_h2=(q-S)/2,
        scalar_squared_mass_eigenvalues=np.linalg.eigvalsh(H).tolist())


def run():
    checks=[];origin_max=0.;oneform_max=0.
    for N in (3,5):
        Z,J,_,D=finite.structures(N,1.,2.3232141058)
        DM=np.zeros_like(D);DM[0,4]=DM[4,0]=11
        basis=algebra_basis(N)
        for B in basis:
            comm=DM@B-B@DM
            origin_max=max(origin_max,float(np.max(np.abs(comm))))
            assert np.max(np.abs(comm))==0
            omega=D@B-B@D
            for A in basis:
                W=A@omega;Wreal=W+J@W.conj()@J
                block=np.array([Wreal[0,4],Wreal[4,0]])
                oneform_max=max(oneform_max,float(np.max(np.abs(block))))
                assert np.max(np.abs(block))==0
                assert np.max(np.abs(Z@W-W@Z))<1e-12
    checks.append('Majorana_block_commutes_with_entire_existing_real_algebra_basis')
    checks.append('all_real_inner_one_forms_leave_Majorana_block_fixed_even_with_weights')

    rng=np.random.default_rng(543);moment_error=0.
    for _ in range(20):
        q=.2;S=.4;r=2.3;k=.7;h=.6;s=1.1
        c=rng.normal(size=3)+1j*rng.normal(size=3);c*=math.sqrt(S)/np.linalg.norm(c)
        Y=np.column_stack([c,1j*c,np.zeros(3)])/math.sqrt(2)
        K=np.diag([k,k,0.]);D=neutral_dirac(Y,K,h,s)
        a=3*q+r*S;b=3*q*q+r*S*S;ck=2*r*k*k;dk=2*r*k**4;ek=r*k*k*S
        got2=r*np.trace(D@D).real+12*q*h*h
        got4=r*np.trace(np.linalg.matrix_power(D,4)).real+12*q*q*h**4
        expected2=4*a*h*h+2*ck*s*s
        expected4=4*b*h**4+8*ek*h*h*s*s+2*dk*s**4
        moment_error=max(moment_error,abs(got2-expected2),abs(got4-expected4))
        assert np.isclose(got2,expected2,rtol=1e-14) and np.isclose(got4,expected4,rtol=1e-14)
    checks.append('complex_protected_pair_weighted_spectral_moments_from_full_neutral_matrix')

    r=F(7,3);q=F(3,10);S=F(1,5);T=3*q+r*S;k=F(2,5)
    b=3*q*q+r*S*S;ck=2*r*k*k;dk=2*r*k**4;ek=r*k*k*S
    B=1/T
    lh,p,ls=couplings(q,S,r,T)
    assert 2*ek/ck==p and 2*dk/(B*ck*ck)==ls
    assert lh-p*p/ls==3*q*q/T
    assert lh*ls-p*p==3*q*q/r>0
    checks.append('canonical_portal_self_coupling_and_exact_Schur_complement')

    for q,S,r in ((F(3,10),F(1,5),F(7,3)),(F(1,2),F(1,10),F(2)),
                  (F(1,3),F(1,4),F(5,2))):
        T=3*q+r*S;lh,p,ls=couplings(q,S,r,T)
        h2=2/q;s2=2*r/T*(1-S/q)
        assert lh*h2+p*s2==2 and p*h2+ls*s2==2
        assert s2/h2==r*(q-S)/T<=r/3
        assert (T/(2*r))*s2/h2==(q-S)/2<=T/6
        eig=vacuum(float(q),float(S),float(r),float(T))['scalar_squared_mass_eigenvalues']
        assert min(eig)>0 and np.allclose(eig,sorted([4.,float(12*(q-S)/T)]),rtol=1e-13)
        assert 2*q<=2*T/3
    checks.append('common_mass_coexistence_vacuum_and_nonhierarchical_scale_bounds')

    # Equal point has zero singlet curvature; reversed inequality has no positive coexistence s².
    for q,S in ((F(1,3),F(1,3)),(F(1,4),F(1,3))):
        r=F(7,3);T=3*q+r*S;lh,p,ls=couplings(q,S,r,T)
        formal_s2=2*r/T*(1-S/q)
        assert formal_s2<=0
        curvature=-2+p*(2/lh)
        assert (curvature==0)==(q==S)
        assert curvature>=0
    # q=0 is a flat angular family, not an isolated coexistence minimum.
    q=F(0);S=F(1,3);r=F(7,3);T=r*S;lh,p,ls=couplings(q,S,r,T)
    assert lh==p==ls and lh*ls-p*p==0
    flat_ratios=[]
    for h2 in (F(1,10),F(1,10000),F(1,10**8)):
        s2=2/S-h2
        assert lh*h2+p*s2==2 and p*h2+ls*s2==2
        flat_ratios.append(float((T/(2*r))*s2/h2))
    assert flat_ratios[-1]>1e7
    checks.append('equal_boundary_zero_mass_reversed_boundary_and_flat_top_zero_exceptions')

    r=matching.BOUND['r'];T=matching.BOUND['T'];w=T/(3+r)
    rows=[]
    for ratio in (1.1,1.5,3.):
        S=T/(3*ratio+r);q=ratio*S;row=vacuum(q,S,r,T)
        row.update(top_to_neutrino_squared_ratio=ratio,q=q,S=S)
        assert row['lambda_effective']>3*w*w/T
        rows.append(row)
    equal=dict(q=w,S=w,original_lambda=w,
        formal_effective_lambda=3*w*w/T,formal_negative_shift=r*w*w/T,
        singlet_vev=0.,singlet_curvature=0.,heavy_threshold_valid=False)
    checks.append('historical_common_weight_joint_branch_and_balanced_invalid_threshold')

    # Splitting the two mass coefficients removes the hierarchy bound: a scope control.
    S=T/(4.5+r);q=1.5*S;lh,p,ls=couplings(q,S,r,T)
    h2=2.;s2=2e8;mh=lh*h2+p*s2;ms=p*h2+ls*s2
    assert abs(mh-ms)>1e6 and s2/h2>r/3
    assert lh*ls-p*p>0
    control=dict(independent_mass_coefficients=[mh,ms],squared_vev_ratio=s2/h2,
                 common_mass_condition_satisfied=False)
    checks.append('independent_mass_inputs_allow_hierarchy_and_delimit_the_bare_obstruction')
    deps=('research_note_532.md','shared_trace_unimodularity.py','research_note_533.md',
          'joint_yukawa_higgs_matching_results.json','research_note_542.md',
          'unified_physics_condition_ledger_542.md')
    return dict(round=543,tests_run=len(checks),failures=0,errors=0,checks=checks,
        origin_commutator_max=origin_max,Majorana_one_form_block_max=oneform_max,
        weighted_moment_max_error=moment_error,historical_weight_ratio=r,T=T,
        equal_boundary=equal,coexistence_examples=rows,
        general_bounds=dict(squared_vev_ratio_upper=r/3,Majorana_squared_over_h2_upper=T/6,
                            strict_coexistence_effective_lambda_lower=3*w*w/T),
        split_mass_counterexample=control,
        zero_top_flat_family_Majorana_squared_ratios=flat_ratios,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(singlet_is_explicit_extra_candidate=True,generated_by_existing_inner_fluctuations=False,
            same_scale_bare_a4_common_mass_potential_only=True,
            scale_bounds_require_nonzero_top_and_strict_radial_coexistence=True,
            scalar_threshold_requires_actual_heavy_mode=True,
            quantum_running_hierarchy_excluded=False,full_RG_or_observable_fit=False,
            full_unification_completed=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','equal_boundary','general_bounds')},ensure_ascii=False))
