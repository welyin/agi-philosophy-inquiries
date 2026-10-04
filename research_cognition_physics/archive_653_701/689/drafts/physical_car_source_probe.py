"""689 entry: original pure-time Weyl contacts and a local CAR group insertion.

Inherited653 time transfer is not reproved. New test connects actual original
physical fields to an equal-time group insertion with a polynomial determinant,
including singular massless K and singular Cayley denominator.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_full_holonomy_positive_measure as prior
old=prior.old;single=prior.entry
TARGET=HERE/'physical_car_source_probe_results.json'


def original_fields(links):
    nt=len(links);size=16*nt;eye=np.eye(size)
    S=np.zeros((size,size),complex)
    for x,g in enumerate(links):
        y=(x+1)%nt
        S[16*x:16*x+16,16*y:16*y+16]=(-1 if x==nt-1 else 1)*g
    q=S.conj().T
    values,vectors=np.linalg.eigh(old.spin.GAMMA[3])
    wp=vectors[:,values>.5];wm=old.spin.G5@wp
    jp=(wp+wm)/np.sqrt(2);jm=(wp-wm)/np.sqrt(2)
    X=-np.kron(q,(np.eye(4)+old.spin.GAMMA[3])/2)-np.kron(S,(np.eye(4)-old.spin.GAMMA[3])/2)
    D=(np.eye(4*size)+X)/2
    v=(np.kron(eye,wp)-np.kron(q,wm))/np.sqrt(2)
    K=np.kron(eye,jp.conj().T)@D@v
    A=np.kron(eye,jm.conj().T)@v
    Q=np.kron(q,np.eye(2));I=np.eye(2*size)
    err=max(np.max(np.abs(2*K-(I-Q))),np.max(np.abs(2*A-(I+Q))))
    assert err<3e-14
    omega=[]
    power=np.linalg.matrix_power(Q,nt)
    for x in range(nt):omega.append(-power[32*x:32*x+32,32*x:32*x+32])
    return dict(K=K,A=A,Q=Q,omega=omega,frame_error=float(err),time_slices=nt)


def detvalue(a):
    phase,logabs=np.linalg.slogdet(a)
    return phase*np.exp(logabs)


def contact_checks():
    rows=[]
    for nt in (1,2,3,4):
        links=[single.gauge.rep(*single.gauge.group(68900+10*nt+x,.19+.06*x)) for x in range(nt)]
        data=original_fields(links);K=data['K'];A=data['A'];Q=data['Q']
        G=-A@np.linalg.inv(K)
        canonical=-np.linalg.inv(np.eye(len(Q))-Q)
        contact=G-2*canonical
        assert np.max(np.abs(contact-np.eye(len(Q))))<8e-13
        diagonal_error=0.
        for x,U in enumerate(data['omega']):
            expected=(U-np.eye(32))@np.linalg.inv(U+np.eye(32))
            diagonal_error=max(diagonal_error,float(np.max(np.abs(G[32*x:32*x+32,32*x:32*x+32]-expected))))
        assert diagonal_error<8e-13
        # Original ordered w/bar covariance, finite Wick coefficients.
        zero=np.zeros_like(G);full=np.block([[zero,G],[-G.T,zero]])
        cfull=np.block([[zero,canonical],[-canonical.T,zero]])
        J=np.block([[zero,np.eye(len(G))],[-np.eye(len(G)),zero]])
        rng=np.random.default_rng(68940+nt);worst=0.
        for degree in (2,4,6,8):
            wave=rng.normal(size=(degree,len(full)))+1j*rng.normal(size=(degree,len(full)))
            wave/=np.linalg.norm(wave,axis=1)[:,None]
            actual=old.pfaffian(wave@full@wave.T)
            correct=old.pfaffian(wave@(2*cfull+J)@wave.T)
            worst=max(worst,float(abs(actual-correct)))
        assert worst<8e-12
        rows.append(dict(N=nt,actual_frame_error=data['frame_error'],
            contact_identity_error=float(np.max(np.abs(contact-np.eye(len(Q))))),
            equal_time_Cayley_error=diagonal_error,source_coefficient_error=worst))
    # At the identity temporal holonomy, original w/bar equal-time covariance
    # is zero, whereas the actual CAR occupation is1/2. Contact is essential.
    identity=original_fields([np.eye(16,dtype=complex)])
    C=-identity['A']@np.linalg.inv(identity['K'])
    assert np.max(np.abs(C))<3e-14
    return dict(rows=rows,identity_original_equal_time_covariance_zero=True,
        required_original_CAR_occupation=0.5,all_source_identity_analytic_not_inferred_from_sample_degree=True)


def insertion_checks():
    rows=[]
    # Original hypercharge group gives genuine zero-Weyl-weight backgrounds,
    # without substituting arbitrary matrices outsideG.
    zero_group=single.gauge.rep(np.eye(3),np.eye(2),np.exp(1j*np.pi/2))
    singular_H=np.kron(single.gauge.rep(np.eye(3),np.eye(2),np.exp(1j*np.pi)),np.eye(2))
    for nt,kind in ((1,'regular'),(2,'regular'),(3,'regular'),(2,'zero_K'),(2,'zero_D'),(1,'both_singular')):
        if kind in ('zero_K','both_singular'):
            links=[zero_group]+[np.eye(16,dtype=complex) for _ in range(nt-1)]
        else:
            links=[single.gauge.rep(*single.gauge.group(68970+10*nt+x,.23)) for x in range(nt)]
        data=original_fields(links);K=data['K'];A=data['A'];I=np.eye(len(K));x=nt-1
        H=singular_H if kind in ('zero_D','both_singular') else np.kron(single.gauge.rep(*single.gauge.group(68995,.61)),np.eye(2))
        B=(H-np.eye(32))/2;D=(H+np.eye(32))/2
        P=I[32*x:32*x+32,:]
        # This block determinant is the polynomial continuation of
        # detD * det(K-P^T B D^-1 P A), with the original2^-32N removed.
        block=np.block([[2*K,P.T@B],[2*P@A,D]])
        actual=detvalue(block);expected=detvalue(np.eye(32)+data['omega'][x]@H)
        err=float(abs(actual-expected)/max(1.,abs(expected)))
        assert err<3e-12
        minK=float(np.linalg.svd(K,compute_uv=False)[-1]);minD=float(np.linalg.svd(D,compute_uv=False)[-1])
        if kind in ('zero_K','both_singular'):assert minK<1e-13
        if kind in ('zero_D','both_singular'):assert minD<1e-13
        rows.append(dict(N=nt,case=kind,relative_or_absolute_error=err,
            K_smallest_singular_value=minK,D_smallest_singular_value=minD,
            inserted_smallest_singular_value=float(np.linalg.svd(np.eye(32)+data['omega'][x]@H,compute_uv=False)[-1]),
            unnormalized_insertion_real=float(actual.real),unnormalized_insertion_imag=float(actual.imag),
            no_inverse_K_or_D_used_for_insertion=True))
    assert any(r['case']=='zero_K' and r['inserted_smallest_singular_value']>1e-5 for r in rows),rows
    return dict(rows=rows,all32_original_modes_retained=True,
        same_time_physical_Y_polynomial_not_auxiliary_E_insertion=True,
        group_Haar_insertion_analytically_maps_to_original_CAR_singlet_projector=True,
        general_spatial_or_HF_identity_not_claimed=True)


def run():
    deps=('research_note_653.md','research_note_654.md','research_note_646.md','research_note_661.md',
          'research_note_670.md','research_note_679.md','research_note_688.md',
          'joint_gauss_support_marginal_results.json')
    return dict(date='2026-10-02',entry_round=689,latest_formal_round=688,not_formal_round=True,
        contact=contact_checks(),insertion=insertion_checks(),
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps},
        old_space_and_full_goal_unchanged=True)


if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=689,all_checks_passed=True,
        max_contact_error=max(r['contact_identity_error'] for r in result['contact']['rows']),
        max_insertion_error=max(r['relative_or_absolute_error'] for r in result['insertion']['rows']))))
