"""684 entry: action covariance, anti-linear square, and a failed doubled repair.
Only the declared original compensating Gaussian sector, not physical Gauss RP.
"""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_local_source_lift as body
import joint_physical_source_reflection as reflection
base=body.base;err=body.err
TARGET=HERE/'compensation_involution_probe_results.json'

def run():
    links,e,phis=body.prior.fixture();mat=base.fixed_matrices(e,phis)
    n=len(mat[2]);g5=mat[1]@mat[1].T-mat[0]@mat[0].T
    perm=np.eye(4)[[2,3,0,1]];ell,_,_=reflection.reflection_matrices(mat,perm)
    t=g5@ell;R=np.kron(np.array([[0.,1.],[1.,0.]]),t)
    assert err(R.imag)==0 and err(R@R.conj()+np.eye(2*n))==0
    covariance=[];freeG=None
    for label,ls in [('nonflat',links),('free',np.tile(np.eye(16,dtype=complex),(2,4,1,1)))]:
        ls_t=ls[:,[2,3,0,1]].copy();ls_t[1]=np.array([u.conj().T for u in ls[1]])
        _,_,_,h,_=base.kernel(ls);_,_,_,ht,_=base.kernel(ls_t)
        k,_,_=body.blocks(g5@h,g5,.2,1);kt,_,_=body.blocks(g5@ht,g5,.2,1)
        G=k.conj().T@k;Gt=kt.conj().T@kt
        defect=err(R.conj().T@Gt@R-G)
        assert defect<1e-12
        covariance.append(dict(background=label,action_covariance_error=defect))
        if label=='free':freeG=G
    phase_rows=[]
    for phase in (0.,.37,np.pi/2):
        rr=np.exp(1j*phase)*R
        residual=err(rr@rr.conj()+np.eye(2*n));assert residual<1e-14
        phase_rows.append(dict(phase=phase,square_minus_identity_error=residual))
    # Two identical independent copies, J real and skew: S is a real involution.
    J=np.array([[0.,1.],[-1.,0.]])
    S=np.kron(J,R);G2=np.kron(np.eye(2),freeG)
    assert err(S@S.conj()-np.eye(4*n))==0
    assert err(S.conj().T@G2@S-G2)<1e-12
    C=np.linalg.inv(freeG)
    pos=np.r_[np.arange(128,256),np.arange(128,256)+n]
    B=(R@C)[np.ix_(pos,pos)]
    # For real Gaussian covariance, Q_ij=E[theta(b_i)b_j]=(S C2)_ij.
    gram=np.block([[np.zeros_like(B),B],[-B,np.zeros_like(B)]])
    assert err(gram-gram.conj().T)<1e-13
    ev,v=np.linalg.eigh((gram+gram.conj().T)/2)
    witness=v[:,0];norm=np.vdot(witness,gram@witness)
    singular=float(np.linalg.norm(B,2))
    assert abs(ev[0]+singular)<1e-13 and norm.real<-1e-4 and abs(norm.imag)<1e-14
    return dict(entry_round=684,latest_formal_round=683,not_formal_round=True,
        action_covariance=covariance,anti_linear_square_is_minus_identity=True,
        phase_cannot_fix_anti_linear_square=phase_rows,
        doubled_repair=dict(real_involution=True,action_invariant=True,
            exact_Gram_structure='[[0,B],[-B,0]], B*= -B',
            B_operator_norm=singular,minimum_Gram_eigenvalue=float(ev[0]),
            witness_value=base.old.cpair(norm),
            witness_vector=[[float(z.real),float(z.imag)] for z in witness]),
        failure_only_for_declared_auxiliary_completion=True,
        original_physical_Gauss_S9_RP_not_decided=True,
        duplicated_sector_not_claimed_as_preserving_physical_determinants_or_sources=True)

if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=684,all_checks_passed=True,
        action_covariance=result['action_covariance'],
        doubled_Gaussian_negative_norm=result['doubled_repair']['witness_value'])))
