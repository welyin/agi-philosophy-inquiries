"""818: full-state parity twirl, explicit logical inputs and second-order noise.
Independent finite identities and original neutral-coefficient diagnostics.
No original continuum covariance, finite-coupling limit or device is simulated.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'817'))
import even_preparation_witness as prep
TARGET=HERE/'preparation_noise_bridge_results.json'
I2=np.eye(2);X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.])
PAULI=[I2,X,Y,Z]

def record_frame(a):
    I=np.eye(a[0].shape[0]);p=I-a[0].conj().T@a[0]
    e10=a[0].conj().T@(a[1]+a[1].conj().T)
    B=I[:,np.where(np.diag(p)>.5)[0]]
    F=np.column_stack((B,e10@B))
    assert np.linalg.norm(F.conj().T@F-I)<1e-12
    return F

def complete_state_twirl():
    d=prep.data();units=d[6]
    a=prep.previous.old.car(6);F=record_frame(a);dim=64
    rng=np.random.default_rng(818);v=rng.normal(size=dim)+1j*rng.normal(size=dim)
    for i in range(dim):
        if i.bit_count()%2:v[i]=0
    v/=np.linalg.norm(v)
    state=np.outer(v,v.conj());twirled=np.zeros((dim,dim),complex)
    for U in units:
        w=np.kron(I2,U)@v;twirled+=np.outer(w,w.conj())/256
    blocks=F.conj().T@twirled@F;shaped=blocks.reshape(2,32,2,32)
    corner=np.einsum('aras->rs',shaped)
    product=np.kron(I2/2,corner)
    factor_error=float(np.linalg.norm(blocks-product))
    assert factor_error<1e-12
    before=F.conj().T@state@F
    before_corner=np.einsum('aras->rs',before.reshape(2,32,2,32))
    before_record=np.einsum('arbr->ab',before.reshape(2,32,2,32))
    reference_defect=float(np.linalg.norm(before-np.kron(before_record,before_corner)))
    assert reference_defect>.1
    rho=(I2+.2*X-.4*Y+.1*Z)/2
    eig,vec=np.linalg.eigh(rho)
    prepared=np.zeros_like(twirled);completeness=np.zeros_like(twirled)
    for weight,vector in zip(eig,vec.T):
        for k in range(2):
            L=np.sqrt(weight)*np.outer(vector,np.eye(2)[:,k])
            K=F@np.kron(L,np.eye(32))@F.conj().T
            prepared+=K@twirled@K.conj().T;completeness+=K.conj().T@K
    reset_error=float(np.linalg.norm(F.conj().T@prepared@F-np.kron(rho,corner)))
    assert reset_error<1e-12 and np.linalg.norm(completeness-np.eye(dim))<1e-12
    parity=np.diag([(-1)**i.bit_count() for i in range(dim)])
    assert abs(np.trace(prepared@parity)-1)<1e-12
    return dict(original_correlated_state_product_defect=reference_defect,
        full_state_product_error_after_twirl=factor_error,
        explicit_reset_error=reset_error,
        all_preparations_preserve_total_parity=True,
        outside_mode_was_included_in_factorization_test=True,
        autonomous_reset_proven=False)

def original_noise_dictionary():
    a,xi,Q,cov=prep.previous.system(4);F=record_frame(a)
    ids=[24,25,30,31,56,57,62,63]
    phi=np.array([0.,.6987151138647895,0.,0.,.545863690226873])
    matrices=[]
    for k in (1,4):
        M=prep.previous.old.vertex.dmass(phi,np.eye(5)[k])
        M=F.conj().T@Q(M[np.ix_(ids,ids)])@F
        matrices.append([np.einsum('ba,arbs->rs',s,M.reshape(2,8,2,8))/2 for s in PAULI])
    B=[sum((np.kron(s,matrices[k][r]) for k,s in enumerate((X,Y))),np.zeros((16,16),complex)) for r in range(4)]
    # The finite T block after the parity-preserving ensemble has trace state.
    env=np.kron(np.diag([1.,0.]),np.eye(8)/8)
    C=np.array([[np.trace(env@a@b) for b in B[1:]] for a in B[1:]])
    V=sum((np.kron(s,b) for s,b in zip(PAULI,B)),np.zeros((32,32),complex))
    joint=np.kron(I2/2,env)
    derivatives=[1j*(np.kron(s,np.eye(16))@V-V@np.kron(s,np.eye(16))) for s in PAULI[1:]]
    actual=np.array([[np.trace(joint@a@b).real for b in derivatives] for a in derivatives])
    R=C.real;predicted=4*(np.trace(R)*np.eye(3)-R)
    error=float(np.linalg.norm(actual-predicted))
    assert error<1e-12 and np.linalg.eigvalsh(actual).min()>1e-3
    def partial(M):return np.einsum('arbr->ab',M.reshape(2,16,2,16))
    def L(rho):
        M=np.kron(rho,env)
        return partial(V@M@V-.5*(V@V@M+M@V@V))
    J=np.zeros((4,4),complex)
    for i in range(2):
        for j in range(2):
            E=np.outer(np.eye(2)[:,i],np.eye(2)[:,j])
            J+=np.kron(E,L(E))/2
    omega=np.array([1.,0.,0.,1.])/np.sqrt(2)
    orth=np.eye(4)-np.outer(omega,omega)
    restriction=orth@J@orth
    noise_eig=np.linalg.eigvalsh(restriction)
    inverse_eig=np.linalg.eigvalsh(-restriction)
    expected=np.r_[0.,np.linalg.eigvalsh(C)]
    assert np.max(abs(noise_eig-expected))<1e-12
    assert inverse_eig.min() < -1e-3
    # This checks the actual finite unitary, not just the double-commutator.
    es,vs=np.linalg.eigh(V);rows=[]
    for epsilon in (.04,.02,.01):
        U=(vs*np.exp(-1j*epsilon*es))@vs.conj().T
        choi=np.zeros((4,4),complex)
        for i in range(2):
            for j in range(2):
                E=np.outer(np.eye(2)[:,i],np.eye(2)[:,j])
                N=partial(U@np.kron(E,env)@U.conj().T)
                choi+=np.kron(E,N)/2
        Fe=float(np.vdot(omega,choi@omega).real)
        rows.append(dict(epsilon=epsilon,entanglement_infidelity_coefficient=(1-Fe)/epsilon**2))
    assert abs(rows[-1]['entanglement_infidelity_coefficient']-np.trace(C).real)<1e-6
    return dict(original_neutral_coefficients_used=True,
        mean_square_Gram_identity_error=error,
        logical_change_Gram_eigenvalues=np.linalg.eigvalsh(actual).tolist(),
        noise_Gram_eigenvalues=np.linalg.eigvalsh(C).tolist(),
        Choi_kernel_positive_eigenvalues=noise_eig.tolist(),
        proposed_inverse_Choi_kernel_minimum=float(inverse_eig.min()),
        predicted_entanglement_infidelity_coefficient=float(np.trace(C).real),
        finite_series_checks=rows,
        original_continuous_Gram_computed=False,
        finite_coupling_original_channel_or_recovery_proven=False)

def run():
    return dict(round=818,all_checks_passed=True,
        full_state_preparation=complete_state_twirl(),
        noise=original_noise_dictionary(),
        scope='Full-state logical preparation identity, ordered noise and formal inverse positivity calibrated independently.',
        original_autonomous_resources_closed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))

