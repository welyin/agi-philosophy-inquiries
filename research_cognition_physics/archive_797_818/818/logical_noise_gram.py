"""818 working: logical Pauli noise Gram in the original neutral CAR coefficients.
Finite boson and product-preparation diagnostics only. No original continuous
physical covariance or autonomous reset is supplied by this calculation.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'815'))
import record_retention_criterion as prev
TARGET=HERE/'logical_noise_gram_results.json'

def partial_env(rho):
    return np.einsum('arbr->ab',rho.reshape(2,16,2,16))

def run():
    a,xi,Q,cov=prev.system(4);eye=np.eye(16)
    e00=eye-a[0].conj().T@a[0];e10=a[0].conj().T@(a[1]+a[1].conj().T)
    columns=np.eye(16)[:,np.where(np.diag(e00)>.5)[0]]
    frame=np.column_stack((columns,e10@columns))
    assert np.linalg.norm(frame.conj().T@frame-np.eye(16))<1e-12
    ids=[24,25,30,31,56,57,62,63]
    take=lambda M:M[np.ix_(ids,ids)]
    phi=np.array([0.,.6987151138647895,0.,0.,.545863690226873])
    H=Q(take(prev.old.vertex.mass(phi)))
    _,vectors=np.linalg.eigh(H);v=frame.conj().T@vectors[:,0]
    reference=np.outer(v,v.conj()).reshape(2,8,2,8)
    sigma=np.einsum('aras->rs',reference)
    rhoA=np.einsum('arbr->ab',reference)
    product=np.kron(rhoA,sigma)
    correlation=float(np.linalg.norm(reference.reshape(16,16)-product))
    assert correlation>1e-3
    I=np.eye(2);X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.])
    pauli=[I,X,Y,Z];bosons=[X,Y]
    operators=[frame.conj().T@Q(take(prev.old.vertex.dmass(phi,np.eye(5)[i])))@frame for i in (1,4)]
    decomposed=[]
    for M in operators:
        blocks=M.reshape(2,8,2,8)
        parts=[np.einsum('ba,arbs->rs',p,blocks)/2 for p in pauli]
        assert np.linalg.norm(sum((np.kron(p,b) for p,b in zip(pauli,parts)),np.zeros((16,16),complex))-M)<1e-12
        decomposed.append(parts)
    # Environment order is finite boson, then the original logical corner.
    B=[sum((np.kron(b,decomposed[i][r]) for i,b in enumerate(bosons)),np.zeros((16,16),complex)) for r in range(4)]
    env=np.kron(np.diag([1.,0.]),sigma)
    C=np.array([[np.trace(env@r@s) for s in B] for r in B])
    assert np.linalg.norm(C-C.conj().T)<1e-12
    assert np.linalg.eigvalsh(C).min()>-1e-12
    assert max(abs(np.trace(env@b)) for b in B)<1e-12
    V=sum((np.kron(p,b) for p,b in zip(pauli,B)),np.zeros((32,32),complex))
    es,vs=np.linalg.eigh(V)
    def full_second(rho):
        R=np.kron(rho,env)
        return partial_env(V@R@V-.5*(V@V@R+R@V@V))
    def noise(rho,reverse=False):
        out=np.zeros((2,2),complex)
        for r in range(1,4):
            for s in range(1,4):
                coefficient=C[r,s] if reverse else C[s,r]
                out+=coefficient*(pauli[r]@rho@pauli[s]-.5*(pauli[s]@pauli[r]@rho+rho@pauli[s]@pauli[r]))
        return out
    h=sum((-C[0,r].imag*pauli[r] for r in range(1,4)),np.zeros((2,2),complex))
    probes=[(I+X)/2,(I+Y)/2,(I+Z)/2,I/2]
    error=wrong=0.
    for rho in probes:
        predicted=noise(rho)-1j*(h@rho-rho@h)
        error=max(error,float(np.linalg.norm(predicted-full_second(rho))))
        wrong=max(wrong,float(np.linalg.norm(noise(rho,True)-1j*(h@rho-rho@h)-full_second(rho))))
    assert error<1e-12
    # A numerical zero here would mean this particular calibration is real;
    # it would not license changing the general ordered formula.
    rows=[];rho=(I+Y)/2;coefficient=full_second(rho)
    for epsilon in (.04,.02,.01,.005):
        U=(vs*np.exp(-1j*epsilon*es))@vs.conj().T
        reduced=partial_env(U@np.kron(rho,env)@U.conj().T)
        defect=float(np.linalg.norm((reduced-rho)/epsilon**2-coefficient))
        rows.append(dict(epsilon=epsilon,second_order_residual=defect))
    assert rows[-1]['second_order_residual']<rows[0]['second_order_residual']/5
    eig=np.linalg.eigvalsh(C[1:,1:])
    assert eig.max()>1e-4
    return dict(round=818,status='working_not_formal',all_checks_passed=True,
        original_neutral_mass_variations_used=True,
        original_correlated_state_product_replacement_distance=correlation,
        logical_noise_Gram_real=C[1:,1:].real.tolist(),
        logical_noise_Gram_imag=C[1:,1:].imag.tolist(),
        logical_noise_Gram_eigenvalues=eig.tolist(),
        full_second_order_dictionary_error=error,
        reversed_order_coefficient_defect=wrong,
        short_time_series_checks=rows,
        pure_environment_cross_Hamiltonian_norm=float(np.linalg.norm(h)),
        original_physical_W_B_computed=False,
        original_continuous_logical_noise_nonzero_proven=False,
        original_autonomous_preparation_or_recovery_proven=False,
        formal_numbered_test_groups_added=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

