"""821 working calibration: one common controlled shift and its logical cost.
The oscillator is a diagnostic of the CCR identity, not a new physical field
or a proven implementation of the original homogeneous compensation.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
TARGET=HERE/'controlled_shift_probe_results.json'
I=np.eye(2);X=np.array([[0.,1],[1,0]],complex)
Y=np.array([[0.,-1j],[1j,0]]);Z=np.diag([1.,-1.])

def unitary(V,epsilon):
    eigen,frame=np.linalg.eigh(V)
    return (frame*np.exp(-1j*epsilon*eigen))@frame.conj().T
def reduce_record(state,N):return np.einsum('arbr->ab',state.reshape(2,N,2,N))

def run():
    N=24
    a=np.zeros((N,N),complex)
    for n in range(1,N):a[n-1,n]=np.sqrt(n)
    q=(a+a.conj().T)/np.sqrt(2);p=(a-a.conj().T)/(1j*np.sqrt(2))
    vac=np.zeros((N,N));vac[0,0]=1
    inputs=[(I+Z)/2,(I-Z)/2,(I+X)/2,(I+.2*X-.4*Y+.1*Z)/2]
    rows=[]
    for eps in (.2,.1,.05):
        U=unitary(np.kron(Z,p),eps);outputs=[]
        mean_error=0.;channel_error=0.;variance_error=0.
        for rho in inputs:
            state=U@np.kron(rho,vac)@U.conj().T;record=reduce_record(state,N)
            predicted=rho.copy();predicted[0,1]*=np.exp(-eps**2);predicted[1,0]*=np.exp(-eps**2)
            shift=np.trace(state@np.kron(I,eps*q)).real
            variance=np.trace(state@np.kron(I,q@q)).real
            mean_error=max(mean_error,abs(shift-eps**2*np.trace(rho@Z).real))
            variance_error=max(variance_error,abs(variance-(.5+eps**2)))
            channel_error=max(channel_error,float(np.linalg.norm(record-predicted)))
            outputs.append(state)
        assert max(mean_error,variance_error,channel_error)<2e-12
        mixed=.37*inputs[0]+.63*inputs[2]
        direct=U@np.kron(mixed,vac)@U.conj().T
        affine_error=float(np.linalg.norm(direct-(.37*outputs[0]+.63*outputs[2])))
        assert affine_error<2e-13
        rows.append(dict(epsilon=eps,physical_mean_shift_error=mean_error,
            logical_dephasing_formula_error=channel_error,
            quantum_second_moment_error=variance_error,
            full_state_affinity_error=affine_error,
            plus_X_coherence_loss=1-np.exp(-eps**2),
            coherence_loss_over_epsilon_squared=(1-np.exp(-eps**2))/eps**2))
    return dict(round=821,formal_round_completed=False,rows=rows,
        single_input_independent_joint_unitary=True,
        physical_mean_shift_order='epsilon^2 = hbar',
        normalized_field_covariance_unchanged_at_leading_order=True,
        logical_content_exactly_undisturbed=False,
        finite_Fock_CCR_claimed=False,
        original_compensation_in_physical_causal_range_proven=False,
        autonomous_compensation_preparation_proven=False,
        formal_test_groups_added=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
