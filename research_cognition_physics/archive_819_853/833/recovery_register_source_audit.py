"""833: the same decoded channel need not have the same register sources.

The finite error/syndrome representation is an exact algebraic counterexample
to inferring source suppression from logical fidelity. It is not a native
decoder constructed from the original graph interaction.
"""
from pathlib import Path
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'recovery_register_source_audit_results.json'


def model():
    i=np.eye(2);x=np.array([[0.,1.],[1.,0.]]);z=np.diag([1.,-1.])
    x1=np.kron(x,i);x2=np.kron(i,x);a=.7
    h=np.kron(i,np.kron(x1,z))+np.kron(i,np.kron(x2,x))+a*np.kron(z,np.kron(x1@x2,i))
    w=np.zeros((2,4,2,2))
    for logical in range(2):w[logical,0,0,logical]=1
    w=w.reshape(16,2)
    decoder=np.eye(16)
    for logical in range(2):
        for env in range(2):decoder[logical*8+3*2+env,logical*8+3*2+env]=(-1)**logical
    pointer=np.kron(i,np.kron(x2,i))
    return h,w,decoder,pointer,z,a


def matrix_derivative(h,w,d,o,n):
    from math import comb
    return sum(comb(n,k)*(1j)**k*(-1j)**(n-k)*
               (w.conj().T@np.linalg.matrix_power(h,k)@d.conj().T@o@d@np.linalg.matrix_power(h,n-k)@w)
               for k in range(n+1))


def decoded_choi(vectors):
    # Stinespring isometry: L, syndrome, original environment, input.
    psi=vectors.reshape(2,4,2,2).transpose(3,0,1,2).reshape(4,8)/np.sqrt(2)
    return psi@psi.conj().T


def run():
    h,w,d,o,z,a=model()
    derivatives=[matrix_derivative(h,w,d,o,n) for n in range(3)]
    assert np.linalg.norm(derivatives[0])<1e-14
    assert np.linalg.norm(derivatives[1])<1e-14
    logical_second=(np.trace(z@derivatives[2])/2).real
    assert abs(logical_second+4*a)<1e-13
    # This quadratic syndrome source has an order-t^2 logical coefficient,
    # although the decoded logical entangled-state infidelity is order t^4.
    expected_input_difference_coefficient=-4*a
    eig,vec=np.linalg.eigh(h)
    bell=np.eye(2).reshape(-1,order='F')/np.sqrt(2)
    rows=[]
    for t in (.01,.02,.04,.08):
        raw=vec@(np.exp(-1j*t*eig)[:,None]*(vec.conj().T@w))
        decoded=d@raw
        observable=decoded.conj().T@o@decoded
        difference=float((observable[0,0]-observable[1,1]).real)
        choi=decoded_choi(decoded)
        infidelity=float(1-(bell.conj()@choi@bell).real)
        # Copy the syndrome to an additional internal register. Tracing both
        # gives the identical logical channel, but a single register loses
        # its off-diagonal source. This does not construct the copying unitary.
        copied=np.zeros((2,4,4,2,2),complex)
        source=decoded.reshape(2,4,2,2)
        for s in range(4):copied[:,s,s,:,:]=source[:,s,:,:]
        cp_psi=copied.transpose(4,0,1,2,3).reshape(4,32)/np.sqrt(2)
        cp_choi=cp_psi@cp_psi.conj().T
        same_channel=float(np.linalg.norm(cp_choi-choi))
        assert same_channel<1e-14
        x=np.array([[0.,1.],[1.,0.]])
        x2=np.kron(np.eye(2),x)
        copied_vectors=copied.reshape(64,2)
        local_operator=np.kron(np.eye(2),np.kron(x2,np.kron(np.eye(4),np.eye(2))))
        joint_operator=np.kron(np.eye(2),np.kron(x2,np.kron(x2,np.eye(2))))
        local_matrix=copied_vectors.conj().T@local_operator@copied_vectors
        joint_matrix=copied_vectors.conj().T@joint_operator@copied_vectors
        copied_source=float(np.linalg.norm(local_matrix))
        assert copied_source==0 and np.linalg.norm(joint_matrix-observable)<1e-14
        assert abs(difference/t**2-expected_input_difference_coefficient)<.04
        assert infidelity<1.1*t**4
        rows.append(dict(time=t,register_source_input_difference=difference,
                         source_difference_over_t2=difference/t**2,
                         decoded_entangled_state_infidelity=infidelity,
                         same_decoded_channel_after_copy_residual=same_channel,
                         one_register_off_diagonal_source_after_copy=copied_source))
    return dict(round=833,all_checks_passed=True,fresh_test_groups=1,
                pointer_second_derivative_logical_Z_coefficient=float(logical_second),
                register_source_input_difference_t2_coefficient=expected_input_difference_coefficient,
                rows=rows,
                pointer_X_can_be_represented_as_dual_rail_quadratic_hopping=True,
                copying_requires_additional_internal_preparation_and_coupling=True,
                copying_full_joint_source_audit_completed=False,
                native_decoder_or_original_Yukawa_coefficient_not_computed=True,
                scope='Two finite Stinespring record realizations have the same decoded channel but different low-degree pointer observables. The first can show logical dependence at order t^2 while logical infidelity is order t^4. Copying hides that single-register coherence but leaves the full joint source and resource audit open; it is not a native autonomous repair.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
