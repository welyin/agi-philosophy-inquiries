"""Finite output recovery test with a faithful mixed gauge factor.

Computes Kraus products of the full enlarged output channel, not just its
logical marginal. A clean flag is a counterexample to dropping faithfulness.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'mixed_output_recovery_boundary_results.json'
sys.path.insert(0,str(HERE.parent/'818'))
import preparation_noise_bridge as old

def run():
    a,_,_,_=old.prep.previous.system(4);F=old.record_frame(a)
    g=[]
    for c in a:g.extend([c+c.conj().T,1j*(c-c.conj().T)])
    H=sum((.5j*t*g[i]@g[i+4] for i,t in enumerate((1.,1.5,.75,1.25))),np.zeros((16,16),complex))
    HI=.125j*g[0]@g[5]-.0625j*g[2]@g[7]
    H=F.conj().T@H@F;HI=F.conj().T@HI@F
    X=old.X;Y=old.Y;eye=np.eye(16);tau=np.eye(8)/8
    V=np.kron(X,H)+np.kron(Y,HI)
    vals,vec=np.linalg.eigh(V)
    embeddings=[np.kron(np.eye(2),np.eye(8)[:,i:i+1])/np.sqrt(8) for i in range(8)]
    def defect(kraus):
        return max(float(np.linalg.norm(k.conj().T@l-np.trace(k.conj().T@l)*np.eye(2)/2)) for k in kraus for l in kraus)
    zero_defect=defect(embeddings);assert zero_defect<1e-14
    rows=[]
    for eps in (.02,.01,.005):
        U=(vec*np.exp(-1j*eps*vals))@vec.conj().T
        blocks=U.reshape(2,16,2,16)
        kraus=[blocks[b,:,0,:]@e for b in range(2) for e in embeddings]
        complete=sum((k.conj().T@k for k in kraus),np.zeros((2,2),complex))
        err=float(np.linalg.norm(complete-np.eye(2)))
        d=defect(kraus);assert err<1e-12 and d>eps*.01
        rows.append(dict(epsilon=eps,full_output_KL_traceless_defect=d,defect_divided_by_epsilon=d/eps,completeness_error=err))
    # Compute the second-order decoded Choi kernel using the actual joint V.
    def decoded_second(rho):
        initial=np.kron(np.diag([1.,0.]),np.kron(rho,tau))
        d=V@initial@V-.5*(V@V@initial+initial@V@V)
        full=np.einsum('aras->rs',d.reshape(2,16,2,16))
        return np.einsum('arbr->ab',full.reshape(2,8,2,8))
    choi=sum((np.kron(np.outer(np.eye(2)[:,i],np.eye(2)[:,j]),decoded_second(np.outer(np.eye(2)[:,i],np.eye(2)[:,j])))/2 for i in range(2) for j in range(2)),np.zeros((4,4),complex))
    omega=np.array([1.,0.,0.,1.])/np.sqrt(2);p=np.eye(4)-np.outer(omega,omega)
    kernel=p@choi@p;ev=np.linalg.eigvalsh(kernel)
    assert ev.min()>-1e-12 and ev.max()>.1
    # A pure flagged input admits a zeroth-order decoder unequal to trace.
    prob=.2;e0=np.array([[1.],[0.]]);e1=np.array([[0.],[1.]])
    flagged=[np.sqrt(1-prob)*np.kron(np.eye(2),e0),np.sqrt(prob)*np.kron(X,e1)]
    recovery=[np.kron(np.eye(2),e0.T),np.kron(X,e1.T)]
    err=0.
    for i in range(2):
        for j in range(2):
            E=np.outer(np.eye(2)[:,i],np.eye(2)[:,j]);out=sum((k@E@k.conj().T for k in flagged),np.zeros((4,4),complex))
            restored=sum((r@out@r.conj().T for r in recovery),np.zeros((2,2),complex))
            err=max(err,float(np.linalg.norm(restored-E)))
    assert err<1e-12 and defect(flagged)<1e-12
    return dict(round=828,all_checks_passed=True,logical_dimension=2,enlarged_fermion_output_dimension=16,
        gauge_state_min_eigenvalue=.125,zero_order_KL_defect=zero_defect,full_output_rows=rows,
        decoded_second_order_Choi_kernel_eigenvalues=ev.tolist(),
        pure_flag_counterexample=dict(initial_flag_rank=1,flag_dimension=2,complete_channel_recovery_error=err,
            full_output_KL_defect=defect(flagged)),
        original_continuum_numerically_computed=False,
        scope='Faithful mixed output contract and vanishing relative first-order channel; not arbitrary enlarged outputs or encodings.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
