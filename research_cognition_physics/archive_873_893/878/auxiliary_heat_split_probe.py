"""Working878: auxiliary trace-class damping alone is not a regional split criterion."""
from pathlib import Path
import argparse,itertools,json,math
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'auxiliary_heat_split_probe_results.json'
def run():
    paulis=[np.eye(2,dtype=complex),np.array([[0,1],[1,0]],complex),
        np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]
    rows=[]
    for n in (1,2,3):
        words=list(itertools.product(range(4),repeat=n));dim=2**n
        mats=[]
        for word in words:
            a=np.ones((1,1),complex)
            for k in word:a=np.kron(a,paulis[k])
            mats.append(a)
        unitary=np.stack([a.reshape(-1,order='F')/math.sqrt(dim) for a in mats],axis=1)
        gram=float(np.max(np.abs(unitary.conj().T@unitary-np.eye(dim*dim))))
        x=np.kron(paulis[1],np.eye(2**(n-1)))
        L=unitary.conj().T@np.kron(np.eye(dim),x)@unitary
        R=unitary.conj().T@np.kron(x.T,np.eye(dim))@unitary
        height=np.array([max([i+1 for i,k in enumerate(w) if k]+[0]) for w in words])
        K=np.diag(height**2)
        comm=L@R-R@L
        nested=(K@L-L@K)@R-R@(K@L-L@K)
        omega=np.eye(dim*dim)[:,0]
        obstruction=float(np.linalg.norm(nested@omega))
        direct=float(np.max(np.abs(nested@omega+2*omega)))
        assert gram<1e-14 and np.max(np.abs(comm))<1e-14
        assert direct<1e-14 and abs(obstruction-2)<1e-14
        heats=[]
        for beta in (.3,.7,1.2):
            exact=float(np.exp(-beta*height**2).sum())
            counted=1+sum(3*4**(ell-1)*math.exp(-beta*ell*ell) for ell in range(1,n+1))
            assert abs(exact-counted)<2e-13
            heats.append(dict(beta=beta,matrix_trace=exact,degeneracy_formula=counted))
        rows.append(dict(spins=n,standard_Hilbert_dimension=dim*dim,gram_error=gram,
            left_right_commutator=float(np.max(np.abs(comm))),
            auxiliary_flow_nonlocal_double_commutator_on_vacuum=obstruction,
            exact_minus_two_vector_residual=direct,heat_checks=heats))
    tails=[]
    for beta in (.3,.7,1.2):
        m=16
        # ratio t_(ell+1)/t_ell=4*exp[-beta*(2ell+1)] decreases.
        total=1+sum(3*4**(ell-1)*math.exp(-beta*ell*ell) for ell in range(1,m+1))
        nextterm=3*4**m*math.exp(-beta*(m+1)**2)
        ratio=4*math.exp(-beta*(2*(m+1)+1))
        upper=nextterm/(1-ratio)
        assert upper>0 and ratio<1
        tails.append(dict(beta=beta,partial_trace=total,analytic_remaining_tail_upper=upper))
    return dict(date='2026-10-06',working_round=878,formal_reports=877,
        cumulative_numbered_groups=3662,fresh_numbered_groups=0,all_checks_passed=True,
        scope='Finite left/right Pauli checks of the explicit tracial infinite tensor-product counterexample; no claim about failure of the original physical E net.',
        finite_level_checks=rows,infinite_heat_trace_bounds=tails,
        abstract_no_split_proof_is_analytic=True,
        actual_E_split_or_failure_proved=False,full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
