"""815: test the state-specific and algebraic record-retention criteria.
Original neutral coefficients are retained in one calibration; all finite boson
and preparation choices are diagnostic, not an evaluation of the continuum W.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'803'))
import quasifree_pairing_probe as old
TARGET=HERE/'record_retention_criterion_results.json'

def system(n):
    a=old.car(n);xi=a+[b.conj().T for b in a];dim=2**n
    def Q(k):
        return sum((xi[i].conj().T@xi[j]*k[i,j]/2 for i in range(2*n) for j in range(2*n)),np.zeros((dim,dim),complex))
    def covariance(v):
        return np.array([[np.vdot(v,a@b.conj().T@v) for b in xi] for a in xi])
    return a,xi,Q,covariance

def original_kernel_check():
    car,xi,Q,cov=system(4);ids=[24,25,30,31,56,57,62,63]
    take=lambda a:a[np.ix_(ids,ids)]
    phi=np.array([0.,.6987151138647895,0.,0.,.545863690226873])
    later=phi*np.array([1.,.8,1.,1.,1.1])
    _,v=np.linalg.eigh(Q(take(old.vertex.mass(phi))));state=v[:,0]
    e,v=np.linalg.eigh(Q(take(old.vertex.mass(later))))
    state=(v*np.exp(-.63j*e))@v.conj().T@state
    P=cov(state);N=np.eye(8)-P
    f=(np.eye(8)[:,2]+1j*np.eye(8)[:,7])/np.sqrt(2)
    cf=np.r_[f[4:].conj(),f[:4].conj()]
    D=-np.outer(f,f.conj())+np.outer(cf,cf.conj())
    keys=[take(old.vertex.dmass(later,np.eye(5)[i])) for i in (1,4)]
    L=[1j*(D@k-k@D) for k in keys]
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
    B=[X,Y];vac=np.array([1.,0.]);W=np.array([[np.vdot(vac,a@b@vac) for b in B] for a in B])
    change=sum((np.kron(b,Q(l)) for b,l in zip(B,L)),np.zeros((32,32),complex))
    vector=change@np.kron(vac,state)
    rows=[]
    for j in range(2):
        Lj=sum((W[j,a]*L[a] for a in range(2)),np.zeros((8,8),complex))
        operator=Q(Lj);contracted=operator@state
        direct=np.kron((vac.conj()@B[j]).reshape(1,2),np.eye(16))@vector
        mean=np.vdot(state,contracted);trace_mean=np.trace(N@Lj)/2
        pair=np.linalg.norm(P@Lj@N,'fro')**2/2
        predicted=abs(mean)**2+pair
        measured=np.linalg.norm(contracted)**2
        reverse=np.linalg.norm(N@Lj@P,'fro')**2/2
        row=dict(test=j,contraction_error=float(np.max(abs(direct-contracted))),
            mean_trace_error=float(abs(mean-trace_mean)),square_error=float(abs(predicted-measured)),
            mean_sector=float(abs(mean)**2),pair_sector=float(pair),measured=float(measured),
            wrong_order_pair_defect=float(abs(reverse-pair)),
            commutator_rank=int(np.linalg.matrix_rank(Lj,tol=1e-12)))
        assert max(row[k] for k in ('contraction_error','mean_trace_error','square_error'))<1e-12
        assert row['commutator_rank']<=4
        rows.append(row)
    assert max(r['wrong_order_pair_defect'] for r in rows)>1e-5
    return dict(cases=rows,change_square=float(np.linalg.norm(vector)**2),
        finite_boson_Gram_rank=int(np.linalg.matrix_rank(W)),
        uses_original_neutral_mass_coefficients=True,original_continuum_W_evaluated=False)

def null_distinctions():
    a,xi,Q,cov=system(2);eye=np.eye(4);p=eye-a[0].conj().T@a[0]
    h=np.array([[0,1],[1,0]],complex)
    K=np.block([[h,np.zeros((2,2))],[np.zeros((2,2)),-h.T]])
    D=np.diag([-1.,0.,1.,0.]);L=1j*(D@K-K@D)
    operator=1j*(p@Q(K)-Q(K)@p)
    vacuum=np.eye(4)[:,0];occupied=np.eye(4)[:,1]
    P=cov(vacuum);mean=np.vdot(vacuum,operator@vacuum)
    pair=np.linalg.norm(P@L@(np.eye(4)-P),'fro')**2/2
    zero=float(np.linalg.norm(operator@vacuum)**2)
    nonzero=float(np.linalg.norm(operator@occupied)**2)
    assert zero<1e-14 and pair<1e-14 and abs(mean)<1e-14 and nonzero>.5
    X=np.array([[0,1],[1,0]],complex)
    redundant=np.kron(X,operator)+np.kron(X,-operator)
    assert np.linalg.norm(redundant)==0
    return dict(state_null=dict(mean_sector=float(abs(mean)**2),pair_sector=float(pair),
            actual_square=zero,other_state_square=nonzero,operator_norm=float(np.linalg.norm(operator,2))),
        redundant_field_components=dict(total_operator_norm=float(np.linalg.norm(redundant)),
            wrong_sum_of_individual_squares=float(2*np.linalg.norm(operator@occupied)**2)),
        first_order_is_not_all_orders=dict(first_order_generator_zero=True,
            second_order_generator_commutator_norm=float(np.linalg.norm(operator,2))))

def run():
    return dict(round=815,all_checks_passed=True,original_kernel=original_kernel_check(),
        null_distinctions=null_distinctions(),
        scope='State-vector, vacuum/pair-sector and algebraic-null criteria calibrated independently.',
        original_record_retention_time_computed=False,
        all_original_protected_records_classified=False,
        finite_coupling_or_autonomous_storage_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists(),'Do not overwrite frozen evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

