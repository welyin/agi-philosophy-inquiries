"""877: finite diagnostic of collared type-I compression, not a PDE/split proof."""
from pathlib import Path
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'collared_region_compression_results.json'

def tnorm(x):
    return float(np.linalg.svd(x,compute_uv=False).sum())
def partial(x,a,b,which):
    z=x.reshape(a,b,a,b)
    return np.einsum('ibjb->ij',z) if which=='local' else np.einsum('aiaj->ij',z)
def maps(r,m=6,d=2,blocks=2):
    loc=m*d*blocks
    P=np.zeros((loc,loc),complex);ks=[]
    for a in range(blocks):
        base=a*d*m
        for j in range(r):
            for i in range(d):P[base+i*m+j,base+i*m+j]=1
        for j in range(r,m):
            K=np.zeros_like(P)
            for i in range(d):K[base+i*m,base+i*m+j]=1
            ks.append(K)
    return P,[P]+ks
def apply(rho,ks,outside):
    return sum(np.kron(k,np.eye(outside))@rho@np.kron(k.conj().T,np.eye(outside)) for k in ks)
def run():
    rng=np.random.default_rng(877)
    m=6;d=2;blocks=2;loc=m*d*blocks;outside=3;dim=loc*outside
    X=rng.normal(size=(dim,5))+1j*rng.normal(size=(dim,5))
    rho=X@X.conj().T;rho/=np.trace(rho)
    records=[]
    for a in range(blocks):
        for i in range(d):
            for j in range(d):
                A=np.zeros((loc,loc),complex)
                A[a*d*m+i*m:a*d*m+(i+1)*m,a*d*m+j*m:a*d*m+(j+1)*m]=np.eye(m)
                records.append(A)
    Y=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim))
    H=(Y+Y.conj().T)/2;w,V=np.linalg.eigh(H)
    U=(V*np.exp(-.19j*w))@V.conj().T
    final=np.diag(np.arange(dim)%3==0).astype(complex)
    M=U.conj().T@final@U
    rows=[]
    for r in (1,2,4,6):
        P,ks=maps(r);state=apply(rho,ks,outside)
        defect=max(tnorm(sum(k.conj().T@A@k for k in ks)-A) for A in records)
        joint=0.
        for A in records:
            for i in range(outside):
                for j in range(outside):
                    B=np.zeros((outside,outside),complex);B[i,j]=1
                    joint=max(joint,abs(np.trace((state-rho)@np.kron(A,B))))
        eps=float(np.trace(rho@np.kron(np.eye(loc)-P,np.eye(outside))).real)
        dist=tnorm(state-rho)
        tail=abs(np.trace(state@np.kron(np.eye(loc)-P,np.eye(outside))))
        outside_defect=tnorm(partial(state,loc,outside,'remote')-partial(rho,loc,outside,'remote'))
        complete=tnorm(sum(k.conj().T@k for k in ks)-np.eye(loc))
        event_error=abs(np.trace((state-rho)@M))
        assert defect<1e-12 and joint<1e-12 and outside_defect<1e-12
        assert tail<1e-12 and complete<1e-12 and abs(np.trace(state)-1)<1e-12
        assert np.linalg.eigvalsh(state).min()>-1e-12
        assert dist<=2*np.sqrt(eps)+eps+1e-12 and event_error<=dist/2+1e-12
        rows.append(dict(local_capacity=4*r,tail_probability=eps,state_norm_error=dist,
            bound=2*np.sqrt(eps)+eps,record_algebra_error=defect,
            record_outside_correlation_error=float(joint),outside_marginal_error=outside_defect,
            output_support_error=float(tail),trace_preservation_error=complete,
            fixed_later_bounded_event_error=float(event_error)))
    P,ks=maps(2)
    Ptot=np.kron(P,np.eye(outside))
    post=Ptot@rho@Ptot;post/=np.trace(post)
    post_defect=tnorm(partial(post,loc,outside,'remote')-partial(rho,loc,outside,'remote'))
    wrong=[P]
    omitted=np.flatnonzero(np.diag(P).real==0)
    for j in omitted:
        k=np.zeros_like(P);k[0,j]=1;wrong.append(k)
    record_loss=max(tnorm(sum(k.conj().T@A@k for k in wrong)-A) for A in records)
    assert post_defect>1e-3 and record_loss>1
    # A genuinely infinite Schmidt-rank input has epsilon=q**r.
    q=.6
    infinite=[dict(rank=r,tail=q**r,exact_norm_error=float(np.sqrt(4*q**r-3*q**(2*r))+q**r),
                   bound=float(2*np.sqrt(q**r)+q**r)) for r in (1,2,4,8,16)]
    # Independent finite matrix verification of the same closed formula.
    N=12
    probs=(1-q)*q**np.arange(N);probs/=probs.sum()
    psi=np.zeros(N*N,complex);psi[np.arange(N)*(N+1)]=np.sqrt(probs)
    rho2=np.outer(psi,psi.conj());checks=[]
    for r in (1,3,7):
        p=np.diag(np.arange(N)<r).astype(complex);operators=[p]
        for n in range(r,N):
            k=np.zeros((N,N),complex);k[0,n]=1;operators.append(k)
        y=apply(rho2,operators,N);epsilon=float(probs[r:].sum())
        expected=np.sqrt(4*epsilon-3*epsilon**2)+epsilon
        error=abs(tnorm(y-rho2)-expected)
        assert error<2e-12
        checks.append(dict(rank=r,formula_error=error))
    return dict(round=877,date='2026-10-06',formal_reports=877,
        cumulative_numbered_groups=3662,fresh_numbered_groups=1,all_checks_passed=True,
        scope='Matrix diagnostic of the normal collared-factor map; not a numerical proof of continuum splitting, energy control, or autonomous dynamics.',
        finite_center_and_record_checks=rows,
        postselection_changes_outside_marginal=post_defect,
        resetting_record_instead_of_only_multiplicity_defect=record_loss,
        infinite_Schmidt_input=infinite,independent_trace_distance_checks=checks,
        original_free_neutral_receiver_net_split_applicability_analytic=True,
        original_full_mixed_gravity_gauge_net_split_proved=False,
        original_unbounded_source_or_spatial_refinement_bridge_proved=False,
        full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
