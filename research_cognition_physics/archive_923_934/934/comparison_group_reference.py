"""Finite comparison data, collective symmetry, and an internal reference witness.

The Clifford 3-design fact is established literature, not a new discovery.
The project test connects it to the prior relation-carrier/actual-direction gap.
"""
from pathlib import Path
import argparse, hashlib, itertools, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TARGET=HERE/'comparison_group_reference_results.json'
I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)

def kron_all(xs):
    out=np.array([[1.]],complex)
    for a in xs:out=np.kron(out,a)
    return out

def lift(a,n):return kron_all([a]*n)
def norm(a):return float(np.linalg.norm(a))
def scalar(z):
    assert abs(complex(z).imag)<2e-11
    return float(complex(z).real)

def binary_octahedral():
    b=np.eye(4);qs=[]
    for a in b:
        qs.extend([a,-a])
    for signs in itertools.product((-1,1),repeat=4):qs.append(np.array(signs)/2)
    for a,c in itertools.combinations(range(4),2):
        for s,t in itertools.product((-1,1),repeat=2):qs.append((s*b[a]+t*b[c])/np.sqrt(2))
    return [q[0]*I-1j*(q[1]*X+q[2]*Y+q[3]*Z) for q in qs]

def permutation_matrix(n,p):
    d=2**n;out=np.zeros((d,d),complex)
    for j in range(d):
        bits=[(j>>(n-1-k))&1 for k in range(n)]
        k=sum(bits[p[a]]<<(n-1-a) for a in range(n))
        out[k,j]=1
    return out

def schur_projector(n):
    perms=[permutation_matrix(n,p) for p in itertools.permutations(range(n))]
    mat=np.column_stack([a.reshape(-1,order='F') for a in perms])
    u,s,_=np.linalg.svd(mat,full_matrices=False)
    rank=int(np.count_nonzero(s>1e-10));basis=u[:,:rank]
    return basis@basis.conj().T,rank,sum(perms)/math.factorial(n)

def collective(a,n):
    return sum(kron_all([a/2 if j==k else I for j in range(n)]) for k in range(n))

def haar_state_twirl(rho,n):
    """Euler Haar integral: exact z dephasings and polynomial Gauss quadrature.

    After both z averages, beta integrands have degree <= n in cos(beta).
    n+1 Legendre nodes integrate that finite polynomial space exactly.
    This evaluates a finite density matrix; it is not an infinite apparatus.
    """
    counts=np.array([j.bit_count() for j in range(2**n)])
    mask=counts[:,None]==counts[None,:]
    mid=rho*mask
    eig,v=np.linalg.eigh(collective(Y,n))
    xs,ws=np.polynomial.legendre.leggauss(n+1)
    out=np.zeros_like(rho)
    for x,w in zip(xs,ws):
        u=(v*np.exp(-1j*np.arccos(x)*eig))@v.conj().T
        out+=(w/2)*(u@mid@u.conj().T)*mask
    return (out+out.conj().T)/2

def run():
    group=binary_octahedral();assert len(group)==48
    stack=np.stack(group)
    closure=max(min(norm(a@b-c) for c in group) for a in group for b in group)
    assert closure<2e-14
    rows=[]
    for n in range(1,5):
        blocks=[lift(u,n) for u in group]
        pg=sum(np.kron(u.conj(),u) for u in blocks)/48
        ph,rank,psym=schur_projector(n)
        expected=(2*4**n+12*2**n+16)/48
        dist=norm(pg-ph)
        assert abs(np.trace(pg).real-expected)<2e-12
        assert rank==[1,2,5,14][n-1]
        assert abs(dist-(0 if n<4 else 1))<2e-12
        rows.append(dict(copies=n,finite_commutant_dimension=round(expected),continuous_commutant_dimension=rank,
                         full_channel_frobenius_difference=dist,finite_projector_error=norm(pg@pg-pg)))
    n=4;z4=np.zeros(16,complex);z4[0]=1
    seed=np.outer(z4,z4.conj())
    finite=sum(lift(u,4)@seed@lift(u,4).conj().T for u in group)/48
    continuous=psym/5
    q=lift(X,4)+lift(Y,4)+lift(Z,4)
    effect=(np.eye(16)+q)/4
    assert norm(effect@effect-effect)<1e-13
    calibrated=[scalar(np.trace(effect@rho)) for rho in (finite,continuous)]
    assert np.allclose(calibrated,[.5,.4],atol=2e-14,rtol=0)
    purity=scalar(np.trace(finite@finite));overlap=scalar(np.trace(continuous@finite))
    assert abs(purity-5/24)<1e-14 and abs(overlap-1/5)<1e-14
    # Keep the four-qubit reference inside the same object, erase only a common axis.
    omega_f=haar_state_twirl(np.kron(finite,finite),8)
    omega_h=haar_state_twirl(np.kron(continuous,finite),8)
    swap=permutation_matrix(8,tuple(range(4,8))+tuple(range(4)))
    read=(np.eye(256)+swap)/2
    internal=[]
    for rho in (omega_f,omega_h):
        a=rho.reshape(16,16,16,16)
        left=np.einsum('abcb->ac',a);right=np.einsum('abad->bd',a)
        covariance=max(norm(collective(p,8)@rho-rho@collective(p,8)) for p in (X,Y,Z))
        # Full pointer density, controlled block-SWAP, then the same pointer + readout.
        initial=np.kron(np.ones((2,2))/2,rho)
        perm=np.argmax(swap,axis=0)
        idx=np.concatenate([np.arange(256),256+perm])
        final=initial[np.ix_(idx,idx)]
        pointer_p=scalar(np.trace(final[:256,:256]+final[256:,256:]+final[:256,256:]+final[256:,:256])/2)
        direct=scalar(np.trace(read@rho))
        item=dict(probability=direct,full_pointer_probability=pointer_p,
                  common_frame_commutator_error=covariance,trace_error=abs(scalar(np.trace(rho))-1),
                  smallest_eigenvalue=float(np.linalg.eigvalsh(rho)[0]),
                  source_marginal_error=norm(left-continuous),reference_marginal_error=norm(right-continuous))
        assert covariance<1e-12 and abs(pointer_p-direct)<1e-13
        assert max(item['source_marginal_error'],item['reference_marginal_error'])<1e-12
        assert item['smallest_eigenvalue']>-1e-12 and item['trace_error']<1e-12
        internal.append(item)
    assert abs(internal[0]['probability']-29/48)<1e-12
    assert abs(internal[1]['probability']-3/5)<1e-12
    gap=internal[0]['probability']-internal[1]['probability']
    assert abs(gap-1/240)<1e-12
    assert norm(omega_h-np.kron(continuous,continuous))<1e-12
    common_invariance=max(norm(swap@collective(p,8)-collective(p,8)@swap) for p in (X,Y,Z))
    assert common_invariance<1e-12
    # A continuous reorientation detects the extra fourth-order record.
    u=np.diag(np.exp(np.array([-1j,1j])*np.pi/8))
    q_change=norm(lift(u,4)@q@lift(u,4).conj().T-q)
    assert q_change>1
    inputs=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/role_algebra_reuse_review.md',
            HERE/'drafts/operation_hypotheses_screening.md',HERE.parent/'research_note_933.md',
            ROOT/'research_cognition_physics/archive_429_466/research_note_430.md',
            ROOT/'research_cognition_physics/archive_429_466/research_note_461.md']
    return dict(round=934,date='2026-10-07',all_scientific_checks_passed=True,
       mature_three_design_result_reused=True,one_bounded_group_and_internal_reference_test=True,
       group_size=48,closure_error=closure,moment_rows=rows,
       four_copy_calibrated_probabilities=calibrated,four_copy_purity=purity,
       fourth_record_continuous_rotation_difference=q_change,
       internal_reference_rows=internal,internal_probability_gap=gap,
       internal_swap_common_frame_invariance_error=common_invariance,
       continuous_joint_factorization_error=norm(omega_h-np.kron(continuous,continuous)),
       error_example=dict(each_total_probability_error=.0005,guaranteed_gap=1/240-.001),
       physical_comparison_permissions_and_preparations_remain_inputs=True,
       all_six_protocols_realized=False,space_or_gauge_dynamics_derived=False,
       full_goal_completed=False,candidate_detail_optimization_stopped=True,
       source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
