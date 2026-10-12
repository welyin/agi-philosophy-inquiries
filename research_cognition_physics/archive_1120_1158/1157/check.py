"""1157 finite certificates. Default recomputes saved results without writing.
Algebraic and circuit checks only; no claim that matrices certify actual apparatus.
"""
from pathlib import Path
import argparse
import json
import numpy as np

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
H = (X + Z) / np.sqrt(2)


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def maxabs(a):
    return float(np.max(np.abs(a)))


def weyl(d):
    return np.roll(np.eye(d, dtype=complex), 1, axis=0), np.diag(np.exp(2j*np.pi*np.arange(d)/d))


def adjoint(kraus, b):
    return sum(k.conj().T @ b @ k for k in kraus)


def trace_b(k, da_out, da_in, db):
    return np.einsum('abcb->ac', k.reshape(da_out, db, da_in, db))


def embed(u, sites, n):
    rest = [j for j in range(n) if j not in sites]
    order = list(sites) + rest
    inverse = np.argsort(order)
    raw = np.kron(u, np.eye(2**len(rest))).reshape([2]*(2*n))
    return raw.transpose(list(inverse)+list(inverse+n)).reshape(2**n, 2**n)


def cnot():
    return np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]],complex)


def trace_distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2)


def exact_certificate():
    # A:2 -> A':3, two Kraus operators; B consists of a qubit and a qutrit.
    rng = np.random.default_rng(1157)
    q, _ = np.linalg.qr(rng.normal(size=(6,2))+1j*rng.normal(size=(6,2)))
    local = [q[:3], q[3:]]
    mix = np.array([[1,1j],[1j,1]],complex)/np.sqrt(2)
    local = [sum(mix[i,j]*local[j] for j in range(2)) for i in range(2)]
    kraus = [np.kron(k,np.eye(6)) for k in local]
    x3,z3 = weyl(3)
    gens = [np.kron(X,np.eye(3)),np.kron(Z,np.eye(3)),np.kron(I,x3),np.kron(I,z3)]
    tp = maxabs(sum(k.conj().T@k for k in kraus)-np.eye(12))
    defect = max(opnorm(adjoint(kraus,np.kron(np.eye(3),g))-np.kron(I,g)) for g in gens)
    factors = [trace_b(k,3,2,6)/6 for k in kraus]
    factor_error = max(maxabs(k-np.kron(l,np.eye(6))) for k,l in zip(kraus,factors))
    # A common Stinespring keeps the same Kraus environment, not just the marginal.
    full = np.vstack(kraus)
    reconstructed = np.vstack([np.kron(l,np.eye(6)) for l in factors])
    env_error = maxabs(full-reconstructed)
    all_units_error = 0.
    for i in range(6):
        for j in range(6):
            b=np.zeros((6,6),complex); b[i,j]=1
            all_units_error=max(all_units_error,maxabs(adjoint(kraus,np.kron(np.eye(3),b))-np.kron(I,b)))
    assert max(tp,defect,factor_error,env_error,all_units_error)<2e-12
    return dict(input_local_dim=2,output_local_dim=3,remote_dims=[2,3],generators=4,
                tp_residual=tp,generator_residual=defect,kraus_factor_residual=factor_error,
                kept_environment_residual=env_error,all_remote_matrix_units_residual=all_units_error)


def weak_condition_failure():
    u=np.diag([1,1,1,-1]).astype(complex)
    bz=np.kron(I,Z); bx=np.kron(I,X)
    faithful=maxabs(u@(np.eye(4)/4)@u.conj().T-np.eye(4)/4)
    zdef=opnorm(u.conj().T@bz@u-bz)
    xdef=opnorm(u.conj().T@bx@u-bx)
    fixed_remote_ready=0.
    for i in range(2):
        for j in range(2):
            a=np.zeros((2,2),complex); a[i,j]=1
            out=u@np.kron(a,I/2)@u.conj().T
            b=np.einsum('abad->bd',out.reshape(2,2,2,2))
            fixed_remote_ready=max(fixed_remote_ready,maxabs(b-(I/2 if i==j else np.zeros((2,2)))))
    plus=np.array([1,1],complex)/np.sqrt(2)
    state=np.kron(np.array([0,1]),plus)
    changed=u@state
    initial_x=float(np.vdot(state,bx@state).real)
    final_x=float(np.vdot(changed,bx@changed).real)
    assert max(faithful,zdef,fixed_remote_ready)<1e-13 and abs(xdef-2)<1e-13
    assert abs(initial_x-1)<1e-13 and abs(final_x+1)<1e-13
    return dict(faithful_ready_invariance=faithful,diagonal_readout_invariance=zdef,
                remote_maximally_mixed_ready_residual=fixed_remote_ready,
                full_X_generator_defect=xdef,remote_X_before_after=[initial_x,final_x],
                scope='Fixed ready and commuting readouts fail; full-state fixed marginal was NOT weakened to this.')


def joint_network():
    # Axes q,a,b,rA,rB. q is unknown; a,b share a common-past Bell state.
    n=5
    prep=np.zeros((32,2),complex)
    for q in range(2):
        for a in range(2):
            prep[16*q+8*a+4*a,q]=1/np.sqrt(2)
    ua=embed(cnot(),[1,3],n)@embed(cnot(),[0,1],n)
    ub=embed(cnot(),[2,4],n)@embed(H,[2],n)
    # C is a declared common-future event, so both records may be used there.
    uc=embed(cnot(),[4,0],n)@embed(np.diag([1,1,1,-1]),[3,4],n)
    va=uc@ub@ua@prep
    vb=uc@ua@ub@prep
    reorder=maxabs(va-vb)
    isometry=maxabs(va.conj().T@va-np.eye(2))
    comm=maxabs(ua@ub-ub@ua)
    # Unknown q entangled with a passive R. All records and q remain outputs.
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    pa=np.kron(va,I)@bell; pb=np.kron(vb,I)@bell
    ref_error=trace_distance(np.outer(pa,pa.conj()),np.outer(pb,pb.conj()))
    # Before convergence, Bell entanglement across the two laboratories survives
    # their coherent records. A=q,a,rA; B=b,rB.
    p=(ub@ua@prep[:,0]).reshape([2]*5).transpose(0,1,3,2,4).reshape(8,4)
    rho=np.outer(p.reshape(-1),p.reshape(-1).conj()).reshape(8,4,8,4)
    pt=rho.transpose(0,3,2,1).reshape(32,32)
    negativity=float((np.sum(abs(np.linalg.eigvalsh(pt)))-1)/2)
    # Event positions use c=1; all transport edges are future timelike.
    events={'P':(-4.,0.),'A':(0.,-2.),'B':(0.,2.),'C':(4.,0.)}
    margins=[]
    for x,y in [('P','A'),('P','B'),('A','C'),('B','C')]:
        tx,xx=events[x]; ty,xy=events[y]
        margins.append((ty-tx)-abs(xy-xx))
    assert max(reorder,isometry,comm,ref_error)<2e-12
    assert abs(negativity-.5)<2e-12 and min(margins)>0
    return dict(reordering_residual=reorder,isometry_residual=isometry,
                spacelike_event_commutator=comm,reference_trace_distance=ref_error,
                full_record_entanglement_negativity=negativity,min_causal_margin=min(margins),
                scope='Logical finite network certificate; actual clocks, control, and fourP require separate material evidence.')


def finite_error():
    rows=[]
    for theta in [1e-6,1e-4,1e-2]:
        u=np.diag(np.exp(-1j*theta*np.array([1,-1,-1,1])))
        ex=opnorm(u.conj().T@np.kron(I,X)@u-np.kron(I,X))
        ez=opnorm(u.conj().T@np.kron(I,Z)@u-np.kron(I,Z))
        eta=np.sqrt(2*ex)+np.sqrt(2*ez)
        twirl=sum(np.kron(I,g.conj().T)@u@np.kron(I,g) for g in [I,X,Z,X@Z])/4
        s=trace_b(twirl,2,2,2)/2
        left,_,right=np.linalg.svd(s)
        vlocal=np.kron(left@right,I)
        dist=opnorm(u-vlocal)
        assert eta<1 and opnorm(u-twirl)<=eta+1e-12 and dist<=2*eta+1e-12
        # Exact half-diamond distance to the identity is sin(theta) in this range.
        truth=float(np.sin(theta))
        rows.append(dict(theta=theta,epsilon_X=ex,epsilon_Z=ez,eta=float(eta),
                         isometry_distance=dist,half_diamond_to_identity=truth,
                         certified_half_diamond_bound=float(min(1,2*eta))))
    # Two independently declared local implementation errors; retained R included.
    alpha,beta=.03,.04
    rz=np.diag(np.exp(-1j*alpha*np.array([1,-1])))
    rx=np.cos(beta)*I-1j*np.sin(beta)*X
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    actual=np.kron(rx@rz,I)@bell
    td=trace_distance(np.outer(actual,actual.conj()),np.outer(bell,bell.conj()))
    bound=float(np.sin(alpha)+np.sin(beta))
    assert td<=bound+1e-12
    return dict(twirl_certificates=rows,composed_reference_trace_distance=td,
                telescoping_bound=bound,
                scope='Existence error bound alone does not certify physical access to the approximating map.')


def check():
    return dict(schema='round1157_checks_v1',status='PASS',check_groups=4,
                scope='Finite algebraic and network checks; no new empirical calibration or actual-apparatus certification.',
                exact_certificate=exact_certificate(),weak_condition_failure=weak_condition_failure(),
                joint_network=joint_network(),finite_error=finite_error())


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float): assert abs(a-b)<1e-11,(a,b)
    else: assert a==b,(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write-once',action='store_true')
    args=parser.parse_args(); result=check(); path=Path(__file__).with_name('results.json')
    if args.write_once:
        if path.exists(): raise FileExistsError(path)
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else: compare(result,json.loads(path.read_text(encoding='utf-8')))
    print(json.dumps(result,ensure_ascii=False,indent=2))
