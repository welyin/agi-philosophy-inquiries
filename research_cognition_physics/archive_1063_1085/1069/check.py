"""Finite witnesses for round 1069. Default recomputes and compares, --write creates results."""
from pathlib import Path
import argparse,json
import numpy as np

def run():
    I=np.eye(2,dtype=complex); X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex); Z=np.diag([1,-1]).astype(complex)
    ket=[np.array([1,0],complex),np.array([0,1],complex)]
    p=[np.outer(v,v.conj()) for v in ket]
    qv=[np.array([3,4],complex)/5,np.array([-4,3],complex)/5]
    q=[np.outer(v,v.conj()) for v in qv]
    P=[np.kron(v,I) for v in p]; Q=[np.kron(v,I) for v in q]
    Vs=[X,Z]; Ms=[np.kron(q[s],Vs[s]) for s in range(2)]
    residuals=[]; counts={}
    def norm(a):return float(np.linalg.norm(a))
    def check(a,b):residuals.append(norm(a-b));assert residuals[-1]<2e-12
    def ptr_k(a):
        return np.trace(a.reshape(2,2,2,2),axis1=1,axis2=3)
    def real(x):
        assert abs(complex(x).imag)<1e-12
        return float(complex(x).real)
    def chan(mats,rho):return sum((m@rho@m.conj().T for m in mats),np.zeros_like(rho))
    for s in range(2):
        check(Ms[s].conj().T@Ms[s],Q[s]);check(Q[s]@Ms[s],Ms[s])
        for i in range(4):
            for j in range(4):
                e=np.zeros((4,4),complex);e[i,j]=1
                check(ptr_k(Ms[s]@e@Ms[s].conj().T),q[s]@ptr_k(e)@q[s])
    counts["instrument_matrix_units"]=32
    rates=[]
    for r in range(2):
        row=[]
        for s in range(2):
            U=np.outer(ket[r],qv[s].conj())+np.outer(ket[1-r],qv[1-s].conj())
            R=np.kron(U,Vs[s].conj().T)
            c=np.vdot(qv[s],ket[r])
            check(R.conj().T@R,np.eye(4))
            check(R@Ms[s]@P[r],c*P[r])
            row.append(real(c.conjugate()*c))
        rates.append(row)
    counts["reference_safe_recovery_operator_identities"]=4
    # Non-Luders memory change, followed by faithful recovery when r,s are known.
    before=np.kron(p[0],p[0])
    out=Ms[0]@before@Ms[0].conj().T/rates[0][0]
    memory_after=np.trace(out.reshape(2,2,2,2),axis1=0,axis2=2)
    check(memory_after,p[1])
    # Same separable mixed relation supplies both phase families.
    sing=np.array([0,1,-1,0],complex)/np.sqrt(2)
    w=.25;omega=(1-w)*np.eye(4)/4+w*np.outer(sing,sing.conj())
    decomp=.25*np.kron(I/2,I/2)
    for pauli in [X,Y,Z]:
        for sign in [-1,1]:
            decomp+=.125*np.kron((I+sign*pauli)/2,(I-sign*pauli)/2)
    check(omega,decomp)
    zeta=.6+.8j
    for basis in [p,q]:
        U=zeta*basis[0]+basis[1]
        UU=np.kron(U,U)
        check(UU@omega@UU.conj().T,omega)
        C=np.kron(basis[0],basis[1])@omega@np.kron(basis[1],basis[0])
        check(norm(C),w/2)
        # Scalar Schur coefficients preserve every in-sector operator.
        A0=P[0]+.6*P[1];A1=.8*P[1]
        for r in range(2):
            for i in range(2):
                for j in range(2):
                    e=np.zeros((2,2),complex);e[i,j]=1
                    rho=np.kron(p[r],e)
                    check(chan([A0,A1],rho),rho)
    dK=[p[0]+.6*p[1],.8*p[1]]
    degraded=chan([np.kron(a,b) for a in dK for b in dK],omega)
    check(degraded[1,2],-.125*.6**2)
    degradation=norm(degraded-omega)
    check(chan([np.kron(a,b) for a in dK for b in dK],np.eye(4)/4),np.eye(4)/4)
    # Retain private hardware but resource need not have full hardware support.
    total_resource=np.kron(omega,np.outer(np.array([1,0,0,0]),np.array([1,0,0,0])))
    assert np.linalg.matrix_rank(total_resource,tol=1e-10)==4
    # Two non-parallel phase generators; actual generated reading distinguishes Y.
    phase=1j*p[0]+p[1]
    rotated=phase.conj().T@q[0]@phase
    effects=[p[0],p[1],q[0],rotated]
    coeff=lambda a:np.array([real(np.trace(a@b))/2 for b in [I,X,Y,Z]])
    rank=int(np.linalg.matrix_rank(np.array([coeff(a) for a in effects])))
    diff_rank=int(np.linalg.matrix_rank(np.array([coeff(a-effects[1])[1:] for a in effects])))
    assert(rank,diff_rank)==(4,3)
    ygap=real(np.trace(rotated@Y));check(ygap,24/25)
    lie_vectors=np.array([[0,0,1],[24/25,0,-7/25],[0,24/25,0]])
    assert np.linalg.matrix_rank(lie_vectors)==3
    # Omitting recoverability: two genuine task sectors have distinct rates.
    multiP=np.kron(np.eye(2),p[0])
    qa=np.outer(np.array([3,4])/5,np.array([3,4])/5)
    qb=np.outer(np.array([4,3])/5,np.array([4,3])/5)
    multiQ=np.block([[qa,np.zeros((2,2))],[np.zeros((2,2)),qb]])
    recovery_failure_rates=[real(multiQ[0,0]),real(multiQ[2,2])]
    assert abs(recovery_failure_rates[0]-recovery_failure_rates[1])>.27
    # No repeatability: SWAP after Q is recoverable but leaks K into next direction read.
    S=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)
    bad=S@Q[0];bad_probs=[]
    for memory in p:
        out=bad@np.kron(p[0],memory)@bad.conj().T
        bad_probs.append(real(np.trace(P[0]@out)))
    check(bad_probs,np.array([9/25,0]))
    # Hidden fine record: coarse binary effect is sharp, individual effect reads K.
    fine=np.kron(q[0],p[0])
    projection=np.kron(ptr_k(fine)/2,I)
    fine_gap=norm(fine-projection)
    assert fine_gap>.7
    # Recovery may be extended maliciously outside its stated history domain.
    U=np.outer(ket[0],qv[0].conj())+np.outer(ket[1],qv[1].conj())
    controlled=np.kron(p[0],I)+np.kron(p[1],Z)
    extraR=controlled@np.kron(U,I)
    check(extraR@Q[0]@P[0],(3/5)*P[0])
    plus=(ket[0]+ket[1])/np.sqrt(2)
    start=U.conj().T@plus
    extra_probs=[]
    read=np.kron(np.outer(plus,plus.conj()),I)
    for memory in p:
        rho=np.kron(np.outer(start,start.conj()),memory)
        out=extraR@rho@extraR.conj().T
        extra_probs.append(real(np.trace(read@out)))
    check(extra_probs,np.array([1.,0.]))
    return {"round":1069,"passed":True,"checks":counts,"rates":rates,
      "resource_eigenvalues":np.linalg.eigvalsh(omega).tolist(),
      "resource_purity":real(np.trace(omega@omega)),
      "resource_separable_decomposition_verified":True,
      "resource_hardware_dimension":16,"resource_hardware_rank":4,
      "resource_cross_block_norm":w/2,"dephasing_cross_before":-.125,
      "dephasing_cross_after":real(degraded[1,2]),"dephasing_relation_defect":degradation,
      "effect_rank":rank,"traceless_difference_rank":diff_rank,"generated_Y_reading_gap":ygap,
      "without_recovery_rates":recovery_failure_rates,
      "without_repeatability_next_joint_probabilities":bad_probs,
      "hidden_fine_record_distance_to_task_algebra":fine_gap,
      "standalone_recovery_next_probabilities":extra_probs,
      "matrix_assertions":len(residuals),"max_formula_residual":max(residuals),
      "spatial_dimension_proved":False,
      "scope":"Finite witnesses; general CP, continuous control and task statements are in proof.md."}

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert abs(a-b)<=2e-12
    else:assert a==b,(a,b)

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args();result=run();p=Path(__file__).with_name("results.json")
    if args.write:
        with p.open("x",encoding="utf-8") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,json.loads(p.read_text(encoding="utf-8")))
    print(json.dumps(result,ensure_ascii=False))
