"""1156 finite matrix certificates; no physical locality inferred from a gate list.
Run --write-once only for first result. Default recomputes and checks saved values.
"""
from pathlib import Path
import argparse, json
import numpy as np

def gate_perm(n, fn):
    U=np.zeros((2**n,2**n),complex)
    for j in range(2**n):
        bits=[(j>>(n-1-k))&1 for k in range(n)]
        out=fn(bits.copy()); i=sum(v<<(n-1-k) for k,v in enumerate(out))
        U[i,j]=1
    return U
def flip(bits, ctrls, target):
    if all(bits[c] for c in ctrls): bits[target]^=1
    return bits
def check():
    tol=2e-13
    V=np.zeros((32,4),complex)
    for x in range(2):
        for y in range(2):
            for r in range(2):
                a,b,e=r,r^(x*y),4*x+2*y+r
                V[(2*a+b)*8+e,2*x+y]=1/np.sqrt(2)
    U=np.eye(32,dtype=complex)
    for c,t in [(0,2),(1,3),(2,0),(3,1)]:
        U=gate_perm(5,lambda z,c=c,t=t:flip(z,[c],t))@U
    Had=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
    U=np.kron(np.eye(16),Had)@U
    for c,t in [(4,0),(4,1)]:
        U=gate_perm(5,lambda z,c=c,t=t:flip(z,[c],t))@U
    U=gate_perm(5,lambda z:flip(z,[2,3],1))@U
    unity=float(np.max(abs(U.conj().T@U-np.eye(32))))
    dilation=float(np.max(abs(U[:,[0,8,16,24]]-V)))
    # A non-scalar spectator clock checks tensor identities, not old-clock reconstruction.
    Hc=np.diag([1,1.1]); H=np.kron(np.eye(32),Hc)+6*np.eye(64)
    W=np.kron(U,np.eye(2))
    p0=H; pz=np.zeros_like(H)
    gamma,kappa=1.25,.75
    p0p=gamma*p0+kappa*pz; pzp=kappa*p0+gamma*pz
    charge=max(float(np.max(abs(W@p-p@W))) for p in [p0,pz,p0p,pzp])
    assert max(unity,dilation,charge)<tol

    # Complete matrix-unit checks imply the identities also with arbitrary passive R.
    ns=0.
    for i in range(4):
        for j in range(4):
            z=np.outer(V[:,i],V[:,j].conj()).reshape(2,2,8,2,2,8)
            target=np.eye(2)/2 if i==j else np.zeros((2,2))
            ns=max(ns,float(np.max(abs(np.einsum('abecbe->ac',z)-target))),
                   float(np.max(abs(np.einsum('abeade->bd',z)-target))))
    # Full E_x prematurely accessible: normal coherent input |+>, B input y.
    premature=[]
    for y in range(2):
        psi=(V[:,y]+V[:,2+y])/np.sqrt(2)
        z=psi.reshape(2,2,2,2,2)
        ex=np.einsum('abxyz,abwyz->xw',z,z.conj())
        plus=np.ones(2)/np.sqrt(2)
        premature.append(float(np.real(plus@ex@plus)))
    bad_ready_B1=[]
    for x in range(2):
        j=(2*x+1)*8
        psi=(U[:,j]+U[:,j+1])/np.sqrt(2)
        prob=np.sum(abs(psi.reshape(2,2,8))**2,axis=(0,2))
        bad_ready_B1.append(float(prob[1]))
    assert abs(bad_ready_B1[0])<tol and abs(bad_ready_B1[1]-1)<tol
    # Both semicausal decompositions preserve the entire E output.
    left=np.zeros_like(V); right=np.zeros_like(V)
    for x in range(2):
        for y in range(2):
            for r in range(2):
                left[(2*r+(r^(x*y)))*8+4*x+2*y+r,2*x+y]+=1/np.sqrt(2)
            for b in range(2):
                a=b^(x*y)
                right[(2*a+b)*8+4*x+2*y+a,2*x+y]+=1/np.sqrt(2)
    factor=max(float(np.max(abs(V-left))),float(np.max(abs(V-right))))
    assert max(ns,factor)<tol
    assert abs(premature[0]-1)<tol and abs(premature[1]-.5)<tol

    probs=np.empty((2,2,2,2))
    for x in range(2):
        for y in range(2):
            probs[x,y]=np.sum(abs(V[:,2*x+y].reshape(2,2,8))**2,axis=2)
    correl=np.einsum('xyab,ab->xy',probs,np.array([[1,-1],[-1,1]]))
    chsh=float(correl[0,0]+correl[0,1]+correl[1,0]-correl[1,1])
    bound=(2-np.sqrt(2))/4
    # The standard quantum-optimal unbiased CHSH distribution saturates TV bound.
    quantum=np.empty_like(probs)
    for x in range(2):
        for y in range(2):
            for a in range(2):
                for b in range(2):
                    quantum[x,y,a,b]=(1+(-1)**(a+b+x*y)/np.sqrt(2))/4
    tv=.5*np.sum(abs(probs-quantum),axis=(2,3))
    X=np.array([[0,1],[1,0]],complex); Z=np.diag([1,-1])
    B0=(Z+X)/np.sqrt(2); B1=(Z-X)/np.sqrt(2)
    S=np.kron(Z,B0+B1)+np.kron(X,B0-B1)
    tsirelson=float(np.linalg.norm(S,2))
    assert abs(chsh-4)<tol and np.max(abs(tv-bound))<tol
    assert abs(tsirelson-2*np.sqrt(2))<tol

    # Interleaving: A side sees box1 then box2; B sees box2 then box1.
    # Geometry: sites z=+-2; local input/output events t=0,1,2,3; all cross
    # laboratory events spacelike. E ports open only at a common future t=6.
    vt=V.reshape(2,2,8,2,2)
    net=np.einsum('abexd,cdfay->cbefxy',vt,vt).reshape(256,4)
    seq=np.zeros_like(net)
    # Kept records of the first local outputs, made by ordinary local CNOTs.
    records=np.zeros((1024,4),complex) # axes a,c,d,b,e,f
    for x in range(2):
        for y in range(2):
            for a in range(2):
                for d in range(2):
                    b=a^(x*d); c=d^(a*y)
                    e=4*x+2*d+a; f=4*a+2*y+c
                    seq[((2*c+b)*8+e)*8+f,2*x+y]+=.5
                    idx=(((((a*2+c)*2+d)*2+b)*8+e)*8+f)
                    records[idx,2*x+y]+=.5
    contraction=float(np.max(abs(net-seq)))
    norm=float(np.max(abs(net.conj().T@net-np.eye(4))))
    norm_record=float(np.max(abs(records.conj().T@records-np.eye(4))))
    cross_ns=0.
    for i in range(4):
        for j in range(4):
            z=np.outer(records[:,i],records[:,j].conj()).reshape(4,4,64,4,4,64)
            target=np.eye(4)/4 if i==j else np.zeros((4,4))
            cross_ns=max(cross_ns,
                float(np.max(abs(np.einsum('abecbe->ac',z)-target))),
                float(np.max(abs(np.einsum('abeade->bd',z)-target))))
    assert max(contraction,norm,norm_record,cross_ns)<tol
    return dict(schema="round1156_checks_v1",status="PASS",check_groups=4,
      scope="Finite certificates support analytic full-type/cut proofs; no empirical data, no local construction of W.",
      complete_boundary=dict(unitary_residual=unity,ready_dilation_residual=dilation,
                             spectator_clock_fourP_residual=charge,
                             old_actual_clock="reused analytically; the 2-level clock here is diagnostic only"),
      complete_causality=dict(matrix_unit_residual=ns,dual_factorization_residual=factor,
                             prematurely_accessible_E_x_plus_probabilities=premature,
                             unadmitted_plus_ready_B_one_probabilities=bad_ready_B1,
                             early_E_scope="countertest to an unadopted enlarged permission, not legal early output"),
      finite_separation=dict(PR_CHSH=chsh,quantum_CHSH_norm=tsirelson,
                             per_context_TV_lower_bound=float(bound),
                             saturating_quantum_TV=tv.tolist()),
      interleaved_two_calls=dict(contraction_residual=contraction,isometry_residual=norm,
           kept_local_records_isometry_residual=norm_record,
           all_matrix_units_local_marginal_residual=cross_ns,
           full_E_outputs_retained=True,no_physical_spacelike_memory_wire_claim=True))
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    elif isinstance(a,float): assert abs(a-b)<1e-11
    else: assert a==b
if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--write-once",action="store_true")
    args=ap.parse_args(); result=check(); out=Path(__file__).with_name("results.json")
    if args.write_once:
        with out.open("x",encoding="utf-8") as f: json.dump(result,f,ensure_ascii=False,indent=2); f.write("\n")
    else: compare(result,json.loads(out.read_text(encoding="utf-8")))
    print(json.dumps(result,ensure_ascii=False,indent=2))
