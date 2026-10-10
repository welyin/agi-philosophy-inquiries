"""Round 1075: deterministic numerical checks; prints JSON, writes nothing."""
import json
import numpy as np
rng=np.random.default_rng(1075)
I=np.eye(2,dtype=complex)
Z=np.diag([1.,-1.]).astype(complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
S=np.eye(4,dtype=complex)[[0,2,1,3]]
Ps=(np.eye(4)+S)/2
Pa=(np.eye(4)-S)/2
singlet=np.array([0,1,-1,0],complex)/np.sqrt(2)
counts=0
residuals={}
def close(name,a,b,tol=4e-12):
    global counts
    counts+=1
    err=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    residuals[name]=max(residuals.get(name,0.),err)
    assert err<tol,(name,err)
def le(name,a,b,tol=4e-12):
    global counts
    counts+=1
    err=float(max(0.,float(a)-float(b)))
    residuals[name]=max(residuals.get(name,0.),err)
    assert err<tol,(name,a,b)
def pure(v):
    v=np.asarray(v,complex);v=v/np.linalg.norm(v)
    return np.outer(v,v.conj())
def state(n):
    a=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    a=a@a.conj().T
    return a/np.trace(a)
def vec():
    a=rng.normal(size=2)+1j*rng.normal(size=2)
    return a/np.linalg.norm(a)
def ptr2(w):
    return np.trace(w.reshape(2,2,2,2),axis1=1,axis2=3)
def ptr1(w):
    return np.trace(w.reshape(2,2,2,2),axis1=0,axis2=2)
def td(a,b):
    return float(np.sum(np.linalg.svd(a-b,compute_uv=False))/2)
def ks(tau):
    out=[np.kron(I,I[b:b+1])@Ps for b in range(2)]
    vals,basis=np.linalg.eigh(tau)
    out.extend(np.sqrt(max(0.,v))*np.outer(basis[:,j],singlet.conj())
               for j,v in enumerate(vals))
    return out
def apply(k,w):
    return sum(a@w@a.conj().T for a in k)
def choi(k):
    # input index first, then output, from matrix units
    return sum(np.kron(np.outer(np.eye(4)[i],np.eye(4)[j]),
               apply(k,np.outer(np.eye(4)[i],np.eye(4)[j])))
               for i in range(4) for j in range(4))
paulis=[(I+sign*a)/2 for a in [X,Y,Z] for sign in [1,-1]]
rho01=np.diag([0.,1.,0.,0.]).astype(complex)
rho10=np.diag([0.,0.,1.,0.]).astype(complex)
A=Ps@rho01@Ps
identity=sum(sign*np.kron(q,q) for sign,q in zip([1,1,1,1,-1,-1],paulis))/4
close("six_preparation_identity",A,identity)
close("paired_input_identity",(rho01+rho10)/2,A+Pa/2)
taus=[I/2,pure([1,0]),pure([1,1j])]+[state(2) for _ in range(6)]
for tau in taus:
    k=ks(tau)
    close("trace_preservation",sum(a.conj().T@a for a in k),np.eye(4))
    le("choi_complete_positivity",0.,np.linalg.eigvalsh(choi(k)).min())
    for i in range(4):
        for j in range(4):
            w=np.outer(np.eye(4)[i],np.eye(4)[j])
            explicit=(ptr1(w)+ptr2(w))/2+np.trace(Pa@w)*(tau-I/2)
            close("all_matrix_units_classification",apply(k,w),explicit)
    for _ in range(12):
        u,v=vec(),vec();r,t=pure(u),pure(v);w=np.kron(r,t)
        c=float(abs(np.vdot(u,v)));gap=(1-c)**2/4
        out=apply(k,w)
        close("pure_agreement",apply(k,np.kron(r,r)),r)
        close("role_exchange",out,apply(k,S@w@S))
        asym=ptr2(Ps@w@Ps)
        close("forced_spectrum",np.linalg.eigvalsh(asym),[(1-c)**2/4,(1+c)**2/4])
        closest=1-float(np.linalg.eigvalsh(out)[-1])
        le("all_pure_target_lower_bound",gap,closest)
        for _ in range(3):
            le("sample_pure_targets",gap,td(out,pure(vec())))
        vals,us=np.linalg.eigh(asym)
        q=pure(us[:,-1])
        close("pointwise_sharpness",td(apply(ks(q),w),q),gap)
        meanout=apply(ks(I/2),w)
        close("covariant_closest_pure",1-np.linalg.eigvalsh(meanout)[-1],(1-c)/2)
    # passive entangled reference (dimension 3) complete extension
    joint=state(12)
    ext=sum(np.kron(a,np.eye(3))@joint@np.kron(a,np.eye(3)).conj().T for a in k)
    le("reference_positivity",0.,np.linalg.eigvalsh(ext).min())
    close("reference_trace",np.trace(ext),1.)
# Pure-state projectors span all Hermitian operators on symmetric subspace.
T=np.array([[1,0,0],[0,1/np.sqrt(2),0],[0,1/np.sqrt(2),0],[0,0,1]],complex)
rows=[]
for _ in range(20):
    q=pure(vec())
    m=T.conj().T@np.kron(q,q)@T
    rows.append(np.r_[np.diag(m).real,m[np.triu_indices(3,1)].real,m[np.triu_indices(3,1)].imag])
close("projector_span_rank",np.linalg.matrix_rank(rows),9)
# Eight-input error criterion, with role asymmetry explicitly present.
base=ks(pure([1,0]))
discard=[np.kron(I,I[b:b+1]) for b in range(2)]
robust=[]
for alpha,beta in [(0.,0.),(.01,0.),(.04,.03),(.12,.10)]:
    reset=[np.outer(np.array([1,0]),np.eye(4)[j]) for j in range(4)]
    k=[np.sqrt(1-alpha-beta)*a for a in base]+[np.sqrt(alpha)*a for a in reset]+[np.sqrt(beta)*a for a in discard]
    delta=max(td(apply(k,np.kron(q,q)),q) for q in paulis)
    kappa=td(apply(k,rho01),apply(k,rho10))
    bound=max(0.,.25-1.5*delta-.5*kappa)
    minimum=1-np.linalg.eigvalsh(apply(k,rho01))[-1]
    le("finite_eight_input_bound",bound,minimum)
    robust.append(dict(reset_weight=alpha,role_bias_weight=beta,delta=delta,kappa=kappa,bound=bound,nearest_pure_error=float(minimum)))
# Directly witness symmetry need not be covariance.
h=(X+Z)/np.sqrt(2)
w=np.kron((I+X)/2,(I+Y)/2)
covgap=td(apply(base,np.kron(h,h)@w@np.kron(h,h).conj().T),h@apply(base,w)@h.conj().T)
le("role_symmetry_not_covariance",.01,covgap)
# Dropping roles, agreement, or purity gives the three stated witnesses.
close("first_party_pure",ptr2(rho01),pure([1,0]))
close("first_party_breaks_roles",td(ptr2(rho01),ptr2(rho10)),1.)
close("mean_orthogonal_mixed",apply(ks(I/2),rho01),I/2)
close("exact_orthogonal_quarter",1-np.linalg.eigvalsh(apply(base,rho01))[-1],.25)
print(json.dumps(dict(round=1075,passed=True,assertion_calls=counts,
 residuals=residuals,max_residual=max(residuals.values()),
 finite_eight_input_examples=robust,noncovariance_gap=covgap,
 empirical_results=0,spatial_generation_goal_completed=False),indent=2))
