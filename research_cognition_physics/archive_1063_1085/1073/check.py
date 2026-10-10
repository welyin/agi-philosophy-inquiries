"""Round 1073: finite Kraus/Choi and reference checks; no SDP or empirical claims."""
import json
import numpy as np

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1,-1]).astype(complex)
P=np.diag([1,0]).astype(complex)
rng=np.random.default_rng(1073)
checks=0
max_equality_residual=0.
def eq(a,b,tol=3e-11):
    global checks,max_equality_residual
    v=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    assert v<=tol,(v,a,b)
    max_equality_residual=max(max_equality_residual,v)
    checks+=1
def le(a,b,tol=3e-11):
    global checks
    assert float(a)<=float(b)+tol,(a,b)
    checks+=1
def psd(a):
    le(-np.linalg.eigvalsh((a+a.conj().T)/2)[0],0)
def trace_norm(a):
    return float(np.linalg.svd(a,compute_uv=False).sum())
def unitar():
    v=rng.normal(size=3);v/=np.linalg.norm(v); t=rng.uniform(-3,3)
    return np.cos(t)*I-1j*np.sin(t)*(v[0]*X+v[1]*Y+v[2]*Z)
def sqrt_psd(a):
    v,w=np.linalg.eigh(a)
    return (w*np.sqrt(np.maximum(v,0)))@w.conj().T
def adj(ks,a):
    return sum((k.conj().T@a@k for k in ks),np.zeros_like(I))
def out(ks,rho,ref=1):
    return sum((np.kron(k,np.eye(ref))@rho@np.kron(k,np.eye(ref)).conj().T for k in ks),
               np.zeros_like(rho))
def rho_rand(d):
    a=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d));r=a@a.conj().T
    return r/np.trace(r)
def choi(ks):
    return sum((np.outer(k.reshape(-1,order="F"),k.reshape(-1,order="F").conj())
                for k in ks),np.zeros((4,4),complex))
def reset_ks(effect,projector):
    ev,vec=np.linalg.eigh(effect)
    _,pv=np.linalg.eigh(projector);v=pv[:,-1]
    return [np.sqrt(max(0,ev[j]))*np.outer(v,vec[:,j].conj()) for j in range(2)]
def bound(eps):
    q=eps/(1-eps)
    return min(1.,float(np.sqrt(q)+q/2+eps))
def instrument(eps,rotation):
    # Uniformly stable reset part plus coherent Kraus terms, all resources explicit.
    p=rotation@P@rotation.conj().T
    e=eps/4
    effect=e*I+(1-2*e)*p
    t=0 if eps==0 else (eps-e)/(1-2*e)
    result=[]
    for es,ps in [(effect,p),(I-effect,I-p)]:
        root=sqrt_psd(es)
        u,v=unitar(),unitar()
        ks=[np.sqrt(1-t)*k for k in reset_ks(es,ps)]
        ks += [np.sqrt(t/3)*u@root,np.sqrt(2*t/3)*v@root]
        result.append(ks)
    return result,p,effect
rows=[]
for eps in [0.,1e-4,.01,.1,.25,.49]:
    ins,p,effect=instrument(eps,unitar())
    q=eps/(1-eps)
    eq(adj(ins[0],I),effect);eq(adj(ins[1],I),I-effect)
    eq(adj(ins[0],I)+adj(ins[1],I),I)
    for ks,es,ps in zip(ins,[effect,I-effect],[p,I-p]):
        psd(choi(ks))
        psd(eps*es-adj(ks,I-es)) # all-input conditional inequality
    evals=np.linalg.eigvalsh(effect)
    le(evals[0],eps);le(1-evals[1],eps)
    le(abs(np.trace(effect).real-1),eps)
    le(np.linalg.norm(effect-p,2),eps)
    leaks=[adj(ins[0],I-p),adj(ins[1],p)]
    psd(q*I-leaks[0]-leaks[1])
    for ks,es,ps,loss in zip(ins,[effect,I-effect],[p,I-p],leaks):
        compressed=[ps@k for k in ks]
        mp=reset_ks(es,ps)
        cp_loss=choi(mp)-choi(compressed)
        psd(cp_loss)
        eq(cp_loss,choi(reset_ks(loss,ps)))
    largest_sample_error=0.
    for ref in [1,2,3]:
        for j in range(4):
            rho=rho_rand(2*ref)
            actual=[out(ks,rho,ref) for ks in ins]
            ideal=[out([ps],rho,ref) for ps in [p,I-p]]
            err=sum(trace_norm(a-b) for a,b in zip(actual,ideal))/2
            largest_sample_error=max(largest_sample_error,err)
            le(err,bound(eps))
            le(sum(np.trace(np.kron(I-ps,np.eye(ref))@a).real
                   for ps,a in zip([p,I-p],actual)),q)
    rows.append({"epsilon":eps,"half_diamond_upper":bound(eps),
                 "largest_sample_reference_distance":largest_sample_error})

# Exact square-root order example: Choi-block/ref direct computation.
sqrt_rows=[]
for eps in [0.,1e-6,1e-4,.01,.1,.49]:
    u=np.array([[np.sqrt(1-eps),-np.sqrt(eps)],
                [np.sqrt(eps),np.sqrt(1-eps)]],complex)
    ins=[[u@P],[u@(I-P)]]
    eq(adj(ins[0],I),P);eq(adj(ins[1],I),I-P)
    eq(adj(ins[0],I-P),eps*P)
    eq(adj(ins[1],P),eps*(I-P))
    for ref in [1,2,3]:
        r=rho_rand(2*ref)
        d=sum(trace_norm(out(k,r,ref)-out([ps],r,ref))
              for k,ps in zip(ins,[P,I-P]))/2
        eq(d,np.sqrt(eps))
    sqrt_rows.append({"epsilon":eps,"exact_half_diamond":float(np.sqrt(eps))})

# Bias threshold sharpness and per-branch repeatability, exact two-setting example.
threshold_rows=[]
for eps in [0.,.01,.1,.49]:
    e=(1-eps)*P;f=P+eps*(I-P)
    eq(f-e,eps*I)
    eq((f-e)-np.trace(f-e)*I/2,np.zeros((2,2)))
    for effect in [e,f]:
        ins=[reset_ks(effect,P),reset_ks(I-effect,I-P)]
        for ks,es in zip(ins,[effect,I-effect]):
            psd(eps*es-adj(ks,I-es))
    threshold_rows.append({"epsilon":eps,"effect_gap":float(np.linalg.norm(f-e,2)),
                           "traceless_gap":0.})

# A finite adaptive tree: same feedback rules in both comparisons.
eps=.0001
menus=[instrument(eps,unitar()) for _ in range(4)]
initial=rho_rand(4)
actual={():initial};ideal={():initial};tree=[]
for n in range(1,4):
    na={};ni={}
    for history in actual:
        idx=(sum(history)+n-1)%4
        ins,p,_=menus[idx]
        for s in [0,1]:
            key=history+(s,)
            na[key]=out(ins[s],actual[history],2)
            ni[key]=out([p if s==0 else I-p],ideal[history],2)
    actual,ideal=na,ni
    eq(sum(np.trace(x).real for x in actual.values()),1.)
    eq(sum(np.trace(x).real for x in ideal.values()),1.)
    d=sum(trace_norm(actual[k]-ideal[k]) for k in actual)/2
    le(d,min(1.,n*bound(eps)))
    tree.append({"calls":n,"sample_full_record_distance":d,
                 "uniform_upper":min(1.,n*bound(eps))})

# Active record memory counterexample: same total instrument, no D-only closure.
m=np.eye(3,dtype=complex)
memory_ins=[]
for s in [1,2]:
    memory_ins.append([np.kron(I,np.outer(m[:,s],m[:,0].conj()))/np.sqrt(2),
                       np.kron(I,np.outer(m[:,s],m[:,s].conj()))])
effects=[sum(k.conj().T@k for k in ks) for ks in memory_ins]
eq(sum(effects),np.eye(6))
init=np.kron(rho_rand(2),np.outer(m[:,0],m[:,0]))
for s,ks in enumerate(memory_ins):
    post=sum(k@init@k.conj().T for k in ks)
    eq(np.trace(post),.5)
    eq(np.trace(effects[1-s]@post),0.)
    eq(effects[s].reshape(2,3,2,3)[:,0,:,0],I/2)

# Explicit positive direction-family calibration, not generation of its S2 domain.
eps=.01
n=np.array([1.,2.,3.]);n/=np.linalg.norm(n)
en=eps*I+(1-2*eps)*(I+n[0]*X+n[1]*Y+n[2]*Z)/2
em=eps*I+(1-2*eps)*(I-n[0]*X-n[1]*Y-n[2]*Z)/2
eq(np.linalg.norm(en-em,2),1-2*eps)
eq(np.trace(en),1)
le(eps,np.linalg.norm(en-em,2))

# Count-zero preliminary channel-family screen: fixed readouts have rank five.
paulis=[X,Y,Z]
basis=[]
for a,b,c,d,e in np.eye(5):
    basis.append(np.array([[a,c,d],[c,b,e],[d,e,-a-b]]))
jac=np.array([[s[0,0],s[1,1],s[0,1],s[0,2],s[1,2]] for s in basis]).T/2
eq(np.linalg.matrix_rank(jac),5)
s=sum(v*x for v,x in zip([.01,.02,.03,.015,-.02],basis))
js=(np.eye(4)+sum(s[i,j]*np.kron(paulis[i],paulis[j].T)
                  for i in range(3) for j in range(3)))/4
psd(js);eq(np.trace(js),1)
le((1-3*np.linalg.norm(s))/4,np.linalg.eigvalsh(js)[0])
print(json.dumps({"round":1073,"passed":True,"checks":checks,
 "max_equality_residual":max_equality_residual,
 "instrument_cases":rows,"sqrt_order_exact_cases":sqrt_rows,
 "sharp_bias_threshold":threshold_rows,"adaptive_tree":tree,
 "memory_counterexample":{"first_effect":"I/2","apparent_repeat_error":0,
                            "D_only_closed_instrument":False},
 "count_zero_channel_screen_rank":5,
 "finite_budget_example":{"epsilon":.0001,"calls":20,
                         "upper":min(1,20*bound(.0001))},
 "universal_claims_by":"analytic proof, not random samples or SDP",
 "actual_spatial_dimension_derived":False},ensure_ascii=False,indent=2))
