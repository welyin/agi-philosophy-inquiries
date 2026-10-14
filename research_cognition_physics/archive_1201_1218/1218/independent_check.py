"""Independent shared-balance composition checks; no author/old code imports.

Exact deletion certificates use integer permutation columns. Floating matrices
check general coherent gates and a passive reference, not a uniform theorem.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
OUT=HERE/'independent_results.json'

def require(value,label):
    if not bool(value):raise AssertionError(label)

def norm(a):return float(np.linalg.norm(a,2))
def distance(a,b):return float(np.linalg.svd(a-b,compute_uv=False).sum()/2)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def local_balance(U,N,which):
    """Basis order A,B,C. Every full two-dimensional energy shell is kept."""
    size=4*(N+1)
    W=np.eye(size,dtype=complex)
    def index(a,b,c):return (2*a+b)*(N+1)+c
    for other in (0,1):
        for e in range(1,N+1):
            ids=([index(0,other,e),index(1,other,e-1)] if which==0
                 else [index(other,0,e),index(other,1,e-1)])
            W[np.ix_(ids,ids)]=U
    total=np.array([a+b+c for a in (0,1) for b in (0,1) for c in range(N+1)])
    require(norm(W.conj().T@W-np.eye(size))<5e-12,'full unitary')
    require(np.max(abs((total[:,None]-total[None,:])*W))<1e-12,'full total energy')
    return W,total

def isometry(W,b):return W@np.kron(np.eye(4),np.asarray(b)[:,None])
def material(V,rho,N):
    out=(V@rho@V.conj().T).reshape(4,N+1,4,N+1)
    return np.trace(out,axis1=1,axis2=3)
def balance(V,rho,N):
    return np.trace((V@rho@V.conj().T).reshape(4,N+1,4,N+1),axis1=0,axis2=2)
def one_marginal(rho,which):
    shaped=rho.reshape(2,2,2,2)
    return np.trace(shaped,axis1=1-which,axis2=3-which)

def lifted_columns(U,N,b):
    energies=[0,1,1,2]
    V=np.zeros((4*(N+1),4),complex)
    for j,e in enumerate(energies):
        for c,weight in enumerate(b):
            if abs(weight)<1e-15:continue
            for i,f in enumerate(energies):
                if abs(U[i,j])<1e-15:continue
                target=c+e-f
                require(0<=target<=N,'ideal path terminal support')
                V[i*(N+1)+target,j]+=weight*U[i,j]
    return V

def flip_permutation(N,which):
    ans=[]
    for a in (0,1):
        for b in (0,1):
            for c in range(N+1):
                bits=[a,b];s=bits[which]
                if (s==0 and c>0) or (s==1 and c<N):
                    bits[which]=1-s;c2=c+(1 if s==1 else -1)
                else:c2=c
                ans.append((2*bits[0]+bits[1])*(N+1)+c2)
    require(sorted(ans)==list(range(4*(N+1))),'complete permutation')
    return ans

def exact_maxmixed_output(N,support,order=(0,1)):
    # All amplitudes +1/sqrt(len(support)); initial material exactly I4/4.
    perms=[flip_permutation(N,k) for k in order]
    V=np.zeros((4*(N+1),4),dtype=object)
    for j in range(4):
        for c in support:
            out=j*(N+1)+c
            for p in perms:out=p[out]
            V[out,j]+=1
    gram=V@V.T
    den=4*len(support)
    rho=[[F(sum(int(gram[i*(N+1)+c,j*(N+1)+c]) for c in range(N+1)),den)
          for j in range(4)] for i in range(4)]
    rC=[F(sum(int(gram[j*(N+1)+c,j*(N+1)+c]) for j in range(4)),den)
        for c in range(N+1)]
    return rho,rC

def entropy(r):
    e=np.linalg.eigvalsh((r+r.conj().T)/2)
    require(min(e)>-5e-12,'state positivity')
    return -sum(float(x)*math.log(float(x)) for x in e if x>1e-14)

def compute():
    groups={};X=np.array([[0,1],[1,0]],complex)
    deletions=[]
    for N in (4,6):
        support=[1,N-1]
        for which in (0,1):
            r,_=exact_maxmixed_output(N,support,(which,))
            require(r==[[F(int(i==j),4) for j in range(4)] for i in range(4)],'single preparation calibration')
        rho,rC=exact_maxmixed_output(N,support)
        target=[F(1,8),F(3,8),F(3,8),F(1,8)]
        require(rho==[[target[i] if i==j else F(0) for j in range(4)] for i in range(4)],'exact joint output')
        defect=sum(abs(target[i]-F(1,4)) for i in range(4))/2
        EA=target[2]+target[3];EB=target[1]+target[3]
        EC=sum(c*p for c,p in enumerate(rC))
        require(defect==F(1,4) and EA==EB==F(1,2) and EC==F(N,2),'exact averages and defect')
        W0,_=local_balance(X,N,0);W1,_=local_balance(X,N,1)
        b=np.zeros(N+1);b[support]=1/np.sqrt(2)
        V=isometry(W1@W0,b)
        out=material(V,np.eye(4)/4,N);bat=balance(V,np.eye(4)/4,N)
        require(norm(out-np.diag([float(x) for x in target]))<1e-12,'exact vs matrix')
        require(norm(one_marginal(out,0)-np.eye(2)/2)<1e-12 and norm(one_marginal(out,1)-np.eye(2)/2)<1e-12,'both marginal checks miss defect')
        mid=balance(isometry(W0,b),np.eye(4)/4,N)
        edge=float(mid[0,0].real+mid[N,N].real)
        require(abs(edge-.5)<1e-12,'intermediate boundary weight')
        totalout=V@(np.eye(4)/4)@V.conj().T
        mutual=entropy(out)+entropy(bat)-entropy(totalout)
        require(mutual>0,'carrier correlations retained')
        safe=np.zeros(N+1);safe[N//2]=1
        safeout=material(isometry(W1@W0,safe),np.eye(4)/4,N)
        require(norm(safeout-np.eye(4)/4)<1e-12,'same mean safe preparation')
        deletions.append(dict(N=N,initial_support=support,
            joint_diagonal_exact=[str(x) for x in target],normalized_unital_defect_exact=str(defect),
            mean_A_exact=str(EA),mean_B_exact=str(EB),mean_C_exact=str(EC),
            both_single_original_preparations_unital=True,both_output_marginals_maxmixed=True,
            intermediate_boundary_weight=edge,joint_entropy_drop=math.log(4)-entropy(out),
            final_carrier_entropy=entropy(bat),final_material_carrier_mutual_information=mutual))
    groups['same_preparation_individual_checks_do_not_compose']=dict(passed=True,examples=deletions)

    # Safe shared use is a joint lift even for distinct coherent local gates.
    N=6;b=np.zeros(N+1,complex);b[2:5]=[1,1j,-1];b/=np.sqrt(3)
    H=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
    R=np.array([[np.sqrt(3)/2,-1j/2],[-1j/2,np.sqrt(3)/2]],complex)
    WA,total=local_balance(H,N,0);WB,_=local_balance(R,N,1)
    U=np.kron(np.eye(2),R)@np.kron(H,np.eye(2))
    V=isometry(WB@WA,b);Videal=lifted_columns(U,N,b)
    e=norm(V-Videal)
    require(e<2e-12 and norm(V.conj().T@V-np.eye(4))<2e-12,'safe complete isometry identity')
    require(norm(material(V,np.eye(4),N)-np.eye(4))<2e-12,'joint unital safe')
    groups['safe_full_input_lift_composition']=dict(passed=True,N=N,initial_support=[2,3,4],
        isometry_identity_error=e,all_intermediate_columns_inside_local_safe_window=True,
        endpoint_carrier_not_assumed_independent_or_returned=True)

    # Passive R check on an arbitrary entangled ABR input, no per-step reset.
    dR=3;z=np.arange(1,13)+1j*np.array([0,1,-1,2,0,1,-2,1,3,0,-1,2]);z=z/np.linalg.norm(z)
    input_joint=np.einsum('jr,c->jcr',z.reshape(4,dR),b).reshape(-1)
    full=np.kron(WB@WA,np.eye(dR))@input_joint
    full=full.reshape(4,N+1,dR)
    direct=np.einsum('icr,jcs->irjs',full,full.conj()).reshape(4*dR,4*dR)
    K=V.reshape(4,N+1,4).transpose(1,0,2)
    zz=np.outer(z,z.conj());expected=sum(np.kron(k,np.eye(dR))@zz@np.kron(k.conj().T,np.eye(dR)) for k in K)
    require(norm(direct-expected)<2e-12,'passive unknown reference')
    groups['passive_reference_and_no_midcourse_repreparation']=dict(passed=True,
        reference_dimension=dR,full_joint_output_error=norm(direct-expected))

    # A shared phase reference generates correlations even in the safe class.
    b0=np.zeros(N+1);b0[3]=1
    HA,_=local_balance(H,N,0);HB,_=local_balance(H,N,1)
    shared=material(isometry(HB@HA,b0),np.diag([1,0,0,0]),N)
    exact=np.eye(4)/4;exact[1,2]=exact[2,1]=.25
    require(norm(shared-exact)<2e-12,'shared phase correlation')
    require(abs(distance(shared,np.eye(4)/4)-.25)<2e-12,'fresh independent replacement changes process')
    groups['safe_joint_unital_does_not_mean_fresh_product_channel']=dict(passed=True,
        shared_coherence_01_10_exact='1/4',distance_from_fresh_product_exact='1/4',
        non_CP_divisibility_claimed=False)

    # Actual intermediate QND copy, never used by the second gate.
    first=isometry(HA,b0).reshape(4,N+1,4)
    rec=np.zeros((4,N+1,2,4),complex)
    for s in range(4):rec[s,:,s//2,:]=first[s]
    final=(np.kron(HB,np.eye(2))@rec.reshape(4*(N+1)*2,4)).reshape(4,N+1,2,4)
    KK=final.reshape(4,2*(N+1),4).transpose(1,0,2)
    unital=sum(k@k.conj().T for k in KK)
    tp=sum(k.conj().T@k for k in KK)
    require(norm(unital-np.eye(4))<2e-12 and norm(tp-np.eye(4))<2e-12,'no-feedback QND composition')
    rho=np.diag([.4,.3,.2,.1]);joint=final.reshape(4*(N+1)*2,4)@rho@final.reshape(4*(N+1)*2,4).conj().T
    rK=np.trace(np.trace(joint.reshape(4,N+1,2,4,N+1,2),axis1=0,axis2=3),axis1=0,axis2=2)
    require(entropy(rK)>.1,'QND blank is consumed')
    groups['midcourse_energy_QND_without_feedback']=dict(passed=True,
        full_material_unital_error=norm(unital-np.eye(4)),TP_error=norm(tp-np.eye(4)),
        final_record_entropy=entropy(rK),record_blank_reuse_claimed=False)

    # Error telescoping must act on the full ideal reachable input subspace.
    theta=.01;delta=.007
    D1=np.diag([np.exp(1j*theta),np.exp(-1j*theta)])
    D2=np.diag([np.exp(-1j*delta),np.exp(1j*delta)])
    AA,_=local_balance(D1@H,N,0);BB,_=local_balance(D2@R,N,1)
    initial=np.kron(np.eye(4),b[:,None])
    a1=norm((AA-WA)@initial);a2=norm((BB-WB)@WA@initial)
    end=norm((BB@AA-WB@WA)@initial)
    require(end<=a1+a2+1e-12,'telescoping complete input bound')
    groups['full_reachable_subspace_error_budget']=dict(passed=True,
        first_error=a1,second_ideal_reachable_error=a2,endpoint_error=end,
        summed_bound=a1+a2,per_column_or_reinitialized_bound_used=False)

    sources=['archive_1201_/1215/proof.md','archive_1201_/1216/proof.md']
    result=dict(round=1218,status='INDEPENDENT_PASS_SCOPED_NOT_ACCEPTANCE',
        all_checks_passed=True,check_groups=len(groups),groups=groups,
        source_sha256={p:sha(BASE/p) for p in sources},independent_code_sha256=sha(Path(__file__)),
        scope=dict(whole_material_domain=True,work_role_adopted=False,
            actual_native_QQ_controller=False,new_clock=False,all_time_safe_reuse=False,
            no_universal_thermodynamic_or_third_law_failure_claim=True))
    return result

def compare(a,b,path=()):
    require(type(a)==type(b),'type '+str(path))
    if isinstance(a,dict):
        require(a.keys()==b.keys(),'keys '+str(path))
        for k in a:compare(a[k],b[k],path+(k,))
    elif isinstance(a,list):
        require(len(a)==len(b),'length '+str(path))
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+(i,))
    elif isinstance(a,float):
        require(math.isfinite(a) and math.isfinite(b) and abs(a-b)<1e-11,'diagnostic '+str(path))
    else:require(a==b,'exact '+str(path))

def main():
    if '--write' in sys.argv and OUT.exists():raise SystemExit('Refuse to overwrite results')
    out=compute()
    if '--write' in sys.argv:
        with OUT.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    compare(out,json.loads(OUT.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',all_checks_passed=True,results_match=True,check_groups=out['check_groups'])))

if __name__=='__main__':main()
