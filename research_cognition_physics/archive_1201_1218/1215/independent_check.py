"""Independent finite balance-window checks. No author/old program imports."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
OUT=HERE/'independent_results.json'
SOURCES=['archive_370_428/research_note_421.md',
         'archive_1201_/1209/proof.md','archive_1201_/1214/proof.md']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def trdist(a,b):return float(np.sum(np.linalg.svd(a-b,compute_uv=False))/2)
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)):return round(float(x),12)
    if isinstance(x,np.integer):return int(x)
    return x

def shells(energies,L,select):
    """Complete finite unitary, including partial lowest/highest shells."""
    energies=list(energies);d=len(energies);W=np.eye(d*(L+1),dtype=complex)
    Q=np.array([e+n for e in energies for n in range(L+1)])
    for E in range(max(energies)+L+1):
        ids=[(i,E-e) for i,e in enumerate(energies) if 0<=E-e<=L]
        if len(ids)!=d:continue
        ind=[i*(L+1)+n for i,n in ids]
        W[np.ix_(ind,ind)]=select(E)
    assert norm(W.conj().T@W-np.eye(len(W)))<2e-12
    assert np.max(np.abs((Q[:,None]-Q[None,:])*W))<2e-12
    return W,Q

def kraus(W,b,d):
    V=(W@np.kron(np.eye(d),b[:,None])).reshape(d,len(b),d)
    return V.transpose(1,0,2)
def phi(K,x):return np.einsum('aij,jk,alk->il',K,x,K.conj())
def residuals(K):
    d=K.shape[1];I=np.eye(d)
    return dict(TP=norm(np.einsum('aji,ajk->ik',K.conj(),K)-I),
                unital=norm(phi(K,I)-I))
def entropy(rho):
    vals=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    assert min(vals)>-2e-12
    return float(-sum(x*np.log(x) for x in vals if x>1e-14))

def run():
    groups={};H=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
    d=3;L=12;e=np.array([0,1,3]);omega=np.exp(2j*np.pi/3)
    U=np.array([[omega**(i*j) for j in range(d)] for i in range(d)])/np.sqrt(d)
    W,Q=shells(e,L,lambda E:U)
    b=np.zeros(L+1,complex);b[5:8]=[1,1j,-1];b/=np.sqrt(3)
    K=kraus(W,b,d);res=residuals(K)
    assert res['TP']<2e-12 and res['unital']<2e-12
    energy_rows=[]
    for j in range(d):
        out=(W@np.kron(np.eye(d)[:,j],b)).reshape(d,L+1)
        deS=float(sum(e*np.sum(abs(out)**2,axis=1))-e[j])
        deC=float(sum(np.arange(L+1)*np.sum(abs(out)**2,axis=0))-
                  sum(np.arange(L+1)*abs(b)**2))
        assert abs(deS+deC)<2e-12
        energy_rows.append(dict(input=j,material_delta=deS,balance_delta=deC))
    groups['finite_positive_ladder_complete_boundaries']=dict(
        dimension=d*(L+1),energies=e.tolist(),support=[5,6,7],**res,
        energy_rows=energy_rows,full_boundary_unitarity=True)

    # Joint SR output, independently using full W tensor I_R and Kraus action.
    r=2;v=np.array([1,2j,-1,3,1j,2],complex);v/=np.linalg.norm(v)
    rho=np.outer(v,v.conj())
    initial=np.einsum('ir,n->inr',v.reshape(d,r),b).reshape(-1)
    final=(np.kron(W,np.eye(r))@initial).reshape(d,L+1,r)
    direct=np.einsum('inr,kns->irks',final,final.conj()).reshape(d*r,d*r)
    expected=sum(np.kron(k,np.eye(r))@rho@np.kron(k.conj().T,np.eye(r)) for k in K)
    joint_err=norm(direct-expected);assert joint_err<2e-12
    groups['unknown_reference_preserved']=dict(reference_dimension=r,joint_output_residual=joint_err)

    # Deletion with identical transition probabilities but unequal shell phases.
    L2=7;m=3;b2=np.zeros(L2+1);b2[m:m+2]=1/np.sqrt(2)
    W0,_=shells([0,1],L2,lambda E:H)
    W1,_=shells([0,1],L2,lambda E:-H if E==m+2 else H)
    a=phi(kraus(W0,b2,2),np.eye(2)/2)
    z=phi(kraus(W1,b2,2),np.eye(2)/2)
    exact=np.array([[.5,.25],[.25,.5]])
    assert np.max(abs(abs(W0)**2-abs(W1)**2))<1e-15
    assert norm(a-np.eye(2)/2)<2e-12 and norm(z-exact)<2e-12
    # Coherence: sigma[m,m+1] * (H00 H10 - H01 H11), divided by d.
    coherence=F(1,2)*(F(1,2)+F(1,2))/2
    assert coherence==F(1,4)
    groups['probability_only_contract_fails']=dict(
        identical_full_transition_probabilities=True,
        normalized_output_exact=[['1/2','1/4'],['1/4','1/2']],
        normalized_unital_defect_exact='1/4',
        output_eigenvalues_exact=['1/4','3/4'],diagnostic_error=norm(z-exact))

    # Energy-only address: a finite pure balance absorbs unknown S information.
    X=np.array([[0,1],[1,0]],complex);bb=np.eye(L2+1)[:,m]
    Wr,_=shells([0,1],L2,lambda E:X if E==m+1 else np.eye(2))
    Kr=kraus(Wr,bb,2);reset=phi(Kr,np.eye(2)/2)
    assert norm(reset-np.diag([1,0]))<1e-14
    bell=np.eye(2).reshape(-1)/np.sqrt(2)
    vi=np.einsum('ir,n->inr',bell.reshape(2,2),bb).reshape(-1)
    vo=(np.kron(Wr,np.eye(2))@vi).reshape(2,L2+1,2)
    assert norm(vo[1])<1e-14
    embedding=np.zeros((L2+1,2));embedding[m,0]=1/np.sqrt(2);embedding[m+1,1]=1/np.sqrt(2)
    assert norm(vo[0]-embedding)<1e-14
    groups['balance_address_reset']=dict(normalized_unital_defect_exact='1/2',
        output_material_ground=True,unknown_R_transferred_to_balance=True,
        reference_residual=norm(vo[0]-embedding),additional_record_dimension=1)

    # Terminal QND is harmless for unitality; feedback can spend the pure record.
    # Logical X=(S,K), e=(0,0,1,1), K carries zero energy.
    e4=[0,0,1,1];copy=np.zeros((4,4));feed=np.zeros((4,4))
    for s in range(2):
        for k in range(2):
            copy[2*s+(k^s),2*s+k]=1
            feed[2*(s^k)+k,2*s+k]=1
    Wf,_=shells(e4,L2,lambda E:feed)
    Wc=np.kron(copy,np.eye(L2+1));whole=Wf@Wc
    Kwhole=kraus(whole,bb,4);whole_res=residuals(Kwhole)
    assert whole_res['unital']<1e-14 and whole_res['TP']<1e-14
    rhoSK=np.kron(np.eye(2)/2,np.diag([1,0]))
    out=phi(Kwhole,rhoSK)
    target=np.trace(out.reshape(2,2,2,2),axis1=1,axis2=3)
    assert norm(target-np.diag([1,0]))<1e-14
    copy_out=copy@rhoSK@copy.T
    assert norm(np.trace(copy_out.reshape(2,2,2,2),axis1=1,axis2=3)-np.eye(2)/2)<1e-14
    groups['QND_record_then_feedback_scope']=dict(
        complete_SK_channel_unital=True,whole_material_residuals=whole_res,
        terminal_QND_preserves_material_maxmixed=True,
        feedback_with_initial_pure_K_resets_S=True,
        full_initial_SK_min_eigenvalue=0,
        statement='The same balance-independent whole SK unitary can consume K purity; S alone is not the material domain.')

    # Genuine two-block first-hit permutation, lift and full terminal record.
    labels=[(i,j) for i in range(3) for j in range(3)]
    T=[]
    for s in labels:
        t=list(s)
        for p,x in enumerate(s):
            if x!=1:t[p]=2-x;break
        T.append(labels.index(tuple(t)))
    P=np.zeros((9,9));P[T,np.arange(9)]=1
    energies=[sum(s) for s in labels];L3=18
    W3,_=shells(energies,L3,lambda E:P)
    b3=np.zeros(L3+1);b3[6:13]=1/np.sqrt(7)
    V=(W3@np.kron(np.eye(9),b3[:,None])).reshape(9,L3+1,9)
    event=[]
    for s in labels:
        event.append(next((int(x==0) for x in s if x!=1),0))
    Vrec=np.zeros((9,L3+1,2,9),complex)
    for i in range(9):Vrec[i,:,event[i],:]=V[i]
    K3=Vrec.reshape(9,2*(L3+1),9).transpose(1,0,2)
    rr=residuals(K3);assert rr['unital']<3e-12 and rr['TP']<3e-12
    probs=[F(1,2),F(3,10),F(1,5)]
    product=[probs[i]*probs[j] for i,j in labels]
    before=sum(p for p,q in zip(product,event) if q)
    after=sum(product[j] for j in range(9) if event[T[j]])
    assert before==F(13,20) and after==F(13,50)
    out=phi(K3,np.diag([float(p) for p in product]))
    expected=np.diag(np.array([float(p) for p in product])[T])
    assert norm(out-expected)<3e-12
    groups['two_block_lift_with_terminal_record']=dict(material_dimension=9,
        first_hit_before_exact=str(before),first_hit_after_exact=str(after),
        event_drift_exact=str(before-after),all_boundary_shells_retained=True,
        full_record_included=True,**rr,
        claim='Small independent membership check; the all-size 1209 theorem is reused analytically.')

    # One explicit positive autonomous endpoint for the reset deletion.
    # All Q eigenvalues integers. W_r is a Hermitian involution.
    Qr=np.diag(np.array([e+n for e in [0,1] for n in range(L2+1)],float))
    G=(np.eye(len(Wr))-Wr)/4
    parent=Qr+G
    vals,vec=np.linalg.eigh(parent)
    endpoint=(vec*np.exp(-2j*np.pi*vals))@vec.conj().T
    assert min(vals)>-1e-12 and norm(endpoint-Wr)<2e-12
    assert norm(Qr@G-G@Qr)<1e-14
    groups['autonomous_endpoint_scope']=dict(positive_parent=True,time='2*pi',
        endpoint_residual=norm(endpoint-Wr),constant_interaction_retained=True,
        physical_clock_or_release_certified=False)
    return clean(dict(schema='1215-independent-v1',status='PASS',all_checks_passed=True,
        check_groups=len(groups),groups=groups,
        source_sha256={p:sha(BASE/p) for p in SOURCES},
        independent_code_sha256=sha(Path(__file__)),
        not_certified=['No new cognitive axiom adopted.',
         'No general unital equals mixed-unitary converse.',
         'No arbitrary initial/boundary battery or permanent cold output.',
         'No native QQ controller implementation or four-law completion.']))

def main():
    result=run()
    if '--write' in sys.argv:
        if OUT.exists():raise SystemExit('Refuse to overwrite existing results')
        OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(OUT.read_text(encoding='utf-8'))
        assert saved==result,'Saved result differs'
    print(json.dumps(dict(status='PASS',all_checks_passed=True,results_match=True,
                         check_groups=result['check_groups'])))

if __name__=='__main__':main()
