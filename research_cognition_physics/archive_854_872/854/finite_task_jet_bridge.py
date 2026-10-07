"""854: finite-order process compression, with history and source jets.

This ladder diagnostic checks a two-stage common process. It is not a
numerical construction of the continuum QFT or of the old graph limit.
"""
from pathlib import Path
import argparse,itertools,json,math
import numpy as np
HERE=Path(__file__).resolve().parent;TARGET=HERE/'finite_task_jet_bridge_results.json'
D=6
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.]).astype(complex)

def pmul(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            k=tuple(x+y for x,y in zip(ka,kb))
            if k[0]<=D:out[k]=out.get(k,0)+va@vb
    return out

def exp_series(h):
    dim=next(iter(h.values())).shape[0]
    term={(0,0,0):np.eye(dim,dtype=complex)};out=dict(term)
    for n in range(1,D+1):
        term={k:1j*v/n for k,v in pmul(term,h).items()}
        for k,v in term.items():out[k]=out.get(k,0)+v
    return out

def exact_u(h,t,u,v):
    m=sum((t**k[0]*u**k[1]*v**k[2]*a for k,a in h.items()),np.zeros_like(next(iter(h.values()))))
    e,w=np.linalg.eigh(m);return (w*np.exp(1j*e))@w.conj().T

def run():
    ne=12;nr=D+1;dim=2*ne
    q=np.zeros((ne,ne),complex)
    for n in range(ne-1):q[n,n+1]=q[n+1,n]=n+1
    number=np.diag(np.arange(ne,dtype=float));eye=np.eye(ne,dtype=complex)
    h1={(1,0,0):np.kron(X,q)/3,(2,0,0):np.kron(Z,number)/5,
        (2,1,0):np.kron(Z,q)/7,(6,2,0):np.kron(Y,eye)/11}
    h2={(1,0,0):np.kron(Z,q)/4,(2,0,0):np.kron(Y,q@q)/6,
        (2,0,1):np.kron(X,q)/9,(6,1,1):np.kron(Z,eye)/13}
    assert all(np.max(abs(a-a.conj().T))==0 for h in (h1,h2) for a in h.values())
    # q^j |0> reaches level j with coefficient j!: this generates the
    # exact word-closed environment through degree D, not just the input.
    v=np.eye(ne)[:,0].astype(complex);cyclic=[]
    for n in range(D+1):
        assert abs(v[n]-math.factorial(n))<1e-12 and max(abs(v[n+1:]),default=0)==0
        cyclic.append(v.copy());v=q@v
    assert np.linalg.matrix_rank(np.stack(cyclic,axis=1))==D+1
    w0=np.kron(I,np.eye(ne)[:,[0]])
    initial=(I+.2*X+.3*Y+.4*Z)/2
    assert np.linalg.eigvalsh(initial).min()>0
    full1,full2=exp_series(h1),exp_series(h2)
    signs=(-1,1)
    def make_projector(m,size,sign):return np.kron((I+sign*m)/2,np.eye(size))
    def histories(s1,s2,size,w):
        out={}
        for a,b in itertools.product(signs,repeat=2):
            k1=make_projector(Z,size,a);k2=make_projector(X,size,b)
            stage1={k:k1@v for k,v in s1.items()}
            v=pmul(s2,stage1);out[a,b]={k:k2@m@w for k,m in v.items()}
        return out
    full_hist=histories(full1,full2,ne,w0)
    # Preserve the whole code matrix algebra by a tensor-compatible fixed P.
    def compress(n):
        emb=np.kron(I,np.eye(ne)[:,:n]);proj=emb@emb.conj().T
        assert max(np.max(abs(proj@np.kron(a,eye)-np.kron(a,eye)@proj)) for a in (X,Y,Z))==0
        hs=[{k:emb.conj().T@v@emb for k,v in h.items()} for h in (h1,h2)]
        ss=[exp_series(h) for h in hs]
        ww=emb.conj().T@w0
        hist=histories(*ss,n,ww)
        embedded={r:{k:emb@v for k,v in a.items()} for r,a in hist.items()}
        return emb,hs,ss,ww,hist,embedded
    emb,hs,ss,ww,compressed,embedded=compress(nr)
    _,bad_hs,bad_ss,bad_w,bad_hist,bad_embedded=compress(1)
    def amplitude_error(ref,other):
        return max(float(np.max(abs(v-other[r].get(k,np.zeros_like(v))))) for r,p in ref.items() for k,v in p.items())
    amp_error=amplitude_error(full_hist,embedded);bad_amp=amplitude_error(full_hist,bad_embedded)
    assert amp_error<2e-12 and bad_amp>.01
    def gram_jets(hist):
        out={}
        for outcome,p in hist.items():
            g={}
            for ka,a in p.items():
                for kb,b in p.items():
                    n=ka[0]+kb[0]
                    if n>D:continue
                    k=(n,ka[1],ka[2],kb[1],kb[2])
                    g[k]=g.get(k,0)+a.conj().T@b
            out[outcome]=g
        return out
    gf,gc,gb=map(gram_jets,(full_hist,compressed,bad_hist))
    def gram_error(a,b,predicate=lambda k:True):
        return max(float(np.max(abs(v-b[r].get(k,np.zeros_like(v))))) for r,g in a.items() for k,v in g.items() if predicate(k))
    gram_err=gram_error(gf,gc);source_err=gram_error(gf,gc,lambda k:sum(k[1:])>0)
    bad_source=gram_error(gf,gb,lambda k:sum(k[1:])>0)
    assert max(gram_err,source_err)<3e-12 and bad_source>.001
    mixed=[]
    for r,g in gf.items():
        for k,m in g.items():
            if sum(k[1:])==2 and sum(k[1:3]) and sum(k[3:]) and np.max(abs(m))>1e-6:
                mixed.append(dict(outcome=list(r),degree_and_sources=list(k),input_expectation_real=float(np.trace(initial@m).real),input_expectation_imag=float(np.trace(initial@m).imag)))
    assert mixed
    # Physical k-source extraction after z=t^3 u carries t^{-k}.
    # For N=4, K=2, the t^6 u^2 term cannot be discarded at order four.
    contact=h1[6,2,0]
    contact_expectation=complex(np.trace(initial@(w0.conj().T@(2j*contact)@w0)))
    assert abs(contact_expectation)>0
    # Discarding the same environment after the first actual record is a
    # different process. Compare all source-free history coefficients.
    sourcefree1={k[0]:v for k,v in full1.items() if k[1:]==(0,0)}
    sourcefree2={k[0]:v for k,v in full2.items() if k[1:]==(0,0)}
    def partial_trace(a):return np.einsum('aibi->ab',a.reshape(2,ne,2,ne))
    vac=np.zeros((ne,ne),complex);vac[0,0]=1
    reset_differences=[]
    for a,b in itertools.product(signs,repeat=2):
        k1=make_projector(Z,ne,a);k2=make_projector(X,ne,b)
        rhos={}
        for n in range(D+1):
            dens=np.zeros((dim,dim),complex)
            for l,ul in sourcefree1.items():
                if n-l in sourcefree1:
                    ur=sourcefree1[n-l];dens+=k1@ul@w0@initial@w0.conj().T@ur.conj().T@k1
            rhos[n]=np.kron(partial_trace(dens),vac)
        for n in range(D+1):
            actual=complex(np.trace(initial@gf[a,b].get((n,0,0,0,0),np.zeros((2,2)))))
            reset=0j
            for i,ui in sourcefree2.items():
                for j,uj in sourcefree2.items():
                    k=n-i-j
                    if k>=0:reset+=np.trace(k2@ui@rhos[k]@uj.conj().T@k2)
            if abs(actual-reset)>1e-10:reset_differences.append(dict(outcome=[a,b],degree=n,difference_real=float((actual-reset).real),difference_imag=float((actual-reset).imag)))
    assert reset_differences
    first_reset=min(x['degree'] for x in reset_differences)
    # The reduced process is exactly unitary and normalized at finite values;
    # this is a property of Q_eff, not convergence of the continuum series.
    finite_rows=[]
    for t in (.15,.35,.7):
        u,v=.4,-.3;us=[exact_u(h,t,u,v) for h in hs]
        probabilities=[]
        for a,b in itertools.product(signs,repeat=2):
            vv=make_projector(X,nr,b)@us[1]@make_projector(Z,nr,a)@us[0]@ww
            probabilities.append(float(np.trace(vv@initial@vv.conj().T).real))
        unit=max(float(np.max(abs(v.conj().T@v-np.eye(2*nr)))) for v in us)
        assert min(probabilities)>0 and abs(sum(probabilities)-1)<2e-13 and unit<2e-13
        finite_rows.append(dict(epsilon=t,sources=[u,v],probabilities=probabilities,normalization_error=abs(sum(probabilities)-1),unitarity_error=unit))
    return dict(round=854,all_checks_passed=True,fresh_test_groups=1,
        diagnostic='two-stage noncommuting ladder process with common memory and bilateral source jets; not the original continuum or graph Hamiltonian',
        max_epsilon_degree=D,reference_dimension=dim,word_closed_dimension=2*nr,input_only_dimension=2,
        amplitude_coefficient_max_error=amp_error,history_gram_coefficient_max_error=gram_err,
        bilateral_source_coefficient_max_error=source_err,input_only_amplitude_error=bad_amp,input_only_source_error=bad_source,
        nonzero_mixed_source_examples=mixed[:4],
        source_filter=dict(target_epsilon_order=4,max_insertions=2,required_series_order=6,normalized_contact_lost_at_order_four=[contact_expectation.real,contact_expectation.imag]),
        reset_first_wrong_epsilon_degree=first_reset,reset_examples=reset_differences[:4],finite_reduced_process_rows=finite_rows,
        all_input_density_matrices_matched_through_declared_order=True,
        original_continuum_remainder_bound_or_old_graph_mapping_claimed=False,
        projected_generators_local_or_universal_Hamiltonian_claimed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
