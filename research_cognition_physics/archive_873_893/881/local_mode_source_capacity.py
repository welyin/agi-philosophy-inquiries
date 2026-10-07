"""881: correlated Gaussian polynomial/source tasks under local capacity cutoffs.
The infinite reference moments are evaluated by Wick's theorem, never by a
large finite reference cutoff. This calibrates the capacity theorem; it is
not a numerical solution of the original curved mixed quantum field theory.
"""
from pathlib import Path
from functools import lru_cache
from fractions import Fraction
from math import comb
import argparse,json
import numpy as np
TARGET=Path(__file__).with_name('local_mode_source_capacity_results.json')
Q=Fraction(1,4)
I=np.eye(4,dtype=complex)
z=np.diag([1.,-1.]); c=np.array([[0.,1.],[0.,0.]])
cA=np.kron(c,np.eye(2)); cB=np.kron(z,c)
nA=cA.conj().T@cA; nB=cB.conj().T@cB
ZA=I-2*nA; ZB=I-2*nB
hop=cA.conj().T@cB+cB.conj().T@cA
def poly(*terms):return list(terms)
def mul(a,b):
    return [(ca*cb,wa+wb,ma@mb) for ca,wa,ma in a for cb,wb,mb in b]
def adj(a):return [(complex(c).conjugate(),w[::-1],m.conj().T) for c,w,m in a]
ONE=poly((1,(),I))
A=poly((1,(0,),ZA))
B=poly((1,(3,),ZB))
G=poly((1,(1,1),ZA),(Fraction(1,3),(2,),ZB))
C=poly((1,(0,2),hop))
WORDS=[ONE,A,B,G,C,mul(A,G),mul(G,A),mul(B,G),mul(G,B),
       mul(G,C),mul(C,G),mul(poly((1,(),nA)),G),
       mul(poly((1,(),nB)),A),mul(poly((1,(),hop)),C)]
LABELS=['1','A','B','G','C','AG','GA','BG','GB','GC','CG','nA G','nB A','hop C']
v=float((1+Q)/(2*(1-Q))); correlation=float(Q**Fraction(1,2)/(1-Q))
COV=np.diag([v]*4).astype(complex)
COV[0,1]=COV[2,3]=.5j; COV[1,0]=COV[3,2]=-.5j
COV[0,2]=COV[2,0]=correlation; COV[1,3]=COV[3,1]=-correlation
@lru_cache(None)
def wick(w):
    if not w:return 1.+0j
    if len(w)%2:return 0j
    return sum(COV[w[0],w[j]]*wick(w[1:j]+w[j+1:]) for j in range(1,len(w)))
def expect(p):
    return sum((complex(c)*wick(w)*m for c,w,m in p),np.zeros_like(I))
def infinite_gram():
    return np.block([[expect(mul(adj(a),b)) for b in WORDS] for a in WORDS])
def field(vec,k):
    axis=k//2
    out=np.zeros_like(vec); low=[slice(None)]*vec.ndim; high=low.copy()
    low[axis]=slice(0,-1); high[axis]=slice(1,None)
    shape=[1]*vec.ndim; shape[axis]=vec.shape[axis]-1
    weights=np.sqrt(np.arange(1,vec.shape[axis])/2).reshape(shape)
    if k%2==0:
        out[tuple(low)]+=weights*vec[tuple(high)]
        out[tuple(high)]+=weights*vec[tuple(low)]
    else:
        out[tuple(low)]+=-1j*weights*vec[tuple(high)]
        out[tuple(high)]+=1j*weights*vec[tuple(low)]
    return out
def apply(p,vec):
    """P p(q,p) P: pad inside each polynomial to preserve boundary contacts."""
    r=vec.shape[0]-1; d=max((len(w) for _,w,_ in p),default=0)
    padded=np.zeros((r+d+1,r+d+1)+vec.shape[2:],complex)
    padded[:r+1,:r+1]=vec
    out=np.zeros_like(vec)
    for coeff,w,m in p:
        term=padded
        for k in reversed(w):term=field(term,k)
        term=np.einsum('ij,abjk->abik',m,term,optimize=True)
        out+=complex(coeff)*term[:r+1,:r+1]
    return out
# Products must insert the local cutoff between coefficient letters, not only
# compress the already multiplied polynomial at the very end.
SEQUENCES=[(),(A,),(B,),(G,),(C,),(A,G),(G,A),(B,G),(G,B),
           (G,C),(C,G),(poly((1,(),nA)),G),
           (poly((1,(),nB)),A),(poly((1,(),hop)),C)]
def columns(r,vacuum=False):
    psi=np.zeros((r+1,r+1),complex)
    if vacuum:psi[0,0]=float(Q**(r+1))**.5
    else:
        for n in range(r+1):psi[n,n]=np.sqrt(1-float(Q))*float(Q)**(n/2)
    base=psi[:,:,None,None]*I[None,None,:,:]
    cols=[]
    for seq in SEQUENCES:
        out=base
        for p in reversed(seq):out=apply(p,out)
        cols.append(out.reshape(-1,4))
    return np.concatenate(cols,axis=1)
def geom_moments(k):
    s=[1/(1-Q)]
    for d in range(1,k+1):
        s.append(Q/(1-Q)*sum(Fraction(comb(d,j))*s[j] for j in range(d)))
    return s
def moment_tail(start,power):
    s=geom_moments(power)
    return (1-Q)*Q**start*sum(Fraction(comb(power,j))*(1+2*start)**(power-j)*2**j*s[j] for j in range(power+1))
def weighted_trace(r,s=3):
    # Exact infinite three-dimensional support representation after weighting.
    eps=float(Q**(r+1)); total=float(moment_tail(0,2*s)); tail=float(moment_tail(r+1,2*s))
    a=np.sqrt(1-float(Q)); b=np.sqrt(max(0,total-tail-a*a)); g=np.sqrt(tail)
    delta=np.array([[eps,0,-a*g],[0,0,-b*g],[-a*g,-b*g,-tail]])
    norm=float(np.abs(np.linalg.eigvalsh(delta)).sum())
    bound=2*np.sqrt(total*tail)+tail
    assert norm<=bound+1e-10
    return dict(weight_power=s,total_moment=total,tail_moment=tail,
                exact_weighted_trace_norm=norm,proved_upper_bound=float(bound))
def matrix(p,r):
    d=4*(r+1)**2
    return apply(p,np.eye(d,dtype=complex).reshape(r+1,r+1,4,d)).reshape(d,d)
def finite_process():
    r=3; aa,bb,gg,cc=[matrix(p,r) for p in (A,B,G,C)]
    herm=max(float(np.max(np.abs(x-x.conj().T))) for x in (aa,bb,gg,cc))
    generator=.21*aa+.21**2*bb+.21**3*.31*gg+.21**6*.31**2*cc
    eig,u=np.linalg.eigh(generator); unitary=(u*np.exp(1j*eig))@u.conj().T
    # Same state and full two-CAR record algebra; finite reference is correlated.
    boson=np.zeros((r+1)**2,complex)
    for n in range(r+1):boson[n*(r+2)]=np.sqrt(1-float(Q))*float(Q)**(n/2)
    rhoB=np.outer(boson,boson.conj());rhoB[0,0]+=float(Q**(r+1))
    record=np.array([0,1,1j,0],complex)/np.sqrt(2)
    rho=np.kron(rhoB,np.outer(record,record.conj()))
    out=unitary@rho@unitary.conj().T
    probs=[float(np.trace(out[j::4,j::4]).real) for j in range(4)]
    return dict(capacity=r,dimension=len(generator),Hermitian_defect=herm,
                local_even_commutator=float(np.max(np.abs(aa@bb-bb@aa))),
                unitarity_defect=float(np.max(np.abs(unitary.conj().T@unitary-np.eye(len(unitary))))),
                record_probabilities=probs,probability_sum_defect=abs(sum(probs)-1),
                parity_defect=float(np.max(np.abs(generator@np.kron(np.eye((r+1)**2),ZA@ZB)
                                                 -np.kron(np.eye((r+1)**2),ZA@ZB)@generator))))
def run():
    gram=infinite_gram()
    assert np.max(np.abs(gram-gram.conj().T))<1e-12
    assert np.min(np.linalg.eigvalsh(gram))>-1e-10
    exact_hessian=2j*expect(C)-expect(mul(G,G))
    rows=[]
    for r in (2,4,8,12,20,28):
        v0=columns(r);v1=columns(r,True)
        approx=v0.conj().T@v0+v1.conj().T@v1
        hess=2j*approx[0:4,16:20]-approx[12:16,12:16]
        err=float(np.linalg.norm(approx-gram,2))
        rows.append(dict(capacity=r,total_finite_dimension=4*(r+1)**2,
                         bosonic_tail_probability=float(Q**(r+1)),
                         full_input_word_Gram_error=err,
                         mixed_source_Hessian_error=float(np.linalg.norm(hess-exact_hessian,2)),
                         weighted_state=weighted_trace(r)))
    assert rows[-1]['full_input_word_Gram_error']<1e-9
    assert rows[-1]['mixed_source_Hessian_error']<1e-10
    assert all(b['full_input_word_Gram_error']<a['full_input_word_Gram_error'] for a,b in zip(rows,rows[1:]))
    process=finite_process()
    assert max(process[k] for k in ('Hermitian_defect','local_even_commutator','unitarity_defect','probability_sum_defect','parity_defect'))<1e-12
    # CAR and omitted-contact checks include cross-region odd x odd, total even.
    car=max(np.max(np.abs(x@y+y@x)) for x,y in ((cA,cB),(cA,cB.conj().T)))
    assert car==0 and np.max(np.abs(hop@(ZA@ZB)-(ZA@ZB)@hop))==0
    return dict(round=881,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3666,
        argument_scope='Finite selected physical modes and all prescribed finite source/history coefficients admit correlated local-capacity approximants. No full-region factorization, continuum CCR in finite dimension, actual-Q identification or finite-coupling continuum convergence is claimed.',
        reference='Infinite Schmidt-rank two-mode squeezed Gaussian, q=1/4; exact infinite Wick moments',
        local_state_map='Product of number-cutoff plus vacuum-reset CPTP maps; no regional marginal product',
        finite_word_labels=LABELS,full_input_matrix_dimension=4,word_Gram_dimension=len(gram),
        maximum_Gram_polynomial_degree=max(len(w) for a in WORDS for b in WORDS for _,w,_ in mul(adj(a),b)),
        source_family='L(t,u)=t A+t^2 B+t^3 u G+t^6 u^2 C; [t^6] d_u^2 U=2i C-G^2',
        exact_qA_qB=correlation,discard_reference_correlation_error=correlation,
        omitted_source_contact_operator_norm=float(np.linalg.norm(2j*expect(C),2)),
        original_CAR_defect=float(car),infinite_Gram_min_eigenvalue=float(np.linalg.eigvalsh(gram)[0]),
        capacity_results=rows,finite_process=process,
        original_curved_QFT_numerically_solved=False,auxiliary_number_is_physical_energy=False,
        source_bounds_uniform_in_mode_resolution=False,original_QE_bridge_complete=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
