"""Read-only calibration: old reciprocal source feeding the adopted comparator.
No new numbered round; no external package installation; no plot.
"""
from fractions import Fraction as F
from itertools import permutations, combinations
import json
import numpy as np

def run():
    checks=0
    def require(condition):
        nonlocal checks
        assert bool(condition)
        checks+=1
    s7=lambda x:x-x**3/6+x**5/120-x**7/5040
    ell=min(s7(F(57,40)),s7(F(63,40)))
    B=F(6711,81920)
    delta=F(4,3)*ell**2-1-2*B
    require(delta>F(7,50))
    require(delta==F(297109212947865250239692843,2104533975040000000000000000))
    eta=F(1,32)
    margins=(F(7,50)/32,F(7,50)/32-eta/16,F(7,50)*eta/8)
    require(margins==(F(7,1600),F(31,12800),F(7,12800)))
    require(min(margins)==F(7,12800))
    error=F(7,25600)
    require(all(x-error>0 for x in margins))

    # The general proof is in screen.md. These are exact finite formula checks.
    for n in (4,6,8):
        for w in ([F(1)]*n,[F(i+1,n) for i in range(n)]):
            total=sum(w)
            row=[wi*(total-wi) for wi in w]
            for i,j in combinations(range(n),2):
                require(row[i]-row[j]==(w[i]-w[j])*(total-w[i]-w[j]))
                require(total-w[i]-w[j]>0)
            if len(set(w))==1: require(len(set(row))==1)
            else: require(len(set(row))>1)

    # Explicit primitive quantum instrument, all six branches retained.
    I2=np.eye(2); I8=np.eye(8)
    pauli=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]])
    rho=lambda v:(I2+np.einsum('a,aij->ij',v,pauli))/2
    def swap(n,i,j):
        out=np.zeros((2**n,2**n))
        for col in range(2**n):
            bits=[(col>>(n-1-k))&1 for k in range(n)]
            bits[i],bits[j]=bits[j],bits[i]
            row=sum(bit<<(n-1-k) for k,bit in enumerate(bits))
            out[row,col]=1
        return out
    SA,SB=swap(3,0,1),swap(3,0,2)
    P=np.diag([1.,0.]); Q=I2-P
    AP,AM=(I8+SA)/2,(I8-SA)/2
    BP,BM=(I8+SB)/2,(I8-SB)/2
    ks=[np.kron(P,x)/np.sqrt(2) for x in (AP,AM,BP,BM)]
    ks += [np.kron(Q,I8)/np.sqrt(2),np.kron(Q,I8)/np.sqrt(2)]
    completeness=sum(k.conj().T@k for k in ks)
    eplus=sum(ks[k].conj().T@ks[k] for k in (0,3,4))
    expect=np.eye(16)/2+np.kron(P,(SA-SB)/4)
    matrix_error=max(float(np.max(np.abs(completeness-np.eye(16)))),float(np.max(np.abs(eplus-expect))))
    require(matrix_error<1e-13)
    ps=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)
    r=np.array([.5,-.5,float(eta)])
    require(float(r@r)<1)
    primitive_errors=[]
    for h in (F(0),F(1,3),F(1)):
        for i,j in combinations(range(4),2):
            density=np.kron(np.diag([float(h),1-float(h)]),np.kron(rho(r),np.kron(rho(ps[i]),rho(ps[j]))))
            actual=float(np.trace(eplus@density).real)
            target=.5+float(h)*float(r@(ps[i]-ps[j]))/8
            primitive_errors.append(abs(actual-target))
            require(abs(actual-target)<1e-13)

    # Recompute the original 445 invariant permutation sector without altering it.
    states=list(permutations(range(4))); lookup={p:k for k,p in enumerate(states)}
    edges=list(combinations(range(4),2)); count=len(states)
    N=np.array([[int(p[i]==j and p[j]==i) for p in states] for i,j in edges])
    require(np.all(N.sum(axis=0)<=2))
    for i in range(4):
        incident=[k for k,(a,b) in enumerate(edges) if i in (a,b)]
        require(np.all(N[incident].sum(axis=0)<=1))
    h0=np.diag(4-2*N.sum(axis=0)).astype(float)
    gamma=np.zeros((count,count))
    for col,p in enumerate(states):
        for i,j in edges:
            q=list(p); q[i],q[j]=q[j],q[i]
            gamma[lookup[tuple(q)],col]+=1
    eps=1/128
    vals,vecs=np.linalg.eigh(h0+eps*gamma)
    initial=lookup[(1,0,3,2)]
    coeff=vecs[initial].conj()
    spectral=[]
    for tau in (1.9,2.,2.1):
        psi=vecs@(np.exp(-1j*tau/eps**2*vals)*coeff)
        prob=np.abs(psi)**2
        h=N@prob
        a,b=h[edges.index((0,1))],h[edges.index((0,2))]
        symmetry=max(abs(h[k]-(a if e in ((0,1),(2,3)) else b)) for k,e in enumerate(edges))
        require(symmetry<1e-9)
        require(b-a>float(delta))
        require(a+b<=1+1e-12)
        contrasts=((b-a)/32,(b-a)/32-(a+b)*float(eta)/16,b*float(eta)/8)
        require(all(x>float(y) for x,y in zip(contrasts,margins)))
        require(abs(np.sum(prob)-1)<1e-12)
        spectral.append({'tau':tau,'old_pair':round(float(a),10),'cross_pair':round(float(b),10),
             'cycle_contrasts':[round(float(x),10) for x in contrasts]})
    return {'classification':'admission_only_count_zero','new_numbered_round':False,
       'new_scientific_calibrations':0,'passed':True,'assertions':checks,
       'exact_delta_lower_bound':str(delta),'exact_cycle_lower_bounds':[str(x) for x in margins],
       'per_report_error_allowance':str(error),'instrument_branch_count':len(ks),
       'quantum_instrument_max_residual':float(f'{max(matrix_error,max(primitive_errors)):.4g}'),
       'original_445_recomputation':spectral,
       'source_gate_is_additional_instrument_permission':True,
       'no_postselection':True,'reciprocal_record_equals_physical_adjacency':False,
       'spatial_dimension_derived':False}

if __name__=='__main__':
    print(json.dumps(run(),ensure_ascii=False,indent=2))
