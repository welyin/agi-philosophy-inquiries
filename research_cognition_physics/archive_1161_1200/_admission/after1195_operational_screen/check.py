from pathlib import Path
from fractions import Fraction
import argparse,json
import numpy as np

def rank_exact(a):
    a=[[Fraction(int(x)) for x in row] for row in a]
    r=0
    for col in range(len(a[0])):
        pivot=next((j for j in range(r,len(a)) if a[j][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        d=a[r][col];a[r]=[v/d for v in a[r]]
        for j in range(len(a)):
            if j!=r:
                d=a[j][col];a[j]=[v-d*w for v,w in zip(a[j],a[r])]
        r+=1
        if r==len(a):break
    return r

def compute():
    I=np.eye(3,dtype=np.int64);D=np.diag([0,1,2]);Q=D@D
    A=np.array([[0,1,0],[1,0,1],[0,1,0]],dtype=np.int64)
    hs=[10*I+A+Q,10*I+A-Q,10*I+D-2*A]
    H=sum(hs);g=[3,1,2];J=sum(w*h for w,h in zip(g,hs))
    comm=lambda x,y:x@y-y@x
    n2=lambda x:int(np.sum(x*x))
    bounds=[min(int(h[i,i]-sum(abs(int(h[i,j])) for j in range(3) if j!=i)) for i in range(3)) for h in hs]
    assert all(np.array_equal(h,h.T) for h in hs) and min(bounds)>0
    assert np.array_equal(H,30*I+D) and np.array_equal(J,60*I+2*D+2*Q)
    pair=[comm(hs[a],hs[b]) for a,b in [(0,1),(0,2),(1,2)]]
    actual=[comm(H,h) for h in hs]
    assert all(n2(x)>0 for x in pair+actual)
    L=np.stack([x.ravel() for x in actual],axis=1)
    assert rank_exact(np.stack([x.ravel() for x in hs],axis=1))==3
    assert rank_exact(L)==1 and np.array_equal(L@g,np.zeros(9,dtype=np.int64))
    assert np.array_equal(actual[0],actual[1]) and np.array_equal(actual[2],-2*actual[0])
    assert n2(comm(H,J))==0
    second=int(J[2,2]-2*J[1,1]+J[0,0]);assert second==4
    rhs=sum((g[b]-g[a])*comm(hs[a],hs[b]) for a,b in [(0,1),(0,2),(1,2)])
    assert np.array_equal(rhs,comm(H,J))
    return {"status":"PASS_WITH_SCOPE","check_groups":3,
      "hamiltonian_terms":[h.tolist() for h in hs],
      "strict_positive_gershgorin_bounds":bounds,"term_span_rank":3,
      "H":H.tolist(),"source_weights":g,"J":J.tolist(),
      "pair_commutator_frobenius_squared":[n2(x) for x in pair],
      "actual_exchange_frobenius_squared":[n2(x) for x in actual],
      "weight_map_rank":rank_exact(L),"weight_kernel_dimension":2,
      "kernel_equation":"g1 + g2 - 2*g3 = 0",
      "source_commutator_squared":n2(comm(H,J)),
      "source_spectral_second_difference":second,
      "scope":"Finite Hermitian source-class counterexample only; no spacetime locality, physical monopole instrument, adopted axiom, or GR counterexample."}

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true");args=parser.parse_args()
    out=compute();p=Path(__file__).with_name("results.json")
    if args.write:
        with p.open("x",encoding="utf-8") as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:assert json.loads(p.read_text(encoding="utf-8"))==out
    print(json.dumps(out,ensure_ascii=False))
