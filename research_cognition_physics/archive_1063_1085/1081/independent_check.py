"""Exact independent direction-network audit; no imports from author checker."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json, sys

def rank(a):
    a=[[F(x) for x in row] for row in a]
    r=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
        if pivot is None: continue
        a[r],a[pivot]=a[pivot],a[r]
        q=a[r][col]; a[r]=[x/q for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                q=a[i][col];a[i]=[x-q*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r

def delta(a,b):return tuple(F(x)-F(y) for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def scaled(s,a):return tuple(F(s)*x for x in a)
def edge(D,i,j):return D[i,j] if i<j else neg(D[j,i])
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def direction_rank(D,n):
    rows=[]
    for (i,j),d in D.items():
        c=((0,-d[2],d[1]),(d[2],0,-d[0]),(-d[1],d[0],0))
        for axis in range(3):
            row=[F(0)]*(3*(n-1))
            for node,sign in ((i,1),(j,-1)):
                if node:
                    for k in range(3):row[3*(node-1)+k]+=sign*c[axis][k]
            rows.append(row)
    return rank(rows)

def run():
    D={(0,1):(-1,0,0),(0,2):(0,-1,0),(0,3):(-2,-3,0),
       (1,2):(1,-1,0),(1,3):(-1,-3,0),(2,3):(-2,F(-5,2),0)}
    triangles=[((0,1,2),(1,1,1)),((0,1,3),(1,1,1)),
               ((0,2,3),(1,2,2)),((1,2,3),(7,8,9))]
    for (i,j,k),w in triangles:
        vecs=(edge(D,i,j),edge(D,j,k),edge(D,k,i))
        assert all(sum(w[t]*vecs[t][s] for t in range(3))==0 for s in range(3))
        assert min(w)>0
    rank2=direction_rank(D,4)
    assert rank2==9 and rank(list(D.values()))==2
    max_norm=max(sum(F(x)**2 for x in d) for d in D.values())/64
    assert max_norm==F(13,64)<F(1,4)
    # All cross-block edges A -> B point e1; only one internal edge per block.
    A={0,2}; B={1,3}; e1=(1,0,0);e2=(0,1,0);e3=(0,0,1)
    S={(i,j):(e2 if (i,j)==(0,2) else e3 if (i,j)==(1,3)
        else e1 if i in A else neg(e1)) for i,j in combinations(range(4),2)}
    assert rank(list(S.values()))==3 and direction_rank(S,4)==8
    for i in A:
        for j in B:assert edge(S,i,j)==e1
    dot=lambda r,d:sum(a*b for a,b in zip(r,d))
    assert dot(e2,edge(S,0,1))==dot(e2,edge(S,1,2))==0
    assert dot(e2,edge(S,0,2))==1
    # Distinct method: solve all cross-product direction constraints simultaneously.
    V=[(0,0,0),(1,0,0),(0,1,0),(0,0,1),(2,0,0),(0,0,0)]
    gains={(i,j):F(i+j+2,j-i+1) for i,j in combinations(range(6),2)}
    P={ij:scaled(a,delta(V[ij[0]],V[ij[1]])) for ij,a in gains.items()}
    assert rank(list(P.values()))==3 and direction_rank(P,6)==14
    for i,j,k in combinations(range(6),3):
        vecs=(edge(P,i,j),edge(P,j,k),edge(P,k,i))
        w=(1/gains[i,j],1/gains[j,k],1/gains[i,k])
        assert all(sum(w[t]*vecs[t][s] for t in range(3))==0 for s in range(3))
    assert P[0,5]==(0,0,0)
    # Continuous 1D polynomial curve has three independent difference vectors.
    C=[(t,t*t,t*t*t) for t in (0,1,2,3)]
    assert rank([delta(C[i],C[0]) for i in (1,2,3)])==3
    return dict(round=1081,passed=True,arithmetic="Fraction exact",assertion_groups=5,
      rank2_counterexample=dict(span=2,anchored_direction_constraint_rank=rank2,
        unknown_coordinates=9,positive_triangle_certificates=4,max_effect_norm_squared=str(max_norm)),
      strict_order_counterexample=dict(span=3,anchored_direction_constraint_rank=8,
        positive_scale_obstruction="cross-block edges force each block to coincide; internal nonzero edge cannot have positive difference gain",
        tie_witness="r=e2: 0~1, 1~2, 0>2",all_r_strict_acyclic_reason="uniform cross-block orientation; blocks of size two"),
      positive_network=dict(nodes=6,span=3,unknown_coordinates=15,constraint_rank=14,
        parallel_stars=True,zero_edge=True,positive_triangle_certificates=20),
      dimension_boundary="curve t -> (t,t^2,t^3) has affine span 3 and dimension 1",
      scope="exact finite certificates and implementation audit; general theorem rests on proof.md; no actual spatial dimension conclusion")

if __name__=="__main__":
    result=run();p=Path(__file__).with_name("independent_results.json")
    if "--write" in sys.argv:
        assert not p.exists();p.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    elif "--check" in sys.argv:assert result==json.loads(p.read_text(encoding="utf-8"))
    print(json.dumps(result,ensure_ascii=False,indent=2))
