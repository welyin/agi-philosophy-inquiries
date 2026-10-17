"""Independent 1273 check: exact four-bit algebra; no author module imports.

One independent group. The rational matrix diagnostics are not a construction
of new physical stress, controls, clocks or geometry. Analytic quantifiers are
reviewed separately. Default mode is read-only; first creation is explicit.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent

class Q:
    """Exact Gaussian rational."""
    __slots__ = ("r", "i")
    def __init__(self, r=0, i=0):
        self.r, self.i = F(r), F(i)
    def __add__(self, y):
        y = asq(y)
        return Q(self.r+y.r, self.i+y.i)
    __radd__ = __add__
    def __neg__(self):
        return Q(-self.r, -self.i)
    def __sub__(self, y):
        return self + (-asq(y))
    def __mul__(self, y):
        y = asq(y)
        return Q(self.r*y.r-self.i*y.i, self.r*y.i+self.i*y.r)
    __rmul__ = __mul__
    def conj(self):
        return Q(self.r, -self.i)
    def norm2(self):
        return self.r*self.r+self.i*self.i
    def __eq__(self, y):
        y=asq(y)
        return self.r == y.r and self.i == y.i
    def __pow__(self, n):
        assert n >= 0
        ans=Q(1)
        for _ in range(n):
            ans=ans*self
        return ans

def asq(x):
    return x if isinstance(x, Q) else Q(x)

def zeros(n):
    return [[Q() for _ in range(n)] for _ in range(n)]

def eye(n):
    a=zeros(n)
    for i in range(n):
        a[i][i]=Q(1)
    return a

def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), Q())
             for j in range(len(b[0]))] for i in range(len(a))]

def mv(a,v):
    return [sum((x*y for x,y in zip(row,v)),Q()) for row in a]

def dagger(a):
    return [[a[j][i].conj() for j in range(len(a))]
            for i in range(len(a[0]))]

def rank_exact(rows):
    a=[[F(v) for v in row] for row in rows if any(row)]
    if not a:
        return 0
    row=0
    for col in range(len(a[0])):
        candidates=[j for j in range(row,len(a)) if a[j][col]]
        if not candidates:
            continue
        j=candidates[0]
        a[j],a[row]=a[row],a[j]
        z=a[row][col]
        a[row]=[v/z for v in a[row]]
        for j in range(len(a)):
            if j != row:
                z=a[j][col]
                a[j]=[v-z*w for v,w in zip(a[j],a[row])]
        row+=1
        if row == len(a):
            break
    return row

BITS=list(product((0,1),repeat=4))
INDEX={x:j for j,x in enumerate(BITS)}
ENERGY=[sum(x) for x in BITS]
SECTOR=[j for j,e in enumerate(ENERGY) if e==2]
EDGES=((0,2),(1,3),(0,1))

def swap(i,j):
    a=zeros(16)
    for k,b in enumerate(BITS):
        c=list(b)
        c[i],c[j]=c[j],c[i]
        a[INDEX[tuple(c)]][k]=Q(1)
    return a

def half_swap(i,j):
    s=swap(i,j)
    return [[Q(F(1,2),F(-1,2))*(int(r==c))+
             Q(F(1,2),F(1,2))*s[r][c] for c in range(16)]
            for r in range(16)]

def apply_heat(v,u):
    out=[Q() for _ in range(32)]
    for o in range(2):
        part=mv(u,v[o*16:(o+1)*16])
        out[o*16:(o+1)*16]=part
    return out

def geometry(n,step):
    f=n*n+3*n+7
    if step==0:
        z=Q(F(3,5),F(4,5))**f
        rot=[[Q(F(3,5)),Q(F(-4,5))],[Q(F(4,5)),Q(F(3,5))]]
    else:
        z=Q(F(5,13),F(12,13))**f
        rot=[[Q(F(5,13)),Q(F(-12,13))],[Q(F(12,13)),Q(F(5,13))]]
    return mm(rot,[[z,Q()],[Q(),z.conj()]])

def apply_geometry(v,step):
    out=[Q() for _ in range(32)]
    for n in range(16):
        z=mv(geometry(ENERGY[n],step),[v[n],v[16+n]])
        out[n],out[16+n]=z
    return out

def compute():
    assert len(SECTOR)==6
    edge_records=[]
    all_rows=[]
    for i,j in EDGES:
        s=swap(i,j)
        assert mm(dagger(s),s)==eye(16)
        rows=set()
        for col,b in enumerate(BITS):
            row=next(r for r in range(16) if s[r][col]==1)
            assert ENERGY[row]==ENERGY[col]
            if col in SECTOR:
                diff=tuple(BITS[row][k]-b[k] for k in range(4))
                rows.add(diff)
                all_rows.append(diff)
        nonzero=sorted(v for v in rows if any(v))
        assert len(nonzero)==2
        edge_records.append({"edge":[i,j],"sector_response_rows":[list(v) for v in nonzero]})
    rank=rank_exact(all_rows)
    assert rank==3
    # The surviving coefficient vector is constant: rank three plus this vector.
    assert all(sum(row)==0 for row in all_rows)
    u=half_swap(0,1)
    assert mm(dagger(u),u)==eye(16)
    assert mm(u,u)==swap(0,1)
    for r in range(16):
        for c in range(16):
            if u[r][c] != 0:
                assert ENERGY[r]==ENERGY[c]

    # This is the reduced four-energy-bit density after tracing the retained
    # orthogonal seed labels, not a replacement for the complete preparation.
    weights={}
    for a,b in product((0,1),repeat=2):
        weights[INDEX[(a,b,1-a,1-b)]]=(F(1,4) if a else F(3,4))*(F(2,5) if b else F(3,5))
    assert sum(weights.values())==1
    before=[sum((p*BITS[j][k] for j,p in weights.items()),F()) for k in range(4)]
    after_prob=[sum((p*u[i][j].norm2() for j,p in weights.items()),F()) for i in range(16)]
    after=[sum((after_prob[j]*BITS[j][k] for j in range(16)),F()) for k in range(4)]
    delta=[b-a for a,b in zip(before,after)]
    assert before==[F(1,4),F(2,5),F(3,4),F(3,5)]
    assert delta==[F(3,40),F(-3,40),F(0),F(0)]
    for c in ((1,1,1,1),(1,0,0,0),(1,-1,0,0),(2,3,5,7)):
        assert sum((F(x)*y for x,y in zip(c,delta)),F())==F(3,40)*(c[0]-c[1])
    # Exact sector factorization with two noncommuting geometry steps and
    # noncommuting heat steps. This is an algebraic diagnostic, not a source model.
    g=mm(geometry(2,1),geometry(2,0))
    assert mm(geometry(2,1),geometry(2,0)) != mm(geometry(2,0),geometry(2,1))
    programs=((u,swap(0,2)),(eye(16),half_swap(1,3)))
    checked_columns=0
    for first,second in programs:
        h=mm(second,first)
        for o in (0,1):
            for n in SECTOR:
                v=[Q() for _ in range(32)]
                v[o*16+n]=Q(1)
                v=apply_heat(apply_geometry(v,0),first)
                v=apply_heat(apply_geometry(v,1),second)
                expected=[g[oo][o]*h[nn][n] for oo in range(2) for nn in range(16)]
                assert v==expected
                checked_columns+=1
    # Analytic preparation/perturbation budget from two applications of Duhamel.
    d1=d2=F(1,10**9)
    eta1,eta2=F(1,10000),F(1,20000)
    budget=d1+d2+eta1+eta2
    assert budget==F(75001,500000000)
    return {
        "round":1273,
        "verdict":"PASS_INDEPENDENT_SCOPED",
        "independent_groups":1,
        "author_code_imported":False,
        "method":"exact Gaussian rationals and four-energy-bit swaps; no numerical quadrature",
        "basis_dimension":16,
        "sector_dimension":6,
        "spectator_degeneracy":960,
        "inherited_full_sector_dimension":6*960,
        "exchange_checks":edge_records,
        "derived_additive_constraint_rank":rank,
        "additive_kernel_dimension":4-rank,
        "half_swap_unitary_and_square_identity":True,
        "first_exchange_K1_before":[str(x) for x in before],
        "first_exchange_K1_after":[str(x) for x in after],
        "first_exchange_K1_delta":[str(x) for x in delta],
        "first_exchange_K0_delta":["0"]*4,
        "noncommuting_sector_factorization_exact_columns":checked_columns,
        "unknown_reference_scope":"exact operator identity extends with I_R; no new preparation permission",
        "duhamel_budget_arithmetic":str(budget),
        "physical_new_source_or_geometry_certified":False,
        "old_3_over_40_not_transported_to_noncentral_coupling":True,
        "code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write","--write-once",dest="write",action="store_true")
    args=parser.parse_args()
    result=compute()
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf-8",newline="\n") as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
        mode="FIRST_WRITE"
    else:
        assert result==json.loads(target.read_text(encoding="utf-8")),"saved result mismatch"
        mode="READ_ONLY_MATCH"
    print(json.dumps({"round":1273,"verdict":result["verdict"],"mode":mode,
                      "exact_factorization_columns":result["noncommuting_sector_factorization_exact_columns"],
                      "additive_rank":result["derived_additive_constraint_rank"],
                      "independent_groups":1},ensure_ascii=False))
