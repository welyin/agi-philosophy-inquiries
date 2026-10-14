"""Unnumbered EPRL boundary probe. Exact CG sums in a small radical field.

No SymPy dependency or imported spin-foam formula. NumPy spectra are diagnostic.
Default compares the saved JSON without writing; --write creates it once.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product, permutations
import json
from math import factorial, gcd, sqrt
from pathlib import Path
import os
for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[name] = "1"
import numpy as np


class R:
    """Exact real sum q[d] sqrt(d), d positive and squarefree."""
    __slots__ = ("q",)
    def __init__(self, q=0):
        self.q = {d: Q(v) for d, v in q.items() if v} if isinstance(q, dict) else ({1: Q(q)} if q else {})
    def __add__(self, other):
        other = other if isinstance(other, R) else R(other)
        out = dict(self.q)
        for d, v in other.q.items():
            out[d] = out.get(d, Q(0)) + v
        return R(out)
    __radd__ = __add__
    def __neg__(self):
        return R({d: -v for d,v in self.q.items()})
    def __sub__(self, other):
        return self + (-other if isinstance(other,R) else -R(other))
    def __rsub__(self, other):
        return R(other) - self
    def __mul__(self, other):
        other = other if isinstance(other, R) else R(other)
        out = {}
        for d,a in self.q.items():
            for e,b in other.q.items():
                g = gcd(d,e)
                f = d*e//(g*g)
                out[f] = out.get(f,Q(0)) + a*b*g
        return R(out)
    __rmul__ = __mul__
    def __eq__(self, other):
        other = other if isinstance(other,R) else R(other)
        return self.q == other.q
    def __bool__(self):
        return bool(self.q)
    def value(self):
        return sum(float(v)*sqrt(d) for d,v in self.q.items())
    def data(self):
        return {str(d): f"{v.numerator}/{v.denominator}" for d,v in sorted(self.q.items())}
    def __repr__(self):
        return " + ".join(f"({v})" + (f"*sqrt({d})" if d!=1 else "") for d,v in sorted(self.q.items())) or "0"


def root(q):
    q = Q(q)
    assert q >= 0
    if not q:
        return R()
    n, outside, inside, p = q.numerator*q.denominator, 1, 1, 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        outside *= p**(e//2)
        if e%2:
            inside *= p
        p += 1
    inside *= n
    ans = R({inside: Q(outside,q.denominator)})
    assert ans*ans == R(q)
    return ans


def integer(q):
    q = Q(q)
    assert q.denominator == 1
    return q.numerator


@lru_cache(None)
def cg(j1,m1,j2,m2,j,m):
    if m1+m2 != m or abs(m)>j or j>j1+j2 or j<abs(j1-j2):
        return R()
    facts = [j+j1-j2, j-j1+j2, j1+j2-j, j1+j2+j+1,
             j+m,j-m,j1-m1,j1+m1,j2-m2,j2+m2]
    if any(x<0 or Q(x).denominator!=1 for x in facts):
        return R()
    f = [factorial(integer(x)) for x in facts]
    rad = Q((2*j+1)*f[0]*f[1]*f[2],f[3])
    for v in f[4:]:
        rad *= v
    sm = Q(0)
    for z in range(integer(j1+j2+j)+2):
        args = [z,j1+j2-j-z,j1-m1-z,j2+m2-z,
                j-j2+m1+z,j-j1-m2+z]
        if any(x<0 or Q(x).denominator!=1 for x in args):
            continue
        den = 1
        for x in args:
            den *= factorial(integer(x))
        sm += Q((-1)**z,den)
    return root(rad)*sm


def magnetic(j):
    return [j-i for i in range(integer(2*j)+1)]


@lru_cache(None)
def invariant(j,a,ms):
    if sum(ms)!=0:
        return R()
    m=ms[0]+ms[1]
    if abs(m)>a:
        return R()
    return cg(j,ms[0],j,ms[1],a,m)*cg(j,ms[2],j,ms[3],a,-m)*((-1)**integer(a-m))*root(Q(1,2*a+1))


def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),R())


def mm(a,b):
    return [[dot(row,col) for col in zip(*b)] for row in a]


def transpose(a):
    return [list(x) for x in zip(*a)]


def eye(n):
    return [[R(int(i==j)) for j in range(n)] for i in range(n)]


def subtract(a,b):
    return [[x-y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]


def determinant(a):
    n=len(a)
    out=R()
    for perm in permutations(range(n)):
        inv=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=R((-1)**inv)
        for i,j in enumerate(perm):
            term=term*a[i][j]
        out=out+term
    return out


def encode(a):
    return [[x.data() for x in row] for row in a]


def diagnostic(a):
    return np.array([[x.value() for x in row] for row in a],dtype=float)


def calculate():
    j, jp, jm = Q(3,2), Q(1), Q(1,2)
    assert (Q(1)+Q(1,3))*j/2==jp and (Q(1)-Q(1,3))*j/2==jm
    m0=list(product(magnetic(j),repeat=4))
    mp=[x for x in product(magnetic(jp),repeat=4) if sum(x)==0]
    mn=[x for x in product(magnetic(jm),repeat=4) if sum(x)==0]
    basis=[[invariant(j,Q(a),ms) for ms in m0] for a in range(4)]
    assert mm(basis,transpose(basis))==eye(4)
    plus=[[invariant(jp,Q(a),ms) for ms in mp] for a in range(3)]
    minus=[[invariant(jm,Q(b),ms) for ms in mn] for b in range(2)]
    assert mm(plus,transpose(plus))==eye(3)
    assert mm(minus,transpose(minus))==eye(2)
    ym=[[(cg(jp,x,jm,y,j,m)) for m in magnetic(j)]
        for x,y in product(magnetic(jp),magnetic(jm))]
    assert mm(transpose(ym),ym)==eye(4)
    rows=[(a,b) for a in range(3) for b in range(2)]
    fusion=[[R() for _ in range(4)] for _ in rows]
    for pi,pms in enumerate(mp):
        for ni,nms in enumerate(mn):
            ms=tuple(x+y for x,y in zip(pms,nms))
            yy=R(1)
            for x,y,m in zip(pms,nms,ms):
                yy=yy*cg(jp,x,jm,y,j,m)
            if not yy:
                continue
            for row,(a,b) in enumerate(rows):
                weight=plus[a][pi]*minus[b][ni]*yy
                if weight:
                    for col in range(4):
                        fusion[row][col]=fusion[row][col]+weight*invariant(j,Q(col),ms)
    gram=mm(transpose(fusion),fusion)
    det=determinant(gram)
    assert det
    assert det.q.keys()=={1} and det.q[1]>0
    assert any(gram[i][j] for i in range(4) for j in range(4) if i!=j)
    q12=[[R(i*(i+1) if i==j2 else 0) for j2 in range(4)] for i in range(4)]
    permutation23=[[sum((basis[i][k]*invariant(j,Q(a),(ms[0],ms[2],ms[1],ms[3]))
                               for k,ms in enumerate(m0)),R()) for a in range(4)] for i in range(4)]
    assert mm(permutation23,transpose(permutation23))==eye(4)
    q13=mm(mm(permutation23,q12),transpose(permutation23))
    c12=subtract(mm(gram,q12),mm(q12,gram))
    c13=subtract(mm(gram,q13),mm(q13,gram))
    assert any(x for row in c12 for x in row)
    assert any(x for row in c13 for x in row)
    edges=[(i,k) for i in range(4) for k in range(i+1,4) if q13[i][k]]
    reached={0}
    for _ in range(4):
        reached.update(k for i,k in edges if i in reached)
        reached.update(i for i,k in edges if k in reached)
    assert len(reached)==4
    # q12 has four distinct eigenvalues, so its commutant is diagonal.
    # Connected offdiagonal graph of q13 forces all four diagonal entries equal.
    g=diagnostic(gram)
    vals,vec=np.linalg.eigh(g)
    f=diagnostic(fusion)
    polar=f@((vec/np.sqrt(vals))@vec.T)
    raw_normalized_overlap=g/np.sqrt(np.outer(np.diag(g),np.diag(g)))
    return {
        "status":"UNNUMBERED_PREVIEW", "gamma":"1/3", "k_four_legs":["3/2"]*4,
        "j_plus":"1/1", "j_minus":"1/2", "input_dimension":4, "output_dimension":6,
        "basis":"(12)->a,(34)->a, coupled to singlet; Condon-Shortley CG convention",
        "fusion_row_order_a_plus_a_minus":[list(x) for x in rows],
        "fusion_column_order_a":[0,1,2,3],
        "radical_encoding":"entry dictionary d -> rational coefficient of sqrt(d)",
        "fusion_exact":encode(fusion), "gram_exact":encode(gram),
        "gram_determinant_exact":det.data(),
        "single_leg_Y_isometry_exact":True,"three_invariant_bases_orthonormal_exact":True,
        "gram_non_scalar_exact":True,"gram_rank_exact":4,
        "q12_exact":encode(q12),"swap23_exact":encode(permutation23),"q13_exact":encode(q13),
        "commutator_G_q12_exact":encode(c12),"commutator_G_q13_exact":encode(c13),
        "q13_nonzero_offdiagonal_graph": [list(x) for x in edges],
        "joint_complex_commutant_dimension_exact":1,
        "commutant_proof":"q12 simple spectrum + connected q13 offdiagonal graph",
        "diagnostic_only":{
            "gram_matrix":g.tolist(),"gram_eigenvalues":vals.tolist(),
            "polar_isometry_residual":float(np.linalg.norm(polar.T@polar-np.eye(4))),
            "raw_normalized_image_overlaps":raw_normalized_overlap.tolist(),
            "raw_success_probability_range":vals[[0,-1]].tolist(),
            "G_q12_commutator_norm":float(np.linalg.norm(diagnostic(c12),2)),
            "G_q13_commutator_norm":float(np.linalg.norm(diagnostic(c13),2))
        },
        "scope":{
            "mature_nonisometry_source":"https://arxiv.org/html/0912.0540v2",
            "corrected_injectivity_source":"https://arxiv.org/html/1109.5023",
            "admissibility":"sum j_plus=4 integer; both target invariant spaces nonzero",
            "not_a_CPTP_spacetime_vertex":True,"no_actual_clock_or_source_backreaction_claim":True,
            "no_installed_dependency":True,"no_quantum_gravity_failure_claim":True
        }
    }


def compare(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(compare(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return isinstance(b,list) and len(a)==len(b) and all(compare(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return isinstance(b,(int,float)) and abs(a-b)<=2e-13*max(1,abs(a),abs(b))
    return a==b


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    args=ap.parse_args()
    result=calculate()
    target=Path(__file__).with_name("probe_results.json")
    if args.write:
        with target.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
        label="WRITE_ONCE_PASS"
    else:
        assert compare(result,json.loads(target.read_text(encoding="utf-8")))
        label="READ_ONLY_MATCH"
    print(label)
    print("G =", result["gram_exact"])
    print("eigenvalues(diagnostic) =",result["diagnostic_only"]["gram_eigenvalues"])
    print("joint commutant dimension =",result["joint_complex_commutant_dimension_exact"])


if __name__=="__main__":
    main()
