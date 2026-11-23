"""Permission-independent one-loop trace check, sharp four-family balance and FUTA mass lemma."""
import argparse
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

import numpy as np


def su5_invariants(labels):
    # Partition row lengths, Weyl's formula and C2 with T(5)=1/2.
    rows = [sum(labels[i:]) for i in range(4)] + [0]
    boxes = sum(rows)
    dim = Q(1)
    for i, j in itertools.combinations(range(5), 2):
        dim *= Q(rows[i] - rows[j] + j - i, j - i)
    cas = (5 * boxes + sum(r * (r + 1 - 2*(i+1)) for i,r in enumerate(rows))
           - Q(boxes * boxes, 5)) / 2
    return dim, cas, dim * cas / 24


def sharp_gap():
    # Four symmetric-family tensors coupled to two 5's; antisym pairs below.
    a = (48 / 95)**.5
    b = (36 / 95)**.5
    up = a * np.array([np.eye(4),np.diag([1,1j,-1,-1j])])
    # Four of the five matchings of K6. Index 5 is the fixed circle vertex.
    lam = np.zeros((4,6,6),complex)
    for k in range(4):
        pairs = [(5,k)]+[((k+i)%5,(k-i)%5) for i in (1,2)]
        for i,j in pairs:
            lam[k,i,j],lam[k,j,i] = b,-b
    cu = 3*np.einsum("aij,bij->ab",up,up.conj())
    cd = 4*np.einsum("iAB,iCB->AC",lam,lam.conj())
    ct = 3*np.einsum("aij,akj->ik",up,up.conj())+np.einsum("iAB,kAB->ik",lam,lam.conj())
    ratios = []
    for block,gauge in [(cu,24/5),(cd,24/5),(ct,36/5)]:
        ratio = np.real(np.diag(block))/gauge
        assert np.allclose(block, np.eye(len(block))*ratio[0]*gauge)
        ratios.append(float(ratio[0]-1))
    assert np.allclose(ratios,[5/19,5/19,-5/19])
    return ratios


def mass_test():
    # Four toy GUT-scale masses in one common unit; real Householder rotations
    # mixing family and fourth pairs. A small holomorphic mass can be restored.
    def house(x):
        x=np.array(x,dtype=float)
        x/=np.linalg.norm(x)
        return np.eye(4)-2*np.outer(x,x)
    left,right=house([1,2,-1,5]),house([3,-2,4,1])
    heavy=np.array([5.,2.,.8])
    delta=3.7
    a,b=abs(left[-1,-1]),abs(right[-1,-1])
    effective={}
    for mu in (0.,.007):
        md=left @ np.diag([*heavy,mu]) @ right.T
        mt=md.copy()
        mt[-1,-1]+=delta
        measured=abs(np.linalg.det(mt))/np.prod(heavy)
        primary=abs(delta)*a*b
        bound=mu*(1+abs(delta)/min(heavy)*((1-a*a)*(1-b*b))**.5)
        assert abs(measured-primary)<bound+1e-12
        effective[str(mu)]={"triplet_threshold":measured,"rank_three_value":primary,
                            "finite_mu_abs_error_bound":bound}
    fup=(8/5)*sum(left[:-1,-1]**2)
    fdn=(6/5)*sum(right[:-1,-1]**2)
    assert np.isclose(delta*((1-fup/(8/5))*(1-fdn/(6/5)))**.5,delta*a*b)
    return effective


def calculation():
    indices={}
    for i,j in itertools.combinations_with_replacement(range(4),2):
        lab=[0]*4
        lab[i]+=1
        lab[j]+=1
        dim,cas,t=su5_invariants(lab)
        indices["".join(map(str,lab))]=[int(dim),float(cas),float(t)]
    assert min(v[-1] for v in indices.values()) == 3.5
    coef=np.array([[2.,1,0],[1,0,1],[0,2,1]])
    rhs={n: np.linalg.solve(coef,[72*n,24*(10-2*n),24*(10-n)]).tolist()
         for n in range(1,6)}
    assert all(np.isclose(vals[-1],24*(10-3*n)) for n,vals in rhs.items())
    return {"twice_fundamental_weights_dim_C2_T":indices,
            "colour_weighted_up_down_adjoint_traces":rhs,
            "sharp_four_family_relative_errors_5_bar5_10":sharp_gap(),
            "FUTA_example":mass_test()}


def rounded(v):
    if isinstance(v,dict): return {str(k):rounded(x) for k,x in v.items()}
    if isinstance(v,list): return [rounded(x) for x in v]
    return round(float(v),10)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    flag=cli.add_mutually_exclusive_group()
    flag.add_argument("--check",action="store_true")
    flag.add_argument("--write-new",action="store_true")
    args=cli.parse_args()
    path=Path(__file__).with_name("independent_results.json")
    result=rounded(calculation())
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            f.write(json.dumps(result,indent=2)+"\n")
    elif args.check:
        assert json.loads(path.read_text(encoding="utf-8"))==result
        print("1479 independent: all cubic traces, sharp four-family defect and colour threshold verified")
    else:
        print(json.dumps(result,indent=2))
