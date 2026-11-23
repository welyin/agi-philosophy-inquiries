"""Check finite-SU(5) traces, isolated-family cuts and light-doublet budgets."""
import argparse
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def named_models(which):
    # Coefficients are squared couplings in units g^2, in the conventional normalization.
    u=np.zeros((4,3,3)); d=np.zeros_like(u)
    for i in range(3):
        factor=.5 if which=="B" and i>0 else 1
        u[i,i,i]=math.sqrt(1.6*factor)
        d[i,i,i]=math.sqrt(1.2*factor)
    if which=="A":
        f=np.diag([0,0,0,1])
    else:
        f=np.diag([0,math.sqrt(.5),math.sqrt(.5),0])
        u[3,1,2]=u[3,2,1]=math.sqrt(.8)
        d[3,1,2]=d[3,2,1]=math.sqrt(.6)
    up=np.einsum("aij,bij->ab",u,u)
    dn=np.einsum("aij,bij->ab",d,d)
    errors=[3*up+4.8*f.T@f-4.8*np.eye(4),
            4*dn+4.8*f@f.T-4.8*np.eye(4),
            4*np.einsum("aji,ajk->ik",d,d)-4.8*np.eye(3),
            3*np.einsum("aij,akj->ik",u,u)+2*np.einsum("aij,akj->ik",d,d)-7.2*np.eye(3)]
    assert max(np.max(abs(m)) for m in errors)<1e-12
    assert np.isclose(np.sum(f*f)+21/5*(15/7),10)
    assert np.allclose(f@f.T,np.diag(np.diag(f@f.T)))
    return dict(f_norm2=float(np.sum(f*f)),adjoint_rank=int(np.linalg.matrix_rank(f)),
                up_frobenius2=float(np.sum(u*u)),down_frobenius2=float(np.sum(d*d)),
                f_singular_squared=sorted(np.diag(f@f.T).tolist(),reverse=True))


def run():
    cfg=read(HERE/"inputs.json")
    data=read(HERE/cfg["sources"]["weak_tree"])
    expected=read(HERE/cfg["sources"]["weak_tree_verification"])["direct_edges"]["min_amplitudes"]
    tree={}
    for edge,central,errors in data["edges"]:
        tree[edge]=(central-2*sum(errors))/1.02
    assert all(math.isclose(x,expected[k]) for k,x in tree.items())
    cuts={}
    for row,col in itertools.product("uct","dsb"):
        rowleak=math.sqrt(sum(a*a for edge,a in tree.items() if edge[0]==row and edge[1]!=col))
        colleak=math.sqrt(sum(a*a for edge,a in tree.items() if edge[1]==col and edge[0]!=row))
        cuts[row+col]=max(rowleak,colleak)
    winner=min(cuts,key=cuts.get)
    assert winner=="tb" and math.isclose(cuts[winner],tree["cb"])
    # If all directions are isolated, the best permutation yields the worst per-edge/row lower bound.
    diagonal={''.join(p):max(cuts[i+j] for i,j in zip("uct",p)) for p in itertools.permutations("dsb")}
    field_counts=[]
    for N in range(1,6):
        n=10-2*N
        field_counts.append(dict(N=N,n=n,required_f_norm2_in_g2=n-N,
                                 trace_compatible=n>=N,adjoint_self_squared_in_g2=float(F(5*N,7))))
    ## The 3 independent trace residuals cancel all U,L exactly for arbitrary integers.
    for N,n,U,L,f in [(3,4,F(24,5),F(18,5),1),(4,2,F(5),F(5),0),(2,6,2,3,4)]:
        r10=3*U+2*L-F(36,5)*N
        r5=3*U+F(24,5)*(f-n)
        rb=8*L+F(24,5)*(f-N-n)
        assert r5-r10+rb/4==6*(f-n+N)
    w=cfg["projection_witness_in_reduced_units"]
    a,b=np.array(w["up"]),np.array(w["down"])
    sb2=w["sin_beta_squared"]
    roots=[]
    for vector,weight in [(a,sb2),(b,1-sb2)]:
        light=vector/math.sqrt(weight)
        assert light@light<1
        h=np.append(light,math.sqrt(1-light@light))
        assert np.isclose(h@h,1)
        roots.append(h.tolist())
    return dict(count_table=field_counts,models={k:named_models(k) for k in "AB"},
                individual_cut_operator_lower=cuts,all_isolated_permutation_lower=diagonal,
                most_permissive_isolated_pair=winner,best_block_rearrangement=cuts[winner],
                best_all_isolated_lower=min(diagonal.values()),
                minimum_any_common_active_threshold_rotation_radians=math.asin(cuts[winner]),
                N4_common_trace_residual_per_field_min=float(F(8,5)),
                reduced_light_Higgs_directions=roots,combined_projected_budget=float(a@a+b@b))


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return math.isclose(a,b,rel_tol=3e-12,abs_tol=1e-13)
    return a==b


if __name__=="__main__":
    args=argparse.ArgumentParser()
    mode=args.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    opt=args.parse_args()
    ans=run()
    dest=HERE/"results.json"
    if opt.write:
        with dest.open("x",encoding="utf-8") as f:
            json.dump(ans,f,indent=2)
            f.write("\n")
    else:
        assert same(ans,read(dest))
        print("1479: finite SU(5) traces, A/B tensors, nine cuts and projections passed.")
