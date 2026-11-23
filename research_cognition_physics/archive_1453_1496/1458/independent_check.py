"""Enumerate the weights of the 15 and all vertices of the measured one-loop section."""
from itertools import combinations, combinations_with_replacement, product
import argparse
import json
from math import exp, log, pi, sqrt
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
DATA=json.loads((HERE/"../1369/inputs.json").read_text(encoding="utf-8-sig"))
MZ=DATA["mu_ref_GeV"]
SM=np.array([4.1,-19/6,-7.])


def weights():
    # (Y, T3, one SU(3) Cartan with fundamental index 1/2).
    five=np.array([[-1/3,0,.5],[-1/3,0,-.5],[-1/3,0,0],[.5,.5,0],[.5,-.5,0]])
    by_weak=np.zeros((3,3))
    for i,j in combinations_with_replacement(range(5),2):
        charge=five[i]+five[j]
        by_weak[int(i>=3)+int(j>=3)] += charge**2 * np.array([.6,1.,1.])/3
    return by_weak[::-1]  # weak triplet, leptoquark, sextet


def vertex_solutions(inv, real_triplet, light_pair=False):
    blocks=list(weights())
    if real_triplet:
        blocks.append(np.array([0.,1/3,0.]))
    columns=np.array([SM]+blocks).T
    # Unknowns are t, followed by u_j in [0,t]. Two measured differences fix a plane.
    A=np.diff(columns,axis=0)
    b=np.diff(inv)
    n=len(blocks)
    base=[]
    if light_pair:
        ell=log(1000/MZ)/(2*pi)
        for j in (1,2):
            row=np.zeros(n+1)
            row[0]=1
            row[j]=-1
            base.append((row,ell))
    hits=[]
    for chosen in combinations(range(1,n+1),n-1-len(base)):
        for upper in product([False,True],repeat=len(chosen)):
            rows=base.copy()
            for j,up in zip(chosen,upper):
                row=np.zeros(n+1)
                row[j]=1
                row[0]=-float(up)
                rows.append((row,0.))
            eq=np.vstack([A]+[r[0] for r in rows])
            targets=np.r_[b,[r[1] for r in rows]]
            if np.linalg.matrix_rank(eq)!=len(targets):
                continue
            v=np.linalg.solve(eq,targets)
            t=v[0]
            u=v[1:]
            if t<=0 or min(u)<-1.e-9 or max(u)>t+1.e-9:
                continue
            unified=inv-columns@v
            assert min(unified)>0 and np.ptp(unified)<1.e-9
            hits.append(MZ*exp(2*pi*t))
    return hits


def solve():
    names=["g1_GUT","g2","g3"]
    boxes=[(DATA[k][0]-2*DATA[k][1],DATA[k][0]+2*DATA[k][1]) for k in names]
    out={}
    for label,points in [("central",[[DATA[k][0] for k in names]]),("two_width",list(product(*boxes)))]:
        total=[]
        combined=[]
        for g in points:
            inv=4*pi/np.square(g)
            assert not vertex_solutions(inv,False)
            a1,a2,a3=inv
            x,y=a1-a2,a2-a3
            budget=25*x-32*y
            analytic_range=[(x+y-.3*budget/17)/11.4,(x+y+.6*budget/23)/11.4]
            all_m=vertex_solutions(inv,True)
            assert np.allclose([min(all_m),max(all_m)],MZ*np.exp(2*pi*np.array(analytic_range)))
            total.extend(all_m)
            subset=vertex_solutions(inv,True,light_pair=True)
            assert subset
            combined.extend(subset)
        out[label]={"unification_scale_any_three_plus_real_triplet_GeV":[min(total),max(total)],
                    "scale_if_triplet_and_leptoquark_are_both_1TeV_GeV":[min(combined),max(combined)]}
    out["complex_scalar_group_weights"]=weights().round(10).tolist()
    return json.loads(json.dumps(out),parse_float=lambda x:float(f"{float(x):.10g}"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    results=solve()
    destination=HERE/"independent_results.json"
    if args.write_new:
        with destination.open("x",encoding="utf-8") as fp:
            json.dump(results,fp,indent=2)
            fp.write("\n")
    else:
        assert results==json.loads(destination.read_text(encoding="utf-8-sig"))
        print("1458: direct weights, threshold vertices and both mass windows verified")
