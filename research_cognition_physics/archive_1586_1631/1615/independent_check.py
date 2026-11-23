"""Independent D14 character and physical qqqL permutation checks."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def representation(k,reflected,harmonic=1,twist=False):
    a=2*math.pi*k*harmonic/14
    plane=np.array([[math.cos(a),-math.sin(a)],[math.sin(a),math.cos(a)]])
    if reflected:
        plane=plane@np.diag([1.,-1.])
    g=np.eye(3)
    g[1:,1:]=plane
    return -g if twist and reflected else g


def flavour_relations():
    def idx(i,j,k):
        return 9*i+3*j+k
    rows=[]
    for i,j,k in itertools.product(range(3),repeat=3):
        row=np.zeros(27)
        for pos,sign in [((i,j,k),1),((j,i,k),1),((k,i,j),-1),((k,j,i),-1)]:
            row[idx(*pos)]+=sign
        rows.append(row)
    return np.array(rows)


def area(x,y):
    return x[1]*y[2]-x[2]*y[1]


def mixed(x,y,z):
    return x[0]*area(y,z)-z[0]*area(x,y)


def run():
    groups=[[representation(k,ref,twist=twist)
             for ref in (False,True) for k in range(14)] for twist in (False,True)]
    constraints=flavour_relations()
    assert np.linalg.matrix_rank(constraints,tol=1e-9)==8
    ans={}
    for label,group in zip(("ordinary","twisted"),groups):
        characters=[]
        cube=np.zeros((27,27))
        for g in group:
            x,y,z=(np.trace(np.linalg.matrix_power(g,j)) for j in (1,2,3))
            characters.append([(x**3+3*x*y+2*z)/6,(x**3-z)/3,(x**3-3*x*y+2*z)/6])
            cube+=np.kron(g,np.kron(g,g))/28
        schur=np.rint(np.mean(characters,axis=0)).astype(int)
        # Orthogonal representatives of the operator quotient satisfy R c=0.
        rank=27-np.linalg.matrix_rank(np.vstack([cube-np.eye(27),constraints]),tol=1e-9)
        assert rank==sum(schur)
        ans[label]={"symmetric_mixed_alternating":schur.tolist(),
                    "duql_three_trivial_leptons":int(round(np.trace(cube)))*3,
                    "qqql_three_trivial_leptons":int(rank)*3}
    assert ans["ordinary"]["qqql_three_trivial_leptons"]==9
    assert ans["twisted"]["qqql_three_trivial_leptons"]==6
    e=np.eye(3)
    coeff=np.array([mixed(e[i],e[j],e[k]) for i,j,k in itertools.product(range(3),repeat=3)])
    assert np.max(abs(constraints@coeff))<1e-12
    assert np.max(abs(cube@coeff-coeff))<1e-12
    varied=[]
    for j in range(1,7):
        total=0.
        for ref in (False,True):
            for k in range(14):
                F=representation(k,ref)
                L=representation(k,ref,j)
                total+=np.trace(F)**3*np.trace(L)
        varied.append(int(round(total/28)))
    assert varied==[10,7,5,4,4,4]
    obs=json.loads((HERE.parent/"1381/inputs.json").read_text(encoding="utf-8-sig"))["ud_beta"]
    bound=[]
    for ud in (obs[0],obs[0]-2*obs[1]):
        s=math.sqrt(1-ud**2)
        u=np.array([1.,1.,0.])/math.sqrt(2)
        d=ud*u+s*e[2]
        null_u=e[1]
        null_d=ud*e[1]+s*e[2]
        assert math.isclose(float(mixed(u,u,d)),s/2,rel_tol=3e-15)
        assert abs(mixed(null_u,null_u,null_d))<1e-14
        bound.append({"ud":ud,"sharp_ordered_mixed_factor_upper":s/2,
                      "squared":s*s/4,"same_CKM_null":float(mixed(null_u,null_u,null_d))})
    ans["qqql_source_strength_examples"]=bound
    ans["authored_duql_when_lepton_doublet_varies"]=varied
    return ans


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    else:
        assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-13)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    path=HERE/"independent_results.json"
    out=run()
    if args.check:
        compare(out,json.loads(path.read_text()))
        print("1615 independent OK: group, Fierz quotient, mass-axis suppression and blind direction")
    else:
        with path.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
