"""Direct charged and neutral vacuum-polarization sum, with U(3) outer bounds."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
SCALE=3/(16*math.pi**2*246.22**2)
UP=[.0022,1.27,172.5]


def load(n,fn):
    return json.loads((HERE.parent/str(n)/fn).read_text(encoding="utf-8-sig"))


def theta(m,n):
    if n==m: return 0.,0.
    if m*n==0: return m*m+n*n,0.
    x,y=m*m,n*n
    t=(x-y)/(x+y)
    if abs(t)<.1:
        a=(x+y)*sum(2*t**(2*k)/((2*k-1)*(2*k+1)) for k in range(1,12))
        b=4*m*n*sum(t**(2*k)/(2*k+1) for k in range(1,12))
    else:
        # Logarithms of true masses avoid atanh rounding for a large hierarchy.
        l=math.log(m/n)
        a=x+y-4*x*y/(x-y)*l
        b=4*m*n*(l/t-1)
    return a,b


def evaluate(L,R,m,n):
    p,q=theta(m,n)
    return (abs(L)**2+abs(R)**2)*p+2*(L.conjugate()*R).real*q


def polar(a):
    u,s,vh=np.linalg.svd(a)
    return u@vh


def ref_ckm():
    p=load(1635,"inputs.json")["witness"]["ordinary_rotation_s12_s23_s13_delta"]
    s12,s23,s13,delta=p
    c12,c23,c13=np.sqrt(1-np.array(p[:3])**2)
    return np.array([[1,0,0],[0,c23,s23],[0,-s23,c23]]) @ np.array(
        [[c13,0,s13*np.exp(-1j*delta)],[0,1,0],[-s13*np.exp(1j*delta),0,c13]]) @ np.array(
        [[c12,s12,0],[-s12,c12,0],[0,0,1]])


def unitary_bounds(matrices):
    # A common Hermitian summand contributes its trace for every unitary.
    # The remaining row sum is bounded by two independent variational inequalities.
    best=-math.inf
    perms=list(itertools.permutations(range(3)))
    trial=min(sum(h[p[i],p[i]].real for i,h in enumerate(matrices)) for p in perms)
    for baseline in matrices:
        delta=[(h-baseline+(h-baseline).conj().T)/2 for h in matrices]
        trace=baseline.trace().real
        row_bound=trace+sum(np.linalg.eigvalsh(h)[0] for h in delta)
        off=[h-np.diag(np.diag(h)) for h in delta]
        diagonal=min(sum(h[p[i],p[i]].real for i,h in enumerate(delta)) for p in perms)
        doubly_stochastic=trace+diagonal+sum(np.linalg.eigvalsh(h)[0] for h in off)
        best=max(best,row_bound,doubly_stochastic)
    # Explicit outward guard in GeV^2, much larger than binary64 scaling times epsilon.
    return [best-1e-5,trial+1e-5]


def result(raw,md):
    arr=np.asarray(raw)
    mat=arr[:,:,0]+1j*arr[:,:,1]
    left,m,vh=np.linalg.svd(mat)
    order=np.argsort(m)
    m=m[order];left=left[:,order];right=vh.conj().T[:,order]
    K=left[:3]; L=left[4]; R=right[4]
    # 2T3 for the -1/3 sector. Identity terms have vanishing off-diagonal entries.
    NL=-K.conj().T@K+np.outer(L.conjugate(),L)
    NR=np.outer(R.conjugate(),R)
    exotic=sum(evaluate(l,r,z,md) for l,r,z in zip(L,R,m))
    neutral=sum(evaluate(NL[i,j],NR[i,j],m[i],m[j]) for i,j in itertools.combinations(range(5),2))
    if abs(NL[1,2])>0: assert (NL[1,2]/NR[1,2]).real>0
    v0=polar(K[:,:3]); Y=v0.conj().T@K
    assert np.max(abs(Y[:,:3]-Y[:,:3].conj().T))<3e-12
    hats=[np.array([theta(up,d)[0] for d in m]) for up in UP]
    h=[(Y*row)@Y.conj().T for row in hats]
    hSM=[np.diag(row[:3]) for row in hats]
    ordinary=sum(np.sum(abs(K[i])**2*hats[i]) for i in range(3))
    ordinary_SM=sum(np.sum(abs(v0[i])**2*hats[i][:3]) for i in range(3))
    fixed_ref=sum(np.sum(abs(ref_ckm()[i])**2*hats[i][:3]) for i in range(3))
    delta=ordinary-ordinary_SM+exotic-neutral
    same=np.array(unitary_bounds([x-y for x,y in zip(h,hSM)]))+exotic-neutral
    absolute=np.array(unitary_bounds(h))+exotic-neutral-fixed_ref
    return {"M_GeV":md,"masses":m.tolist(),"ordinary_W_GeV2":ordinary,"exotic_W_GeV2":exotic,
            "neutral_Z_GeV2":neutral,"SM_W_same_GeV2":ordinary_SM,"SM_fixed_ref_GeV2":fixed_ref,
            "Delta_rho_at_ordinary_rotation":SCALE*delta,
            "Delta_rho_at_fixed_ref":SCALE*(ordinary+exotic-neutral-fixed_ref),
            "T":127.95*SCALE*delta, "inf_Delta_rho_same_U_bracket":(same*SCALE).tolist(),
            "inf_Delta_rho_fixed_Uref_bracket":(absolute*SCALE).tolist()}


def electroweak(row):
    c2=.7688;alpha=1/127.95
    reference=[80.3545,.0059,.231534,.000062]
    m,w=80.3602,.0099
    residual={}
    for n in (1,2):
        m0,e0,s0,ds=reference
        xmin=((m-n*w)/(m0+n*e0))**2-1+2*((row[0]-n*row[1])-(s0+n*ds))
        xmax=((m+n*w)/(m0-n*e0))**2-1+2*((row[0]+n*row[1])-(s0-n*ds))
        residual[str(n)]={"X":[xmin,xmax],"Delta_rho_if_U_and_vertex_zero":[xmin/c2,xmax/c2]}
    return residual


def rounded(x):
    if isinstance(x,dict):return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [rounded(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.10g}")
    return x


def all_data():
    source=load(1635,"results.json")
    out={name:result(source[name]["weak_basis_mass_GeV"],M)
        for name,M in (("common_witness",1600),("high_mass_witness",2600))}
    solved=load(1640,"independent_results.json")
    out["retuned_tree_solutions"]=[result(row["weak_basis_mass_re_im"],row["M"]) for row in solved["witnesses"]]
    out["primary_1640_tree_solutions"]=[result(row["weak_basis_mass_re_im"],row["M"]) for row in solved["primary_calibration_SVD"]]
    for label, val in (("published_CMS",[.23152,.00031]),("preliminary_CMS",[.23157,.00031])):
        ranges=electroweak(val)
        for n in (1,2):
            ent=ranges[str(n)]
            for key in ("common_witness","high_mass_witness"):
                rho=out[key]["inf_Delta_rho_fixed_Uref_bracket"][0]
                mismatch=max(0,.7688*rho-ent["X"][1])
                ent[key+"_minimal_extra_abs_X"]=mismatch
                ent[key+"_required_max_U_when_no_other_residual"]=(ent["X"][1]-.7688*rho)*(4*(1-.7688))*127.95
        out[label]=ranges
    return rounded(out)


if __name__=="__main__":
    cmd=argparse.ArgumentParser();cmd.add_argument("--check",action="store_true");args=cmd.parse_args()
    target=HERE/"independent_results.json";data=all_data()
    if args.check:
        assert json.loads(target.read_text(encoding="utf-8"))==data
        print("1642 independent neutral/charged sum, two variational bounds, two CMS inputs: PASS")
    else:
        target.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
        for k in ("common_witness","high_mass_witness","published_CMS"):print(k,data[k])
