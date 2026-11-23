"""A finite fourth spectator, two named preparation laws, and four-dimensional vacuum matching."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent
RAD=runpy.run_path(str(HERE/"../1751/check.py"))
REF=runpy.run_path(str(HERE/"../1748/check.py"))
R=4000.
ETA=.001
LAMBDA=1e-6


def pair_response(x,w):
    pairs=list(itertools.combinations(range(4),2))
    weights=np.array([w[i]*w[j] for i,j in pairs])
    prob=weights/weights.sum()
    out=np.array([sum(p*(u-x[i])*(u-x[j]) for p,(i,j) in zip(prob,pairs)) for u in x])
    return out,prob


def direct_trace(x,w):
    """Construct the six orthonormal antisymmetric states and evaluate the source."""
    pairs=list(itertools.combinations(range(4),2))
    vectors=np.zeros((16,len(pairs)))
    for k,(i,j) in enumerate(pairs):
        vectors[4*i+j,k]=1/math.sqrt(2)
        vectors[4*j+i,k]=-1/math.sqrt(2)
    density=(vectors*np.array([w[i]*w[j] for i,j in pairs]))@vectors.T
    density/=density.trace()
    X=np.diag(x)
    diff=np.diag(np.kron(X,np.eye(16))-np.kron(np.eye(4),np.kron(X,np.eye(4))))*np.diag(
        np.kron(X,np.eye(16))-np.kron(np.eye(16),X))
    return diff.reshape(4,16)@density.diagonal()


def radial_four(j,eta):
    states={}
    for ss in itertools.product((-1.,1.),repeat=4):
        s=np.array(ss)
        # delta around the two outer roots; avoids subtracting unit entries
        # to form the very small angular curvatures of opposite-sign modes.
        delta=np.zeros(4)
        for _ in range(6):
            delta-=(2*delta+3*s*delta**2+delta**3-eta*j)/(2+6*s*delta+3*delta**2)
        phi=s+delta
        assert np.max(abs(2*delta+3*s*delta**2+delta**3-eta*j))<1e-15
        rad=3*phi**2-1
        pair=[]
        for i,k in itertools.combinations(range(4),2):
            h=phi[i]**2+phi[k]**2+phi[i]*phi[k]-1 if s[i]==s[k] else eta*(j[i]-j[k])/(phi[i]-phi[k])
            pair.extend([h,h])
        h=np.r_[rad,pair]
        if min(h)<=0:
            continue
        tree=(eta*j/phi)**2/4-eta*j*phi
        cw=np.sum(h*h*(np.log(h/2)-1.5))/(64*math.pi**2)
        tag="".join("+" if t>0 else "-" for t in s)
        states[tag]=(tree,cw,phi,float(min(h)))
    assert len(states)==5
    tref,cref,phi,minimum=states["+-++"]
    delta={k:dict(tree=float(sum(t-tref)),CW=float(c-cref)) for k,(t,c,_,_) in states.items()}
    assert all(row["tree"]>=0 for row in delta.values())
    caps=[row["tree"]/(-row["CW"]) for row in delta.values() if row["CW"]<0]
    limit=min(caps)
    def winner(lam):
        return min(delta,key=lambda k:delta[k]["tree"]+lam*delta[k]["CW"])
    assert winner(limit/10)=="+-++"
    return dict(vacua_relative_to_alternating=delta,
                all_positive_relative_at_old_lambda=delta["++++"]["tree"]+LAMBDA*delta["++++"]["CW"],
                first_relative_matching_crossing_lambda=limit,winner_at_old_lambda=winner(LAMBDA),
                winner_at_one_tenth_crossing=winner(limit/10),smallest_tree_Hessian=minimum,
                alternating_tree_field=phi)


def data(h):
    a=(h[1]-h[0])/(h[2]-h[0])
    x=np.array([0.,a,1.,R])
    rows={}
    for p in (1,2):
        w=(1+x)**-p
        J,prob=pair_response(x,w)
        np.testing.assert_allclose(J,direct_trace(x,w),atol=2e-14)
        J3=pair_response(x,np.r_[w[:3],0.])[0][:3]
        frac=sum(v for v,(i,k) in zip(prob,itertools.combinations(range(4),2)) if k==3)
        rows[str(p)]=dict(raw_selection=J,relative_to_no_fourth_for_light=J[:3]/J3,
                           signs=np.sign(J),pair_contains_fourth=frac,
                           smallest_pair_probability=prob.min(),population_weight_condition=max(prob)/min(prob))
    acrit=a*(1-a)/(1-a-4*a)
    # t=(R-a)/(1+R)^2=a_crit. The larger root is the only crossing above R=1.
    rmin=(1-2*acrit+math.sqrt(1-4*acrit*(1+a)))/(2*acrit)
    independent_uniform_caps=[a*(1-a)/((R-a)*(1-2*a)),(1-a)/((R-1)*(2-a))]
    assert rows["1"]["signs"][1]>0 and np.all(rows["2"]["signs"]==[1,-1,1,1]) and R>rmin>1
    J=rows["2"]["raw_selection"]
    return dict(middle=a,exact_large_R_crossing_for_rate_preparation=rmin,
                uniform_light_weight_allowed_fourth=min(independent_uniform_caps),
                common_alignment_tests=rows,vacuum=radial_four(J/max(abs(J)),ETA))


def solve(step):
    bg=REF["OLD"]["mother"](step)
    Y=bg["yW"][2,:,:3]
    d=bg["yW"][1,:,:3]@bg["right_quark_to_GJ"]
    labs=[("quark",d),("lepton",Y)]
    sections={}
    for key,y in labs:
        h=np.linalg.svd(y,compute_uv=False)[::-1]**2
        sections[key]=data(h)
    sections["common_uniform_cap"]=min(sections[k]["uniform_light_weight_allowed_fourth"] for k,_ in labs)
    L=np.linalg.svd(Y)[0][:,::-1]
    X,_,_,_=REF["terms"](bg,0,(1,1,1))
    true=REF["terms"](bg,0,(1,-1,1))[0]
    # Fourth line is a separately protected/no-mixing direct summand of the light material port.
    lepton=sections["lepton"]
    runtime={}
    def response(sign):
        XS=X@((L*sign)@L.conj().T)
        vals,vec=np.linalg.eigh(np.eye(3)+XS.conj().T@XS)
        F=(vec/np.sqrt(vals))@vec.conj().T
        result=RAD["observables"](bg,F)
        return dict(charged_pulls=result["charged_pulls"],angles=result["angles"],clocks=result["solar_center_clocks"],
                    clock_pulls=result["clock_pulls"],quark_radius=result["quark_mass_radius"],
                    FC68_90=result["FC68_FC90"],passed=result["passed"],X_minus_1748_norm=np.linalg.norm(XS-true,2))
    for p,term in lepton["common_alignment_tests"].items():
        runtime["polar_p"+p]=response(term["signs"][:3])
    vac=lepton["vacuum"]
    runtime["alternating_tree_finite_4d"]=response(vac["alternating_tree_field"][:3])
    j=lepton["common_alignment_tests"]["2"]["raw_selection"].copy()
    j/=max(abs(j))
    winner=vac["winner_at_old_lambda"]
    shift=np.array([1. if s=="+" else -1. for s in winner])
    for k in range(5):
        shift-=(shift**3-shift-ETA*j)/(3*shift**2-1)
    runtime["old_lambda_lowest_of_tree_vacua"]=response(shift[:3])
    assert runtime["polar_p2"]["passed"] and not runtime["polar_p1"]["passed"]
    assert runtime["alternating_tree_finite_4d"]["passed"] and not runtime["old_lambda_lowest_of_tree_vacua"]["passed"]
    return sections,runtime


def clean(o):
    if isinstance(o,dict):return {k:clean(v) for k,v in o.items()}
    if isinstance(o,(np.ndarray,list,tuple)):return [clean(v) for v in o]
    if isinstance(o,(np.bool_,bool)):return bool(o)
    if isinstance(o,(float,np.floating)):return float(f"{o:.10g}")
    return o


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    cfg=p.parse_args()
    fine,flows=solve(.125)
    coarse,other=solve(.25)
    delta=max(abs(flows[k]["clocks"][0]-other[k]["clocks"][0]) for k in flows)
    assert delta<1e-9
    result=clean(dict(fourth_response_R=R,two_sources=fine,known_light_material=flows,step_halving_max_dm31=delta))
    path=HERE/"independent_results.json"
    if cfg.check:
        assert result==json.loads(path.read_bytes())
        print("1753 full-support antisymmetric source, suppression law, four-field vacuum and light matter checked")
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")

