"""Four nucleon searches with an arbitrary complex one-direction loss."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent
ORDER=("p_pi","p_K","n_pi","n_K")
REF=1e15


def load(path):
    return json.loads(path.read_text(encoding="utf8"))


def input_rows():
    now=load(BASE/"inputs.json")
    old=load(BASE/now["previous_nucleons"])
    t12,t13=load(BASE/now["angles"])["high_scale"]["angles_theta12_theta13_theta23_delta"][:2]
    ud=math.cos(t12[0])*math.cos(t13[0])
    us=math.sin(t12[0])*math.cos(t13[0])
    pion=now["pion_W0_GeV2"]
    us_d,ud_s,ds_u=(old["W0_GeV2"][k] for k in ("us_d","ud_s","ds_u"))
    masses=old["mass_GeV"]
    def factor(m,mes,tau):
        width=m*(1-(mes/m)**2)**2/(32*math.pi)
        return math.sqrt(width*tau*now["year_seconds"]/now["hbar_GeV_s"])/REF**2
    # S-K 2510.26232, both single-pion lifetime limits given in the same paper.
    ks=np.array([factor(masses["p"],now["pion_mass_GeV"],now["pion_lifetime_years"]),
                 factor(masses["p"],masses["K+"],old["proton"]["partial_lifetime_lower_yr"]),
                 factor(masses["n"],.1349768,1.4e33)/math.sqrt(2),
                 factor(masses["n"],masses["K0"],old["neutron"]["partial_lifetime_lower_yr"])])
    rows=np.array([[ud*pion[0],0],[us*ud_s[0],ud*us_d[0]],
                   [ud*pion[0],0],[-us*ud_s[0],ud*ds_u[0]]])*ks[:,None]
    def box(w):
        t=t12[0]+w*t12[1]
        c13min=math.cos(t13[0]+w*t13[1])
        c13max=math.cos(t13[0]-w*t13[1])
        low=lambda data:abs(data[0])-w*sum(data[1:])
        high=lambda data:abs(data[0])+w*sum(data[1:])
        A=ks[2]*math.cos(t)*c13min*low(pion)
        B=ks[1]*math.sin(t)*c13max*high(ud_s)
        C=ks[1]*math.cos(t)*c13min*low(us_d)
        assert min(A,B,C)>0
        return REF*math.sqrt(A*C/math.hypot(A+B,C))
    return rows,{str(w):box(w) for w in (0,1,2)}


def optimum(rows,z):
    x,y=rows.T
    base=(1+z*z)*(x*x+y*y)/2
    dirs=np.column_stack([x*x-y*y,2*x*y])*(1-z*z)/2
    candidates=[np.zeros(2)]
    for g in dirs:
        if np.linalg.norm(g)>1e-14:
            candidates.append(g/np.linalg.norm(g))
    for i,j in itertools.combinations(range(len(rows)),2):
        v=dirs[i]-dirs[j]
        length=np.linalg.norm(v)
        if length<1e-14:
            continue
        unit=v/length
        k=(base[i]-base[j])/length
        if abs(k)<=1+1e-13:
            b=math.sqrt(max(0,1-k*k))*np.array([-unit[1],unit[0]])
            candidates.extend([k*unit+b,k*unit-b])
    for i,j,k in itertools.combinations(range(len(rows)),3):
        v=np.stack([dirs[i]-dirs[j],dirs[i]-dirs[k]])
        if abs(np.linalg.det(v))>1e-13:
            t=np.linalg.solve(v,[base[i]-base[j],base[i]-base[k]])
            if np.linalg.norm(t)<=1+1e-13:
                candidates.append(t)
    q,p=min(((float(max(base-dirs@u)),u) for u in candidates),key=lambda x:x[0])
    # A density on a real plane always equals Re(nn†) for a normalized complex n.
    r=np.array([[(1+p[0])/2,p[1]/2],[p[1]/2,(1-p[0])/2]])
    ev,vec=np.linalg.eigh(r)
    n=vec@np.array([math.sqrt(max(ev[0],0)),1j*math.sqrt(max(ev[1],0))])
    g=np.eye(2)-(1-z*z)*np.outer(n,n.conj())
    rates=np.einsum("ij,jk,ik->i",rows,g,rows).real
    np.testing.assert_allclose(rates,base-dirs@p,rtol=2e-12,atol=1e-12)
    assert np.linalg.norm(p)<=1+2e-13
    return {"scale_GeV":REF*q**.25,"max_abs_Wilson_GeV_minus2":1e-30/math.sqrt(q),
            "rates_for_C_1e_minus30":rates.tolist(),"minimax":q,"bloch":p.tolist()}


def audit():
    rows,box=input_rows()
    out={}
    for z in (0.,.05923,.5,1.):
        p=optimum(rows[:2],z)
        all_four=optimum(rows,z)
        assert all_four["scale_GeV"]>=p["scale_GeV"]
        out[str(z)]={"proton_only":p,"four_searches":all_four,
                     "scale_improvement":all_four["scale_GeV"]/p["scale_GeV"]}
    np.testing.assert_allclose(out["0.0"]["four_searches"]["scale_GeV"],box["0"],rtol=3e-14)
    return {"channel_order":list(ORDER),"rows_for_C_1e_minus30":rows.tolist(),
            "isoscalar_pion_response_weight_ratio_n_to_p":float((rows[2,0]/rows[0,0])**2),
            "loss_floor_examples":out,"four_channel_robust_scale_lower_GeV":box}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
