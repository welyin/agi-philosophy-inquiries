"""Phase-independent off-diagonal test for the triplet-dominant SO(10) mass sum rule."""
import argparse
from itertools import combinations, product
import json
from math import sqrt
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))

def box(pair):
    return pair[0]-2*pair[1],pair[0]+2*pair[-1]

def tau_enclosure(a,b,c):
    A,B,C = a[1],b[1],c[1]
    a,b,c = a[0],b[0],c[0]
    p1=(max(0.,sqrt(a*c)-sqrt((1-a)*(1-c)*B))**2,
        min(1.,(sqrt(A*C)+sqrt((1-a)*(1-c)*B))**2))
    p2=(max(0.,sqrt((1-A)*c)-sqrt(A*(1-c)*B))**2,
        min(1.,(sqrt((1-a)*C)+sqrt(A*(1-c)*B))**2))
    return [p1,p2,((1-C)*(1-B),(1-c)*(1-b))]

def max_purity(bounds):
    # Convex maximum on the intersection of the probability simplex with a box.
    values=[]
    for i,j in combinations(range(3),2):
        k=3-i-j
        for u,v in product(bounds[i],bounds[j]):
            w=1-u-v
            if bounds[k][0]-1.e-13 <= w <= bounds[k][1]+1.e-13:
                values.append(u*u+v*v+w*w)
    assert values
    return max(values)

def subboxes(interval,count):
    cuts=np.linspace(*interval,count+1)
    return list(zip(cuts[:-1],cuts[1:]))

def solve():
    cfg=read("inputs.json")
    quarks=read(cfg["high_quarks"])["high_scale"]["down"]
    leptons=read(cfg["high_leptons"])["scales"]["1e16"][3:]
    osc=read(cfg["oscillation"])
    juno=read(cfg["juno"])["juno_207_days"]
    d,s,b=[box(v)[1] for v in quarks]
    e,mu=[box(v)[1] for v in leptons[:2]]
    tau=box(leptons[2])[0]
    dnorm=sqrt(d*d+s*s+b*b)
    out={}
    reactor=tuple((1-sqrt(1-v))/2 for v in box(osc["DayaBay_amplitude"]))
    for order in ["NO","IO"]:
        solar=box(cfg["NO_solar_sin2"] if order=="NO" else osc[order]["solar_sine_squared"])
        atm=box(osc[order]["atmospheric_sine_squared"])
        base=tau_enclosure(solar,reactor,atm)
        rough=max_purity(base)
        best=max(max_purity(tau_enclosure(aa,bb,cc))
                 for aa,bb,cc in product(subboxes(solar,8),subboxes(reactor,8),subboxes(atm,8)))
        floor=(tau-mu)*sqrt(1-best)
        rough_floor=(tau-mu)*sqrt(1-rough)
        assert rough_floor > dnorm and floor>=rough_floor-1.e-13
        # Second certificate: min singular value of D-E never exceeds s+mu.
        mS=box(juno["solar_21"])
        mA=box(juno["normal_atmospheric_31" if order=="NO" else "inverted_abs_31"])
        if order=="IO":
            ztau=sqrt(base[2][0])
            ze=sqrt(reactor[1])
            shape=sqrt((1-reactor[1])*mA[0]/(mA[0]+mS[1]))
        else:
            ztau=sqrt(base[0][0])
            ze=sqrt((1-solar[0])*(1-reactor[0]))
            shape=sqrt(reactor[0]+mS[0]/mA[1]*solar[0]*(1-reactor[1]))
        def wedge_gap(err):
            min_w=ztau*(tau-s-mu-err)-err
            if min_w <= 0:
                return -np.inf
            row_max=b*(s+ze*(e+s+mu+err)+err)/min_w
            return (tau-b-err)*shape-e-err-row_max
        left,right=0.,tau-b
        if wedge_gap(0)>0:
            for _ in range(80):
                mid=(left+right)/2
                if wedge_gap(mid)>0:
                    left=mid
                else:
                    right=mid
        out[order]={"tau_probability_rectangle":base,"max_purity_coarse":rough,
            "max_purity_512_cells":float(best),"coarse_F_norm_floor":rough_floor,
            "matched_D_F_norm_upper":dnorm,"F_norm_floor":float(floor),
            "min_additive_F_error":float(floor-dnorm),
            "min_fraction_of_tau_mu_gap_additive_F":float((floor-dnorm)/(tau-mu)),
            "necessary_bottom_min_given_other_down_high":float(sqrt(floor**2-d*d-s*s)),
            "electron_wedge_operator_error_at_least":float(left)}
    # Check all covariant identities on a complex unitary, non-aligned exact mass pair.
    rng=np.random.default_rng(1453)
    u=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    E=np.diag([.001,.1,1.])
    n=u.conj()@np.diag([.02,.07,.21])@u.conj().T
    D=E+(.2+.7j)*n
    L=u.T@E@u
    transformed=u.T@D@u
    assert np.linalg.norm((transformed-L)-np.diag(np.diag(transformed-L)))<1.e-14
    off=L-np.diag(np.diag(L))
    p=abs(u)**2
    direct=np.linalg.norm(E)**2-np.sum((np.diag(E)@p)**2)
    variance=sum((E[i,i]-E[j,j])**2*(p[i]@p[j]) for i,j in combinations(range(3),2))
    assert np.isclose(direct,variance) and np.linalg.norm(off)**2>=variance-1.e-14
    return out

def quantize(obj):
    if isinstance(obj,dict):
        return {k:quantize(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [quantize(v) for v in obj]
    return round(float(obj),12)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    flags=p.add_mutually_exclusive_group()
    flags.add_argument("--write-new",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=quantize(solve())
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    elif args.check:
        assert result==read("results.json")
        print("1453 tau off-diagonal and electron exterior-power certificates verified")
    else:
        print(json.dumps(result,indent=2))
