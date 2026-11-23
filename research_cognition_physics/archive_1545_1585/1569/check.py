"""Topological endpoint phases and parity-resolved grouped trapping."""
import argparse
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent


def endpoint(m):
    # Phases in turns. At the end of the lifted color path each triplet gets m/3.
    color = F(m,3)
    su5 = [(color - F(m,3)) % 1] * 3 + [F(m)%1,F(0)]
    su4 = [(color+F(m,6))%1]*3 + [(-F(m,2))%1]
    su2 = [(F(m,2))%1,(-F(m,2))%1]
    assert all(x == 0 for x in su5)
    assert len(set(su4+su2)) == 1
    local = [(color+F(2*m,3))%1, (color-F(m,3))%1, F(-m)%1]
    assert not any(local)
    return {"PS_loop_class": int(2*su4[0]), "PS_smooth": su4[0]==0,
            "SU5_smooth_sector": True,"color_center_mod3": m%3}


def uniform_p0(n, lam, single):
    terms=[]
    for k in range(n+1):
        flip=single*k/n+(1-single)*2*k*(n-k)/(n*(n-1))
        terms.append(math.comb(n,k)*math.exp(-2*lam*flip))
    return math.fsum(terms)/(2**n)


def upper(n, fraction_single, miss, alpha):
    needed = (alpha-miss)/(1-miss)
    rank = n if fraction_single>0 else n-1
    if needed<=2**(-rank):
        return None
    lo,hi=0.,1.
    while uniform_p0(n,hi,fraction_single)>needed:
        hi*=2
    for _ in range(80):
        m=(lo+hi)/2
        if uniform_p0(n,m,fraction_single)>needed: lo=m
        else: hi=m
    return hi


def enumerate_small():
    # Poisson-thinned independent mark streams. Exact XOR convolution, independently of characters.
    n=4
    lam=2.3
    single=0.4
    marks={1<<i:single/n for i in range(n)}
    for i,j in itertools.combinations(range(n),2):
        marks[(1<<i)|(1<<j)] = (1-single)/math.comb(n,2)
    dist=[0.]*(1<<n)
    dist[0]=1.
    for h,p in marks.items():
        odd=-math.expm1(-2*lam*p)/2
        dist=[x*(1-odd)+dist[j^h]*odd for j,x in enumerate(dist)]
    assert abs(sum(dist)-1)<2e-14
    assert abs(dist[0]-uniform_p0(n,lam,single))<2e-14
    return dist[0]


def compute(cfg):
    n=cfg["tube_groups"]
    alpha=cfg["test_tail"]
    assert cfg["rings"]==n*cfg["rings_per_group"]
    models={}
    for name,single in cfg["normalized_visible_layouts"].items():
        row={}
        for miss in cfg["conditional_any_nonzero_parity_miss"]:
            cap=upper(n,single,miss,alpha)
            row[str(miss)]={"accepted_mean_cap":cap,
            "effective_visible_cross_section_nb":1000*cap/cfg["luminosity_inverse_microbarn"]}
        models[name]=row
    exact_single=-(n/2)*math.log(2*alpha**(1/n)-1)
    assert abs(models["single_uniform"]["0.0"]["accepted_mean_cap"]-exact_single)<2e-13
    return {"endpoint_classes_m_minus6_to6":{str(i):endpoint(i) for i in range(-6,7)},
            "accepted_event_limits":models,"small_non_character_crosscheck_p0":enumerate_small(),
            "fewest_independent_parity_modes_for_finite_95pct_cap":math.floor(-math.log2(alpha))+1,
            "single_unknown_concentration": {str(m):upper(m,1,0,alpha) for m in (2,3,4,5,10,57)},
            "always_same_tube_pair_p0":1,"fixed_two_distinct_tubes_high_rate_floor":0.5}


def rounded(x):
    if isinstance(x, dict):return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,float):return float(format(x,".12g"))
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--write-new",action="store_true")
    g.add_argument("--check",action="store_true")
    args=p.parse_args()
    cfg=json.loads((HERE/"inputs.json").read_text(encoding="utf-8-sig"))
    result=rounded(compute(cfg))
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    elif result!=json.loads(path.read_text(encoding="utf-8-sig")):
        raise SystemExit("Grouped magnetic-response certificate differs")
    print(result["accepted_event_limits"])
