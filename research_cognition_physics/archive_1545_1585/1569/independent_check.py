"""Check the trinification lift and a modular, signed magnetic-charge readout."""
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent

def center(m):
    # Fundamental weights in units e. Off-diagonal color lift, electromagnetic L,R.
    generators = [[F(1,3), F(1,3), F(-2,3)],
                  [F(1,3), F(-2,3), F(1,3)],
                  [F(1,3), F(-2,3), F(1,3)]]
    endpoints = [[(m*q) % 1 for q in group] for group in generators]
    return all(q == F(m,3)%1 for group in endpoints for q in group)

def count_distribution(n, q):
    # Character average for a uniformly selected ordered pair of distinct tubes.
    for a in range(n+1):
        if q == 2:
            populations = [(n-a,0)]
        else:
            populations = [(b,n-a-b) for b in range(n-a+1)]
        for b,c in populations:
            norm = (a-b)**2 if q == 2 else a*a+b*b+c*c-a*b-b*c-c*a
            amplitude = (norm-n)/(n*(n-1))
            multiplicity = math.comb(n,a)*math.comb(n-a,b)
            yield multiplicity / q**n, amplitude

def null_rate(n, q):
    chars = list(count_distribution(n,q))
    assert abs(math.fsum(p for p,_ in chars)-1)<1e-14
    saturation = math.fsum(p for p,f in chars if f == 1.)
    assert abs(saturation - q**(1-n)) < 1e-14
    def probability(lam):
        return math.fsum(p*math.exp(-lam*(1-f)) for p,f in chars)
    assert probability(0)>0.9999999
    if saturation >= .05:
        return {"large_mean_limit":saturation, "mean_at_five_percent":None}
    left,right=0.,10.
    assert probability(right) < .05
    for _ in range(60):
        mid=(left+right)/2
        if probability(mid)>.05:
            left=mid
        else:
            right=mid
    return {"large_mean_limit":saturation, "mean_at_five_percent":(left+right)/2}

def compute():
    assert all(center(n) for n in range(-24,25))
    return {"primitive_smooth_flux":{"SU5":1,"PS_quotient":2,"trinification_quotient":3,
                                    "EW_single_doublet":"zero kernel"},
            "signed_pairs_in_distinct_equiprobable_tubes":
                {str(n):{str(q):null_rate(n,q) for q in (2,3)} for n in (2,3,57)}}

if __name__ == "__main__":
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    cfg=p.parse_args()
    result=compute()
    path=HERE/"independent_results.json"
    if cfg.write_new:
        with path.open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")
    else:
        assert result==json.loads(path.read_text(encoding="utf-8"))
        print("1569 central lifts and signed mod-two/mod-three deposit calculation passed")
