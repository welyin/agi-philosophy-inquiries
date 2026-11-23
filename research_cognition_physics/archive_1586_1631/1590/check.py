"""Eliminate one common radiative muon mass; cross-check integral and published functions."""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent


def load(name):
    return json.loads((ROOT/name).read_text(encoding="utf-8-sig"))


def exact_functions(fermion,m1,m2):
    with localcontext() as ctx:
        ctx.prec=60
        mass=Decimal(fermion)
        x1,x2=[(Decimal(m)/mass)**2 for m in (m1,m2)]
        def mass_fn(x):
            return Decimal(1) if x==1 else x*x.ln()/(x-1)
        def dipole_fn(x,Q):
            if x==1:
                return -Decimal(2)/3+Q
            d=1-x
            return ((3*x-1)/d**2+2*x*x*x.ln()/d**3+
                    2*Q*(1/d+x*x.ln()/d**2))
        diff=mass_fn(x2)-mass_fn(x1)
        return float(diff), [float((dipole_fn(x1,q)-dipole_fn(x2,q))/diff/mass**2)
                             for q in (0,1)]


def integrate(fermion,m1,m2):
    # t=1-exp(-q) exposes small scalar masses without a near-endpoint mesh loss.
    step=45/200000
    q=(np.arange(200000)+.5)*step
    one_t=np.exp(-q)
    t=-np.expm1(-q)
    delta1=one_t*fermion**2+t*m1**2
    delta2=one_t*fermion**2+t*m2**2
    dmass=m2**2-m1**2
    denominator=np.sum(one_t*np.log1p(t*dmass/delta1)/dmass)
    return [float(np.sum(2*one_t*t*one_t*(charge-one_t)/(delta1*delta2))/denominator)
            for charge in (0,1)]


def coefficient_checks():
    # Coefficient-level audit of the two expansions; the proof covers all orders.
    for sector in ("small_scalar","large_scalar"):
        old_pos=old_neg=Fraction(0)
        for n in range(200):
            if sector=="small_scalar":
                den=Fraction(1,n+2)
                pos=Fraction(2*(n+1),(n+3)*(n+4))
                neg=Fraction(4*(n+1),(n+2)*(n+3)*(n+4))
            else:
                den=Fraction(1,(n+1)*(n+2))
                pos=Fraction(4*(n+1),(n+2)*(n+3)*(n+4))
                neg=Fraction(2*(n+1),(n+3)*(n+4))
            assert pos>=den/3 and pos-old_pos<=den/3 and neg-old_neg<=den/3
            old_pos,old_neg=pos,neg


def calculate():
    cfg=load("inputs.json")
    m=cfg["muon_mass_GeV"]
    scale=cfg["scale_in_a_mu"]
    cases=[]
    for entry in cfg["benchmarks"]:
        fermion,m1,m2=entry["masses_GeV"]
        dmass,coeff=exact_functions(fermion,m1,m2)
        integral=integrate(fermion,m1,m2)
        assert all(math.isclose(x,y,rel_tol=3e-7) for x,y in zip(coeff,integral))
        low=1/(3*max(fermion,m2)**2)
        high=1/(3*min(fermion,m1)**2)
        assert low<=coeff[1]<=high and 0<-coeff[0]<=high
        product=32*math.pi**2*m/(fermion*entry["sin_2theta"]*dmass)
        cases.append({"label":entry["label"],"Q0_and_Q1_dipoles_1e11":[m*m*k/scale for k in coeff],
                      "absolute_Q0_times_Mmax_squared":abs(coeff[0])*max(fermion,m2)**2,
                      "absolute_Yukawa_product_needed":product,
                      "mu_mix_GeV":entry["sin_2theta"]*(m2**2-m1**2)/(math.sqrt(2)*cfg["v_GeV"]),
                      "degenerate_positive_kernel_floor_GeVm2":low})
    assert cases[-1]["absolute_Q0_times_Mmax_squared"]<1/3
    a,ea=cfg["experiment_world"]
    s,es=cfg["WP20_theory"]
    box={}
    wp20=[a-s, math.hypot(ea,es)]
    for name,(delta,error) in (("WP25",cfg["WP25_difference_as_published"]),("historical_WP20",wp20)):
        lo=delta-cfg["widths"]*error
        hi=delta+cfg["widths"]*error
        tests={}
        for eps in cfg["relative_response_cases"]:
            item={"positive_largest_mass_min_GeV":m*math.sqrt((1-eps)/(3*hi*scale))}
            if lo>0:
                item["positive_smallest_mass_max_GeV"]=m*math.sqrt((1+eps)/(3*lo*scale))
                item["negative_possible_with_sign_preserved"]=False
            else:
                r=cfg["bounded_ratio_example"]
                item["negative_largest_mass_min_if_compressed_GeV"]=m*math.sqrt((1-eps)/(3*(-lo)*scale))/r
            tests[str(eps)]=item
        box[name]={"difference_and_error":[delta,error],"working_interval":[lo,hi],"mass_tests":tests}
    coefficient_checks()
    return {"benchmarks":cases,"theory_branches":box}


def same(a,b):
    if isinstance(a,dict):
        return set(a)==set(b) and all(same(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    if isinstance(a,str) or isinstance(a,bool):
        return a==b
    return math.isclose(a,b,rel_tol=5e-11,abs_tol=1e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    output=calculate()
    if args.check:
        if not same(output,load("results.json")):
            raise SystemExit("1590 results changed")
        print("1590: mass elimination, dipole integral and spectral bounds verified")
    else:
        (ROOT/"results.json").write_text(json.dumps(output,indent=2)+"\n",encoding="utf-8")
