"""Rational all-scale envelopes when the sterile singular values track the top."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
KLO=F("78.95683520")
KHI=F("78.95683522")
ROUND=10**13


def load(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"),parse_float=F)


def ceil_grid(x):
    return F(-((-x.numerator*ROUND)//x.denominator),ROUND)


def bound_mass(row):
    return row[0]+2*row[1]


def lower_gauge_squared(value,b,left,right):
    # All input values and interval endpoints are rational.
    t=left if b>0 else right
    K=KHI if b>0 else KLO
    return 1/(1/value**2-b*t/K)


def envelope(r):
    ys=load("../1386/inputs.json")["yukawa"]
    gs=load("../1369/inputs.json")
    gauge=[gs[k][0]-2*gs[k][1] for k in ("g1_GUT","g2","g3")]
    bs=[F(41,10),F(-19,6),F(-7)]
    U=bound_mass(ys["up"][-1])**2
    S=3*sum(bound_mass(row)**2 for row in ys["up"][:2]+ys["down"])
    S+=sum(bound_mass(row)**2 for row in ys["charged_leptons"])
    n=16384
    T=F("39.44")
    dt=T/n
    p=F("1.002")
    maximum=U
    scalar_cap=S
    min_gap=F(100)
    for j in range(n):
        g=[lower_gauge_squared(v,b,j*dt,(j+1)*dt) for v,b in zip(gauge,bs)]
        G=F(17,20)*g[0]+F(9,4)*g[1]+8*g[2]
        min_gap=min(min_gap,F(2,5)*g[0]+8*g[2])
        ub,sb=U*p,S*p
        top_driver=(F(9,2)+r)*ub+sb-G
        xu=dt*top_driver/(KLO if top_driver>=0 else KHI)
        xs=dt*((3+r)*ub+F(5,2)*sb)/KLO
        assert max(xu,xs)<1
        new_u=ceil_grid(U/(1-xu))
        new_s=ceil_grid(S/(1-xs))
        assert max(U,new_u)<ub and new_s<sb
        maximum=max(maximum,U,new_u)
        scalar_cap=max(scalar_cap,new_s)
        U,S=new_u,new_s
    assert maximum<4*F("3.14159") and min_gap>F(3,2)
    return {"r":str(r),"intervals":n,"max_top_yukawa":math.sqrt(float(maximum)),
            "max_other_weighted_trace":float(scalar_cap),"color_gauge_gap_floor":float(min_gap)}


def run():
    assert math.log(1.22e19/91.1876)<39.44
    return {"rank_bound":3,"window_log_ceiling_rational":"39.44",
            "integer_barrier_results":[envelope(F(3)),envelope(F(16,5))],
            "necessary_max_neutral_over_top_somewhere":math.sqrt(16/15)}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    out=HERE/"independent_results.json"
    result=run()
    if args.write_new:
        with out.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2)
            handle.write("\n")
    else:
        assert result==json.loads(out.read_text())
    print(json.dumps(result))
