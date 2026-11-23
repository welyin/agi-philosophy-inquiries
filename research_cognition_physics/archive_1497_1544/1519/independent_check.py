"""Boundary-current exchange and threshold matching in a finite interval tower."""
import argparse
import json
import math
from pathlib import Path

OUT = Path(__file__).with_name("independent_results.json")


def mode(j,rho):
    x0=(j+.5)*math.pi
    if rho==0:
        return x0,2.
    # tan(x-x0)=rho/x avoids cancellation at a pole of cot(x).
    lo,hi=0.,math.pi/2
    for _ in range(60):
        mid=(lo+hi)/2
        if (x0+mid)*math.tan(mid) < rho:
            lo=mid
        else:
            hi=mid
    shift=(lo+hi)/2
    x=x0+shift
    boundary_norm=math.cos(shift)**2/(.5+math.sin(2*shift)/(4*x))
    return x,boundary_norm


def calculate():
    out={}
    for rho in (0,1,4,9):
        full=1/(1+rho)
        out[str(rho)]={}
        for n in (2,8,32):
            modes=[mode(j,rho) for j in range(n)]
            green=sum(w/(x*x) for x,w in modes)
            logdet=2*sum(math.log(x/((j+.5)*math.pi)) for j,(x,w) in enumerate(modes))
            tail_bound=2/(math.pi**2*n)
            assert -1e-13<=full-green<=tail_bound+1e-13
            assert -1e-13<=math.log1p(rho)-logdet<=rho*tail_bound+1e-13
            if rho>0:
                h=1e-5
                deriv=sum((2*math.log(mode(j,rho+h)[0]/mode(j,rho-h)[0]))/(2*h) for j in range(n))
                assert abs(deriv-green)<2e-8
            out[str(rho)][str(n)]={
                "dimensionless_exchange_sum":green,
                "exchange_omission_upper":tail_bound,
                "log_squared_mass_product":logdet,
                "untruncated_exchange":full,
                "untruncated_log_squared_product":math.log1p(rho)
            }
    return out


def verify(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            verify(a[k],b[k])
    else:
        assert math.isclose(a,b,rel_tol=2e-12,abs_tol=1e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    m=parser.add_mutually_exclusive_group(required=True)
    m.add_argument("--check",action="store_true")
    m.add_argument("--write-new",action="store_true")
    opt=parser.parse_args()
    result=calculate()
    if opt.write_new:
        with OUT.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        verify(result,json.loads(OUT.read_text(encoding="utf8")))
    print("Verified Robin spectral thresholds, exchange, and finite tower tail bounds.")
