"""Certify full-box monotonicity and construct a physical two-plane saturation."""

from argparse import ArgumentParser
from itertools import product
from math import prod, sqrt
from pathlib import Path
import json
import numpy as np

HERE = Path(__file__).resolve().parent


def terms(values):
    u1,u2,u3,d1,d2,d3=values
    g=prod(values)**(1/3)
    return g, [(g/(u2*d3))**2,(g/(u3*d2))**2,(u2/u3)**2,(d2/d3)**2]


def achieve(values):
    u1,u2,u3,d1,d2,d3=values
    g,(a,b,c,d)=terms(values)
    p=(u2*d2/g)**2*(1-a)*(1-b)/((1-c)*(1-d))
    assert 0 < p < 1 and u1*d1 < g
    V=np.eye(3)
    V[1:,1:]=[[sqrt(p),sqrt(1-p)],[-sqrt(1-p),sqrt(p)]]
    s=np.linalg.svd(np.diag(values[:3]) @ V @ np.diag(values[3:]),compute_uv=False)
    assert abs(s[1]/g-1)<1e-11
    return p


def results():
    archive=HERE.parent
    inputs=json.loads((archive/"1386"/"inputs.json").read_text(encoding="utf-8"))["yukawa"]
    data=np.array(inputs["up"]+inputs["down"])
    center,error=data.T
    # All four correction factors are exact monomials; extrema at box corners.
    all_z=[terms(center+3*np.array(sign)*error)[1] for sign in product((-1,1),repeat=6)]
    zhi=np.max(all_z,axis=0)
    deriv=2*sum(zhi/(1-zhi))
    assert max(zhi)<1 and deriv < 2/3
    signs=np.array([-1,1,-1,-1,1,-1])
    def best(width):
        return achieve(center+width*signs*error)
    cms=0.955
    out={"three_width_bound_on_log_correction_gradient":round(deriv,10),
         "sharp_upper_vs_width":[round(best(w),10) for w in range(4)]}
    for r in (1.,1.05,1.2):
        required=cms/(cms+r*(1-cms))
        left,right=2.,3.
        assert best(left)<required<best(right)
        for _ in range(50):
            mid=(left+right)/2
            if best(mid)<required:
                left=mid
            else:
                right=mid
        out["critical_width_relative_weight_"+str(r)]=round(right,9)
    p=best(2)
    weight=2.4
    R=weight*p/(1-p+weight*p)
    assert R>cms
    out["counterexample"]={"p_exact_shape":round(p,10),"relative_weight":weight,"R":round(R,10)}
    return out


if __name__=="__main__":
    parser=ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    obj=results()
    target=HERE/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf8") as file:
            json.dump(obj,file,indent=2)
            file.write("\n")
    else:
        assert obj==json.loads(target.read_text(encoding="utf8"))
    print(json.dumps(obj))
