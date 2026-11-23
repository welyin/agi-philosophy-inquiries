"""Phase-resolved amplitude, and a distinct CMS raw count for equal moduli."""
import argparse
import cmath
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())


def poisson_upper(n):
    def cdf(mu):
        return math.fsum(math.exp(-mu+k*math.log(mu)-math.lgamma(k+1)) for k in range(n+1))
    lo, hi=float(n),2*n+25.
    for unused in range(80):
        mid=(lo+hi)/2
        if cdf(mid)>.05:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def minimum(phase,p,q,plus,minus):
    e=cmath.exp(1j*phase)
    a=(p+q*e)/2
    b=(p-q*e)/2
    # Independent extreme positive Gram matrix by a single complex channel.
    u=-plus*(b/abs(b))*(a.conjugate()/abs(a))
    v=minus
    w=(u+v)/2
    z=(u-v)/2
    actual=abs(p*w+q*e*z)**2
    floor=(abs(b)*minus-abs(a)*plus)**2
    assert math.isclose(floor,actual,abs_tol=1e-10)
    return floor


def compute():
    atlas=CFG["VBF_WH_same_SRminus"]
    ua=poisson_upper(atlas["raw"])
    uc=poisson_upper(130)
    ans={"ATLAS_raw_total_upper":ua,"CMS_raw_total_upper":uc,"boxes":{}}
    for j in (0,1,2):
        p=math.sqrt(.95-j*.20-.02)
        q=math.sqrt(.88-j*.20-.10)
        plus=math.sqrt(2.93+j*.35)
        minus=math.sqrt(285-j*45)
        lo,hi=math.pi/2,math.pi
        if minimum(lo,p,q,plus,minus)>ua:
            boundary=lo
        else:
            for unused in range(75):
                mid=(lo+hi)/2
                if minimum(mid,p,q,plus,minus)>ua:
                    hi=mid
                else:
                    lo=mid
            boundary=(lo+hi)/2
        floor=minimum(math.pi,p,q,plus,minus)
        cms=366-j*(3+68)
        ans["boxes"][str(j)]={
          "ATLAS_opposite_floor_count":floor,
          "ATLAS_negative_half_phase_halfwidth_about_pi_degrees":180-math.degrees(boundary),
          "CMS_equal_moduli_threshold":math.sqrt(uc/cms),
          "CMS_equal_moduli_signal_floor_at_associated_minimum":cms*p*p,
          "CMS_equal_moduli_lower_exceeds_raw_upper":cms*p*p>uc}
    return ans


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    dest=HERE/"independent_results.json"
    out=compute()
    if args.write:
        with dest.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        compare(out,json.loads(dest.read_bytes()))
    print(json.dumps(out))
