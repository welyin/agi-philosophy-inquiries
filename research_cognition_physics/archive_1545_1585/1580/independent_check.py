"""Independent 1580 Gram-source and neutral-scalar flavour checks."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
CMS_EM = 4.4e-5
RATE_TRANSFER = 1.01/.99
ATMOS_23 = {"NO": (.435,.585), "IO": (.440,.584)}


def read(n):
    return json.loads((STAGE/str(n)/"inputs.json").read_text(encoding="utf-8-sig"))


def eps_interval():
    charged = read(1386)["yukawa"]["charged_leptons"]
    e, mu, tau = [(a-2*b,a+2*b) for a,b in charged]
    return ((mu[0]-e[1])/(tau[1]-e[1]), (mu[1]-e[0])/(tau[0]-e[0]))


def oscillation_bounds(order, box):
    _, sl, sh = box["s12"]
    _, zl, zh = box["s13"]
    tl, th = ATMOS_23[order]
    solar = box["solar"][2]
    atmosphere = box["atm_abs"][1]
    r = solar/(atmosphere if order == "NO" else atmosphere-solar)
    solar_mix = math.sqrt(sh*(1-sh))
    k = r*solar_mix / (math.sqrt(zl)*(1-r*sh))
    tan23 = math.sqrt(th/(1-th))
    assert k*tan23 < 1
    em_over_et = (tan23+k)/(1-k*tan23)
    common = math.sqrt(zh)*(1+r*sh)
    em_max = math.sqrt(1-zl)*(common*math.sqrt(th)+r*solar_mix*math.sqrt(1-tl))
    et_max = math.sqrt(1-zl)*(common*math.sqrt(1-tl)+r*solar_mix*math.sqrt(th))
    mt_min = (1-zh)*math.sqrt(min(tl*(1-tl),th*(1-th)))-r/2
    eps = eps_interval()[1]
    et_em = (eps*em_over_et)**2*RATE_TRANSFER
    mt_em = (eps*em_max/((1-eps)*mt_min))**2*RATE_TRANSFER
    mt_et = (et_max/((1-eps)*mt_min))**2*RATE_TRANSFER
    return {"r21_over_abs_31_upper":r, "F_em_over_et_upper":em_over_et,
            "F_em_over_abs_31_upper":em_max,"F_et_over_abs_31_upper":et_max,
            "F_mt_over_abs_31_lower":mt_min,
            "rate_et_over_em_upper":et_em, "rate_mt_over_em_upper":mt_em,
            "rate_mt_over_et_upper":mt_et, "B_et_from_CMS_em":et_em*CMS_EM,
            "B_mt_from_CMS_em":mt_em*CMS_EM}


def mixing(x,y,z,phase):
    s1,s2,s3 = np.sqrt([x,y,z])
    c1,c2,c3 = np.sqrt([1-x,1-y,1-z])
    h = np.exp(1j*phase)
    return np.array([[c1*c3,s1*c3,s3/h],
                     [-s1*c2-c1*s2*s3*h,c1*c2-s1*s2*s3*h,s2*c3],
                     [s1*s2-c1*c2*s3*h,-c1*s2-s1*c2*s3*h,c2*c3]])


def matrix_test(order, box, bounds):
    eps = sum(eps_interval())/2
    e = np.diag([0,eps,1])
    for delta in [0.,math.pi,0.83,4.18,5.21]:
        u = mixing(box["s12"][0],sum(ATMOS_23[order])/2,box["s13"][0],delta)
        r = box["solar"][0]/box["atm_abs"][0]
        spec = [0.,r,1.] if order=="NO" else [1-r,1.,0.]
        f = (u*np.array(spec)) @ u.conj().T
        c = e@f-f@e
        d = c.conj().T@c
        for i,j,k in ((0,1,2),(0,2,1),(1,2,0)):
            assert np.isclose(d[i,j],(e[k,k]-e[i,i])*(e[k,k]-e[j,j])*f[i,k]*f[k,j])
        for numerator,denominator,key in [((0,2),(0,1),"et_over_em"),
                                           ((1,2),(0,1),"mt_over_em"),
                                           ((1,2),(0,2),"mt_over_et")]:
            assert RATE_TRANSFER*abs(d[numerator]/d[denominator])**2 <= bounds["rate_"+key+"_upper"]
        weak_j = (u[0,0]*u[1,1]*u[0,1].conj()*u[1,0].conj()).imag
        invariant = np.trace(c@d).imag
        predicted = 6*eps*(1-eps)*abs(np.prod([spec[j]-spec[i] for i,j in itertools.combinations(range(3),2)]))*abs(weak_j)
        assert np.isclose(abs(invariant),predicted,atol=1e-15)
        k0 = 1-eps
        other = d+k0*((f[0,0]+f[1,1]).real*f-f@f+np.eye(3))
        assert np.min(np.linalg.eigvalsh(other))>=-1e-12
        assert abs(other[0,1])<1e-14
        assert abs(np.trace(c@other).imag-invariant)<1e-15
        if delta == 4.18:
            saved = {"delta_rad":delta, "weak_J":float(weak_j),"signed_source":float(invariant),
                     "D_em_et_mt_abs":[float(abs(d[i,j])) for i,j in [(0,1),(0,2),(1,2)]],
                     "Dprime_em_et_mt_abs":[float(abs(other[i,j])) if abs(other[i,j])>1e-13 else 0.
                                              for i,j in [(0,1),(0,2),(1,2)]],
                     "Dprime_smallest_eigenvalue":float(np.linalg.eigvalsh(other)[0])}
    return saved


def necklaces():
    mismatch = {}
    for length in range(1,7):
        bad = []
        for word in itertools.product("EF",repeat=length):
            w = "".join(word)
            if w[::-1] not in [w[i:]+w[:i] for i in range(length)]:
                bad.append(w)
        mismatch[length] = len(bad)
    assert [mismatch[n] for n in range(1,6)] == [0]*5 and mismatch[6] > 0
    return mismatch


def clean(obj):
    if isinstance(obj,float):
        return float(f"{obj:.12g}")
    if isinstance(obj,dict):
        return {str(k):clean(v) for k,v in obj.items()}
    if isinstance(obj,list) or isinstance(obj,tuple):
        return [clean(x) for x in obj]
    return obj


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write",action="store_true")
    args = parser.parse_args()
    cap = read(1370)
    bounds = {o:oscillation_bounds(o,cap[o]) for o in ATMOS_23}
    records = clean({"round":1580, "words_not_cyclically_equal_to_their_reverse":necklaces(),
               "charged_gap_fraction_bounds":eps_interval(), "endpoint_transfer":RATE_TRANSFER,
               "CMS_em_H125_limit":CMS_EM,"oscillation_enclosures":bounds,
               "same_source_competitors":{o:matrix_test(o,cap[o],bounds[o]) for o in ATMOS_23}})
    path = HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as fh:
            json.dump(records,fh,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)
            fh.write("\n")
    elif records != json.loads(path.read_text(encoding="utf-8")):
        raise SystemExit("independent result mismatch")
    print("1580 independent check passed")


if __name__ == "__main__":
    main()
