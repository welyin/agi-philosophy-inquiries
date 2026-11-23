"""Independently invert beta moments before testing the mass trace."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def box(row):
    a, low, high = row
    return a-2*low, a+2*high


def possible(q_box, p_box):
    """Interval cover in the measured moment and gap variables; no root scans."""
    saved = {}
    for label, phase in (("unit", 1.), ("cube", -.5)):
        q = q_box[1]
        kmax = math.floor(math.sqrt(3*(q+1/q-2)/(4*(1+phase))))
        for k in range(kmax+1):
            length = 3+2*k*k
            branches = set()
            qs = np.linspace(*q_box, 129)
            ps = np.linspace(*p_box, 25)
            for a,b in zip(qs,qs[1:]):
                for c,d in zip(ps,ps[1:]):
                    lo = length*a / (6*(d*(b*b-1)+1/3))
                    hi = length*b / (6*(c*(a*a-1)+1/3))
                    min_sum = 4*lo + length**2/hi + 8*k*k*phase
                    max_sum = 4*hi + length**2/lo + 8*k*k*phase
                    if max_sum < 6*(a+1/a)-1e-10 or min_sum > 6*(b+1/b)+1e-10:
                        continue
                    if hi >= length/2-1e-10:
                        branches.add("large")
                    if lo <= length/2+1e-10:
                        branches.add("small")
            if branches:
                saved[f"{label}_{k}"] = sorted(branches)
    return saved


def inspect(k, q, cosine, solar):
    length = 3+2*k*k
    trace = 6*(q+1/q)-8*k*k*cosine
    r = (trace + math.sqrt(trace*trace-16*length*length))/8
    m_b = math.sqrt(solar*q/(6*r))
    a = np.array([0.,1.,1.])
    b = np.array([1.,1.+k,k-1.])
    m = m_b*(r*np.outer(a,a)+np.exp(1j*math.acos(cosine))*np.outer(b,b))
    m1,m2,m3 = sorted(np.linalg.svd(m,compute_uv=False))
    assert m1 < 1e-14 and abs(m2**2 - solar) < 1e-15 and abs(m3/m2-q)<1e-10
    assert np.isclose(np.linalg.norm(m[0]),abs(m[0,0])*math.sqrt(length))
    return {"m_bb_meV": 1000*abs(m[0,0]), "m_beta_meV": 1000*np.linalg.norm(m[0])}


def compute():
    inp=json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig"))
    solar=box(inp["KamLAND_dm21_eV2"])
    atm=box(inp["T2K_dm32_NO_eV2"])
    effective=box(inp["MINOS_effective_atm_eV2"])
    p=[(1-math.sqrt(1-A))/2 for A in box(inp["DayaBay_sin2_2theta13"])]
    t2k=[math.sqrt(1+atm[0]/solar[1]), math.sqrt(1+atm[1]/solar[0])]
    minos=[math.sqrt(effective[0]/solar[1]),math.sqrt(1+effective[1]/solar[0])]
    branches_t2k=possible(t2k,p)
    branches_minos=possible(minos,p)
    assert branches_t2k=={"cube_2":["large"]}
    assert branches_minos=={"cube_2":["large"],"unit_1":["large"]}
    mbeta=[1e3*math.sqrt(s/3+t*v) for s,t,v in zip(solar,atm,p)]
    q=math.sqrt(1+inp["T2K_dm32_NO_eV2"][0]/inp["KamLAND_dm21_eV2"][0])
    solar_needed=[1-2/(3*(1-x)) for x in reversed(p)]
    assert solar_needed[0]>box(inp["JUNO_sin2_theta12_NO"])[1]
    q_alt = 5.23
    # For this near-corner frequency ratio use an actual simultaneous solar input.
    # The representative spectrum printed below keeps its same explicit solar center.
    alt_solar = solar[1]
    dm31_alt = q_alt*q_alt*alt_solar
    assert effective[0] <= dm31_alt and 0 <= dm31_alt-effective[0] <= alt_solar
    return {"T2K_branches":branches_t2k,"MINOS_branches":branches_minos,
            "T2K_primary_witness":inspect(2,q,-.5,inp["KamLAND_dm21_eV2"][0]),
            "MINOS_competitor_spectral_witness":inspect(1,5.23,1.,inp["KamLAND_dm21_eV2"][0]),
            "primary_beta_meV_outer":mbeta,
            "primary_mbb_meV_outer":[v/math.sqrt(11) for v in mbeta],
            "normalized_electron_row_norm_to_ee":math.sqrt(11),
            "solar_needed":solar_needed}


def equal(x,y):
    if isinstance(x,dict):
        return x.keys()==y.keys() and all(equal(x[k],y[k]) for k in x)
    if isinstance(x,list):
        return len(x)==len(y) and all(equal(a,b) for a,b in zip(x,y))
    if isinstance(x,str):
        return x==y
    return math.isclose(x,y,rel_tol=1e-11,abs_tol=1e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    options=parser.add_mutually_exclusive_group()
    options.add_argument("--write",action="store_true")
    options.add_argument("--check",action="store_true")
    result=compute()
    destination=ROOT/"independent_results.json"
    if parser.parse_args().write:
        with destination.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        assert equal(result,json.loads(destination.read_text(encoding="utf-8")))
    print(json.dumps(result,indent=2))
