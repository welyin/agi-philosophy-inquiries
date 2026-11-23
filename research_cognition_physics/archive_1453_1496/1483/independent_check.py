"""Independent complex-Gram and two-body threshold checks, masses in GeV."""
from pathlib import Path
from itertools import product
import argparse
import json
import math
import numpy as np

DIR=Path(__file__).resolve().parent
DATA=json.loads((DIR/"inputs.json").read_text(encoding="utf-8-sig"))


def momentum(parent,kaon,invisible):
    if invisible>=parent-kaon:
        return 0.
    return math.sqrt((parent**2-(kaon+invisible)**2)*(parent**2-(kaon-invisible)**2))/(2*parent)


def rate(phases,p,d):
    pz,nz=phases
    return nz/pz*(np.vdot(d,d)/np.vdot(p,p)).real


def run():
    m=DATA["mass_GeV"]
    phases=[4*momentum(m[s],m[k],0)**2/m[s] for s,k in [("p","K+"),("n","K0")]]
    a,b,c=[DATA["W0_GeV2"][z][0] for z in ("us_d","ud_s","ds_u")]
    t=math.tan(DATA["CKM_matching"]["angle_rad"][0])
    x=np.array([a,t*b])
    y=np.array([c,-t*b])
    common=rate(phases,x,y)
    # Haar rotation of complete invisible flavor basis and arbitrary matched relative phases.
    rng=np.random.default_rng(1483)
    z,_=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    rows=np.exp(1j*np.array([.74,-1.17]))[:,None]*z[:2]
    assert abs(rate(phases,x@rows,y@rows)-common)<1e-12
    m0=DATA["illustrative_common_Gram_floor"]
    def best(s):
        r=np.outer(y,y)/(y@y)-s*np.outer(x,x)/(x@x)
        ev=np.linalg.eigvalsh(r)
        return ev[-1]+m0*ev[0]
    lo,hi=1.,1/m0
    for _ in range(65):
        mid=(lo+hi)/2
        if best(mid)>0:
            lo=mid
        else:
            hi=mid
    dmin=(m["n"]-m["K0"])*1e3
    dmax=(m["p"]-m["K+"])*1e3
    mass=.443
    free_p=momentum(m["p"],m["K+"],mass)
    free_n=momentum(m["n"],m["K0"],mass)
    assert dmin<1000*mass<dmax and free_n==0 and free_p>0
    mK=m["K+"]
    electron=.000510999
    muon=.1056584
    assert electron+mass<mK<muon+mass
    return {"ideal_ratio":common,"max_common_Gram_ratio":hi*common,
            "min_common_Gram_ratio":common/hi,"open_only_for_proton_band_MeV":[dmin,dmax],
            "kaon_p_in_free_proton_decay_MeV":{"massless":momentum(m["p"],mK,0)*1e3,
            "443_MeV":free_p*1e3},"threshold_if_electron_kaon_open_MeV":1e3*(mK-electron),
            "threshold_if_muon_kaon_open_MeV":1e3*(mK-muon)}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    mode=p.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=p.parse_args()
    res=run(); path=DIR/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(res,stream,indent=2)
            stream.write("\n")
    else:
        saved=json.loads(path.read_text(encoding="utf-8-sig"))
        def eq(a,b):
            if isinstance(a,dict): return a.keys()==b.keys() and all(eq(a[k],b[k]) for k in a)
            return np.allclose(a,b,rtol=1e-12,atol=1e-12)
        assert eq(saved,res)
        print("1483 independent width, Gram and threshold checks passed.")
