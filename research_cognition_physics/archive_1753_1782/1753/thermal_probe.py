"""Can one temperature sweep discriminate an extra participating level in a pair-exchange probe?"""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
CORE=runpy.run_path(str(HERE/"check.py"))


def thermal_signed_response(h,energy,beta):
    pairs=list(itertools.combinations(range(len(h)),2))
    en=np.array([energy[j]+energy[k] for j,k in pairs])
    weights=np.exp(-beta*(en-min(en)))
    return np.dot(weights,CORE["response"](h,1,2))/sum(weights)


def unique_beta(h,energy):
    # Work on the m positive gaps to avoid division by tiny thermal probabilities.
    a=h[1]-h[0]
    gaps=h[2:]-h[1]
    de=energy[2:]-energy[0]
    assert min(de)>0
    def excess(beta):
        v=gaps*np.exp(-beta*de)
        return sum(x*y for x,y in itertools.combinations(v,2))/sum(v)/a-1
    assert excess(0)>0
    lo,hi=0.,1.
    while excess(hi)>0:
        hi*=2
    for _ in range(80):
        mid=(lo+hi)/2
        if excess(mid)>0:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def sector(h):
    x=(h-h[0])/(h[-1]-h[0])
    rows={}
    for r in (1.05,2.,10.,100.):
        xx=np.r_[x,r]
        beta=unique_beta(xx,xx)
        a,b,d=x[1],1-x[1],r-x[1]
        np.testing.assert_allclose(a/d*math.exp(beta*r)+a/b*math.exp(beta),1.,rtol=1e-12)
        below=thermal_signed_response(xx,xx,1.2*beta)
        above=thermal_signed_response(xx,xx,.8*beta)
        assert below<0<above
        assert thermal_signed_response(x,x,.8*beta)<0
        assert thermal_signed_response(x,x,1.2*beta)<0
        rows[str(r)]=dict(beta_star=beta,T_star_in_units_E3_minus_E1=1/beta,
                           normalized_response_at_Tstar_over_1p2=below,
                           normalized_response_at_Tstar_over_p8=above)
    # Five or more lines, including the case where spectator energies differ from the mass source.
    hh=np.r_[x,1.4,3.2,4.7]
    ee=np.array([0.,.02,1.,1.2,1.9,3.])
    beta=unique_beta(hh,ee)
    assert thermal_signed_response(hh,ee,.9*beta)>0>thermal_signed_response(hh,ee,1.1*beta)
    return dict(four_labels=rows,six_label_witness_beta=beta)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    b=CORE["PORT"]["OLD"]["mother"](.125)
    out={}
    for label,index in (("down",1),("charged_leptons",2)):
        eig=np.linalg.svd(b["yW"][index,:,:3],compute_uv=False)[::-1]**2
        out[label]=sector(eig)
    out=json.loads(json.dumps(out),parse_float=lambda x:float(f"{float(x):.11g}"))
    path=HERE/"thermal_probe_results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1753 one-parameter thermal sign test and strict crossing thresholds verified")
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2,allow_nan=False)
            f.write("\n")
