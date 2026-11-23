"""A two-source neutrino mass bound from the spectrum of the very same 24_F."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))

def verify_spectral_inequality():
    rng=np.random.default_rng(1572)
    for _ in range(150):
        u,w=rng.normal(size=(2,3))+1j*rng.normal(size=(2,3))
        u,w=u/np.linalg.norm(u), w/np.linalg.norm(w)
        a,b=np.exp(rng.uniform(-2,2,2))
        mat=a*np.outer(u,u)+b*np.outer(w,w)
        h,l,z=np.linalg.svd(mat,compute_uv=False)
        assert z<1e-12
        assert l/h <= min(a,b)/max(a,b)+1e-12
        assert abs(l*h-a*b*(1-abs(np.vdot(u,w))**2)) < 1e-11
    masses=read("results.json")["center"]["both_real_positive_mass_witness_GeV"]
    t,o,m=[masses[key] for key in ("weak3","color8","dirac32")]
    assert abs((t+o)/(2*m)-1)<1e-10
    singlet=(t+4*m)/5
    assert abs(singlet/read("results.json")["center"]["physical_mass_singlet_GeV"]-1)<1e-10

def compute():
    verify_spectral_inequality()
    oscillation=read("../1370/inputs.json")
    dm=oscillation["NO"]
    ratios={"NO":sqrt(dm["solar"][1]/dm["atm_abs"][2])}
    dm=oscillation["IO"]
    ratios["IO"]=sqrt(1-dm["solar"][2]/dm["atm_abs"][1])
    result={}
    for match,box in read("results.json")["outer_envelopes"].items():
        # Enlarge rounded saved upper endpoint before using it in an inequality.
        x=box["M3_over_M32_envelope"][1]*(1+1e-10)
        heavy_ratio=(4-x)/(5*x)
        result[match]={"singlet_over_triplet_floor":heavy_ratio,"minimum_Yukawa_norm_ratios":{}}
        for eps in (0,.05):
            output={}
            for ordering,observed in ratios.items():
                r=(observed-eps)/(1+eps)
                value=sqrt(heavy_ratio*r)
                clebsch=sqrt(3/5)
                output[ordering]={"solar_atmospheric_singular_ratio_lower":r,
                                  "singlet_over_triplet_Yukawa_minimum":value,
                                  "relative_perturbation_required_if_both_norms_close":
                                      (value-clebsch)/(value+clebsch)}
            result[match]["minimum_Yukawa_norm_ratios"][str(eps)]=output
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    answer=compute()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as f:
            f.write(json.dumps(answer,indent=2)+"\n")
    else:
        assert answer==read("independent_results.json")
        print("1572 adjoint singlet, two-rank spectral inequality and oscillation windows passed")
