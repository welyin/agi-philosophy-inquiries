"""CP lower bound, external tag control and a two-channel CPT check. No author imports."""
import argparse
import cmath
import json
from math import asin, hypot, pi, sqrt
from pathlib import Path

import numpy as np


def rescatter(phase):
    # CP-even final-state unitary acting on both conjugate short-distance sources.
    u=np.array([[1,1j],[1,-1j]])/sqrt(2)
    amp=u@np.array([1.,cmath.exp(1j*phase)])
    conjug=u@np.array([1.,cmath.exp(-1j*phase)])
    g,h=abs(amp)**2,abs(conjug)**2
    assert abs(sum(g)-sum(h))<1.e-13
    return (g-h)/(g+h)


def compute():
    theta=.35
    rot=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    v2=np.diag(np.exp([.2j,.9j]))@rot@np.diag(np.exp([.1j,-.6j]))
    lam=v2[0]*v2[1].conj()
    assert abs(lam[0]+lam[1])<1e-15 and abs((lam[0]*lam[1].conj()).imag)<1e-15
    # Fourier 3x3 has nonzero invariant and an explicitly closed row pair.
    n=np.arange(3)
    v3=np.exp(2j*pi*np.outer(n,n)/3)/sqrt(3)
    lam=v3[0]*v3[1].conj()
    j=(lam[2]*lam[0].conj()).imag
    assert abs(sum(lam))<1e-15 and np.isclose(abs(j),sqrt(3)/18)
    assert abs(j)<=min(abs(lam[2]),abs(lam[0]))*abs(lam[1])
    contract=np.array([[1,1],[1,1j]])/2
    lamc=contract[0]*contract[1].conj()
    assert np.linalg.svd(contract,compute_uv=False)[0]<1
    comp=(lamc[0]*lamc[1].conj()).imag
    assert np.isclose(abs(comp),1/16)
    ads=-.403
    # Observable 2-path amplitude bound is independent of the microscopic CP source.
    r=abs(ads)/(1+sqrt(1-ads**2))
    assert np.isclose(2*r/(1+r*r),abs(ads))
    a_ads=rescatter(asin(-ads))
    prompt=-18.2e-4
    a_d=rescatter(asin(-prompt/2))
    assert np.isclose(a_ads[0],ads) and np.isclose(a_d[0]-a_d[1],prompt)
    window=[(-18.2-4*(3.2+.9))*1e-4,(-18.2+4*(3.2+.9))*1e-4]
    assert max(window)<0
    return {"ADS":{"central":ads,"quadrature_ratio_only":abs(ads)/hypot(.056,.011),
                   "published_null_likelihood_significance":8.0,"central_min_amplitude_ratio":r},
            "pure_prompt_Dstar_tag_deltaACP":{"central":prompt,
                    "four_times_linear_marginal_error_window":window},
            "D_all_tags":{"central":-.00154,"error":.00029,"published_significance":5.3},
            "unitary_three_quartet_abs":abs(j),"nonclosed_two_by_two_quartet_abs":abs(comp),
            "independent_CP_source_two_channel_asymmetries":{"ADS_like":a_ads.tolist(),
                    "charm_prompt_difference_example":a_d.tolist()},
            "minimal_additional_direction_supplies_abs_DFT_row_overlap":abs(lam[1]),
            "same_domain_upper_bound_from_1377_standard_template":8.949266579,
            "upper_two_marginal_box_edge":3.797598741}


def round_numbers(x):
    if isinstance(x,dict):
        return {k:round_numbers(v) for k,v in x.items()}
    if isinstance(x,list):
        return [round_numbers(v) for v in x]
    return float(format(x,".10g")) if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    output=Path(__file__).with_name("independent_results.json")
    answer=round_numbers(compute())
    if args.write:
        with output.open("x",encoding="utf-8") as f:
            json.dump(answer,f,indent=2,ensure_ascii=False)
            f.write("\n")
    else:
        assert answer==json.loads(output.read_text(encoding="utf-8"))
    print("1380 independent unitary, charm-only lower bound and CPT control verified.")
