"""Exact two-family extrema and one-loop QED matching, followed by radial derivatives."""
import argparse
import json
import math
from pathlib import Path
import numpy as np


HERE=Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8-sig"))


def L(q,s):
    b2=q*q+(1-q*q)*s*s
    trace=12*(q+1/q)-143*q/(12*b2)
    ans=(trace+math.sqrt(trace*trace-4))/2
    direct=np.array([[0.,.5*q/math.sqrt(b2)],
        [6*math.sqrt(b2),6*(1-q*q)*s*math.sqrt(1-s*s)/math.sqrt(b2)]])
    assert math.isclose(ans,np.linalg.cond(direct),rel_tol=5e-14)
    return ans


def ratio_qed(R,alpha):
    return R*(1+1.5*alpha/math.pi*math.log(R))


def derive():
    i=read(HERE/"inputs.json")
    mass=read(HERE/i["quark_masses_and_muonium"])
    d,ed=mass["lattice"]["md"]; m,em=mass["lattice"]["ms"]
    pole,ep=mass["muonium"]["pole_ratio"]
    _,s,errors=next(x for x in read(HERE/i["charm"])["edges"] if x[0]=="cd")
    qlo,qhi=(d-2*ed)/(m+2*em),(d+2*ed)/(m-2*em)
    slo,shi=(s-2*sum(errors))/1.02,(s+2*sum(errors))/.98
    assert 0<qlo<qhi<slo<shi<1
    lo,hi=L(qhi,slo),L(qlo,shi)
    central=L(d/m,s)
    inverse=i["inverse_alpha_on_shell"]
    common=ratio_qed(pole,1/inverse)
    t=common/central
    rem=i["differential_transport_remainder"]
    tlo=ratio_qed(pole-2*ep,1/inverse)*rem[0]/hi
    thi=ratio_qed(pole+2*ep,1/inverse)*rem[1]/lo
    weak=read(HERE/i["weak_scale_charged_yukawas"])["yukawa"]["charged_leptons"]
    high=read(HERE/i["SM_running_table"])["scales"]["1e16"]
    rates={}
    diff=2*(1-t)/(1+t)
    for name,ke,km in [("electron_only",3-2*t,1.),("muon_only",1.,3-2/t),
                        ("minimax",1+diff,1-diff)]:
        assert math.isclose((3-ke)/(3-km),t)
        total=1+i["dimuon_SM_BR"]*(km*km-1)+i["cms_dielectron"]["SM_BR_scale"]*(ke*ke-1)
        rates[name]={"kappa_e":ke,"kappa_mu":km,"tree_signal_mu":km*km/total,
                     "electron_tree_signal_scale":i["cms_dielectron"]["SM_BR_scale"]*ke*ke/total}
    # x=(v+h)/v; powers x^1, x^3, x^5, x^7 at x=1.
    from math import prod
    def derivative(terms,order):
        return sum(coeff*prod(power-j for j in range(order)) for power,coeff in terms)
    for terms,expected in [([(1,-1),(3,2),(5,-1)],(0,0,-8)),
                            ([(3,3),(5,-3),(7,1)],(1,1,0))]:
        assert tuple(derivative(terms,j) for j in (0,1,2))==expected
    return {"texture_L_central":central,"texture_L_global":[lo,hi],
            "weak_common_QED_ratio":common,"QED_factor":common/pole,
            "SM_tabulated_ratios_MZ_GUT":[weak[1][0]/weak[0][0],high[4][0]/high[3][0]],
            "T_pole":pole/central,"T_QED":t,"T_box_pureQED":[ratio_qed(pole-2*ep,1/inverse)/hi,
                          ratio_qed(pole+2*ep,1/inverse)/lo],"T_box_with_transport":[tlo,thi],
            "minimax_central":diff,"minimax_lower_inc_transport":2*(1-thi)/(1+thi),
            "absolute_lower_extra_differential_transfer_to_remove":lo/ratio_qed(pole+2*ep,1/inverse)-1,
            "examples":rates,"ee_limit_over_SM_scale":i["cms_dielectron"]["BR_95CL"]/i["cms_dielectron"]["SM_BR_scale"]}


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    return np.allclose(a,b,rtol=1e-11,atol=1e-15)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true"); mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    values=derive()
    dest=HERE/"results.json"
    if args.write:
        with dest.open("x",encoding="utf-8") as f:
            json.dump(values,f,indent=2,allow_nan=False); f.write("\n")
    else:
        assert equal(values,read(dest))
    print(json.dumps(values,indent=2))
