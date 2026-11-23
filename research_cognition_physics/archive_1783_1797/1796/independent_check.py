"""Charged partner checks for the two neutral-degenerate aligned flavour sources."""

import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
V, MW, MZ, MT = 246.22, 80.37, 91.1876, 172.5
LEPTONS = np.array([.000510999, .1056584, 1.77686])
DOWN = np.array([.005, .1, 4.8])
UP = np.array([.002, 1.3, MT])
NODES, WEIGHTS = np.polynomial.legendre.leggauss(64)
UNIT_NODES, UNIT_WEIGHTS = (NODES+1)/2, WEIGHTS/2


def rho_shift(charged, neutral):
    x, y = charged**2, neutral**2
    # Stable at exact and near degeneracy; F(x,y)=(x-y)^2 int_0^1 t(1-t)/(tx+(1-t)y) dt.
    f = (x-y)**2 / (8*math.pi**2*V**2)
    return f * np.dot(UNIT_WEIGHTS,UNIT_NODES*(1-UNIT_NODES)/(UNIT_NODES*x+(1-UNIT_NODES)*y))


def mass_bounds(neutral, rho_cap=.00031+2*.00019):
    endpoints = []
    for sign in (-1,1):
        lo, hi = 0., neutral/2
        for _ in range(70):
            mid = (lo+hi)/2
            if rho_shift(neutral+sign*mid,neutral) < rho_cap:
                lo=mid
            else:
                hi=mid
        endpoints.append(neutral+sign*(lo+hi)/2)
    max_square_split = 4*math.pi*V**2
    quartic = [math.sqrt(neutral**2-max_square_split),math.sqrt(neutral**2+max_square_split)]
    permitted = [max(endpoints[0],quartic[0]),min(endpoints[1],quartic[1])]
    assert permitted[0] > neutral-MW and permitted[1] < neutral+MW
    return {"conditional_rho_only":endpoints, "scalar_partial_wave_only":quartic,
            "combined":permitted, "rho_at_on_shell_W_edges":[rho_shift(neutral+sign*MW,neutral) for sign in (-1,1)]}


def qcd_factor(scale):
    def run(alpha,m1,m2,nf):
        return alpha/(1+alpha*(11-2*nf/3)/(2*math.pi)*math.log(m2/m1))
    start=.1181
    at_top=run(start,MZ,MT,5)
    at_heavy=run(at_top,MT,scale,6)
    return (at_top/start)**(12/23)*(at_heavy/at_top)**(4/7)


def matched_inputs():
    r=runpy.run_path(str(HERE.parent/"1795/check.py"))
    old=r["OLD"]
    high=runpy.run_path(str(old["HERE"]/old["CFG"]["charged_background"]))
    flow,_=high["charged_background"](.1)
    table,ckm=r["sources"]()
    return table,ckm,np.array(flow["central"]("down")),np.array(flow["central"]("charged_leptons"))


def widths(row, ckm, yd, ye, vevB, mass):
    z, down, lepton=row
    t=vevB**2
    denom=vevB*math.sqrt(1-t)
    d=qcd_factor(mass)*(math.sqrt(2)*down-t*np.diag(yd))/denom
    e=(math.sqrt(2)*lepton-t*np.diag(ye))/denom
    lphase=(1-(LEPTONS/mass)**2)**2
    light=np.sum(abs(e)**2,axis=0)*lphase
    quark=ckm@d
    qi2=(UP[:,None]/mass)**2
    qj2=(DOWN[None,:]/mass)**2
    a=1-qi2-qj2
    kin=a*a-4*qi2*qj2
    qfloor=3*abs(quark)**2*kin**1.5/a
    return {"mu_over_tau":light[1]/light[2], "electron_over_mu":light[0]/light[1],
            "mu_over_tau_normalized_to_MZ_mass_ratio":light[1]/light[2]*(ye[2]/ye[1])**2,
            "widths_lnu_GeV":(mass/(16*math.pi)*light).tolist(),
            "unknown_up_yukawa_cb_over_tau_lower":qfloor[1,2]/light[2],
            "unknown_up_yukawa_tb_over_tau_lower":qfloor[2,2]/light[2],
            "unknown_up_yukawa_inclusive_b_over_tau_lower":qfloor[:,2].sum()/light[2]}


def run():
    table,ckm,yd,ye=matched_inputs()
    labels=["A_vanishes_first_edge","A_vanishes_third_edge"]
    output={}
    for scale in (5000.,30000.):
        poles=mass_bounds(scale)
        point={key:widths(table[key],ckm,yd,ye,.95,scale) for key in labels}
        equal={key:widths(table[key],ckm,yd,ye,math.sqrt(.75),scale) for key in labels}
        for k in ["mu_over_tau","electron_over_mu"]:
            assert abs(equal[labels[0]][k]/equal[labels[1]][k]-1)<1e-10
        sep=max(p["mu_over_tau"] for p in point.values())/min(p["mu_over_tau"] for p in point.values())
        output[str(int(scale))]={"charged_mass_window_GeV":poles, "vevB_0p95_charged_at_M":point,
                                "vevB_sqrt_0p75_charged_at_M":equal,"largest_symmetric_fractional_ratio_error_for_disjoint_predictions":(sep-1)/(sep+1),
                                "common_quark_amplitude_QCD_factor":qcd_factor(scale)}
    return {"EW_rho0_conditional_upper_minus_one":.00069,"charged_to_neutral_threshold_closed_both_directions_if_mass_above_GeV":(4*math.pi*V**2+MW**2)/(2*MW),
            "two_neutral_mass_scales_GeV":output}


def cleaned(o):
    if isinstance(o,dict):
        return {k:cleaned(v) for k,v in o.items()}
    if isinstance(o,list):
        return [cleaned(v) for v in o]
    if isinstance(o,(float,np.floating)):
        return float(f"{o:.10g}")
    return o


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--check",action="store_true")
    args=p.parse_args()
    answer=cleaned(run())
    data=json.dumps(answer,indent=2,allow_nan=False)+"\n"
    target=HERE/"independent_results.json"
    if args.check:
        assert json.loads(target.read_bytes())==answer
    else:
        target.write_text(data,encoding="utf-8")
    print("1796 charged partner",hashlib.sha256(data.encode()).hexdigest()[:12])
