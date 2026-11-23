"""Positive dispersion integral and on-shell electroweak doublet production."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def spectral(m, t1, t2, scalar=False, order=96):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    v=(1+nodes)/2
    r=v**3/4 if scalar else v*(3-v*v)/2
    cut=4*m*m
    denom=(cut+t1*(1-v*v))*(cut+t2*(1-v*v))
    return (t2-t1)/(3*math.pi*ALPHA_INV)*sum(weights*cut*v*r/denom)


def h(x,ms):
    return (2*x+math.log(math.hypot(math.exp(x),ms)))/3


def other(eq,m,ms):
    budget=2*h(math.log(eq),ms)-h(math.log(m),ms)
    low,high=math.log(m),math.log(1.e7)
    assert h(low,ms)<=budget<h(high,ms)
    for _ in range(75):
        mid=(low+high)/2
        if h(mid,ms)>budget:
            high=mid
        else:
            low=mid
    return math.exp((low+high)/2)


def weight(m,energy):
    y=4*(m/energy)**2
    return 0. if y>=1 else math.sqrt(1-y)*(1+y/2)


def currents(energy):
    weak=read("../1369/inputs.json")
    y2=.6*weak["g1_GUT"][0]**2
    w2=weak["g2"][0]**2
    sw=y2/(y2+w2)
    # A stipulated Born convention: alpha(0), weak angle from the measured MZ couplings.
    z=(energy**2/(energy**2-MZ**2))/(sw*(1-sw))
    def helicity(q,t):
        return [(-q+e*(t-q*sw)*z)**2 for e in (-.5+sw,sw)]
    return {"sin2_theta":sw,"charged":sum(helicity(-1,-.5))/2,
            "neutral":sum(helicity(0,.5))/2,
            "unit_pb":4*math.pi/(3*ALPHA_INV**2*energy*energy)*3.89379372e8}


def main():
    old=read("../1698/results.json")["shared_witness"]
    rival=read("../1698/independent_results.json")["actual_zero_boundary_solutions"]
    classes=[("30TeV",old["soft_GeV"],old["weak_fermions_GeV"][0])]+[
        (k,row["soft_GeV"],row["gauge"]["equal_weak_fermions_GeV"]) for k,row in rival.items()]
    windows={name:cfg[key] for name,key in [("opal","opal_transfer_GeV2"),
            ("l3_mid","l3_mid_transfer_GeV2"),("l3","l3_transfer_GeV2")]}
    meters={}
    spectra={}
    for name,ms,eq in classes:
        cap=other(eq,MZ,ms)
        spectra[name]={"equal_GeV":eq,"max_allowed_GeV":cap}
        b=math.hypot(ms,eq)
        meters[name]={win:{"two_fermions":2*spectral(eq,*interval),
                          "four_complex_scalars":4*spectral(b,*interval,scalar=True)}
                      for win,interval in windows.items()}
        assert abs(2*h(math.log(eq),ms)-h(math.log(cap),ms)-h(math.log(MZ),ms))<4e-14
    old_results=read("results.json")["witnesses"]
    worst=0.
    for name,v in meters.items():
        for win,pair in v.items():
            parent=old_results[name]["alpha_windows"][win]
            for a,b in [("two_fermions","weak_fermions"),("four_complex_scalars","weak_scalars")]:
                worst=max(worst,abs(pair[a]/parent[b]-1))
    assert worst<2e-9
    mstar=spectra["30TeV"]["max_allowed_GeV"]
    energy=400.
    cr=currents(energy)
    unit=cr["unit_pb"]
    def production(f):
        return {name:unit*f*cr[name] for name in ("charged","neutral")}
    light_min=weight(MZ,energy)+weight(mstar,energy)
    light_max=2*weight(classes[0][2],energy)
    heavy_max=weight(MZ,energy)
    assert all(2*sp["max_allowed_GeV"]>energy for sp in list(spectra.values())[1:])
    return {"spectral_integrals":meters,"relative_difference_with_parent_fourier_integral":worst,
            "weak_spectra":spectra,"exact_threshold_window_GeV":[2*mstar,2*classes[1][2]],
            "Born_400GeV":{**cr,"light_min_weight":light_min,"light_max_weight":light_max,
                           "heavy_max_weight":heavy_max,"minimum_gap_weight":weight(mstar,energy),
                           "light_min_pb":production(light_min),"light_max_pb":production(light_max),
                           "heavy_max_pb":production(heavy_max),"gap_pb":production(weight(mstar,energy))},
            "Born_300GeV":{**currents(300.),"light_equal_weight":2*weight(classes[0][2],300.),
                           "common_maximally_split_weight":weight(MZ,300.)},
            "same_split_just_below_threshold":{k:weight(MZ,2*mstar-.01)+
                    weight(row["max_allowed_GeV"],2*mstar-.01) for k,row in spectra.items()}}


cfg=read("inputs.json")
ALPHA_INV=cfg["inverse_alpha0"]
MZ=read("../1369/inputs.json")["mu_ref_GeV"]
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=json.dumps(main(),ensure_ascii=False,indent=2)+"\n"
    dest=HERE/"independent_results.json"
    if args.check:
        assert json.loads(result)==json.loads(dest.read_bytes())
        print("1701 independent verified")
    else:
        dest.write_text(result,encoding="utf-8")
