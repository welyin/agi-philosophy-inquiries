"""Reproduce the neutral-current inversion and its finite-window holdouts."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def limits(pair, n):
    value, width = pair
    return value-n*width, value+n*width


def invert(mz, va, aa, gf):
    factor = mz * math.sqrt(math.sqrt(2) * gf)
    left, right = factor*(va+aa), factor*(va-aa)
    triple = right-2*left
    charge2 = right*triple
    weak = math.sqrt(-2*(left-right)*triple)
    assert math.isclose(weak**2, charge2+triple**2)
    return dict(left=left, right=right, triple=triple, weak=weak, inverse_alpha=4*math.pi/charge2)


def certificate():
    cfg = load(HERE / "inputs.json")
    prior = load(HERE / cfg["prior_Z_and_W"])
    gf = load(HERE / cfg["muon_normalisation"])["formula"]["Gmu_GeV_minus2"]
    pars = [cfg["Z_mass_GeV"], prior["Z_flavor_specific"]["gV"], prior["Z_flavor_specific"]["gA"]]
    center = invert(*(p[0] for p in pars), gf)
    windows = []
    for n in (1,2):
        corners = [invert(*v, gf) for v in itertools.product(*(limits(p,n) for p in pars))]
        intervals = {key:[min(p[key] for p in corners),max(p[key] for p in corners)] for key in center}
        electromagnetic = limits(cfg["dispersion_inverse_alpha_MZ"],n)
        factor_interval = [electromagnetic[0]/intervals["inverse_alpha"][1],
                           electromagnetic[1]/intervals["inverse_alpha"][0]]
        charged = []
        for width,br,mw in itertools.product(limits(prior["W_width_GeV"],n),limits(prior["BR_e"],n),prior["W_mass_window_GeV"]):
            charged.append(math.sqrt(48*math.pi*width*br/mw))
        assert min(charged)<intervals["weak"][0]<intervals["weak"][1]<max(charged)
        observed = cfg["LEP2_one_parameter_g1Z"]["published_95"]
        needed_ref = [intervals["triple"][0]/observed[1],intervals["triple"][1]/observed[0]]
        windows.append(dict(marginal_widths=n,pole_prediction=intervals,
                            measured_inverse_alpha=electromagnetic,
                            predicted_e_squared_over_external_e_squared=factor_interval,
                            smallest_relative_charge_squared_gap=max(0.,factor_interval[0]-1.,1.-factor_interval[1]),
                            direct_weak=[min(charged),max(charged)], admissible_TGC_reference=needed_ref))
    assert windows[0]["smallest_relative_charge_squared_gap"] > 0
    assert windows[1]["smallest_relative_charge_squared_gap"] == 0
    mz = cfg["Z_mass_GeV"][0]
    mw = prior["nominal_W_mass_GeV"]
    direct = math.sqrt(48*math.pi*prior["W_width_GeV"][0]*prior["BR_e"][0]/mw)
    alpha = cfg["dispersion_inverse_alpha_MZ"][0]
    charge = math.sqrt(4*math.pi/alpha)
    z0 = math.sqrt(direct**2-charge**2)
    tgc_ref_gmu = math.sqrt(4*math.sqrt(2)*gf)*mw**2/mz
    tgc_ref_running_os = charge*(mw/mz)/math.sqrt(1-(mw/mz)**2)
    cosine = cfg["hidden_neutral_rotation_cosine"]
    sine = math.sqrt(1-cosine**2)
    t = cosine*z0
    r0, l0 = charge**2/z0,(charge**2-direct**2/2)/z0
    r, l = center["right"],center["left"]
    xr, xl = (r-cosine*r0)/sine,(l-cosine*l0)/sine
    # X'=-s Z0+c X so that the massive-plane rotation is orthogonal.
    rot = np.array([[cosine,sine],[-sine,cosine]])
    couplings = rot @ np.array([[z0,r0,l0],[0.,xr,xl]])
    assert np.allclose(couplings[0],[t,r,l])
    assert np.allclose([couplings[:,0] @ couplings[:,1],couplings[:,0] @ couplings[:,2]],
                       [charge**2,charge**2-direct**2/2])
    hidden_norm = abs(couplings[1,0])
    needed = abs(charge**2-t*r)/hidden_norm
    assert math.isclose(abs(couplings[1,1]),needed)
    return dict(pole_center=center, direct_weak_center=direct,
                electromagnetic_and_direct_W_prediction= z0,
                mass_ratio_ansatz_at_same_direct_W=direct*mw/mz,
                central_predicted_over_observed_e_squared=alpha/center["inverse_alpha"],
                central_TGC_ratio_if_Gmu_on_shell=center["triple"]/tgc_ref_gmu,
                central_TGC_ratio_if_running_e_on_shell=center["triple"]/tgc_ref_running_os,
                reference_coefficients=[tgc_ref_gmu,tgc_ref_running_os],windows=windows,
                extra_neutral=dict(negative_WWZ=t, fitted_electron_Z=[l,r],
                                   invisible_triple_norm=hidden_norm, needed_eR_hidden_norm=needed,
                                   new_electron_gauge_charges_times_coupling=[xl,xr],
                                   hidden_physical_chiral=[float(couplings[1,2]),float(couplings[1,1])],
                                   signed_interference_product=t*r,
                                   no_extra_neutral_required_eR_T3_for_full_negative_branch=(-z0*r-charge**2)/direct**2))


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(v,b[k]) for k,v in a.items())
    if isinstance(a,(list,tuple)):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return math.isclose(a,b,rel_tol=2e-12,abs_tol=1e-14)
    return a==b


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    choice=parser.add_mutually_exclusive_group()
    choice.add_argument("--write",action="store_true")
    choice.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=certificate()
    target=HERE/"results.json"
    if args.write:
        with target.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        assert equal(result,load(target))
        print("1475: matched neutral inversion, charge and W comparisons, LEP2 transfer and hidden-state saturation passed.")
