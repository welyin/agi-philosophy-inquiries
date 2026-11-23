"""Three-bin threshold and Higgs-background contrasts. Couplings are MS-bar proxies."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PER_TRIPLET = 1/(3*math.pi)


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def shift(q, masses):
    """Inverse alpha change from positive Dirac triplets relative to a low anchor."""
    return -PER_TRIPLET*sum(max(0., math.log(q/m)) for m in masses)


def hat(m, q):
    if not q[0] < m < q[2]:
        return 0.
    if m < q[1]:
        return math.log(m/q[0])/math.log(q[1]/q[0])
    return math.log(q[2]/m)/math.log(q[2]/q[1])


def contrast_coeffs(q):
    gaps = np.diff(np.log(q))
    return np.array([-1/gaps[0], 1/gaps[0]+1/gaps[1], -1/gaps[1]])


def compute():
    inp = read("inputs.json")
    cms = inp["measurement"]["profile_table3"]
    setup = inp["theory"]
    q = np.array(cms["Q_GeV"])
    a, da = np.array(cms["alpha_s"]), np.array(cms["symmetric_total"])
    xlo, xhi = 1/(a+da), 1/(a-da)
    coeff = contrast_coeffs(q)
    assert abs(sum(coeff)) < 1e-13 and abs(coeff @ np.log(q)) < 1e-13
    corner = [coeff @ x for x in itertools.product(*zip(xlo, xhi))]
    n = setup["new_Dirac_triplets"]
    b = (11-2*setup["baseline_active_Dirac_triplets"]/3)/(2*math.pi)
    factor = n*PER_TRIPLET
    offsets = np.log(q/q[0])
    # For every mass individually the offset lies between both slopes.
    anchor = [float(max(xlo-(b-factor)*offsets)), float(min(xhi-b*offsets))]
    assert anchor[0] < anchor[1]
    anchor_witness = sum(anchor)/2
    cases = {}
    for m in setup["mass_cases_GeV"]:
        deltas = np.array([shift(point,[m]*n) for point in q])
        contrast = coeff @ deltas
        assert math.isclose(contrast, factor*hat(m,q), abs_tol=1e-13)
        local = anchor_witness+b*offsets+deltas-deltas[0]
        assert all(xlo <= local) and all(local <= xhi)
        cases[str(m)] = {"hat":hat(m,q), "model_contrast":float(contrast),
                         "anchored_alpha":(1/local).tolist()}
    pair = setup["split_pair_GeV"]
    geom = math.sqrt(math.prod(pair))
    same = [shift(point,pair)-shift(point,[geom,geom]) for point in q]
    assert max(abs(d) for d in same) < 1e-12
    sep = shift(geom,pair)-shift(geom,[geom,geom])
    # A finite smooth mass law exp(kappa*h/v) checks the Higgs derivative.
    h_step = 1e-5
    def across(h, kappas):
        moving = [m*math.exp(k*h) for m,k in zip(pair,kappas)]
        return shift(q[2],moving)-shift(q[0],moving)
    derivatives = {}
    for k in (0,1):
        num = (across(h_step,[k,k])-across(-h_step,[k,k]))/(2*h_step)
        assert math.isclose(num, factor*k, rel_tol=1e-9, abs_tol=1e-13)
        derivatives[str(k)] = {"delta_inverse_alpha_s_at_v": across(0,[k,k]),
                               "d_d_h_over_v": num,
                               "hGG_coefficient_in_units_alpha_s_over_8pi": 4*k/3}
    ideal = factor/(2*sum(abs(coeff)))
    return {
        "contrast_coefficients":coeff.tolist(), "table3_proxy_central":float(coeff @ (1/a)),
        "table3_marginal_box_projection":[min(corner), max(corner)],
        "n2_model_contrast_range":[0.,factor], "universal_anchor_interval":anchor,
        "synthetic_common_anchor":anchor_witness, "cases":cases,
        "max_uniform_inverse_alpha_halfwidth_for_separating_peak_vs_zero":ideal,
        "split_and_common_mass_alias":{"masses_GeV":pair, "common_mass_GeV":geom,
                                     "max_abs_at_registered_Q":max(abs(d) for d in same),
                                     "inverse_alpha_difference_at_common_mass":sep},
        "Higgs_response_for_same_spectrum":derivatives}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=5e-10,abs_tol=2e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check",action="store_true")
    modes.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=compute()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        assert same(result,read("results.json"))
        print("1379: CMS proxy box, triangular kernel, threshold aliases and Higgs derivative verified")
