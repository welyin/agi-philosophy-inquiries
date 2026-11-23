"""Minimum Yukawa norm for signed mass sources and a bottomonium interference holdout."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def sign_budget(c, groups):
    weights = []
    for sign, costs in groups:
        roots = np.sqrt(np.array(costs))
        weight = (1+sign*c)/2
        weights.extend((sign, float(p), cost) for p, cost in zip(weight*roots/sum(roots), costs))
    used = sum(p for sign,p,cost in weights)
    signed = sum(sign*p for sign,p,cost in weights)
    np.testing.assert_allclose([used, signed], [1., c], atol=2e-14)
    gram = sum(p*np.outer([1, sign], [1, sign]) for sign,p,cost in weights)
    np.testing.assert_allclose(gram, [[1.,c],[c,1.]], atol=2e-14)
    return sum(cost/p for sign,p,cost in weights)


def rare_decay():
    # CLEO electronic width in keV (not the hadronic product), Eq. 38/43 of 1505.03870.
    alpha = 1/137.036
    width_keV = 1.354
    mv, mh = 9.4603, 125.10
    norm = 6/alpha * width_keV*1e-6/mv*(1-(mv/mh)**2)**2
    out = {}
    for sign in (-1, 1):
        def ratio(real, imag, r):
            return norm*abs(1-(real+1j*imag)*r-(.0184-.0015j))**2
        central = ratio(.948, .130, sign)
        # Conservative component intervals, measured dielectron-width errors add linearly.
        rs = [sign*.9, sign*1.1]
        real = [1-.0184-a*r for a in (.868, 1.028) for r in rs]
        imag = [.0015-a*r for a in (.092, .168) for r in rs]
        def squared_bounds(x):
            lo, hi = min(x), max(x)
            return (0 if lo<=0<=hi else min(lo*lo,hi*hi), max(lo*lo,hi*hi))
        rl, rh = squared_bounds(real)
        il, ih = squared_bounds(imag)
        out[str(sign)] = {"central_ratio_to_diphoton": central,
                          "conditional_ratio_box": [norm*(1-.048/width_keV)*(rl+il),
                                                    norm*(1+.048/width_keV)*(rh+ih)],
                          "central_BR_if_Bgg_is_0p00227": central*.00227}
    assert out["1"]["conditional_ratio_box"][1] < out["-1"]["conditional_ratio_box"][0]
    return out


def run():
    # Full SM, mu=M_h row of Huang and Zhou Tables 2/3; m_f=y_f*246/sqrt(2).
    masses = np.array([[165.05,.607,.00120],[2.768,.05200,.00261],
                       [1.73369,.102065,.00048449]])
    up, down, leptons = np.sum(masses**2, axis=1)*np.array([3,3,1])*2/246.**2
    examples = {}
    for c in (.95, .99, .999):
        examples[str(c)] = {
            "opposite_up_same_down_and_leptons": sign_budget(c, [(1,[down+leptons]),(-1,[up])]),
            "opposite_down_same_up_and_leptons": sign_budget(c, [(1,[up+leptons]),(-1,[down])]),
            "opposite_down_and_leptons_shared": sign_budget(c, [(1,[up]),(-1,[down+leptons])]),
            "opposite_down_and_leptons_separate": sign_budget(c, [(1,[up]),(-1,[down,leptons])])}
    c, k = .97, -.8
    p = (1-c*c)/(1-2*c*k+k*k)
    residual = np.array([[1.,c],[c,1.]]) - p*np.outer([1.,k],[1.,k])
    eigenvalues = np.linalg.eigvalsh(residual)
    assert abs(eigenvalues[0]) < 1e-14 and eigenvalues[1] > 0
    return {"round":1792, "yukawa_source_weights": [up,down,leptons], "budget_minima": examples,
            "finite_top_trial":{"c_lower":c,"k_up_upper":k,"up_vev_weight_max":p,
                                "trace_lower_from_up_sector_alone":up/p},
            "upsilon1S": rare_decay(), "ATLAS_139fb_UL_assuming_SM_production":.00025}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    body = json.dumps(run(), indent=2, allow_nan=False)+"\n"
    path = HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==json.loads(body)
    else:
        path.write_text(body, encoding="utf8")
    print("1792 independent", hashlib.sha256(body.encode()).hexdigest()[:12])
