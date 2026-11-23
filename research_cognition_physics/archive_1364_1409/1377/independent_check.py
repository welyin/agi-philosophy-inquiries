"""Independent leading-order 1377 check; offline, Python + NumPy."""
import argparse
import json
from math import asin, log, pi, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
MH = 125.0
V, DV = 248.401, .032
# MZ, up/charm/top/down/strange/bottom. From the archived Eq. (4) inputs.
Y = np.array([7.04e-6, .00356, .967, 1.54e-5, .000306, .01630])
DY = np.array([.15e-6, .00006, .004, .02e-5, .000004, .00009])


def form(tau):
    if tau < 1:
        nodes, wg = np.polynomial.legendre.leggauss(32)
        x = (1 + nodes[:, None]) / 2
        y = (1 - x) * (1 + nodes[None, :]) / 2
        weights = wg[:, None] * wg[None, :] * (1 - x) / 4
        return 4 * np.sum(weights * (1 - 4*x*y) / (1 - tau*4*x*y))
    if tau == 1:
        return 2.
    # Product-free branch log well conditioned even for electron-scale masses.
    absorptive = log(tau) + 2*np.log1p(sqrt(1-1/tau)) - 1j*pi
    return 2/tau * (1 - .25*(1-1/tau)*absorptive**2)


def radius_enclosure():
    tau_low = MH**2 / (2*(V+DV)**2*(Y+DY)**2)
    light = np.delete(tau_low, 2)
    assert (light > 10).all()
    light_sum_bound = sum(2/light * (1 + .25*(np.log(4*light)**2 + pi**2)))
    top_tau_high = MH**2 / (2*(V-DV)**2*(Y[2]-DY[2])**2)
    assert 0 < top_tau_high < 1
    # In the triangle 0 <= 4xy <= 1.  Light bound decreases for all tau >= 10.
    top_bound = (4/3) / (1-top_tau_high)
    re_min = 4/3 - light_sum_bound
    abs_max = top_bound + light_sum_bound
    new = 8/3
    return {"sum_abs_light_upper": light_sum_bound, "top_abs_upper": top_bound,
            "SM_real_lower": re_min, "SM_modulus_upper": abs_max,
            "pure_gg_lower": 1 + (2*new*re_min + new**2) / abs_max**2}


def compute():
    tau = MH**2 / (2*V**2*Y**2)
    total = sum(map(form, tau))
    for t in [.001, .16, .4, .7]:
        assert abs(form(t)-2*(t+(t-1)*asin(sqrt(t))**2)/t**2) < 1.e-9
    pure = abs(1+(8/3)/total)**2
    h250 = abs(1+2*form((MH/500)**2)/total)**2
    eta = .49 * .0262 / 1.17
    bound = radius_enclosure()
    # Do not assume the size of the positive bbH contribution from the new model.
    bound["group_over_vbf_lower_with_positive_bbH"] = (1-eta)*bound["pure_gg_lower"]
    # Direct rounded cross sections from Fig. 10, fb; Table 8 SM normalization in fb.
    rect = [(1120-130)/(110+40)*92/1170, (1120+130)/(110-40)*92/1170]
    assert bound["group_over_vbf_lower_with_positive_bbH"] > 5 > rect[1]
    assert abs(total-(1.3277768094+.05304766592j)) < 1e-9
    assert abs(pure-9.037458798)<1e-7 and abs(h250-9.218942666)<1e-7
    y=.12
    z=sqrt(15)*y
    x=sqrt(1-16*y*y)
    weights=np.array([1.,1.,16.])
    vacuum=np.array([x,z,y/4])
    h=np.array([x,0.,4*y])
    single=np.dot(weights*vacuum,h)
    double=np.dot(weights*h,h)
    defect=double-single**2
    assert np.isclose(np.dot(vacuum*weights,vacuum),1) and np.isclose(h@h,1)
    assert np.isclose(single,1) and np.isclose(defect,240*y*y)
    assert np.isclose(defect,np.dot(weights,(h-single*vacuum)**2))
    assert defect >= z*z/(1-z*z)
    return {"m_h_GeV":MH,"SM_amplitude":[total.real,total.imag],
            "central_extra_family_Rgg_infimum":pure,"two_250_GeV_Rgg":h250,
            "central_bbH_fraction_upper":eta,"central_group_lower_with_SM_bbH":(1-eta)*pure+eta,
            "full_input_box":bound,
            "ATLAS_table8_signal_strength_ratio":.96/1.21,
            "ATLAS_figure10_rounded_ratio_one_marginal_width":rect,
            "ATLAS_figure10_ggF_VBF_correlation":-.22,
            "counterexample":{"single_W":single,"double_W":double,"new_kappa":h[1]/z,
                              "geometric_defect":defect,"lower_bound_from_doublet_alone":z*z/(1-z*z)},
            "source":"https://arxiv.org/pdf/2004.03447#page=37"}


def rounded(x):
    if isinstance(x, dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x, list):
        return [rounded(v) for v in x]
    return float(format(x, ".10g")) if isinstance(x,(float,np.floating)) else x


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    p.add_argument("--check",action="store_true",help="compare with the independent result (default)")
    args=p.parse_args()
    result=rounded(compute())
    archive=HERE/"independent_results.json"
    if args.write:
        with archive.open("x",encoding="utf-8") as file:
            json.dump(result,file,indent=2,ensure_ascii=False)
            file.write("\n")
    else:
        assert result == json.loads(archive.read_text(encoding="utf-8"))
    print("1377 independent triangle integrals, whole mass rectangle and ATLAS bounds verified.")
