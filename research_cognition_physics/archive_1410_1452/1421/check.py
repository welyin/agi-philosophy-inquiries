"""Test a common direction using a weak interference readout and a tagged nuclear CC readout."""
import argparse
import json
from math import hypot, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def compute():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    e, est, esy = cfg["elastic_coefficient"]
    c, cst, csy = cfg["ground_state_carbon_in_1e42_cm2"]
    (v, dv), (a, da) = cfg["charm_VA_reference_branch"]
    widths = cfg["marginal_widths"]
    k = 4*cfg["sigma0_over_E_in_1e45_cm2_per_MeV"]
    left, right = abs(v+a)/2, abs(v-a)/2
    # Only the two absolute chiral amplitudes enter the calibration from CHARM II.
    # Relax their correlation to make the wrong-sign exclusion conservative.
    llo = left-widths*(dv+da)/2
    rlo = right-widths*(dv+da)/2
    ehi, clo = e+widths*hypot(est, esy), c-widths*hypot(cst, csy)
    c0hi = max(cfg["published_nuclear_predictions"])
    ratio_hi = ehi/clo
    positive_floor = k*((1+llo)**2+rlo**2/3)/c0hi
    scale_limit = ratio_hi/positive_floor
    p_upper = (ehi/k-llo**2-rlo**2/3)/(1+2*llo)
    p_lower = clo/c0hi
    assert 0 < p_upper < p_lower < 1 and scale_limit < 1
    p_fake = (e/k-left**2-right**2/3)/(1+2*left)
    assert 0 < p_fake < 1
    physical_elastic = k*((1-left)**2+right**2/3)
    assert abs(physical_elastic-e) < hypot(est,esy)
    # General density and coherent chirality matrices, including an all-flavor mirror.
    eye = np.eye(3)
    phase = np.exp(.73j)
    axis = np.array([sqrt(.2)*phase, sqrt(.5), sqrt(.3)])
    projection = np.outer(axis, axis.conj())
    rho = np.diag([.7, .2, .1])
    p = float(np.trace(rho@projection).real)
    total_l = -left*eye + projection
    total_r = right*eye
    alternate_nc = left*eye-2*projection
    for y in (0., .24, .78, .99):
        def rate(l, r):
            kernel = l.conj().T@l + (1-y)**2*(r.conj().T@r) - .014*y*(l.conj().T@r+r.conj().T@l)/2
            return float(np.trace(kernel@rho).real)
        assert np.isclose(rate(total_l, total_r), rate(projection+alternate_nc, -total_r))
    integ = np.trace((total_l@total_l+total_r@total_r/3)@rho).real
    assert np.isclose(integ, left**2+right**2/3 +(1-2*left)*p)
    term = total_l+total_r
    opposite = alternate_nc+projection-total_r
    assert np.allclose(opposite,-term)
    return {"round":1421,"cross_section_units_as_registered":True, "positive_left_min_abs":llo,
            "observed_ratio_upper":ratio_hi, "positive_left_ratio_lower":positive_floor,
            "relative_response_max_for_positive_solution":scale_limit,
            "independent_equal_fractional_error_needed_if_opposite":(1-scale_limit)/(1+scale_limit),
            "positive_left_elastic_only_p_at_center":p_fake, "two_width_absolute_elastic_p_upper":p_upper,
            "two_width_carbon_p_lower":p_lower, "negative_left_p1_elastic_coefficient":physical_elastic,
            "conditional_forward_dressed_mirror":True,"new_SM_input_values":0}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    output=compute()
    if args.check:
        previous=json.loads((HERE/"results.json").read_text(encoding="utf-8"))
        assert output.keys()==previous.keys() and all(np.isclose(output[k],previous[k]) for k in previous)
        print("1421 common flux inequality, charged overlap and full-amplitude mirror verified")
    else:
        print(json.dumps(output,indent=2))
