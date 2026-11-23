"""Finite charged-current rank bounds from W widths, chirality, beta, charm and top."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
BRANCHING = [(0.1083,0.0010),(0.1094,0.0008),(0.1077,0.0021)]  # CMS 2201.07861
WIDTH = (2.202,0.047)  # ATLAS 2403.15085, GeV
W_UNIT_MAX = .260  # one fully open chirality, including pole/rate matching
RATE_RETENTION = .95  # example lower fraction after template recalibration


def calculate():
    br = [x-2*err for x,err in BRANCHING]
    width = WIDTH[0]-2*WIDTH[1]
    totals = np.array(br)*width/W_UNIT_MAX
    h_tau_upper = -.995+2*math.hypot(.010,.003)
    # Product of e and mu left fractions exceeds Q_LL in the registered vector model.
    qll = .955
    tau_left = (1-h_tau_upper)/2
    # Infimum of a*x+b*y, xy>=qll, 0<=x,y<=1. Evaluate the stationary point and boundaries.
    a,b,c = totals
    x = min(1.,max(qll,math.sqrt(qll*b/a)))
    y = qll/x
    left_sum = a*x+b*y+c*tau_left
    weak_norm_upper = 1.02*(1/125)*81/(12*.22)  # common independent gauge-budget example
    assert weak_norm_upper < W_UNIT_MAX
    assert abs(a*x+b*y-min(a*u+b*qll/u for u in np.linspace(qll,1,501)))<1e-7
    ud = .97373-2*.00031
    cs = .9663-2*math.sqrt(sum(x*x for x in [.0053,.0039,.0019,.0040]))
    tb = .988-2*.051  # fixed known decay channels, floating t-row mix; CMS BSM scenario 1
    vec = [ud*math.sqrt(b*RATE_RETENTION), cs*math.sqrt(min(a,b)*RATE_RETENTION)]
    top = .95*tb  # left vertex transport/projection budget
    threshold = sum(vec)+top-2
    cd = .97373
    sd = math.sqrt(1-cd**2)
    v = np.array([[cd,sd,0.],[-sd,cd,0.],[0.,0.,1.]])
    left = np.diag([0.,0.,1.])@v
    right = np.diag([1.,1.,0.])@v
    assert np.linalg.norm(left+right-v)<1e-14 and np.linalg.matrix_rank(left)==1
    assert abs(cd-.9663)<2*math.sqrt(sum(x*x for x in [.0053,.0039,.0019,.0040]))
    return dict(width_lower_GeV=width,branching_lower=br,
                upper_GeV_one_chirality=W_UNIT_MAX, independent_electroweak_budget=weak_norm_upper,
                inclusive_row_floors=totals.tolist(), e_mu_QLL_floor=qll, tau_fraction_lower=tau_left,
                minimizing_e_mu_left_fractions=[x,y], left_trace_lower=left_sum,
                sufficient_common_retention=2/left_sum, example_retention=RATE_RETENTION,
                retained_left_trace_lower=left_sum*RATE_RETENTION,
                unit_width_critical_if_perfect_matching=W_UNIT_MAX*left_sum/2,
                conditional_quark_vector_lower=vec, conditional_top_left_lower=top,
                available_light_right_amplitude_sum_budget=threshold,
                mirror_quark_witness={"abs_diagonals":list(map(float,abs(np.diag(v)))),
                                      "left_rank":1,"right_rank":2})


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    result=json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"independent_results.json"
    if check:
        assert json.loads(dest.read_bytes()) == json.loads(result)
    else:
        dest.write_text(result,encoding="utf-8")
    print("1781 charged currents and chirality",hashlib.sha256(result.encode()).hexdigest()[:12])
