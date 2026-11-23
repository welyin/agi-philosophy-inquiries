"""One-loop uniform bound at an independent hypercharge--weak matching scale."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# GUT-normalized hypercharge. Antusch et al., arXiv:2510.01312v2, table 2, 10^9 GeV.
COUPLINGS = [(0.509051, 0.000035), (0.57639, 0.00002), (0.6872, 0.0008)]
BETA = [4.1, -19/6, -7]
L = 16 * math.pi**2


def evolved_square(g, b, t):
    return 1 / (g**(-2) - 2*b*t/L)


def calculate():
    boxes = [[x-2*s,x+2*s] for x,s in COUPLINGS]
    times = [(a**(-2)-b**(-2))*L/(2*(BETA[0]-BETA[1]))
             for a,b in itertools.product(*boxes[:2])]
    lo, hi = min(times), max(times)
    # Outward bounds on Tr Hu, Tr Hd and Tr He from the same published row.
    initial_trace_max = 3 * ((0.5838+2*0.0048)**2 + 0.003**2+0.0001**2) + 9*0.009**2 + 3*0.011**2
    denominator = 1 - 5*initial_trace_max*hi/L
    assert denominator > 0
    trace_max = initial_trace_max/denominator
    g1_sq_max = evolved_square(boxes[0][1],BETA[0],hi)
    g3_sq_min = evolved_square(boxes[2][0],BETA[2],hi)
    # Uniform envelope of L*d(log determinant ratio)/dt.
    upper_slope = 6*g1_sq_max - 24*g3_sq_min + trace_max/2
    assert upper_slope < 0
    data = json.loads((HERE.parent/"1402"/"inputs.json").read_text(encoding="utf8"))
    center, half = np.array(data["scales"]["1e9"]).T
    start_upper = np.prod((center+2*half)[:3])/np.prod((center-2*half)[3:])
    upper_ratio = start_upper*math.exp(upper_slope*lo/L)
    eta = (1-upper_ratio**(1/3))/(1+upper_ratio**(1/3))
    return dict(gauge_pair="g1_GUT=g2; separate matching condition from g2=g3",
                log_scale_shift_bounds=[lo, hi], scale_GeV_bounds=[1e9*math.exp(t) for t in [lo,hi]],
                trace_initial_upper=initial_trace_max, trace_uniform_upper=trace_max,
                beta_numerator_uniform_upper=upper_slope, starting_volume_ratio_upper=start_upper,
                volume_ratio_upper_at_match=upper_ratio,
                g3_over_matched_g_lower=math.sqrt(g3_sq_min/g1_sq_max),
                minimum_max_fractional_singular_gain_change_for_equality=eta)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.check:
        saved=json.loads((HERE/"independent_results.json").read_text(encoding="utf8"))
        assert saved.keys()==result.keys()
        for k,v in result.items():
            if isinstance(v,str):
                assert v==saved[k]
            else:
                assert np.allclose(v,saved[k],rtol=1e-10,atol=1e-14), k
        print("1403 hypercharge--weak scale, uniform volume bound and correction budget verified")
    else:
        print(json.dumps(result,indent=2,allow_nan=False))
