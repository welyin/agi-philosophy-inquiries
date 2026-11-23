"""Independent flavor-vector and neutral-kaon propagation check."""
import argparse
import cmath
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def avg(values):
    return sum(values) / len(values)


def kaon_rates(background, shift, epsilon):
    amplitude = [z - shift for z in background]
    neutral = [abs(1j*z.imag + epsilon*z.real)**2 / (1+abs(epsilon)**2)
               for z in amplitude]
    return [avg([abs(z)**2 for z in amplitude]), avg(neutral)]


def evaluate():
    # Archive-independent extraction of the intervals quoted by the experiments and the authors.
    plus = (9.6-2*1.8, 9.6+2*1.9)
    charged_pref = (5.173-2*.025)*.997
    neutral_pref = 22.31+2*.13
    gn_upper = plus[1]*neutral_pref/charged_pref
    h = math.pi*.23116*127.9/.225**5
    slots = [complex(-1.30,.39),complex(-1.29,.38),complex(-1.21,.37)]
    mean = avg(slots)
    offsets = [z-mean for z in slots]
    eps = .0023*cmath.exp(.76j)  # stress example; the contraction proof is independent of epsilon
    first = kaon_rates(slots,0,eps)
    mirrored = kaon_rates(slots,2*mean,eps)
    decomp = [x+y for x,y in zip(kaon_rates([mean],0,eps),kaon_rates(offsets,0,eps))]
    for a,b,c in zip(first,mirrored,decomp):
        assert abs(a-b)<2e-15 and abs(a-c)<2e-15
    assert first[1] <= first[0] and gn_upper < 220
    fu=.010042538351926054
    obs=math.sqrt(plus[1]/charged_pref)
    # Conservative independent upper of top and charm terms in a three-angle rectangle.
    st,ca,cb=.24,.05,.005
    bound_top=(st*ca+cb)*(ca+st*cb)*(1.475+.016)/.225**5
    bound_charm=(st+ca*cb)*(.405+.048)/.225
    return {"neutral_max_1e11":gn_upper, "limit_div_max":220/gn_upper,
            "h":h, "broad_W_product":fu*(obs+bound_top+bound_charm)/h,
            "smallest_extra_background":h*3.034414244877883e-7/fu-obs-bound_top-bound_charm,
            "nonuniform_model_rates":first,"indistinguishable_F":[(2*mean/h).real,(2*mean/h).imag],
            "positive_flavor_variance":avg([abs(z)**2 for z in offsets])}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    results=evaluate()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as file:
            json.dump(results,file,indent=2)
            file.write("\n")
    else:
        assert json.loads(path.read_text(encoding="utf-8"))==results
        print("1461 independent OK")
