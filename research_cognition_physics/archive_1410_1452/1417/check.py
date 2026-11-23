"""One-loop canonical gauge matching with arbitrary full scalar-doublet thresholds."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
B = np.array([41/10, -19/6, -7.0])
DOUBLE = np.array([1/10, 1/6, 0.0])
TAU = 2 * math.pi


def colour_prediction(a1, a2, effective):
    q = (23 + effective)*2.5/(109 - effective)
    inverse_colour = a2 - q*(a1-a2)
    return 1/inverse_colour, q


def evaluate(gauge, ref, ceiling, cuts):
    a = 4*math.pi/np.square(gauge)
    x, y = a[0]-a[1], a[1]-a[2]
    t, u = (5*x+2*y)/44, (218*y-115*x)/44
    assert t > 0 and u > 0
    effective = u/t
    count = math.ceil(effective)
    unification = ref*math.exp(TAU*t)
    assert unification < ceiling
    inverse_at_match = a - B*t - DOUBLE*u
    assert np.ptp(inverse_at_match) < 1e-12 and inverse_at_match.min() > 0
    endpoints = [colour_prediction(a[0], a[1], h)[0] for h in [0,1]]
    _, q = colour_prediction(a[0], a[1], 1)
    residual = y-q*x
    spread = residual/(1+q)
    # Equal a1 and a3 matching shifts at the minimum. Actual positive t changes.
    t_split = (x+spread)/(B[0]-B[1]+DOUBLE[0]-DOUBLE[1])
    uv = a-(B+DOUBLE)*t_split
    assert abs(uv[1]-uv[0]-spread) < 1e-12 and abs(uv[2]-uv[0]) < 1e-12
    log_mass_budget = TAU*(count*t-u)
    common_mass = ref*math.exp(log_mass_budget/count)
    assert ref <= common_mass <= unification and 0 < t_split < math.log(ceiling/ref)/TAU
    # Independent direct evolution over unequal threshold segments at same log budget.
    logs = np.linspace(0, 2*log_mass_budget/count, count)
    assert logs[-1] < TAU*t
    knots = sorted(set([0.0, *logs, TAU*t]))
    running = a.copy()
    for left, right in zip(knots, knots[1:]):
        live = np.count_nonzero(logs <= left)
        running -= (right-left)*(B+live*DOUBLE)/TAU
    assert np.max(abs(running-inverse_at_match)) < 1e-12
    return {
        "observed_matched_alpha3":1/a[2], "one_extra_predicted_alpha3_min":endpoints[0],
        "one_extra_predicted_alpha3_max":endpoints[1], "effective_extras":effective,
        "minimum_extra_count":count, "M_X_GeV":unification, "unified_inverse_alpha":inverse_at_match[0],
        "log_mass_budget":log_mass_budget, "minimal_count_geomean_mass_GeV":common_mass,
        "one_extra_min_colour_only_inverse_matching_shift_abs":residual,
        "one_extra_min_full_inverse_matching_spread":spread,
        **{f"minimum_extras_strictly_below_{int(cut)}_GeV_at_min_count":
           max(0,count-math.floor(log_mass_budget/math.log(cut/ref))) for cut in cuts},
    }


def calculate():
    conf = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    inherited = json.loads((HERE/conf["inherited"]).read_text(encoding="utf-8"))
    base = np.array([inherited[k] for k in ["g1_GUT","g2","g3"]])
    args = inherited["mu_ref_GeV"], conf["ceiling_GeV"], conf["counting_scales_GeV"]
    boxes = {}
    for width in conf["marginal_widths"]:
        runs = [evaluate(base[:,0]+width*np.array(sign)*base[:,1],*args)
                for sign in itertools.product([-1,1], repeat=3)]
        boxes[str(width)] = {k:[min(row[k] for row in runs),max(row[k] for row in runs)] for k in runs[0]}
        assert boxes[str(width)]["one_extra_predicted_alpha3_max"][1] < boxes[str(width)]["observed_matched_alpha3"][0]
        assert boxes[str(width)]["minimum_extra_count"] == [7,7]
    return {"round":1417,"central":evaluate(base[:,0],*args),"boxes":boxes}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    actions=parser.add_mutually_exclusive_group()
    actions.add_argument("--check",action="store_true")
    actions.add_argument("--write",action="store_true")
    args=parser.parse_args()
    out=json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"results.json"
    if args.write:
        dest.write_text(out,encoding="utf-8")
    elif json.loads(dest.read_text(encoding="utf-8")) != json.loads(out):
        raise SystemExit("Saved results differ")
    print("1417: calibration, holdout, threshold feasibility, sharp split and unequal-threshold replay verified.")
