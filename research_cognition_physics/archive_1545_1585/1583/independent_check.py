"""Independent light-flavour certificate in the exactly decoupled 10 mass subcase."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def compact(item):
    if isinstance(item, dict):
        return {k: compact(v) for k,v in item.items()}
    if isinstance(item, list):
        return [compact(x) for x in item]
    if isinstance(item, float):
        return float(format(item, ".11g"))
    return item


def smooth_compressions(p, s):
    # Build all strata with the two free principal angles independently set.
    X = np.diag([p, p, p, s, s])
    count = 0
    for y,z in itertools.combinations_with_replacement([0, .13, .5, .81, 1], 2):
        W = np.zeros((5,3))
        W[0,0] = 1
        W[1,1], W[3,1] = math.sqrt(1-y), math.sqrt(y)
        W[2,2], W[4,2] = math.sqrt(1-z), math.sqrt(z)
        assert np.allclose(W.T@W, np.eye(3))
        assert np.allclose(np.linalg.eigvalsh(W.T@X@W), [p, p+(s-p)*y, p+(s-p)*z])
        count += 1
    return count


def sharp_point(p, s, fraction, ckm_floor):
    # Both physical heavy singular values are one, including the limiting points.
    r = fraction
    mass = np.array([[0,0,0,1,0], [math.sqrt(1-r*r),0,0,0,r]], dtype=float)
    _, heavy, right = np.linalg.svd(mass, full_matrices=True)
    assert np.allclose(heavy, [1,1])
    W = right[2:,:].T
    charges = np.linalg.eigvalsh(W.T@np.diag([p,p,p,s,s])@W)
    span = float(charges[-1]-charges[0])
    assert abs(span-(s-p)*(1-r*r)) < 2e-14
    return {"mass_row_small_ratio":r, "light_charge_span":span,
            "one_of_four_left_vertices_floor":span*ckm_floor}


def kaon_bound(span):
    # Outward masses, width and matched f0 for a reproducible historical holdout.
    kaon, pion, total, form = .493, .140, 5.4e-17, .90
    br_limit = 7.3e-11
    ds_left = span/2
    br_times_F2 = kaon**3*(1-(pion/kaon)**2)**3*form**2*ds_left**2/(64*math.pi*total)
    return {"specific_shape_ds_left":ds_left, "BR_times_F_PQ_squared_min_GeV2":br_times_F2,
            "historical_branching_limit":br_limit, "F_PQ_GeV_if_tree_only":math.sqrt(br_times_F2/br_limit),
            "F_PQ_GeV_if_retained_amplitude_90pct":.9*math.sqrt(br_times_F2/br_limit)}


def result():
    ckm = read(HERE.parent/"1498"/"independent_results.json")["unit_spectral_span_unweighted_floor_full_box"]
    p,s = read(HERE/"inputs.json")["charge"]["base"]
    assert (p,s) == (1,3)
    lattice = [(x,y) for x in range(-60,61) for y in range(-200,201) if 13*x-4*y == 1]
    assert all((x-1)%4 == 0 and y == 3+13*((x-1)//4) for x,y in lattice)
    points = {str(r):sharp_point(p,s,r,ckm) for r in (0,.04,.0625,.25,1)}
    kaon = kaon_bound(points["0.04"]["light_charge_span"])
    return compact({"anomalies":{"K_parent":13*p-4*s,"K_QCD":3*(13*p-4*s),
                                  "E_EM":8*(13*p-4*s),"E_over_K_QCD":8/3},
                    "compression_cases_checked":smooth_compressions(p,s), "sharp_points":points,
                    "specific_45_degree_down_point":kaon})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    r = result()
    target = HERE/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(r,stream,indent=2)
            stream.write("\n")
    elif read(target) != r:
        raise SystemExit("Certificate changed")
    print("1583 independent verified")
