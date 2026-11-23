"""Independent full-box one-loop coloured-doublet bound and exact two-scale optimization."""
from itertools import product
import argparse
import json
import math
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
Q1, Q2 = 753.66, 1600.5
CMS = [(Q1, .0898, .0007, .0008), (Q2, .0821, .0008, .0010)]
N = 4
NORMALIZE = 6*math.pi


def read(rel):
    return json.loads((HERE/rel).read_bytes())


def visible_widths(row):
    _, alpha, low, high = row
    return [1/alpha-1/(alpha+high), 1/(alpha-low)-1/alpha]


def coupled_minimum(geom, mz, upper):
    a, b = math.log(Q1/mz), math.log(Q2/mz)
    d = b-a
    upper_log = math.log(upper/mz)
    budget = N*math.log(geom/mz)
    assert 0<a<b<upper_log and 0<budget<N*b
    full = int(budget//b)
    residual = budget-full*b
    logs = [b]*full + [residual] + [0]*(N-full-1)
    minimum = ((N-full)*d-max(0,residual-a))/NORMALIZE
    # Each cell in the arrangement is a box cut by one half-space.
    # A vertex has all coordinates at knots, or at most one off a knot with tight budget.
    knots = (0,a,b,upper_log)
    vertices = [point for point in product(knots,repeat=N) if sum(point)<=budget+1e-12]
    for fixed in product(knots,repeat=N-1):
        last = budget-sum(fixed)
        if 0<=last<=upper_log:
            vertices.append((*fixed,last))
    def evaluate(point):
        return sum(max(0,b-x)-max(0,a-x) for x in point)/NORMALIZE
    assert abs(evaluate(logs)-minimum)<1e-14
    assert abs(min(map(evaluate,vertices))-minimum)<1e-14
    assert abs(max(map(evaluate,vertices))-N*d/NORMALIZE)<1e-14
    return {"min_increment":minimum, "max_increment":N*d/NORMALIZE,
            "minimizer_doublets_GeV":[mz*math.exp(x) for x in sorted(logs)],
            "cumulative_min_at_Q1":N*max(0,math.log(Q1/geom))/NORMALIZE,
            "cumulative_min_at_Q2":N*max(0,math.log(Q2/geom))/NORMALIZE,
            "value_from_subtracting_cumulative_minima":N*d/NORMALIZE}


def solve():
    old = runpy.run_path(str(HERE.parent / "1663/independent_check.py"))
    cfg = read("../1663/inputs.json")
    measured = read("../1369/inputs.json")
    mz = measured["mu_ref_GeV"]
    upper = cfg["matching_GeV"]
    L = math.log(upper/mz)
    # Recompute X and Y together for each of the same 16 outer-box corners.
    gains = read("../1663/results.json")["effective_gain_box"]
    g = [measured[k] for k in ("g1_GUT","g2","g3")]
    corners = []
    for signs in product((-1,1),repeat=3):
        inverse = [4*math.pi/(v+2*s*dv)**2 for (v,dv),s in zip(g,signs)]
        for gain in gains:
            X,Y = old["solve"](inverse,gain,L)
            corners.append((X,Y,signs,gain))
    X,Y,signs,gain = min(corners,key=lambda c:c[0])
    geom = upper*math.exp(-X/N)
    assert abs(geom-read("../1663/independent_results.json")["four_B_order_statistic_upper_bounds_GeV"][0])<1e-6
    central_gain = read("../1663/independent_results.json")["physical_central_gain"]
    g0 = [4*math.pi/v**2 for v,_ in g]
    X0,Y0 = old["solve"](g0,central_gain,L)
    results = {}
    for label,x,y in (("outer_box",X,Y),("central_point",X0,Y0)):
        geo = upper*math.exp(-x/N)
        ans = coupled_minimum(geo,mz,upper)
        # Use the complete log budget to realize both threshold sums by a triplet at MU
        # and a common sextet far above Q2. Heavy coloured fermions can also exceed Q2.
        masses = ans["minimizer_doublets_GeV"]
        sext = geo*math.exp(y/N)
        assert abs(sum(math.log(upper/m) for m in masses)-x) < 1e-11
        assert abs(sum(math.log(sext/m) for m in masses)-y) < 1e-11
        assert Q2 < sext <= upper
        widths = visible_widths(CMS[1])
        results[label] = {"geo_max_GeV":geo, "partner_masses_Delta_S_GeV":[upper,sext],
                          "high_Q_min_over_inverse_error_down":ans["cumulative_min_at_Q2"]/widths[0],
                          **ans}
    return {"reference":"https://arxiv.org/html/2412.16665v2",
            "table":"https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-24-007/CMS-SMP-24-007_Table_003.pdf",
            "CMS_flavour_number":5, "CMS_RGE_order":5,
            "inverse_widths_down_up_by_Q":{str(row[0]):visible_widths(row) for row in CMS},
            "same_corner":{"X":X,"Y":Y,"g123_signs":list(signs),"k":gain},
            "CMS_inverse_two_centers_difference":1/CMS[1][1]-1/CMS[0][1],
            "scenarios":results}


def stable(v):
    if isinstance(v,dict):
        return {k:stable(x) for k,x in v.items()}
    if isinstance(v,list): return [stable(x) for x in v]
    if isinstance(v,float): return float(f"{v:.11g}")
    return v


if __name__=="__main__":
    arg = argparse.ArgumentParser()
    arg.add_argument("--check",action="store_true")
    args=arg.parse_args()
    out = stable(solve())
    if args.check:
        assert out==read("independent_results.json")
        print("1665 independent PASS: coupled corners, exact increment extrema, both realizable spectra, CMS scheme")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as f:
            f.write(json.dumps(out,indent=2)+"\n")
