"""A flavour-only discriminator for a correction-only and a full-mass triplet."""
import argparse
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
POLE = 2000.
# The fixed, small-reactor-angle branch and two printed JUNO marginal widths.
def read(file):
    return json.loads((HERE/file).read_bytes())


def mul(a,b):
    values = [x*y for x,y in itertools.product(a,b)]
    return min(values), max(values)


def div(a,b):
    assert b[0]>0
    return mul(a,(1/b[1],1/b[0]))


def oscillation_box():
    reactor, solar, _ = read("../1659/results.json")["joint_windows"]
    x = [(1-math.sqrt(1-i))/2 for i in reactor]
    clocks = read("../1668/inputs.json")
    d, lo, hi = clocks["juno_NO_dm21"]
    d21 = (d-2*lo,d+2*hi)
    d, lo, hi = clocks["juno_NO_dm31"]
    d31 = (d-2*lo,d+2*hi)
    return x,solar,div(d21,d31),d31


def pole_error():
    # Below the new scalar's 2 TeV matching scale the archived EFT has a
    # single light Higgs and no active sterile columns. Singular masses
    # at the charged endpoint have y_tau<.011; the .02 norm cap is
    # preserved by the differential bound verified by R1729 below 30 TeV.
    leg_log = 1.5*.02**2*math.log(POLE/91.1876)/(16*math.pi**2)
    return math.expm1(4*leg_log)


def ordering_requirement(arm, error):
    x, s, r, gap = oscillation_box()
    # Each inclusive pole lepton fraction can differ from its weak-scale
    # value by ≤ decay transfer. Conservatively allow it for BOTH models.
    e = error+2*pole_error()
    p = [arm[0]-e,arm[0]+e]
    difference = arm[1]-arm[2]
    ratio = div((difference-2*e, difference+2*e), (1-3*p[1],1-3*p[0]))
    e2 = (s[0]*(1-x[1]),s[1]*(1-x[0]))
    den_ratio = (1-3*x[1]+r[0]*(1-3*e2[1]),
                 1-3*x[0]+r[1]*(1-3*e2[0]))
    c = (1-s[1]*(1+x[1]),1-s[0]*(1+x[0]))
    den = (1-x[1]-r[1]*c[1],1-x[0]-r[0]*c[0])
    num = mul(ratio,den_ratio)
    cp = 2*r[1]*math.sqrt(s[1]*(1-s[1])*x[1])
    lo,hi = div((num[0]-cp,num[1]+cp),den)
    t = [(1+lo)/2,(1+hi)/2]
    lower = p[0]*(1+r[0])-x[1]-r[1]*e2[1]
    upper = p[1]*(1+r[1])-x[0]-r[0]*e2[0]
    xmass = div((lower,upper),(1-3*p[1],1-3*p[0]))
    m0 = [math.sqrt(i*j) for i,j in zip(xmass,gap)]
    return t,m0


def signed_score(p):
    return p[1]-p[2]+(1/3-p[0])/4


def no_score_proof():
    x,s,r,_ = oscillation_box()
    tmin = read("../1659/results.json")["joint_windows"][2][0]
    e2_min = s[0]*(1-x[1])
    low3 = (2*tmin-1)*(1-x[0])+(1/3-x[1])/4
    low2 = -11/12 + .75*e2_min
    margin = low3+r[1]*low2
    assert margin>0
    return margin


def inverted_electron_test():
    x,solar,_,_ = oscillation_box()
    # JUNO IO dm31=-2.482(+.025,-.026)e-3 eV^2. Use its smallest
    # two-width magnitude for an upper bound on dm21/abs(dm31).
    d,_,upper = read("../1668/inputs.json")["juno_NO_dm21"]
    rmax = (d+2*upper)/(.002482-2*.000025)
    e2 = solar[0]*(1-x[1])
    margin = 1/3-x[1] + rmax*(e2-1/3)
    assert margin>0
    return margin


def rank_fractions(b):
    h = np.outer(b,b)
    pair = np.abs(h)**2
    p = pair.sum(axis=1)/pair.sum()
    np.testing.assert_allclose(p,abs(b)**2/np.vdot(b,b).real)
    # Unordered same-sign pair on the diagonal, both cross permutations off it.
    rates = np.triu(pair + pair.T,1)+np.diag(np.diag(pair))
    rates /= rates.sum()
    np.testing.assert_allclose((rates.sum(axis=1)+rates.sum(axis=0))/2,p)
    return p


def candidate(b):
    p = rank_fractions(b)
    angle,m0 = ordering_requirement(p,0.)
    limit = .4355 # allow the published T2K .436 to round downward
    low,high = 0.,.045
    assert ordering_requirement(p,high)[0][1]>limit
    for _ in range(42):
        mid = (low+high)/2
        if ordering_requirement(p,mid)[0][1]<limit:
            low=mid
        else:
            high=mid
    assert angle[1]<limit
    # Score has range 2. Bounds use an independent event, not independent
    # lepton legs; the average of its two charges is one bounded observation.
    gap = -signed_score(p)-4*pole_error()
    assert gap>0
    return {"rank_one_leg_fractions_MZ":p.tolist(),"rank_one_score_MZ":signed_score(p),
            "full_source_required_theta23_sq":angle,
            "full_source_required_lightest_mass_NO_eV_from_e_alone":m0,
            "guaranteed_score_separation_at_pole_upper_octant_box":gap,
            "each_leg_fraction_error_radius_vs_T2K_rounding_floor":low,
            "ideal_independent_signal_event_count_Hoeffding_95_one_sided":math.ceil(8*math.log(20)/gap**2)}


def run():
    src = runpy.run_path(str(HERE/"check.py"))
    cases = {}
    for case in ("central_open","registered_open"):
        point = src["endpoint"](case,.20)
        b, c = point["b"], point["C"]
        row = candidate(b)
        solutions = []
        for label, roots in (("historical",src["FROZEN"]),
                             ("two_clocks",src["clocks"](c,b,point["n"]))):
            for rho in roots:
                mass = c + rho*np.outer(b,b)
                squares = np.abs(mass)**2
                full = squares.sum(axis=1)/squares.sum()
                solutions.append({"gain_selection":label, "rho":rho,
                                  "full_source_MZ_probabilities":full.tolist(),
                                  "full_source_MZ_score":signed_score(full)})
        row["same_full_mass_matrix_comparisons"]=solutions
        cases[case] = row
    return {"pole_GeV":POLE,"each_leg_MZ_to_pole_probability_error":pole_error(),
            "NOvA_branch_analytic_positive_numerator_lower":no_score_proof(),
            "inverted_order_electron_excess_over_one_third_numerator_lower":inverted_electron_test(),
            "source_historical_T2K_FC": "https://arxiv.org/pdf/2303.03222#page=35",
            "source_JUNO": "https://arxiv.org/abs/2609.37895",
            "source_doubly_charged_widths":"https://arxiv.org/abs/0805.3536",
            "points":cases}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(run()),parse_float=lambda z:float(f"{float(z):.11g}"))
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==result
        print("1731 verified the inclusive triplet discriminator")
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
