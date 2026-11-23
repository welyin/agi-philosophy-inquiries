"""Rebuild the one-loop slope theorem, threshold witness, and finite-scale audits."""
import argparse
from fractions import Fraction as F
from itertools import combinations_with_replacement
import json
from math import log, pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
IN = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))


def weyl(d3, d2, y):
    return (F(2, 5)*d3*d2*y*y,
            F(d3, 3) if d2 == 2 else F(0),
            F(d2, 3) if d3 == 3 else F(0))


def block(d3, y):
    fields = [(d3, 2, y), (d3, 1, -y-F(1, 2)), (d3, 1, -y+F(1, 2))]
    return tuple(map(sum, zip(*(weyl(*field) for field in fields))))


def excess(b):
    return b[0]-F(9, 5)*b[1]+F(4, 5)*b[2]


def fractions(r):
    return list(map(str, r))


def slope_check():
    checks = 0
    for k in range(-8, 9):
        for sign in [-1, 1]:
            y = F(sign, 6)+k
            assert excess(block(3, y)) == F(8, 5)*(3*k*k+sign*k)
            checks += 1
        assert excess(block(1, F(-1, 2)+k)) == F(8, 5)*k*(k-1)
        checks += 1

    # For total cost <= 1, nonnegative costs imply these finite lists exhaust all blocks.
    quarks = [(s, k) for s in (-1, 1) for k in (-1, 0, 1)
              if 3*k*k+s*k <= 2]
    leptons = [(l, l*(l-1)//2) for l in range(-1, 3)]
    examples = {0: [], 1: []}
    all_min = set()
    for qq in combinations_with_replacement(quarks, 3):
        sum_y3 = sum(F(s, 2)+3*k for s, k in qq)
        dq = sum((3*k*k+s*k)//2 for s, k in qq)
        for ll in combinations_with_replacement(leptons, 3):
            if sum_y3+sum(F(-1, 2)+l for l, _ in ll):
                continue
            cost = dq+sum(d for _, d in ll)
            if cost not in examples:
                continue
            examples[cost].append({"Q_sign_k": qq, "L_l": [l for l, _ in ll]})
            if cost == 0:
                mirrors = sum(s < 0 for s, k in qq)
                assert mirrors == sum(l == 1 for l, d in ll)
                all_min.add(mirrors)
    sm = tuple(3*(a+b) for a,b in zip(block(3, F(1, 6)), block(1, F(-1, 2))))
    assert sorted(all_min) == list(range(4)) and sm == (4,4,4)
    assert examples[1] and any(row["L_l"] == [-1, 0, 1] for row in examples[1])
    gap = F(16, 5)
    return {"checked_block_identities": checks, "SM_matter_B123": fractions(sm),
            "SM_total_b123": fractions([a+F(b) for a,b in zip(sm, IN["conventions"]["bosons_and_Higgs"])]),
            "ordinary_minus_mirror_families_at_SM_slopes": [3-2*m for m in sorted(all_min)],
            "zero_and_first_gap_all_unlabelled_configurations": examples,
            "first_gap_b1": str(gap), "first_gap_inverse_slope_per_decade": float(gap)*log(10)/(2*pi),
            "isolated_one_family_positive_gap": str(F(64,5)),
            "three_equal_copies_positive_gap": str(F(192,5))}


def threshold_check():
    d = IN["split_witness"]
    states = {"ordinary_pair": [(2,F(1,2),d["doublet_log_fraction"]),(1,F(1),d["standard_singlet_one"])],
              "exotic_pair": [(2,F(3,2),d["doublet_log_fraction"]),(1,F(1),d["exotic_singlet_one"]),
                              (1,F(2),d["exotic_singlet_two"])]}
    def integrated(st, end):
        vec = [F(0)]*3
        for d2, charge, expo in st:
            weight = max(F(0), end-(1-F(expo)))
            for i, bi in enumerate(weyl(1,d2,charge)):
                vec[i] += 2*bi*weight
        return vec
    full = [integrated(x, F(1)) for x in states.values()]
    mid = [integrated(x, F(1,2)) for x in states.values()]
    assert full[0] == full[1] and mid[0][0]-mid[1][0] == F(8,25)
    return {"same_two_endpoint_mean_B123": fractions(full[0]),
            "half_interval_integrals_in_full_log_units": list(map(fractions, mid)),
            "absolute_inverse_alpha1_midpoint_difference": float(F(8,25))*log(d["Q1_over_Q0"])/(2*pi),
            "masses_over_Q0": {name:[d["Q1_over_Q0"]**(1-float(F(ex))) for _,_,ex in st] for name,st in states.items()}}


def z_test():
    cfg = IN["conventions"]
    smin, smax = cfg["weak_mixing_sin_squared_bounds"]
    _, upper = cfg["common_Z_normalization_bounds"]
    eps = cfg["lepton_specific_absolute_budget"]
    values = IN["lep_sld"]["absolute_gL"]
    lower = {f:a-2*e for f,(a,e) in values.items()}
    lower_neutral = {f:(v-eps)/upper for f,v in lower.items()}
    assert min(lower_neutral.values()) > smax
    plower = {f:2*(smin+v) for f,v in lower_neutral.items()}
    assert 2 < sum(plower.values()) < 3
    return {"two_marginal_width_projector_diagonal_lower":plower,
            "sum_lower":sum(plower.values()),"required_integer_left_doublet_rank":3,
            "absolute_budget_max_for_strict_sign_test":min(lower.values())-upper*smax}


def opal_check(order):
    cfg = IN["opal"]
    a, b = cfg["spacelike_Q_squared_GeV2"]
    mu, *errors = cfg["measurement"]
    total = sqrt(sum(x*x for x in errors))
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x = (nodes+1)/2
    u = x*(1-x)
    def signal(m):
        return float(cfg["alpha_zero"]*4/pi * np.dot(weights, u*np.log1p((b-a)*u/(m*m+a*u)))/cfg["unit"])
    def minimum_slope(m):
        return float(cfg["alpha_zero"]*4/pi * np.dot(weights, u*(a*u)/(m*m+a*u))*log(b/a)/cfg["unit"])
    def mass_edge(maximum, probe=signal):
        if probe(0) <= maximum:
            return 0.
        lo, hi = 0., 100.
        assert probe(hi)<maximum
        for _ in range(65):
            mid = (lo+hi)/2
            if probe(mid) > maximum:
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2
    predictions = {}
    for err in cfg["SM_errors_for_correlated_and_uncorrelated_hadronic_endpoints"]:
        # Add independent envelopes outward; no confidence interval is assigned.
        tolerance = mu-cfg["SM_prediction"]+2*(total+err)
        predictions[str(err)] = {"max_extra_in_1e-5": tolerance,
                                 "necessary_mass_GeV":mass_edge(tolerance),
                                 "positive_slope_weight_envelope_GeV":mass_edge(tolerance,minimum_slope)}
    positive_only = {str(n): mass_edge(mu+n*total-cfg["known_leptonic_contribution"]) for n in (1,2)}
    return {"combined_measurement_error_in_1e-5":total,"SM_envelopes_two_widths":predictions,
            "positive_hadrons_only_mass_lower_one_two_widths_GeV":positive_only,
            "charge_two_signal_in_1e-5":{str(m):signal(m) for m in [0,1,2,10,100]}}


def calculate():
    result = {"round":1780,"slopes":slope_check(),"thresholds":threshold_check(),
              "Z":z_test(),"OPAL":opal_check(64)}
    opal112 = opal_check(112)
    for key, field in result["OPAL"]["SM_envelopes_two_widths"].items():
        assert abs(field["necessary_mass_GeV"]-opal112["SM_envelopes_two_widths"][key]["necessary_mass_GeV"])<1e-10
    return result


def main():
    opt = argparse.ArgumentParser()
    opt.add_argument("--check", action="store_true")
    out = HERE/"results.json"
    result = json.loads(json.dumps(calculate()))
    if opt.parse_args().check:
        saved = json.loads(out.read_text(encoding="utf-8"))
        assert saved == result
        print("round 1780 verified: exact G6 slopes, thresholds, Z chirality, OPAL charge-two envelope")
    else:
        out.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n",encoding="utf-8")


if __name__ == "__main__":
    main()
