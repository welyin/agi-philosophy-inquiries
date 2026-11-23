"""One-loop physical X pole, nuclear normalization, and held-out strong coupling."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
PI = math.pi
W = np.array([5, -3, -2])
C = np.array([5, 3, 2])
SM = np.array([41/10, -19/6, -7])
D = np.array([2/5, 2/3, 0])
T = np.array([4/15, 0, 2/3])
HC = T/4
S3 = np.array([0, 1/3, 0])
S8 = np.array([0, 0, 1/2])


def read(name):
    return json.loads((HERE/name).read_bytes())


CFG = read("inputs.json")
EW = read(CFG["low_gauge"])
STRONG = read(CFG["lattice_strong"])
FIVES = read(CFG["fixed_split_fives"])
FLAVOUR = read(CFG["flavour_screen"])
NUC = read(CFG["nucleon"])
MZ = EW["mu_ref_GeV"]


def inverse_alpha(g):
    return 4*PI/g**2


def physical_gauges(tag, shifts=(0, 0, 0, 0)):
    def shifted(pair, shift):
        return pair[0]+pair[1]*shift
    a1 = inverse_alpha(shifted(EW["g1_GUT"], shifts[0]))
    a2 = inverse_alpha(shifted(EW["g2"], shifts[1]))
    if tag == "common_MZ":
        a3 = inverse_alpha(shifted(EW["g3"], shifts[2]))
    else:
        five = shifted(STRONG["strong_lattice_5f"]["ALPHA_2026"], shifts[2])
        top = shifted(STRONG["top_MS_matching_GeV"], shifts[3])
        a3 = 1/five+math.log(top/MZ)/(3*PI)
    return np.array([a1, a2, a3])


def maximize(avec, floor):
    log_q = 2*PI*float(W@avec)/44
    x = (log_q+1-math.log(floor/MZ))/21
    q = MZ*math.exp(log_q)
    return log_q, x, q


def run_gauge(avec, mx, fields):
    poles = sorted({MZ, mx, *(mass for mass, _ in fields if MZ < mass < mx)})
    gauges = avec.copy()
    integrals = np.zeros(3)
    for lo, hi in zip(poles, poles[1:]):
        beta = SM.copy()
        for pole, weights in fields:
            if pole <= lo:
                beta += weights
        length = math.log(hi/lo)/(2*PI)
        z = beta*length/gauges
        assert np.all(z < 1)
        for i, bi in enumerate(beta):
            integrals[i] += -math.log1p(-z[i])/bi if bi else length/gauges[i]
        gauges -= beta*length
    return gauges+C/(12*PI), math.exp(float(np.array([23/20, 9/4, 2])@integrals))


def low_strong_holdout(s, ksum, a1, a2):
    # No strong coupling or D/T scale is used in this conditional inversion.
    lam = (s/ksum)**.25
    def u1_upper(t):
        return a1+C[0]/(12*PI)-SM[0]*t/(2*PI)
    def allowed(t):
        return MZ*math.exp(t)*math.sqrt(u1_upper(t)/(4*PI)) >= lam
    lo, hi = 0., 35.
    assert not allowed(lo) and allowed(hi)
    assert u1_upper(hi) > SM[0]/(4*PI)
    for _ in range(100):
        mid = (lo+hi)/2
        if allowed(mid):
            hi = mid
        else:
            lo = mid
    maximum = (5*a1-3*a2)/2-(42*lo-2)/(4*PI)
    return {"Lambda_GeV": lam, "MX_min_GeV": MZ*math.exp(lo),
            "a3_max": maximum, "alpha3_min": 1/maximum if maximum > 0 else None}


def central(tag, ksum):
    a = physical_gauges(tag)
    base = FIVES[tag]["two_split_complete_pairs"]
    flavour = FLAVOUR[tag]
    a0 = base["inverse_alphas_at_matching"][0]
    strength = flavour["ten_pion_e_mu_strength_min"]
    scale = (strength/ksum)**.25
    assert math.isclose(scale, flavour["Lambda_geV_lower"]["one_width_shape_box"], rel_tol=2e-12)
    fields = [(m, wgt) for masses, wgt in [(base["M_D_GeV"], D), (base["M_T_GeV"], T)]
              for m in masses]
    floor_results = []
    for floor in CFG["adjoint_pole_floor_GeV"]:
        t, x, q = maximize(a, floor)
        mx = q*math.exp(x)
        u = math.log(mx/floor)
        assert math.isclose(q, base["matching_scale_GeV"], rel_tol=2e-11)
        assert abs(44*x-2*u-2) < 2e-13
        a5 = a0+3/(12*PI)-((SM[1]+4/3)*x+u/3)/(2*PI)
        matched, a2short = run_gauge(a, mx, fields+[(floor, S3), (floor, S8), (q, HC)])
        assert max(abs(matched-a5)) < 2e-9
        g5 = math.sqrt(4*PI/a5)
        cap = mx/g5
        assert a5 > 25/(8*PI)
        floor_results.append({"floor_GeV": floor, "MX_over_Q_max": math.exp(x),
                              "MX_max_GeV": mx, "m_HC_at_optimum_GeV": q, "a5": a5,
                              "g5": g5, "MX_over_g5_joint_max_GeV": cap,
                              "A2_full_max_to_survive": (cap/scale)**2,
                              "S_max_at_A2_1": cap**4*ksum,
                              "A2short_gauge_only_at_joint_optimum": a2short})
    corner_mx = [qcorner*math.exp(xcorner)
                 for shifts in itertools.product([-2, 2], repeat=4)
                 for _, xcorner, qcorner in [maximize(physical_gauges(tag, shifts), MZ)]]
    r = flavour["MX_Q_lower_if_A2_1_same_g"]
    t = math.log(q/MZ)
    extra = (42*math.log(r)-2*t-2)/(2*PI)
    return {"low_a": a.tolist(), "Q_GeV": q, "flavour_strength_at_1811_gate": strength,
            "conditional_nuclear_scale_GeV": scale, "threshold_optima": floor_results,
            "two_width_box_MX_upper_range_GeV": [min(corner_mx), max(corner_mx)],
            "old_g_reference_r": r, "required_extra_projected_inverse_coupling": extra,
            "minimum_max_abs_inverse_coupling_mismatch": extra/10,
            "single_hypercharge_shift_required_if_fixed_r": -extra/5,
            "required_downward_change_in_low_a3_at_fixed_r": extra/2}


def tidy(obj):
    if isinstance(obj, (float, np.floating)):
        return float(f"{obj:.12g}")
    if isinstance(obj, dict):
        return {k: tidy(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [tidy(i) for i in obj]
    return obj


def run():
    for field, n in [(SM, 44), (C, 12), (D, 0), (T, 0), (HC, 0), (S3, -1), (S8, -1),
                     ([-1, -3, 2], 0), ([7, 9, 4], 0), ([1, 0, 0], 5)]:
        assert abs(W@field-n) < 1e-12
    pion = runpy.run_path(str(HERE/CFG["pion_weights"]))
    factors = [pion["response"](i, True)[0] for i in range(2)]
    ksum = sum(NUC["hbar_GeV_s"]/(NUC["seconds_per_year"]*tau*factor)
               for factor, tau in zip(factors, NUC["lifetime_lower_years_90"]))
    channel = {tag: central(tag, ksum) for tag in FLAVOUR if tag in FIVES}
    # Largest allowed a3 is attained by largest a1, smallest a2:
    # increasing a1 also decreases the minimum mass implied by the nuclear bound.
    best = physical_gauges("common_MZ", (-2, 2, 0, 0))
    nominal = physical_gauges("common_MZ")
    tests = {str(s): {"central_EW": low_strong_holdout(s, ksum, *nominal[:2]),
                       "two_width_EW_most_permissive": low_strong_holdout(s, ksum, *best[:2])}
             for s in CFG["independently_specified_visible_strength_tests"]}
    return tidy({"one_width_pion_K": factors, "sum_width_limit_over_K": ksum,
                 "gauge_and_nuclear": channel, "independent_visible_strength_to_strong": tests})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    content = (json.dumps(run(), indent=2, ensure_ascii=False)+"\n").encode()
    target = HERE/"results.json"
    if args.check:
        assert json.loads(target.read_bytes()) == json.loads(content)
    else:
        target.write_bytes(content)
    print(hashlib.sha256(content).hexdigest()[:12])
