"""Neutral-B oscillation, dimuon branching fractions, and lifetime moments."""
import argparse
import json
from math import isclose, log, pi, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "inputs.json").read_text(encoding="utf-8-sig"))


def box(row, k):
    """Center, symmetric width, or separate down/up statistical and systematic widths."""
    if len(row) == 2:
        low = high = row[1]
    else:
        low, high = row[1] + row[3], row[2] + row[4]
    return [row[0] - k * low, row[0] + k * high]


def positive_ratio(numerators, denominators):
    low = high = 1.
    for left, right in numerators:
        low, high = low * left, high * right
    for left, right in denominators:
        assert left > 0
        low, high = low / right, high / left
    return [low, high]


def S0(x):
    return (4*x - 11*x**2 + x**3) / (4*(1-x)**2) - 3*x**3*log(x) / (2*(1-x)**3)


def calculate():
    c = DATA["const_GeV"]
    sd = DATA["short_distance"]
    osc, tau, bag = DATA["osc_per_ps"], DATA["tau_ps"], DATA["bag"]
    cms = DATA["CMS_s_in_1e9"]
    run3, combined = DATA["LHCb_run3_2024_25"], DATA["LHCb_run123"]
    beta = {q: sqrt(1 - (2*c["m_mu"]/c["M_B" + q])**2) for q in ("s", "d")}
    normalization = 3 * (c["GF"] * c["MW"] * c["m_mu"])**2 / pi**3
    K = normalization * sd["CA_at5GeV"][0]**2 / (sd["S0_published"][0]*sd["etaB"][0])
    phase_ratio = beta["d"] / beta["s"]
    ratio = phase_ratio * bag["s_over_d"][0] * tau["d"][0]/tau["s_heavy"][0] * osc["d"][0]/osc["s"][0]
    absolute = {q: K * beta[q] * osc[q][0] * tau[clock][0] / bag[q][0]
                for q,clock in (("s","s_heavy"),("d","d"))}
    assert isclose(absolute["s"]*1e9,3.81,abs_tol=.005)
    assert isclose(absolute["d"]*1e10,1.031,abs_tol=.001)
    assert abs(S0(sd["xt"][0])-sd["S0_published"][0]) < .001
    assert abs((sd["mtMS"][0]/c["MW"])**2 - sd["xt"][0]) < .001
    results = {}
    delta = DATA["local_correction_trial"]
    for k in (1,2):
        r = positive_ratio([[phase_ratio,phase_ratio],box(bag["s_over_d"],k),
                            box(tau["d"],k),box(osc["d"],k)],
                           [box(tau["s_heavy"],k),box(osc["s"],k)])
        source = box(cms,k)
        pred = [10 * a*b for a,b in zip(source,r)]
        direct = box(combined["d_over_s"],k)
        shape = box(run3["d_in_1e10"],k)
        rho_cross = positive_ratio([shape],[pred])
        rho_direct = positive_ratio([direct],[r])
        triple_floor = DATA["competition_double_ratio"] * r[0] * delta[0]
        assert triple_floor > direct[1]
        # A joint distribution is not needed for a deterministic enclosure of positive ratios.
        results[str(k)] = {"d_over_s_theory":r,"d_prediction_in_1e10":pred,
                           "d_Run3_box_in_1e10":shape,"LHCb_combined_ratio_box":direct,
                           "cross_sample_rate_modifier":rho_cross,"LHCb_ratio_rate_modifier":rho_direct,
                           "threefold_rate_ratio_floor_with_trial_correction":triple_floor,
                           "LHCb_ratio_upper":direct[1]}
    clock = DATA["LHCb_JpsiKK_2011_12_ps_inv"]
    width, splitting = clock["Gamma_s"][0], clock["deltaGamma_s"][0]
    heavy, light = 1/(width-splitting/2), 1/(width+splitting/2)
    shapes = []
    for weight in (0,.25,.5,.75,1):
        mean = weight*heavy + (1-weight)*light
        eff = (weight*heavy**2+(1-weight)*light**2)/mean
        assert isclose(1/mean,(heavy+light-eff)/(heavy*light),rel_tol=1e-12)
        shapes.append({"heavy_weight":weight,"integral_clock_ps":mean,
                       "effective_lifetime_ps":eff,"same_s_rate_d_prediction_factor":heavy/mean})
    return {"short_distance_S0_recalculated":S0(sd["xt"][0]),"K_SM":K,"phase_space_ratio_d_s":phase_ratio,
            "absolute_SM_s_and_d":absolute,"d_over_s_theory":ratio,
            "s_calibrated_K":cms[0]*1e-9*bag["s"][0]/(osc["s"][0]*tau["s_heavy"][0]*beta["s"]),
            "CMS_s_to_d_in_1e10":10*cms[0]*ratio,"Run3_cross_experiment_modifier":run3["d_in_1e10"][0]/(10*cms[0]*ratio),
            "combined_direct_modifier":combined["d_over_s"][0]/ratio,"boxes":results,
            "independent_clock_moment_examples":shapes}


def agrees(a,b):
    if isinstance(a,dict):
        return a.keys() == b.keys() and all(agrees(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(agrees(x,y) for x,y in zip(a,b))
    return isclose(a,b,rel_tol=1e-10,abs_tol=1e-23)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    opts = parser.add_mutually_exclusive_group(required=True)
    opts.add_argument("--check",action="store_true")
    opts.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    result = calculate()
    path = HERE / "results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as dest:
            json.dump(result,dest,indent=2,allow_nan=False)
            dest.write("\n")
    else:
        assert agrees(result,json.loads(path.read_text(encoding="utf-8-sig")))
        print("1404 oscillation, absolute/double-ratio and two-moment checks passed")
