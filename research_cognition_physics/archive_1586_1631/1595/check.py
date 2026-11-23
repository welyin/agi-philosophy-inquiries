"""Joint LFV, tau signal and oscillation tests of a common Higgs mass direction."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8-sig"))


def null_column(order, a, b, c, delta):
    s12, s13, s23 = map(math.sqrt,(a,b,c))
    c12, c13, c23 = map(math.sqrt,(1-a,1-b,1-c))
    phase = np.exp(1j*delta)
    if order == "NO":
        return np.array([c12*c13,-s12*c23-c12*s23*s13*phase,
                         s12*s23-c12*c23*s13*phase])
    return np.array([s13/phase,c13*s23,c13*c23])


def two_body(masses, mh, i, j):
    # Exact coefficient of |Yij|^2+|Yji|^2 and maximum interference.
    u, v = (masses[i]/mh)**2,(masses[j]/mh)**2
    phase = math.sqrt((1-u-v)**2-4*u*v)
    return phase*(1-u-v), phase*4*math.sqrt(u*v)


def partial(y, masses, mh, i, j):
    if i == j:
        beta = math.sqrt(1-4*(masses[i]/mh)**2)
        return beta**3*y[i,i].real**2 + beta*y[i,i].imag**2
    factor, interference = two_body(masses,mh,i,j)
    return (factor*(abs(y[i,j])**2+abs(y[j,i])**2)
            -interference*(y[i,j]*y[j,i]).real)


def blind_example(cfg, angles, atmospheric, limits):
    # Real rank-one third source, zero coupling to the measured h.
    n = null_column("NO",angles["NO"]["s12"][1],angles["NO"]["s13"][2],atmospheric["NO"][1],0).real
    ex = cfg["example"]
    masses = cfg["charged_masses_GeV"]
    d = np.diag(masses)
    a = ex["normalized_etau"]
    b = (n[2]-a*n[0])/abs(n[1])
    t = np.array([a,-b,0])/n[2]
    third = np.outer(t,n@d)
    rest = d-third
    assert abs(n@t-1)<2e-15 and np.linalg.matrix_rank(third)==1
    assert np.linalg.matrix_rank(rest)==2 and np.linalg.norm(n@rest)<1e-14
    # Each of two loop sources has unit Higgs/mass gain; the third has zero gain.
    y = rest/masses[2]
    tau_br = ex["true_tau_signal_pb"]/ex["sigma_fid_reference_pb"]
    rates = {key:tau_br*partial(y,masses,cfg["mh_GeV"],i,j)/partial(y,masses,cfg["mh_GeV"],2,2)
             for key,i,j in [("et",0,2),("mt",1,2),("em",0,1)]}
    assert rates["et"]<limits["atlas"]["MC_marginal95"]["et"]
    assert rates["mt"]<limits["atlas"]["MC_marginal95"]["mt"]
    assert rates["em"]<limits["cms"]["em_upper95"]
    return {"n":n.tolist(),"third_left_vector":t.tolist(),"Y_times_v_over_mtau":y.tolist(),
            "Btautau_selected":tau_br,"predicted_B":rates}


def calculate():
    cfg = read("inputs.json")
    angles = read(cfg["angles12_13"])
    atmospheric = read(cfg["angles23"])["theta23_sin2"]
    limits = read(cfg["lfv"])
    data = cfg["tau_measurement"]
    low = data["central_pb"]-data["widths"]*math.hypot(data["stat_pb"],data["syst_down_pb"])
    sigma, error = data["SM_reference_times_BR_pb"]
    fid_upper = (sigma+data["widths"]*error)/data["SM_tau_BR_rounded_lower"]
    kin = []
    for i in [0,1]:
        a,b = two_body(cfg["charged_masses_GeV"],cfg["mh_GeV"],i,2)
        kin.append(a-b/2)
    floor = cfg["two_body_transfer_floor"]
    assert floor<min(kin)
    all_bounds = {}
    for q in cfg["relative_response"]:
        e, mu = [math.sqrt(lim*fid_upper/(q*floor*low))
                 for lim in (limits["atlas"]["MC_marginal95"][key] for key in ["et","mt"])]
        by_order = {}
        for order in ["NO","IO"]:
            a = angles[order]["s12"][1]
            b = angles[order]["s13"][2]
            c = atmospheric[order][1 if order=="NO" else 2]
            n = abs(null_column(order,a,b,c,0))
            cap = (e*n[0]+mu*n[1])/n[2]
            # Critical relative yield transfer below which a zero third response can reenter.
            ans = {"normalized_tau_column_LFV_amplitude_caps":[e,mu],"ratio_radius":cap,
                   "q_required_for_strict_blind_exclusion":q*cap**2,
                   "unmodeled_fraction_needed_if_eta_zero":max(0.,1-cap),
                   "wrong_unit_sign_excluded":bool(cap<2)}
            if cap<1:
                ans.update({"phase_degrees_max":math.degrees(math.asin(cap)),
                            "kappa_modulus_range":[1/(1+cap),1/(1-cap)]})
            by_order[order] = ans
        all_bounds[str(q)] = by_order
    return {"tau_lower_pb":low,"fid_reference_upper_pb":fid_upper,
            "kinematic_width_lower_factors":kin,"bounds":all_bounds,
            "normal_blind_effective_witness":blind_example(cfg,angles,atmospheric,limits)}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(val,b[k]) for k,val in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=3e-11,abs_tol=1e-16)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=calculate()
    if args.check:
        if not same(result,read("results.json")):
            raise SystemExit("mismatch in tau, LFV or null-column checks")
        print("1595: two mass orders, relative complex tau disk, finite mass widths and rank-one blind witness agree")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
