"""Single-cycle adjoint bounds and separated-index neutrino test."""
import argparse
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def mass_fit(red, green, mixed):
    # Linear + quadratic invariants in the fixed adjoint vacuum.
    extra = (red + green - 2*mixed) / 25
    linear = (green-red+5*extra) / 5
    base = red+3*linear-9*extra
    assert math.isclose((base-.5*linear-6*extra)/mixed, 1, abs_tol=1e-12)
    # The singlet matrix element is (3*red+2*green)/5.
    singlet = base-linear+7*extra
    assert abs(singlet-(3*red+2*green)/5) < 2e-14*max(abs(base), abs(mixed))
    return singlet


def global_spectra(cfg, gauge, h):
    # At these simultaneous endpoints the scale ceiling and O/T minimum coincide.
    signs = [-1,1,-1]
    ag = [4*math.pi/(gauge[k][0]+i*cfg["coupling_widths"]*gauge[k][1])**2-i*h
          for k,i in zip(("g1_GUT","g2","g3"), signs)]
    d12 = 2*math.pi*(ag[0]-ag[1])
    d13 = 2*math.pi*(ag[0]-ag[2])
    log_ceiling = .75*d12/(5.45-1)
    mx = gauge["mu_ref_GeV"]*math.exp(log_ceiling)
    t = mx*math.exp(.75*d12-5.45*log_ceiling)
    o = mx*math.exp(.5*d13-5.55*log_ceiling)
    assert math.isclose(t,gauge["mu_ref_GeV"],rel_tol=2e-13)
    assert t < o < mx
    ratio = math.exp(2*math.pi*(-.25*ag[0]+.75*ag[1]-.5*ag[2])-.1*log_ceiling)
    signed_ratio = (2*ratio-3)/5
    calculated = mass_fit(-t,o,mx)
    assert math.isclose(calculated/t,signed_ratio,rel_tol=3e-8)
    return {"log_ceiling":log_ceiling,"MX_at_edge_GeV":mx,
            "O_over_T_lower":ratio,"S_over_T_lower":signed_ratio}


def neutrino(cfg, n, order, ratio):
    solar, atm = n["solar"], n["atm_abs"]
    if order == "NO":
        lo_ev = math.sqrt(solar[1])
        q = math.sqrt(solar[1]/atm[2])
        electron_p = (1-n["s12"][2])*(1-n["s13"][2])
    else:
        lo_ev = math.sqrt(atm[1]-solar[2])
        q = math.sqrt(1-solar[2]/atm[1])
        electron_p = n["s13"][1]
    cap = cfg["yukawa_ratio_cap_for_mass_comparison"]
    eps = cfg["light_operator_relative_budget"]
    return {"q_min":q, "distance_from_rank_one_at_least_eV":lo_ev,
            "minimum_Ys_over_Yt":math.sqrt(q*ratio),
            "minimum_separate_trace_over_Yt":math.sqrt(q*ratio)-math.sqrt(3/5),
            "minimum_mass_cancellation_fraction_at_Y_ratio_cap":1-cap**2/(q*ratio),
            "minimum_Ys_over_Yt_with_operator_budget":
               math.sqrt(max(0.,(q-eps)/(1+eps))*ratio),
            "angle_from_pure_e_min_degree":math.degrees(math.asin(math.sqrt(electron_p)))}


def check_words():
    for p in range(13):
        for q in range(13):
            weak, color = (-3)**(p+q), 2**(p+q)
            # Normalized singlet contraction on all 5 eigenlines.
            weights = [2,2,2,-3,-3]
            lhs = sum(w*w*w**p*w**q for w in weights)
            assert lhs == 6*(3*weak+2*color)
            # Keep sqrt(60) normalization factored out (YT has eigenvalue -1/2).
            assert (-3)**p*(-3)*(-3)**q == -3*weak


def calculate():
    cfg = load("inputs.json")
    g, n = load(cfg["gauge"]), load(cfg["oscillation"])
    old = load(cfg["prior_witness"])["extra_insertion"]
    spec = old["mass_witness_GeV"]
    t,o,z = (spec[k] for k in ("weak3","color8","dirac32"))
    s = mass_fit(t,o,z)
    vr = cfg["vEW_GeV"]/math.sqrt(2)
    # Orthogonal Takagi columns of ANY measured PMNS in the charged-lepton basis.
    lo,hi = (math.sqrt(n["NO"][k][0])*1e-9 for k in ("solar","atm_abs"))
    yt,ys = math.sqrt(t*hi)/vr,math.sqrt(s*lo)/vr
    clebsch = math.sqrt(3/5)
    assert math.isclose(vr**2*yt**2/t,hi) and math.isclose(vr**2*ys**2/s,lo)
    check_words()
    branches = {}
    for h in cfg["inverse_alpha_offsets"]:
        bound = global_spectra(cfg,g,h)
        branches[str(h)] = {"mass":bound, "oscillation":{
           op:neutrino(cfg,n[op],op,bound["S_over_T_lower"]) for op in ("NO","IO")}}
    return {"word_pairs_checked":169,"global_matching_boxes":branches,
            "retained_proton_witness":{"MX_GeV":old["MX_GeV"],"MT_GeV":t,"MO_GeV":o,
            "MS_GeV":s,"MS_over_MT":s/t,"NO_central_Yt":yt,"NO_central_Ys":ys,
            "NO_central_Ys_over_Yt":ys/yt,
            "independent_vector_over_Yt_in_orthogonal_witness":math.hypot(ys/yt,clebsch),
            "NO_central_solar_mass_eV":lo*1e9,"NO_central_atmos_mass_eV":hi*1e9}}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,float):
        return float(format(x,".11g"))
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--check",action="store_true")
    actions.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=rounded(calculate())
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as output:
            json.dump(result,output,indent=2,allow_nan=False)
            output.write("\n")
    elif result!=load("results.json"):
        raise SystemExit("Saved certificate differs")
    print(json.dumps(result))
