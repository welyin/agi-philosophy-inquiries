"""Additional chiral generations in the soft single- and double-Higgs gluon operators."""
import argparse
import json
import math
from pathlib import Path

BASE = Path(__file__).resolve().parent

def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def response(beta, mixing):
    s, c = math.sin(mixing), math.cos(mixing)
    up = s+c/math.tan(beta)
    down = s-c*math.tan(beta)
    return s,up,down

def solve():
    cfg = read(BASE/"inputs.json")
    vev = read(BASE/cfg["weak"])["v246_GeV"]
    mh = read(BASE/cfg["higgs"])["atlas"]["mass_GeV"][0]
    ymax = math.sqrt(4*math.pi*cfg["perturbative_yukawa_squared_over_four_pi"])
    error = 0.
    for n in range(1,6):
        for beta in (.13,.46,.88,1.21,1.45):
            for mixing in (-1.46,-.61,.25,.96,1.44):
                s,ku,kd = response(beta,mixing)
                l = n*(ku+kd)
                q = n*(ku**2+kd**2)
                predicted = 2*n-2*s*l+l*l/n
                err = abs(q-predicted)/(1+q)
                error = max(error,err)
                assert err < 1e-13
                # Polynomial whose integer root is the generation count.
                assert abs(2*n*n-(q+2*s*l)*n+l*l) < 1e-10
    examples = {}
    for s in cfg["higgs_vector_examples"]:
        beta = math.atan(math.sqrt((1+s)/(1-s)))
        mixing = 2*beta-math.pi/2
        actual,up,down = response(beta,mixing)
        assert abs(actual-s)<1e-13 and abs(up-1)+abs(down+1)<1e-12
        cap = ymax*vev*math.sqrt(1-s)/2
        assert abs(cap-ymax*vev*math.cos(beta)/math.sqrt(2)) < 1e-10
        examples[str(s)] = {"tan_beta":math.tan(beta),"mass_cap_GeV":cap,
                             "min_HH_s_hat_over_four_M_cap_squared":(mh/cap)**2}
    limits = {}
    for k in cfg["expansion_ratio_examples"]:
        # At pair threshold s_hat=4*m_h**2, even before higher invariant mass bins.
        limits[str(k)] = 1-4*mh**2/(k*ymax**2*vev**2)
    assert limits["1"] < .95
    window = cfg["atlas"]["gghh"]
    return dict(vev_GeV=vev, benchmark_higgs_mass_GeV=mh, assumed_yukawa_ceiling=ymax,
                maximum_rounded_algebra_error=round(error,12),
                one_extra_pair_at_exact_soft_single_blind={"c_ggh":0.,"c_gghh":-2/3,
                    "nominal_distance_to_fixed_ATLAS_interval":window[0]+2/3},
                min_Q_at_one_family_given_abs_L_budget={str(e):2-2*e+e*e for e in cfg["soft_one_family_deviation_budget"]},
                exact_blind_examples=examples, necessary_abs_kappaV_ceiling_for_requested_expansion=limits,
                matching="Kinematic ratio below 1 is necessary for a threshold expansion; it does not bound the experimental transfer error.")

if __name__ == "__main__":
    arg = argparse.ArgumentParser()
    opt = arg.add_mutually_exclusive_group(required=True)
    opt.add_argument("--write-new", action="store_true")
    opt.add_argument("--check", action="store_true")
    args = arg.parse_args()
    path = BASE/"results.json"
    result = solve()
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    else:
        assert result == read(path)
        print("1496: doublet projection, integer count and heavy expansion ceilings verified.")
