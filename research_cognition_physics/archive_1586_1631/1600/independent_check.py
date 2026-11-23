"""Independent two-body weak- and quark-flavor check for round 1600."""
import argparse
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent


def certificate():
    # KLOE, MuLan and FNAL/MILC. Masses are the kinematic working values in GeV.
    kaon, pion, muon, electron = .493677, .13957039, .1056584, .00051099895
    tau_k, dt_k = 12.347e-9, .030e-9
    tau_mu, dt_mu = 2196980.3e-12, 2.2e-12
    form, dform = .9704, .0032
    kinematics = 2*(muon/kaon)**3*((1-(electron/muon)**2)/(1-(pion/kaon)**2))**3
    ratio = kinematics*tau_mu/(tau_k*form**2)
    # Relative rate transfer covers matching/radiative and isospin corrections only
    # when the same low-energy chiral operator survives. No detector refit.
    ratio_low = .9*kinematics*(tau_mu-2*dt_mu)/((tau_k+2*dt_k)*(form+2*dform)**2)
    ratio_high = 1.1*kinematics*(tau_mu+2*dt_mu)/((tau_k-2*dt_k)*(form-2*dform)**2)
    kaon_limits = {"NA62_public_reanalysis_2016_24":1.4e-11, "E787_E949":7.3e-11}
    twist = 58e-6
    transfers = {}
    for key, bound in kaon_limits.items():
        transfers[key] = {"central_mu_conditional_ceiling":ratio*bound,
                          "envelope_mu_conditional_ceiling": ratio_high*bound,
                          "direct_mu_to_conditional_ceiling":twist/(ratio_high*bound)}
    prior = json.loads((HERE.parent/"1598"/"results.json").read_bytes())
    blind = {}
    for branch, outcomes in prior["central_gamma_offsets_degrees"].items():
        modes = {}
        for tag, record in outcomes["tags"].items():
            p = record["down_flavor_fractions"]
            r = [math.sqrt(x) for x in p]
            # Distance to the union of the two coordinate hyperplanes where e-mu vanishes.
            theta = math.asin(min(r[:2]))
            drop = 0 if r[0]<r[1] else 1
            w = r.copy()
            w[drop] = 0.
            norm = math.hypot(*w)
            w = [x/norm for x in w]
            assert math.isclose(sum(x*y for x,y in zip(w,r)),math.cos(theta),rel_tol=1e-13)
            assert w[0]*w[1]==0 and r[0]*r[1]>0
            off = r[0]*r[1]
            s = math.sin(math.radians(.2))
            modes[tag] = {"first_mu_blind_angle_degrees":math.degrees(theta),
                           "tag_phase_independent_ratio_at_0p2deg":[ratio_low*max(0,1-s/off)**2,
                                                                   ratio_high*(1+s/off)**2]}
        blind[branch]=modes
    assert ratio_low < ratio < ratio_high
    return {"BR_mu_per_BR_K_central":ratio, "registered_two_width_and_10pct":[ratio_low,ratio_high],
            "TWIST_90pct_Aminus1_endpoint_standard_normalization":twist,
            "muon_inference":transfers, "alignment_orientations":blind}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=2e-11,abs_tol=1e-22)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=certificate()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        if not same(result,json.loads(path.read_bytes())):
            raise SystemExit("1600 independent result differs")
        print("1600 independent two-body widths, two experimental transfers and six exact misalignment witnesses verified.")
