"""1419: independent chiral norm with proton/Cs carried into both LEP helicities."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(file):
    return json.loads((HERE / file).read_text(encoding="utf-8"))


def distance_from_zero(a, b):
    lo, hi = sorted((a, b))
    return max(lo, -hi, 0.)


def compute():
    cfg = read("inputs.json")
    inputs = cfg["low_energy"]
    lep = cfg["LEP"]
    pv = read(inputs["PV"])
    pi = read(inputs["pion"])["PIENU"]
    gf = read(inputs["GF"])["const_GeV"]["GF"]
    v2 = 1 / (math.sqrt(2)*gf)
    wlow = inputs["sin2_low_illustrative_tree_1418"]
    w = lep["hard_inputs"]["sin2_weak"]
    e2 = 4*math.pi/lep["hard_inputs"]["alpha_em_inverse"]
    smin = (lep["sqrt_s_GeV"]*lep["sqrt_sprime_over_s_cut"])**2
    mz = lep["hard_inputs"]["MZ_GeV"]
    p = pv["proton"]
    cs = pv["caesium"]
    dp = p["pv_electron_scattering_only"][0]-p["SM"][0]
    dc = cs["value_and_errors_experiment_atom"][0]-cs["SM"][0]
    wp = 2*math.hypot(p["pv_electron_scattering_only"][1],p["SM"][1])
    wc = 2*math.hypot(*cs["value_and_errors_experiment_atom"][1:],cs["SM"][1])
    c1 = []
    for i,j in itertools.product((-1,1), repeat=2):
        # 2u+d = P ; 188u+211d = Cs; determinant 234.
        P, Cs = -(dp+i*wp)/2, -(dc+j*wc)/2
        u = (211*P-Cs)/234
        c1.append((u,P-2*u))
    seen, ep, sp = pi["exp_in_units_1e_4"]
    theory, et = pi["same_paper_SM"]
    rad = 2*math.hypot(ep,sp)
    a = math.sqrt((seen-rad)/(theory+2*et))
    b = math.sqrt((seen+rad)/(theory-2*et))
    measurement = lep["sigma_had_pb_and_errors"]
    upper = measurement[0]+2*math.hypot(*measurement[1:])
    unit = 3*smin*lep["normalization_pb_per_GeV_minus2"]/(96*math.pi)
    rows = []
    for cc in [(a-1,b-1),(-1-b,-1-a)]:
        sq = []
        for chirality in range(2):
            for i,(charge,iso) in enumerate(((2/3,.5),(-1/3,-.5))):
                flip = -2*(iso-2*wlow*charge)
                samples = [flip-2*iso*d-(vec[i] if chirality else 0)
                           for d in cc for vec in c1]
                minimum = distance_from_zero(min(samples),max(samples))
                ge = w if chirality else w-.5
                sm_max = (2*e2*abs(charge)/smin
                          + abs(e2*ge*(iso-2*w*charge)/(w*(1-w)))/(smin-mz*mz))
                sq.append(max(0.,2*minimum/v2-sm_max)**2)
        floors = [unit*sum(sq[:2]),unit*sum(sq[2:])]
        rows.append({"left_and_right_pb_floors":floors,
                     "combined_pb_floor":sum(floors),
                     "minimum_positive_transfer_weight_for_exclusion":upper/sum(floors)})
    return {"round":1419,"two_width_LEP_upper_pb":upper,"positive_and_reflected_CC":rows}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    group = p.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = p.parse_args()
    result = compute()
    out = HERE / "independent_results.json"
    if args.write:
        with out.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    else:
        saved = json.loads(out.read_text(encoding="utf-8"))
        for old,new in zip(saved["positive_and_reflected_CC"], result["positive_and_reflected_CC"], strict=True):
            assert all(math.isclose(a,b,rel_tol=1e-10) for a,b in zip(old["left_and_right_pb_floors"],new["left_and_right_pb_floors"], strict=True))
            for key in ("combined_pb_floor","minimum_positive_transfer_weight_for_exclusion"):
                assert math.isclose(old[key],new[key],rel_tol=1e-10)
        assert saved["round"]==result["round"]
        assert math.isclose(saved["two_width_LEP_upper_pb"],result["two_width_LEP_upper_pb"],rel_tol=1e-10)
    print(json.dumps(result))
