"""Two CP-even Higgs poles with one top Yukawa direction. No fit to unseen states."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(file):
    return json.loads((HERE / file).read_text(encoding="utf-8-sig"))


def interval(row, widths):
    value, low, high = row
    return value-widths*low, value+widths*high


def angle_at_top(target, z, lower, upper):
    def top(angle):
        return math.cos(angle)+z*math.sin(angle)
    if target <= top(lower):
        return lower
    if target >= top(upper):
        return upper
    for _ in range(75):
        mid = (lower+upper)/2
        if top(mid) < target:
            lower = mid
        else:
            upper = mid
    return (lower+upper)/2


def constraints(p, data, with_transfer):
    mult = float(with_transfer)
    transfer = data["transfers"]
    z_mz = math.sqrt(1/p**2-1)
    z = z_mz * (1-mult*(1-transfer["z_ratio_floor"]))
    widths = data["widths"]
    vector_min = interval(data["ATLAS"]["vector"],widths)[0]-mult*transfer["absolute_light_vector"]
    t_lo, t_hi = interval(data["ATLAS"]["top"],widths)
    t_lo -= mult*transfer["absolute_light_top"]
    t_hi += mult*transfer["absolute_light_top"]
    opening = math.acos(vector_min)
    assert opening < math.atan(z)
    lo = angle_at_top(t_lo,z,-opening,0)
    hi = angle_at_top(t_hi,z,0,opening)
    assert abs(math.cos(hi)+z*math.sin(hi)-t_hi)<1e-10
    top_H = math.sqrt(1+z*z-t_hi*t_hi)
    c_min, c_max = interval(data["CMS"]["tg"],widths)
    g_min, g_max = interval(data["CMS"]["Zg"],widths)
    ratio_min, ratio_max = c_min/g_max, c_max/g_min
    ratio_min -= mult*transfer["absolute_CMS_ratio"]
    ratio_max += mult*transfer["absolute_CMS_ratio"]
    assert 0<ratio_min<ratio_max<1+z*z and z*z>1/8
    dr = ratio_max-1
    cms_top = (z*z-dr)/math.sqrt(z*z+dr*dr)
    min_vec = math.cos(max(abs(lo),abs(hi)))
    return {"z_min":z,"ATLAS_top_interval":[t_lo,t_hi],"ATLAS_vector_floor_input":vector_min,
            "ATLAS_angle_degrees":[math.degrees(lo),math.degrees(hi)],
            "ATLAS_vector_min_geometry":min_vec,"ATLAS_abs_top_H_floor":top_H,
            "ATLAS_Htt_squared_floor":top_H**2,
            "CMS_abs_ratio_box":[ratio_min,ratio_max],"CMS_abs_top_H_floor":cms_top,
            "CMS_Htt_squared_floor":cms_top**2}


def pole_examples(data):
    cfg = data["aligned_examples"]
    gauge = load(data["gauge"])
    initial = [gauge[key][0] for key in ("g1_GUT","g2","g3")]
    b = [21/5,-3,-7]
    c = [17/20,9/4,8]
    K = 8*math.pi**2
    rows = []
    # Central two-doublet example, top-dominant, no extra signed response.
    for cutoff, mass in zip(cfg["cutoffs_GeV"], cfg["masses_GeV"]):
        count = 65536
        step = math.log(cutoff/gauge["mu_ref_GeV"])/count
        def kernel(j):
            t = (j+.5)*step
            return math.exp(sum(c_i/b_i*math.log1p(-b_i*g_i**2*t/K)
                                for b_i,c_i,g_i in zip(b,c,initial)))
        J = (4.5/K)*step*math.fsum(kernel(j) for j in range(count))
        comp = 1/math.sqrt(J)
        p = cfg["top_effective_low"]/comp
        z = math.sqrt(1/p**2-1)
        rows.append({"cutoff_GeV":cutoff,"Y_comp_MZ":comp,"p":p,
                     "light_V_top":[1.,1.],"other_V_top":[0.,z],
                     "other_neutral_charged_mass_GeV":mass,
                     "h_charged_trilinear_tree":0.,
                     "lambda_v":cfg["mh_GeV"]**2/(2*cfg["v_GeV"]**2)})
    return rows


def calculate():
    data = load("inputs.json")
    p = load(data["previous"])["second_doublet"]["vev_share_upper"]
    rows = pole_examples(data)
    assert all(r["p"]<p and r["lambda_v"]>0 for r in rows)
    assert rows[0]["other_V_top"][1]>rows[1]["other_V_top"][1]
    return {"bare_projection_cap":p,"nominal":constraints(p,data,False),
            "registered_transfers":constraints(p,data,True),
            "alignment_witnesses":rows,
            "test_example":{"a_interval":[.94,.96],"t_interval":[.96,1.04],
                            "max_abs_t_minus_a":.10,
                            "min_required_difference":constraints(p,data,True)["z_min"]*math.sqrt(1-.96**2)}}


def equal(a,b):
    if isinstance(a,dict):
        return set(a)==set(b) and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.check:
        if not equal(result,load("results.json")):
            raise SystemExit("1588: archived results differ")
        print("1588: two independent data projections and bounded alignment witnesses verified")
    else:
        (HERE/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
