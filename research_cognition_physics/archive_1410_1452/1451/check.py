"""Match a horizontal gauge anomaly using colored fermions and actual QCD running."""
import argparse
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(file):
    return json.loads((HERE/file).read_text(encoding="utf8"))


def inventory(charges):
    # numbers: color dimension, weak dimension, 6Y, twice Dynkin indices for color and weak
    reps = {"Q":(3,2,1,1,1), "uc":(3,1,-4,1,0), "dc":(3,1,2,1,0),
            "L":(1,2,-3,0,1),"ec":(1,1,6,0,0)}
    light = Counter()
    heavy = Counter()
    weights = 0
    def record(dest, info, x):
        nc,nw,y,tc,tw=info
        dest.update({"SU3_squared_F":nw*tc*x, "SU2_squared_F":nc*tw*x,
                     "Y_squared_F":nc*nw*y*y*x, "Y_F_squared":nc*nw*y*x*x,
                     "F_cubed":nc*nw*x**3, "gravity_F":nc*nw*x})
    for name, values in charges.items():
        nc,nw,y,tc,tw = reps[name]
        for x in values:
            record(light,reps[name],x)
            for k in range(1,x+1):
                record(heavy,reps[name],k-1)
                record(heavy,(nc,nw,-y,tc,tw),-k)
                weights += nw*tc
    assert all(light[key]+heavy[key]==0 for key in light)
    assert light["SU3_squared_F"] == weights
    return {"light":dict(light),"heavy":dict(heavy),"W":weights}


def analyze():
    cfg=read("inputs.json")
    data=read(cfg["mass_boxes"])
    e=cfg["illustrative_epsilon"]
    coeff=cfg["coefficient_budget_B"]
    output={}
    for scale in ("MZ","1e16"):
        lu,ld=data[scale]["log_volume_ranges"][:2]
        net=lu[0]+ld[0]-6*math.log(coeff)
        output[scale]={"epsilon_upper_if_QCD_AF":math.exp(-net/10),
                       "min_integer_D_at_illustrative_epsilon":math.ceil(net/-math.log(e))}
    example=inventory(cfg["quark_charge_control"])
    p=cfg["quark_charge_control"]
    low=read(cfg["light_yukawas_and_ckm"])
    sing=np.concatenate([np.array(low["yukawa"][sec])[:,0]/e**np.array(p[name])
                         for sec,name in [("up","uc"),("down","dc"),("charged_leptons","ec")]])
    assert min(sing)>1/coeff and max(sing)<coeff
    example["bare_core_singular_extrema"]=[min(sing),max(sing)]
    w=example["W"]
    bnew=float(Fraction(7)-Fraction(2,3)*w)
    sample=cfg["cms_R32"]
    q1,q2=sample["bins_GeV"]
    low_a,high_a=sample["alpha_centers"]
    dl,du=sample["working_symmetric_errors"]
    err=sample["difference_in_inverse_alpha_transfer"]
    interval=math.log(q2/q1)
    baseline=7/(2*math.pi)*interval
    screened=bnew/(2*math.pi)*interval
    bins=[]
    for n in (1,2):
        lo=1/(high_a+n*du) - 1/(low_a-n*dl)
        hi=1/(high_a-n*du) - 1/(low_a+n*dl)
        candidate=q2*math.exp(-3*math.pi*(baseline+err-lo)/w)
        bins.append({"box_width":n,"inverse_difference":[lo,hi],
                     "common_threshold_lower_if_above_first_bin":candidate if candidate>q1 else None,
                     "minimum_transfer_if_active_over_full_interval":max(lo-screened,screened-hi,0.)})
    uv=cfg["high_window"]
    end=uv["UV_GeV"]
    max_alpha=uv["alpha_control"]
    hi_results=[]
    for threshold in uv["illustrative_common_thresholds_GeV"]:
        at=1/high_a+7/(2*math.pi)*math.log(threshold/q2)
        at_end=at+bnew/(2*math.pi)*math.log(end/threshold)
        crosses=threshold*math.exp((at-1/max_alpha)*2*math.pi/(-bnew))
        hi_results.append({"common_threshold_GeV":threshold,"alpha_control_crossing_GeV":crosses,
                           "alpha_at_end_if_inside_control":1/at_end if at_end>=1/max_alpha else None})
    assert bins[0]["minimum_transfer_if_active_over_full_interval"]>err
    assert bins[1]["common_threshold_lower_if_above_first_bin"] is None
    return {"mass_only":output,"telescoping_anomalies":example,"one_loop_b3_above_all":bnew,
            "R32_SM_difference":baseline,"R32_if_all_active_difference":screened,
            "R32_sensitivity":bins,"UV_examples":hi_results}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=analyze()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        assert result==read("results.json")
    print("1451 checked",result["one_loop_b3_above_all"])
