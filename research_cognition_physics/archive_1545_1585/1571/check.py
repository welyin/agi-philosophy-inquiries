"""Intersect the global charge lattice with the weak-angle data using independent inverse equations."""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))

def level(a, c, angle, beta):
    by,bw,bc=beta
    t=(a*angle-c)/(bw-bc)
    return (a*(1-angle)-by*t)/(c-bc*t)

def angle(a,c,k,beta):
    by,bw,bc=beta
    t=(a-(1+k)*c)/(by+bw-(1+k)*bc)
    return (c+(bw-bc)*t)/a

def min_winding(k):
    half_norm=F(7,12)+k/4
    n=1
    while (n*n*half_norm).denominator != 1:
        n+=1
    return n

def compute():
    cfg=read("inputs.json")
    inp=read(cfg["upstream"])
    a,da=inp["A_em_MSbar"]
    aerr=2*da+inp["extra_A_for_scheme"]+inp["extra_A_for_alpha_s"]
    a_box=[a-aerr,a+aerr]
    g,dg=inp["alpha_s_lattice"]
    c_box=[1/(g+2*dg),1/(g-2*dg)]
    s,ds=cfg["SLD_polarized_leptonic_angle"]
    sld=[s-2*ds,s+2*ds]
    shift=cfg["effective_to_MSbar_and_other_residual_sensitivity"]
    beta=inp["SM_beta_Y_2_3"]
    k_box=[level(ai,ci,u,beta) for ai,ci,u in
           product(a_box,c_box,[sld[0]-shift,sld[1]+shift])]
    assert F(11,9)<min(k_box)<max(k_box)<F(4,3)
    levels={}
    for k in (F(11,9),F(4,3),F(5,3),F(5,4),F(32,25)):
        values=[angle(ai,ci,float(k),beta) for ai,ci in product(a_box,c_box)]
        distance=max(0,min(values)-sld[1],sld[0]-max(values))
        n=min_winding(k)
        quot=[p for p in cfg["standard_quotients"] if (6//p)%n==0]
        levels[str(k)]={"minimum_winding_from_norm":n,"unexcluded_quotients_by_norm":quot,
                        "one_loop_MSbar_angle":[min(values),max(values)],
                        "minimum_total_angle_correction_to_SLD_two_widths":distance}
    up=read(cfg["upstream_results"])["boxes"]["2"]["allowed_k"]
    return {"SLD_only_continuous_level_box":[min(k_box),max(k_box)],
            "SLD_two_widths_effective_angle":sld,"levels":levels,
            "prior_low_energy_inverse_level_boxes":up}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=compute()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as dest:
            dest.write(json.dumps(result,indent=2)+"\n")
    else:
        assert result==read("results.json")
        print("1571 continuous pole window and all integral hypercharge levels checked")
