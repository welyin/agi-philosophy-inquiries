"""Test a coherent two-tag filter, calibration tolerance, and the matched axis."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT=Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT/name).read_text(encoding="utf-8-sig"))


def two_width(row):
    w=2*sum(row[1:])
    return row[0]-w,row[0]+w


def solve():
    cfg=read("inputs.json")
    q=read(cfg["quark_inputs"])
    b=read(cfg["belle"])["BelleII_Vcb_BGL"]
    ratio=q["ub_over_cb"]
    cb=b[0]
    def angle(v,r):
        return math.asin(v*math.hypot(1,r))
    theta=angle(cb,ratio[0])
    theta_box=[angle(tv,rv) for tv,rv in zip(two_width(b),two_width(ratio))]
    m=math.ceil(math.pi/(4*theta))
    cos_m=math.cos(2*m*theta)
    p=-cos_m/(1-cos_m)
    assert math.cos(2*(m-1)*theta)>0 and 0<=p<=1
    residual=max(abs(p+(1-p)*math.cos(2*m*angle_i)) for angle_i in theta_box)
    assert 0<2*m*theta_box[0] and 2*m*theta_box[1]<math.pi
    # Changing to a general complex basis leaves the logical subspace and its partner intact.
    rng=np.random.default_rng(1490)
    U=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    t=U[:,0]
    b=U[:,0]*math.cos(theta)+U[:,1]*math.sin(theta)
    n=U[:,2]
    identity=np.eye(3)
    def ray(x):
        return np.outer(x,x.conj())
    unitary=(identity-2*ray(t))@(identity-2*ray(b))
    powers=np.linalg.matrix_power(unitary,m)
    K=p*identity+(1-p)*(powers+powers.conj().T)/2
    assert np.linalg.norm(K-ray(n))<1e-12
    vec=rng.normal(size=(3,2))+1j*rng.normal(size=(3,2))
    vec/=np.linalg.norm(vec)
    assert np.linalg.norm(K@vec-n[:,None]*(n.conj()@vec)[None,:])<1e-13
    low,high=theta_box
    c=math.cos(theta)
    alt=math.ceil((1+math.log(.01)/math.log(c))/2)
    true_prob=np.linalg.norm(n.conj()@vec)**2
    return {"theta_degrees":math.degrees(theta),"working_theta_degrees":list(map(math.degrees,theta_box)),
            "max_power_per_path":m,"number_of_reflections":2*m,"zero_power_path_probability":p,
            "worst_filter_residual_at_central_setting":residual,
            "sufficient_theta_calibration_deg_for_1pct":math.degrees(.01/(2*m*(1-p))),
            "alternating_measure_cycles_for_1pct":alt,
            "random_entangled_success_probability":true_prob}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=solve()
    path=ROOT/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        assert result==read(path.name)
        print("1490 independent coherent filters, unknown partner and finite precision checks passed.")
