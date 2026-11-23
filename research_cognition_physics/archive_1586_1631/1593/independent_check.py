"""Finite spacelike electromagnetic response across both L3 large-angle windows."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads(name.read_text(encoding="utf-8-sig"))


def integration_measure(order):
    x,w=np.polynomial.legendre.leggauss(order)
    u=(x+1)/2
    return u*(1-u), w*x**2/(8*math.pi)


def inverse_alpha_step(mass,lower,upper,order=64):
    k,measure=integration_measure(order)
    return float(np.dot(measure, np.log1p(k*(upper-lower)/(mass**2+k*lower))))


def finite_pair(doublet,singlet,mix_sq,lower,upper,order=64,extra_diagonal=0.):
    # mix_sq is the off-diagonal mass-squared entry, not a probability.
    A=doublet**2+extra_diagonal+mix_sq**2/singlet**2
    B=singlet**2
    trace=A+B
    determinant=(doublet**2+extra_diagonal)*B
    large=(trace+math.hypot(A-B,2*mix_sq))/2
    small=determinant/large
    k,measure=integration_measure(order)
    x,y=k*lower,k*upper
    step=np.log1p((y-x)*(trace+x+y)/(determinant+trace*x+x*x))
    poles=sum(inverse_alpha_step(math.sqrt(m),lower,upper,order) for m in (small,large))
    direct=float(np.dot(measure,step))
    assert math.isclose(poles,direct,rel_tol=2.e-12,abs_tol=1.e-19)
    return direct


def compute():
    data=read(HERE.parent/"1592"/"results.json")
    weak=read(HERE.parent/"1369"/"inputs.json")["mu_ref_GeV"]
    central=data["central"]
    caps=[central["seven_doublet_geomean_cap_GeV"],
          data["two_width_working_box"]["seven_doublet_geomean_cap_GeV"][1]]
    s=central["seven_scale_max_GeV"]
    vev=246.22
    induced_diagonal=vev**2/2
    # S at the scale of meeting leaves the unbroken hypercharge-singlet
    # running length zero. A real trilinear of 0.1 S still mixes the two fields.
    mixing=.1*s*vev/math.sqrt(2)
    result={}
    for label,qlo,qhi in [("L3_2000_large_angle",12.25,3434.),
                          ("L3_2005_large_angle",1800.,21600.)]:
        assert weak**2>math.sqrt(qlo*qhi)/4
        floors=[7*inverse_alpha_step(cap,qlo,qhi) for cap in caps]
        if label=="L3_2000_large_angle":
            assert 7*math.log(qhi/qlo)/(12*math.pi)<1.29
        original=7*finite_pair(caps[0],s,mixing,qlo,qhi)
        changed=7*finite_pair(caps[0],s,mixing,qlo,qhi,extra_diagonal=induced_diagonal)
        assert original>=floors[0]*(1-1e-13) and changed<floors[0]
        # Test genuine EWSB pole splittings without a small-mixing assumption.
        vs=[finite_pair(190.,450.,b*190*450,qlo,qhi) for b in (0,.3,1,4,30)]
        assert all(p<q for p,q in zip(vs,vs[1:]))
        for cap,v in zip(caps,floors):
            assert math.isclose(v,7*inverse_alpha_step(cap,qlo,qhi,96),rel_tol=2e-13)
        result[label]={"spacelike_endpoints_GeV2":[qlo,qhi],"seven_pairs_floor_central_and_across_box":floors,
                       "seven_pairs_ceiling_arbitrary_positive_square_mixing":7*math.log(qhi/qlo)/(12*math.pi),
                       "ceiling_if_all_fourteen_poles_above_MZ":14*inverse_alpha_step(weak,qlo,qhi),
                       "seven_200_GeV_doublets_only":7*inverse_alpha_step(200.,qlo,qhi),
                       "same_unbroken_with_unity_diagonal_portal":changed,
                       "scaled_pair_splittings":vs}
    result["L3_2000_large_angle"]["reported"]=[3.80,.61,1.14,1.29]
    result["L3_2005_large_angle"]["reported_template_C"]=[1.05,.07,.14]
    return result


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(v,w) for v,w in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-18)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=compute()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as out:
            json.dump(data,out,indent=2)
            out.write("\n")
    else:
        assert same(data,read(path))
    print(json.dumps(data))
