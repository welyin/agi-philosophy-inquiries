"""Test one invariant triplet per sector against two measured left-handed weak edges."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def archived(n):
    return json.loads((HERE.parent / str(n) / "inputs.json").read_text(encoding="utf-8-sig"))


def yukawa(a,b,v):
    return np.array([[a*v[0],b*v[2],b*v[1]], [b*v[2],a*v[1],b*v[0]], [b*v[1],b*v[0],a*v[2]]])


def check_bases():
    w = np.exp(2j*math.pi/3)
    f = w**np.outer(np.arange(3),np.arange(3))/math.sqrt(3)
    bases = [np.eye(3)] + [np.diag([1,1,w**s])@f for s in range(3)]
    a,b = 1+.3j,-.4+.7j
    d = abs(a)**2+2*abs(b)**2
    weights = [abs(a)**2/d]+[abs(a+2*b*w**(-s))**2/(3*d) for s in range(3)]
    for i,u in enumerate(bases):
        for other in bases[:i]:
            np.testing.assert_allclose(abs(other.conj().T@u)**2, np.ones((3,3))/3)
        for v in ([1,0,0],[0,1,0],[0,0,1],[1+.2j,-3j, 2-.4j]):
            matrix = u.conj().T@yukawa(a,b,v)@u.conj()
            fraction = np.linalg.norm(np.diag(matrix))**2/np.linalg.norm(matrix)**2
            np.testing.assert_allclose(fraction,weights[i])
    np.testing.assert_allclose(sum(weights),2)
    # Two independent VEV components in a triplet can produce physical weak CP.
    up,down = yukawa(2,1,[1,2,0]),yukawa(2,1,[1,1j,0])
    h,k=up@up.conj().T,down@down.conj().T
    determinant = np.linalg.det(h@k-k@h)
    assert abs(determinant.imag)>1
    return dict(example_four_weights=weights,commutator_det_imag=round(determinant.imag))


def cert():
    u = np.array(archived(1527)["one_TeV_DR"]["up"])
    d = np.array(archived(1526)["down_yukawa_dr_1TeV"])
    weak = archived(1476)
    edges = {name:(value,errors) for name,value,errors in weak["edges"]}
    transfer = weak["edge_work_box"]["amplitude_transfer"]
    lower = {ab:(edges[ab][0]-2*sum(edges[ab][1]))/(1+transfer) for ab in ["cb","tb"]}
    intervals = math.asin(lower["cb"]), math.acos(lower["tb"])
    distances = [intervals[0],math.pi/2-intervals[1],math.acos(1/math.sqrt(3))-intervals[1]]
    assert min(distances) > 0
    rho,light = [],[]
    for data in (u,d):
        ratios = (data[:2,0]+2*data[:2,1])/(data[2,0]-2*data[2,1])
        rho.append(ratios[1])
        light.append(np.linalg.norm(ratios))

    def angle_budget(correction):
        errors = [ell+2*correction for ell in light]
        assert max((ell+correction)/(1-correction) for ell in light)<(1-1/math.sqrt(2))/math.sqrt(3)
        assert all(2*e+r<1 for e,r in zip(errors,rho))
        return sum(math.asin(e/(1-r-e)) for e,r in zip(errors,rho))

    lower_eps,upper_eps=0.,.03
    for _ in range(75):
        trial=(lower_eps+upper_eps)/2
        if angle_budget(trial)<min(distances):
            lower_eps=trial
        else:
            upper_eps=trial
    angle=angle_budget(0)
    return dict(group=check_bases(), data=dict(rho_ceiling=rho,light_over_heavy_ceiling=light,
                weak_amplitudes_lower=lower, physical_heavy_angle_range=intervals),
                max_unbroken_total_turn_rad=angle,
                same_or_orthogonal_wrong_edge_upper=math.sin(angle),
                unbiased_heavy_overlap_upper=math.cos(math.acos(1/math.sqrt(3))-angle),
                minimum_correction_frobenius_over_heavy=lower_eps,
                input_to_light_triplet_distance_floor=(1-1/math.sqrt(2))/math.sqrt(3))


def equivalent(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equivalent(v,b[k]) for k,v in a.items())
    return np.allclose(a,b,rtol=1e-11,atol=1e-13)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=cert()
    path=HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,allow_nan=False,indent=2)+"\n")
    else:
        assert equivalent(result,json.loads(path.read_text(encoding="utf-8")))
    print("1528 independent triplet and weak-current tests passed")
