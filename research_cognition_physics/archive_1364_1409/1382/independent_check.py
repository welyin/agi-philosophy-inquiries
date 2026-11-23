"""Check anomaly budgets, global charge lattices and the reflected photon amplitude."""
import argparse
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def spectrum(z):
    charged = [(n+1)/2 for n in z]+[(n-1)/2 for n in z]
    left_hypercharges = [u for n in z for u in (n/2,n/2,-(n+1)/2,-(n-1)/2)]
    assert sum(left_hypercharges)==0
    assert sum(u**3 for u in left_hypercharges)==-3*sum(z)/4
    return sum(u*u for u in charged)


def quotient_orders(z):
    # List all four subgroup orders under which every color-singlet Weyl field descends.
    fields=[(t,q) for n in z for t,q in ((1,3*n),(0,-3*(n+1)),(0,-3*(n-1)))]
    assert all(q.denominator==1 for _,q in fields)
    return [order for order in (1,2,3,6)
            if all((3*t+q)*(6//order)%6==0 for t,q in fields)]


def evaluate():
    baseline=json.loads((HERE.parent/"1373/results.json").read_text(encoding="utf-8"))
    amp=complex(*baseline["SM_gamma_amplitude"])
    mh=json.loads((HERE.parent/"1373/inputs.json").read_text(encoding="utf-8"))["higgs_mass_GeV"]
    # S<=23 bounds each |z|<=sqrt(46-1) and D<=4 when |z|>=3.
    candidate_z=[Q(z) for z in (-5,-3,3,5)]
    extremal={}
    for sign in (1,-1):
        candidates=[(spectrum(s),s) for n in (1,3)
                    for s in combinations_with_replacement(candidate_z,n) if sum(s)==-sign]
        s,z=min(candidates)
        assert s==23 and quotient_orders(z)==[1,2,3,6]
        ordinary=Q(5,3) if sign==1 else Q(1)
        ratio=abs(amp+4*(s+ordinary)/3)**2/abs(amp)**2
        extremal[str(sign)]={"minimal_exotic_S":int(s),"z":[int(n) for n in z],"ratio_floor":ratio}
    z=[Q(1,3)]*3+[Q(7,3),-Q(7,3)]
    assert sum(z)==1 and len(z)%2==1
    assert quotient_orders(z)==[1,2]
    S=spectrum(z)+1
    assert S==Q(82,9)
    # Positive triangle formula; invert in tau without using the analytic loop function.
    nodes,weights=np.polynomial.legendre.leggauss(72)
    x=(1+nodes[:,None])/2
    v=(1+nodes[None,:])/2
    wt=weights[:,None]*weights[None,:]/4
    k=4*x*(1-x)*v
    def loop(tau):
        return 4*np.sum(wt*(1-x)*(1-k)/(1-tau*k))
    goal=-2*amp.real/float(S)
    assert 4/3<goal<2
    low,high=0.,1.
    for _ in range(72):
        t=(low+high)/2
        if loop(t)<goal:
            low=t
        else:
            high=t
    t=(low+high)/2
    reflected=abs(amp+float(S)*loop(t))**2/abs(amp)**2
    np.testing.assert_allclose(reflected,1,atol=1e-12)
    return {"integer_cases":extremal, "fractional_Z2":{"z":[str(n) for n in z],
            "quotient_orders":quotient_orders(z),"total_charge_squares":str(S),
            "reflected_target_fermion_loop":goal,"one_common_mass_GeV":mh/(2*sqrt(t)),
            "photon_width_ratio":reflected}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    value=evaluate()
    path=HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(value,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        saved=json.loads(path.read_text(encoding="utf-8"))
        for sign,case in value["integer_cases"].items():
            assert case["z"]==saved["integer_cases"][sign]["z"]
            assert case["minimal_exotic_S"]==saved["integer_cases"][sign]["minimal_exotic_S"]
            np.testing.assert_allclose(case["ratio_floor"],saved["integer_cases"][sign]["ratio_floor"])
        for k,v in value["fractional_Z2"].items():
            if isinstance(v,float):
                np.testing.assert_allclose(v,saved["fractional_Z2"][k],rtol=1e-11)
            else:
                assert v==saved["fractional_Z2"][k]
        print("1382: anomaly minimum, representation descent and equal-width example verified.")
