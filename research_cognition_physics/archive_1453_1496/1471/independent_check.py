"""Test the complete finite lattice and the exterior-energy test for a second photon."""
import argparse
from fractions import Fraction
from itertools import product
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def solve():
    p = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    forms = p["global_forms"]
    cases = 0
    for k1,k2,t1,t2,n in product(range(-4,7),range(-2,3),range(3),range(3),range(7)):
        allowed = []
        for order in forms:
            subgroup = {(6*i//order)%6 for i in range(order)}
            electric = all((h*(k+2*t))%6==0 for h in subgroup for k,t in ((k1,t1),(k2,t2)))
            if electric and n%6 in subgroup:
                allowed.append(order)
        required = 6//math.gcd(6,n)
        compatible = math.gcd(6,k1+2*t1,k2+2*t2)
        assert allowed == [x for x in forms if x%required==0 and compatible%x==0]
        cases += 1
    # Exact rational data for the real, canonically normalized two-photon example.
    n = p["positive_test"]["magnetic_n"]
    q = Fraction(p["positive_test"]["electric_k"],6)
    inv_alpha = Fraction(137036,1000)
    relative = Fraction(1,25)
    budget = 350
    def exterior_times_R(m):
        return inv_alpha*(n*n+(m-n*q)**2/relative)/8
    field = {str(m):float(exterior_times_R(m)) for m in range(5)}
    # Integer m outside this interval is farther from nq and hence costs more.
    lo = math.ceil(float(n*q)-math.sqrt(float((8*budget/inv_alpha-n*n)*relative)))
    hi = math.floor(float(n*q)+math.sqrt(float((8*budget/inv_alpha-n*n)*relative)))
    assert (lo,hi)==(2,2) and exterior_times_R(lo)<=budget
    assert all(exterior_times_R(m)>budget for m in [lo-1,hi+1])
    hidden_charge_floor = (n*q)**2/(8*budget/inv_alpha-n*n)
    return {"lattice_cases":cases,"alpha_inverse_illustration":float(inv_alpha),
            "hidden_electric_gram_cap":float(relative),"exterior_ER_budget":budget,
            "magnetic_energy_floor_times_R":field,"allowed_hidden_winding":[lo,hi],
            "gram_lower_if_hidden_winding_zero":float(hidden_charge_floor)}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=solve()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as out:
            json.dump(data,out,indent=2)
            out.write("\n")
    else:
        assert data==json.loads(path.read_text(encoding="utf-8"))
        print("1471: center subgroups and the two-photon energy obstruction verified")
