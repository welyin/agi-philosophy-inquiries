"""Exact fifth-order character test and unrestricted-rank sufficiency witnesses."""
import argparse
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent


def elementary(values, degree):
    if degree == 0:
        return 1
    return sum(math.prod(subset) for subset in itertools.combinations(values,degree))


def annihilating_polynomial(a):
    b = 1-a
    s1 = -2-a
    s2 = 9-b*b-2*b*s1
    lhs = (b-a)*(a*(b+a)+(b+3*a)*s1+2*s2)+10
    assert lhs == 16*(a-1)**2*(a+2)
    assert s1*s1-2*s2 == (a+2)*(7*a-10)
    return lhs


def parent_moments(fundamental):
    """Tr H^j − Tr(−H)^j on wedge^2 F + (N−4) Fbar."""
    assert sum(fundamental)==0
    charges = [x+y for x,y in itertools.combinations(fundamental,2)]
    charges += [-q for q in fundamental]*(len(fundamental)-4)
    return [sum(x**j-(-x)**j for x in charges) for j in range(1,6)]


def expansion_check(a,b,spectators,color,weak,t):
    """Compute e5 directly from diagonal entries, then from independent factor expansion."""
    assert sum(color) == 0 and sum(weak) == 0 and 3*a+2*b+sum(spectators)==0
    p,q,r = elementary(color,2), elementary(color,3), elementary(weak,2)
    y=[a]*3+[b]*2+spectators
    S=[elementary(spectators,j) for j in range(4)]
    coeff_prt = a+S[1]
    coeff_qt2 = b*b + 2*b*S[1] + S[2]
    coeff_pt3 = a*b*b+(b*b+2*a*b)*S[1]+(a+2*b)*S[2]+S[3]
    coeff_rt3 = a**3+3*a*a*S[1]+3*a*S[2]+S[3]
    actual = [x+t*a for x in color]+[x+t*b for x in weak]+[t*x for x in spectators]
    reduced=q*r + p*r*t*coeff_prt + q*t*t*coeff_qt2
    reduced+=(p*coeff_pt3+r*coeff_rt3)*t**3+elementary(y,5)*t**5
    assert reduced == elementary(actual,5)
    return actual


def run():
    # A polynomial of degree at most three is fixed by four evaluations.
    node_values=[annihilating_polynomial(F(i)) for i in (-3,-1,0,2)]
    assert node_values == [-256,64,32,64]
    r = []
    for a in (-2,1):
        b=1-a
        s1=-2-a
        s2=12-3*a*a
        r.append({"a":a,"b":b,"spectator_sum":s1,"spectator_sum_of_squares":s1*s1-2*s2})
    color, weak = [1,3,-4], [2,-2]
    tests=0
    for a,b,left,t in itertools.product((-2,0,F(2,3)),(1,3),(-3,-1,0),(0,1,F(-2,3))):
        s=[left,-3*a-2*b-left]
        expansion_check(a,b,s,color,weak,t)
        tests+=1
    canonical=[v-2 for v in color]+[v+3 for v in weak]
    sm=parent_moments(canonical)
    validity=[]
    for n in (5,6,7,14,37):
        row=parent_moments(canonical+[0]*(n-5))
        assert row == sm == [0,0,0,0,-120*math.prod(canonical)]
        validity.append(n)
    exotic=expansion_check(-2,1,[1,3],color,weak,1)
    exotic_moments=parent_moments(exotic)
    assert exotic_moments[:4] == sm[:4] and exotic_moments[4] != sm[4]
    return {"coefficient_identity_nodes":[-3,-1,0,2],
            "cast_integer_node_values":list(map(int,node_values)),
            "possible_roots_and_positivity":r,"expanded_cartan_tests":tests,
            "finite_positive_parent_dimension_checks":validity,
            "canonical_chiral_character_moments":sm,
            "noncanonical_SU7_chiral_character_moments":exotic_moments}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--check",action="store_true")
    flags.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=run()
    out=HERE/"independent_results.json"
    if args.write_new:
        with out.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    elif result!=json.loads(out.read_text(encoding="utf-8-sig")):
        raise SystemExit("1573 independent fifth-order identity changed")
    print("1573 exact fifth-order necessity and anomaly-free realizability verified")
