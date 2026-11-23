"""Reproduce the anomaly lattice with vectorlike matter and test finite resolution."""
import argparse
import json
from fractions import Fraction as F
from itertools import product
from pathlib import Path

BASE = Path(__file__).resolve().parent
GENERATORS = {1:(0,0),2:(0,1),3:(1,0),6:(-2,3)}
TARGET = F(5,3)
WEAK = F(3,2)


def index(c, w, y):
    """One fundamental or singlet in each factor, with axial PQ difference +1."""
    return (int(c==3)*w, int(w==2)*c, c*w*y*y)


def matter_lattice(q):
    y3, y2 = GENERATORS[q]
    slots = [(3,1,y3),(1,2,y2),(1,1,q)]
    assert all((2*int(c==3)+3*int(w==2)+y)%q==0 for c,w,y in slots)
    return tuple(index(*slot) for slot in slots)


def coefficient(q, n3, n2, n1):
    return tuple(sum(n*column[i] for n,column in zip((n3,n2,n1),matter_lattice(q)))
                 for i in range(3))


def counts(q, triple):
    k,h,b = triple
    v = matter_lattice(q)
    rest, remainder = divmod(b-k*v[0][2]-h*v[1][2],q*q)
    return None if remainder else (k,h,rest)


def photon(t):
    k,h,b=t
    return F(b+18*h,18*k)


def joint_nearest(q,k):
    h=k*WEAK.numerator//WEAK.denominator
    base=coefficient(q,k,h,0)
    real_n=(18*k*TARGET-base[2]-18*h)/(q*q)
    n0=real_n.numerator//real_n.denominator
    candidates=[(abs(photon(triple)-TARGET),triple) for n in (n0,n0+1)
                for triple in [coefficient(q,k,h,n)]]
    delta, triple=min(candidates)
    expected={1:0,2:F(1,9*k),3:0,6:F(1,k)}[q] if k%2 else 0
    assert delta == expected and abs(F(h,k)-WEAK)==F(k%2,2*k)
    assert coefficient(q,*counts(q,triple))==triple
    return {"K":triple, "signed_messenger_counts": counts(q,triple),
            "photon":str(photon(triple)), "weak":str(F(h,k)),
            "abs_photon_error":str(delta), "abs_weak_error":str(abs(F(h,k)-WEAK))}


def audit():
    for q in GENERATORS:
        for k,h,b in product(range(-5,6), range(-5,6), range(-43,44)):
            ns=counts(q,(k,h,b))
            # The published spin-bundle condition, independent of the three
            # triangular basis columns used for the explicit completion.
            topological=(b%q==0 and (24*k+18*h+b)%(q*q)==0)
            assert topological == (ns is not None)
            if ns is not None:
                assert coefficient(q,*ns)==(k,h,b)
        for k in range(1,60):
            joint_nearest(q,k)
    real_quark=index(3,2,1)
    assert real_quark==(2,3,6) and photon(real_quark)==TARGET
    two_pairs=tuple(x+y for x,y in zip(index(3,1,-2), index(1,2,0)))
    assert two_pairs==(1,1,12) and photon(two_pairs)==TARGET
    for q in GENERATORS:
        assert (counts(q,two_pairs) is not None)==(q in (1,3))
    three_lepton_shift=(0,3,-54)
    after=tuple(-x+y for x,y in zip(real_quark,three_lepton_shift))
    assert after==(-2,0,-60) and photon(after)==TARGET
    assert all(counts(q,after) is not None for q in GENERATORS)
    return {"gauge_basis_in_K3_K2_K1_order":{str(q):matter_lattice(q) for q in GENERATORS},
            "two_messenger_species":two_pairs,"lepton_current_shift_if_allowed":after,
            "both_ratios_nearest":{str(k):{str(q):joint_nearest(q,k) for q in GENERATORS}
                                   for k in (1,2,5,101)}}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
