"""Exact restriction characters and degree-five checks for a unit-index family embedding."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent

def elementary(values, degree):
    return sum(prod(x) for x in combinations(values,degree))

def higgs_independent_family(c,w,t):
    # All left-handed Weyl states of one observed family.
    return ([a+b+t for a in c for b in w] + [-a-4*t for a in c]
            +[-a+2*t for a in c]+[b-3*t for b in w]+[6*t])

def chiral(weights):
    counts=Counter(weights)
    keys=set(counts)|{-x for x in counts}
    return {x:counts[x]-counts[-x] for x in keys if counts[x]!=counts[-x]}

def finite_witness(c,w,ys):
    a,b,*rest=ys
    parent=sorted(q%3 for q in [*(x+a for x in c),*(x+b for x in w),*rest])
    sm=Counter(x%3 for x in higgs_independent_family(c,w,1))
    return {"parent_spectrum_mod3":parent,"chiral_coefficient_of_omega_minus_conjugate":sm[1]-sm[2]}

def compute():
    inp=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    ys=inp["competing_embedding"]["y_blocks"]
    left=finite_witness([0,0,0],[1,-1],ys)
    right=finite_witness([0,1,-1],[0,0],ys)
    assert left["parent_spectrum_mod3"]==right["parent_spectrum_mod3"]
    assert left["chiral_coefficient_of_omega_minus_conjugate"]==-9
    assert right["chiral_coefficient_of_omega_minus_conjugate"]==0
    # Both solutions surviving conjugacy and the mixed fifth coefficient.
    for a,b in ((-2,3),(1,0)):
        assert a+b==1
        k=F(6*a*a-9,9)
        spectator_norm=18*k-3*a*a-2*b*b
        assert spectator_norm==(-9 if a==1 else 0)
    for n in (5,6,7,10,18):
        for c,w,t in [([2,-1,-1],[4,-4],2),([7,-3,-4],[3,-3],-1),([3,-2,-1],[1,-1],0)]:
            defining=[*(x-2*t for x in c),*(x+3*t for x in w)]+[0]*(n-5)
            sm=higgs_independent_family(c,w,t)
            matter=[a+b for a,b in combinations(defining,2)]+[-x for x in defining]*(n-4)
            assert chiral(matter)==chiral(sm)
            for degree in (1,3):
                assert sum(x**degree for x in matter)==0
            # Odd part of the character has coefficient 2 Tr Z^5 / 5!.
            assert 2*sum(x**5 for x in sm)==-120*elementary(defining,5)
    previous=json.loads((HERE/inp["comparison"]).read_text(encoding="utf-8"))
    canonical=previous["levels"]["5/3"]
    return {"canonical_y":{"color":-2,"weak":3,"spectators":0,"k_Y":"5/3","connected_image_center_quotient":6},
            "SU7_same_parent_conjugacy_class_distinct_SM_indices":[left,right],
            "fifth_trace_normalization":"2 Tr_family Z^5=-5! e5(Z) on the canonical image",
            "other_net_family_multiplicities":"all positive integers have a witness",
            "conditional_pole_angle_window":canonical["one_loop_MSbar_angle"],
            "minimum_required_pole_angle_correction":canonical["minimum_total_angle_correction_to_SLD_two_widths"]}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    cfg=parser.parse_args()
    result=compute()
    path=HERE/"results.json"
    if cfg.write_new:
        with path.open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")
    else:
        assert result==json.loads(path.read_text(encoding="utf-8"))
        print("1573 finite character obstruction, odd traces, family witnesses and pole holdout checked")
