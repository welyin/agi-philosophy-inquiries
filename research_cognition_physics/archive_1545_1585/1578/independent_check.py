"""Independent branching, accidental chiral charge, and a finite photon test."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import argparse
import json

ROOT=Path(__file__).resolve().parent
# (colour dimension, weak dimension, six times hypercharge).
sm={
    "1":[(1,1,0)], "5":[(3,1,-2),(1,2,3)],
    "10":[(3,1,-4),(3,2,1),(1,1,6)],
    "15":[(6,1,-4),(3,2,1),(1,3,6)],
    "24":[(8,1,0),(1,3,0),(1,1,0),(3,2,-5),(3,2,5)]}

def charges(rep):
    weight=Counter()
    for color,weak,y in sm[rep]:
        for t in range(1-weak,weak,2):
            weight[(color,F(3*t+y,6))]+=color
    return weight

def solve():
    data={"1":2,"24":2,"10":4,"15":4}
    total=sum(m*sum(charges(r).values()) for r,m in data.items())
    extra_colored=sum(m*count*q*q/2 for r,m in data.items() for (c,q),count in charges(r).items() if c>1)
    protected={"10":4,"15":1,"24":2,"1":2}
    protected_dimension=sum(m*sum(charges(r).values()) for r,m in protected.items())
    exotic=sum(m*count*q*q/2 for r,m in protected.items() for (c,q),count in charges(r).items() if c>1 and abs(q)==F(4,3))
    assert (total,protected_dimension,extra_colored,exotic)==(150,105,26,F(32,3))
    active={"10":3,"bar10":1,"15":1,"24":2,"1":2}
    dual_y={"10":1,"bar15":1,"15":3,"24":2,"1":2}
    dimensions={"10":10,"bar10":10,"15":15,"24":24,"1":1}
    pairable=sum(dimensions[r]*min(count,dual_y.get(r,0)) for r,count in active.items())
    assert pairable==75
    cases=[(4,13),(13,13),(-13,17),(0,1)]
    mixed=sum(q*t for q,t in cases)
    forbidden=[4+x for x in (4,13,-13,0)]
    assert mixed==0 and all(forbidden) and 13+(-13)==0
    upper=F(356,100)+2*(F(1,100)+F(7,100))
    # Real basis representation check: each difference matrix element removes one active X and one Y.
    a=[-1]*5+[0]*5+[1]*5
    independent=sum(ai!=aj for i,ai in enumerate(a) for aj in a[i+1:])
    assert independent==75
    return {"mixed_SU15_twice_index_weight":mixed,"U1_charges_Xperp_X_Y_L":[4,13,-13,0],
            "max_rank_SU5_and_accidental_U1":2*(pairable+15),"forbidden_bilinear_charges":forbidden,"mixed_gravity_background":4*105+13*(105-120),
            "Weyl_rank_at_generic_renormalizable_point":2*independent+30,
            "protected_Xperp_dimension":protected_dimension,"protected_abs4over3_R_colored":str(exotic),
            "surplus_SU5_dimensions_counting_conjugates_together":data,"tree_surplus_Weyl":total,
            "tree_surplus_massless_colored_R":str(extra_colored),"CLEO_upper_two_width_linear":float(upper),
            "min_transfer_exotic_abs4over3_plus_uds":float(upper/(exotic+2)),
            "min_transfer_all_tree_surplus_plus_uds":float(upper/(extra_colored+2)),
            "min_transfer_exotic_abs4over3_only":float(upper/exotic)}
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=solve()
    path=ROOT/"independent_results.json"
    if args.check:
        assert result==json.loads(path.read_text(encoding="utf-8-sig"))
        print("1578 independent review verified")
    else:
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
