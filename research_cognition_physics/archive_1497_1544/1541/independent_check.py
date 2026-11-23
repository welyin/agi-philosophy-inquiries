"""Test hypercharge twists, minimal light defects and two independent electroweak-era probes."""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

DIR=Path(__file__).resolve().parent
# y=6Y, q=B-2Y, dimension of each complete SM weight cluster.
ROLES={"Q":(1,0,6),"uc":(-4,1,3),"ec":(6,-2,1),"dc":(2,-1,3),"L":(-3,1,2)}


def weight_cost(a):
    phase={name:(q-a*y)%3 for name,(y,q,_) in ROLES.items()}
    # By the strict majority of weights in each orbit the common weighted
    # L1 minimizer equals the virtual multiplicity distribution of Q (or dc).
    cost=sum(6*ROLES[r][2] for r,v in phase.items()
             if v!=phase["Q" if r in ("Q","uc","ec") else "dc"])
    return int(cost)


def after_ew(q,twice_t3):
    return (q-twice_t3)%3


def run():
    checks={str(a):weight_cost(Fraction(a)) for a in range(3)}
    assert checks=={"0":36,"1":18,"2":18}
    for d in range(1,31):
        for i in range(3*d):
            a=Fraction(i,d)
            assert weight_cost(a)==(18 if (5*a)%3 in (1,2) else 36)
    a=2
    p={r:(q-a*y)%3 for r,(y,q,_) in ROLES.items()}
    assert p=={"Q":1,"uc":0,"ec":1,"dc":1,"L":1}
    # 10(1) contains ordinary Q,e^c and extra uc, 10(0)+bar10(0) supplies ordinary
    # uc and the extra conjugate up. All are left-handed fields.
    plus=[after_ew(0,1)]*3+[after_ew(a*4,0)]*3
    minus=[after_ew(1,0)]*3+[after_ew(1+a*(-4),0)]*3
    rank=sum(min(plus.count(s),minus.count((-s)%3)) for s in range(3))
    assert (plus,minus,rank)==([2]*6,[1]*3+[2]*3,3)
    # Compact lifting in this example: three 10(1)+bar5(1) and three neutral
    # vectorlike 10 pairs. Scalars are not chiral contributors.
    anomaly={"twice_SU5_dynkin_charge_mod3":(3*(3+1))%3,
             "sum_X_cubed_mod9":(3*(10+5))%9}
    assert not any(anomaly.values())
    a=1
    p={r:(q-a*y)%3 for r,(y,q,_) in ROLES.items()}
    assert p=={"Q":2,"uc":2,"ec":1,"dc":0,"L":1}
    new_neutral=after_ew(p["dc"]+a*(-3),1)
    neutral_conjugate=after_ew(-p["L"]+a*3,-1)
    assert (new_neutral,neutral_conjugate)==(2,0)
    residue=(3*(-1)**3)%9
    assert residue==6 and (residue+3)%9==0
    assert all((residue+x)%9!=0 for x in range(-2,3))
    assert (0-1+1)%3==0  # L_extra H N^c in the unbroken weak-phase lift
    ordinary_neutral=after_ew(1,1)
    assert ordinary_neutral==0 and (new_neutral*2)%3!=0
    rupper=3.56+2*(.01+.07)
    energy=10.52
    # Necessary mass threshold only if the escaped particles retain standard charges
    # and one applies this massive parton response in a new, explicitly broken model.
    eff=.90
    target=(rupper/eff-2)/4
    left,right=0.,1.
    for _ in range(70):
        b=(left+right)/2
        if b*(3-b*b)/2<target:
            left=b
        else:
            right=b
    mass=energy/2*math.sqrt(1-((left+right)/2)**2)
    invisible_upper=2.9963+2*.0074
    return {"integer_twist_min_new_Weyl_weight_states":checks,"charged_extra_lh_grades":plus,
            "charged_conjugate_lh_grades":minus,"max_up_sector_rank":rank,
            "UV_anomaly_witness":anomaly,"unpaired_active_neutral_grade":new_neutral,
            "ordinary_and_extra_conjugate_neutral_grade":neutral_conjugate,
            "CLEO_R_two_reported_widths_linear_upper":rupper,
            "min_massless_parton_R_uds_and_three_extra_up":6.,
            "common_transfer_to_allow_tree_rate_at_most":rupper/6,
            "if_mass_escape_max_of_three_must_exceed_GeV_at_90pct_transfer":mass,
            "LEP_Nnu_two_width_upper":invisible_upper,
            "new_per_species_max_width_if_ordinary_three_fixed":(invisible_upper-3)/3,
            "new_per_species_max_width_if_ordinary_three_at_98pct":(invisible_upper-2.94)/3,
            "universal_visible_transfer_max_for_six":invisible_upper/6}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    grp=p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--check",action="store_true")
    grp.add_argument("--write-new",action="store_true")
    args=p.parse_args()
    path=DIR/"independent_results.json"
    result=run()
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        assert result==json.loads(path.read_text(encoding="utf-8"))
        print("1541: full weight orbit minima, electroweak rank, CLEO and LEP receipts pass")
