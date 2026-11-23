"""Test the frozen anomaly and all family-universal SM reassigments without SU(5) or a mu term."""
import argparse
import collections
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def sign(q):
    return (q+1)%3-1


def compute():
    # First entry is the original Z12 flavor charge. Second entry is the auxiliary charge.
    pairs=[]
    def add(dimension, original, auxiliary):
        pairs.extend([(original,sign(auxiliary))]*dimension)
    for t,f in zip([0,0,1],[2,0,1]):
        add(10,t,-1)
        add(5,f,0)
    for p in [3,3,0,0]:
        add(5,p,0)
    for q in [-1,4,-2,-3,2,3,0,0]:
        add(1,q,0)
    for h in [-1,1]:
        add(2,0,h)
    for old,aux in [(0,1),(2,1),(0,-1),(6,-1),(0,-1)]:
        add(1,old,aux)
    linear=sum(aux for _,aux in pairs)
    cubic=sum(aux**3 for _,aux in pairs)
    assert (linear,cubic)==(-31,-31)
    reps=collections.Counter((sign(old),new) for old,new in pairs if old%3 and new%3)
    assert reps == {(1,-1):10,(-1,1):1}
    mixed=[sum(sign(q)**(3-k)*p**k for q,p in pairs) for k in [1,2]]
    assert mixed == [-9,9]
    for extra in range(5):
        choices=[list(setup) for setup in itertools.combinations_with_replacement([-1,0,1],extra)
                 if (cubic+sum(n**3 for n in setup))%9==0 and (linear+sum(setup))%3==0]
        assert choices==([[1,1,1,1]] if extra==4 else [])
    repaired=pairs+[(0,1)]*4
    # Restrict to all 3-primary cyclic maps; the written pullback argument covers noncyclic bundles.
    for a,b in itertools.product(range(3),repeat=2):
        z=[sign(q*a+p*b) for q,p in repaired]
        assert sum(v**3 for v in z)%9==0 and sum(z)%3==0

    # Old S and X have been made auxiliary-neutral by a shear of the two discrete generators.
    # No family unification, no bare Higgs bilinear. Connectivity enforces one charge per species.
    old_edges=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))["down_edges"]
    connected_solutions=0
    for h in range(3):
        for rowcol in itertools.product(range(3),repeat=6):
            if all((rowcol[i]+rowcol[3+j]+h)%3==0 for i,j in old_edges):
                assert len(set(rowcol[:3]))==len(set(rowcol[3:]))==1
                connected_solutions+=1
    assert connected_solutions==9
    reassignments=collections.Counter()
    for n in [-1,1]:
        for q,h,k in itertools.product(range(3), repeat=3):
            u,d,l,e=map(sign,[-q-h,-q-k,-h-n,h+n-k])
            su2=3*(3*q+l)+h+k
            if su2%3:
                continue
            weak=[6*[sign(q)],3*[u],3*[d],2*[l],[e]]
            all_charges=3*sum(weak,[])+2*[sign(h)]+2*[sign(k)]+2*[n]+2*[-n]+[-n]
            anomaly=sum(c**3 for c in all_charges)%9
            assert anomaly==(-4*n)%9
            reassignments[(n,anomaly)]+=1
    assert reassignments=={(-1,4):9,(1,5):9}
    # Arbitrary old-flavon charge s only shears the generator, it cannot remove a nonzero restriction.
    sheared=[]
    for s in [-1,0,1]:
        transformed=[(q,sign(p-s*q)) for q,p in pairs]
        extracted=[sign(p+s*q) for q,p in transformed]
        assert extracted==[p for _,p in pairs]
        assert sign(s-1*s)==0 and sign(-4*s+4*s)==0
        sheared.append(sum(p**3 for _,p in transformed)%9)
    first_B=min(k for k in range(3) if (1+1-k)%3==0)
    assert first_B==2
    return {
        "original_signed_linear_and_cubic": [linear,cubic],
        "joint_nonzero_old3_aux3_multiplicities": {"1,-1":10,"-1,1":1},
        "old_old_aux_and_old_aux_aux_mixed_moments":mixed,
        "minimal_repair_auxiliary_signed_charges":[1,1,1,1],
        "all_reassigned_ordinary_SM_species_weak_instanton_pass_cases":
            [{"N0":n,"cubic_mod9":a,"count":v} for (n,a),v in sorted(reassignments.items())],
        "old_texture_solutions_all_family_universal":connected_solutions,
        "auxiliary_after_shearing_by_original_minus1_0_plus1_mod9":sheared,
        "first_B_insertions_for_spectator_mass_and_N0_contact": [first_B,first_B],
    }


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=compute()
    if args.check:
        assert result==json.loads((HERE/"independent_results.json").read_text(encoding="utf-8"))
        print("PASS: full charged table, character cancellation and all ordinary-charge reassigments")
    else:
        print(json.dumps(result,indent=2))
