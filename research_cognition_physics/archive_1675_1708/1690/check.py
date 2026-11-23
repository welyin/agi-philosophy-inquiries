"""Exact indices for the specified SU(5)-extending fermions, including center flux."""
import argparse
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def rep_weights(roots, rep):
    if rep=="10":
        return [i+j for i,j in combinations(roots,2)]
    if rep in ("5","bar5"):
        sign=1 if rep=="5" else -1
        return [sign*x for x in roots]
    return [0]


def indices(a,b,rep):
    return sum(x*y for x,y in zip(rep_weights(a,rep),rep_weights(b,rep)))


def split_fluxes():
    for prefix in product((-1,0,1),repeat=4):
        yield (*prefix,-sum(prefix))


def calculate():
    cfg=read("inputs.json")
    spec=read(cfg["spectrum"])
    previous=read(cfg["previous"])
    med=read(cfg["mediator"])
    fields=[("10",q) for q in spec["su5_generation_charges_10"]]
    fields.extend(("bar5",q) for q in spec["su5_generation_charges_5bar"])
    heavy=med["flavor_charge"]
    fields.extend((rep,heavy[name]) for rep,name in
                  [("5","P"),("bar5","Pbar"),("5","T"),("bar5","Tbar")])
    fields.extend(("1",heavy[name]) for name in spec["neutral_charged_under_flavor"])
    a=cfg["test_cartan"]
    quadratic={}
    cubic={}
    for rep in ("1","5","bar5","10"):
        weights=rep_weights(a,rep)
        quadratic[rep]=Fraction(sum(x*x for x in weights),sum(x*x for x in a))
        cubic[rep]=Fraction(sum(x**3 for x in weights),sum(x**3 for x in a))
    assert quadratic["10"]==3 and cubic["10"]==1 and cubic["bar5"]==-1
    assert sum(cubic[rep] for rep,_ in fields)==0
    coefficient=sum(q*quadratic[rep] for rep,q in fields)
    assert coefficient==12
    a1=sum(q*len(rep_weights(a,rep)) for rep,q in fields)
    a3=sum(q**3*len(rep_weights(a,rep)) for rep,q in fields)
    assert (a1,a3)==(58,388)
    old_pairs=previous["singlet_solutions_by_multiplicity"]["2"]
    assert all((a1+sum(pair))%6==0 and (a3+sum(q**3 for q in pair))%36==0
               for pair in old_pairs)
    for p,h in product(range(12),repeat=2):
        assert (coefficient*p*h)%12==0
    # Split bundles on S^2 x S^2, including bundles unavailable on the product cover.
    total=0
    admissible={str(n):0 for n in cfg["global_forms"]}
    for a,b in product(split_fluxes(),repeat=2):
        fundamental=indices(a,b,"5")
        actual=sum(q*indices(a,b,rep) for rep,q in fields)
        assert actual==coefficient*fundamental
        assert actual%12==0
        for n in cfg["global_forms"]:
            multiple=6//n
            if sum(a[3:])%multiple==0 and sum(b[3:])%multiple==0:
                admissible[str(n)]+=1
        total+=1
    flux=[1,0,0,-1,0]
    weak_c1=(sum(flux[3:]),sum(flux[3:]))
    k3=-Fraction(2,3)
    k2=-Fraction(1,2)
    kY=-Fraction(5,6)
    k_total=k3+k2+kY
    index=sum(q*indices(flux,flux,rep) for rep,q in fields)
    assert k_total==-indices(flux,flux,"5")==-2 and index==24
    example={"fundamental_c2":int(k_total),"signed_chiral_index_weighted_by_flavor":index,
             "projective_SU3_and_SU2_and_abelian_c2":[str(s) for s in (k3,k2,kY)],
             "weak_determinant_chern_flux":list(weak_c1)}
    competitor={}
    for j in cfg["global_forms"]:
        k=6//j
        charge_index=2*k**2
        competitor[str(j)]={"mixed_exponent_mod12_for_t1_m1":charge_index%12,
                             "allowed":charge_index%12==0}
        for t,m in product(range(12),repeat=2):
            assert ((t*m*m*charge_index)%12==0)==((t*m*m*k*k)%6==0)
    return {"local_SU5_cubic_sum":0,"mixed_index_coefficient":int(coefficient),
            "mixed_bordism_generator_order":12,"passing_neutral_pairs":old_pairs,
            "before_singlets_A1_A3":[a1,a3],"flux_backgrounds_checked":total,
            "counts_liftable_to_each_SM_quotient_on_test_surface":admissible,
            "explicit_nonliftable_example":example,"added_hypercharge_pair_test":competitor}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write",action="store_true")
    opts=cli.parse_args()
    answer=calculate()
    if opts.write:
        (HERE/"results.json").write_text(json.dumps(answer,indent=2)+"\n",encoding="utf-8")
    else:
        assert answer==read("results.json")
        print("PASS: SU(5) and mixed indices, six neutral pairs and 6561 split quotient bundles")
