"""Reproduce the baryon selection, projected Ward identities and twisted Weyl index."""
import argparse
from fractions import Fraction
from itertools import combinations, product
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def cyclic_check(max_order):
    admitted = protected = 0
    for n in range(2,max_order+1):
        for q,l,h in product(range(n),repeat=3):
            d = (3*q+l)%n
            if (2*(l+h))%n or 3*d%n:
                continue
            admitted += 1
            if not d:
                continue
            protected += 1
            assert n//math.gcd(n,d)==3
            u, e, down = -q-h,-l+h,-q+h
            assert (u+e-2*q+d)%n==0 and (q-e-l+down)%n==0
            for db,dl in product(range(-3,4),repeat=2):
                if (dl-db)%2==0:
                    assert ((3*q-h)*db+(l+h)*dl-d*db)%n==0
    return {"solutions_scanned":admitted,"protecting_solutions":protected}


def index_cost(lam, families):
    # Each phase is stored in turns, without any restriction to cube roots.
    # Uniform phases common on a single Weyl orbit can be dropped.
    third=Fraction(1,3)
    sector={"u":(third-lam)%1,"Q":Fraction(),"e":(third+lam)%1,
            "d":(third+lam)%1,"L":Fraction()}
    answer=[]
    for orbit in ({"u":3,"Q":6,"e":1},{"d":3,"L":2}):
        distance=0
        for angle in {sector[k] for k in orbit}:
            values={k:families*int(sector[k]==angle) for k in orbit}
            costs=[sum(dim*abs(values[k]-c) for k,dim in orbit.items())
                   for c in range(families+1)]
            distance+=min(costs)
        answer.append(distance)
    return answer


def mass_example(data):
    charges=data["X_mod3"]
    f={key:(x+data["six_Y"][key])%3 for key,x in charges.items()}
    assert [f[key] for key in ["Q","u","e","d","L","H"]]==[1,0,1,1,1,2]
    adj=data["adjoint_eigenvalues"]
    # Adjoint acts on the antisymmetric two-index representation.
    m=[adj[i]+adj[j]-2*adj[0] for i,j in combinations(range(5),2)]
    assert sorted(m)==[-10]+[-5]*6+[0]*3
    return {"commuting_phases_mod3":f,"10_pair_masses_over_yV":sorted(m),
            "extra_light_Weyl_components":2*m.count(0)*data["families"]}


def ward_example():
    z=np.exp(2j*np.pi/3)
    neutral=np.diag([z,z.conjugate()])
    d=np.array([1,0],complex)
    l=np.array([0,1],complex)
    f=np.array([[.25,.02j],[-.02j,.5]])  # positive propagator, units mass^(-2)
    exact=np.diag(np.diag(f))
    comm=np.linalg.norm(neutral@f-f@neutral,2)
    assert np.all(np.linalg.eigvalsh(f)>0) and d.conj()@exact@l==0
    actual=abs(d.conj()@f@l)
    assert abs(comm/math.sqrt(3)-actual)<1e-13
    # A single pole can saturate the minimax defect if its phase is free.
    assert np.allclose([abs(-1-z),abs(-1-z.conjugate())],[1,1])
    return {"two_pole_exact_cross":0,"two_pole_broken_cross":actual,
            "commutator_norm":comm,"saturating_bound":comm/math.sqrt(3),
            "single_phase_minimax_relative_defect":1}


def calculate(data):
    n=data["families"]
    # Generic values cannot produce a new coefficient coincidence with Q or d.
    # Include a nontrivial possible SU(5)-centre lift as an explicit example.
    cases={"one":Fraction(0),"omega":Fraction(1,3),"omega_squared":Fraction(2,3),
           "generic":Fraction(1,17),"centre_extension_example":Fraction(1,15)}
    costs={k:index_cost(v,n) for k,v in cases.items()}
    assert costs["omega"]==[6,12] and costs["omega_squared"]==[18,0]
    assert all(sum(c)==(18 if k in ("omega","omega_squared") else 36)
               for k,c in costs.items())
    # One unknown index column on each of the ten wedge-square and five dual weights.
    # The positive and negative SU(5)-centre classes are pairwise distinct.
    assert len({2%5,-2%5,1%5,-1%5})==4
    exp=data["experiments"]
    mean=exp["Ge76_reported_half_life_yr_90"]/math.log(2)
    assert math.isclose(mean,exp["Ge76_ppp_and_ppn_each_partial_mean_lifetime_yr_90"],rel_tol=.01)
    return {"charge_enumeration":cyclic_check(data["discrete_enumeration_max_modulus"]),
            "weyl_orbits_10_bar5_costs":costs,"least_added_components":6*n,
            "adjoint_example":mass_example(data),"ward_example":ward_example(),
            "Ge76_converted_mean_lifetime_yr":mean}


def close(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(close(a[k],b[k]) for k in a)
    return np.allclose(a,b,rtol=2e-10,atol=1e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    switch=parser.add_mutually_exclusive_group(required=True)
    switch.add_argument("--check",action="store_true")
    switch.add_argument("--write",action="store_true")
    args=parser.parse_args()
    out=calculate(read("inputs.json"))
    if args.write:
        with (HERE/"results.json").open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2,allow_nan=False)
            f.write("\n")
    elif not close(read("results.json"),out):
        raise SystemExit("Results differ")
    else:
        print("1541: selection, Ward identity, continuous twists and sharp chiral inventory verified")
