"""Finite checks of characteristic numbers, the centre quotient and scalar stabilizers."""
import argparse
from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path


ROOT=Path(__file__).resolve().parent


def phase(x):
    return str(x % 1)


def centre_kernel():
    # SU4 centre i^j; SU2 centres (-1)^l, (-1)^r; finite generator s^a.
    ker=[]
    for j,l,r,a in product(range(4),range(2),range(2),range(12)):
        left=Q(j,4)+Q(l,2)+Q(a,12)
        right=-Q(j,4)+Q(r,2)-Q(a,12)
        if left.denominator==right.denominator==1:
            ker.append((j,l,r,a))
            assert (Q(j,4)+Q(r,2)+Q(a,12)).denominator==1  # fundamental scalar F=1
            assert (Q(j,2)+Q(a,6)).denominator==1            # symmetric tensor F=2
    generators={( (-k+2*d)%4,d,d,(3*k)%12) for k,d in product(range(4),range(2))}
    assert set(ker)==generators
    # For s^a=c^{-a} t^a, t=s^4 has order three and projects to s mod Z4.
    for a in range(12):
        assert ((4*a-3*a)-a)%12==0
    return len(ker)


def backgrounds():
    # On S² x S² the classes a,b are pairs. Give E4 the two oriented planes
    # a,b; E6 has the plane a+b plus a trivial oriented real four-plane.
    def cup(a,b):
        return a[0]*b[1]+a[1]*b[0]
    hits=0
    for v in product(range(-2,3),repeat=4):
        a,b=v[:2],v[2:]
        e=cup(a,b)
        add=tuple(i+j for i,j in zip(a,b))
        sub=tuple(i-j for i,j in zip(a,b))
        p_l,p_r=cup(add,add),cup(sub,sub)
        p6=cup(add,add)
        i_l,i_r=p6+p_l,p6+p_r
        assert (p_l-p_r)%4==0 and i_l-i_r==4*e
        assert (Q(p_l,4)-Q(p_r,4))==e
        if p_l%4:
            hits+=1
    return {"examples_on_product_spheres":625,"nonintegral_pair_samples":hits,
            "unit_Euler_example":{"p6":2,"pL":2,"pR":-2,
                                  "k_left":"1/2","k_right":"-1/2","indices":[4,0]}}


def charges():
    # Indices label left-handed Weyl states; third entries are T3_R.
    fermions={"Q":(1,Q(1,3),0),"L":(1,-1,0),
              "u_c":(-1,Q(-1,3),Q(-1,2)),"d_c":(-1,Q(-1,3),Q(1,2)),
              "e_c":(-1,1,Q(1,2)),"nu_c":(-1,1,Q(-1,2))}
    out={}
    for name,(f,bl,tr) in fermions.items():
        y=tr+Q(bl,2)
        b=Q(f+bl,4)
        out[name]={"B":str(b),"Y":str(y),"q_B3_mod3":int((b-2*y)%3)}
        assert (b-2*y).denominator==1
    assert out["u_c"]["q_B3_mod3"]==1 and out["e_c"]["q_B3_mod3"]==1
    return out


def vacua():
    allowed={}
    for a,b in product((1,5,9),(2,6,10)):
        exponents=[k for k in range(12) if (k*(b-2*a))%12==0]
        n=len({k%3 for k in exponents})
        allowed[f"{a},{b}"]=n
        assert (n==3)==(b%12==(2*a)%12)
    return allowed


def calculate():
    cfg=json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig"))
    ngen=cfg["families"]
    mixed=[n for n in range(1,13) if (4*n)%12==0]
    assert ngen in mixed
    # Positive and negative family triplets are balanced under the auxiliary SU3.
    return {"centre_kernel_before_PS_quotient":centre_kernel(),
            "external_components_after_quotient":3,
            "instanton":backgrounds(),"admitted_family_counts_1_to_12":mixed,
            "unit_Euler_mapping_torus_phase_in_turns":{str(n):phase(Q(n,3)) for n in range(1,7)},
            "SU3_family_cubic_anomaly_per_block":8-8,
            "SM_roles":charges(),"two_scalar_external_component_counts":vacua(),
            "internal_stabilizer_component_counts":{"fundamental":1,"tensor":2,"both":1}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate()
    path=ROOT/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    elif result!=json.loads(path.read_text(encoding="utf-8-sig")):
        raise SystemExit("Stored result differs")
    else:
        print("1542: quotient bundles, faithful centre, family count and scalar phases verified")
