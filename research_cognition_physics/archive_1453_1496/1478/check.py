"""Exact arithmetic for the remaining vector charge and its finite-lattice test."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SCALE=10**21


def window(entry,k):
    c,s=F(entry["value"]),F(entry["width"])
    return c-k*s,c+k*s


def pair(xs):
    return list(map(str,xs))


def ceil(x):
    return -(-x.numerator//x.denominator)


def verify_charges(n):
    # Positive electron unit; W sets the common electroweak ladder interval.
    w=1-n
    u,d=(2-n)/3,(-1+2*n)/3
    e,nu=-1,-n
    assert 2*u+d == 1 and u+2*d == n and u-d == nu-e == w
    assert n==1+e-nu
    assert 2*(1-n)+2*e == -2*n
    assert 2*(1-n)+2*e-2*nu == 0
    h=w/2
    q,l=u-h,nu-h
    assert 3*q+l==0 and 2*q-u-d==0 and 2*l-e-nu==0
    assert 3*(2*q**3-u**3-d**3)+2*l**3-e**3-nu**3==0


def solve():
    p=json.loads((ROOT/"inputs.json").read_text(encoding="utf-8"))
    Z,A=p["sf6"]["Z"],p["sf6"]["A"]
    k=F(A-Z,A)
    rows={}
    for width in (1,2):
        neutron=window(p["neutron"],width)
        sf=window(p["sf6"],width)
        pred=(k*neutron[0],k*neutron[1])
        agreement=(max(neutron[0],sf[0]/k),min(neutron[1],sf[1]/k))
        assert agreement==neutron and sf[0]<=pred[0]<=pred[1]<=sf[1]
        nue=(F(A,Z)*(sf[0]-neutron[1]),F(A,Z)*(sf[1]-neutron[0]))
        dpe=(nue[0]+neutron[1],nue[1]+neutron[0])
        signs={}
        for label,edge,sign in (("negative",-agreement[0],-1),("positive",agreement[1],1)):
            minimal=ceil(SCALE/edge)
            # A particular cover where all quark charges also have integral weights.
            N=minimal+(2*sign-minimal)%3
            delta=F(sign,N)
            assert agreement[0]<=SCALE*delta<=agreement[1]
            verify_charges(delta)
            up,down=N*(2-delta)/3,N*(-1+2*delta)/3
            weak_step=N-sign
            for top,bottom in ((up,down),(-sign,-N),(weak_step,0)):
                assert top%1==bottom%1==0 and top-bottom==weak_step
                doublet_weight=top+bottom
                assert doublet_weight%2==weak_step%2
                assert (2*top)%2==(2*bottom)%2==0
            signs[label]={"min_colorless_electron_index":minimal,
                          "explicit_common_integral_index":N,"neutron_charge":str(delta),
                          "predicted_sf6_per_nucleon":str(k*delta)}
        rows[str(width)]={"universal_neutron_times_1e21":pair(agreement),
                          "sf6_prediction_from_neutron_times_1e21":pair(pred),
                          "free_beta_nue_times_1e21":pair(nue),
                          "free_p_plus_e_times_1e21":pair(dpe),"lattice":signs}
    for n in (F(0),F(-1,10**22),F(1,10**22)): verify_charges(n)
    return rows


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=cli.parse_args()
    results=solve()
    out=ROOT/"results.json"
    if args.write_new:
        with out.open("x",encoding="utf-8") as f:
            json.dump(results,f,indent=2)
            f.write("\n")
    else:
        assert results==json.loads(out.read_text(encoding="utf-8"))
        print("1478: charge conservation, mixed and cubic anomalies, working boxes and compact counterexamples verified")
