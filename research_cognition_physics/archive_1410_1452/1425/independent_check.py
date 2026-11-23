"""Distance from the physical scalar projector to the rate-allowed cone."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def run():
    src=json.loads((HERE.parent/"1386"/"inputs.json").read_text(encoding="utf-8-sig"))
    e,mu,tau=src["yukawa"]["charged_leptons"]
    e=e[0]-2*e[1]
    mu=mu[0]+2*mu[1]
    tau=tau[0]-2*tau[1]
    mass_factor=(tau*tau+e*e)/(tau*tau+mu*mu)
    amp,sig=src["DayaBay_amplitude"]
    x=(1-math.sqrt(1-amp-2*sig))/2
    atm=src["NO"]["atmospheric_sine_squared"]
    y=atm[0]+2*atm[2]
    old=[math.sqrt(1-x),math.sqrt(x*y),math.sqrt(x*(1-y))]
    def distance(ratio):
        k=math.sqrt(ratio*mass_factor*.99/1.01)
        normal=[k/math.hypot(k,1),-1/math.hypot(k,1),0]
        dot=sum(a*b for a,b in zip(normal,old))
        closest=[(a-dot*b)/math.sqrt(1-dot**2) for a,b in zip(old,normal)]
        assert math.isclose(closest[1],k*closest[0])
        return {"min_mu_over_e":k, "min_angle_degrees":math.degrees(math.asin(dot)),
                "closest_moduli":closest}
    # Check exact reconstruction from all three strictly positive rates, allowing phases.
    f=[math.sqrt(.29),math.sqrt(.37)*cmath.exp(.83j),math.sqrt(.34)*cmath.exp(-1.34j)]
    rij={(i,j):2.16*abs(f[i]*f[j].conjugate())**2 for i,j in itertools.combinations(range(3),2)}
    em,et,mt=rij[0,1],rij[0,2],rij[1,2]
    unnorm=[em*et,em*mt,et*mt]
    inferred=[val/sum(unnorm) for val in unnorm]
    assert all(math.isclose(p,abs(z)**2) for p,z in zip(inferred,f))
    return {"two_width_sign_only":distance(1.),"two_width_sign_and_et_upper":distance(1.25),
            "three_positive_mass_divided_rates":list(rij.values()),"inferred_port_weights":inferred}


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    out=run()
    target=HERE/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    else:
        assert equal(out,json.loads(target.read_text(encoding="utf8")))
        print("1425 direct cone and three-rate inverse verified")
