"""Euler-operator certificates and Lorentz-scalar finite cross sections."""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def value(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def moments(powers):
    # Lagrange weights evaluate polynomials at p=1 on the K power nodes.
    weights=[math.prod(Fraction(1-other,p-other) for other in powers if other!=p) for p in powers]
    derivatives=[sum(w*math.prod(p-i for i in range(k)) for w,p in zip(weights,powers))
                 for k in range(len(powers)+1)]
    target=[int(k<2) for k in range(len(powers)+1)]
    assert derivatives[:-1]==target[:-1]
    gap=derivatives[-1]-target[-1]
    assert gap==-math.prod(1-p for p in powers)
    return int(gap)


def minkowski(a,b):
    return a[...,0]*b[...,0]-np.sum(a[...,1:]*b[...,1:],axis=-1)


def rate(energy,m,h,v,kappa,eta,conversion):
    E=energy/2
    pin,pout=math.sqrt(E**2-m**2),math.sqrt(E**2-h**2)
    p1=np.array([E,0,0,pin])
    p2=np.array([E,0,0,-pin])
    # Vectorized quadrature in polar angle, endpoint-safe independent of author's gamma trace.
    theta=np.linspace(0,np.pi,8193)
    q=np.stack([np.full_like(theta,E),pout*np.sin(theta),np.zeros_like(theta),pout*np.cos(theta)],axis=1)
    q2=np.array([2*E,0,0,0])-q
    a=np.full_like(theta,m/v**2*(eta+kappa*3*h**2/(energy**2-h**2)))
    b=np.zeros_like(q)
    for outgoing in (q,q2):
        internal=p1-outgoing
        coeff=(m*kappa/v)**2/(minkowski(internal,internal)-m*m)
        a+=coeff*m
        b+=coeff[:,None]*internal
    spin=a*a*(minkowski(p2,p1)-m*m)+2*a*m*minkowski(p2-p1,b)
    spin+=2*minkowski(p1,b)*minkowski(p2,b)-(minkowski(p1,p2)+m*m)*minkowski(b,b)
    weighted=spin*np.sin(theta)
    steps=np.diff(theta)
    integral=float(np.sum(steps*(weighted[:-1]+weighted[1:])/2))
    return conversion*integral*pout/pin/(64*math.pi*energy**2)


def run():
    src=value(HERE/"inputs.json")
    muin=value(HERE/src["shared"]["muon_mass_and_higgs_mass"])
    m=muin["masses"]["muon_MeV"]/1000.
    mh=muin["cms"]["higgs_mass_GeV"]
    v=(math.sqrt(2)*value(HERE/src["shared"]["fermi_constant"])["GF"])**-.5
    rates={}
    authors=value(HERE/"results.json")
    for E in src["muon_collider_prospective_energies_GeV"]:
        rates[str(E)]={}
        for key,k,e in [("linear",1,0),("cubic",3,6),("two_aligned",1,-8),("three_aligned",1,0)]:
            result=rate(E,m,mh,v,k,e,src["conversion_fb_per_GeV_minus2"])
            assert math.isclose(result,authors["unpolarized_muon_direct_hh_born_fb"][str(E)][key],rel_tol=1e-6)
            rates[str(E)][key]=result
    nearest={str(k):moments(list(range(3,3+2*k,2))) for k in range(1,9)}
    wider=moments([3,7,11,15])
    return {"first_new_derivative_at_closest_odd_powers":nearest,"first_new_derivative_3_7_11_15":wider,
            "direct_muon_pair_to_hh_crosscheck_fb":rates}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=run()
    dest=HERE/"independent_results.json"
    if args.write_new:
        with dest.open("x",encoding="utf-8") as handle:
            json.dump(data,handle,indent=2);handle.write("\n")
    elif args.check:
        def close(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(close(a[k],b[k]) for k in a)
            return math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-14)
        assert close(data,value(dest))
        print("1436 independent Lorentz trace and exact rational derivative theorem verified")
    else:
        print(json.dumps(data,indent=2))
