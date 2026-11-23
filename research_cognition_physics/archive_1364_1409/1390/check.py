"""Finite-area Dalitz inversion and dimension-six/eight potential competitors."""
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def potential_coefficients():
    # h/v is the expansion variable, all potentials divided by m^2 v^2.
    y = np.polynomial.Polynomial([0.,1.,0.5])
    def derivatives(polynomial):
        return [round(float(polynomial.deriv(n)(0)),11) for n in (2,3,4)]
    sm = y**2/2
    dim6 = sm+y**3/3
    # f^4 log(f^2) series to order four in h/v.
    phi=np.polynomial.Polynomial([1.,1.])
    log=np.polynomial.Polynomial([0.]+[2*(-1)**(n+1)/n for n in range(1,5)])
    gw=phi**4*(log-0.5)/8
    return {k:derivatives(p) for k,p in {"SM":sm,"GW":gw,"dim6":dim6,"dim8":dim6+y**4/6}.items()}


def moments(energy, center, halfwidth, order):
    nodes, weights=np.polynomial.legendre.leggauss(order)
    q = center+halfwidth*nodes[:,None]
    w = halfwidth*nodes[None,:]
    # Use a lower bound throughout the box on the maximum allowed |w|.
    qlo, qhi=center-halfwidth,center+halfwidth
    assert 4 < qlo and qhi < (math.sqrt(energy)-1)**2
    lam=(energy-qhi-1)**2-4*qhi
    assert math.sqrt(lam*qlo*(qlo-4))/(2*qhi)>halfwidth
    average=(energy+3-q)/2
    shape=1/(q-1)+1/(average+w-1)+1/(average-w-1)
    measure=np.outer(weights,weights)/4
    return float((measure*shape).sum()),float((measure*shape**2).sum())


def rate(a,b_eff,mu,nu):
    return b_eff**2 + 2*b_eff*a*a*mu + a**4*nu


def recover(a, pair, r):
    mu1,nu1=pair[0]
    mu2,nu2=pair[1]
    return ((r[1]-r[0])-a**4*(nu2-nu1))/(2*a*a*(mu2-mu1))


def certificate():
    spec=json.loads((HERE/"inputs.json").read_bytes())
    potentials=potential_coefficients()
    assert potentials=={"SM":[1,3,3],"GW":[1,5,11],"dim6":[1,5,15],"dim8":[1,5,19]}
    by_energy={}
    for energy in spec["s_in_mass_squared"]:
        centers=[(energy+3)/3,spec["asymmetric_pair_mass_squared"]]
        halfwidth=0.05
        pair=[moments(energy,q,halfwidth,12) for q in centers]
        assert np.allclose(pair,[moments(energy,q,halfwidth,20) for q in centers],rtol=0,atol=2e-15)
        point_shape=[1/(q-1)+2/((energy+3-q)/2-1) for q in centers]
        trials={}
        for label,a,b,d in [("SM",3.,3.,0.),("GW",5.,11.,0.),("dim6",5.,15.,0.),
                            ("dim6_with_contact",5.,15.,float(F(-4,15)))]:
            b_eff=b+d*(energy-1)
            rates=[rate(a,b_eff,mu,nu) for mu,nu in pair]
            inferred=recover(a,pair,rates)
            assert abs(inferred-b_eff)<1e-11
            factor=sum(rates)/(2*a*a*abs(pair[1][0]-pair[0][0]))
            trials[label]={"rates":rates,"effective_quartic":inferred,
                           "one_percent_worst_abs_error":spec["relative_reduced_rate_test"]*factor}
        by_energy[str(energy)]={"point_shape":point_shape,"moments":pair,"labels":trials}
    e1,e2=spec["s_in_mass_squared"]
    e_values=[by_energy[str(e)]["labels"]["dim6_with_contact"]["effective_quartic"] for e in (e1,e2)]
    contact=(e_values[1]-e_values[0])/(e2-e1)
    quartic=e_values[0]-(e1-1)*contact
    assert math.isclose(contact,-4/15) and math.isclose(quartic,15)
    return {"m2_cubic_quartic_units":potentials,
            "hh_elastic_threshold_units":{k:5*a*a/3-b for k,(_,a,b) in potentials.items()},
            "mu_source_at_reduced_threshold":{k:[a/3,(b+a*a)/8] for k,(_,a,b) in potentials.items()},
            "bin_halfwidth_mass_squared":0.05,"dalitz_bins":by_energy,
            "source_contact_recovered_at_two_energies":{"quartic":quartic,"contact":contact},
            "both_cubic_templates_inside_combined_interval":all(spec["hh_combined_ATLAS_CMS"]["interval"][0]<x<spec["hh_combined_ATLAS_CMS"]["interval"][1] for x in (1,5/3))}


def agrees(a,b):
    if isinstance(a,dict):
        return isinstance(b,dict) and a.keys()==b.keys() and all(agrees(v,b[k]) for k,v in a.items())
    if isinstance(a,(tuple,list)):
        return isinstance(b,(tuple,list)) and len(a)==len(b) and all(agrees(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12) if isinstance(a,float) else a==b


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    result=certificate()
    target=HERE/"results.json"
    if p.parse_args().write:
        with target.open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2)
            fp.write("\n")
    else:
        assert agrees(result,json.loads(target.read_bytes()))
        print("1390: invariant potential derivatives and finite Dalitz inversions verified")
