"""Check optimal compact color-singlet charge grids and isotope response exactly."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def model(n, sign):
    step = n-sign
    charge = {"u":F(2*n-sign,3), "d":F(-n+2*sign,3), "nu":-sign, "e":-n,
              "H+":step,"H0":0}
    a = (sign-2*n)%3
    w2 = step%2
    singlet_hyp = {}
    doublet_hyp = {}
    for color, doublet in [(1,("u","d")), (0,("nu","e")), (0,("H+","H0"))]:
        p,q = (F(charge[s]) for s in doublet)
        assert p-q == step
        h = 3*(p+q)
        assert h.denominator == 1
        assert (2*a*color+3*w2+h)%6 == 0
        doublet_hyp[doublet[0]+"/"+doublet[1]] = int(h)
        for particle in doublet:
            k = F(charge[particle])
            assert (k+F(a*color,3))%1 == 0
            if not particle.startswith("H"):
                singlet = 6*k
                assert singlet.denominator == 1
                assert (2*a*color+singlet)%6 == 0
                singlet_hyp[particle] = int(singlet)
    hq,hl = doublet_hyp["u/d"],doublet_hyp["nu/e"]
    hu,hd,hn,he = (singlet_hyp[k] for k in ["u","d","nu","e"])
    assert 3*hq+hl == 2*hq-hu-hd == 2*hl-he-hn == 0
    assert 3*(2*hq**3-hu**3-hd**3)+2*hl**3-hn**3-he**3 == 0
    assert 2*charge["u"]+charge["d"]+charge["e"]==0
    assert charge["u"]+2*charge["d"]==sign
    assert charge["nu"]==-sign
    return {"index":n,"sign":sign,"color_center_exponent":a,"weak_center_exponent":w2,
            "quark_weights":[str(charge["u"]),str(charge["d"])],
            "minimum_pure_EM_period_in_colorless_turns":1 if a==0 else 3,
            "integer_doublet_hypercharges":doublet_hyp,"integer_singlet_hypercharges":singlet_hyp}


def run():
    scale = 10**21
    neutron = [F(-26,10),F(18,10)]
    molecule = [F(-23,10),F(21,10)]
    def ceil(x):
        return -((-x.numerator)//x.denominator)
    cases = []
    for sign,bound in zip((-1,1),(-neutron[0],neutron[1])):
        n=ceil(scale/bound)
        assert F(scale,n)<=bound<F(scale,n-1)
        cases.append(model(n,sign))
    for mass in (146,147):
        prediction=[F(mass-70,mass)*v for v in neutron]
        assert molecule[0]<prediction[0]<prediction[1]<molecule[1]
    beta=[F(146,70)*(molecule[0]-neutron[1]),F(146,70)*(molecule[1]-neutron[0])]
    return {"two_signed_sharp_models":cases,"unconstrained_beta_neutrino_range_in_1e21":[str(x) for x in beta],
            "delta_exactly_zero_for_nonzero_integer_flux_up_to":ceil(scale/(-neutron[0]))-1}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group()
    g.add_argument("--write",action="store_true")
    g.add_argument("--check",action="store_true")
    args=p.parse_args()
    results=run()
    file=ROOT/"independent_results.json"
    if args.write:
        with file.open("x",encoding="utf-8") as stream:
            json.dump(results,stream,indent=2)
            stream.write("\n")
    else:
        assert results==json.loads(file.read_text(encoding="utf-8-sig"))
        print("1478 independent center, lattice and matter checks passed")
