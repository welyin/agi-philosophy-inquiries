"""Subtracted photon--family polarization from the charged mirror spectrum."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def integral(q, doublet, singlet, n):
    # Difference of logarithms, evaluated without a cancellation at small momentum.
    nodes, weights = np.polynomial.legendre.leggauss(n)
    x = (nodes+1)/2
    k = x*(1-x)
    arg = (1/doublet**2-1/singlet**2)*q**2*k/(1+q**2*k/singlet**2)
    return float(np.dot(weights/2,k*np.log1p(arg)))


def run():
    mass = read(HERE/"../1622/results.json")["witness_L_masses_GeV"]
    obs = read(HERE/"../1619/inputs.json")
    e2 = 4*math.pi/obs["alpha_inverse"]
    out = {}
    for r in [.5,1.,2.]:
        rows = {}
        asym = e2*math.log(r)/(12*math.pi**2)
        for q in [10.,200.,1000.,1e5]:
            integ = np.array([integral(q,m,r*m,96) for m in mass])
            for i,m in enumerate(mass):
                assert math.isclose(integ[i],integral(q,m,r*m,128),abs_tol=4e-9)
            # T_11 has charge 1/sqrt(2); relative electromagnetic vertex to g_F/sqrt(2).
            relative = e2*integ/(4*math.pi**2)
            central = sum(integ)/(4*math.pi**2*math.sqrt(6))
            assert all(abs(a)<=abs(asym)+1e-13 for a in relative)
            assert all(z*math.log(r)>=0 for z in relative)
            rows[str(q)]={"relative_electron_vector_vertex_shift":float(relative[0]),
                          "central_inverse_kinetic_shift":central,
                          "heavy_flavour_vector_vertex_shifts":relative[1:].tolist()}
        edges=[abs(rows[str(q)]["relative_electron_vector_vertex_shift"]) for q in [10.,200.,1000.,1e5]]
        assert edges==sorted(edges)
        out[str(r)]={"by_spacelike_Q_GeV":rows,"asymptotic_relative_vertex_shift":asym,
                     "minimum_max_at_0_and_200":edges[1]/2}
    return out


def compare(old,new):
    if isinstance(old,dict):
        return old.keys()==new.keys() and all(compare(old[k],new[k]) for k in old)
    if isinstance(old,list):
        return len(old)==len(new) and all(compare(x,y) for x,y in zip(old,new))
    return math.isclose(old,new,rel_tol=1e-8,abs_tol=1e-13)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    val=run()
    path=HERE/"independent_results.json"
    if cli.parse_args().check:
        assert compare(read(path),val)
        print("1624 PASS: photon-family loop, signed subtraction, paired threshold and asymptote")
    else:
        with path.open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(val,indent=2)+"\n")
