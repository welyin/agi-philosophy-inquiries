"""Replay charge embedding, B-L selection and neutrino-scale finite-defect bounds."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def run():
    cfg=read(HERE/"inputs.json")
    source=read(HERE/cfg["oscillation_input"])["juno_207_days"]
    v=cfg["vev_upper_GeV"]
    branch_masses={}
    spectra={}
    s,s_err=source["solar_21"]
    for order,key in (("NO","normal_atmospheric_31"),("IO","inverted_abs_31")):
        a,a_down,_=source[key]
        n=cfg["working_widths"]
        solar=s-n*s_err
        atm=a-n*a_down
        off=np.array([0.,solar,atm]) if order=="NO" else np.array([atm,atm+solar,0.])
        masses=np.sqrt(off)
        heavy,middle=np.sort(masses)[::-1][:2]
        norm1=float(np.sum(masses))
        y=math.sqrt(2)*heavy*1e-9/v
        eta=2*heavy*1e-9/v**2
        branch_masses[order]=masses
        spectra[order]=dict(lightest_zero_in_eV=masses.tolist(),min_spectral_mass_eV=heavy,
                            min_rank_one_distance_eV=middle,min_mass_sum_eV=norm1,
                            min_Dirac_largest_yukawa=y, min_Weinberg_norm_GeVinv=eta,
                            max_inverse_Weinberg_norm_GeV=1/eta,
                            min_16pi2_trace_gammaS_Dirac=10*float(masses@masses)*1e-18/v**2,
                            min_16pi2_trace_gammaS_seesaw={str(m0):10*m0*norm1*1e-9/v**2
                                                         for m0 in cfg["min_heavy_mass_examples_GeV"]})
    # A genuinely complex light Takagi basis. The two degenerate columns saturate
    # the nuclear-norm bound with no Majorana-phase or flavour-basis restriction.
    U=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/math.sqrt(3)
    for order,masses in branch_masses.items():
        diag=masses*1e-9
        source_mass=1e10
        nz=diag>0
        Y=1j*(math.sqrt(2*source_mass)/v)*U[:,nz]*np.sqrt(diag[nz])
        lo=-v**2/(2*source_mass)*(Y@Y.T)
        target=(U*diag)@U.T
        assert np.allclose(lo,target,atol=1e-25,rtol=2e-12)
        assert math.isclose(5*np.linalg.norm(Y,"fro")**2,
            spectra[order]["min_16pi2_trace_gammaS_seesaw"][str(cfg["min_heavy_mass_examples_GeV"][0])])
    components={"Q":(1,F(1,6),F(1,3)),"uc":(1,F(-2,3),F(-1,3)),
                "ec":(1,1,1),"dc":(-3,F(1,3),F(-1,3)),"L":(-3,F(-1,2),-1),
                "Hu":(-2,F(1,2),0),"Hd":(2,F(-1,2),0)}
    for x,y,bl in components.values():
        assert x+4*y==5*bl
    X={"T":1,"F":-3,"H":-2,"Hb":2,"S":0}
    monomials=["T T H","T F Hb","Hb S H","Hb H","S S S","S S"]
    assert all(sum(X[f] for f in mon.split())==0 for mon in monomials)
    mixed=3*(F(3,2)-3*F(1,2))+4*(-2*F(1,2)+2*F(1,2))
    assert mixed==0 and 3*F(-1,3)+2*F(1,2)==0
    return dict(hypercharge_fundamental=["-1/3"]*3+["1/2"]*2, total_b5_three_families_even_adj_four_pairs=0,
                extra_dynkin_index={"odd_adjoint":5,"15_and_conjugate":7,"true_gauge_singlets":0},
                budget_after_three_families_and_one_adjoint=4,
                mixed_SU5_squared_X=str(mixed),X_of_Weinberg_operator=2*(X["F"]+X["H"]),
                low_energy_mass_work_boxes=spectra)


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    if isinstance(a,float):
        return math.isclose(a,b,rel_tol=2e-12,abs_tol=0.)
    return a==b


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    group=cli.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write",action="store_true")
    args=cli.parse_args()
    result=run()
    dst=HERE/"results.json"
    if args.write:
        with dst.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        assert same(result,read(dst))
        print("1481: SU(5) charges, Dirac/Majorana bounds and rank-two equality witnesses verified.")
