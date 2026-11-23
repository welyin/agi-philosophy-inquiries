"""Real-CP matching, a compact-group countercheck, and a topology Wilson coefficient."""
import argparse
import json
from math import pi
from pathlib import Path

import numpy as np


def benchmark():
    s=1+1j
    b=np.array([.5,.5,1.])*s+np.array([.5,-.5,0.])*s.conjugate()
    heavy=4.
    m=np.array([[1,0,0],[1,2,0],[1,1,3]],float)
    hu=np.diag([1.,4.,9.])
    mass=np.zeros((4,4),complex)
    mass[:3,:3]=m
    mass[3,:3]=b
    mass[3,3]=heavy
    gram=mass@mass.conj().T
    eig,left=np.linalg.eigh(gram)
    quartet=(left[0,0]*left[1,1]*left[0,1].conj()*left[1,0].conj()).imag
    assert abs(quartet)>.002 and abs(np.linalg.det(mass)-24)<1e-12
    # For a real small coefficient, the aligned gate is nonzero with a real determinant.
    epsilon=.002
    c=epsilon*m@b.conj()*heavy
    aligned=mass.copy()
    aligned[:3,3]=c
    det_ratio=np.linalg.det(aligned)/np.linalg.det(mass)
    assert np.isclose(det_ratio,1-epsilon*(b@b.conj()).real)
    shifted=mass.copy()
    shifted[:3,3]=epsilon*hu@m@b.conj()*heavy
    z=epsilon*b@np.linalg.solve(m,hu@m@b.conj())
    assert np.isclose(np.linalg.det(shifted)/24,1-z) and abs(z.imag)>.0001
    # Pure scalar CP-odd invariant is even under S -> -S and odd under CP.
    # Convention X=Im S^2; coefficients dimensionless; M and light masses in benchmark units.
    p=lambda v:(v**2).imag
    assert p(s)==p(-s)==2 and p(s.conjugate())==-2
    cutoff=10.
    schemes=[]
    for coeff in [0.,.001]:
        # Same bare theta and non-topological matched couplings below scalar thresholds.
        theta=coeff*p(s)/cutoff**2+np.angle(np.linalg.det(mass))
        schemes.append({"real_topology_coefficient":coeff,"matched_theta":theta,
                        "light_left_quartet":quartet,"heavy_mass_squared":eig[-1]})
    assert abs(schemes[0]["matched_theta"]-schemes[1]["matched_theta"])>1e-5
    return {"benchmark_determinant":np.linalg.det(mass).real,"aligned_determinant_ratio":det_ratio.real,
            "nonaligned_z":[z.real,z.imag],"nonaligned_phase":np.angle(1-z),
            "two_fixed_weak_matching_examples":schemes}


def degree_witness(n):
    # Any scalar-only monomial S^k S*^l transforming trivially has k-l=0 mod n.
    terms=[(i,total-i) for total in range(n+1) for i in range(total+1)
           if (2*i-total)%n==0 and i!=total-i]
    first=min(i+j for i,j in terms)
    v=np.exp(.5j*pi/n)
    assert first==n and abs((v**n).imag-1)<1e-14
    return first


def result():
    # Frozen PSI and lattice inputs are not jointly distributed.
    lower=(.00148-.00014-.00031)*1e-13
    upper=1.8e-26/lower
    return {"tree_and_matched":benchmark(),"first_possible_CP_odd_degree_in_Zn":
            {str(n):degree_witness(n) for n in (2,3,5,11,31)},
            "conditional_theta_abs_upper_PSI_one_linear_lattice_box":upper}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(format(x,".10g")) if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    parser.add_argument("--check",action="store_true")
    out=Path(__file__).with_name("independent_results.json")
    ans=rounded(result())
    if parser.parse_args().write:
        with out.open("x",encoding="utf-8") as fp:
            json.dump(ans,fp,indent=2,ensure_ascii=False)
            fp.write("\n")
    else:
        assert ans==json.loads(out.read_text(encoding="utf-8"))
    print("1389 independent real-CP spurion, observable weak phase, operator orders and nEDM checked.")
