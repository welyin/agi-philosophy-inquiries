"""Independent finite scalar responses and simultaneous real-axis deformations."""
import argparse
import itertools
import json
from math import atan, isclose, sqrt
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def normalized(a):
    return a/np.linalg.norm(a)


def chiral_mass_basis(y):
    # Recover right singular kets independently from the left Hermitian eigensystem.
    square,L=np.linalg.eigh(y@y.conj().T)
    masses=np.sqrt(square)
    R=(y.conj().T@L)/masses
    np.testing.assert_allclose(R.conj().T@R,np.eye(3),atol=2e-10)
    np.testing.assert_allclose(L.conj().T@y@R,np.diag(masses),atol=2e-12)
    return masses,L,R


def flavor_fraction(coupling,mass,parent_mass):
    scores=np.zeros((3,3))
    for i,j in itertools.product(range(3),repeat=2):
        a,b=mass[i]/parent_mass,mass[j]/parent_mass
        if a+b>=1:
            continue
        gij,gji=coupling[i,j],coupling[j,i]
        phase=sqrt((1-(a+b)**2)*(1-(a-b)**2))
        spin=(1-a*a-b*b)*(abs(gij)**2+abs(gji)**2)-4*a*b*(gij*gji).real
        scores[i,j]=phase*spin
    return float((scores.sum()-np.trace(scores))/scores.sum())


def calculate():
    cfg=read("inputs.json")
    inputs=cfg["witness_in_units_of_yb"]
    A=inputs["A"]
    d=np.array(inputs["direct"])
    u=normalized(np.array(inputs["u"]))
    seed=np.array(inputs["v_seed"])
    constraint=normalized(u/d)
    v=normalized(seed-constraint*np.dot(seed,constraint))
    X=-1j*A*np.outer(u,v)
    D=np.diag(d)
    y=D+X
    m,L,R=chiral_mass_basis(y)
    w=d*v
    # Vandermonde determinant in the three pairs (u_i,w_i), up to ordering sign.
    minor=np.linalg.det(np.column_stack([u**2,u*w,w**2]))
    h=y@y.conj().T
    I=(h[0,1]*h[1,2]*h[2,0]).imag
    assert isclose(abs(I),A**3*abs(minor),rel_tol=1e-10)
    J=(L[0,0]*L[1,1]*L[0,1].conj()*L[1,0].conj()).imag
    assert isclose(abs(I),abs(J)*np.prod([m[j]**2-m[i]**2 for i,j in ((0,1),(0,2),(1,2))]),rel_tol=1e-10)
    q=-A*u*v/d
    assert abs(q.sum())<1e-12 and isclose(sum(abs(q)),2*max(abs(q)),rel_tol=1e-12)
    eps=cfg["inverse_direct_independent_fraction"]
    corners=[np.angle(np.linalg.det(np.diag(d/(1+eps*np.array(s)))+X))
             for s in itertools.product([-1,1],repeat=3)]
    worst=atan(eps*sum(abs(q)))
    np.testing.assert_allclose([min(corners),max(corners)],[-worst,worst],atol=2e-14)

    G=L.conj().T@X@R
    signed=np.diag(G)/m
    np.testing.assert_allclose(signed.sum(),0,atol=2e-10)
    assert isclose(np.dot(m*m,signed.real),A*A,rel_tol=1e-10)
    step=1e-5
    left=np.linalg.svd(D+(1-step)*X,compute_uv=False)[::-1]
    right=np.linalg.svd(D+(1+step)*X,compute_uv=False)[::-1]
    np.testing.assert_allclose((np.log(right)-np.log(left))/(2*step),signed.real,atol=4e-8)
    strength=abs(G/A)**2
    fraction=1-np.trace(strength)
    assert isclose(strength.sum(),1)
    assert isclose(fraction,flavor_fraction(G,0*m,1.),rel_tol=1e-13)
    # Same scalar reference and physical point, but an extra even X^2 direct interaction.
    aligned=L.conj().T@(D+X)@R
    assert abs(flavor_fraction(aligned,0*m,1.))<1e-14
    for x in (-1.2,-.5,.5,1.,1.2):
        assert abs(np.angle(np.linalg.det((1+x*x)/2*D+x*X)))<3e-14
    weak=read(cfg["weak_matching_inputs"])
    yb=weak["yukawa"]["down"][-1][0]
    masses_GeV=m*yb*246.22/sqrt(2)
    finite_R=flavor_fraction(G,masses_GeV,91.1876)
    return {"simultaneous_five_percent_max_theta_rad":worst,
            "inverse_direct_signed_coefficients":q.tolist(),"J":float(J),
            "relative_mass_slopes":signed.real.tolist(),"relative_diagonal_phase_slopes":signed.imag.tolist(),
            "light_scalar_linear_portal_chiral_fractions":strength.tolist(),
            "different_flavor_fraction_open_massless":float(fraction),
            "tree_level_91GeV_different_flavor_fraction":finite_R,
            "countermodel_different_flavor_fraction":0.}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    return np.allclose(a,b,rtol=2e-7,atol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        assert same(result,read("independent_results.json"))
        print("1414 mass determinant, flavor response and exact deformation bounds verified")
