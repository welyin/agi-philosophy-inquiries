"""Independent determinant permutation, mass-pole interpolation, and finite stress tests."""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def complex_array(values):
    a=np.asarray(values,dtype=float)
    return a[...,0]+1j*a[...,1]


def permutation_det(mat):
    value=0j
    for perm in itertools.permutations(range(len(mat))):
        parity=sum(perm[i]>perm[j] for i,j in itertools.combinations(range(len(mat)),2))
        term=math.prod(mat[i,col] for i,col in enumerate(perm))
        value+=(-1)**parity*term
    return value


def check_mass(w):
    M=complex_array(w["weak_basis_mass_GeV"])
    det=permutation_det(M)
    # Independent of either Schur inversion or differentiating slogdet.
    cof=(-1)**(3+4)*permutation_det(np.delete(np.delete(M,3,axis=0),4,axis=1))
    k=cof/det
    assert abs(k.imag)<3e-15
    u,m,vh=np.linalg.svd(M)
    order=np.argsort(m)
    m=m[order]; u=u[:,order]; v=vh.conj().T[:,order]
    pole=(v[-1,:]*u[-2,:].conj())/m
    assert abs(k-sum(pole))<2e-12
    # Only b and two heavy poles participate in the aligned singlet row.
    k_light=pole[2]
    closed=k_light*np.prod(1-(m[2]/m[3:])**2)
    assert abs(closed-k)<3e-12
    for E in (0.01+.03j, -0.2+1j, -1.25/k, (0.3+.56j)/k):
        shifted=M.copy()
        shifted[3,4]+=E
        ratio=permutation_det(shifted)/det
        assert abs(ratio-(1+k*E))<1e-11

    # Tangent point to the circle of all complex perturbations, and a real sign flip.
    ang=1e-10
    r=math.sin(ang)
    z=r*np.exp(1j*(math.pi/2+ang))
    assert math.isclose(np.angle(1+z),ang,rel_tol=1.e-12)
    assert abs(abs(np.angle(1-1.25))-math.pi)<1e-12
    # Scale *all* one-Higgs entries, including the reverse entry, at fixed hard masses.
    E=(.02+.04j)/k
    shifted=M.copy(); shifted[3,4]=E
    hard=np.zeros_like(M)
    hard[3,3]=M[3,3]; hard[4,4]=M[4,4]
    def phase_at(logv):
        mat=hard+math.exp(logv)*(shifted-hard)
        return np.angle(permutation_det(mat)/det)
    h=1e-5
    slope=(phase_at(h)-phase_at(-h))/(2*h)
    expected=(2*k*E/(1+k*E)).imag
    assert abs(slope-expected)<1e-8

    # Finite, complex, determinant-positive wave function factors provide a control example.
    rng=np.random.default_rng(1637)
    H=rng.normal(size=(5,5))+1j*rng.normal(size=(5,5))
    ev,V=np.linalg.eigh((H+H.conj().T)*.02)
    L=(V*np.exp(ev))@V.conj().T
    assert abs(np.angle(permutation_det(L@M)/det))<1e-11
    return {"abs_inverse_entry_per_GeV":float(f"{abs(k):.12g}"),
            "abs_light_only_per_GeV":float(f"{abs(k_light):.12g}"),
            "relative_heavy_pole_shift":float(f"{abs(closed/k_light)-1:.9g}"),
            "all_phase_radius_GeV_at_1em10":float(f"{r/abs(k):.12g}"),
            "background_phase_derivative_test":float(f"{expected:.9g}")}


def run():
    rows=read("../1635/results.json")
    return {k:check_mass(rows[k]) for k in ("common_witness","high_mass_witness")}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=run()
    if args.check:
        assert result==read("independent_results.json")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")
    print("1637 independent determinant / pole residues / soft phase PASS")
