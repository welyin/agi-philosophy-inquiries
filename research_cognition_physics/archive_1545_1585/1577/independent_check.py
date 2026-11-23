"""Use charged and neutral SU(5) tensors directly in the parent 15."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
N = 5


def block(i, j, x):
    t = np.zeros((3*N, 3*N), dtype=complex)
    t[i*N:(i+1)*N,j*N:(j+1)*N] = x
    if i != j:
        t[j*N:(j+1)*N,i*N:(i+1)*N] = x.conj().T
    return t


def gram(operators, adj, sym):
    basis = [t/np.linalg.norm(t) for t in operators]
    assert np.allclose([[np.vdot(x,y).real for x in basis] for y in basis], np.eye(len(basis)))
    columns = [np.r_[(t@adj-adj@t).ravel(), (t@sym+sym@t.T).ravel()] for t in basis]
    return np.array([[np.vdot(x,y).real for y in columns] for x in columns])


def round_spectrum(values):
    return [round(float(x), 10) for x in values]


def mass_point(a, b1, b2):
    eye = np.eye(N)
    adj = a*(block(2,2,eye)-block(0,0,eye))
    sym = b1*block(0,1,eye) + b2*block(2,1,eye)
    f = np.diag([1.,-1.,0.,0.,0.])
    anti = np.zeros((N,N))
    anti[0,1], anti[1,0] = 1, -1
    spectra = []
    for tensor in (anti, f):
        for phase in (1,1j):
            cm = gram([block(i,1,phase*tensor) for i in (0,2)], adj,sym)
            vals = np.linalg.eigvalsh(cm)
            if phase == 1:
                spectra.append(vals)
            else:
                assert np.allclose(vals, spectra[-1])
    aa, bb = a*a, b1*b1+b2*b2
    assert np.allclose(spectra[0], [aa,aa+bb])
    assert np.allclose(spectra[1], [aa+bb,aa+4*bb])
    real = [block(i,i,f) for i in (0,1,2)] + [block(0,2,f)]
    vals = np.linalg.eigvalsh(gram(real+[block(0,2,1j*f)], adj,sym))
    assert np.isclose(vals[0], 0)
    im = np.argmin(abs(vals-(4*aa+bb)))
    assert np.isclose(vals[im],4*aa+bb)
    roots = np.delete(vals,[0,im])
    assert np.isclose(sum(roots),4*aa+5*bb)
    assert np.isclose(sum(x*y for x,y in itertools.combinations(roots,2)),4*bb**2+16*aa*bb)
    ratio = np.prod(roots)/(12*aa*bb**2)
    assert np.isclose(ratio,(2*b1*b2/bb)**2)
    return {"inputs_g_equals_one": [a,b1,b2], "ten": round_spectrum(spectra[0]),
            "fifteen": round_spectrum(spectra[1]), "adjoint_real_squared_masses": round_spectrum(roots),
            "adjoint_imag_squared_mass": float(4*aa+bb), "bridge_product": round(float(ratio),10),
            "fractional_asymmetry": round(abs(b1*b1-b2*b2)/bb,10)}


def epsilon_witness():
    eps = np.zeros([5]*5)
    for p in itertools.permutations(range(5)):
        inv = sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))
        eps[p] = (-1)**inv
    q = np.eye(5)[:,0]
    x, y = np.zeros((5,5)), np.zeros((5,5))
    x[1,2], x[2,1], y[3,4], y[4,3] = 1,-1,1,-1
    def overlap(q,x,y):
        return np.einsum("ijklm,i,jk,lm",eps,q,x,y)
    initial = overlap(q,x,y)
    rng = np.random.default_rng(1577)
    u, _ = np.linalg.qr(rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)))
    u[:,0] /= np.linalg.det(u)
    assert np.isclose(initial,4.) and np.isclose(overlap(u@q,u@x@u.T,u@y@u.T),initial)
    assert np.isclose(overlap(q,x+y,x+y),8.)
    assert all((1+k)%5 for k in (0,2,-2))
    assert (1+2+2)%5 == 0
    return {"singlet_contraction_distinct_states":int(initial),
            "singlet_contraction_identical_superposed_states":8, "minimum_broken_vectors":2}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write",action="store_true")
    args = parser.parse_args()
    points = [mass_point(*p) for p in [(1.,1.,1.),(1.1,1.,1.),(2.,1.,1.),(1.,1.,2.),(1.,2.,1.)]]
    assert points[3]["adjoint_real_squared_masses"] == points[4]["adjoint_real_squared_masses"]
    data = {"round":1577, "mass_squared":points, "static_charge":epsilon_witness(),
            "embedding_index":3, "pi3_coset_order":3, "g5_over_g15":1/math.sqrt(3)}
    path = HERE / "independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as fh:
            json.dump(data,fh,indent=2,sort_keys=True,allow_nan=False)
            fh.write("\n")
    elif data != json.loads(path.read_text(encoding="utf-8")):
        raise SystemExit("independent result mismatch")
    print("1577 independent check passed")


if __name__ == "__main__":
    main()
