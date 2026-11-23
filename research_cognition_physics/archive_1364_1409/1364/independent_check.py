"""Exact root-count witnesses for the tree-level single-singlet mass completion."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np


def mul(x, y):
    corners = [a*b for a in x for b in y]
    return min(corners), max(corners)


def inv(x):
    assert x[0]*x[1] > 0
    return 1/x[1], 1/x[0]


def plus(x, y):
    return x[0]+y[0], x[1]+y[1]


def scale(x, factor):
    return mul(x, (factor, factor))


def value(coefficients, x):
    acc = (Fraction(0), Fraction(0))
    for a in reversed(coefficients):
        acc = plus(mul(acc, x), (a, a))
    return acc


def isolate(polynomial, lo, hi):
    lo, hi = Fraction(lo), Fraction(hi)
    a, b = value(polynomial,(lo,lo))[0], value(polynomial,(hi,hi))[0]
    assert a*b < 0
    for _ in range(56):
        half = (lo+hi)/2
        y = value(polynomial,(half,half))[0]
        if y*a < 0:
            hi = half
        else:
            assert y != 0
            lo = half
    return lo, hi


def decimal_box(box, places=8):
    factor = 10**places
    return [float((box[0]*factor).__floor__()/factor), float((box[1]*factor).__ceil__()/factor)]


def integer_data(mass):
    """All real and imaginary matrix operations here use a doubled integer matrix."""
    x, y = mass.real.astype(np.int64), mass.imag.astype(np.int64)
    block = np.block([[x, -y],[y, x]])
    H = block@block.T
    powers = [np.linalg.matrix_power(H, j) for j in range(5)]
    assert max(np.abs(p).max() for p in powers) < 10**8
    traces = [int(np.trace(p[:4, :4])) for p in powers]
    ch = [1]  # Descending characteristic-polynomial coefficients, via Newton's identities.
    for j in range(1,5):
        numerator = -sum(ch[j-k]*traces[k] for k in range(1,j+1))
        assert numerator % j == 0
        ch.append(numerator//j)
    offdiag = {}
    real_edges = []
    for i,j in [(0,1),(1,2),(2,0)]:
        real_edges.append(int(H[i,j]))
        assert H[4+i,j] == 0
        # For off-diagonals N of adj(sI-H), N=H^3+(s+c1)H^2+(s^2+c1*s+c2)H.
        poly = [powers[3]+ch[1]*powers[2]+ch[2]*H, powers[2]+ch[1]*H, H]
        offdiag[i,j] = ([int(p[i,j]) for p in poly], [int(p[4+i,j]) for p in poly])
    return ch[::-1], offdiag, real_edges


def root_certificate(mass, brackets):
    polynomial, edges, a = integer_data(mass)
    roots = [isolate(polynomial,*pair) for pair in brackets]
    # Four disjoint sign changes exhaust the roots of the monic quartic.
    derivative = [j*polynomial[j] for j in range(1,5)]
    def projected(s):
        recip = inv(value(derivative, s))
        return [(mul(value(r,s),recip), mul(value(i,s),recip)) for r,i in edges.values()]
    blocks = projected(roots[-1])
    linear = (Fraction(0),Fraction(0))
    quadratic = linear
    for i in range(3):
        j,k = (i+1)%3,(i+2)%3
        linear = plus(linear, scale(blocks[i][1], a[j]*a[k]))
        quadratic = plus(quadratic, scale(plus(mul(blocks[j][0],blocks[k][1]),
                                              mul(blocks[j][1],blocks[k][0])), a[i]))
    next_root = plus(roots[-1],scale(mul(linear,inv(quadratic)), -1))
    # J = Im(P_mass0[0,1] * conjugate(P_mass1[0,1])) from the light spectral projectors.
    real0, imag0 = projected(roots[0])[0]
    real1, imag1 = projected(roots[1])[0]
    quartet = plus(mul(imag0,real1), scale(mul(real0,imag1),-1))
    assert quartet[0]*quartet[1] > 0
    return ({"quartic_ascending": polynomial, "spectral_boxes": [decimal_box(z) for z in roots],
             "other_root_box": decimal_box(next_root), "light_quartet_box": decimal_box(quartet),
             "F_at_true_plus_1_box": decimal_box(plus(linear,quadratic)),
             "F_at_true_plus_4_box": decimal_box(plus(scale(linear,4),scale(quadratic,16)))},
            float(sum(next_root)/2))


def polar_fit(h, x):
    masses, basis = np.linalg.eigh(h)
    new = h + (x-masses[-1])*np.outer(basis[:,-1],basis[:,-1].conj())
    # These two examples have nonzero tree edges 01,02. In general use a spanning forest.
    phases = np.exp(1j*np.array([0,np.angle(new[0,1]),np.angle(new[0,2]),0]))
    new = phases[:,None]*new*phases[None,:].conj()
    assert abs(new[:3,:3].imag).max() < 2.e-11
    chol = np.linalg.cholesky(new[:3,:3].real)
    b = np.linalg.solve(chol, new[:3,3]).conj()
    m0 = np.sqrt(new[3,3].real - np.vdot(b,b).real)
    reconstructed = np.zeros((4,4),complex)
    reconstructed[:3,:3],reconstructed[3,:3],reconstructed[3,3] = chol, b, m0
    assert np.linalg.norm(new - reconstructed@reconstructed.conj().T) < 4.e-11
    assert abs(np.angle(np.linalg.det(reconstructed))) < 1.e-12
    return float(m0)


def run():
    cases = [
        ("registered", [[1,0,0],[1,2,0],[1,1,3]], [1,1j,1+1j], [(0,1),(3,4),(9,10),(23,24)]),
        ("two_physical_roots", [[1,0,0],[-1,2,0],[1,-1,3]], [-1+2j,-2+1j,-2+2j],
         [(0,1),(2,3),(10,11),(37,38)]),
    ]
    data = {}
    for key, m, b, brackets in cases:
        mass = np.zeros((4,4),complex)
        mass[:3,:3], mass[3,:3], mass[3,3] = m, b, 4
        H = mass@mass.conj().T
        eigen, basis = np.linalg.eigh(H)
        weak = basis[:3,:3]
        deficit = np.eye(3) - weak@weak.conj().T
        assert np.allclose(deficit, np.outer(basis[:3,-1],basis[:3,-1].conj()))
        report, other = root_certificate(mass, brackets)
        accepted = [eigen[-1]] + ([other] if other > eigen[2] else [])
        assert len(accepted) == (1 if key=="registered" else 2)
        report["reconstructed_positive_M0"] = [polar_fit(H,x) for x in accepted]
        if key=="registered":
            for shift in [1,4]:
                upper = H[:3,:3]+shift*deficit
                observed = (upper[0,1]*upper[1,2]*upper[2,0]).imag
                lower,higher = report[f"F_at_true_plus_{shift}_box"]
                assert lower <= observed <= higher < 0
        angle = .37
        right = mass @ np.diag([np.exp(1j*angle),1.,1.,1.])
        assert np.allclose(right@right.conj().T,H)
        assert np.isclose(np.angle(np.linalg.det(right)/np.linalg.det(mass)),angle)
        c = np.array([.021j,.013,-.044j])
        perturbed = mass.copy()
        perturbed[:3,3] = c
        z = np.asarray(b) @ np.linalg.solve(np.asarray(m),c) / 4
        assert np.allclose(np.linalg.det(perturbed)/np.linalg.det(mass),1-z)
        report["reverse_column_phase_example"] = float(np.angle(1-z))
        data[key] = report
    return data


def compare(actual,saved):
    if isinstance(actual,dict):
        assert actual.keys() == saved.keys()
        for key,item in actual.items():
            compare(item,saved[key])
    else:
        assert np.allclose(actual,saved,atol=1.e-12,rtol=1.e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args = parser.parse_args()
    result = run()
    target = Path(__file__).with_name("independent_results.json")
    if args.write:
        with target.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
    else:
        compare(result,json.loads(target.read_text(encoding="utf-8")))
    print("PASS", {k:(v["other_root_box"], v["light_quartet_box"]) for k,v in result.items()})
