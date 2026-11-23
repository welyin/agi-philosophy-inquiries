"""Check the two-current projection, kaon inequality, and an exact dark ds channel."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def spectrum(mass):
    _, singular, left_dagger = np.linalg.svd(mass.conj().T)
    left = left_dagger.conj().T[:, np.argsort(singular)]
    p = np.outer(left[-1].conj(), left[-1])
    z = left.conj().T @ np.diag([-1., -1., -1., 0.]) @ left
    np.testing.assert_allclose(z + np.eye(4), p, atol=1e-13)
    j = (left[0, 0]*left[1, 1]*left[0, 1].conj()*left[1, 0].conj()).imag
    return sorted(singular**2), float(abs(p[0, 1])), float(j)


def mass_matrix(m, c, b, t):
    return np.vstack([np.column_stack([m, c]), np.r_[t*b, 4*t]])


def calculate():
    old_limit = json.loads((HERE/"../1387/results.json").read_text("utf-8"))["projection_abs_upper"]
    couplings = [1.1e12, 1.6e12]
    kaon = {"FV_2016_22_2016_24_GeV": couplings,
            "crossings_GeV": [old_limit*x/2 for x in couplings],
            "joint_P_upper_for_fa_1e6_and_1e9_GeV": [min(old_limit, 2*f/couplings[-1]) for f in (1e6, 1e9)]}

    marker = np.array([1., 2., 2.])/3
    orthogonal = 2*np.outer(marker, marker)-np.eye(3)
    m = orthogonal @ np.diag([1., 2., 3.])
    c = orthogonal[:, 2]/2
    b = np.array([0, 1j, 1])
    zero_case = []
    for t in (1., 2., 11., 50.):
        full = mass_matrix(m, c, b, t)
        masses, p, j = spectrum(full)
        dark = np.r_[orthogonal[:, 0], 0.]
        np.testing.assert_allclose(full @ full.conj().T @ dark, dark, atol=1e-12)
        np.testing.assert_allclose(np.linalg.det(full), 23*t, atol=1e-12)
        assert abs(masses[0]-1) < 1e-11 and p < 1e-12 and abs(j) > .01
        zero_case.append({"t": t, "mass_squares": [float(i) for i in masses],
                          "weak_quartet": j, "direct_ds_exact": 0})
    w = np.array([0., -2j, 5.])
    eff = orthogonal @ (np.diag([1., 4., 9.25]) - np.outer(w, w.conj())/18) @ orthogonal.T
    # For these integral coefficients 2916 H_eff has Gaussian-integer entries.
    scaled = 2916*eff
    integer_matrix = np.rint(scaled.real).astype(int)+1j*np.rint(scaled.imag).astype(int)
    np.testing.assert_allclose(scaled, integer_matrix, atol=1e-10)
    triple = integer_matrix[0,1]*integer_matrix[1,2]*integer_matrix[2,0]
    invariant = Fraction(round(triple.imag), 2916**3)
    assert invariant != 0
    light, left = np.linalg.eigh(eff)
    limit_j = float((left[0,0]*left[1,1]*left[0,1].conj()*left[1,0].conj()).imag)
    assert min(np.diff(light)) > 1 and abs(limit_j) > .01

    # A different texture with a nonzero, decaying tree signal (the main benchmark).
    m0 = np.array([[1.,0.,0.],[1.,2.,0.],[1.,1.,3.]])
    b0 = np.array([1., 1j, 1+1j])
    decoupling = []
    for t in (1., 30., 300.):
        _, p, j = spectrum(mass_matrix(m0, np.zeros(3), b0, t))
        decoupling.append({"t": t, "t2_P_upper_abs": t*t*p, "weak_quartet": j})
    assert abs(decoupling[-1]["t2_P_upper_abs"]-.00178475426) < 1e-8
    assert abs(decoupling[-1]["weak_quartet"]-.00203369668) < 5e-9
    return {"kaon": kaon, "exact_tree_zero": zero_case,
            "zero_example_limit_mass_squares": light.tolist(),
            "zero_example_limit_weak_quartet": limit_j,
            "zero_example_limit_CP_odd_triplet_exact": str(invariant),
            "nonzero_decoupling": decoupling}


def verify(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            verify(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            verify(x, y)
    elif isinstance(a, (int, float)):
        np.testing.assert_allclose(a, b, rtol=3e-9, atol=2e-12)
    else:
        assert a == b


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    target = HERE/"independent_results.json"
    if args.write:
        with target.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        verify(result, json.loads(target.read_text("utf-8")))
        print("1393 independent verification passed")
