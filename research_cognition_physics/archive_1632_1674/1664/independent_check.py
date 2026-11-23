"""Independent checks on positive 15 mass blocks, curvature and the scalar vacuum."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def interpolate(delta, mixed, sextet):
    cross = (delta + sextet - 2*mixed)/25
    linear = (sextet - mixed)/5 - 2*cross
    common = mixed + linear + 6*cross
    return common, linear, cross


def evaluate(coeff, x, y):
    common, linear, cross = coeff
    return common + (x+y)*linear + x*y*cross


def logdet(block):
    return np.linalg.slogdet(block)[1]


def positive_matrices_test():
    # Complex positive Hermitian blocks whose common diagonalization is impossible.
    factors = [np.array([[1, 2j], [.5, 1]], dtype=complex),
               np.array([[.3j, .1], [-.2, .4j]], dtype=complex),
               np.array([[2, 1], [1j, 3]], dtype=complex)]
    delta, mixed, sextet = [z@z.conj().T + .1*np.eye(2) for z in factors]
    assert np.linalg.norm(delta@sextet - sextet@delta) > .1
    coeff = interpolate(delta, mixed, sextet)
    # Test on the normalized symmetric two-index space, including a complex-rotated adjoint.
    rng = np.random.default_rng(1664)
    unitary = np.linalg.qr(rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)))[0]
    sigma = (unitary*np.array([2,2,2,-3,-3]))@unitary.conj().T
    pair_basis = []
    eye = np.eye(5)
    edges = [(i,j) for i in range(5) for j in range(i,5)]
    for i,j in edges:
        xy = np.kron(eye[i],eye[j]) + np.kron(eye[j],eye[i])
        pair_basis.append(xy/np.linalg.norm(xy))
    V = np.column_stack(pair_basis)
    w = V.conj().T@np.kron(unitary,unitary)@V
    one = V.T@(np.kron(sigma,eye)+np.kron(eye,sigma))@V
    two = V.T@np.kron(sigma,sigma)@V
    a,b,c = coeff
    projected = np.kron(np.eye(15),a) + np.kron(one,b) + np.kron(two,c)
    reference = np.zeros_like(projected)
    for k,(i,j) in enumerate(edges):
        target = sextet if j<3 else (delta if i>=3 else mixed)
        reference[2*k:2*k+2,2*k:2*k+2] = target
    U = np.kron(w,np.eye(2))
    np.testing.assert_allclose(U.conj().T@projected@U,reference,atol=2e-13)
    x = (logdet(delta)-logdet(mixed))/2
    y = (logdet(sextet)-logdet(mixed))/2
    L = np.linalg.cholesky(mixed)
    inv = np.linalg.inv(L)
    curvature = np.linalg.eigvalsh(inv@(delta+sextet-2*mixed)@inv.conj().T)
    amplitude_bound = np.exp(x)+np.exp(y)-2  # n = 2
    assert curvature[-1] >= amplitude_bound
    midpoint = (delta+sextet)/2
    mid_sum = (logdet(delta)+logdet(sextet))/2 - logdet(midpoint)
    assert mid_sum < 0
    # Test the analytic global quartic completion even far outside the chosen vacuum orbit.
    K = .4*max(np.linalg.norm(b+u*c,2) for u in (-3,2)) + np.linalg.norm(c,2)/2
    for scale in (.01,1,10,1e3):
        for _ in range(50):
            eig = rng.normal(size=5)*scale
            eig -= eig.mean()
            W = np.sum(((eig+3)*(eig-2))**2)
            minima = [np.linalg.eigvalsh(evaluate(coeff,eig[i],eig[j]))[0] for i,j in edges]
            assert min(minima) >= -2*K*math.sqrt(W)-1e-12
    return {"adjoint_15_multiplicities_D_B_S": [3,6,6],
            "single_leg_midpoint_X_plus_Y": mid_sum,
            "quadratic_curvature_upper_eigenvalue": curvature[-1],
            "curvature_trace_AMGM_lower_bound": amplitude_bound,
            "sufficient_quartic_K": K}


def compute():
    prior = read("../1663/results.json")
    both = read("../1663/independent_results.json")
    X = min(prior["sum_log_triplet_over_doublet"][0],both["whole_box_X"][0])-1e-5
    Y = min(prior["sum_log_sextet_over_doublet"][0],both["whole_box_Z"][0])-1e-5
    lower_bound = math.exp(X/2)+math.exp(Y/2)-2
    central = prior["central_witness"]
    # Exact rational reconstruction is essential when a squared doublet mass lies 27 digits below Delta.
    scale = Fraction(10**16)
    delta = (Fraction(str(central["common_triplet_mass_GeV"]))/scale)**2
    sextet = (Fraction(str(central["common_sextet_mass_GeV"]))/scale)**2
    mixed = [(Fraction(str(v))/scale)**2 for v in central["doublet_masses_GeV"]]
    coeff = [interpolate(delta,v,sextet) for v in mixed]
    for v, abc in zip(mixed,coeff):
        assert [evaluate(abc,*e) for e in [(-3,-3),(-3,2),(2,2)]] == [delta,v,sextet]
    linear = np.array([float(row[1]) for row in coeff])
    cross = np.array([float(row[2]) for row in coeff])
    K = .4*max(max(abs(linear+u*cross)) for u in (-3,2))+max(abs(cross))/2
    actual = float(max((delta+sextet-2*v)/v for v in mixed))
    assert actual > lower_bound
    return {"four_copy_full_box_curvature_lower_bound": lower_bound,
            "four_copy_full_box_minimum_generalized_B_over_endpoint_mean_upper_bound": 2/(2+lower_bound),
            "central_witness_exact_rational_mass_squared_reconstruction": True,
            "central_witness_curvature_upper_eigenvalue": actual,
            "central_quartic_K_in_units_MU": K, "central_allowed_example_lambda": 1,
            "central_lambda_minus_K_squared": 1-K*K,
            "noncommuting_test": positive_matrices_test()}


def stable(x):
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.10g}")
    return x


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args = cli.parse_args()
    ans = stable(compute())
    if args.check:
        assert ans == read("independent_results.json")
        print("1664 independent PASS: no commuting assumption, 15-space projection, exact benchmark, global quartic bound")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as f:
            f.write(json.dumps(ans,indent=2)+"\n")
