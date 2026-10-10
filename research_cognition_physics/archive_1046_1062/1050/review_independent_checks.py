"""Independent 1050 review: rational bounds, direct Slater quadrature and CAR.

No import of author science.  This is not a Coulomb time-evolution simulation.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import itertools
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
RESULT = HERE / "review_independent_results.json"


def lift_one_body(matrix):
    basis = [n for n in range(16) if n.bit_count() == 2]
    index = {n: k for k, n in enumerate(basis)}
    result = np.zeros((6, 6), complex)
    for col, bits in enumerate(basis):
        for p, q in itertools.product(range(4), repeat=2):
            if not (bits >> q) & 1:
                continue
            sign = (-1) ** ((bits & ((1 << q) - 1)).bit_count())
            temp = bits ^ (1 << q)
            if (temp >> p) & 1:
                continue
            sign *= (-1) ** ((temp & ((1 << p) - 1)).bit_count())
            result[index[temp | (1 << p)], col] += sign * matrix[p, q]
    return result


def slater_review(R):
    # Prolate coordinates, with Laguerre variable z=R*(xi-1).
    # Direct overlap and derivative of the right orbital in real position.
    z, wz = np.polynomial.laguerre.laggauss(16)
    eta, we = np.polynomial.legendre.leggauss(12)
    xi = 1 + z[:, None] / R
    eta = eta[None, :]
    jac = xi**2 - eta**2
    weights = wz[:, None] * we[None, :]
    factor = R**2 * math.exp(-R) / 4
    overlap_quad = factor * float(np.sum(weights * jac))
    derivative_quad = factor * float(np.sum(
        weights * jac * (1-xi*eta)/(xi-eta)))
    cross_x_quad = factor * float(np.sum(weights * jac * R*xi*eta/2))
    S = math.exp(-R) * (1+R+R**2/3)
    Sp = -math.exp(-R) * (R+R**2)/3
    assert abs(overlap_quad-S) < 3e-14
    assert abs(derivative_quad+Sp) < 3e-14
    assert abs(cross_x_quad) < 3e-14
    G = np.array([[1., S], [S, 1.]])
    vals, vecs = np.linalg.eigh(G)
    lowdin = (vecs / np.sqrt(vals)) @ vecs.T
    x = lowdin @ np.diag([-R/2, R/2]) @ lowdin
    p_raw = np.array([[0, -1j*derivative_quad],
                      [1j*derivative_quad, 0]])
    p = lowdin @ p_raw @ lowdin
    ell = R/math.sqrt(1-S*S)
    assert np.linalg.norm(x-np.diag([-ell/2, ell/2])) < 5e-14
    p_exact = np.array([[0, 1j*Sp], [-1j*Sp, 0]])/math.sqrt(1-S*S)
    assert np.linalg.norm(p-p_exact) < 5e-14
    x6 = lift_one_body(np.kron(x, np.eye(2)))
    j6 = lift_one_body(np.kron(p, np.eye(2)))
    Q = lift_one_body(np.diag([-.5, -.5, .5, .5]))
    D = Q @ Q
    assert np.linalg.norm(x6-ell*Q) < 8e-14
    assert np.linalg.norm(D@D-D) < 1e-14
    assert np.trace(D) == 2
    assert abs(np.linalg.norm(j6, 2)-2*abs(Sp)/math.sqrt(1-S*S)) < 8e-14
    controlled_phase = np.diag(np.exp(-1j*math.pi*np.diag(Q)))
    K0 = (np.eye(6)+controlled_phase)/2
    K1 = (np.eye(6)-controlled_phase)/2
    assert np.linalg.norm(K0-(np.eye(6)-D)) < 1e-14
    assert np.linalg.norm(K1-D) < 1e-14
    return {"R": R, "overlap": S,
            "direct_quadrature_overlap_error": abs(overlap_quad-S),
            "direct_momentum_derivative_error": abs(derivative_quad+Sp),
            "ell": ell, "physical_j_compressed_norm": float(np.linalg.norm(j6,2)),
            "full_six_dimensional_Luders_identity_passed": True}


def run():
    R = F(1000)
    # e^R > R^3/3!, without an exponential float underflow.
    S_upper = 6*(1+R+R**2/3)/R**3
    assert S_upper < F(1,100)
    s = F(1,100)
    b = F(20,7)
    assert b*b > 8/(1-s)
    # H bound includes the positive nuclear constant 1/R.
    assert (10-1/R)**2 > 98/(1-s)
    E = F(42)
    J = F(37,2)
    assert J*J > 8*E
    T = R/F(137)
    L_upper = F(22,7)/R
    eps = L_upper*(b+T*J/2)
    assert eps == F(74239,335650)
    assert eps < F(9,40)
    assert 1-2*eps > F(11,20)
    # Independent finite integration of the chosen sin² envelope.
    u, wu = np.polynomial.legendre.leggauss(32)
    times = (u+1)*float(T)/2
    beta = 2*math.pi/(float(R*T))*np.sin(math.pi*times/float(T))**2
    area = float(wu@beta)*float(T)/2
    first_moment = float(wu@(times*beta))*float(T)/2
    assert abs(area-math.pi/float(R)) < 1e-15
    assert abs(first_moment-float(T)*area/2) < 1e-14
    return {
        "round": 1050, "review_only": True, "all_passed": True,
        "does_not_import_author_science": True,
        "does_not_simulate_full_Coulomb_propagation": True,
        "rational_certificate": {
            "overlap_upper_without_underflow": str(S_upper),
            "b_upper": str(b), "H_norm_with_nuclear_constant_upper": "10",
            "E_upper": str(E), "J_upper": str(J), "T": str(T),
            "epsilon_upper": str(eps), "epsilon_decimal": float(eps),
            "readout_probability_contrast_lower": str(1-2*eps),
            "current_first_moment_error_upper": str(eps*(2*J+4*L_upper))},
        "symmetric_pulse": {"area": area, "first_moment": first_moment},
        "direct_Slater_CAR_checks": [slater_review(R) for R in (.7, 2., 5.)],
        "scope": "Initial charge measurement followed by exact natural W_T; NR controlled dipole rights only",
    }


def compare(left, right):
    if isinstance(left, dict):
        assert left.keys() == right.keys()
        for key in left:
            compare(left[key], right[key])
    elif isinstance(left, list):
        assert len(left) == len(right)
        for a, b in zip(left, right):
            compare(a, b)
    elif isinstance(left, float):
        assert math.isclose(left, right, rel_tol=1e-11, abs_tol=1e-13)
    else:
        assert left == right


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with RESULT.open("x", encoding="utf-8") as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write("\n")
    else:
        compare(result, json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"round": 1050, "independent_review_checks_passed": True,
                      "epsilon_upper": result["rational_certificate"]["epsilon_upper"],
                      "mode": "exclusive_write" if args.write else "read_only"}))


if __name__ == "__main__":
    main()
