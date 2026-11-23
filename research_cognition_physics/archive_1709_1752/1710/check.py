"""Reconstruct the CP basis, then test which gauge-kinetic jets fix the strong phase."""
import argparse
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


CFG = read("inputs.json")
FLAVOR = read(CFG["color_sector"])


def cp_matching():
    frozen = read(CFG["joint_seed"])["predeclared_endpoints"][CFG["candidate"]]
    model = runpy.run_path(str(HERE / CFG["charged_graph"]))
    model["full"]["cfg"]["direct_semileptonic_moduli"]["Vub"][0] = frozen["Vub"]
    root = model["sim"]["graph"](frozen["k"], frozen["yd_widths"], frozen["port"], .2)[1]
    np.testing.assert_allclose(root["gamma_MZ_deg"], frozen["quark"]["gamma_deg"], rtol=0, atol=1e-9)
    block = root["down_texture_10PeV"]
    reduced = np.array(block["real"]) + 1j * np.array(block["imag"])
    phi = math.pi * FLAVOR["scalar_angle_over_pi"]
    omega = np.exp(4j * phi)
    long = reduced[1, 2].imag / omega.imag
    short = reduced[1, 2].real - long * omega.real
    amplitudes = [reduced[0, 1].real, reduced[1, 0].real, reduced[1, 1].real,
                  reduced[2, 1].real, reduced[2, 2].real, short, long]
    powers = np.array(read(CFG["down_powers"])["powers_A_B_C_E_F_R_D"])
    a, b, c, e, f, u, v = amplitudes * np.exp(1j * powers * phi)
    down = np.array([[0, a, 0], [b, c, u + v], [0, e, f]])
    up = np.diag(root["Y_up_at_10PeV"]) @ np.diag(np.exp(1j * phi * np.array([0, 0, 2])))
    left = np.diag(np.exp(1j * phi * np.array([0, 0, 1])))
    right = np.diag(np.exp(1j * phi * np.array([2, 0, 1])))
    np.testing.assert_allclose(down, left @ reduced @ right, rtol=2e-15, atol=1e-18)
    y, y_ref = np.linalg.svd(down, compute_uv=False), np.linalg.svd(reduced, compute_uv=False)
    np.testing.assert_allclose(y, y_ref, rtol=2e-13)
    def weak_cp(matrix):
        h = matrix @ matrix.conj().T
        return (h[0, 1] * h[1, 2] * h[2, 0]).imag
    np.testing.assert_allclose(weak_cp(down), weak_cp(reduced), rtol=2e-13)
    assert abs(weak_cp(down)) > 1e-21
    s = FLAVOR["scalar_magnitude_over_cutoff"] * np.exp(1j * phi)
    color = FLAVOR["cutoff_TeV"] * np.array([[0, s**3], [-s**3, 1.]])
    colored_masses = np.linalg.svd(color, compute_uv=False)
    np.testing.assert_allclose(colored_masses, read(CFG["heavy_validation"])["conditional_masses_TeV_descending"])
    light = np.linalg.det(up) * np.linalg.det(down)
    heavy = np.linalg.det(color)
    total = light * heavy
    assert abs(np.angle(total)) < 1e-14 and total.real > 0
    return {"gamma_MZ_deg": root["gamma_MZ_deg"], "cp_basis_signed_weak_edges": amplitudes,
            "upper_weak_cp_triple": weak_cp(down), "light_determinant_phase_mod2pi": np.angle(light) % (2*math.pi),
            "heavy_determinant_phase_mod2pi": np.angle(heavy) % (2*math.pi),
            "full_determinant_phase_rad": np.angle(total), "color_Dirac_masses_TeV": colored_masses.tolist()}


def first_operators():
    fields = {"S": 1, "X": 0, "A2": 2, "A3": 3}
    q = FLAVOR["flavor_charge"]
    out = {}
    for degree in range(1, 4):
        rows = []
        for term in combinations_with_replacement(fields, degree):
            if sum(q[f] for f in term) % 12 == 0 and sum(fields[f] for f in term) % 6:
                rows.append("*".join(term))
        out[str(degree)] = rows
    assert out == {"1": [], "2": [], "3": ["S*X*A3", "X*A2*A2"]}
    return out


def bounds(radius, x, jet_order, lift_all_charges):
    # z = I in the one-generator class and z = S^4 on the full frozen-X valley.
    z = radius * np.exp(-2j * math.pi / 3)
    b = abs(z.imag)
    factor = Polynomial([radius**2, radius, 1])
    p = (factor**jet_order).integ()
    p.coef[0] -= p(z).real
    # In the full charge ring each z^j is lifted as z^j X^(j mod 3).
    lifted = [c / x**(j % 3) if lift_all_charges else c for j, c in enumerate(p.coef)]
    budget = sum(abs(c) for c in lifted)
    normalized = p / budget
    recovered = sum(c * z**i * (x**(i % 3) if lift_all_charges else 1.) for i, c in enumerate(lifted))
    np.testing.assert_allclose(recovered, p(z), atol=1e-18)
    imag_shift = abs(normalized(z).imag)
    assert abs(normalized(z).real) < 1e-16
    for order in range(1, jet_order+1):
        assert abs(normalized.deriv(order)(z)) < 1e-14
    assert abs(normalized.deriv(jet_order+1)(z)) > 0
    integral = math.prod(2*j/(2*j+1) for j in range(1, jet_order+1))
    np.testing.assert_allclose(imag_shift, b**(2*jet_order+1) * integral / budget, rtol=1e-13)
    rho = CFG["cauchy_radius"]
    imax = math.ceil(rho/(1-rho))
    max_derivative = imax * rho**(imax-1)
    cap = max_derivative * b * integral * (b/(rho-radius))**(2*jet_order)
    assert cap >= imag_shift
    # Integral of endpoint Hermite interpolation on the vertical conjugacy segment.
    n = jet_order
    weights = [math.factorial(2*n-1-j)*math.factorial(n)
               / (math.factorial(2*n)*math.factorial(j+1)*math.factorial(n-1-j))
               for j in range(n)]
    for degree in range(2*n):
        def endpoint_derivative(endpoint, j):
            if j > degree:
                return 0.
            return math.factorial(degree)/math.factorial(degree-j)*endpoint**(degree-j)
        approx = sum(w * (endpoint_derivative(0, j) + (-1)**j*endpoint_derivative(1, j))
                     for j, w in enumerate(weights))
        assert abs(approx-1/(degree+1)) < 1e-14
    remainder = b * (2*b)**(2*n) / (math.comb(2*n, n)*(1-radius)**(2*n+1))
    error_gains = [2*b * w * (2*b)**j for j,w in enumerate(weights)]
    return {"local_jet_order": jet_order, "witness_shift_abs_rad": imag_shift,
            "bound_for_all_unit_budget_differences_rad": cap,
            "finite_jet_errors_gains_first_to_nth": error_gains,
            "finite_jet_unit_budget_remainder_rad": remainder,
            "finite_jet_upper_at_each_error_1e_minus9": remainder + 1e-9*sum(error_gains),
            "witness_first_distinguishing_derivative_modulus": abs(normalized.deriv(jet_order+1)(z)),
            "witness_real_coeffs_in_registered_monomial_order": [c / budget for c in lifted]}


def compute():
    eps, x = FLAVOR["scalar_magnitude_over_cutoff"], FLAVOR["real_X_over_cutoff"]
    limit = read(CFG["single_source_edm"])["conditional_edm"]["theta_limit_rad"]
    results = {}
    for whole, key in ((False, "single_I"), (True, "whole_ring_along_radial_valley")):
        radius = eps**4 if whole else x * eps**4
        results[key] = [bounds(radius, x, n, whole) for n in range(5)]
        order = 4 if whole else 3
        assert results[key][order-1]["witness_shift_abs_rad"] > limit
        assert 2 * results[key][order]["bound_for_all_unit_budget_differences_rad"] < limit
    leading = x * eps**4 * math.sqrt(3) / 2
    return {"actual_1708_matching": cp_matching(), "phase_nonreal_holomorphic_terms": first_operators(),
            "single_operator_theta_per_coefficient_rad": leading, "single_operator_bound": limit/leading,
            "single_source_EDM_working_bound_rad": limit, **results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = compute()
    destination = HERE / "results.json"
    if args.check:
        assert json.loads(destination.read_bytes()) == results
        print("1710 CP matching, charge invariants and finite-derivative bounds verified")
    else:
        destination.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
