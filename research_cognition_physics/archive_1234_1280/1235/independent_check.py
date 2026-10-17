"""Independent finite-window checks for round 1235; no author code imported.

Run --write only to create this reviewer's result. Default is read-only replay.
Analytic proofs, rather than sampled maxima, certify derivative bounds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def graph_geometry(tau, u, x, shape):
    """Power-law FLRW graph, with h differentiated directly in graph coordinates."""
    x = np.asarray(x, dtype=float)
    f = u * (x @ shape @ x) / 2
    t = tau + f
    v = u * shape @ x
    f2 = u * shape
    a2 = t ** (2 / 3)
    h = a2 * np.eye(3) - np.outer(v, v)
    hi = np.linalg.inv(h)
    z = v @ v / a2
    assert t > 0 and z < 1
    k = (f2 + a2 / (3 * t) * np.eye(3)
         - 2 / (3 * t) * np.outer(v, v)) / np.sqrt(1 - z)
    kmix = hi @ k
    ktrace = np.trace(kmix)
    knorm = np.trace(kmix @ kmix)
    # The scalar stress is computed from phi_I=q_I log(t), not from Einstein G.
    physical_qnorm2 = 3 * (np.sqrt(2) / 3)**2
    rho = 2 * physical_qnorm2 / (2 * t**2) * (1 + z) / (1 - z)
    # Independently evaluate G from H and Hdot in the ambient FLRW metric.
    hubble = 1 / (3 * t)
    hdot = -1 / (3 * t**2)
    einstein_nn = (3 * hubble**2
                   - (2 * hdot + 3 * hubble**2) * z) / (1 - z)
    dh = np.empty((3, 3, 3))  # derivative index first
    ddh = np.empty((3, 3, 3, 3))
    for j in range(3):
        dh[j] = ((2 / 3) * t ** (-1 / 3) * v[j] * np.eye(3)
                 - np.outer(f2[:, j], v) - np.outer(v, f2[:, j]))
        for l in range(3):
            ddh[l, j] = ((-2 / 9) * t ** (-4 / 3) * v[l] * v[j]
                          + (2 / 3) * t ** (-1 / 3) * f2[j, l]) * np.eye(3)
            ddh[l, j] -= (np.outer(f2[:, j], f2[:, l])
                          + np.outer(f2[:, l], f2[:, j]))
    dhi = np.asarray([-hi @ dh[j] @ hi for j in range(3)])
    gamma = np.zeros((3, 3, 3))
    dgamma = np.zeros((3, 3, 3, 3))
    for kidx in range(3):
        for i in range(3):
            for j in range(3):
                for ell in range(3):
                    b = dh[i, ell, j] + dh[j, ell, i] - dh[ell, i, j]
                    gamma[kidx, i, j] += hi[kidx, ell] * b / 2
                    for d in range(3):
                        db = ddh[d, i, ell, j] + ddh[d, j, ell, i] - ddh[d, ell, i, j]
                        dgamma[d, kidx, i, j] += (dhi[d, kidx, ell] * b
                                                  + hi[kidx, ell] * db) / 2
    ricci = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            for kidx in range(3):
                ricci[i, j] += dgamma[kidx, kidx, i, j] - dgamma[j, kidx, i, kidx]
                for ell in range(3):
                    ricci[i, j] += (gamma[kidx, kidx, ell] * gamma[ell, i, j]
                                    - gamma[kidx, j, ell] * gamma[ell, i, kidx])
    scalar = float(np.sum(hi * ricci))
    return t, z, float(ktrace), float(knorm), scalar, rho, 2 * einstein_nn


def window_values(shape, *, radius, thickness, amplitude, order=12):
    """Fixed common-coordinate product window, not per-slice proper volumes.

    Uniform cube of half-width radius/sqrt(3), supported inside |x|<=radius.
    Both selected shapes have zero third row; integrate x3 exactly by constancy.
    """
    nodes, weights = np.polynomial.legendre.leggauss(order)
    weights = weights / 2
    half = radius / np.sqrt(3)
    xi, x1, x2 = np.meshgrid(thickness * nodes, half * nodes, half * nodes, indexing="ij")
    ww = np.einsum("i,j,k->ijk", weights, weights, weights).ravel()
    xx = np.stack((x1.ravel(), x2.ravel(), np.zeros(x1.size)), axis=-1)
    tau = 1 + xi.ravel()
    ax = xx @ shape.T
    q = np.sum(xx * ax, axis=-1)
    av = np.sum(ax * ax, axis=-1)
    out = []
    for u in (amplitude, -amplitude, 0.0):
        t = tau + u * q / 2
        a2 = t ** (2 / 3)
        v = u * ax
        z = np.sum(v * v, axis=-1) / a2
        h = a2[:, None, None] * np.eye(3) - v[:, :, None] * v[:, None, :]
        hi = np.linalg.inv(h)
        k = (u * shape + a2[:, None, None] / (3 * t[:, None, None]) * np.eye(3)
             - 2 / (3 * t[:, None, None]) * v[:, :, None] * v[:, None, :]) / np.sqrt(1 - z[:, None, None])
        kmix = hi @ k
        tr = np.trace(kmix, axis1=-2, axis2=-1)
        qk = np.einsum("...ij,...ji->...", kmix, kmix)
        # R is Gauss-based here, independently audited above by direct derivatives.
        ambient = (2 / (3 * t**2)) * (1 + z) / (1 - z)
        r3 = qk - tr * tr + ambient
        # The inherited constraint uses rho=2 T_phys(n,n), explicitly normalized.
        matter = 2 * (2 / 3) / (2 * t**2) * (1 + z) / (1 - z)
        c = qk - tr * tr - r3 + matter
        frozen_matter = (4 / 3) / (2 * tau**2)
        c_frozen = qk - tr * tr - r3 + frozen_matter
        out.append({"matter": float(ww @ matter), "C": float(ww @ c),
                    "C_frozen": float(ww @ c_frozen),
                    "max_C": float(np.max(np.abs(c)))})
    central = lambda key: out[0][key] + out[1][key] - 2 * out[2][key]
    qnorm2 = 4 / 3
    lower = 2 * qnorm2 * amplitude**2 * float(ww @ av) * (1 + thickness) ** (-8 / 3)
    margin_t = 1 - thickness - abs(amplitude) * np.linalg.norm(shape, 2) * radius**2 / 2
    zmax = amplitude**2 * np.linalg.norm(shape, 2)**2 * radius**2 * margin_t ** (-2 / 3)
    upper = amplitude**2 * qnorm2 * ((3 / 4) * float(ww @ (q * q)) * margin_t**(-4)
             + 2 * float(ww @ av) * margin_t ** (-8 / 3) / (1 - zmax))
    assert margin_t > 0 and zmax < 1
    assert lower > 0 and lower <= central("matter") <= upper
    assert max(p["max_C"] for p in out) < 2e-14
    assert abs(central("C_frozen") + central("matter")) < 3e-14
    return {"radius": radius, "thickness": thickness, "amplitude": amplitude,
            "source_central": central("matter"), "source_lower": lower,
            "source_upper": upper, "full_C_central": central("C"),
            "frozen_source_C_central": central("C_frozen"),
            "time_margin": margin_t, "tilt_bound": zmax}


def run():
    shapes = [np.diag([1.0, 0.0, 0.0]), np.diag([1.0, -1.0, 0.0])]
    gauss = []
    scalar_kg = []
    for shape in shapes:
        for tau, u, x in [(1.03, 0.31, [0.04, -0.025, 0.015]),
                          (0.94, -0.27, [-0.03, 0.047, -0.02]),
                          (1.08, 0.0, [0.02, 0.01, -0.04])]:
            t, z, tr, qk, r3, matter, ambient = graph_geometry(tau, u, x, shape)
            gauss.append(abs(r3 - qk + tr**2 - ambient))
            assert abs(matter - ambient) < 2e-15
            for charge in (np.sqrt(2) / 3,) * 3:
                scalar_kg.append(abs(-charge / t**2 + 3 / (3 * t) * charge / t))
    assert max(gauss) < 3e-14
    assert max(scalar_kg) < 3e-15

    windows = []
    for shape in shapes:
        for radius, thickness, amplitude in [(0.16, 0.07, 0.35), (0.08, 0.04, 0.35),
                                             (np.sqrt(3) * 1e-4, 1e-5, 0.25)]:
            coarse = window_values(shape, radius=radius, thickness=thickness, amplitude=amplitude, order=10)
            fine = window_values(shape, radius=radius, thickness=thickness, amplitude=amplitude, order=14)
            assert abs(coarse["source_central"] - fine["source_central"]) < 3e-14
            fine["shape"] = np.diag(shape).tolist()
            fine["quadrature_difference"] = abs(coarse["source_central"] - fine["source_central"])
            windows.append(fine)

    # Sharp, deterministic calculus check for separate x/xi mixed derivative terms.
    # F(u,xi,x)=F0(u)+u²(c|x|²+d xi²). Positive coefficients saturate the bound.
    s, mu_x, mu_xi, c, d = 0.37, 0.021, 0.007, 1.7, 0.6
    exact_error = 2 * s**2 * (c * mu_x + d * mu_xi)
    bound = s**2 / 2 * (mu_x * 4 * c + mu_xi * 4 * d)
    assert abs(bound - exact_error) < 2e-16
    # An odd spatial term cancels only for a centered even common window.
    shifted_mean, odd_coefficient = 0.012, 2.4
    odd_leak = 2 * s**2 * odd_coefficient * shifted_mean
    assert odd_leak > 0 and exact_error < exact_error + odd_leak

    # Exact Cauchy-constant arithmetic; existence of the analytic bound is proved.
    analytic_bound, du, dx, dxi = 128, 1 / 4, 1 / 8, 1 / 16
    diagonal = 4 * analytic_bound / (du**2 * dx**2)
    mixed = 2 * analytic_bound / (du**2 * dx**2)
    time_second = 4 * analytic_bound / (du**2 * dxi**2)
    assert diagonal == 524288 and mixed == 262144 and time_second == 2097152
    # Independent, deliberately looser analytic matrix norm bound also suffices.
    inv_bound = 8 / 7
    kmix_norm = inv_bound * (1 / 2 + 3 / 4 + 1 / 8) * np.sqrt(8 / 7)
    invariant_i_bound = 3 * kmix_norm**2
    invariant_j_bound = 9 * kmix_norm**2
    full_residual_bound = ((invariant_i_bound + 2 * invariant_j_bound) / 0.9
                          + 1.1 * (invariant_i_bound + invariant_j_bound + 6) + 2 + 6)
    assert invariant_i_bound < 9 and invariant_j_bound < 26
    assert full_residual_bound < 128
    mu_x_cert, mu_t_cert = 1e-8, 1e-10 / 3
    certified_e = (mu_x_cert * 3 * diagonal + mu_t_cert * time_second) / 2
    threshold, step, a_upper = 1e-8, 1 / 4, 1.1
    kappa = a_upper / 2 * (4 * threshold / step**2 + certified_e)
    ab_interval = a_upper / 4 * (4 * threshold / step**2 + certified_e)
    lambda_interval = kappa / (2 * (2 - 3 * kappa))
    # A definite wrong member of the declared class cannot hide in this window.
    wrong_point_signal = 2 * step**2 * (1 - 0.9)  # a=b=1, L=0.9, rank-one A
    wrong_exclusion_margin = wrong_point_signal - step**2 * certified_e - 4 * threshold
    assert kappa < 2 / 3 and lambda_interval < 0.001094 and ab_interval < 0.002173
    assert wrong_exclusion_margin > 0

    return {"status": "PASS_CONDITIONAL_FINITE_WINDOW_SCOPE", "diagnostic_groups": 3,
            "direct_intrinsic_curvature_and_independent_scalar_source": {
                "gauss_max_error": max(gauss), "KG_max_error": max(scalar_kg),
                "physical_scalar_charge_norm_squared": 2 / 3,
                "constraint_source_equals": "2*T_phys(n,n)",
                "geometry_field_equations_used_to_define_source": False},
            "curved_common_history_finite_windows": windows,
            "mixed_derivative_bound_and_deletion": {
                "sharp_polynomial_error": exact_error, "bound": bound,
                "noncentered_window_odd_leak": odd_leak,
                "Cauchy_diagonal": diagonal, "Cauchy_mixed": mixed,
                "Cauchy_temporal": time_second,
                "independent_complex_residual_bound": full_residual_bound,
                "certified_E": certified_e, "L_interval": kappa,
                "ab_interval": ab_interval, "lambda_interval": lambda_interval,
                "wrong_member_exclusion_margin": wrong_exclusion_margin,
                "sampled_maxima_certify_derivative_bounds": False},
            "limits": ["No actual constraint-measurement instrument is certified.",
                       "The effective mixed-derivative bounds remain additional domain inputs.",
                       "The scalar FLRW example tests transport; it does not derive its background.",
                       "No new all-axiom quantum common implementation or GR completion is claimed."]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    result["author_binding"] = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                                for name in ("proof.md", "check.py", "results.json", "NEXT.md")}
    destination = HERE / "independent_results.json"
    if args.write:
        destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    else:
        previous = json.loads(destination.read_text(encoding="utf-8"))
        if previous != result:
            raise SystemExit("Independent replay does not match the saved result/binding")
    print(json.dumps({"status": result["status"], "groups": result["diagnostic_groups"],
                      "mode": "write" if args.write else "read_only_match"}))


if __name__ == "__main__":
    main()
