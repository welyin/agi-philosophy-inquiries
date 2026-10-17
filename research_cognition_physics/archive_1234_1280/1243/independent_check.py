"""Independent 1243 audit: rational finite-sheet and weak-trace checks.

No author code is imported. Fractions certify the displayed arithmetic;
geometric identities and quantifiers are proved in independent_review.md.
The lapse example is an interface calibration, not an implemented instrument.
Default execution only compares the saved result. --write creates it once.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def poly_add(*items):
    out = [F(0)] * max(map(len, items), default=1)
    for p in items:
        for i, v in enumerate(p):
            out[i] += v
    return out


def poly_scale(p, c):
    return [c * x for x in p]


def poly_mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def derivative(p):
    return [F(i) * p[i] for i in range(1, len(p))] or [F(0)]


def value(p, x):
    return sum((v * x ** i for i, v in enumerate(p)), F(0))


def integral(p, a, b):
    return sum((v * (b ** (i + 1) - a ** (i + 1)) / (i + 1)
                for i, v in enumerate(p)), F(0))


def pack(q):
    q = F(q)
    return {"rational": str(q), "decimal_diagnostic": float(q)}


def certified_budget(ell):
    delta, b0, B = F(1, 10**6), F(1, 10**6), F(1, 100)
    Dsigma, Lr, Lf = F(1, 10**4), F(1, 1000), F(1, 1000)
    dr, df, cover = F(1, 10**5), F(1, 10**5), F(1, 10**5)
    Rstar = Fstar = F(1, 100)
    z = b0 * ell + B * ell**2 / 2
    assert 0 <= z < 1
    M = 1 / (1 - z)  # exp(z) <= sum z^n = 1/(1-z)
    N = b0**2 / 2 + b0*B*ell/3 + (Dsigma+B**2/2)*ell**2/6
    parts = {
        "fixed_readout_noise": 2*delta/ell**2,
        "terminal_expansion": 2*b0/ell,
        "transverse_transport": dr+df,
        "area_kernel_derivative": Lr*ell/3,
        "heat_kernel_derivative": 2*Lf*ell/3,
        "area_jacobian": (M-1)*(Rstar+Fstar),
        "shear_expansion": M*N,
        "finite_cover": cover,
    }
    zeta = sum(parts.values(), F(0))
    epsS, W = 9*zeta, F(53, 4)
    trace_bound = epsS*W + F(1, 10**5)
    return {
        "ell": pack(ell), "exp_argument": pack(z), "exp_upper": pack(M),
        "parts": {k: pack(v) for k, v in parts.items()},
        "zeta": pack(zeta), "epsilon_S": pack(epsS),
        "W1_upper": pack(W), "trace_difference_bound": pack(trace_bound),
        "full_equation_average_bound": pack(epsS+2*trace_bound),
        "gap_to_conditional_lower_005": pack(F(1, 20)-trace_bound),
    }


def finite_windows():
    ell = F(1, 10)
    area_weight, heat_weight = [ell, -F(1)], [F(0), F(1)]
    area_mass = integral(area_weight, F(0), ell)
    heat_mass = integral(heat_weight, F(0), ell)
    area_x = integral(poly_mul(area_weight, [F(0), F(1)]), F(0), ell)/area_mass
    area_x2 = integral(poly_mul(area_weight, [F(0), F(0), F(1)]), F(0), ell)/area_mass
    heat_x = integral(poly_mul(heat_weight, [F(0), F(1)]), F(0), ell)/heat_mass
    assert area_mass == heat_mass == ell**2/2
    assert (area_x, area_x2, heat_x) == (ell/3, ell**2/6, 2*ell/3)

    # Isotropic Rosen scale w(lambda); F=w^2. Raychaudhuri is checked
    # by polynomial identities, with a nonzero terminal theta retained.
    w = [F(1), F(1, 1000), -F(1, 100)]
    wp, wpp = derivative(w), derivative(derivative(w))
    area_J = poly_mul(w, w)
    F_r = poly_scale(poly_mul(w, wpp), -F(2))
    F_half_theta2 = poly_scale(poly_mul(wp, wp), F(2))
    d_Ftheta = derivative(poly_scale(poly_mul(w, wp), F(2)))
    ray_rhs = poly_add(F_half_theta2, poly_scale(F_r, -F(1)))
    assert poly_add(d_Ftheta, poly_scale(ray_rhs, -F(1))) == [F(0)]*len(d_Ftheta)
    k_area, k_heat = [ell, F(1)], [F(0), -F(1)]
    end = 2*ell*w[1]
    ricci_term = integral(poly_mul(k_area, F_r), -ell, F(0))
    expansion_term = -integral(poly_mul(k_area, F_half_theta2), -ell, F(0))
    delta_area = value(area_J, F(0))-value(area_J, -ell)
    assert delta_area == end+ricci_term+expansion_term
    wrong_kernel_ricci = integral(poly_mul(k_heat, F_r), -ell, F(0))
    assert wrong_kernel_ricci != ricci_term
    assert end != 0 and expansion_term != 0

    main = certified_budget(ell)
    small = certified_budget(ell/10)
    main_bound = F(main["trace_difference_bound"]["rational"])
    small_bound = F(small["trace_difference_bound"]["rational"])
    assert main_bound < F(42, 1000)
    assert F(1, 20)-main_bound > F(8, 1000)
    assert small_bound > main_bound
    assert F(321, 4) < 81  # inherited algebra, not a new reconstruction theorem
    return {
        "status": "PASS",
        "kernel_moments": {
            "equal_mass": pack(area_mass), "area_mean_x": pack(area_x),
            "area_mean_x2": pack(area_x2), "heat_mean_x": pack(heat_x),
        },
        "independent_Rosen_polynomial_identity": {
            "area_difference": pack(delta_area), "terminal": pack(end),
            "ricci_integral": pack(ricci_term), "expansion_integral": pack(expansion_term),
            "wrong_heat_kernel_ricci": pack(wrong_kernel_ricci),
            "kernel_swap_error": pack(wrong_kernel_ricci-ricci_term),
            "exact_residual": pack(delta_area-end-ricci_term-expansion_term),
        },
        "declared_budget": main, "ten_times_shorter_window": small,
        "scope": "Finite data and regularity inputs are conditional; no instrument realization is certified.",
    }


def weak_trace():
    h, d = F(1, 4), F(1, 2)
    q = [F(15, 16)/h, F(0), -F(30, 16)/h**3, F(0), F(15, 16)/h**5]
    qp = derivative(q)
    assert integral(q, -h, h) == 1
    assert value(q, h) == value(q, -h) == 0
    assert value(qp, h) == value(qp, -h) == 0
    q_var = integral(qp, -h, F(0))-integral(qp, F(0), h)
    assert q_var == F(15, 8)/h
    # CDF difference U = integral_0^d q(t-s) ds times transverse Q.
    # Fubini yields exact mass and first moment without quadrature.
    U_mass, U_first = d, d**2/2
    derivative_L1_upper = F(2)+3*d*q_var
    assert derivative_L1_upper == F(53, 4)

    # Independent nonconstant-volume calibration: g=(-N(t)^2,1,1,1),
    # N=1+nu*t, rho=N, V^0=U/N, S00=3N^3, Sii=N, Lambda=-3N.
    # T00=0, Tii=-4N/kappa is conserved. The pure lapse metric is flat,
    # but its Christoffel term is nonzero and essential to the pairing.
    nu = F(1, 5)
    lapse_max = 1+nu*(d+h)
    lhs = 3*nu*d
    density_derivative_pairing = 6*nu*U_mass
    christoffel_pairing = -3*nu*U_mass
    rhs = density_derivative_pairing+christoffel_pairing
    assert rhs == lhs == F(3, 10)
    wrong_no_connection = density_derivative_pairing
    assert wrong_no_connection != lhs
    lapse_W_upper = 2+lapse_max*3*d*q_var
    assert lapse_W_upper == F(239, 16)
    # |S|_e=sqrt(12)N <4N, with the same positive auxiliary norm.
    assert abs(lhs) < 4*lapse_max*lapse_W_upper

    # Reconstruct the coordinate Christoffel symbols and covariant
    # derivatives independently, retaining density in the final pairing.
    tensor_rows = []
    for t in (-h, F(0), d, d+h):
        N = 1+nu*t
        metric, inverse = [-N*N, F(1), F(1), F(1)], [-1/N**2, F(1), F(1), F(1)]
        dg = [[[F(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
        dg[0][0][0] = -2*N*nu
        Gamma = [[[inverse[a]*(dg[b][a][c]+dg[c][a][b]-dg[a][b][c])/2
                   for c in range(4)] for b in range(4)] for a in range(4)]
        S = [[F(0) for _ in range(4)] for _ in range(4)]
        T = [[F(0) for _ in range(4)] for _ in range(4)]
        dS = [[[F(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
        dT = [[[F(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
        S[0][0], dS[0][0][0] = 3*N**3, 9*N*N*nu
        for i in (1, 2, 3):
            S[i][i], dS[0][i][i] = N, nu
            T[i][i], dT[0][i][i] = -4*N, -4*nu

        def divergence(tensor, dtensor):
            return [sum((inverse[a]*(dtensor[a][a][b]
                        -sum((Gamma[c][a][a]*tensor[c][b]
                              +Gamma[c][a][b]*tensor[a][c] for c in range(4)), F(0)))
                        for a in range(4)), F(0)) for b in range(4)]

        divS, divT = divergence(S, dS), divergence(T, dT)
        assert sum((inverse[a]*S[a][a] for a in range(4)), F(0)) == 0
        assert divS == [-3*nu, 0, 0, 0] and divT == [0]*4
        coefficients = []
        for U, Ut in ((F(1), F(0)), (F(0), F(1))):
            V = [U/N, F(0), F(0), F(0)]
            dV = [[F(0) for _ in range(4)] for _ in range(4)]
            dV[0][0] = Ut/N-U*nu/N**2
            nabla = [[inverse[a]*(dV[a][b]
                       +sum((Gamma[b][a][c]*V[c] for c in range(4)), F(0)))
                      for b in range(4)] for a in range(4)]
            pairing_density = N*sum((S[a][b]*nabla[a][b]
                                     for a in range(4) for b in range(4)), F(0))
            coefficients.append(pairing_density)
        assert coefficients == [F(0), -3*N]
        tensor_rows.append({"t": pack(t), "lapse": pack(N),
                            "divS": [pack(v) for v in divS],
                            "divT": [pack(v) for v in divT],
                            "density_pairing_U_and_Uprime_coefficients": [pack(v) for v in coefficients]})

    # Conservation deletion: S=0, Lambda=t, T=Lambda*g/kappa.
    # kappa <j,V> = integral U = d, so Lambda_a-Lambda_b=-d.
    deletion_lhs, kappa_j_pair = -d, U_mass
    assert deletion_lhs == -kappa_j_pair and deletion_lhs != 0

    # Conservation alone allows high-frequency trace gradients. For
    # S00=3f,Sii=f,Lambda=-3f, E00=0 and spatial Eii=4f;
    # f=alpha sin(Nt), so j=0, sup |S|_e^2=12 alpha^2.
    alpha = F(1, 1000)
    gradient_values = [-3*alpha*n for n in (100, 10000)]
    assert gradient_values[1]/gradient_values[0] == 100
    return {
        "status": "PASS",
        "compact_window": {
            "normalization": pack(integral(q, -h, h)), "q_total_variation": pack(q_var),
            "U_L1": pack(U_mass), "U_first_moment": pack(U_first),
            "flat_derivative_sum_upper": pack(derivative_L1_upper),
        },
        "nonconstant_volume_calibration": {
            "lapse_slope": pack(nu), "lapse_max": pack(lapse_max),
            "Lambda_window_difference": pack(lhs),
            "density_partial_derivative_contribution": pack(density_derivative_pairing),
            "connection_contribution": pack(christoffel_pairing),
            "full_weak_pairing": pack(rhs), "omitted_connection_wrong_pairing": pack(wrong_no_connection),
            "W1_auxiliary_norm_upper": pack(lapse_W_upper), "exact_identity_residual": pack(lhs-rhs),
            "coordinate_tensor_diagnostics": tensor_rows,
        },
        "conservation_deletion": {
            "all_null_residuals": pack(0), "Lambda_window_difference": pack(deletion_lhs),
            "kappa_source_divergence_pairing": pack(kappa_j_pair),
        },
        "no_point_gradient_bound": {
            "fixed_S_sup_norm_squared": pack(12*alpha**2),
            "Lambda_derivatives_at_zero": [pack(v) for v in gradient_values],
            "divergence_defect": pack(0),
        },
        "scope": "Tensor interface calibrations, not physical material realizations; no claim of pointwise Lambda constancy.",
    }


def compute():
    names = ["admission.md", "proof.md", "check.py", "results.json", "NEXT.md"]
    return {
        "round": 1243, "audit": "independent", "scientific_groups": 2,
        "status": "PASS", "author_sha256": {
            name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names
        },
        "independent_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "groups": {"finite_windows": finite_windows(), "weak_trace": weak_trace()},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = HERE / "independent_results.json"
    result = compute()
    if args.write:
        with output.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("CREATED_ONCE PASS independent scientific groups=2")
    else:
        saved = json.loads(output.read_text(encoding="utf-8"))
        if saved != result:
            raise AssertionError("Saved independent results or bound author files differ")
        print("READ_ONLY_MATCH PASS independent scientific groups=2")
    print("finite_trace_upper", result["groups"]["finite_windows"]["declared_budget"]["trace_difference_bound"]["rational"])


if __name__ == "__main__":
    main()
