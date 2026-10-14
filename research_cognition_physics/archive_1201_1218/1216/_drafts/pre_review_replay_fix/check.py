"""1216: fixed two-time non-CP certificate. Default is read-only replay.

All strict decisions use integer/Fraction certificates. NumPy only proposes
trial inverses and a Choi witness. No temporal scan or clock simulation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
D = 2**80
PROPOSAL_D = 2**40
VECTOR_D = 2**30
N = 16


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def zeros(n):
    return np.zeros((n, n), dtype=object)


def eye_int(n, scale=1):
    return np.eye(n, dtype=object)*scale


def round_div(x, denominator):
    return (int(x) + denominator//2)//denominator


def round_matrix(a, denominator):
    return np.array([round_div(x, denominator) for x in a.flat], dtype=object).reshape(a.shape)


def fixed_fraction_matrix(a, scale=D):
    return np.array([round_div(F(x).numerator*scale, F(x).denominator)
                     for x in a.flat], dtype=object).reshape(a.shape)


def raw_product(a, b):
    ar, ai = a
    br, bi = b
    return ar@br-ai@bi, ar@bi+ai@br


def product(a, b, extra_divisor=1):
    r, i = raw_product(a, b)
    return round_matrix(r, D*extra_divisor), round_matrix(i, D*extra_divisor)


def dagger(a):
    return a[0].T, -a[1].T


def to_float(a, scale=D):
    return np.array(a[0], dtype=float)/scale + 1j*np.array(a[1], dtype=float)/scale


def frobenius_squared(a, denominator):
    return F(sum(int(x)**2 for x in a[0].flat)+sum(int(x)**2 for x in a[1].flat),
             denominator**2)


def proposal(a, scale):
    return (np.array([[int(round(float(x.real)*scale)) for x in row] for row in a], dtype=object),
            np.array([[int(round(float(x.imag)*scale)) for x in row] for row in a], dtype=object))


def encode(a):
    return {"real": [[int(x) for x in row] for row in a[0]],
            "imag": [[int(x) for x in row] for row in a[1]]}


def certified_exponential(h_r, t):
    # A is exactly anti-Hermitian. Scaling is fixed at 12 for both times.
    real = zeros(16)
    imag = fixed_fraction_matrix(-F(t, 4096)*h_r)
    require(np.array_equal(imag, imag.T), "anti-Hermitian fixed input")
    a_norm_bound = F(max(sum(abs(int(x)) for x in row) for row in imag), D)
    require(a_norm_bound <= F(1, 2), "scaled input norm")
    a = real, imag
    term = eye_int(16, D), zeros(16)
    result = term[0].copy(), term[1].copy()
    degree = 30
    for k in range(1, degree+1):
        term = product(term, a, k)
        result = result[0]+term[0], result[1]+term[1]
    tail = F(2, 2**31*math.factorial(31))
    error = F(2*16*degree, D) + tail
    bounds = [error]
    for _ in range(12):
        require(error < F(1, 1000), "squaring error hypothesis")
        result = product(result, result)
        error = F(2001, 1000)*error + F(16, D)
        bounds.append(error)
    # h_r differs only in symmetric sqrt(3)/80 entries; H_AB has two copies.
    original_h_error = F(2, 80*D)
    phase_error = t*original_h_error + F(4096*16, D)
    total_error = error + phase_error
    require(total_error < F(1, 10**12), "exponential finite certificate")
    return result, {
        "t": t, "scaled_A_norm_bound": str(a_norm_bound),
        "Taylor_degree": degree, "squarings": 12,
        "Taylor_tail_bound": str(tail),
        "squaring_error_bounds": [str(x) for x in bounds],
        "Hamiltonian_and_input_rounding_phase_bound": str(phase_error),
        "total_operator_error_bound": str(total_error),
        "certified_operator_error_less_than": "1/1000000000000",
    }


def reduced_channel(u, b_fixed):
    cols_r, cols_i = [], []
    for j in range(4):
        for i in range(4):
            inp = zeros(16)
            inp[4*i:4*(i+1), 4*j:4*(j+1)] = b_fixed
            mid = product(u, (inp, zeros(16)))
            out = product(mid, dagger(u))
            red_r, red_i = zeros(4), zeros(4)
            for a in range(4):
                for c in range(4):
                    red_r[a, c] = sum(out[0][4*a+k, 4*c+k] for k in range(4))
                    red_i[a, c] = sum(out[1][4*a+k, 4*c+k] for k in range(4))
            cols_r.append(red_r.reshape(-1, order="F"))
            cols_i.append(red_i.reshape(-1, order="F"))
    return np.stack(cols_r, axis=1), np.stack(cols_i, axis=1)


def choi_numerator(v):
    # Returned numerator is NOT yet divided by d=4.
    out_r, out_i = zeros(16), zeros(16)
    for i in range(4):
        for j in range(4):
            for a in range(4):
                for b in range(4):
                    out_r[4*i+a, 4*j+b] = v[0][a+4*b, i+4*j]
                    out_i[4*i+a, 4*j+b] = v[1][a+4*b, i+4*j]
    return out_r, out_i


def compute():
    paths = [
        "archive_370_428/research_note_390.md",
        "archive_1201_/1202/results.json",
        "archive_1201_/1203/proof.md",
        "archive_1201_/1203/results.json",
        "archive_1201_/1204/proof.md",
        "archive_1201_/1214/proof.md",
    ]
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    hashes["archive_1201_/1216/check.py"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    cold = json.loads((ROOT/paths[1]).read_text(encoding="utf-8"))
    p = [F(x) for x in cold["one_marginal_rational"]]
    require(sum(p) == 1 and min(p) > 0, "cold source")
    prep = cold["preparation"]
    rootden = int(prep["root_denominator"])
    boxes = [(F(a, rootden), F(b, rootden)) for a, b in prep["root_boxes"]]
    mids = [(a+b)/2 for a, b in boxes]
    root_error = F(1, 10**18)
    minsep = min(boxes[j][0]-boxes[i][1] for i in range(4) for j in range(i+1, 4))
    require(minsep > F(38, 1000), "global minimum spectral separation")
    require(boxes[-1][1]-boxes[0][0] < F(1094, 1000), "spectral width")
    sqrt_num = math.isqrt(3*D*D)
    require(sqrt_num**2 <= 3*D*D < (sqrt_num+1)**2, "sqrt integer bracket")
    sqrt3 = F(sqrt_num, D)
    hr = np.array([[F(0), F(-1, 5), F(0), -sqrt3/80],
                   [F(-1, 5), F(1), F(0), F(0)],
                   [F(0), F(0), F(1), F(0)],
                   [-sqrt3/80, F(0), F(0), F(-1, 40)]], dtype=object)
    ident4 = np.eye(4, dtype=object)
    qmat = zeros(4)
    qmat[1, 2] = qmat[2, 1] = 1
    H = np.kron(hr, ident4)+np.kron(ident4, hr)+F(1, 10)*np.kron(qmat, qmat)
    require(max(sum(abs(x) for x in row) for row in H) <= F(5, 2), "H norm")
    groups = {}
    groups["fixed_original_process_and_history"] = {
        "pass": True, "v": 1, "physical_times": [301, 601], "QQ_durations": [300, 600],
        "copy_scalar_support": [603, 605], "B": "same certified cold r_C at t=0",
        "independent_initial_factors": "unknown A/R_A, cold B/factory R_B, normal C, blank K",
        "intermediate_environment_reset": False, "time_scan": False,
        "history_reused": "390 finite-generator and tan(gt) results; not a new Markov theorem",
    }

    bhat = zeros(4)
    for i in range(4):
        proj = ident4.copy()
        for j in range(4):
            if i != j:
                proj = proj @ ((hr-mids[j]*ident4)/(mids[i]-mids[j]))
        bhat += p[i]*proj
    b_fixed = fixed_fraction_matrix(bhat)
    gap = F(38, 1000)
    width = F(1094, 1000)
    h_error = F(1, 80*D)
    factor_norm_bound = (width+h_error+root_error)/(gap-2*root_error)
    factor_error = (h_error+root_error)/(gap-2*root_error) \
        + width*(2*root_error)/(gap*(gap-2*root_error))
    require(factor_norm_bound < 30 and factor_error < F(2, 10**15), "projector factor bounds")
    b_error = 4*3*30**2*F(2, 10**15) + F(16, D)
    require(b_error < F(3, 10**11), "reference trace error")
    groups["native_reference_polynomial_certificate"] = {
        "pass": True, "source_populations": [str(x) for x in p],
        "sqrt3_lower": str(sqrt3), "sqrt3_upper": str(F(sqrt_num+1, D)),
        "root_midpoints": [str(x) for x in mids],
        "factor_norm_bound": str(factor_norm_bound), "factor_error_bound": str(factor_error),
        "reference_trace_error_bound": str(b_error), "reference_fixed_denominator": D,
        "reference_fixed_real": [[int(x) for x in row] for row in b_fixed],
    }

    us, certs = certified_exponential(H, 300)
    ut, certt = certified_exponential(H, 600)
    groups["two_finite_integer_exponential_certificates"] = {
        "pass": True, "fixed_denominator": D, "certificates": [certs, certt],
        "floating_diagonalization_used_for_exponential": False,
    }
    shat, that = reduced_channel(us, b_fixed), reduced_channel(ut, b_fixed)
    alpha = F(1, 10**12)
    b_safe = F(3, 10**11)
    column_error = 2*alpha+alpha**2+b_safe*(1+alpha)**2+F(256, D)*(2+alpha)
    channel_error = 4*column_error
    q = F(1, 10**8)
    require(channel_error < q, "Liouville certificate")
    groups["reduced_channel_error_certificate"] = {
        "pass": True, "matrix_unit_output_trace_error_bound": str(column_error),
        "Liouville_operator_error_bound": str(channel_error),
        "q_used": str(q), "vec": "column-major: row=a+4*b, column=i+4*j",
        "exact_dyadic_scale": D,
    }

    sf, tf = to_float(shat), to_float(that)
    xf = np.linalg.inv(sf)
    yf = tf@xf
    x, y = proposal(xf, PROPOSAL_D), proposal(yf, PROPOSAL_D)
    xn2 = frobenius_squared(x, PROPOSAL_D)
    yn2 = frobenius_squared(y, PROPOSAL_D)
    xs = raw_product(x, shat)
    rx = xs[0]-eye_int(16, PROPOSAL_D*D), xs[1]
    sx = raw_product(shat, x)
    rx_right = sx[0]-eye_int(16, PROPOSAL_D*D), sx[1]
    ys = raw_product(y, shat)
    rv = ys[0]-PROPOSAL_D*that[0], ys[1]-PROPOSAL_D*that[1]
    rx2 = frobenius_squared(rx, PROPOSAL_D*D)
    rxright2 = frobenius_squared(rx_right, PROPOSAL_D*D)
    rv2 = frobenius_squared(rv, PROPOSAL_D*D)
    r = F(1, 10**8)
    require(xn2 < 36 and yn2 < 16, "proposal norm certificates")
    require(max(rx2, rxright2, rv2) < r*r, "exact residual certificates")
    groups["exact_inverse_and_intermediate_residuals"] = {
        "pass": True, "proposal_denominator": PROPOSAL_D,
        "X_Frobenius_squared": str(xn2), "Y_Frobenius_squared": str(yn2),
        "X_Shat_minus_I_Frobenius_squared": str(rx2),
        "Shat_X_minus_I_Frobenius_squared": str(rxright2),
        "Y_Shat_minus_That_Frobenius_squared": str(rv2),
        "r_used": str(r), "X_norm_upper": 6, "Y_norm_upper": 4,
    }

    jn = choi_numerator(y)
    jf = to_float(jn, 4*PROPOSAL_D)
    evals, evecs = np.linalg.eigh((jf+jf.conj().T)/2)
    zf = evecs[:, 0]
    # Fix an irrelevant global phase for reproducible proposed coefficients.
    pivot = int(np.argmax(np.abs(zf)))
    zf *= np.exp(-1j*np.angle(zf[pivot]))
    zr = np.array([int(round(float(a.real)*VECTOR_D)) for a in zf], dtype=object)
    zi = np.array([int(round(float(a.imag)*VECTOR_D)) for a in zf], dtype=object)
    z_norm_num = sum(int(a)**2 for a in zr)+sum(int(a)**2 for a in zi)
    jzr = jn[0]@zr-jn[1]@zi
    jzi = jn[0]@zi+jn[1]@zr
    quadratic_real_num = sum(int(a)*int(b) for a, b in zip(zr, jzr)) \
        + sum(int(a)*int(b) for a, b in zip(zi, jzi))
    rayleigh = F(quadratic_real_num, 4*PROPOSAL_D*z_norm_num)
    require(z_norm_num > 0 and rayleigh < F(-3, 50), "strict rational negative direction")
    groups["rational_normalized_Choi_witness"] = {
        "pass": True, "Choi_row": "i*4+a", "Choi_column": "j*4+b",
        "Choi_normalization": "1/4", "vector_denominator": VECTOR_D,
        "vector_real": [int(a) for a in zr], "vector_imag": [int(a) for a in zi],
        "vector_norm_squared": str(F(z_norm_num, VECTOR_D**2)),
        "exact_real_Rayleigh": str(rayleigh), "strict_Rayleigh_upper": "-3/50",
        "diagnostic_normalized_Choi_min_eigenvalue": float(evals[0]),
        "diagnostic_smin_Shat": float(np.linalg.svd(sf, compute_uv=False)[-1]),
        "diagnostic_X_Frobenius": float(np.sqrt(float(xn2))),
        "diagnostic_Y_Frobenius": float(np.sqrt(float(yn2))),
        "floating_eigenvalue_is_strict_certificate": False,
    }

    eps = F("0.000109379990")
    xi = 4*eps+q
    actual_inverse_residual = r+6*xi
    actual_choi_error = 8*(r+5*xi)
    require(actual_inverse_residual < F(1, 4), "actual inverse exists")
    require(actual_choi_error < F(9, 500), "actual Choi transport")
    actual_rayleigh_upper = rayleigh+actual_choi_error
    require(actual_rayleigh_upper < F(-21, 500), "strict actual non-CP margin")
    groups["same_clock_transport_and_actual_non_CP"] = {
        "pass": True, "epsilon_each_half_diamond": str(eps),
        "Liouville_actual_to_fixed_bound_xi": str(xi),
        "actual_inverse_residual_bound": str(actual_inverse_residual),
        "actual_inverse_operator_norm_upper": 8,
        "normalized_Choi_actual_to_Y_bound": str(actual_choi_error),
        "rational_actual_Rayleigh_upper": str(actual_rayleigh_upper),
        "certified_actual_Rayleigh_less_than": "-21/500",
        "both_comparators_before_copy": True, "actual_copy_tails_retained": True,
    }

    groups["scope_and_no_new_machine"] = {
        "pass": True, "certified": "A-only family is not CP-divisible at t_s=301,t_t=601",
        "A_pair_trace_distance_backflow_claim": False,
        "full_AB_non_CP_divisibility_claim": False,
        "any_time_local_generator_impossible_claim": False,
        "finite_integrated_positive_GKLS_interface_qualified": False,
        "intermediate_measurement_instrument_built": False,
        "third_law_or_thermodynamics_failure_claim": False,
        "new_cognitive_axioms": 0,
    }
    require(len(groups) == 8 and all(g["pass"] for g in groups.values()), "check groups")
    return {
        "round": 1216, "status": "AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE",
        "all_checks_passed": True, "check_groups": 8, "source_sha256": hashes,
        "groups": groups,
        "exact_matrix_certificate": {
            "S_That_denominator": D, "X_Y_denominator": PROPOSAL_D,
            "Shat": encode(shat), "That": encode(that), "X": encode(x), "Y": encode(y),
            "ordering": "column vec; normalized Choi input-major then output",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    target = HERE/"results.json"
    if args.write and target.exists():
        raise SystemExit("Refusing --write: results.json already exists; use read-only replay.")
    result = compute()
    if args.write:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
    stored = json.loads(target.read_text(encoding="utf-8"))
    require(result == stored, "saved results mismatch; no file overwritten")
    print(json.dumps({"status": "PASS", "all_checks_passed": True,
                      "results_match": True, "check_groups": 8,
                      "mode": "write_and_replay" if args.write else "read_only_replay"}))


if __name__ == "__main__":
    main()
