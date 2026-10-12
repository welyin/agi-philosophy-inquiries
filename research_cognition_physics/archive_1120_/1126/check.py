"""1126: compact CAR loading, spatial receiving, and physical-field budgets.

Default is read-only recomputation. --write-results registers results.json.
Finite CAR matrices test algebra only, not a cutoff of the continuum Dirac
vacuum or a microscopic realization of the declared regional gate windows.
The compact-envelope budgets use exact rational polynomial integrals.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.json"
TOL = 3e-10


def residual(x, label, tol=TOL):
    value = float(np.linalg.norm(x))
    if not math.isfinite(value) or value > tol:
        raise AssertionError(f"{label}: {value}")
    return value


def annihilators(n):
    dim = 1 << n
    ans = []
    for j in range(n):
        c = np.zeros((dim, dim), complex)
        for state in range(dim):
            if state & (1 << j):
                sign = (-1) ** ((state & ((1 << j) - 1)).bit_count())
                c[state ^ (1 << j), state] = sign
        ans.append(c)
    return ans


def hop_swap(c, d):
    # The exact generator has spectrum -1,0,1 for normalized orthogonal CAR
    # modes. Thus exp(-i*pi*G/2)=I-iG-G^2, with no exponentiation cutoff.
    gen = d.conj().T @ c + c.conj().T @ d
    return np.eye(c.shape[0]) - 1j * gen - gen @ gen


def second_quantize(mat, cs):
    out = np.zeros_like(cs[0])
    for j in range(len(cs)):
        for k in range(len(cs)):
            if mat[j, k] != 0:
                out += mat[j, k] * (cs[j].conj().T @ cs[k])
    return out


def half_trace(a):
    return float(np.sum(abs(np.linalg.eigvalsh((a + a.conj().T) / 2))) / 2)


def receiver_unit(out, i, j):
    # Receiver occupies the last/high two bits. All other modes are retained
    # in out and traced only for this separately declared receiver channel.
    return out[:, i].reshape(4, -1) @ out[:, j].reshape(4, -1).conj().T


def receiver_reference(out):
    # Apply the complete isometry to one half of normalized Bell input.
    wave = out.reshape(4, -1, 2).transpose(0, 2, 1).reshape(8, -1) / math.sqrt(2)
    return wave @ wave.conj().T


def group_full_car_channel():
    # d1,d2,u1,v1,u2,v2,r1,r2. The v modes are environmental modes.
    # u is the EXACT transported source mode in this algebraic frame.
    cs = annihilators(8)
    d1, d2, u1, v1, u2, v2, r1, r2 = cs
    ident = np.eye(256)
    vac = ident[:, 0]
    load = hop_swap(u2, d2) @ hop_swap(u1, d1)
    ideal_read = hop_swap(u2, r2) @ hop_swap(u1, r1)
    ideal = ideal_read @ load
    r = {
        "load_unitarity": residual(load.conj().T @ load - ident, "source swap"),
        "ideal_receiver_unitarity": residual(ideal.conj().T @ ideal - ident, "ideal chain"),
        "exact_receiver_r1": residual(ideal.conj().T @ r1 @ ideal + d1, "r1 out=-d1"),
        "exact_receiver_r2": residual(ideal.conj().T @ r2 @ ideal + d2, "r2 out=-d2"),
    }
    # A full, normal, even-parity two-particle Slater state. Its local u
    # marginals are occupied with probabilities p1,p2 and correlated with v.
    p1, p2 = 0.35, 0.62
    create1 = math.sqrt(p1) * u1.conj().T + math.sqrt(1 - p1) * v1.conj().T
    create2 = math.sqrt(p2) * u2.conj().T + math.sqrt(1 - p2) * v2.conj().T
    env = create1 @ create2 @ vac
    source = np.column_stack([d1.conj().T @ env, d2.conj().T @ env])
    r["full_input_isometry"] = residual(source.conj().T @ source - np.eye(2), "source isometry")
    full_ideal = ideal @ source
    code = np.eye(4)[:, [1, 2]]
    matrix_units = 0.0
    for i in range(2):
        for j in range(2):
            target = np.outer(code[:, i], code[:, j])
            matrix_units = max(matrix_units, residual(receiver_unit(full_ideal, i, j) - target,
                                                      "ideal receiver matrix unit"))
    r["all_receiver_matrix_units"] = matrix_units
    # The same environmental vector factors from both logical columns.
    environment = full_ideal[:, 0].reshape(4, -1)[1]
    target_full = np.column_stack([np.kron(code[:, i], environment) for i in range(2)])
    r["complement_independent_of_unknown_Q"] = residual(full_ideal - target_full,
                                                       "full environment factor")
    ideal_bell = receiver_reference(full_ideal)
    bell = np.zeros(8, complex)
    bell[1 * 2 + 0] = 1 / math.sqrt(2)
    bell[2 * 2 + 1] = 1 / math.sqrt(2)
    r["Bell_reference_recovery"] = residual(ideal_bell - np.outer(bell, bell.conj()), "Bell recovery")

    angle = 0.07
    g1 = math.cos(angle) * u1 + math.sin(angle) * v1
    g2 = math.cos(angle) * u2 + math.sin(angle) * v2
    actual_read = hop_swap(g2, r2) @ hop_swap(g1, r1)
    actual = actual_read @ load
    full_actual = actual @ source
    mode_distance = 2 * abs(math.sin(angle / 2))
    analytical_bound = math.pi * mode_distance
    unitary_distance = float(np.linalg.norm(actual_read - ideal_read, 2))
    stinespring_distance = float(np.linalg.norm(full_actual - full_ideal, 2))
    bell_error = half_trace(receiver_reference(full_actual) - ideal_bell)
    if max(unitary_distance, stinespring_distance, bell_error) > analytical_bound + TOL:
        raise AssertionError("Stinespring/Duhamel bound failed")
    rng = np.random.default_rng(1126)
    qr = rng.normal(size=(2, 5)) + 1j * rng.normal(size=(2, 5))
    qr /= np.linalg.norm(qr)
    reference_vector_error = float(np.linalg.norm((full_actual - full_ideal) @ qr))
    if reference_vector_error > stinespring_distance + TOL:
        raise AssertionError("passive reference bound")
    recv_mixed = sum(receiver_unit(full_actual, i, i) for i in range(2)) / 2
    failure_probability = float((recv_mixed[0, 0] + recv_mixed[3, 3]).real)

    # Necessary opposite test: an orthogonal receiver has no source channel.
    # This also prevents identifying global invertibility with local access.
    miss = hop_swap(v2, r2) @ hop_swap(v1, r1) @ load @ source
    miss_r = max(residual(receiver_unit(miss, 0, 0) - receiver_unit(miss, 1, 1), "missed diagonal"),
                 residual(receiver_unit(miss, 0, 1), "missed coherence"),
                 residual(receiver_unit(miss, 1, 0), "missed conjugate coherence"))
    return {
        "scope": "eight-mode exact CAR algebra, with correlated occupied environment; not continuum Dirac truncation",
        "dimension": 256,
        "local_field_occupations": [p1, p2],
        "receiver_mismatch_angle": angle,
        "mode_distance": mode_distance,
        "analytic_half_diamond_upper_bound": analytical_bound,
        "whole_unitary_distance": unitary_distance,
        "complete_stinespring_distance": stinespring_distance,
        "Bell_output_half_trace_distance": bell_error,
        "five_dimensional_reference_vector_error": reference_vector_error,
        "receiver_vacuum_or_double_occupation_probability": failure_probability,
        "orthogonal_receiver_constant_channel_residual": miss_r,
        "residuals": r,
    }


def derivative(p):
    return [i * p[i] for i in range(1, len(p))]


def polynomial_norm_squared(p):
    return sum((a * b * Fraction(2, i + j + 1)
                for i, a in enumerate(p) for j, b in enumerate(p) if (i + j) % 2 == 0), Fraction(0))


def exact_envelope():
    p = [Fraction(0)] * 9
    for k in range(5):
        p[2 * k] = Fraction((-1) ** k * math.comb(4, k))
    i0 = polynomial_norm_squared(p)
    d1 = polynomial_norm_squared(derivative(p)) / i0
    d2 = polynomial_norm_squared(derivative(derivative(p))) / i0
    if (i0, d1, d2) != (Fraction(65536, 109395), Fraction(34, 7), Fraction(459, 7)):
        raise AssertionError("exact compact polynomial integrals")
    return i0, d1, d2


def dirac_matrices():
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1., -1.]).astype(complex)
    zero = np.zeros((2, 2), complex)
    alphas = [np.block([[zero, s], [s, zero]]) for s in (sx, sy, sz)]
    beta = np.diag([1., 1., -1., -1.]).astype(complex)
    chi = np.array([1., 0., 1., 0.], complex) / math.sqrt(2)
    return alphas, beta, chi


def group_continuous_budget():
    i0, d1_exact, d2_exact = exact_envelope()
    d1, d2 = float(d1_exact), float(d2_exact)
    m, r, a, travel = 1., 2., 3., 12.
    a_half_sq = m * m + 2 * d1
    a_norm_sq = m ** 4 + 4 * m * m * d1 + 2 * d2 + 2 * d1 * d1
    constant = 2 * math.sqrt(a_half_sq) + travel * math.sqrt(a_norm_sq) + 4 * math.sqrt(d1)
    carrier = 57975
    error_bound = math.pi * constant / carrier
    if not error_bound < 0.01:
        raise AssertionError("fixed finite carrier does not meet declared threshold")
    per_gate = {
        "energy": 2 * math.sqrt(carrier ** 2 + 2 * (m * m + 3 * d1)),
        "p_x": 2 * math.sqrt(2 * d1),
        "p_y": 2 * math.sqrt(2 * d1),
        "p_z": 2 * math.sqrt(carrier ** 2 + 2 * d1),
    }
    alphas, beta, chi = dirac_matrices()
    point_res = residual(alphas[2] @ chi - chi, "right-moving spinor")
    # Independent pointwise continuum symbol diagnostics, at moderate K.
    # They are NOT quadrature estimates of the norm or certified constants.
    test_k, test_t = 8., 1.3
    points = [(0., 0., 8.), (2., -3., 4.), (1., 0., 13.), (0., 2., -2.), (2., 3., 0.)]
    rows = []
    for values in points:
        p = np.array(values)
        h = sum(p[j] * alphas[j] for j in range(3)) + m * beta
        omega = math.sqrt(float(p @ p) + m * m)
        pminus = (np.eye(4) - h / omega) / 2
        neg = float(np.linalg.norm(pminus @ chi) ** 2)
        formula = (1 - p[2] / omega) / 2
        point_res = max(point_res, residual(neg - formula, "negative-spectrum weight"),
                        residual(h @ h - omega ** 2 * np.eye(4), "Dirac dispersion"))
        evolved = math.cos(omega * test_t) * chi - 1j * math.sin(omega * test_t) * h @ chi / omega
        translated = np.exp(-1j * p[2] * test_t) * chi
        difference = float(np.linalg.norm(evolved - translated))
        good = bool(p[2] >= test_k / 2)
        transverse = float(p[0] ** 2 + p[1] ** 2 + m * m)
        bound = (2 * math.sqrt(transverse) + test_t * transverse) / test_k if good else 2.
        if difference > bound + TOL:
            raise AssertionError("Dirac good/bad region point bound")
        rows.append({"momentum": list(values), "negative_spectral_weight": neg,
                     "mode_symbol_difference": difference, "good_region": good, "analytic_point_bound": bound})
    return {
        "scope": "exact rational H2 compact-envelope budgets; point values only diagnose continuum symbols",
        "envelope_is_C_infinity": False,
        "I0_exact": str(i0), "d1_exact": str(d1_exact), "d2_exact": str(d2_exact),
        "mass": m, "source_and_receiver_radius": r, "window_half_duration": a,
        "center_separation_and_center_time": travel, "total_time": travel + 2 * a,
        "boundary_spatial_envelope_radius": r + a,
        "spatial_envelope_separation_margin": travel - 2 * (r + a),
        "envelope_inside_ball_margin": r - math.sqrt(3),
        "propagation_constant_C": constant, "fixed_carrier_K": carrier,
        "mode_error_upper_bound": constant / carrier,
        "recovery_half_diamond_upper_bound": error_bound,
        "per_gate_field_four_momentum_defect_norms": per_gate,
        "complete_process_field_four_momentum_upper_bounds": {k: 2 * v for k, v in per_gate.items()},
        "field_mean_energy_upper_bound": 2 * per_gate["energy"],
        "field_second_energy_moment_upper_bound": (2 * per_gate["energy"]) ** 2,
        "pointwise_symbol_residual": point_res,
        "pointwise_diagnostics": rows,
    }


def group_four_momentum():
    cs = annihilators(6)
    identity = np.eye(64)
    # Per flavor order f,d,orthogonal derivative mode. Free compressed G
    # need not preserve these three modes; the DIFFERENCE does close here.
    swap = hop_swap(cs[3], cs[4]) @ hop_swap(cs[0], cs[1])
    v3 = np.array([[0, -1j, 0], [-1j, 0, 0], [0, 0, 1]], complex)
    _, d1_exact, _ = exact_envelope()
    d1 = float(d1_exact)
    cases = [("energy", 3., math.sqrt(1 + 3 * d1)),
             ("p_z", 3., math.sqrt(d1)), ("p_x", 0., math.sqrt(d1)),
             ("p_y", 0., math.sqrt(d1))]
    rows = []
    for label, mu, sigma in cases:
        g3 = np.array([[mu, 0, sigma], [0, 0, 0], [sigma, 0, -0.6]], complex)
        b3 = v3.conj().T @ g3 @ v3 - g3
        expected = np.array([[-mu, 0, -sigma], [0, mu, 1j * sigma],
                             [-sigma, -1j * sigma, 0]], complex)
        g = np.kron(np.eye(2), g3)
        b = np.kron(np.eye(2), b3)
        qg, qb = second_quantize(g, cs), second_quantize(b, cs)
        direct = swap.conj().T @ qg @ swap - qg
        lam = math.sqrt(mu * mu + 2 * sigma * sigma)
        exact_norm = 2 * lam
        norm = float(max(abs(np.linalg.eigvalsh(qb))))
        vals, vecs = np.linalg.eigh(g)
        covariance = (vecs * (vals < -1e-12)) @ vecs.conj().T
        trace_shift = float(np.trace(covariance @ b).real)
        normal = qb - trace_shift * identity
        res = max(residual(b3 - expected, "three-mode physical difference"),
                  residual(np.linalg.eigvalsh(b3) - [-lam, 0, lam], "three-mode spectrum"),
                  residual(direct - qb, "whole CAR physical difference"),
                  residual(norm - exact_norm, "whole CAR defect norm"),
                  residual(normal + trace_shift * identity - direct, "normal ordering scalar kept"))
        rows.append({"component": label, "mu": mu, "sigma": sigma,
                     "one_particle_difference_norm": lam,
                     "whole_CAR_defect_norm": norm, "normal_order_scalar": trace_shift,
                     "residual": res})
    return {"scope": "six-mode algebra diagnostic; compression difference, not continuum spectrum truncation",
            "dimension": 64, "cases": rows}


def calculate():
    return {"round": 1126, "scientific_calibration_increment": 0, "groups_passed": 3,
            "complete_CAR_receiver_channel": group_full_car_channel(),
            "continuous_compact_budget": group_continuous_budget(),
            "complete_field_four_momentum": group_four_momentum(),
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def compare(current, saved, path="root"):
    if isinstance(current, dict):
        if not isinstance(saved, dict) or current.keys() != saved.keys():
            raise AssertionError(f"result keys differ at {path}")
        for key in current:
            compare(current[key], saved[key], f"{path}.{key}")
    elif isinstance(current, list):
        if not isinstance(saved, list) or len(current) != len(saved):
            raise AssertionError(f"result list differs at {path}")
        for i, (a, b) in enumerate(zip(current, saved)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(current, (float, int)) and not isinstance(current, bool):
        if not isinstance(saved, (float, int)) or not math.isclose(current, saved, rel_tol=2e-9, abs_tol=3e-10):
            raise AssertionError(f"numeric result differs at {path}: {current} versus {saved}")
    elif current != saved:
        raise AssertionError(f"result differs at {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_results:
        RESULTS.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        status = "written"
    else:
        if not RESULTS.exists():
            raise SystemExit("No saved results; first registration requires --write-results.")
        compare(result, json.loads(RESULTS.read_text(encoding="utf-8")))
        status = "matched read-only"
    budget = result["continuous_compact_budget"]
    print(json.dumps({"round": 1126, "groups_passed": 3, "saved_results": status,
                      "recovery_half_diamond_upper_bound": budget["recovery_half_diamond_upper_bound"],
                      "field_mean_energy_upper_bound": budget["field_mean_energy_upper_bound"],
                      "finite_CAR_Bell_error": result["complete_CAR_receiver_channel"]["Bell_output_half_trace_distance"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
