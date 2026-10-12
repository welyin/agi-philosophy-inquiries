"""1125: support audit of the actual 1111 control reused at a distant source.

Default: read-only recomputation and comparison. --write-results writes only
this round's results.json. Continuous support and energy-moment claims are
proved in proof.md. The oscillator is a Weyl/phase diagnostic, not a spatial
KG cutoff or a physical implementation of the controller or measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / "archive_1106_1119" / "1111" / "check.py"
RESULTS = HERE / "results.json"
TOL = 2e-10


def residual(x, label, tol=TOL):
    value = float(np.linalg.norm(x))
    if not math.isfinite(value) or value > tol:
        raise AssertionError(f"{label}: {value}")
    return value


def exp_minus(h, t=1.0):
    vals, vecs = np.linalg.eigh(h)
    return (vecs * np.exp(-1j * t * vals)) @ vecs.conj().T


def pure_distance(v, w):
    # Both arguments are normalized; use the absolute overlap formula and
    # clamp floating-point roundoff at unit overlap.
    overlap = np.vdot(v, w)
    return math.sqrt(max(0.0, 1.0 - min(1.0, abs(overlap) ** 2)))


def group_complete_control():
    spec = importlib.util.spec_from_file_location("round1111_reused_1125", OLD)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    m = old.model()
    hd, gamma, hc = m["hout"], m["gamma"], m["hc"]
    pc = np.kron(np.eye(8), np.kron(np.diag([0., 0., 1.]), np.eye(3)))
    pd = gamma @ np.kron(np.eye(4), pc) @ gamma.conj().T
    delta = math.pi / 8
    source_basis = m["phase"] * m["ready"]
    clock_h = hc + m["baseline"] * np.eye(4)
    ws = gamma @ np.kron(exp_minus(clock_h, math.pi - delta), np.eye(72)) @ source_basis
    wf = gamma @ np.kron(exp_minus(clock_h, math.pi), np.eye(72)) @ source_basis
    # This independently calls the old uninterrupted physical evolution.
    ws_old = m["w"](math.pi - delta)
    r = {
        "old_actual_input_match": residual(ws - ws_old, "original W_s"),
        "gamma_intertwiner": residual(hd @ gamma - gamma @ m["hin"], "Gamma intertwiner"),
        "projector": residual(pd @ pd - pd, "actual complete projector"),
        "conserved_control": residual(hd @ pd - pd @ hd, "complete control commutes H_D"),
        "source_eigenbranches": residual(pd @ ws - ws @ np.diag([0., 1.]), "reachable source branches"),
        "source_isometry": residual(ws.conj().T @ ws - np.eye(2), "source isometry"),
    }

    # Independent bounded matrices approximate q,p only on these low-energy
    # vectors. No finite matrix satisfies the CCR on its whole Hilbert space.
    cutoff = 40
    lower = np.zeros((cutoff, cutoff), complex)
    for n in range(1, cutoff):
        lower[n - 1, n] = math.sqrt(n)
    q = (lower + lower.conj().T) / math.sqrt(2)
    p = (lower - lower.conj().T) / (1j * math.sqrt(2))
    vac = np.eye(cutoff, dtype=complex)[:, 0]
    a, b = 0.8, 0.6
    c, d = 0.4, 1.2
    chi = a * b / 2
    kappa = a * d - b * c
    # q pulse first, p pulse second; the BCH phase is +ab/2.
    u_time = exp_minus(b * p) @ exp_minus(a * q)
    u_weyl = exp_minus(a * q + b * p)
    v_time = u_time @ vac
    v_weyl = u_weyl @ vac
    wf_probe = exp_minus(c * q + d * p, -1)
    effect = (np.eye(cutoff) - (wf_probe + wf_probe.conj().T) / 2) / 2
    r["ordered_source_with_phase"] = residual(v_time - np.exp(1j * chi) * v_weyl, "source phase on vacuum")
    nu = math.exp(-(c * c + d * d) / 4)
    probabilities = [float(np.vdot(v, effect @ v).real) for v in (vac, v_time)]
    expected = [(1 - nu) / 2, (1 - nu * math.cos(kappa)) / 2]
    r["weyl_local_probe_formula"] = residual(np.array(probabilities) - expected, "Weyl probe probabilities")
    r["vacuum_characteristic"] = residual(np.vdot(vac, wf_probe @ vac) - nu, "normal Gaussian characteristic")
    r["source_normalization"] = residual(np.vdot(v_time, v_time) - 1, "normalization")
    top_probability = float(np.sum(abs(v_time[-4:]) ** 2))
    if top_probability > 1e-20:
        raise AssertionError("oscillator diagnostic reaches cutoff")

    # Preserve all 288 device components; no bare receiver replacement.
    joint = np.column_stack([np.kron(wf[:, 0], vac), np.kron(wf[:, 1], v_time)])
    no_phase = np.column_stack([np.kron(wf[:, 0], vac), np.kron(wf[:, 1], v_weyl)])
    gram = joint.conj().T @ joint
    r["complete_source_isometry"] = residual(gram - np.eye(2), "whole D-field source isometry")
    rng = np.random.default_rng(1125)
    source_reference = rng.normal(size=(2, 3)) + 1j * rng.normal(size=(2, 3))
    source_reference /= np.linalg.norm(source_reference)
    out_reference = joint @ source_reference
    r["passive_reference_decode"] = residual(joint.conj().T @ out_reference - source_reference,
                                               "full source-reference decode")
    matrix_unit_residual = 0.0
    for i in range(2):
        for j in range(2):
            unit = np.zeros((2, 2), complex)
            unit[i, j] = 1
            matrix_unit_residual = max(matrix_unit_residual,
                                       residual(gram @ unit @ gram - unit, "complete matrix units"))
    r["all_source_matrix_units"] = matrix_unit_residual
    plus = np.ones(2) / math.sqrt(2)
    phase_error = pure_distance(joint @ plus, no_phase @ plus)
    r["phase_deletion_distance_formula"] = residual(phase_error - abs(math.sin(chi / 2)),
                                                    "deleted conditional phase is observable")
    return {
        "scope": "original 288-dimensional D and source; finite oscillator algebra diagnostic, not a KG spatial truncation",
        "device_dimension": 288,
        "oscillator_dimension": cutoff,
        "source_relative_phase": chi,
        "probe_symplectic_pairing": kappa,
        "probe_probabilities": probabilities,
        "probe_signal_difference": probabilities[1] - probabilities[0],
        "phase_deleted_complete_output_half_trace_distance": phase_error,
        "oscillator_top_four_level_probability": top_probability,
        "residuals": r,
    }


def bump(x):
    x = np.asarray(x)
    return np.exp(-1.0 / (1.0 - x * x))


def j1_over_z_from_series(z_squared):
    # Entire series in z^2, including its continuous value 1/2 at z=0.
    term = np.full_like(np.asarray(z_squared, dtype=float), 0.5)
    total = term.copy()
    for k in range(1, 9):
        term = -term * z_squared / (4 * k * (k + 1))
        total += term
    return total


def compact_pairing_diagnostic(order):
    # Normalized radial bumps at the same center and normalized time bumps.
    # The entire product support is strictly timelike, so no lightcone delta
    # is sampled or dropped. Integrate the smooth interior kernel only.
    nodes, weights = leggauss(order)
    tw = weights * bump(nodes)
    tw /= tw.sum()
    radius = 0.02 * (nodes + 1) / 2
    rw = weights * radius ** 2 * bump((nodes + 1) / 2)
    rw /= rw.sum()
    angles, aw = leggauss(12)
    aw /= 2
    rr = radius[:, None, None]
    ss = radius[None, :, None]
    distance_squared = rr * rr + ss * ss - 2 * rr * ss * angles[None, None, :]
    spatial_weights = rw[:, None, None] * rw[None, :, None] * aw[None, None, :]
    answer = 0.0
    for i, source_node in enumerate(nodes):
        for j, probe_node in enumerate(nodes):
            dt = 0.2 + 0.01 * (probe_node - source_node)
            tau_squared = dt * dt - distance_squared
            if float(tau_squared.min()) <= 0:
                raise AssertionError("test supports must be strictly timelike")
            interior = j1_over_z_from_series(tau_squared) / (4 * math.pi)
            answer += tw[i] * tw[j] * float(np.sum(spatial_weights * interior))
    return answer


def group_continuous_support():
    ts = 23 / 15 + math.pi - math.pi / 8
    delta = math.pi / 8
    spatial_radius, time_radius = 0.02, 0.01
    dt_min, dt_max = 0.18, 0.22
    r_max = 2 * spatial_radius
    tau_min = math.sqrt(dt_min ** 2 - r_max ** 2)
    origin_spacelike_margin = math.sqrt(37) - spatial_radius - (ts + 0.26)
    window_remaining = delta - 0.26
    near_D_spacelike_margin = 6 - 2 * spatial_radius - dt_max
    if min(origin_spacelike_margin, window_remaining, near_D_spacelike_margin, tau_min) <= 0:
        raise AssertionError("strict geometry margin missing")
    kernel_lower = (1 - dt_max ** 2 / 8) / (8 * math.pi)
    kernel_upper = 1 / (8 * math.pi)
    kappa_lower, kappa_upper = 8 * math.pi ** 2 * kernel_lower, math.pi
    signal_lower = math.exp(-0.5) * math.cos(math.pi * dt_max ** 2 / 16) ** 2
    if signal_lower <= 0.6064:
        raise AssertionError("analytic signal bound")

    # This grid is a diagnostic of the series; the alternating-series proof
    # in proof.md, not the grid, establishes the whole continuum interval.
    z = np.linspace(tau_min, dt_max, 301)
    values = j1_over_z_from_series(z * z) / (4 * math.pi)
    if float(values.min()) < kernel_lower or float(values.max()) > kernel_upper:
        raise AssertionError("Bessel analytic enclosure")
    pairing12 = compact_pairing_diagnostic(12)
    pairing24 = compact_pairing_diagnostic(24)
    quad_change = abs(pairing24 - pairing12)
    if quad_change > 5e-7:
        raise AssertionError("smooth compact-bump diagnostic not resolved")
    if not kernel_lower <= pairing24 <= kernel_upper:
        raise AssertionError("integrated pairing outside analytic enclosure")

    # Spatial b is a normalized 3D radial smooth bump. Quadrature estimates
    # eta and source strength for inspection; the theorem defines them by
    # exact integrals and does not depend on these estimates.
    nodes, weights = leggauss(96)
    u, wu = (nodes + 1) / 2, weights / 2
    one = 4 * math.pi * spatial_radius ** 3 * float(np.sum(wu * u ** 2 * bump(u)))
    two = 4 * math.pi * spatial_radius ** 3 * float(np.sum(wu * u ** 2 * bump(u) ** 2))
    b_l2_squared = two / (one * one)
    eta = math.sqrt(2 / b_l2_squared)
    source_strength = 8 * math.pi ** 2 / eta
    var_bound = eta * eta * b_l2_squared / 2
    residual(var_bound - 1, "vacuum variance upper bound normalization")
    readout_rows = []
    for error in (0.0, 0.05, 0.35):
        readout_rows.append({"uniform_probability_error": error,
                             "signal_lower_bound": max(0.0, signal_lower - 2 * error),
                             "causal_common_prediction_error_lower_bound": max(0.0, signal_lower / 2 - error)})
    return {
        "scope": "one fixed smooth compact KG source/probe; continuum assertions proved analytically",
        "field_mass": 1.0,
        "original_window_start": ts,
        "original_window_duration": delta,
        "remote_center": [1.0, 6.0, 0.0],
        "spatial_support_radius": spatial_radius,
        "time_support_halfwidth": time_radius,
        "source_to_probe_timelike_tau_bounds": [tau_min, dt_max],
        "source_origin_to_latest_probe_spacelike_margin": origin_spacelike_margin,
        "window_remaining_after_probe": window_remaining,
        "near_D_source_to_probe_spacelike_margin": near_D_spacelike_margin,
        "interior_kernel_absolute_bounds": [kernel_lower, kernel_upper],
        "absolute_symplectic_pairing_bounds": [kappa_lower, kappa_upper],
        "normalized_pairing_quadrature_12": pairing12,
        "normalized_pairing_quadrature_24": pairing24,
        "quadrature_change_diagnostic_only": quad_change,
        "absolute_kappa_quadrature_diagnostic": 8 * math.pi ** 2 * pairing24,
        "b_l2_squared_quadrature_estimate": b_l2_squared,
        "eta_quadrature_estimate": eta,
        "source_strength_quadrature_estimate": source_strength,
        "vacuum_variance_analytic_upper_bound": var_bound,
        "vacuum_characteristic_lower_bound": math.exp(-0.5),
        "analytic_signal_lower_bound": signal_lower,
        "analytic_causal_error_lower_bound": signal_lower / 2,
        "imperfect_readout_bounds": readout_rows,
    }


def group_local_controls():
    # Factor order: left local control, right local control, passive Z, probe.
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], complex)
    zero, one = I[:, 0], I[:, 1]
    p0, p1 = np.outer(zero, zero), np.outer(one, one)
    theta = 0.7
    right_probe = np.kron(p0, np.kron(I, I)) + np.kron(p1, np.kron(I, exp_minus(X, theta)))
    evolution = np.kron(I, right_probe)
    rng = np.random.default_rng(112503)
    v = rng.normal(size=8) + 1j * rng.normal(size=8)
    v /= np.linalg.norm(v)
    full = np.kron(v, zero)
    initial = np.outer(full, full.conj())

    def outside_left(rho):
        return np.trace(rho.reshape(2, 8, 2, 8), axis1=0, axis2=2)

    baseline = evolution @ initial @ evolution.conj().T
    operations = {
        "local_flip": [X],
        "local_dephasing": [p0, p1],
        "local_reset": [p0, np.outer(zero, one)],
    }
    rows = []
    for name, kraus in operations.items():
        residual(sum(k.conj().T @ k for k in kraus) - I, "local CPTP completeness")
        changed = np.zeros_like(initial)
        for k in kraus:
            ext = np.kron(k, np.eye(8))
            changed += ext @ initial @ ext.conj().T
        final = evolution @ changed @ evolution.conj().T
        rows.append({"operation": name,
                     "complete_outside_left_state_residual": residual(outside_left(final) - outside_left(baseline),
                                                                        "local operations cannot change pre-correlated remote control")})

    # Cat encoding preserves arbitrary Q and reference but is an already
    # distributed input, not a free encoding of a fresh distant unknown Q.
    cat = np.column_stack([np.kron(zero, zero), np.kron(one, one)])
    encoded_pd_left = np.kron(p1, I)
    encoded_pd_right = np.kron(I, p1)
    code_r = max(residual(encoded_pd_left @ cat - cat @ p1, "left code action"),
                 residual(encoded_pd_right @ cat - cat @ p1, "right code action"),
                 residual(cat.conj().T @ cat - I, "complete cat isometry"))
    mother_difference = float(np.linalg.norm(encoded_pd_left - encoded_pd_right, 2))
    residual(mother_difference - 1, "same encoded action does not identify mother operators")
    locally_flipped_zero = np.kron(X, I) @ cat[:, 0]
    remote_one_probability = float(np.vdot(locally_flipped_zero, encoded_pd_right @ locally_flipped_zero).real)
    residual(remote_one_probability, "local fresh flip does not update remote copy")
    # The analytic near-D KG comparison has exactly zero commutator by
    # support, not a small number produced by a Fourier cutoff.
    gap = 6 - 0.04 - 0.22
    if gap <= 0:
        raise AssertionError("near-D source is not spacelike to remote probe")
    return {
        "scope": "strict near-D KG support and complete local-control algebra, not an implemented relativistic apparatus",
        "near_D_source_spacelike_margin": gap,
        "near_D_source_probe_pairing_exact": 0,
        "near_D_source_probe_signal_exact": 0,
        "precorrelated_local_operation_checks": rows,
        "cat_code_isometry_and_branch_residual": code_r,
        "distinct_local_control_operators_norm_difference": mother_difference,
        "remote_copy_one_probability_after_only_local_flip": remote_one_probability,
    }


def calculate():
    return {
        "round": 1125,
        "scientific_calibration_increment": 0,
        "groups_passed": 3,
        "complete_control_and_phase": group_complete_control(),
        "continuous_KG_support_and_signal": group_continuous_support(),
        "local_support_and_preshared_controls": group_local_controls(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reused_1111_source_sha256": hashlib.sha256(OLD.read_bytes()).hexdigest(),
    }


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
        if not isinstance(saved, (float, int)) or not math.isclose(current, saved, rel_tol=2e-9, abs_tol=2e-10):
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
    print(json.dumps({"round": 1125, "groups_passed": 3, "saved_results": status,
                      "analytic_signal_lower_bound": result["continuous_KG_support_and_signal"]["analytic_signal_lower_bound"],
                      "analytic_causal_error_lower_bound": result["continuous_KG_support_and_signal"]["analytic_causal_error_lower_bound"],
                      "phase_deleted_complete_output_distance": result["complete_control_and_phase"]["phase_deleted_complete_output_half_trace_distance"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
