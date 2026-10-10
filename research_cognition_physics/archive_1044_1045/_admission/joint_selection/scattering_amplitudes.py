"""Two finite-mass particles: bounded barrier, complete outgoing instrument.

Numerics check formulas, not the analytic all-g / finite-time proof.
Only NumPy and Python standard library. Default is read-only verification;
--write creates the colocated results file.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

LO, HI = 0.99, 1.01
WIDTH = HI - LO


def amplitudes(k, g):
    k = np.asarray(k, dtype=float)
    q = np.sqrt(g - k * k)
    denominator = np.cosh(q) + 1j * (g - 2 * k * k) * np.sinh(q) / (2 * k * q)
    t = np.exp(-1j * k) / denominator
    r = -1j * g * np.sinh(q) / (2 * k * q * denominator)
    return r, t


def match_boundary(k, g):
    """Independent four-equation continuity solve; unknowns r,A,B,t."""
    q = np.sqrt(g - k * k)
    ep, em, ek = np.exp(q), np.exp(-q), np.exp(1j * k)
    matrix = np.array([
        [1, -1, -1, 0],
        [-1j * k, -q, q, 0],
        [0, ep, em, -ek],
        [0, q * ep, -q * em, -1j * k * ek],
    ], dtype=complex)
    return np.linalg.solve(matrix, np.array([-1, -1j * k, 0, 0], dtype=complex))


def quadrature(n):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    k = (LO + HI) / 2 + WIDTH * nodes / 2
    weights *= WIDTH / 2
    f = np.sqrt(8 / (3 * WIDTH)) * np.sin(np.pi * (k - LO) / WIDTH) ** 2
    return k, weights, f


def run(n=256):
    k, weights, f = quadrature(n)
    psi = np.sqrt(weights) * f
    distribution = psi * psi
    normalization = float(np.sum(distribution))
    rows = []
    all_boundary_errors = []
    for g in (2.0, 3.0, 4.0):
        r, t = amplitudes(k, g)
        reflectance = np.abs(r) ** 2
        q = np.sqrt(g - k * k)
        F = g * g * np.sinh(q) ** 2 / (4 * k * k * q * q)
        dlogF = 2 / g + (q / np.tanh(q) - 1) / (q * q)

        # Row order: internal 0/1, then motion T/R, then quadrature mode.
        WT = np.zeros((4 * n, 2), dtype=complex)
        WR = np.zeros_like(WT)
        WT[:n, 0] = psi
        WT[2 * n:3 * n, 1] = t * psi
        WR[3 * n:4 * n, 1] = r * psi
        completeness = WT.conj().T @ WT + WR.conj().T @ WR
        effectR = WR.conj().T @ WR
        effectT = WT.conj().T @ WT
        outgoing = WT + WR

        # A coherent unknown-input witness with a retained two-level reference.
        bell = np.array([1, 0, 0, 1], complex) / np.sqrt(2)
        outT = np.kron(WT, np.eye(2)) @ bell
        outR = np.kron(WR, np.eye(2)) @ bell
        instrument_trace = float(np.vdot(outT, outT).real + np.vdot(outR, outR).real)

        # COM K mean 0; marginal momentum moments need no COM sampling.
        pR = float(distribution @ reflectance)
        impulse = float(distribution @ (2 * k * reflectance))
        second_impulse_channel = float(distribution @ (4 * k * k * reflectance))
        coherent_t = complex(distribution @ t)
        conditional_mean_k_R = float(distribution @ (k * reflectance)) / pR
        target_mean_out_data1 = -float(distribution @ k) + impulse

        # Reflection flips relative momentum, leaves K and H0 unchanged.
        # Fixed COM support makes the selected local probe PVM equal R/T.
        probe_T_min = -0.1 / 2 + LO
        probe_R_max = 0.1 / 2 - LO
        target_R_min = -0.1 / 2 + LO
        target_T_max = 0.1 / 2 - LO

        sample_errors = []
        for kk in np.linspace(LO, HI, 17):
            rr, tt = amplitudes(kk, g)
            independent = match_boundary(kk, g)
            sample_errors.append(float(max(abs(independent[0] - rr), abs(independent[3] - tt))))
        all_boundary_errors.extend(sample_errors)
        rows.append({
            "g": g,
            "reflection_probability_data1": pR,
            "transmission_probability_data1": 1 - pR,
            "target_mean_impulse_data1": impulse,
            "target_mean_out_data1": target_mean_out_data1,
            "reflection_weighted_4k2": second_impulse_channel,
            "conditional_target_relative_momentum_R": conditional_mean_k_R,
            "transmitted_coherence_overlap": [coherent_t.real, coherent_t.imag],
            "logical_sqrt_instrument_coherence": float(np.sqrt(1 - pR)),
            "sqrt_instrument_phasefree_magnitude_gap": float(np.sqrt(1 - pR) - abs(coherent_t)),
            "completeness_error": float(np.max(np.abs(completeness - np.eye(2)))),
            "effectR_eigenvalues": np.linalg.eigvalsh(effectR).tolist(),
            "effectT_eigenvalues": np.linalg.eigvalsh(effectT).tolist(),
            "outgoing_isometry_error": float(np.max(np.abs(outgoing.conj().T @ outgoing - np.eye(2)))),
            "retained_reference_instrument_trace": instrument_trace,
            "flux_error": float(np.max(np.abs(np.abs(r) ** 2 + np.abs(t) ** 2 - 1))),
            "reflectance_formula_error": float(np.max(np.abs(reflectance - F / (1 + F)))),
            "dlogF_positive_sample_min": float(np.min(dlogF)),
            "boundary_match_error": max(sample_errors),
            "local_probe_T_momentum_lower": probe_T_min,
            "local_probe_R_momentum_upper": probe_R_max,
            "target_R_momentum_lower": target_R_min,
            "target_T_momentum_upper": target_T_max,
        })

    # Explicit interval enclosure; no numerical maximization or continuity argument.
    def bound_F(g, upper):
        qmin, qmax = np.sqrt(g - HI * HI), np.sqrt(g - LO * LO)
        if upper:
            return g * g * np.sinh(qmax) ** 2 / (4 * LO * LO * qmin * qmin)
        return g * g * np.sinh(qmin) ** 2 / (4 * HI * HI * qmax * qmax)

    F2upper, F4lower = bound_F(2, True), bound_F(4, False)
    gap_enclosure = F4lower / (1 + F4lower) - F2upper / (1 + F2upper)
    result = {
        "scope": "Asymptotic full-state instrument and analytic amplitude checks; finite-time certificate supplied separately.",
        "quadrature_order": n,
        "normalization": normalization,
        "rows": rows,
        "g4_minus_g2_probability_gap": rows[2]["reflection_probability_data1"] - rows[0]["reflection_probability_data1"],
        "interval_probability_gap_lower": float(gap_enclosure),
        "max_boundary_match_error": max(all_boundary_errors),
        "free_energy_support_upper": 0.1 ** 2 / 4 + HI ** 2,
        "interaction_operator_norm_upper": 4.0,
    }
    assert abs(normalization - 1) < 2e-13
    assert gap_enclosure > 0.30
    assert max(all_boundary_errors) < 1e-13
    for row in rows:
        for key in ("completeness_error", "outgoing_isometry_error", "flux_error", "reflectance_formula_error"):
            assert row[key] < 3e-13
        assert abs(row["retained_reference_instrument_trace"] - 1) < 3e-13
        assert row["dlogF_positive_sample_min"] > 0
        assert row["sqrt_instrument_phasefree_magnitude_gap"] > 0
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    refinement = run(384)
    delta = max(abs(a["reflection_probability_data1"] - b["reflection_probability_data1"])
                for a, b in zip(result["rows"], refinement["rows"]))
    result["quadrature_refinement_probability_error"] = delta
    assert delta < 2e-13
    output = Path(__file__).with_name("scattering_amplitudes_results.json")
    if args.write:
        output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    else:
        old = json.loads(output.read_text(encoding="utf-8"))
        for row, prev in zip(result["rows"], old["rows"]):
            assert abs(row["reflection_probability_data1"] - prev["reflection_probability_data1"]) < 1e-12
            assert abs(row["target_mean_impulse_data1"] - prev["target_mean_impulse_data1"]) < 1e-12
    print(json.dumps({"passed": True, "normalization": result["normalization"],
                      "probability_gap": result["g4_minus_g2_probability_gap"],
                      "uniform_gap_lower": result["interval_probability_gap_lower"],
                      "boundary_match_error": result["max_boundary_match_error"],
                      "refinement_error": delta}, ensure_ascii=False))


if __name__ == "__main__":
    main()
