"""Admission diagnostic only; identities are inherited from 956 and 1042.

No apparatus implementation or new scientific round is claimed.  Default runs
checks and prints a small JSON report. --write-results saves the same report.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def evolve(h, t):
    e, u = np.linalg.eigh(h)
    return (u * np.exp(-1j * t * e)) @ u.conj().T


def model(d):
    J = 1.0
    U = J * (1 - 2 * d) / d
    v = J * np.sqrt((1 - d) / d) / 2
    # Basis: bare singlet, three triplets, symmetric/antisymmetric doublon.
    h = np.diag([0.0, 0.0, 0.0, 0.0, U, U])
    h[0, 4] = h[4, 0] = -2 * v
    hpos = h + J * np.eye(6)
    D = np.diag([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
    W = np.zeros((6, 4))
    W[0, 0], W[4, 0] = np.sqrt(1 - d), np.sqrt(d)
    W[1:4, 1:4] = np.eye(3)
    P = np.diag([1.0, 0.0, 0.0, 0.0])
    hcode = J * (np.eye(4) - P)
    A = [(np.eye(6) - D) @ W, D @ W]
    energy_after = sum(a.conj().T @ hpos @ a for a in A)
    transfer = energy_after - hcode
    identities = {
        "isometry": opnorm(W.T @ W - np.eye(4)),
        "intertwining": opnorm(hpos @ W - W @ hcode),
        "charge_effect": opnorm(W.T @ D @ W - d * P),
        "measurement_material_energy_change": opnorm(
            transfer - 2 * J * (1 - d) * P
        ),
    }
    # Verify full initial matrix-valued spectral probabilities independently.
    eigenvalues, eigenvectors = np.linalg.eigh(hpos)
    spectral_checks = []
    for e in sorted(set(np.round(eigenvalues, 10))):
        selected = np.abs(eigenvalues - e) < 1e-8
        proj = eigenvectors[:, selected] @ eigenvectors[:, selected].T
        compressed = W.T @ proj @ W
        target = P if abs(e) < 1e-8 else np.eye(4) - P if abs(e - J) < 1e-8 else np.zeros((4, 4))
        spectral_checks.append(opnorm(compressed - target))
    identities["initial_matrix_spectral_measure"] = max(spectral_checks)
    for t in [0.0, 0.3, 1.0]:
        identities[f"evolution_{t}"] = opnorm(
            evolve(hpos, t) @ W - W @ evolve(hcode, t)
        )
    # Keep a passive entangled reference through the complete physical instrument.
    logical_ref = np.zeros(8, dtype=complex)
    logical_ref[0], logical_ref[3] = 1 / np.sqrt(2), 1 / np.sqrt(2)
    rho = np.outer(logical_ref, logical_ref.conj())
    physical = np.kron(W, np.eye(2))
    before = physical @ rho @ physical.conj().T
    after = sum(np.kron(a, np.eye(2)) @ rho @ np.kron(a, np.eye(2)).conj().T for a in A)
    Href = np.kron(hpos, np.eye(2))
    reference_transfer = float(np.trace(Href @ (after - before)).real)
    assert abs(reference_transfer - J * (1 - d)) < 1e-12
    # Two actual D readings with the physical six-state evolution in between.
    t = np.pi / 4
    branch11 = D @ evolve(hpos, t) @ D @ W[:, 0]
    p11 = float(np.vdot(branch11, branch11).real)
    omega = J / d
    p11_formula = d * (1 - 4 * d * (1 - d) * np.sin(omega * t / 2) ** 2)
    assert abs(p11 - p11_formula) < 1e-12
    assert max(identities.values()) < 1e-12
    return {
        "d": d,
        "J": J,
        "U": U,
        "v": float(v),
        "dimensionless_U_over_v": float(U / v),
        "full_positive_spectrum": eigenvalues.tolist(),
        "initial_spectral_weights_on_code": "P_s at 0; I-P_s at J; zero elsewhere",
        "charge_probability_on_singlet": d,
        "material_energy_change_on_singlet": float(2 * J * (1 - d)),
        "material_energy_change_with_passive_reference": reference_transfer,
        "two_readings_p11_t_pi_over_4": p11,
        "identity_residuals": identities,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    rows = [model(1 / 4), model(1 / 8)]
    n, threshold = 256, 48
    # Fixed inference rule: report d=1/4 iff the number of positive readings >48.
    errors = []
    for p in [1 / 4, 1 / 8]:
        pmf = [math.comb(n, k) * p**k * (1 - p) ** (n - k) for k in range(n + 1)]
        errors.append(float(sum(pmf[:threshold + 1]) if p == 1 / 4 else sum(pmf[threshold + 1:])))
        assert abs(sum(pmf) - 1) < 1e-12
    report = {
        "status": "not_admitted_as_new_scientific_round",
        "models": rows,
        "fixed_readout_probability_gap": 1 / 8,
        "two_reading_p11_gap": abs(rows[0]["two_readings_p11_t_pi_over_4"] - rows[1]["two_readings_p11_t_pi_over_4"]),
        "conditional_calibration_diagnostic": {
            "given_permissions": "independent ground preparations and physical D instruments; neither is internally constructed",
            "copies": n,
            "threshold_positive_count": threshold,
            "errors_d_one_quarter_then_one_eighth": errors,
            "mean_total_material_energy_change": [n * r["material_energy_change_on_singlet"] for r in rows],
            "scope": "parameter identification, not natural-coupling selection or apparatus energy accounting",
        },
        "scientific_increment": "none; 956 and 1042 identities reparameterized to reject a weak admission route",
    }
    text = json.dumps(report, ensure_ascii=False, indent=2)
    if args.write_results:
        Path(__file__).with_name("material_same_spectrum_audit_results.json").write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
