"""Fixed two-cell scalar-clock diagnostic; no interval/time-window certificate.

Only numpy is used. First save: --write (exclusive); default recomputes read-only.
All matrices are the explicitly declared finite effective rules, not silently
identified with compression of an interacting infinite oscillator theory.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "diagnostic_results.json"
ID3 = np.eye(3, dtype=complex)


def tensor(parts):
    out = np.array([[1.0 + 0j]])
    for item in parts:
        out = np.kron(out, item)
    return out


def embed(op, index):
    return tensor([op if j == index else ID3 for j in range(6)])


def expect(states, op):
    return np.real(np.einsum("ij,ij->j", states.conj(), op @ states))


def run():
    a = np.diag(np.sqrt([1.0, 2.0]), 1).astype(complex)
    q = (a + a.conj().T) / np.sqrt(2.0)
    p = (a - a.conj().T) / (1j * np.sqrt(2.0))
    # Exact compression of infinite HO quadratics to |0>, |1>, |2>.
    q2 = np.diag([0.5, 1.5, 2.5]).astype(complex)
    q2[0, 2] = q2[2, 0] = np.sqrt(2.0) / 2.0
    p2 = np.diag([0.5, 1.5, 2.5]).astype(complex)
    p2[0, 2] = p2[2, 0] = -np.sqrt(2.0) / 2.0
    eq, uq = np.linalg.eigh(q)
    sinq = (uq * np.sin(eq)) @ uq.conj().T
    qnorm = float(max(abs(eq)))
    Q = [embed(q, j) for j in range(6)]
    P = [embed(p, j) for j in range(6)]
    Q2 = [embed(q2, j) for j in range(6)]
    P2 = [embed(p2, j) for j in range(6)]
    S = [embed(sinq, j * 3) @ Q[3 * j + 1] @ Q[3 * j + 2]
         for j in range(2)]
    mass = [sum((Q2[j] + P2[j]) / 2.0 for j in range(3 * i, 3 * i + 3))
            for i in range(2)]
    difference_square = sum((Q[s + 3] - Q[s]) @ (Q[s + 3] - Q[s])
                            for s in range(3))
    link = difference_square / 20.0
    H0 = mass[0] + mass[1] + link
    H = H0 + (S[0] + S[1]) / 4.0
    source_N = [mass[i] + difference_square / 40.0 + S[i] / 4.0
                for i in range(2)]
    source_logv = [
        sum((Q2[j] - P2[j]) / 2.0 for j in range(3 * i, 3 * i + 3))
        + S[i] / 4.0 for i in range(2)]
    source_logw = link
    # The generator identity is a matrix identity, not inferred from energy data.
    euler_error = float(np.max(np.abs(H - source_N[0] - source_N[1])))
    qdot = [1j * (H @ Q[3 * i] - Q[3 * i] @ H) for i in range(2)]
    qdot_error = max(float(np.max(np.abs(qdot[i] - P[3 * i])))
                     for i in range(2))
    t0 = np.array([1.0, 1j, 0.0]) / np.sqrt(2.0)
    phis = [np.array([1.0, sign, 0.0]) / np.sqrt(2.0) for sign in [1, -1]]
    blank = np.array([1.0, 0.0, 0.0])
    initial = np.column_stack([
        tensor([t0, ph, blank, t0, ph, blank]).ravel() for ph in phis])
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    coeff = eigenvectors.conj().T @ initial
    def state(t):
        return eigenvectors @ (np.exp(-1j * eigenvalues * t)[:, None] * coeff)
    identity = np.eye(3**6)
    effect = (identity + Q[2] / qnorm) / 2.0
    final = state(0.4)
    times = np.linspace(0.0, 0.4, 81)
    clock_data = []
    max_energy_error = 0.0
    initial_energy = expect(initial, H)
    for t in times:
        st = state(float(t))
        max_energy_error = max(max_energy_error,
                               float(max(abs(expect(st, H) - initial_energy))))
        row = {"t": float(t), "cells": []}
        for i in range(2):
            mean = expect(st, Q[3 * i])
            variance = expect(st, Q[3 * i] @ Q[3 * i]) - mean**2
            row["cells"].append({
                "mean": mean.tolist(),
                "mean_derivative": expect(st, qdot[i]).tolist(),
                "variance_of_actual_Q": variance.tolist()})
        clock_data.append(row)
    parity1 = np.diag([1.0, -1.0, 1.0])
    parity = tensor([ID3, parity1, parity1, ID3, parity1, parity1])
    parity_error = float(np.max(np.abs(parity @ H @ parity - H)))
    probs = expect(final, effect)
    source_data = {}
    operators = {
        "dN1": source_N[0], "dN2": source_N[1],
        "dlogv1": source_logv[0], "dlogv2": source_logv[1],
        "dlogw": source_logw,
        "dLambda": S[0] + S[1],
    }
    for label, op in operators.items():
        source_data[label] = {
            "initial_mean": expect(initial, op).tolist(),
            "final_mean": expect(final, op).tolist(),
            "final_second_moment": expect(final, op @ op).tolist(),
        }
    # Clock force retains interaction: i[lambda f(QT)PhiPsi, PT].
    clock_force = [1j * (S[i] @ P[3*i] - P[3*i] @ S[i]) / 4.0
                   for i in range(2)]
    result = {
        "status": "FLOAT_DIAGNOSTIC_ONLY_NOT_A_TIME_WINDOW_CERTIFICATE",
        "dimension": 729,
        "parameters": {"w": 0.1, "lambda": 0.25, "t_star": 0.4,
                       "N": [1, 1], "v": [1, 1], "m": [1, 1, 1]},
        "quantization": {
            "quadratics": "compress infinite HO P^2,Q^2 separately",
            "f": "sin of the finite Q matrix",
            "link": "square of the finite two-cell Q difference",
            "Q_norm": qnorm,
            "Q2_minus_Q_squared_diag": np.real(np.diag(q2-q@q)).tolist(),
        },
        "algebra_diagnostics": {
            "H_sum_dNiH_max_error": euler_error,
            "clock_Qdot_equals_P_max_error": qdot_error,
            "joint_Phi_Psi_parity_max_error": parity_error,
            "eigh_residual_max": float(np.max(np.abs(
                H @ eigenvectors-eigenvectors*eigenvalues))),
            "H_min_eigenvalue_diagnostic": float(eigenvalues[0]),
            "H_minus_three_quarters_H0_min_eigenvalue_diagnostic":
                float(np.linalg.eigvalsh(H-0.75*H0)[0]),
        },
        "record": {
            "probability_plus_minus": probs.tolist(),
            "minus_minus_plus": float(probs[1]-probs[0]),
            "effect_is_explicit_terminal_primitive": True,
        },
        "complete_energy": {
            "initial_plus_minus": initial_energy.tolist(),
            "final_plus_minus": expect(final, H).tolist(),
            "sampled_max_conservation_error": max_energy_error,
            "clock_initial_HO_energy_per_cell": 1.0,
            "all_six_oscillators_and_interaction_retained": True,
        },
        "sources": source_data,
        "clock": {
            "sampled_min_mean_derivative": min(
                min(c["mean_derivative"]) for row in clock_data for c in row["cells"]),
            "sampled_max_variance": max(
                max(c["variance_of_actual_Q"]) for row in clock_data for c in row["cells"]),
            "endpoint": [clock_data[0], clock_data[-1]],
            "interaction_force_at_final": [expect(final, f).tolist()
                                           for f in clock_force],
            "not_a_sharp_monotone_clock_or_proper_time_certificate": True,
        },
        "scope": [
            "One explicit two-cell finite effective quantum rule.",
            "No continuum Hilbert stress, Einstein response, full CC/SR or hardware claim.",
            "dNi,dlogvi,dlogw are exact derivatives of the declared Hamiltonian family.",
            "The two preparations have equal source expectations by parity; "
            "their local Psi effect differs.",
            "Sampling and eigh residuals are diagnostics, not rigorous certificates.",
        ],
        "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    assert abs(initial_energy[0]-5.2) < 1e-12
    assert euler_error < 1e-12 and qdot_error < 1e-12 and parity_error < 1e-12
    assert probs[1] > probs[0]
    return result


def compare(a, b, path=""):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for key in a:
            compare(a[key], b[key], path+"/"+key)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+"/"+str(i))
    elif isinstance(a, float):
        assert abs(a-b) <= 2e-11*(1+abs(a)), (path, a, b)
    else:
        assert a == b, (path, a, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("SAVED_EXCLUSIVELY")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf-8")))
        print("READ_ONLY_MATCH_DIAGNOSTIC")
    print(json.dumps({k: result[k] for k in ["record", "complete_energy", "clock"]},
                     ensure_ascii=False))
