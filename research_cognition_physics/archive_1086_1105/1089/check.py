"""Round 1089: fixed-time parallel signaling and explicitly scoped budgets.

Only small full Hamiltonians are allocated. Large-copy claims use exact
integer inequalities and analytic tensor product formulas, not simulation
of thousands of qubits. Default operation is read-only verification.
"""
from pathlib import Path
from fractions import Fraction
import json
import math
import sys
import numpy as np

HERE = Path(__file__).resolve().parent


def kron_all(items):
    out = np.array([[1.0]], dtype=complex)
    for item in items:
        out = np.kron(out, item)
    return out


def trace_distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a - b))) / 2)


def receiver_state(global_state, copies):
    # Pair ordering S0,R0,S1,R1,...; move all S factors before all R factors.
    psi = global_state.reshape([2] * (2 * copies))
    psi = np.transpose(psi, list(range(0, 2 * copies, 2)) +
                       list(range(1, 2 * copies, 2)))
    psi = psi.reshape(2**copies, 2**copies)
    return psi.T @ psi.conj()


def run():
    exchange = np.zeros((4, 4), dtype=complex)
    exchange[1, 2] = exchange[2, 1] = 1
    qubit0 = np.array([1, 0], dtype=complex)
    qubit1 = np.array([0, 1], dtype=complex)
    matrix_rows = []
    residuals = []
    for angles in ([.07], [.11, .23], [.08, .13, .19]):
        n = len(angles)
        h = np.zeros((4**n, 4**n), dtype=complex)
        for i, a in enumerate(angles):
            h += a * kron_all([exchange if j == i else np.eye(4)
                               for j in range(n)])
        e, v = np.linalg.eigh(h)
        u = (v * np.exp(-1j * e)) @ v.conj().T
        residuals.append(float(np.max(np.abs(u.conj().T @ u - np.eye(4**n)))))
        states = []
        probabilities = []
        for bit in (qubit0, qubit1):
            psi = kron_all([np.kron(bit, qubit0)[:, None]] * n).ravel()
            rho = receiver_state(u @ psi, n)
            states.append(rho)
            probabilities.append(float(1 - rho[0, 0].real))
        predicted = 1 - math.prod(math.cos(a)**2 for a in angles)
        delta = trace_distance(*states)
        error = max(abs(delta - predicted),
                    abs(probabilities[1] - probabilities[0] - predicted))
        residuals.append(error)
        # Exact commuting-tensor norm sum is checked independently by spectrum.
        residuals.append(abs(float(max(abs(e))) - sum(angles)))
        v0 = np.zeros((4**n, 2**n), dtype=complex)
        for source_index in range(2**n):
            bits = [(source_index >> (n-1-i)) & 1 for i in range(n)]
            output_index = sum(bit * 2 * 4**(n-1-i) for i, bit in enumerate(bits))
            v0[output_index, source_index] = 1
        dilation = u @ v0
        compression = kron_all([np.diag([1, math.cos(a)]) for a in angles])
        residuals.append(float(np.max(np.abs(v0.conj().T @ dilation - compression))))
        dilation_distance = float(np.linalg.norm(dilation - v0, 2))
        distance_formula = math.sqrt(2 * (1 - math.prod(math.cos(a) for a in angles)))
        residuals.append(abs(dilation_distance - distance_formula))
        assert dilation_distance <= math.sqrt(sum(a*a for a in angles)) + 1e-13
        # Maximally entangled unknown source/reference probe; receiver retained.
        psi_ref = dilation.reshape([2] * (2*n) + [2**n]) / math.sqrt(2**n)
        psi_ref = np.transpose(psi_ref, list(range(0, 2*n, 2)) +
                               list(range(1, 2*n, 2)) + [2*n]).reshape(2**n, 4**n)
        rho_ref = psi_ref.T @ psi_ref.conj()
        zero_receiver = np.zeros((2**n, 2**n), dtype=complex)
        zero_receiver[0, 0] = 1
        ideal_ref = np.kron(zero_receiver, np.eye(2**n) / 2**n)
        reference_probe_distance = trace_distance(rho_ref, ideal_ref)
        assert reference_probe_distance <= min(1, distance_formula) + 1e-13
        matrix_rows.append({"copies": n, "hilbert_dimension": 4**n,
                            "angles": angles,
                            "contrast": round(delta, 13),
                            "formula": round(predicted, 13),
                            "dilation_operator_distance": round(dilation_distance, 13),
                            "entangled_reference_probe_distance": round(reference_probe_distance, 13),
                            "all_input_half_diamond_upper": round(math.sqrt(sum(a*a for a in angles)), 13)})

    # General-encoding caveat: same physical line, coherent source pair.
    angle = .3
    e, v = np.linalg.eigh(angle * exchange)
    u = (v * np.exp(-1j * e)) @ v.conj().T
    states = []
    for sign in (1, -1):
        source = (qubit0 + sign * qubit1) / math.sqrt(2)
        states.append(receiver_state(u @ np.kron(source, qubit0), 1))
    coherent = trace_distance(*states)
    residuals.append(abs(coherent - math.sin(angle)))
    assert coherent > math.sin(angle)**2

    # p=1/1000, contrast >=9/10. Exact integer comparisons certify minimal N.
    n = 2302
    assert 10 * 999**n <= 1000**n
    assert 10 * 999**(n-1) > 1000**(n-1)
    delta_n = -math.expm1(n * math.log1p(-.001))

    # No sampling proof: the all-partition upper bound is analytic in proof.md.
    budget_rows = []
    partition_count = 0
    for a in (.1, .4, 1.0, 1.5):
        optimum = math.sin(a)**2
        for m in (1, 2, 3, 8, 32):
            equal = -math.expm1(2 * m * math.log(math.cos(a / m)))
            assert equal <= optimum + 1e-13
            budget_rows.append({"total_action": a, "copies": m,
                                "equal_allocation_contrast": round(equal, 13),
                                "concentrated_optimum_for_basis_code": round(optimum, 13),
                                "all_input_half_diamond_upper_equal_allocation": round(min(1, a/math.sqrt(m)), 13)})
        for i in range(11):
            for j in range(11-i):
                parts = [a*i/10, a*j/10, a*(10-i-j)/10]
                contrast = 1 - math.prod(math.cos(x)**2 for x in parts)
                assert contrast <= optimum + 1e-13
                partition_count += 1

    scale_rows = []
    for k in (10, 30, 100, 300, 1000):
        p = 1 / k**2
        scale_rows.append({"k": k, "single_probability": p,
                           "copies": k**2,
                           "contrast": round(-math.expm1(k**2 * math.log1p(-p)), 13),
                           "total_action": round(k**2 * math.asin(1/k), 10)})
    # A finite transport witness has a strict probability-error margin.
    d, eta = Fraction(9, 10), Fraction(1, 100)
    transported_floor = d - 2 * eta
    assert transported_floor == Fraction(22, 25)
    assert max(residuals) < 1e-12
    return {"schema": "round1089_resource_indexed_signaling_v1",
            "passed": True, "fixed_time": 1,
            "scope": "finite quantum signaling protocol and resource quantifiers; not a spacetime model",
            "small_full_matrix_checks": matrix_rows,
            "full_matrix_residual_below_1e_12": True,
            "large_copy_exact_certificate": {"single_p": "1/1000", "target_contrast": "9/10",
                                            "minimal_copies": n, "qubits": 2*n,
                                            "contrast": round(delta_n, 13),
                                            "both_integer_inequalities_verified": True},
            "fixed_budget": budget_rows,
            "three_part_budget_samples": partition_count,
            "coherent_encoding_warning": {"angle": angle,
                                          "basis_contrast": round(math.sin(angle)**2, 13),
                                          "coherent_contrast": round(coherent, 13),
                                          "basis_bound_is_not_general_capacity_bound": True},
            "nonuniform_scale_limit": scale_rows,
            "limit_when_Np_tends_to_1": round(1-math.exp(-1), 13),
            "transport_probability_floor": str(transported_floor),
            "joint_FUCP_CO_model_certified": False,
            "Lorentz_proved_or_disproved": False}


if __name__ == "__main__":
    result = run()
    if "--write" in sys.argv:
        (HERE / "results.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    else:
        assert result == json.loads((HERE / "results.json").read_text(encoding="utf8"))
    print(json.dumps({"round": 1089, "passed": True,
                      "mode": "write" if "--write" in sys.argv else "read_only_compare"}))
