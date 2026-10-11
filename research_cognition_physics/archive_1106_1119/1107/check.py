"""Round 1107: exact reset weights and finite checks of the analytic bound.

Default execution is read-only and compares saved results. --write-results
creates results.json once. No optimal-code or physical-realization claim.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCES = [
    "research_cognition_physics/archive_1086_1105/1098/proof.md",
    "research_cognition_physics/archive_1086_1105/1098/results.json",
    "research_cognition_physics/archive_1086_1105/_admission/after_1097_full_influence.md",
    "research_cognition_physics/archive_1009_1043/research_note_1041.md",
    "research_cognition_physics/archive_1009_1043/1041/relay_information_admission.md",
    "research_cognition_physics/认知时空操作公理与洛伦兹运动学_阶段论文.md",
]


def strings(matrix):
    return [[str(value) for value in row] for row in matrix]


def reset_decomposition(matrix):
    rows, columns = len(matrix), len(matrix[0])
    assert all(matrix[y][x] > 0 for y in range(rows) for x in range(columns))
    assert all(sum(matrix[y][x] for y in range(rows)) == 1 for x in range(columns))
    minima = [min(row) for row in matrix]
    alpha = sum(minima)
    reset = [[value / alpha] * columns for value in minima]
    residual = None if alpha == 1 else [
        [(matrix[y][x] - minima[y]) / (1 - alpha) for x in range(columns)]
        for y in range(rows)]
    if residual is not None:
        assert all(value >= 0 for row in residual for value in row)
        assert all(sum(residual[y][x] for y in range(rows)) == 1 for x in range(columns))
        assert all(alpha * reset[y][x] + (1 - alpha) * residual[y][x] == matrix[y][x]
                   for y in range(rows) for x in range(columns))
    else:
        assert reset == matrix
    return alpha, reset, residual


def group_decompositions():
    cases = {
        "bsc_one_tenth": [[F(9, 10), F(1, 10)], [F(1, 10), F(9, 10)]],
        "biased_binary": [[F(4, 5), F(2, 5)], [F(1, 5), F(3, 5)]],
        "three_output": [[F(1, 2), F(1, 4)], [F(1, 3), F(1, 4)], [F(1, 6), F(1, 2)]],
        "constant": [[F(2, 3), F(2, 3)], [F(1, 3), F(1, 3)]],
    }
    answer = {}
    for name, matrix in cases.items():
        alpha, reset, residual = reset_decomposition(matrix)
        answer[name] = dict(T=strings(matrix), alpha=str(alpha), R=strings(reset),
                            residual=None if residual is None else strings(residual))
    assert answer["bsc_one_tenth"]["alpha"] == "1/5"
    return answer


def bit(function, argument):
    return (function >> argument) & 1


def group_feedback():
    # Fixed A->B, B->A, A->B slots. A has fresh bit s; B outputs f(y1,y3).
    # Local deterministic strategies exhaust this finite classical subfamily.
    outputs = list(itertools.product((0, 1), repeat=3))
    count, reset_independent, decomposition_checked = 0, 0, 0
    max_contrast, minimum_error = F(0), F(1)
    for first in range(4):
        for feedback in range(4):
            for last in range(16):
                for decoder in range(16):
                    probabilities, reset_counts = [], []
                    for s in (0, 1):
                        numerator, reset_numerator = 0, 0
                        for y1, y2, y3 in outputs:
                            x1 = bit(first, s)
                            x2 = bit(feedback, y1)
                            x3 = bit(last, 2 * s + y2)
                            out = bit(decoder, 2 * y1 + y3)
                            if out:
                                numerator += ((9 if y1 == x1 else 1)
                                              * (9 if y2 == x2 else 1)
                                              * (9 if y3 == x3 else 1))
                                reset_numerator += 1
                        probabilities.append(F(numerator, 1000))
                        reset_counts.append(reset_numerator)
                    assert reset_counts[0] == reset_counts[1]
                    reset_independent += 1
                    q, reset_probability = F(1, 125), F(reset_counts[0], 8)
                    for probability in probabilities:
                        remainder = (probability - q * reset_probability) / (1 - q)
                        assert 0 <= remainder <= 1
                    decomposition_checked += 1
                    contrast = abs(probabilities[1] - probabilities[0])
                    average_error = (probabilities[0] + 1 - probabilities[1]) / 2
                    assert average_error >= q / 2
                    max_contrast = max(max_contrast, contrast)
                    minimum_error = min(minimum_error, average_error)
                    count += 1
    assert count == 4096
    return dict(protocols=count, all_reset_independent=reset_independent,
                convex_decomposition_checked=decomposition_checked,
                uses_including_return=3, reset_weight="1/125",
                maximum_classical_contrast=str(max_contrast),
                minimum_average_classical_error=str(minimum_error),
                scope="fixed three-slot deterministic classical subfamily only; general quantum theorem is analytic")


def maximally_entangled(d):
    vector = np.eye(d, dtype=complex).reshape(d * d) / np.sqrt(d)
    return np.outer(vector, vector.conj())


def trace_distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a - b)).sum() / 2)


def group_choi_witness():
    rng = np.random.default_rng(1107)
    rows, maximum_residual = [], 0.0
    for d in (2, 3, 4):
        phi = maximally_entangled(d)
        states = [np.eye(d) / d, np.diag([1] + [0] * (d - 1))]
        for _ in range(8):
            matrix = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
            state = matrix @ matrix.conj().T
            states.append(state / np.trace(state))
        q = F(1, 125)
        bound = q * (1 - F(1, d * d))
        for index, tau in enumerate(states):
            replacement = np.kron(np.eye(d) / d, tau)
            overlap = float(np.trace(phi @ replacement).real)
            residual = abs(overlap - 1 / (d * d))
            maximum_residual = max(maximum_residual, residual)
            assert residual < 1e-13
            output = float(q) * replacement + (1 - float(q)) * phi
            bell_deficit = 1 - float(np.trace(phi @ output).real)
            distance = trace_distance(output, phi)
            assert abs(bell_deficit - float(bound)) < 1e-13
            assert distance + 1e-13 >= float(bound)
            if index == 0:
                assert abs(distance - float(bound)) < 1e-13
        rows.append(dict(d=d, sigma_cases=len(states), reset_weight=str(q),
                         exact_effect_deficit=str(bound),
                         mixed_identity_example_saturates_effect_and_trace_bound=True))
    return dict(rows=rows, maximum_overlap_residual=maximum_residual,
                diamond_optimization_performed=False,
                note="maximally entangled input gives a lower bound; all strategies covered by proof")


def group_teleportation_comparison():
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Z = np.diag([1, -1]).astype(complex)
    phi = maximally_entangled(2)
    paulibells = []
    for U in (I, X, Z, X @ Z):
        op = np.kron(I, U)
        paulibells.append(op @ phi @ op.conj().T)
    rows = []
    for repeats, e in [(1, F(1, 10)), (3, F(7, 250))]:
        N = 2 * repeats
        weights = [(1 - e)**2, e * (1 - e), e * (1 - e), e**2]
        choi = sum(float(w) * state for w, state in zip(weights, paulibells))
        strategy_error = 1 - weights[0]
        q = F(1, 5)**N
        bound = q * F(3, 4)
        remainder = (choi - float(q) * np.eye(4) / 4) / (1 - float(q))
        assert np.linalg.eigvalsh(remainder).min() >= -1e-13
        assert abs(np.trace(remainder) - 1) < 1e-13
        assert abs(trace_distance(choi, phi) - float(strategy_error)) < 1e-13
        assert strategy_error >= bound > 0
        rows.append(dict(repeats_per_label=repeats, total_uses=N,
                         label_error=str(e), old_strategy_error=str(strategy_error),
                         all_protocols_reset_bound=str(bound),
                         positive_remainder_choi_verified=True))
    assert rows[1]["old_strategy_error"] == "3451/62500"
    assert rows[1]["all_protocols_reset_bound"] == "3/62500"
    return dict(rows=rows, old_strategy_not_reoptimized=True,
                known_Pauli_half_diamond_formula_reused=True)


def group_shared_noise_control():
    # Identical BSC marginals, but a common Z violates the independent-slot contract.
    kernel = [[F(0) for _ in range(4)] for _ in range(8)]
    for message in range(4):
        b1, b2 = (message >> 1) & 1, message & 1
        for z, probability in ((0, F(9, 10)), (1, F(1, 10))):
            y0, y1, y2 = z, b1 ^ z, b2 ^ z
            assert (y1 ^ y0, y2 ^ y0) == (b1, b2)
            kernel[(y0 << 2) | (y1 << 1) | y2][message] += probability
    assert all(sum(row[column] for row in kernel) == 1 for column in range(4))
    joint_reset_weight = sum(min(row) for row in kernel)
    assert joint_reset_weight == 0
    return dict(marginal_bsc_error="1/10", marginal_reset_weight="1/5",
                actual_joint_kernel=strings(kernel), joint_common_reset_weight="0",
                two_bits_exactly_recovered_with_uses=3,
                independent_uses_assumption_satisfied=False,
                hypothetical_independent_weight="1/125",
                actual_time_and_resources_certified=False)


def group_finite_budget():
    rows = [dict(d=d, N=N, reset_weight=str(F(1, 5)**N),
                 necessary_error_lower_bound=str(F(1, 5)**N * (1 - F(1, d*d))))
            for d in (2, 3, 4) for N in (0, 1, 2, 6, 10)]
    thresholds = []
    for epsilon in (F(1, 100), F(1, 10000)):
        N = 0
        while F(3, 4) * F(1, 5)**N > epsilon:
            N += 1
        thresholds.append(dict(target_error=str(epsilon), first_N_not_excluded_by_bound=N,
                               sufficiency_or_deadline_feasibility_asserted=False))
    constant_rows = []
    for d in (2, 3, 4):
        bounds = [F(1)**N * (1 - F(1, d*d)) for N in (0, 1, 2, 6, 10)]
        assert len(set(bounds)) == 1 and bounds[0] > 0
        constant_rows.append(dict(d=d, alpha="1", all_N_bound=str(bounds[0]),
                                  repetition_cannot_approach_identity=True))
    return dict(rows=rows, thresholds=thresholds, all_finite_N_bound_positive=True,
                constant_full_support_boundary=constant_rows,
                bsc_uniform_positive_bound_over_all_N=False,
                input_after_deadline_or_overhead_ignored=False,
                physical_N_limit_assumed_without_evidence=False)


def compute():
    audit = json.loads((HERE / "source_audit.json").read_text(encoding="utf8"))
    source_hashes = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in SOURCES}
    assert audit["local_source_sha256"] == source_hashes
    return dict(round=1107, status="PASS_BOUNDED_PROTOCOL_OBSTRUCTION",
                groups={"reset_decompositions":group_decompositions(),
                        "finite_feedback_protocols":group_feedback(),
                        "maximally_entangled_witness":group_choi_witness(),
                        "old_Bell_strategy_comparison":group_teleportation_comparison(),
                        "shared_noise_assumption_control":group_shared_noise_control(),
                        "finite_budget_quantifiers":group_finite_budget()},
                source_sha256=source_hashes,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                new_adopted_cognitive_axioms=0, scientific_calibration_increment=0,
                empirical_data=False, full_SR_joint_countermodel=False,
                arbitrary_precision_excluded=False, goal_complete=False,
                general_quantifiers_proved_analytically=True)


def compare(saved, actual):
    if isinstance(actual, dict):
        assert saved.keys() == actual.keys()
        for key in actual:
            compare(saved[key], actual[key])
    elif isinstance(actual, list):
        assert len(saved) == len(actual)
        for old, new in zip(saved, actual):
            compare(old, new)
    elif isinstance(actual, float):
        assert np.isclose(saved, actual, rtol=1e-12, atol=1e-12), (saved, actual)
    else:
        assert saved == actual, (saved, actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    actual = compute()
    if args.write_results:
        with (HERE / "results.json").open("x", encoding="utf8") as stream:
            json.dump(actual, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(json.loads((HERE / "results.json").read_text(encoding="utf8")), actual)
    print(json.dumps(dict(round=1107, status=actual["status"], groups=len(actual["groups"]),
                          saved_results_match=not args.write_results,
                          scientific_calibration_increment=0)))


if __name__ == "__main__":
    main()
