"""Single-copy corrections on the two-dimensional |00>, |11> code."""

import argparse
import json
from pathlib import Path

import numpy as np


def phase_state(unit):
    return np.array([np.sqrt(3.0), unit], dtype=complex) / 2


def density(vector):
    return np.outer(vector, vector.conj())


def trace_distance(left, right):
    return np.abs(np.linalg.eigvalsh(left - right)).sum() / 2


def joint_encoding(first, second, first_target, second_target):
    # Match the full complex Gram matrix to build an isometry on the code.
    overlap = np.vdot(first, second) / np.vdot(first_target, second_target)
    if abs(overlap) > 1 + 1e-13:
        raise ValueError("These pure targets would increase distinguishability")
    env_first = np.array([1.0, 0.0])
    env_second = np.array([overlap, np.sqrt(max(0.0, 1 - abs(overlap) ** 2))])
    inputs = np.column_stack((first, second))
    targets = np.column_stack((np.kron(first_target, env_first),
                               np.kron(second_target, env_second)))
    encoding = targets @ np.linalg.inv(inputs)
    np.testing.assert_allclose(encoding.conj().T @ encoding, np.eye(2), atol=1e-13)
    return encoding, overlap


def data_output(encoding, state):
    amplitudes = (encoding @ state).reshape(2, 2)
    return amplitudes @ amplitudes.conj().T


def run():
    first, second, third = map(phase_state, (1, -1, 1j))
    inputs = np.column_stack((first, second))
    a, b = np.linalg.solve(inputs, third)
    initial_overlap_sq = abs(np.vdot(first, second)) ** 2
    np.testing.assert_allclose([a, b], [(1 + 1j) / 2, (1 - 1j) / 2])

    # Independent linear-algebra evaluation across the two limiting cases
    # (reset, no compression) and the intermediate pure-target choices.
    angles = np.linspace(0, np.pi, 65)
    for angle in angles:
        target = phase_state(np.exp(1j * angle))
        encoding, _ = joint_encoding(first, second, first, target)
        target_overlap_sq = abs(np.vdot(first, target)) ** 2
        np.testing.assert_allclose(data_output(encoding, first), density(first), atol=1e-13)
        np.testing.assert_allclose(data_output(encoding, second), density(target), atol=1e-13)
        reduced = data_output(encoding, third)
        linear_entropy = 1 - np.trace(reduced @ reduced).real
        predicted = 2 * abs(a * b) ** 2 * (1 - target_overlap_sq) * (
            1 - initial_overlap_sq / target_overlap_sq)
        np.testing.assert_allclose(linear_entropy, predicted, atol=1e-13)

    encoding, overlap = joint_encoding(first, second, first, third)
    reduced = data_output(encoding, third)
    purity = np.trace(reduced @ reduced).real
    np.testing.assert_allclose(overlap, (3 - 1j) / 5)
    np.testing.assert_allclose(purity, 71 / 80)
    # Test coherent recovery with all four matrix units, including off-diagonals.
    for vector in np.eye(2):
        for other in np.eye(2):
            unit = np.outer(vector, other)
            stored = encoding @ unit @ encoding.conj().T
            np.testing.assert_allclose(encoding.conj().T @ stored @ encoding,
                                       unit, atol=1e-13)
    initial_distance = trace_distance(density(first), density(second))
    compressed_distance = trace_distance(data_output(encoding, first),
                                         data_output(encoding, second))
    lower_bound = (initial_distance - compressed_distance) / 2
    np.testing.assert_allclose(lower_bound, (np.sqrt(3) / 2 - np.sqrt(3 / 8)) / 2)
    return {
        "schmidt_weights": ["3/4", "1/4"],
        "calibration_phase_inputs_in_pi_units": ["0", "1", "1/2"],
        "first_two_targets_in_pi_units": ["0", "1/2"],
        "environment_overlap_real_imaginary": ["3/5", "-1/5"],
        "third_output_purity": "71/80",
        "third_output_linear_entropy": "9/80",
        "data_only_return_error_lower_bound": round(float(lower_bound), 12),
        "full_memory_roundtrip_matrix_units": 4,
        "verified_overlap_tradeoff_cases": len(angles),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    target_file = Path(__file__).with_name("three_phase_feedback_results.json")
    if args.check:
        if json.loads(target_file.read_text(encoding="utf-8")) != result:
            raise ValueError("Saved results differ from this run")
    else:
        with target_file.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
