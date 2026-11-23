"""Distributed heralded phase increments, including dephasing on the flag link."""

import argparse
import json
import math
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from joint_action import joint_window


def previous_window_truth_table():
    # Read the toggle into a clean flag and return all three old lanes. The
    # coherent extension in the note is an additional physical assumption.
    checked = 0
    for bits in product((0, 1), repeat=8):
        a, b, *state = bits
        lanes = tuple(zip(state[::2], state[1::2]))
        pulses = (a, b, a, b)
        after = joint_window(lanes, pulses)
        flag = lanes[0][1] ^ after[0][1]
        assert flag == a*b
        assert joint_window(after, pulses) == lanes
        checked += 1
    return checked

def branch_operators(m, sign, visibility):
    """Enumerate control paths by Walsh signs and keep the two flag-environment amplitudes."""
    size = 2**m
    marked = size - 1
    data_mark = np.diag([-1j * sign, 1j * sign])
    branch = []
    for syndrome in range(size):
        signs = [(-1) ** ((syndrome & x).bit_count()) for x in range(size)]
        plain = (sum(signs[:marked]) / size) * np.eye(2)
        marked_path = (signs[marked] / size) * data_mark
        branch.append((plain + visibility * marked_path,
                       math.sqrt(1 - visibility**2) * marked_path))
    return branch


def branch_checks():
    unit = np.array([[0, 0], [1, 0]], dtype=complex)
    checked = 0
    for m in range(2, 8):
        p = 2.0**(-m)
        success = (1 - p)**2 + p*p
        beta = 2 * math.atan2(p, 1 - p)
        for sign in (-1, 1):
            undo_failure = np.diag([1, 1j * sign])
            for visibility in (0, 0.35, 0.8, 1):
                total_effect = np.zeros((2, 2), dtype=complex)
                for w, operators in enumerate(branch_operators(m, sign, visibility)):
                    probability = success if w == 0 else 2 * p*p
                    effect = sum(k.conj().T @ k for k in operators)
                    np.testing.assert_allclose(effect, probability * np.eye(2), atol=1e-13)
                    total_effect += effect
                    cell = sum(k @ unit @ k.conj().T for k in operators) / probability
                    if w:
                        cell = undo_failure @ cell @ undo_failure.conj().T
                        expected = visibility
                    else:
                        expected = math.cos(beta) + 1j * sign * visibility * math.sin(beta)
                    np.testing.assert_allclose(cell, expected * unit, atol=1e-13)
                    # Diagonal Kraus operators: diagonal cells, this unit, and its
                    # adjoint determine all outputs even with an arbitrary partner.
                    checked += 1
                np.testing.assert_allclose(total_effect, np.eye(2), atol=1e-13)
    return checked


def collinear_limit():
    """Run one radian of nominal turning, averaging completed applications."""
    rows = []
    for visibility in (0, 0.5, 0.9, 1):
        limit = np.exp(-(1 - visibility) + 1j * visibility)
        errors = []
        for m in (4, 8, 12, 16):
            p = 2.0**(-m)
            fail = 2 * p * (1 - p)
            completed = (1 - 2*p + 1j * visibility * fail) / (1 - visibility * fail)
            repeated = completed ** (2**(m-1))
            errors.append(float(abs(repeated - limit)))
        assert all(a > b for a, b in zip(errors, errors[1:]))
        assert errors[-1] < 2e-5
        rows.append({"link_visibility": visibility, "limit_modulus": round(abs(limit), 9),
                     "limit_angle_radians": visibility,
                     "absolute_complex_errors": [round(e, 11) for e in errors]})
    return rows


def run():
    rows = []
    for m in (2, 4, 8, 12):
        p = Fraction(1, 2**m)
        success = (1-p)**2 + p**2
        tangent = p / (1-p)
        cos_beta = (1 - tangent**2) / (1 + tangent**2)
        sin_beta = 2*tangent / (1 + tangent**2)
        assert cos_beta**2 + sin_beta**2 == 1
        assert (2**m - 1) * 2 * p**2 == 1 - success
        if m == 2:
            assert (cos_beta, sin_beta) == (Fraction(4, 5), Fraction(3, 5))
        rows.append({"controls": m, "success_probability": str(success),
                     "rotation_radians": round(2 * math.atan(tangent), 10),
                     "expected_attempts": str(1 / success),
                     "failure_after_three_attempts": str((1 - success)**3)})
    return {"previous_window_truth_table_checks": previous_window_truth_table(),
            "individual_syndrome_checks": branch_checks(), "resource_examples": rows,
            "one_nominal_radian": collinear_limit()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    output = Path(__file__).with_name("coherent_joint_increment_results.json")
    if args.check:
        assert json.loads(output.read_text(encoding="utf-8")) == result
    else:
        with output.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
