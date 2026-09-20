"""Role covariance, real control completion, and assisted state exchange.

The real matrix framework remains a candidate assumption. Relabeling a model,
covariance of its capabilities, and executing SWAP are separate requirements.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bipartite_composition import (GENERATORS, GENERATOR_LABELS, HIDDEN,
                                   FULL_BASIS, integer_determinant, visible_summary)
from interaction_noise_obstruction import (EXAMPLE, bracket_table, exact_coefficients,
                                           hamiltonian, independent_rows, joint_flow)
from quantum_interface_audit import (ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                     effective_readout, rotation_unitary, two_rebit_states)


SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=int)
ROLE_PERMUTATION = (1, 0, 4, 5, 2, 3)
PAULIS = {"I": IDENTITY, "X": PAULI_X, "Y": PAULI_Y, "Z": PAULI_Z}
SINGLET = (np.eye(4) - SWAP) / 2
# W source W^T = sign * target. Each W uses two real two-system pi/2 pulses.
ROUTES = {"IIY": ("IIY", 1, ()),
          "XXY": ("YII", 1, ("XIY", "YXI")),
          "YYY": ("YII", -1, ("XIY", "XYI")),
          "ZZY": ("YII", 1, ("XYI", "IXY"))}
TARGET_TERMS = (("IIY", 1), ("XXY", -1), ("YYY", -1), ("ZZY", -1))


def swap_coefficients(coefficients):
    values = exact_coefficients(coefficients)
    return tuple(values[i] for i in ROLE_PERMUTATION)


def generated_lie_rank(generators):
    table = bracket_table()
    basis = independent_rows(((1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0), *generators))
    while True:
        candidates = list(basis)
        for a, b in product(basis, repeat=2):
            candidates.append([sum(a[i] * b[j] * table[i][j][k] for i in range(6) for j in range(6))
                               for k in range(6)])
        following = independent_rows(candidates)
        if len(following) == len(basis):
            return len(basis)
        basis = following


def role_closure_report(coefficients):
    values = exact_coefficients(coefficients)
    nonlocal_present = any(values[2:])
    return {"coefficients_exact": [str(x) for x in values],
            "role_reversed_coefficients_exact": [str(x) for x in swap_coefficients(values)],
            "rank_before_role_closure": generated_lie_rank([values]),
            "rank_after_role_closure": generated_lie_rank([values, swap_coefficients(values)]),
            "analytic_role_closed_rank": 6 if nonlocal_present else 2}


def proper_rotation_to(vector):
    vector = np.asarray(vector, dtype=float)
    vector = vector / np.linalg.norm(vector)
    first = np.eye(len(vector))[0]
    difference = first - vector
    if np.linalg.norm(difference) < 1e-14:
        return np.eye(len(vector))
    axis = difference / np.linalg.norm(difference)
    reflection = np.eye(len(vector)) - 2 * np.outer(axis, axis)
    correction = np.eye(len(vector))
    correction[1, 1] = -1
    # The product is proper and sends e0 to vector. The individual reflections
    # are algebraic construction devices, not assumed available controls.
    return reflection @ correction


def swap_separating_state(proper_gate):
    relative = SWAP.T @ proper_gate
    _, _, plus_rows = np.linalg.svd(relative - np.eye(4))
    _, _, minus_rows = np.linalg.svd(relative + np.eye(4))
    plus, minus = plus_rows[-1], minus_rows[-1]
    minus = minus - np.dot(plus, minus) * plus
    minus /= np.linalg.norm(minus)
    return (plus + minus) / math.sqrt(2)


def readout_alignment(first, second):
    basis = []
    for vector in (first, second, *np.eye(4)):
        current = np.asarray(vector, dtype=float).copy()
        for previous in basis:
            current -= np.dot(previous, current) * previous
        if np.linalg.norm(current) > 1e-12:
            basis.append(current / np.linalg.norm(current))
        if len(basis) == 4:
            break
    frame = np.column_stack(basis)
    if np.linalg.det(frame) < 0:
        frame[:, -1] *= -1
    return frame.T


def pauli_word(word):
    result = np.ones((1, 1), dtype=complex)
    for letter in word:
        result = np.kron(result, PAULIS[letter])
    return result


def multiply_words(left, right):
    table = {("X", "Y"): (1j, "Z"), ("Y", "X"): (-1j, "Z"),
             ("Y", "Z"): (1j, "X"), ("Z", "Y"): (-1j, "X"),
             ("Z", "X"): (1j, "Y"), ("X", "Z"): (-1j, "Y")}
    phase, letters = 1, []
    for a, b in zip(left, right):
        if a == "I":
            letters.append(b)
        elif b == "I":
            letters.append(a)
        elif a == b:
            letters.append("I")
        else:
            value, letter = table[a, b]
            phase *= value
            letters.append(letter)
    return phase, "".join(letters)


def check_route(target):
    source, expected_sign, route = ROUTES[target]
    sign, current = 1, source
    for generator in route:
        phase, current = multiply_words(generator, current)
        if phase not in (1j, -1j):
            raise ValueError("The route must use anticommuting Pauli words.")
        sign *= int((-1j * phase).real)
    return current == target and sign == expected_sign


def pauli_rotation(word, angle):
    matrix = pauli_word(word)
    return math.cos(angle / 2) * np.eye(len(matrix)) - 1j * math.sin(angle / 2) * matrix


def assisted_swap_pulses(parameter):
    """Chronological pulses, each with at most two nonidentity factors."""
    pulses = []
    for target, coefficient in TARGET_TERMS:
        source, sign, route = ROUTES[target]
        pulses.extend((word, -math.pi / 2) for word in reversed(route))
        pulses.append((source, sign * coefficient * parameter / 2))
        pulses.extend((word, math.pi / 2) for word in route)
    return pulses


def compiled_assisted_swap(parameter):
    result = np.eye(8, dtype=complex)
    for word, angle in assisted_swap_pulses(parameter):
        result = pauli_rotation(word, angle) @ result
    return result


def explicit_assisted_path(parameter):
    return np.kron(np.eye(4) - SINGLET, IDENTITY) + np.kron(SINGLET, rotation_unitary(2 * parameter))


def exact_assisted_swap_certificate():
    # At parameter pi all sixteen angles are +/-pi/2. Each numerator
    # I - i*sign*P is an integer real matrix, with common divisor sqrt(2).
    pulses = []
    for target, coefficient in TARGET_TERMS:
        source, sign, route = ROUTES[target]
        pulses.extend((word, -1) for word in reversed(route))
        pulses.append((source, sign * coefficient))
        pulses.extend((word, 1) for word in route)
    numerator = np.eye(8, dtype=np.int64)
    for word, sign in pulses:
        raw = np.eye(8) - 1j * sign * pauli_word(word)
        integer = raw.real.astype(np.int64)
        if np.any(raw.imag != 0) or not np.array_equal(raw.real, integer):
            raise AssertionError("The pulse numerator must be an exact real integer matrix.")
        numerator = integer @ numerator
    denominator = 2**(len(pulses) // 2)
    target = denominator * np.kron(SWAP, np.eye(2, dtype=np.int64))
    verified = len(pulses) % 2 == 0 and np.array_equal(numerator, target)
    return {"verified": bool(verified), "pulse_count": len(pulses),
            "two_system_pulse_count": sum(sum(letter != "I" for letter in word) == 2 for word, _ in pulses),
            "single_system_pulse_count": sum(sum(letter != "I" for letter in word) == 1 for word, _ in pulses),
            "pulse_angle_unit": "pi/2", "chronological_pulses": [[word, sign] for word, sign in pulses],
            "common_denominator": denominator, "integer_product": numerator.tolist(),
            "target": "SWAP_AB tensor I_C", "auxiliary_returned_for_all_joint_inputs": True,
            "scope": "Exact relative to available real single/pair rotations on all three pairs; an arbitrary primitive H may only approximate those pair rotations."}


class RoleSymmetryAndSwapTests(unittest.TestCase):
    def test_role_permutation_matches_all_six_generators_exactly(self):
        for i, generator in enumerate(GENERATORS):
            expected = GENERATORS[ROLE_PERMUTATION[i]]
            np.testing.assert_array_equal(SWAP @ generator @ SWAP.T, expected)
        np.testing.assert_array_equal(SWAP @ HIDDEN @ SWAP.T, HIDDEN)

    def test_role_covariance_is_a_property_of_the_gate_set(self):
        for time in (-.4, .8, 1.7):
            actual = SWAP @ joint_flow(EXAMPLE, time) @ SWAP.T
            expected = joint_flow(swap_coefficients(EXAMPLE), time)
            np.testing.assert_allclose(actual, expected, atol=5e-16)
        self.assertNotEqual(tuple(EXAMPLE), swap_coefficients(EXAMPLE))

    def test_all_nonlocal_patterns_become_six_dimensional_under_role_closure(self):
        for coefficients in product((0, 1), repeat=4):
            report = role_closure_report((Fraction(1, 3), Fraction(-1, 7), *coefficients))
            self.assertEqual(report["rank_after_role_closure"], report["analytic_role_closed_rank"])
        report = role_closure_report((0, 0, 1, 0, 0, 0))
        self.assertEqual(report["rank_before_role_closure"], 4)
        self.assertEqual(report["rank_after_role_closure"], 6)

    def test_role_symmetric_local_control_does_not_require_interaction(self):
        report = role_closure_report((1, 2, 0, 0, 0, 0))
        self.assertEqual(report["rank_after_role_closure"], 2)
        coefficients = (1, 1, 2, 3, 2, 3)
        self.assertEqual(tuple(coefficients), swap_coefficients(coefficients))
        self.assertEqual(role_closure_report(coefficients)["rank_before_role_closure"], 6)

    def test_proper_real_rotations_prepare_arbitrary_pure_states(self):
        rng = np.random.default_rng(50)
        for vector in (*np.eye(4), -np.eye(4)[0], *rng.normal(size=(15, 4))):
            vector = vector / np.linalg.norm(vector)
            gate = proper_rotation_to(vector)
            np.testing.assert_allclose(gate[:, 0], vector, atol=8e-16)
            np.testing.assert_allclose(gate.T @ gate, np.eye(4), atol=2e-15)
            self.assertAlmostEqual(np.linalg.det(gate), 1., places=14)

    def test_mixing_proper_rotation_preparations_fills_the_real_positive_cone(self):
        rng = np.random.default_rng(150)
        for _ in range(8):
            matrix = rng.normal(size=(4, 4))
            density = matrix @ matrix.T
            density /= np.trace(density)
            weights, vectors = np.linalg.eigh(density)
            reconstructed = sum(weight * np.outer(proper_rotation_to(vector)[:, 0], proper_rotation_to(vector)[:, 0])
                                for weight, vector in zip(weights, vectors.T))
            np.testing.assert_allclose(reconstructed, density, atol=1e-15)

    def test_bare_swap_has_the_wrong_orientation_even_up_to_global_sign(self):
        self.assertEqual(integer_determinant(SWAP.tolist()), -1)
        self.assertEqual(integer_determinant((-SWAP).tolist()), -1)
        np.testing.assert_array_equal(SWAP.T @ SWAP, np.eye(4, dtype=int))
        for index in range(4):
            ket = np.eye(4)[index]
            np.testing.assert_array_equal(SWAP @ ket, np.eye(4)[(index % 2) * 2 + index // 2])

    def test_every_sample_proper_gate_has_maximal_uniform_swap_channel_error(self):
        rng = np.random.default_rng(250)
        gates = [np.eye(4), joint_flow(EXAMPLE, .73).real]
        for _ in range(12):
            gate, _ = np.linalg.qr(rng.normal(size=(4, 4)))
            if np.linalg.det(gate) < 0:
                gate[:, -1] *= -1
            gates.append(gate)
        readout = np.kron(IDENTITY, effective_readout(0, 1, 1.))
        for gate in gates:
            self.assertAlmostEqual(np.linalg.norm(gate - SWAP, ord=2), 2., places=14)
            self.assertAlmostEqual(np.linalg.norm(gate + SWAP, ord=2), 2., places=14)
            state = swap_separating_state(gate)
            first, second = SWAP @ state, gate @ state
            self.assertAlmostEqual(np.dot(first, second), 0., places=13)
            alignment = readout_alignment(first, second)
            first, second = alignment @ first, alignment @ second
            difference = (first @ readout @ first - second @ readout @ second).real
            self.assertAlmostEqual(abs(difference), ALPHA, places=13)

    def test_exact_word_routes_use_only_allowed_single_or_pair_generators(self):
        for target, (source, sign, route) in ROUTES.items():
            self.assertTrue(check_route(target))
            conjugator = np.eye(8, dtype=complex)
            for word in route:
                self.assertEqual(word.count("Y"), 1)
                self.assertLessEqual(sum(letter != "I" for letter in word), 2)
                conjugator = pauli_rotation(word, math.pi / 2) @ conjugator
            np.testing.assert_allclose(conjugator @ pauli_word(source) @ conjugator.conj().T,
                                       sign * pauli_word(target), atol=5e-16)

    def test_assisted_path_starts_at_identity_and_ends_at_swap_times_identity(self):
        np.testing.assert_allclose(explicit_assisted_path(0), np.eye(8), atol=0)
        np.testing.assert_allclose(explicit_assisted_path(math.pi), np.kron(SWAP, IDENTITY), atol=2e-16)
        for parameter in (.2, .7, 1.8):
            path = explicit_assisted_path(parameter)
            np.testing.assert_allclose(path.conj().T @ path, np.eye(8), atol=3e-16)
            self.assertAlmostEqual(np.linalg.det(path).real, 1., places=14)

    def test_pair_pulse_circuit_reproduces_the_entire_continuous_assisted_path(self):
        for parameter in (0., -.3, .7, math.pi / 2, math.pi):
            compiled = compiled_assisted_swap(parameter)
            np.testing.assert_allclose(compiled, explicit_assisted_path(parameter), atol=2e-15)
        self.assertEqual(len(assisted_swap_pulses(math.pi)), 16)

    def test_sixteen_pulse_swap_has_an_exact_integer_certificate(self):
        report = exact_assisted_swap_certificate()
        self.assertTrue(report["verified"])
        self.assertEqual(report["common_denominator"], 256)
        self.assertEqual(report["two_system_pulse_count"], 12)
        self.assertEqual(report["single_system_pulse_count"], 4)

    def test_auxiliary_is_returned_even_for_correlated_joint_inputs(self):
        rng = np.random.default_rng(350)
        unitary = compiled_assisted_swap(math.pi)
        ideal = np.kron(SWAP, IDENTITY)
        for _ in range(6):
            matrix = rng.normal(size=(8, 8))
            density = matrix @ matrix.T
            density /= np.trace(density)
            actual = unitary @ density @ unitary.conj().T
            expected = ideal @ density @ ideal.T
            np.testing.assert_allclose(actual, expected, atol=7e-16)
            before_aux = np.einsum("abad->bd", density.reshape(4, 2, 4, 2))
            after_aux = np.einsum("abad->bd", actual.reshape(4, 2, 4, 2))
            np.testing.assert_allclose(before_aux, after_aux, atol=1e-15)

    def test_role_symmetry_and_assisted_swap_leave_the_local_tomography_defect(self):
        plus, minus = two_rebit_states()
        np.testing.assert_allclose(visible_summary(plus), visible_summary(minus), atol=0)
        for density in (plus, minus):
            np.testing.assert_allclose(SWAP @ density @ SWAP, density, atol=0)
        for a, b in product((IDENTITY, PAULI_X, PAULI_Z), repeat=2):
            self.assertEqual(np.trace((plus - minus) @ np.kron(a, b)).real, 0.)
        self.assertEqual(np.trace((plus - minus) @ HIDDEN).real, 2.)
        self.assertEqual(len(FULL_BASIS), 10)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(RoleSymmetryAndSwapTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"generator_order": list(GENERATOR_LABELS),
              "generator_permutation_under_role_exchange": ROLE_PERMUTATION,
              "role_covariance_condition": "S G S = G; availability symmetry within the same model, not merely passive renaming",
              "nonlocal_role_covariant_connected_control_closure": "SO(4), within the real two-system framework",
              "local_role_covariant_alternative": "SO(2) tensor SO(2)",
              "single_family_example": role_closure_report((0, 0, 1, 0, 0, 0)),
              "bare_swap_determinant_exact": integer_determinant(SWAP.tolist()),
              "bare_continuous_swap_on_two_real_systems": False,
              "bare_uniform_channel_swap_error_worst_case_trace_distance": 1,
              "bare_swap_distinction_with_aligned_original_readout_TV": ALPHA,
              "assisted_swap_certificate": exact_assisted_swap_certificate(),
              "additional_assisted_scope": "Three-system real tensor composition; same full real pair controls on AB, AC, BC; no auxiliary reset needed for this circuit",
              "role_covariance_or_assisted_swap_implies_local_tomography": False,
              "joint_positive_cone_generated_in_real_framework": "All real symmetric 4x4 PSD matrices, using proper rotations and classical mixing",
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("role_symmetry_and_swap_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
