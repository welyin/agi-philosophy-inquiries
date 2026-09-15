"""Every nonlocal real two-system flow permits finite readout amplification.

Assumes the signed joint flow, all local rotations, independent ancilla reset,
and the real composite of rounds 47-48. This is not a theorem for all GPTs.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from bipartite_composition import GENERATORS, GENERATOR_LABELS
from certified_intervals import Interval as I, SCALE, ceil_div, pi_interval, sin_interval
from operational_effect_closure import majority_error, majority_certificate, exact_copy_gate
from quantum_interface_audit import (ALPHA, IDENTITY, PAULI_Y, PAULI_Z, PAULI_X,
                                     branch_kraus, rotation_unitary)


PARITY = ((1, 1, 1, 1, 1, 1), (1, 1, 1, 1, -1, -1),
          (1, 1, -1, -1, 1, 1), (1, 1, -1, -1, -1, -1))
FILTER_SIGNS = {"B": (1, 1, -1, -1), "A": (1, -1, 1, -1)}
FAMILY_INDICES = {"B": (2, 3), "A": (4, 5)}
HALF_TURNS = (np.eye(4, dtype=complex), np.kron(-1j * PAULI_Y, IDENTITY),
              np.kron(IDENTITY, -1j * PAULI_Y), np.kron(-1j * PAULI_Y, -1j * PAULI_Y))
EXAMPLE = (Fraction(1, 3), Fraction(-1, 4), Fraction(2, 5), Fraction(-1, 7),
           Fraction(1, 6), Fraction(1, 8))


def exact_coefficients(values):
    if len(values) != 6:
        raise ValueError("Supply six coefficients in YI,IY,YX,YZ,XY,ZY order.")
    if any(not isinstance(x, (int, str, Fraction)) for x in values):
        raise TypeError("Certificates require exact integers or fractions, not floats.")
    return tuple(Fraction(x) for x in values)


def hamiltonian(coefficients):
    return sum(float(x) * p for x, p in zip(coefficients, GENERATORS))


def sectors(coefficients):
    a, b, c, d, e, f = [float(x) for x in coefficients]
    return (a * GENERATORS[0] + e * GENERATORS[4] + f * GENERATORS[5],
            b * GENERATORS[1] + c * GENERATORS[2] + d * GENERATORS[3])


def joint_flow(coefficients, time):
    """Closed Pauli formula; independent eigh exponential is used in tests."""
    result = np.eye(4, dtype=complex)
    for part in sectors(coefficients):
        radius = math.sqrt(max(0., np.trace(part @ part).real / 4))
        if radius:
            result = result @ (math.cos(time * radius / 2) * np.eye(4)
                               - 1j * math.sin(time * radius / 2) * part / radius)
    return result


def family_coefficients(coefficients, input_system):
    if input_system not in FAMILY_INDICES:
        raise ValueError("Input system must be A or B.")
    return tuple(x if i in FAMILY_INDICES[input_system] else 0 for i, x in enumerate(coefficients))


def classification(coefficients):
    values = exact_coefficients(coefficients)
    squares = {role: sum(values[i]**2 for i in indices) for role, indices in FAMILY_INDICES.items()}
    active = [role for role in ("A", "B") if squares[role] > 0]
    return {"active_input_families": active, "generated_lie_dimension": 2 + 2 * len(active),
            "family_norm_squares_exact": {k: str(v) for k, v in squares.items()},
            "both_original_interfaces_stable": not active,
            "scope": "At least one indicated input can be amplified. This does not assert amplification of the opposite input when only one family is active."}


def bracket(left, right):
    # Hermitian-generator bracket corresponding to [-iH/2,-iK/2].
    return -.5j * (left @ right - right @ left)


def bracket_table():
    return [[[int(round(np.trace(p @ bracket(a, b)).real / 4)) for p in GENERATORS]
             for b in GENERATORS] for a in GENERATORS]


def independent_rows(rows):
    basis = []
    for row in rows:
        value = list(map(Fraction, row))
        for pivot, existing in basis:
            factor = value[pivot]
            value = [x - factor * y for x, y in zip(value, existing)]
        pivot = next((i for i, x in enumerate(value) if x), None)
        if pivot is not None:
            divisor = value[pivot]
            basis.append((pivot, [x / divisor for x in value]))
            basis.sort(key=lambda pair: pair[0])
    return [row for _, row in basis]


def exact_lie_rank(coefficients):
    values, table = exact_coefficients(coefficients), bracket_table()
    basis = independent_rows(((1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0), values))
    while True:
        expanded = list(basis)
        for a, b in product(basis, repeat=2):
            expanded.append([sum(a[i] * b[j] * table[i][j][k] for i in range(6) for j in range(6))
                             for k in range(6)])
        following = independent_rows(expanded)
        if len(following) == len(basis):
            return len(basis)
        basis = following


def filtered_flow(coefficients, input_system, time, steps):
    if not isinstance(steps, int) or steps < 1:
        raise ValueError("Use a positive integer number of product-formula steps.")
    if input_system not in FILTER_SIGNS:
        raise ValueError("Input system must be A or B.")
    step = np.eye(4, dtype=complex)
    for q, sign in zip(HALF_TURNS, FILTER_SIGNS[input_system]):
        step = step @ q @ joint_flow(coefficients, sign * time / (4 * steps)) @ q.conj().T
    return np.linalg.matrix_power(step, steps)


def compiled_pointer_gate(coefficients, input_system, steps):
    indices = FAMILY_INDICES[input_system]
    x, z = [float(coefficients[i]) for i in indices]
    radius = math.hypot(x, z)
    if radius == 0:
        raise ValueError("The selected input family is absent.")
    angle = math.atan2(z, x) - math.pi / 2
    if input_system == "B":
        axis = np.kron(IDENTITY, rotation_unitary(angle))
        local = np.kron(rotation_unitary(math.pi / 2), IDENTITY)
    else:
        axis = np.kron(rotation_unitary(angle), IDENTITY)
        local = np.kron(IDENTITY, rotation_unitary(math.pi / 2))
    filtered = filtered_flow(coefficients, input_system, -math.pi / (2 * radius), steps)
    return local @ axis @ filtered @ axis.conj().T


def target_pointer_gate(input_system):
    if input_system == "B":
        return exact_copy_gate()
    if input_system == "A":
        p0, p1 = (IDENTITY + PAULI_Z) / 2, (IDENTITY - PAULI_Z) / 2
        return np.kron(p0, IDENTITY) + np.kron(p1, -1j * PAULI_Y)
    raise ValueError("Input system must be A or B.")


def pointer_branch(density, gate, input_system, outcome):
    p0 = (IDENTITY + PAULI_Z) / 2
    joint = np.kron(p0, density) if input_system == "B" else np.kron(density, p0)
    joint = gate @ joint @ gate.conj().T
    result = np.zeros((4, 4), dtype=complex)
    for k in branch_kraus(0, outcome, 1.):
        extended = np.kron(k, IDENTITY) if input_system == "B" else np.kron(IDENTITY, k)
        result += extended @ joint @ extended.conj().T
    tensor = result.reshape(2, 2, 2, 2)
    return np.einsum("abad->bd", tensor) if input_system == "B" else np.einsum("abcb->ac", tensor)


def majority_probability(density, gate, input_system, cycles=3):
    if cycles < 1 or cycles % 2 == 0:
        raise ValueError("Use a positive odd number of cycles.")
    result = 0.
    for record in product((0, 1), repeat=cycles):
        if sum(record) <= cycles // 2:
            continue
        output = density.copy()
        for outcome in record:
            output = pointer_branch(output, gate, input_system, outcome)
        result += np.trace(output).real
    return result


def majority_effect(gate, input_system, cycles=3):
    pz = majority_probability((IDENTITY + PAULI_Z) / 2, gate, input_system, cycles)
    mz = majority_probability((IDENTITY - PAULI_Z) / 2, gate, input_system, cycles)
    px = majority_probability((IDENTITY + PAULI_X) / 2, gate, input_system, cycles)
    e0, ez = (pz + mz) / 2, (pz - mz) / 2
    return e0 * IDENTITY + ez * PAULI_Z + (px - e0) * PAULI_X


def finite_violation_certificate(coefficients, input_system=None):
    values = exact_coefficients(coefficients)
    norms = {role: sum(values[i]**2 for i in indices) for role, indices in FAMILY_INDICES.items()}
    if input_system is None:
        input_system = max(norms, key=norms.get)
    if input_system not in norms or norms[input_system] == 0:
        raise ValueError("A nonzero interaction family is required.")
    l1 = sum(abs(x) for x in values)
    ratio = l1**2 / norms[input_system]
    alpha = 4 * sin_interval(I.rational(1, 4))
    gain = alpha * (1 - alpha**2) / 4
    gain_lower = Fraction(gain.lo, SCALE)
    # ||compiled C - C|| <= K/n, K=3*pi^2*||H||_1^2/(128*r^2).
    constant = 3 * pi_interval()**2 * ratio / 128
    constant_upper = Fraction(constant.hi, SCALE)
    required = 6 * constant_upper / gain_lower
    steps = ceil_div(required.numerator, required.denominator)
    gate_error = constant_upper / steps
    probability_error = 3 * gate_error
    margin = gain_lower - probability_error
    return {"coefficients_exact": [str(x) for x in values], "input_system": input_system,
            "classification": classification(values), "l1_norm_exact": str(l1),
            "selected_family_norm_squared_exact": str(norms[input_system]),
            "ideal_three_cycle_gain_lower_exact": str(gain_lower),
            "trotter_constant_upper_exact": str(constant_upper),
            "steps_per_compiled_gate": steps, "joint_flow_pulses_per_gate": 4 * steps,
            "ancilla_cycles": 3, "total_joint_flow_pulses": 12 * steps,
            "compiled_gate_error_upper_exact": str(gate_error),
            "full_event_probability_error_upper_exact": str(probability_error),
            "gain_over_original_success_bound_lower_exact": str(margin),
            "gain_over_original_success_bound_lower_diagnostic": float(margin),
            "strict_finite_violation": margin > 0,
            "scope": "Ideal signed control and exact local rotations; analytic approximation bounds, not hardware or floating-point error certification."}


class InteractionNoiseObstructionTests(unittest.TestCase):
    def test_half_turn_parities_and_filters_are_exact(self):
        for q, row in zip(HALF_TURNS, PARITY):
            for p, sign in zip(GENERATORS, row):
                np.testing.assert_array_equal(q @ p @ q.conj().T, sign * p)
        for role in ("A", "B"):
            selected = sum(sign * q @ hamiltonian(EXAMPLE) @ q.conj().T
                           for sign, q in zip(FILTER_SIGNS[role], HALF_TURNS)) / 4
            np.testing.assert_allclose(selected, hamiltonian(family_coefficients(EXAMPLE, role)), atol=2e-16)

    def test_two_three_dimensional_sectors_commute_and_square_to_scalars(self):
        for coefficients in (EXAMPLE, (1, 2, 3, 4, 5, 6), (0, 1, 0, 0, -2, 3)):
            a, b = sectors(coefficients)
            np.testing.assert_allclose(a @ b, b @ a, atol=3e-15)
            for part in (a, b):
                np.testing.assert_allclose(part @ part, np.trace(part @ part).real / 4 * np.eye(4), atol=3e-15)

    def test_signed_closed_flow_matches_independent_spectral_exponential(self):
        for coefficients, time in product((EXAMPLE, (1, 2, 3, 4, 5, 6), (0, 0, 0, 0, 0, 0)), (-1.7, -.03, .4)):
            h = hamiltonian(coefficients)
            eigenvalues, vectors = np.linalg.eigh(h)
            expected = (vectors * np.exp(-.5j * time * eigenvalues)) @ vectors.conj().T
            actual = joint_flow(coefficients, time)
            np.testing.assert_allclose(actual, expected, atol=4e-15)
            np.testing.assert_allclose(joint_flow(coefficients, -time) @ actual, np.eye(4), atol=7e-16)

    def test_integer_brackets_match_the_pauli_matrices(self):
        for i, j in product(range(6), repeat=2):
            np.testing.assert_array_equal(hamiltonian(bracket_table()[i][j]), bracket(GENERATORS[i], GENERATORS[j]))
        for i, j in product((0, 4, 5), (1, 2, 3)):
            np.testing.assert_array_equal(bracket(GENERATORS[i], GENERATORS[j]), np.zeros((4, 4)))

    def test_exact_lie_closure_has_only_dimensions_two_four_or_six(self):
        for c, d, e, f in product((0, 1), repeat=4):
            coefficients = (Fraction(1, 3), Fraction(-2, 7), c, d, e, f)
            self.assertEqual(exact_lie_rank(coefficients), classification(coefficients)["generated_lie_dimension"])
        self.assertEqual(exact_lie_rank((1, 0, Fraction(1, 10**15), 0, 0, 0)), 4)

    def test_filtered_exponentials_converge_with_the_analytic_bound(self):
        l1 = float(sum(abs(x) for x in EXAMPLE))
        for role, time in product(("A", "B"), (-1.2, .8)):
            previous = math.inf
            for steps in (4, 16, 64):
                expected = joint_flow(family_coefficients(EXAMPLE, role), time)
                error = np.linalg.norm(filtered_flow(EXAMPLE, role, time, steps) - expected, ord=2)
                self.assertLessEqual(error, 3 * time**2 * l1**2 / (32 * steps) + 2e-14)
                self.assertLess(error, previous)
                previous = error

    def test_copy_compiler_works_for_both_input_orientations(self):
        for role in ("A", "B"):
            certificate = finite_violation_certificate(EXAMPLE, role)
            gate = compiled_pointer_gate(EXAMPLE, role, certificate["steps_per_compiled_gate"])
            error = np.linalg.norm(gate - target_pointer_gate(role), ord=2)
            self.assertLess(error, float(Fraction(certificate["compiled_gate_error_upper_exact"])))
            np.testing.assert_allclose(gate.conj().T @ gate, np.eye(4), atol=2e-11)

    def test_ideal_pointer_instrument_retains_the_selected_axis(self):
        density = (IDENTITY + .3 * PAULI_Z + .4 * PAULI_X) / 2
        for role, outcome in product(("A", "B"), (0, 1)):
            sign = 2 * outcome - 1
            expected = np.diag([(1 + sign * ALPHA) * density[0, 0] / 2,
                                (1 - sign * ALPHA) * density[1, 1] / 2])
            np.testing.assert_allclose(pointer_branch(density, target_pointer_gate(role), role, outcome), expected, atol=3e-16)

    def test_compiled_complete_three_cycle_effect_obeys_uniform_probability_bound(self):
        error_ideal = majority_error(3, (1 - ALPHA) / 2)
        effect_ideal = (IDENTITY + (1 - 2 * error_ideal) * PAULI_Z) / 2
        for role in ("A", "B"):
            certificate = finite_violation_certificate(EXAMPLE, role)
            gate = compiled_pointer_gate(EXAMPLE, role, certificate["steps_per_compiled_gate"])
            effect = majority_effect(gate, role)
            bound = float(Fraction(certificate["full_event_probability_error_upper_exact"]))
            self.assertLess(np.linalg.norm(effect - effect_ideal, ord=2), bound)
            success = .5 + np.trace(effect @ PAULI_Z).real / 2
            self.assertGreater(success, (1 + ALPHA) / 2)

    def test_rational_certificate_preserves_at_least_half_the_ideal_gain(self):
        old_one = majority_certificate(1)
        old_three = majority_certificate(3)
        for coefficients in (EXAMPLE, (1, 0, Fraction(1, 100), 0, 0, 0), (0, 2, 0, 0, -3, 1)):
            certificate = finite_violation_certificate(coefficients)
            gain = Fraction(certificate["ideal_three_cycle_gain_lower_exact"])
            independent_upper = Fraction(old_one["error_upper_exact"]) - Fraction(old_three["error_lower_exact"])
            independent_lower = Fraction(old_one["error_lower_exact"]) - Fraction(old_three["error_upper_exact"])
            self.assertLess(abs(gain - (independent_upper + independent_lower) / 2), Fraction(1, 10**26))
            self.assertGreaterEqual(Fraction(certificate["gain_over_original_success_bound_lower_exact"]), gain / 2)
            self.assertTrue(certificate["strict_finite_violation"])

    def test_arbitrarily_small_interactions_are_not_mistaken_for_zero(self):
        weak = (1, 0, Fraction(1, 10**12), 0, 0, 0)
        certificate = finite_violation_certificate(weak)
        self.assertGreater(certificate["steps_per_compiled_gate"], 10**24)
        self.assertEqual(certificate["classification"]["generated_lie_dimension"], 4)
        self.assertTrue(certificate["strict_finite_violation"])
        rescaled = tuple(17 * x for x in EXAMPLE)
        self.assertEqual(finite_violation_certificate(EXAMPLE)["steps_per_compiled_gate"],
                         finite_violation_certificate(rescaled)["steps_per_compiled_gate"])

    def test_local_only_and_absent_families_are_rejected(self):
        self.assertTrue(classification((1, 2, 0, 0, 0, 0))["both_original_interfaces_stable"])
        with self.assertRaises(ValueError):
            finite_violation_certificate((1, 2, 0, 0, 0, 0))
        with self.assertRaises(ValueError):
            finite_violation_certificate((0, 0, 1, 0, 0, 0), "A")
        with self.assertRaises(TypeError):
            finite_violation_certificate((.1, 0, 1, 0, 0, 0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(InteractionNoiseObstructionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    certificate = finite_violation_certificate(EXAMPLE)
    gate = compiled_pointer_gate(EXAMPLE, certificate["input_system"], certificate["steps_per_compiled_gate"])
    effect = majority_effect(gate, certificate["input_system"])
    report = {"coefficient_order": list(GENERATOR_LABELS), "half_turn_parity_table": PARITY,
              "family_filter_signs_before_dividing_by_four": FILTER_SIGNS,
              "lie_dimensions": {"no_active_family": 2, "one_active_family": 4, "two_active_families": 6},
              "every_nonzero_interaction_has_finite_noise_bound_violation": True,
              "necessary_scope": ["Real two-rebit composite", "All local SO(2) controls", "Signed flow exp(-itH/2) for arbitrary real t", "Independent pure ancilla reset and original fresh readout", "The active family determines which system serves as input"],
              "universal_all_GPT_no_interaction_theorem": False,
              "quantum_theory_derived_from_cognition": False,
              "exact_finite_target_gate_synthesis_proven": False,
              "finite_error_formula": "||C_n-C|| <= 3*pi^2*L^2/(128*r^2*n)",
              "effect_norm_closure_on_an_active_input": "All real 2x2 effects, using finite protocols of increasing length",
              "finite_example": certificate,
              "finite_example_floating_diagnostics": {
                  "gate_error": float(np.linalg.norm(gate - target_pointer_gate(certificate["input_system"]), ord=2)),
                  "three_cycle_discrimination_success": .5 + float(np.trace(effect @ PAULI_Z).real / 2),
                  "original_success_upper_bound": (1 + ALPHA) / 2},
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("interaction_noise_obstruction_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
