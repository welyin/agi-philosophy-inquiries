"""Actual readout effects, the first-noise bound, and conditional amplification.

Ancilla amplification uses the two-system gate declared in round 47 and an
independent ancilla reset. A separate many-pointer construction is also checked.
None of these candidate additions is derived from cognition.
See research_note_48.md for analytic proofs and the meaning of state dimension.
"""

import argparse
import json
import math
import unittest
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from certified_intervals import Interval as I, SCALE, ceil_div, sin_interval
from quantum_interface_audit import (ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                     apply_branch, branch_kraus, effective_readout,
                                     noise_kraus, rotation_unitary, two_rebit_states)
from bipartite_composition import interaction


def effect_coordinates(effect):
    return np.array([np.trace(effect @ p).real / 2 for p in (IDENTITY, PAULI_Z, PAULI_X)])


def noise_adjoint(effect):
    return sum(k.conj().T @ effect @ k for k in noise_kraus())


def protocol_effect(tree):
    """A leaf is its acceptance probability; a node is (angle, eta, child0, child1)."""
    if not isinstance(tree, tuple):
        if not 0 <= tree <= 1:
            raise ValueError("A leaf acceptance probability must lie in [0,1].")
        return float(tree) * IDENTITY
    angle, eta, *children = tree
    return sum(k.conj().T @ protocol_effect(children[s]) @ k
               for s in (0, 1) for k in branch_kraus(angle, s, eta))


def first_noise_preimage(tree):
    """Remove only the first noise channel; keep all later protocol branches."""
    if not isinstance(tree, tuple):
        return protocol_effect(tree)
    angle, eta, *children = tree
    rotation = rotation_unitary(angle)
    result = np.zeros((2, 2), dtype=complex)
    for s in (0, 1):
        sign = 2 * s - 1
        l = rotation @ np.diag(np.sqrt([(1 + sign * eta) / 2, (1 - sign * eta) / 2])) @ rotation.conj().T
        result += l @ protocol_effect(children[s]) @ l
    return result


def forward_acceptance(density, tree):
    if not isinstance(tree, tuple):
        return float(tree) * np.trace(density).real
    angle, eta, *children = tree
    return sum(forward_acceptance(apply_branch(density, angle, s, eta), children[s]) for s in (0, 1))


def random_protocol(rng, depth):
    if depth == 0:
        return float(rng.choice([0., 1., rng.random()]))
    return (float(rng.uniform(-math.pi, math.pi)), float(rng.uniform(0, 1)),
            random_protocol(rng, depth - 1), random_protocol(rng, depth - 1))


def one_read_implementation(e0, vector):
    """Construct an effect in the exact finite-protocol set using eta=1."""
    vector = np.asarray(vector, dtype=float)
    length = float(np.linalg.norm(vector))
    if not 0 <= e0 <= 1 or length > ALPHA * min(e0, 1 - e0) + 1e-14:
        raise ValueError("Effect lies outside the original local protocol set.")
    angle = math.atan2(vector[1], vector[0]) if length else 0.
    acceptance = (e0 - length / ALPHA, e0 + length / ALPHA)
    return sum(acceptance[s] * effective_readout(angle, s, 1) for s in (0, 1))


def compiled_copy_gate(angle=math.pi):
    """Ancilla is wire A; original input is B. No wire swap is assumed.

    Implements I_A tensor |0><0|_B + R_y(angle)_A tensor |1><1|_B
    using one YX gate and local rotations.
    """
    change_axis = np.kron(IDENTITY, rotation_unitary(-math.pi / 2))
    yz_part = change_axis @ interaction(-angle / 2, "YX") @ change_axis.conj().T
    return np.kron(rotation_unitary(angle / 2), IDENTITY) @ yz_part


def exact_copy_gate():
    p0, p1 = (IDENTITY + PAULI_Z) / 2, (IDENTITY - PAULI_Z) / 2
    return np.kron(IDENTITY, p0) + np.kron(-1j * PAULI_Y, p1)


def reusable_pointer_branch(density, outcome, contrast=1.):
    """Prepare A in |0>, couple A/B, read A, and discard A; keep input B."""
    p0 = (IDENTITY + PAULI_Z) / 2
    joint = np.kron(p0, density)
    gate = compiled_copy_gate()
    joint = gate @ joint @ gate.conj().T
    output = np.zeros((4, 4), dtype=complex)
    for k in branch_kraus(0, outcome, contrast):
        extended = np.kron(k, IDENTITY)
        output += extended @ joint @ extended.conj().T
    return np.einsum("abad->bd", output.reshape(2, 2, 2, 2))


def reusable_majority_probability(density, read_count, contrast=1.):
    majority_error(read_count, .1)
    result = 0.
    for record in product((0, 1), repeat=read_count):
        if sum(record) <= read_count // 2:
            continue
        state = density.copy()
        for outcome in record:
            state = reusable_pointer_branch(state, outcome, contrast)
        result += np.trace(state).real
    return result


def copy_isometry(pointer_count):
    if not isinstance(pointer_count, int) or pointer_count < 0:
        raise ValueError("Use a nonnegative integer number of pointers.")
    encoding = IDENTITY.copy()
    gate = compiled_copy_gate().reshape(2, 2, 2, 2)
    for previous_count in range(pointer_count):
        # Add |0> at the front, then act on the new first wire and the last
        # (input) wire. Existing pointer wires stay between them.
        expanded = np.kron(np.array([[1.], [0.]]), encoding)
        tensor = expanded.reshape(2, 2**previous_count, 2, 2)
        encoding = np.einsum("abcd,ckdj->akbj", gate, tensor).reshape(-1, 2)
    return encoding


def majority_error(pointer_count, error):
    if not isinstance(pointer_count, int) or pointer_count < 1 or pointer_count % 2 == 0:
        raise ValueError("Use a positive odd number of pointers.")
    if not 0 <= error <= 1:
        raise ValueError("An error probability must lie in [0,1].")
    return sum(math.comb(pointer_count, k) * error**k * (1 - error)**(pointer_count - k)
               for k in range((pointer_count + 1) // 2, pointer_count + 1))


def majority_effect_by_circuit(pointer_count, contrast=1.):
    # Small independent full-matrix calculation, not used for large certificates.
    majority_error(pointer_count, .1)  # Validate the count.
    encoding = copy_isometry(pointer_count)
    result = np.zeros((2, 2), dtype=complex)
    for outcomes in product((0, 1), repeat=pointer_count):
        if sum(outcomes) <= pointer_count // 2:
            continue
        effect = np.ones((1, 1), dtype=complex)
        for outcome in outcomes:
            effect = np.kron(effect, effective_readout(0, outcome, contrast))
        effect = np.kron(effect, IDENTITY)  # The original input is not read.
        result += encoding.conj().T @ effect @ encoding
    return result


def majority_certificate(pointer_count, contrast=Fraction(1)):
    eta = Fraction(contrast)
    if not 0 < eta <= 1:
        raise ValueError("Amplification requires a fixed contrast in (0,1].")
    alpha = 4 * sin_interval(I.rational(1, 4))
    error_lower = (1 - eta * Fraction(alpha.hi, SCALE)) / 2
    error_upper = (1 - eta * Fraction(alpha.lo, SCALE)) / 2
    # The binomial upper tail increases with the error probability. Evaluate
    # endpoint polynomials as exact Fractions, then round outward only once.
    lower = majority_error(pointer_count, error_lower)
    upper = majority_error(pointer_count, error_upper)
    enclosure = I(lower.numerator * SCALE // lower.denominator,
                  ceil_div(upper.numerator * SCALE, upper.denominator))
    return {"pointer_count": pointer_count, "contrast_exact": str(eta),
            "error_lower_exact": str(Fraction(enclosure.lo, SCALE)),
            "error_upper_exact": str(Fraction(enclosure.hi, SCALE)),
            "error_diagnostic": (enclosure.lo + enclosure.hi) / (2 * SCALE),
            "effective_visibility_diagnostic": 1 - (enclosure.lo + enclosure.hi) / SCALE,
            "strictly_positive_error": enclosure.lo > 0,
            "resource_interpretation": "m independent pointers, or m reset/read cycles of one reusable pointer",
            "scope": "Requires the declared real two-system gate and independent pure ancilla reset; the many-pointer version additionally assumes scalable real composition."}


def real_state_dimension(system_count):
    if not isinstance(system_count, int) or system_count < 1:
        raise ValueError("Use a positive integer number of real two-level systems.")
    size = 2**system_count
    homogeneous = size * (size + 1) // 2
    return {"systems": system_count, "matrix_size": size,
            "homogeneous_state_dimension": homogeneous,
            "normalized_state_dimension": homogeneous - 1,
            "local_product_statistic_dimension": 3**system_count}


class OperationalEffectClosureTests(unittest.TestCase):
    def test_state_dimension_counts_real_symmetric_matrices_not_spacetime(self):
        for n in (1, 2, 3):
            basis = []
            for labels in product(range(4), repeat=n):
                if labels.count(3) % 2:
                    continue
                operator = np.ones((1, 1), dtype=complex)
                for label in labels:
                    operator = np.kron(operator, (IDENTITY, PAULI_Z, PAULI_X, PAULI_Y)[label])
                np.testing.assert_array_equal(operator.imag, np.zeros_like(operator.real))
                basis.append(operator.reshape(-1))
            self.assertEqual(np.linalg.matrix_rank(basis), real_state_dimension(n)["homogeneous_state_dimension"])
        self.assertEqual([real_state_dimension(n)["homogeneous_state_dimension"] for n in (1, 2, 3)], [3, 10, 36])

    def test_first_noise_factorization_for_adaptive_protocols(self):
        rng = np.random.default_rng(48)
        for depth in range(1, 5):
            for _ in range(6):
                tree = random_protocol(rng, depth)
                effect, preimage = protocol_effect(tree), first_noise_preimage(tree)
                np.testing.assert_allclose(effect, noise_adjoint(preimage), atol=1e-15)
                self.assertGreaterEqual(np.linalg.eigvalsh(preimage).min(), -1e-15)
                self.assertLessEqual(np.linalg.eigvalsh(preimage).max(), 1 + 1e-15)
                e0, ez, ex = effect_coordinates(effect)
                self.assertLessEqual(math.hypot(ez, ex), ALPHA * min(e0, 1 - e0) + 1e-15)

    def test_protocol_effects_match_forward_branch_probabilities(self):
        tree = random_protocol(np.random.default_rng(148), 4)
        effect = protocol_effect(tree)
        for angle in (0., .41, 1.7, math.pi):
            density = (IDENTITY + math.cos(angle) * PAULI_Z + math.sin(angle) * PAULI_X) / 2
            self.assertAlmostEqual(forward_acceptance(density, tree), np.trace(density @ effect).real, places=14)
        complement = IDENTITY - effect
        e0, ez, ex = effect_coordinates(complement)
        self.assertLessEqual(math.hypot(ez, ex), ALPHA * min(e0, 1 - e0) + 1e-15)

    def test_every_effect_in_the_shrunken_set_has_a_one_read_construction(self):
        for e0, fraction, angle in product((0., .1, .5, .85, 1.), (0., .4, 1.), (0., .7, 2.1)):
            vector = ALPHA * min(e0, 1 - e0) * fraction * np.array([math.cos(angle), math.sin(angle)])
            actual = one_read_implementation(e0, vector)
            expected = e0 * IDENTITY + vector[0] * PAULI_Z + vector[1] * PAULI_X
            np.testing.assert_allclose(actual, expected, atol=3e-16)
        with self.assertRaises(ValueError):
            one_read_implementation(.5, [.5, 0])

    def test_repeated_readout_does_not_beat_single_input_discrimination_bound(self):
        plus, minus = (IDENTITY + PAULI_Z) / 2, (IDENTITY - PAULI_Z) / 2
        for adaptive in (False, True):
            distributions = []
            for density in (plus, minus):
                probabilities = []
                for record in product((0, 1), repeat=5):
                    state = density.copy()
                    for step, outcome in enumerate(record):
                        angle = .19 * sum(record[:step]) if adaptive else 0.
                        state = apply_branch(state, angle, outcome, .65 if adaptive else 1.)
                    probabilities.append(np.trace(state).real)
                distributions.append(np.array(probabilities))
            variation = np.abs(distributions[0] - distributions[1]).sum() / 2
            self.assertLessEqual(variation, ALPHA + 2e-15)
            if not adaptive:
                self.assertAlmostEqual(variation, ALPHA, places=14)

    def test_controlled_copy_is_compiled_from_the_previous_gate_and_local_rotations(self):
        p0, p1 = (IDENTITY + PAULI_Z) / 2, (IDENTITY - PAULI_Z) / 2
        for angle in (0., .21, 1.4, math.pi, -math.pi):
            gate = compiled_copy_gate(angle)
            expected = np.kron(IDENTITY, p0) + np.kron(rotation_unitary(angle), p1)
            np.testing.assert_allclose(gate, expected, atol=4e-16)
            np.testing.assert_allclose(compiled_copy_gate(-angle) @ gate, np.eye(4), atol=5e-16)
            np.testing.assert_array_equal(gate.imag, np.zeros((4, 4)))
        np.testing.assert_allclose(compiled_copy_gate(), exact_copy_gate(), atol=4e-16)

    def test_fanout_copies_only_the_selected_basis_information(self):
        for count in (1, 3, 5):
            encoding = copy_isometry(count)
            expected = np.zeros_like(encoding)
            expected[0, 0], expected[-1, 1] = 1, 1
            np.testing.assert_allclose(encoding, expected, atol=2e-15)
            np.testing.assert_allclose(encoding.conj().T @ encoding, IDENTITY, atol=3e-15)
        encoding = copy_isometry(1)
        input_plus_x = (IDENTITY + PAULI_X) / 2
        joint = (encoding @ input_plus_x @ encoding.conj().T).reshape(2, 2, 2, 2)
        pointer = np.einsum("abcb->ac", joint)
        np.testing.assert_allclose(pointer, IDENTITY / 2, atol=5e-16)
        self.assertGreater(np.linalg.norm(pointer - input_plus_x), .7)

    def test_reusable_pointer_instrument_preserves_the_read_axis_populations(self):
        for angle, eta, outcome in product((0., .7, math.pi / 2), (0., .15, 1.), (0, 1)):
            density = (IDENTITY + math.cos(angle) * PAULI_Z + math.sin(angle) * PAULI_X) / 2
            sign, beta = 2 * outcome - 1, ALPHA * eta
            expected = np.diag([(1 + sign * beta) * density[0, 0] / 2,
                                (1 - sign * beta) * density[1, 1] / 2])
            np.testing.assert_allclose(reusable_pointer_branch(density, outcome, eta), expected, atol=6e-16)
        density = (IDENTITY + PAULI_Z) / 2
        indirect = sum(reusable_pointer_branch(density, s, 1) for s in (0, 1))
        direct = sum(apply_branch(density, 0, s, 1) for s in (0, 1))
        self.assertAlmostEqual(np.trace(indirect @ PAULI_Z).real, 1.)
        self.assertAlmostEqual(np.trace(direct @ PAULI_Z).real, ALPHA)

    def test_one_resettable_pointer_matches_many_independent_pointers(self):
        for count, eta, angle in product((1, 3, 5), (.15, 1.), (0., .7, math.pi)):
            density = (IDENTITY + math.cos(angle) * PAULI_Z + math.sin(angle) * PAULI_X) / 2
            error = majority_error(count, (1 - ALPHA * eta) / 2)
            expected = (1 + (1 - 2 * error) * math.cos(angle)) / 2
            self.assertAlmostEqual(reusable_majority_probability(density, count, eta), expected, places=13)

    def test_full_circuit_effect_matches_the_binomial_prediction(self):
        for count, eta in product((1, 3, 5), (.15, .85, 1.)):
            error = majority_error(count, (1 - ALPHA * eta) / 2)
            expected = (IDENTITY + (1 - 2 * error) * PAULI_Z) / 2
            actual = majority_effect_by_circuit(count, eta)
            np.testing.assert_allclose(actual, expected, atol=3e-15)

    def test_rational_certificates_give_positive_finite_errors_and_strict_improvements(self):
        previous = Fraction(1)
        for count in (1, 3, 5, 7, 9):
            certificate = majority_certificate(count)
            lower, upper = Fraction(certificate["error_lower_exact"]), Fraction(certificate["error_upper_exact"])
            self.assertGreater(lower, 0)
            self.assertLess(upper, previous)
            self.assertLess(upper - lower, Fraction(1, 10**26))
            diagnostic = majority_error(count, (1 - ALPHA) / 2)
            self.assertAlmostEqual(float((lower + upper) / 2), diagnostic, places=16)
            previous = lower
        self.assertLess(Fraction(majority_certificate(9)["error_upper_exact"]), Fraction(1, 10**9))
        for count in (0, 2, -1):
            with self.assertRaises(ValueError):
                majority_certificate(count)

    def test_exponential_bound_and_uniform_ideal_effect_approximation(self):
        for count, eta in product((1, 3, 9, 31), (.1, .5, 1.)):
            beta = ALPHA * eta
            error = majority_error(count, (1 - beta) / 2)
            self.assertLessEqual(error, (1 - beta**2)**(count / 2) + 1e-16)
        error = majority_error(3, (1 - ALPHA) / 2)
        measured = majority_effect_by_circuit(3)
        for e0, length, angle in ((.5, .5, .7), (.2, .15, -1.2), (.8, .2, .3)):
            rotation = rotation_unitary(angle)
            noisy_effect = rotation @ measured @ rotation.conj().T
            actual = (e0 + length) * noisy_effect + (e0 - length) * (IDENTITY - noisy_effect)
            ideal = e0 * IDENTITY + length * (math.cos(angle) * PAULI_Z + math.sin(angle) * PAULI_X)
            self.assertLessEqual(np.linalg.norm(actual - ideal, ord=2), error + 2e-15)

    def test_every_strict_effect_is_exact_with_sufficient_finite_amplification(self):
        measured = majority_effect_by_circuit(5)
        visibility = 1 - 2 * majority_error(5, (1 - ALPHA) / 2)
        for e0, relative_length, angle in product((.1, .5, .8), (.9, .99999), (0., .6)):
            length = relative_length * min(e0, 1 - e0)
            acceptance_minus, acceptance_plus = e0 - length / visibility, e0 + length / visibility
            self.assertGreater(acceptance_minus, 0)
            self.assertLess(acceptance_plus, 1)
            rotation = rotation_unitary(angle)
            rotated = rotation @ measured @ rotation.conj().T
            actual = acceptance_plus * rotated + acceptance_minus * (IDENTITY - rotated)
            expected = e0 * IDENTITY + length * (math.cos(angle) * PAULI_Z + math.sin(angle) * PAULI_X)
            np.testing.assert_allclose(actual, expected, atol=3e-15)

    def test_every_finite_pointer_record_has_positive_probability(self):
        count, eta = 3, 1.
        encoding = copy_isometry(count)
        minimum = ((1 - ALPHA) / 2)**count
        for outcomes in product((0, 1), repeat=count):
            effect = np.ones((1, 1), dtype=complex)
            for outcome in outcomes:
                effect = np.kron(effect, effective_readout(0, outcome, eta))
            pulled = encoding.conj().T @ np.kron(effect, IDENTITY) @ encoding
            self.assertGreaterEqual(np.linalg.eigvalsh(pulled).min(), minimum * (1 - 1e-12))

    def test_sharp_local_effects_still_do_not_reveal_the_joint_yy_coordinate(self):
        difference = two_rebit_states()[0] - two_rebit_states()[1]
        for a, b in product((0., .37, 1.2, math.pi / 2), repeat=2):
            sharp_a = (IDENTITY + math.cos(a) * PAULI_Z + math.sin(a) * PAULI_X) / 2
            sharp_b = (IDENTITY + math.cos(b) * PAULI_Z + math.sin(b) * PAULI_X) / 2
            self.assertAlmostEqual(np.trace(difference @ np.kron(sharp_a, sharp_b)).real, 0., places=15)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(OperationalEffectClosureTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "finite_single_input_effect_set": "E=e0*I+ez*Z+ex*X; |e| <= alpha*min(e0,1-e0)",
        "local_protocol_optimal_antipodal_discrimination_success": (1 + ALPHA) / 2,
        "local_protocol_closure_contains_sharp_effects": False,
        "copy_gate_compiled_from_round_47_YX_and_local_rotations": True,
        "minimal_extra_capability": "One independently resettable pure ancilla, repeated uses of the round-47 two-system gate, and classical record processing",
        "reusable_pointer_selective_map": "R_s(rho)=sum_{b=0,1} [1+s*alpha*eta*(-1)^b]/2 * P_b*rho*P_b",
        "many_pointer_alternative_extra_assumptions": ["Arbitrarily many independently prepared pure real ancillas", "The same joint gate on each fresh ancilla and original input", "Consistent real tensor composition and independent local readout channels"],
        "finite_ancilla_effect_set": "{0,I} union {real E: 0<E<I}; conditional on gate and reset capabilities",
        "ancilla_effect_norm_closure": "All real 2x2 effects 0 <= E <= I; conditional on gate and reset capabilities",
        "finite_majority_protocol_is_exactly_sharp": False,
        "uniform_measurement_error_formula": "delta_m=sum_{k>m/2} binom(m,k)*epsilon^k*(1-epsilon)^(m-k), epsilon=(1-alpha*eta)/2",
        "exponential_error_bound": "delta_m <= (1-(alpha*eta)^2)^(m/2), fixed eta>0, odd m",
        "majority_certificates": [majority_certificate(m) for m in (1, 3, 5, 7, 9)],
        "real_state_dimensions": [real_state_dimension(n) for n in (1, 2, 3)],
        "local_tomography_follows_from_sharpening": False,
        "spacetime_dimension_or_string_theory_derived": False,
        "quantum_theory_derived_from_cognition": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("operational_effect_closure_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
