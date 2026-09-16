"""Round 135: sharp deterministic common-interface error for n independent qubits.

The exact bound optimizes all CPTP channels on the fixed standard encoded input.
Complex orientation projectors are used only to prove a real dual inequality.
The construction uses real Kraus maps; its general gate compilation is separate.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import unittest

import numpy as np

from approximate_joining_optimality import integer_dual, six_probes
from approximate_subject_joining import apply_kraus, logical_newcomer, trace_distance
from common_orientation_structure import encode_state, normalizer_flip
from complex_control_from_reference import real_lift
from encoded_composition_audit import independent_encoding, product_state
from flagged_subject_joining import joining_kraus, orientation_flags, random_state
from quantum_interface_audit import PAULI_X, PAULI_Y, PAULI_Z


def validate_count(count):
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise ValueError("A positive integer subject count is required.")
    return count


def optimal_support(count):
    n = validate_count(count)
    return sum((F(math.comb(n, k), 2**n)*F(2, 3)**min(k, n-k)
                for k in range(n+1)), F(0))


def optimal_error(count):
    return 1-optimal_support(count)


def sequential_error(count):
    return 1-F(5, 6)**(validate_count(count)-1)


def orientation_plan(flag):
    if not flag or flag[0] != 0 or any(bit not in (0, 1) for bit in flag):
        raise ValueError("Use canonical relative flags with first bit zero.")
    flip = int(sum(flag) > len(flag)/2)
    return flip, tuple(i for i, bit in enumerate(flag) if bit != flip)


def repaired_state(state):
    """Conditional repair of a conjugated qubit, not the averaged two-party map."""
    return (state + np.trace(state)*np.eye(2))/3


@lru_cache(maxsize=4)
def multisubject_kraus(count):
    n = validate_count(count)
    base = joining_kraus((2,)*n)
    paulis = (np.eye(2), PAULI_X.real, PAULI_Z.real)
    result = []
    for flag, intake in base.items():
        flip, wrong = orientation_plan(flag)
        gauge = normalizer_flip(2**n) if flip else np.eye(2**(n+1))
        for outcomes in product(range(3), repeat=len(wrong)):
            operators = [np.eye(2) for _ in range(n)]
            for site, outcome in zip(wrong, outcomes):
                operators[site] = paulis[outcome]/math.sqrt(3)
            repair = np.kron(np.eye(2), product_state(operators).real)
            result.append(repair @ gauge @ intake)
    return tuple(result)


def predicted_logical(states):
    n = len(states)
    result = np.zeros((2**n, 2**n), dtype=complex)
    for flag in orientation_flags(n):
        _, wrong = orientation_plan(flag)
        result += product_state([repaired_state(s) if i in wrong else s for i, s in enumerate(states)])
    return result/2**(n-1)


def y_projector(bit):
    return (np.eye(2) + (-1)**bit*PAULI_Y)/2


def reorder_qubits(matrix, order):
    n = len(order)
    return matrix.reshape((2,)*(2*n)).transpose(tuple(order)+tuple(n+i for i in order)).reshape(matrix.shape)


def moment_matrices():
    """Six-state exact second moments, with input then output for each target."""
    identity = np.eye(4, dtype=np.int64)
    omega = np.array([1, 0, 0, 1], dtype=np.int64)
    swap = np.array([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,1]], dtype=np.int64)
    return identity + np.outer(omega, omega), identity + swap  # 6 A, 6 B


def local_spectral_certificate():
    a, b = moment_matrices()
    identity = np.eye(4, dtype=np.int64)
    # A=(I+|Omega><Omega|)/6 >= 0. B and both upper-bound slacks
    # obey symmetric projector polynomials, so all eigenvalues are 0 or 2.
    checks = (b, 3*identity-a, 2*identity-b)
    for matrix in checks:
        if not np.array_equal(matrix, matrix.T) or np.any(matrix@matrix-2*matrix):
            raise ArithmeticError("Exact local PSD certificate failed.")
    return {"matched_moment_norm_exact": "1/2", "mismatched_moment_norm_exact": "1/3",
            "integer_projector_polynomials_verified": True,
            "requires_float_eigenvalues_for_proof": False}


def score_matrix_from_moments(count):
    """Audit helper; general proof uses local bounds instead of growing this array."""
    n = validate_count(count)
    a, b = (matrix/6 for matrix in moment_matrices())
    size = 2**(3*n+1)
    grouped = np.zeros((size, size), dtype=complex)
    for bits in product((0, 1), repeat=n):
        reference = product_state([y_projector(bit).T for bit in bits])
        for out_bit in (0, 1):
            moments = product_state([a if bit == out_bit else b for bit in bits])
            grouped += np.kron(np.kron(reference, y_projector(out_bit)), moments)/2**n
    order = (tuple(range(n)) + tuple(n+1+2*i for i in range(n))
             + (n,) + tuple(n+2+2*i for i in range(n)))
    result = reorder_qubits(grouped, order)
    if np.max(abs(result.imag)) > 1e-14:
        raise ArithmeticError("The paired orientation score must be real.")
    return result.real


def row(count):
    n = validate_count(count)
    return {"subjects": n, "optimal_support_exact": str(optimal_support(n)),
            "optimal_joint_trace_error_exact": str(optimal_error(n)),
            "specified_sequential_error_exact": str(sequential_error(n)),
            "error_improvement_exact": str(sequential_error(n)-optimal_error(n)),
            "joint_error_is_a_local_error": False}


class MultisubjectJoiningOptimumTests(unittest.TestCase):
    def test_exact_local_certificate_and_actual_six_state_moments(self):
        self.assertTrue(local_spectral_certificate()["integer_projector_polynomials_verified"])
        a, b = moment_matrices()
        np.testing.assert_array_equal(sum(np.kron(s.T, s) for s in six_probes()), a)
        np.testing.assert_array_equal(sum(np.kron(s, s) for s in six_probes()), b)

    def test_complexification_is_only_an_exact_decomposition_of_real_inputs(self):
        rng = np.random.default_rng(135)
        states = [random_state(rng, 2) for _ in range(3)]
        decomposed = np.zeros((64, 64), dtype=complex)
        for bits in product((0, 1), repeat=3):
            refs = product_state([y_projector(bit) for bit in bits])
            targets = product_state([s.conj() if bit else s for s, bit in zip(states, bits)])
            decomposed += np.kron(refs, targets)/8
        np.testing.assert_allclose(decomposed, independent_encoding(states), atol=1e-16)

    def test_factorized_score_reproduces_previous_full_integer_certificate(self):
        actual = score_matrix_from_moments(2)
        np.testing.assert_allclose(actual, integer_dual()[1]/576, atol=2e-17)

    def test_actual_real_maps_are_trace_preserving_for_one_to_four_subjects(self):
        for n in range(1, 5):
            ks = multisubject_kraus(n)
            np.testing.assert_allclose(sum(k.T@k for k in ks), np.eye(4**n), atol=2e-15)
            self.assertTrue(all(np.isrealobj(k) for k in ks))
        self.assertEqual(len(multisubject_kraus(3)), 10)

    def test_unknown_mixed_input_output_formula(self):
        rng = np.random.default_rng(235)
        for n in (2, 3, 4):
            states = [random_state(rng, 2) for _ in range(n)]
            output = apply_kraus(multisubject_kraus(n), independent_encoding(states))
            np.testing.assert_allclose(output, encode_state(predicted_logical(states)), atol=5e-16)
            self.assertLessEqual(trace_distance(output, encode_state(product_state(states))), float(optimal_error(n))+1e-14)

    def test_all_216_three_party_probes_and_random_four_party_pure_states_attain_bound(self):
        for states in product(six_probes(), repeat=3):
            output = apply_kraus(multisubject_kraus(3), independent_encoding(states))
            target = product_state(states)
            self.assertAlmostEqual(trace_distance(output, encode_state(target)), .25)
            self.assertAlmostEqual(np.trace(output@real_lift(target)), .75)
        rng = np.random.default_rng(335)
        for _ in range(4):
            states = []
            for _ in range(4):
                vector = rng.normal(size=2)+1j*rng.normal(size=2)
                vector /= np.linalg.norm(vector)
                states.append(np.outer(vector, vector.conj()))
            output = apply_kraus(multisubject_kraus(4), independent_encoding(states))
            self.assertAlmostEqual(trace_distance(output, encode_state(product_state(states))), 3/8)

    def test_retained_kraus_environment_is_reversible_including_external_correlations(self):
        v = np.vstack(multisubject_kraus(3))
        np.testing.assert_allclose(v.T@v, np.eye(64), atol=2e-15)
        rng = np.random.default_rng(435)
        joint = rng.normal(size=(64, 3))
        joint /= np.linalg.norm(joint)
        np.testing.assert_allclose(v.T@(v@joint), joint, atol=4e-16)

    def test_exact_table_and_large_n_joint_precision_limit(self):
        self.assertEqual([optimal_error(n) for n in range(1, 7)],
                         [F(0), F(1,6), F(1,4), F(3,8), F(65,144), F(469,864)])
        self.assertEqual(sequential_error(3)-optimal_error(3), F(1,18))
        for n in range(1, 81):
            self.assertLessEqual(optimal_error(n), sequential_error(n))
            self.assertLessEqual(optimal_support(n), 2*F(5,6)**n)
        for invalid in (0, -1, 2.5, True):
            with self.assertRaises(ValueError):
                optimal_support(invalid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(MultisubjectJoiningOptimumTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 135,
        "model": "One independent standard real encoding of each arbitrary unknown qubit",
        "optimized_channels": "All deterministic real CPTP channels, also bounded for complex CPTP on these inputs",
        "optimal_support_formula": "2^(-n) sum_k binomial(n,k) (2/3)^min(k,n-k)",
        "worst_joint_error_formula": "1 - optimal_support(n)",
        "upper_bound_on_support": "2 (5/6)^n",
        "proof": "Exact local second moments and tensor-product operator bounds give a Choi dual",
        "local_certificate": local_spectral_certificate(),
        "old_host_exact_preservation_imposed_on_global_optimum": False,
        "all_environments_can_be_retained": True,
        "general_n_native_gate_compilation_completed": False,
        "rows": [row(n) for n in range(1, 13)],
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("multisubject_joining_optimum_results.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
