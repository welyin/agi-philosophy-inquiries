"""Round 126: sharp one-copy universal exact joining bound for standard encodings.

Analytic proof: mixed orientation sectors require a forbidden linear conjugation
on at least one arbitrary pure input. Every success Kraus operator vanishes on
those sectors. Only 2 of 2**n equally weighted sectors can contribute.

An independent integer certificate covers two qubits with 16 product probes:
2048 linear constraints on 128 real Kraus entries have exact rank 126.
"""

import argparse
import json
import unittest
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np

from common_orientation_structure import encode_state, orientation
from complex_control_from_reference import real_lift
from encoded_composition_audit import (
    aligned_isometry, aligned_success_branch, independent_encoding, product_state,
)
from flagged_subject_joining import joining_kraus, random_state


def exact_success_bound(count):
    if not isinstance(count, int) or count < 1:
        raise ValueError("A positive integer subject count is required.")
    return Fraction(1, 2 ** (count - 1))


def qubit_probes():
    """Exactly dyadic projectors: |0>, |1>, |+X>, |+Y>."""
    return (
        np.diag([1., 0.]).astype(complex),
        np.diag([0., 1.]).astype(complex),
        np.array([[.5, .5], [.5, .5]], dtype=complex),
        np.array([[.5, -.5j], [.5j, .5]], dtype=complex),
    )


@lru_cache(maxsize=2)
def support_constraints(include_y=True):
    """Integer equations for (I-P_out) K P_in=0, row-major vectorization.

All operands and intermediate values are small dyadic numbers, exactly
representable in binary64. Multiplication by 16 is checked exactly before the
integer conversion. The subsequent rank certificate uses integer arithmetic.
"""
    probes = qubit_probes() if include_y else qubit_probes()[:3]
    rows = []
    for first, second in product(probes, repeat=2):
        pin = 4 * independent_encoding((first, second))
        pout = real_lift(np.kron(first, second))
        rows.append(np.kron(np.eye(8) - pout, pin.T))
    scaled = 16 * np.vstack(rows)
    if not np.array_equal(scaled, np.rint(scaled)):
        raise ArithmeticError("Noninteger coefficient in the exact certificate.")
    return scaled.astype(np.int64)


def rank_mod_1009(matrix):
    """1009 is prime. A nonzero minor here is also nonzero over the rationals."""
    modulus = 1009
    work = np.array(matrix, dtype=np.int64) % modulus
    rank = 0
    for col in range(work.shape[1]):
        candidates = np.flatnonzero(work[rank:, col])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        work[[rank, pivot]] = work[[pivot, rank]]
        work[rank] = (work[rank] * pow(int(work[rank, col]), -1, modulus)) % modulus
        factors = work[rank + 1:, col].copy()
        work[rank + 1:] = (work[rank + 1:] - factors[:, None] * work[rank]) % modulus
        rank += 1
        if rank == work.shape[0]:
            break
    return rank


def integer_kernel_basis():
    base = np.kron(np.array([[1, 0, 0, -1], [0, 1, 1, 0]], dtype=np.int64),
                   np.eye(4, dtype=np.int64))
    j = orientation(4).astype(np.int64)
    return base, j @ base


@lru_cache(maxsize=1)
def rank_certificate():
    equations = support_constraints()
    basis = integer_kernel_basis()
    if any(np.any(equations @ k.ravel()) for k in basis):
        raise ArithmeticError("The proposed kernel vectors do not solve the equations.")
    # Orthogonal nonzero integer vectors prove two independent kernel directions.
    gram = np.array([[int(np.sum(a * b)) for b in basis] for a in basis])
    if not np.array_equal(gram, 16 * np.eye(2, dtype=np.int64)):
        raise ArithmeticError("Independent kernel certificate failed.")
    rank = rank_mod_1009(equations)
    if rank != 126:
        raise ArithmeticError("The finite-field lower rank certificate failed.")
    return {
        "product_preparations": 16,
        "integer_equations": int(equations.shape[0]),
        "real_unknown_kraus_entries": int(equations.shape[1]),
        "prime_modulus": 1009,
        "rank_mod_prime": rank,
        "exact_rational_rank": 126,
        "exact_kernel_dimension": 2,
        "kernel_basis": "C, J_out C; C=sqrt(2) K_plus is integer",
        "kernel_gram_exact": gram.tolist(),
        "probability_bound_exact": "1/2",
        "rank_certificate_is_not_a_floating_svd": True,
    }


def sector_vector(signs):
    return product_state([np.array([1., 1j * sign]) / np.sqrt(2) for sign in signs]).ravel()


def sector_expansion(states):
    count = len(states)
    dimension = 2 ** count * int(np.prod([len(rho) for rho in states]))
    result = np.zeros((dimension, dimension), dtype=complex)
    for signs in product((1, -1), repeat=count):
        vector = sector_vector(signs)
        logical = product_state([rho if sign == 1 else rho.conj()
                                 for rho, sign in zip(states, signs)])
        result += np.kron(np.outer(vector, vector.conj()), logical) / 2 ** count
    return result


class UniversalJoiningBoundTests(unittest.TestCase):
    def test_complex_sector_expansion_is_only_algebra_and_matches_real_input(self):
        rng = np.random.default_rng(126)
        for dimensions in ((2, 2), (2, 3), (2, 2, 2)):
            states = [random_state(rng, d) for d in dimensions]
            np.testing.assert_allclose(sector_expansion(states), independent_encoding(states), atol=2e-16)

    def test_exact_integer_rank_certificate_exhausts_every_qubit_success_kraus(self):
        certificate = rank_certificate()
        self.assertEqual(certificate["exact_kernel_dimension"], 2)
        self.assertEqual(certificate["exact_rational_rank"], 126)

    def test_support_constraints_match_direct_matrix_action(self):
        rng = np.random.default_rng(226)
        k = rng.integers(-3, 4, size=(8, 16))
        expected = []
        for first, second in product(qubit_probes(), repeat=2):
            pin = 4 * independent_encoding((first, second))
            pout = real_lift(np.kron(first, second))
            expected.extend((16 * (np.eye(8) - pout) @ k @ pin).ravel())
        np.testing.assert_array_equal(support_constraints() @ k.ravel(), expected)

    def test_kernel_annihilates_both_misaligned_complex_sectors(self):
        identity = np.eye(4)
        for k in integer_kernel_basis():
            for signs in ((1, -1), (-1, 1)):
                embedding = np.kron(sector_vector(signs)[:, None], identity)
                np.testing.assert_allclose(k @ embedding, 0., atol=0.)

    def test_multi_kraus_success_probability_and_trace_nonincrease(self):
        rng = np.random.default_rng(326)
        states = [random_state(rng, 2), random_state(rng, 2)]
        source = independent_encoding(states)
        base, phase = (k / np.sqrt(2) for k in integer_kernel_basis())
        ks = [.3 * base + .4 * phase, .2 * base - .1 * phase]
        weight = .3 ** 2 + .4 ** 2 + .2 ** 2 + .1 ** 2
        output = sum(k @ source @ k.T for k in ks)
        np.testing.assert_allclose(output, weight / 2 * encode_state(product_state(states)), atol=1e-16)
        self.assertLessEqual(np.linalg.eigvalsh(sum(k.T @ k for k in ks)).max(), 1.)
        self.assertAlmostEqual(np.trace(output), .15)
        excessive = 1.01 * base
        self.assertGreater(np.linalg.eigvalsh(excessive.T @ excessive).max(), 1.)

    def test_filter_attains_general_bound_including_unequal_dimensions(self):
        rng = np.random.default_rng(426)
        for dimensions in ((2,), (2, 2), (2, 3), (2, 2, 2), (2, 2, 2, 2)):
            states = [random_state(rng, d) for d in dimensions]
            actual = aligned_success_branch(states)
            bound = exact_success_bound(len(states))
            self.assertAlmostEqual(np.trace(actual), float(bound))
            np.testing.assert_allclose(actual, float(bound) * encode_state(product_state(states)), atol=2e-16)

    def test_aligned_sector_projector_has_only_two_reference_dimensions(self):
        for count in range(2, 6):
            v = aligned_isometry(count)
            projector = v @ v.T
            self.assertAlmostEqual(np.trace(projector), 2.)
            for signs in product((1, -1), repeat=count):
                vector = sector_vector(signs)
                expected = vector if len(set(signs)) == 1 else np.zeros_like(vector)
                np.testing.assert_allclose(projector @ vector, expected, atol=6e-16)

    def test_omitting_complex_probes_allows_deterministic_real_input_join(self):
        equations = support_constraints(False)
        ks = joining_kraus((2, 2))
        for k in ks.values():
            np.testing.assert_allclose(equations @ k.ravel(), 0., atol=2e-15)
        np.testing.assert_allclose(sum(k.T @ k for k in ks.values()), np.eye(16), atol=5e-16)
        self.assertLess(rank_mod_1009(equations), 126)
        for first, second in product(qubit_probes()[:3], repeat=2):
            source = independent_encoding((first, second))
            output = sum(k @ source @ k.T for k in ks.values())
            np.testing.assert_allclose(output, encode_state(np.kron(first, second)), atol=3e-16)

    def test_resource_count_is_independent_orientation_clusters_not_all_old_members(self):
        self.assertEqual(exact_success_bound(2), Fraction(1, 2))
        self.assertEqual(exact_success_bound(5), Fraction(1, 16))
        self.assertEqual(exact_success_bound(1000).denominator.bit_length(), 1000)
        for invalid in (0, -1, 1.5):
            with self.assertRaises(ValueError):
                exact_success_bound(invalid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(UniversalJoiningBoundTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 126,
        "assumptions": {
            "standard_encoding": "E(rho)=real_lift(rho)/2",
            "independent_inputs": True,
            "one_copy_per_unknown_state": True,
            "all_complex_density_inputs_in_each_dimension_at_least_two": True,
            "allowed_success_map": "Arbitrary real completely positive trace-nonincreasing map",
            "fixed_input_independent_auxiliaries_and_finite_adaptation_included": True,
            "success_target": "Exactly E(tensor_i rho_i), after allowed decoding and correction",
        },
        "sharp_probability_bound_exact": "2^(1-n)",
        "one_new_subject_into_one_existing_unknown_cluster_exact": "1/2",
        "deterministic_universal_exact_join_exists_under_these_assumptions": False,
        "proof": "Every successful Kraus operator kills all mixed orientation sectors",
        "finite_qubit_certificate": rank_certificate(),
        "rows": [{"independent_clusters": n, "success_upper_bound_exact": str(exact_success_bound(n))}
                 for n in (1, 2, 3, 4, 8, 16)],
        "arbitrary_protocol_retrying_the_same_input_is_covered": True,
        "unknown_input_correlated_resource_or_extra_copies_covered": False,
        "approximate_joining_optimum_solved": False,
        "arbitrary_alternative_encodings_covered": False,
        "gate_communication_or_time_lower_bound_claimed": False,
        "attainment_uses_ideal_real_control_and_sector_readout": True,
        "attainment_compiled_into_old_noisy_primitives": False,
        "quantum_theory_derived_from_cognition": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("universal_joining_bound_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
