"""Real joint-state tomography from correlated records of blocks of size <= 2.

Bilocal tomography is a known property of real quantum theory. Here it is
implemented with the project's noisy readout and explicit real pair gate.
Pair marginals alone are deliberately distinguished from block joint records.
"""

import argparse
import json
import math
import unittest
from itertools import combinations, product
from pathlib import Path

import numpy as np

from bipartite_composition import interaction
from quantum_interface_audit import ALPHA, IDENTITY, effective_readout
from role_symmetry_and_swap import pauli_word


def real_words(n):
    return tuple("".join(word) for word in product("IXYZ", repeat=n) if word.count("Y") % 2 == 0)


def dimensions(n):
    total = (4**n + 2**n) // 2
    return {"systems": n, "real_homogeneous": total, "real_normalized": total - 1,
            "single_block_homogeneous": 3**n, "missing_from_single_blocks": total - 3**n,
            "all_pair_marginals_homogeneous": 1 + 2*n + 5*math.comb(n, 2)}


def blocks_for_word(word):
    positions = [i for i, letter in enumerate(word) if letter == "Y"]
    if len(positions) % 2:
        raise ValueError("A real symmetric Pauli word must have an even number of Y factors.")
    blocks = [(tuple(positions[j:j+2]), "YY") for j in range(0, len(positions), 2)]
    blocks.extend(((i,), letter) for i, letter in enumerate(word) if letter in "XZ")
    return tuple(blocks)


def embed_operator(operator, sites, n):
    sites = tuple(sites)
    if len(set(sites)) != len(sites) or any(i < 0 or i >= n for i in sites):
        raise ValueError("Sites must be distinct and in range.")
    if operator.shape != (2**len(sites), 2**len(sites)):
        raise ValueError("Operator dimension does not match its sites.")
    result = np.zeros((2**n, 2**n), dtype=complex)
    rest = tuple(i for i in range(n) if i not in sites)
    for row, column in product(range(2**n), repeat=2):
        rb = tuple((row >> (n - 1 - i)) & 1 for i in range(n))
        cb = tuple((column >> (n - 1 - i)) & 1 for i in range(n))
        if all(rb[i] == cb[i] for i in rest):
            r = sum(rb[i] << (len(sites) - 1 - j) for j, i in enumerate(sites))
            c = sum(cb[i] << (len(sites) - 1 - j) for j, i in enumerate(sites))
            result[row, column] = operator[r, c]
    return result


def block_effect(label, sign, contrast=1.):
    if sign not in (-1, 1):
        raise ValueError("Use a sign of -1 or +1.")
    outcome = (sign + 1) // 2
    if label == "YY":
        gate = interaction(math.pi / 2, "YX")
        return gate.conj().T @ np.kron(IDENTITY, effective_readout(0, outcome, contrast)) @ gate
    if label not in ("X", "Z"):
        raise ValueError("Unknown block observable.")
    return effective_readout(math.pi / 2 if label == "X" else 0, outcome, contrast)


def record_effect(word, signs, contrast=1.):
    blocks = blocks_for_word(word)
    if len(signs) != len(blocks):
        raise ValueError("One recorded sign is required for every nontrivial block.")
    effect = np.eye(2**len(word), dtype=complex)
    for (sites, label), sign in zip(blocks, signs):
        effect = effect @ embed_operator(block_effect(label, sign, contrast), sites, len(word))
    return effect


def record_distribution(density, word, contrast=1.):
    return {signs: float(np.trace(density @ record_effect(word, signs, contrast)).real)
            for signs in product((-1, 1), repeat=len(blocks_for_word(word)))}


def moment_from_records(distribution, word, contrast=1.):
    attenuation = (ALPHA * contrast)**len(blocks_for_word(word))
    if attenuation == 0:
        raise ValueError("Nonzero contrast is required for tomography.")
    return sum(math.prod(signs) * value for signs, value in distribution.items()) / attenuation


def reconstruct_from_moments(moments):
    n = len(next(iter(moments)))
    return sum(value * pauli_word(word) for word, value in moments.items()) / 2**n


def pair_marginal_words(n):
    return tuple(word for word in real_words(n) if sum(letter != "I" for letter in word) <= 2)


def partition_words(n, partition):
    sites = [i for block in partition for i in block]
    if sorted(sites) != list(range(n)):
        raise ValueError("A partition must cover every site exactly once.")
    return tuple(word for word in real_words(n)
                 if all(sum(word[i] == "Y" for i in block) % 2 == 0 for block in partition))


def split_record_counterexample(word="YYX"):
    operator = pauli_word(word)
    return (np.eye(len(operator)) + operator) / len(operator), (np.eye(len(operator)) - operator) / len(operator)


class BilocalRecordTomographyTests(unittest.TestCase):
    def test_exact_parity_counts_and_two_body_marginal_counts(self):
        for n in range(1, 8):
            self.assertEqual(len(real_words(n)), dimensions(n)["real_homogeneous"])
            self.assertEqual(len(pair_marginal_words(n)), dimensions(n)["all_pair_marginals_homogeneous"])
        self.assertEqual([dimensions(3)[key] for key in ("real_homogeneous", "single_block_homogeneous", "all_pair_marginals_homogeneous")], [36, 27, 22])

    def test_even_y_words_form_an_exact_real_orthogonal_matrix_basis(self):
        for n in range(1, 5):
            matrices = np.array([pauli_word(word) for word in real_words(n)])
            np.testing.assert_array_equal(matrices.imag, 0)
            np.testing.assert_array_equal(matrices, matrices.transpose(0, 2, 1))
            flat = matrices.reshape(len(matrices), -1)
            np.testing.assert_array_equal(flat.conj() @ flat.T, 2**n * np.eye(len(matrices)))

    def test_every_even_word_factors_into_disjoint_real_single_or_pair_blocks(self):
        for n in range(1, 7):
            for word in real_words(n):
                reconstructed = ["I"] * n
                for sites, label in blocks_for_word(word):
                    self.assertLessEqual(len(sites), 2)
                    for site, letter in zip(sites, label):
                        self.assertEqual(reconstructed[site], "I")
                        reconstructed[site] = letter
                self.assertEqual("".join(reconstructed), word)

    def test_pair_readout_is_an_existing_gate_followed_by_original_noisy_readout(self):
        for label, sign, eta in product(("X", "Z", "YY"), (-1, 1), (.2, .85, 1.)):
            expected = (np.eye(2**len(label)) + sign * ALPHA * eta * pauli_word(label)) / 2
            np.testing.assert_allclose(block_effect(label, sign, eta), expected, atol=3e-16)

    def test_nonadjacent_pair_embedding_preserves_site_order(self):
        np.testing.assert_array_equal(embed_operator(pauli_word("YY"), (0, 2), 3), pauli_word("YIY"))
        np.testing.assert_array_equal(embed_operator(pauli_word("XZ"), (2, 0), 3), pauli_word("ZIX"))

    def test_entire_block_povms_are_positive_normalized_and_have_expected_parity(self):
        for word in ("III", "XZY", "YYX", "YZY", "YYYY", "XZZX"):
            if word.count("Y") % 2:
                with self.assertRaises(ValueError):
                    blocks_for_word(word)
                continue
            effects = [(signs, record_effect(word, signs, .85))
                       for signs in product((-1, 1), repeat=len(blocks_for_word(word)))]
            for _, effect in effects:
                self.assertGreaterEqual(np.linalg.eigvalsh(effect).min(), -5e-16)
            np.testing.assert_allclose(sum(effect for _, effect in effects), np.eye(2**len(word)), atol=1e-15)
            parity = sum(math.prod(signs) * effect for signs, effect in effects)
            np.testing.assert_allclose(parity, (ALPHA * .85)**len(blocks_for_word(word)) * pauli_word(word), atol=7e-16)

    def test_all_three_system_moments_recover_a_density_matrix_from_actual_joint_records(self):
        rng = np.random.default_rng(51)
        matrix = rng.normal(size=(8, 8))
        density = matrix @ matrix.T
        density /= np.trace(density)
        moments = {word: moment_from_records(record_distribution(density, word, .73), word, .73)
                   for word in real_words(3)}
        np.testing.assert_allclose(reconstruct_from_moments(moments), density, atol=5e-16)

    def test_pair_marginals_and_all_single_block_joint_records_both_miss_yyx(self):
        plus, minus = split_record_counterexample()
        for density in (plus, minus):
            self.assertGreaterEqual(np.linalg.eigvalsh(density).min(), 0)
            self.assertAlmostEqual(np.trace(density).real, 1)
        observed = set(pair_marginal_words(3)) | set(partition_words(3, ((0,), (1,), (2,))))
        for word in observed:
            self.assertEqual(np.trace((plus - minus) @ pauli_word(word)).real, 0)
        self.assertEqual(len(observed), 30)

    def test_joint_block_records_distinguish_when_each_block_marginal_is_identical(self):
        plus, minus = split_record_counterexample()
        distributions = [record_distribution(density, "YYX") for density in (plus, minus)]
        for sign in (-1, 1):
            for column in (0, 1):
                for distribution in distributions:
                    self.assertAlmostEqual(sum(p for signs, p in distribution.items() if signs[column] == sign), .5)
        distance = sum(abs(distributions[0][s] - distributions[1][s]) for s in distributions[0]) / 2
        self.assertAlmostEqual(distance, ALPHA**2, places=14)
        for sign, distribution in zip((1, -1), distributions):
            for (a, b), probability in distribution.items():
                self.assertAlmostEqual(probability, (1 + sign * a * b * ALPHA**2) / 4, places=14)

    def test_a_single_fixed_partition_misses_other_pairings(self):
        partitions = (((0, 1), (2,)), ((0, 2), (1,)), ((1, 2), (0,)))
        spaces = [set(partition_words(3, partition)) for partition in partitions]
        self.assertEqual([len(space) for space in spaces], [30, 30, 30])
        self.assertEqual(len(set.union(*spaces)), 36)
        self.assertNotIn("YXY", spaces[0])

    def test_four_system_pairing_redundancy_does_not_overcount_yyyy(self):
        pairings = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
        self.assertTrue(all("YYYY" in partition_words(4, pairing) for pairing in pairings))
        self.assertEqual(dimensions(4)["real_homogeneous"], 81 + 54 + 1)
        self.assertEqual(len(set.union(*(set(partition_words(4, pairing)) for pairing in pairings))), 136)

    def test_zero_contrast_cannot_be_inverted_and_identity_needs_no_measurement(self):
        with self.assertRaises(ValueError):
            moment_from_records({(-1,): .5, (1,): .5}, "YY", 0)
        self.assertEqual(moment_from_records({(): 1.}, "III", 0), 1.)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(BilocalRecordTomographyTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"scope": "Usual real tensor composition with the declared pair gates and original nonzero noisy readout; not derived from cognition",
              "dimension_table": [dimensions(n) for n in range(1, 9)],
              "three_system_hidden_words": [word for word in real_words(3) if "Y" in word],
              "three_system_pair_marginal_dimension": 22,
              "three_system_single_records_plus_pair_marginals_dimension": 30,
              "single_fixed_pair_partition_dimension": 30,
              "three_pair_partitions_union_dimension": 36,
              "finite_setting_bilocal_tomography": True,
              "pair_marginals_are_sufficient": False,
              "state_counterexample": "rho_+/- = (I_8 +/- YYX)/8",
              "counterexample_block_joint_record_TV": ALPHA**2,
              "word_record_mean": "(alpha*eta)^(number of X/Z factors + number of Y factors/2) * Tr(rho P)",
              "exact_probabilities_are_not_finite_sample_estimates": True,
              "known_theorem_source": "https://arxiv.org/pdf/1005.4870, Sections 1, 4.4",
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("bilocal_record_tomography_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
