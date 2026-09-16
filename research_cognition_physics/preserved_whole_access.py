"""Round 79: an unchanged whole with different local access and message rules.

The preparation always retains both targets and both shared reference systems.
Removing access never deletes a register or applies a state-erasure channel.
The exact message formula concerns the specified two local noisy YY readings.
"""

import argparse
import json
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from differentiation_channel_audit import distributed_whole
from independent_source_alignment import keep_systems
from quantum_interface_audit import ALPHA
from role_symmetry_and_swap import pauli_word


SIGNS = (-1, 1)


def message_channel(rows):
    channel = np.asarray(rows, dtype=float)
    if channel.ndim != 2 or channel.shape[0] != 2 or channel.shape[1] < 1:
        raise ValueError("Two input rows and at least one message symbol required.")
    if not np.all(np.isfinite(channel)) or np.any(channel < 0):
        raise ValueError("Message probabilities must be finite and nonnegative.")
    if not np.allclose(channel.sum(axis=1), 1, atol=1e-14, rtol=0):
        raise ValueError("Each input row must sum to one.")
    return channel


def joint_records(whole, reference_access=True, visibility=ALPHA):
    if not 0 <= visibility <= 1:
        raise ValueError("Visibility must lie in [0,1].")
    # Both effects act on the SAME four-register preparation.
    # XI accesses only the target; YY also accesses the local reference.
    observable = pauli_word("YY" if reference_access else "XI")
    effects = {s: (np.eye(4) + s * visibility * observable) / 2 for s in SIGNS}
    return {(a, c): float(np.trace(whole @ np.kron(effects[a], effects[c])).real)
            for a, c in product(SIGNS, repeat=2)}


def receiver_records(whole, channel, reference_access=True, visibility=ALPHA):
    channel = message_channel(channel)
    joint = joint_records(whole, reference_access, visibility)
    return {(m, c): sum(joint[a, c] * channel[i, m] for i, a in enumerate(SIGNS))
            for m, c in product(range(channel.shape[1]), SIGNS)}


def total_variation(first, second):
    return sum(abs(first[key] - second[key]) for key in first) / 2


def message_distinguishability(channel):
    channel = message_channel(channel)
    return np.abs(channel[1] - channel[0]).sum() / 2


def formula_record(q, channel, visibility=ALPHA):
    channel = message_channel(channel)
    return {(m, c): (channel[0, m] + channel[1, m]
                    + c * visibility**2 * q * (channel[1, m] - channel[0, m])) / 4
            for m, c in product(range(channel.shape[1]), SIGNS)}


class PreservedWholeAccessTests(unittest.TestCase):
    def test_whole_is_real_normalized_and_global_alternatives_remain_distinct(self):
        plus, minus = distributed_whole(1), distributed_whole(-1)
        for whole in (plus, minus):
            np.testing.assert_array_equal(whole.imag, np.zeros_like(whole.imag))
            self.assertAlmostEqual(np.trace(whole).real, 1)
            self.assertGreaterEqual(np.linalg.eigvalsh(whole).min(), -2e-16)
        self.assertAlmostEqual(np.abs(np.linalg.eigvalsh(plus - minus)).sum()/2, 1)
        np.testing.assert_array_equal(plus @ minus, np.zeros((16, 16)))

    def test_each_entire_local_lab_has_identical_marginal_even_with_its_reference(self):
        for q in (-1, 0, 1):
            whole = distributed_whole(q)
            for sites in ((0, 1), (2, 3)):
                np.testing.assert_array_equal(keep_systems(whole, sites, 4), np.eye(4)/4)

    def test_changing_access_does_not_change_or_discard_the_whole(self):
        for q in (-1, 1):
            whole = distributed_whole(q)
            before = whole.copy()
            whole.setflags(write=False)
            for access in (False, True, False):
                receiver_records(whole, np.eye(2), access)
                np.testing.assert_array_equal(whole, before)
            self.assertEqual(whole.shape, (16, 16))

    def test_target_only_real_product_effects_annihilate_the_global_difference(self):
        delta = distributed_whole(1) - distributed_whole(-1)
        # These span all real symmetric single-target effects.
        for left, right in product("IXZ", repeat=2):
            effect = np.kron(pauli_word(left + "I"), pauli_word(right + "I"))
            self.assertEqual(np.trace(delta @ effect), 0)
        # A joint real YY target effect retains the distinction in the whole.
        self.assertEqual(np.trace(delta @ pauli_word("YIYI")).real, 2)

    def test_full_matrix_records_match_the_message_formula_for_general_channels(self):
        rng = np.random.default_rng(79)
        for size in (1, 2, 3, 5):
            channel = rng.random((2, size))
            channel /= channel.sum(axis=1, keepdims=True)
            for q in (-1, -.4, 0, .7, 1):
                actual = receiver_records(distributed_whole(q), channel)
                expected = formula_record(q, channel)
                self.assertAlmostEqual(sum(actual.values()), 1)
                for key in actual:
                    self.assertAlmostEqual(actual[key], expected[key], places=15)

    def test_receiver_distinction_is_visibility_squared_times_message_tv(self):
        channels = ([[1], [1]], np.eye(2), [[.75, .25], [.25, .75]],
                    [[.75, 0, .25], [0, .75, .25]], [[.6, .3, .1], [.2, .2, .6]])
        for channel in channels:
            for visibility in (0, .3, ALPHA, 1):
                first = receiver_records(distributed_whole(1), channel, True, visibility)
                second = receiver_records(distributed_whole(-1), channel, True, visibility)
                self.assertAlmostEqual(total_variation(first, second),
                                       visibility**2 * message_distinguishability(channel), places=15)

    def test_exact_bayes_decision_and_all_binary_postprocessings(self):
        channel = [[.6, .3, .1], [.2, .2, .6]]
        first = receiver_records(distributed_whole(1), channel)
        second = receiver_records(distributed_whole(-1), channel)
        keys = tuple(first)
        minimum = min(sum(second[k] if decision else first[k]
                          for k, decision in zip(keys, decisions))/2
                      for decisions in product((False, True), repeat=len(keys)))
        self.assertAlmostEqual(minimum, (1-total_variation(first, second))/2)
        # Maximum-likelihood rule: sign(c*(W(m|+)-W(m|-))).
        rule_error = 0.
        for (m, c) in keys:
            evidence = c*(channel[1][m] - channel[0][m])
            rule_error += (second[m, c] if evidence > 0 else first[m, c] if evidence < 0
                           else (first[m, c]+second[m, c])/2)/2
        self.assertAlmostEqual(rule_error, minimum)

    def test_no_message_or_no_reference_access_gives_no_receiver_information(self):
        for channel, access in (([[1], [1]], True), (np.eye(2), False)):
            first = receiver_records(distributed_whole(1), channel, access)
            second = receiver_records(distributed_whole(-1), channel, access)
            self.assertEqual(first, second)
        for invalid in ([[1, 1], [0, 1]], [[-.1, 1.1], [0, 1]], [[1]], [[float("nan")], [1]]):
            with self.assertRaises(ValueError):
                message_channel(invalid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(PreservedWholeAccessTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 79,
        "user_clarification": "The whole remains; differentiation changes local access and cooperation.",
        "whole_preparation": "rho_q(targets) tensor tau(shared references), regrouped into two local labs",
        "differentiation_applies_a_state_erasure_channel": False,
        "references_remain_in_the_whole_when_access_is_denied": True,
        "global_state_pair_trace_distance_exact": "1",
        "each_local_lab_state_exact": "I_4/4, independent of q",
        "no_reference_access_even_with_real_local_protocols_and_messages_TV_exact": "0",
        "local_reference_access_without_a_message_receiver_TV_exact": "0",
        "fixed_protocol_with_message_receiver_TV_exact": "beta^2 * TV(W(.|+), W(.|-))",
        "one_noiseless_bit_original_readout_TV_diagnostic": ALPHA**2,
        "one_noiseless_bit_equal_prior_error_diagnostic": (1-ALPHA**2)/2,
        "one_noiseless_bit_internal_old_readouts": 2,
        "fixed_readout_protocol_is_globally_optimal_claimed": False,
        "unchanged_real_whole_is_a_counterexample_to_complex_necessity": True,
        "quantum_theory_derived_from_cognition": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}
    }
    if args.write_results:
        Path(__file__).with_name("preserved_whole_access_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
