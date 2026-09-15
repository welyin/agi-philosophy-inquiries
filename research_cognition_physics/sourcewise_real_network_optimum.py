"""Round 67: optimize every Bob POVM for fixed sourcewise real references.

This is a restricted strategy optimum, not the general real-network bound.
An exact positive-semidefinite dual certificate removes the Bob optimizer.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bell_network_statistics import BELL_SIGNS, SETTINGS, bell_projectors, network_score
from bilocal_record_tomography import embed_operator
from independent_source_alignment import keep_systems, pair_reference
from role_symmetry_and_swap import pauli_word


def one_real_source():
    """Order: outer target, outer reference, middle target, middle reference."""
    return embed_operator(bell_projectors()[0], (0, 2), 4) @ embed_operator(pair_reference(), (1, 3), 4)


def outer_axes():
    return tuple(pauli_word(word) for word in ("XI", "YY", "ZI"))


def pulled_axes():
    return (pauli_word("XI"), -pauli_word("YY"), pauli_word("ZI"))


def conditional_middle(observable, sign, source=None):
    if source is None:
        source = one_real_source()
    return keep_systems(np.kron((np.eye(4)+sign*observable)/2, np.eye(4)) @ source, (2, 3), 4)


def target_only_bell_povm():
    # Bob wire order: B1 target, B1 reference, B2 target, B2 reference.
    return tuple(embed_operator(projector, (0, 2), 4) for projector in bell_projectors())


def score_operators():
    return tuple(pauli_word(word) for word in ("XIXI", "YYYY", "ZIZI"))


def dual_certificate():
    qs = score_operators()
    reference_parity = pauli_word("IYIY")
    dual = 2*np.eye(16)+reference_parity
    scores = tuple(sum(sign*q for sign, q in zip(signs, qs)) for signs in BELL_SIGNS)
    return dual, scores


def score_for_povm(effects):
    _, scores = dual_certificate()
    return 2*math.sqrt(2)/16*sum(np.trace(effect @ score).real for effect, score in zip(effects, scores))


def complete_real_records(effects=None):
    if effects is None:
        effects = target_only_bell_povm()
    axes = outer_axes()
    charlie = tuple((axes[i]+sign*axes[j])/math.sqrt(2) for i, j, sign in SETTINGS)
    left = {(x, a): conditional_middle(obs, a) for x, obs in enumerate(axes) for a in (-1, 1)}
    right = {(z, c): conditional_middle(obs, c) for z, obs in enumerate(charlie) for c in (-1, 1)}
    return {(x,z,a,b,c): float(np.trace(effects[b] @ np.kron(left[x,a], right[z,c])).real)
            for x,z,a,b,c in product(range(3),range(6),(-1,1),range(4),(-1,1))}


class SourcewiseRealNetworkOptimumTests(unittest.TestCase):
    def test_source_is_real_positive_normalized_with_maximally_mixed_middle(self):
        source = one_real_source()
        np.testing.assert_allclose(source.imag, 0, atol=0)
        self.assertGreater(np.linalg.eigvalsh(source).min(), -2e-16)
        self.assertAlmostEqual(np.trace(source).real, 1.)
        np.testing.assert_allclose(keep_systems(source, (2,3), 4), np.eye(4)/4, atol=1e-16)

    def test_outer_axes_are_real_symmetric_pairwise_anticommuting_reflections(self):
        axes = outer_axes()
        for axis in axes:
            np.testing.assert_array_equal(axis.imag, 0)
            np.testing.assert_array_equal(axis, axis.T)
            np.testing.assert_array_equal(axis @ axis, np.eye(4))
        for i,j,_ in SETTINGS:
            np.testing.assert_array_equal(axes[i]@axes[j]+axes[j]@axes[i], np.zeros((4,4)))

    def test_actual_source_contractions_give_the_pulled_back_effects(self):
        for outer, middle in zip(outer_axes(), pulled_axes()):
            for sign in (-1,1):
                np.testing.assert_allclose(conditional_middle(outer,sign), (np.eye(4)+sign*middle)/8, atol=6e-17)

    def test_score_operators_commute_and_their_product_tracks_reference_parity(self):
        qs = score_operators()
        for q in qs:
            np.testing.assert_array_equal(q@q, np.eye(16))
        for i,j,_ in SETTINGS:
            np.testing.assert_array_equal(qs[i]@qs[j], qs[j]@qs[i])
        np.testing.assert_array_equal(qs[0]@qs[1]@qs[2], -pauli_word("IYIY"))

    def test_dual_slacks_satisfy_exact_integer_polynomial_and_are_positive(self):
        dual, scores = dual_certificate()
        for score in scores:
            slack = dual-score
            np.testing.assert_array_equal(slack, slack.T)
            np.testing.assert_array_equal(slack@slack, 4*slack)
            self.assertGreater(np.linalg.eigvalsh(slack).min(), -2e-15)
        self.assertEqual(np.trace(dual), 32)

    def test_target_bell_measurement_saturates_every_complementary_slack(self):
        dual, scores = dual_certificate()
        effects = target_only_bell_povm()
        np.testing.assert_allclose(sum(effects), np.eye(16), atol=2e-16)
        for effect, score in zip(effects, scores):
            np.testing.assert_allclose((dual-score)@effect, 0, atol=4e-16)
        self.assertAlmostEqual(score_for_povm(effects), 4*math.sqrt(2), places=14)

    def test_full_record_calculation_agrees_with_the_linear_score_operators(self):
        records = complete_real_records()
        self.assertGreaterEqual(min(records.values()), 0.)
        self.assertAlmostEqual(network_score(records), score_for_povm(target_only_bell_povm()), places=13)
        for x,z in product(range(3),range(6)):
            self.assertAlmostEqual(sum(records[x,z,a,b,c] for a,b,c in product((-1,1),range(4),(-1,1))),1.)

    def test_random_general_povms_obey_dual_identity(self):
        rng = np.random.default_rng(67)
        dual, scores = dual_certificate()
        for _ in range(4):
            matrices = [rng.normal(size=(16,16)) for _ in range(4)]
            raw = [m@m.T for m in matrices]
            vals, vecs = np.linalg.eigh(sum(raw))
            inverse = (vecs/np.sqrt(vals))@vecs.T
            effects = [inverse@e@inverse for e in raw]
            deficit = sum(np.trace(e@(dual-m)).real for e,m in zip(effects,scores))
            self.assertGreaterEqual(deficit,0.)
            self.assertAlmostEqual(4*math.sqrt(2)-score_for_povm(effects),2*math.sqrt(2)/16*deficit,places=13)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SourcewiseRealNetworkOptimumTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round":67,"scope":"Fixed independent Bell-target times real-reference sources and fixed lifted outer measurements; every Bob POVM",
              "exact_optimum":"4*sqrt(2)","optimum_diagnostic":score_for_povm(target_only_bell_povm()),
              "dual_operator":"2I+Y_R1 Y_R2","exact_slack_identity":"(dual-M_b)^2=4*(dual-M_b)",
              "all_bob_outputs_retained":True,"general_real_network_bound":False,
              "quantum_theory_derived_from_cognition":False,
              "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("sourcewise_real_network_optimum_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
