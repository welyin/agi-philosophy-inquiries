"""Round 63: delayed flag processing recovers Pauli tomography statistics.

There is no physical feedforward before the target measurements. Independent
reference sources prepare a resource for a separate, independent target state.
The target is an additional input, not a two-source Bell experiment.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator, real_words
from independent_source_alignment import conditional_reference, keep_systems, middle_branch
from quantum_interface_audit import ALPHA, two_rebit_states
from role_symmetry_and_swap import pauli_word
from shared_real_reference_tomography import assisted_records


def relabel_c_outcome(word, outcomes, flag):
    labels = [i for i, axis in enumerate(word) if axis != "I"]
    corrected = list(outcomes)
    if word[2] == "Y":
        corrected[labels.index(2)] *= flag
    return tuple(corrected)


def tagged_records(density, word, alignment_visibility=ALPHA):
    if len(word) != 3:
        raise ValueError("Use three target systems A,B,C.")
    return {(flag,)+outcomes: probability/2
            for flag in (-1, 1)
            for outcomes, probability in assisted_records(density, word, conditional_reference(flag, alignment_visibility)).items()}


def corrected_records(density, word, alignment_visibility=ALPHA):
    result = {}
    for (flag, *outcomes), probability in tagged_records(density, word, alignment_visibility).items():
        corrected = relabel_c_outcome(word, outcomes, flag)
        result[corrected] = result.get(corrected, 0.)+probability
    return result


def reconstruct_delayed(density, alignment_visibility=ALPHA):
    if not 0 < alignment_visibility <= 1:
        raise ValueError("Complete reconstruction needs positive alignment visibility.")
    output = np.zeros((8, 8), dtype=complex)
    for word in real_words(3):
        weight = sum(axis != "I" for axis in word)
        factor = alignment_visibility if word[2] == "Y" else 1.
        records = corrected_records(density, word, alignment_visibility)
        moment = sum(math.prod(outcomes)*p for outcomes, p in records.items())/(ALPHA**weight*factor)
        output += moment*pauli_word(word)/8
    return output


def outer_tagged_records(density, alignment_visibility=ALPHA):
    result = {}
    for flag in (-1, 1):
        reference = keep_systems(conditional_reference(flag, alignment_visibility), (0, 2), 3)
        for outcomes, p in assisted_records(density, "YY", reference).items():
            result[(flag,)+outcomes] = p/2
    return result


class DelayedReferenceRecordsTests(unittest.TestCase):
    def test_actual_middle_branches_give_the_declared_joint_record_probabilities(self):
        rho = (np.eye(8)+.3*pauli_word("YYX")+.2*pauli_word("XYY"))/8
        for word in ("YYX", "XYY", "YIY", "XZZ"):
            actual = tagged_records(rho, word)
            for flag in (-1, 1):
                reference = 2*middle_branch(flag)
                direct = assisted_records(rho, word, reference)
                for outcomes, p in direct.items():
                    self.assertAlmostEqual(actual[(flag,)+outcomes], p/2, places=14)

    def test_delayed_relabeling_matches_physical_reference_correction_for_full_pauli_records(self):
        rng = np.random.default_rng(63)
        matrix = rng.normal(size=(8, 8))
        rho = matrix @ matrix.T
        rho /= np.trace(rho)
        for word in real_words(3):
            delayed = corrected_records(rho, word)
            physical = assisted_records(rho, word, conditional_reference(1, ALPHA))
            for key in physical:
                self.assertAlmostEqual(delayed[key], physical[key], places=14)

    def test_all_thirty_six_coordinates_reconstruct_with_delayed_classical_flags(self):
        rng = np.random.default_rng(163)
        matrix = rng.normal(size=(8, 8))
        rho = matrix @ matrix.T
        rho /= np.trace(rho)
        for gamma in (.2, ALPHA):
            np.testing.assert_allclose(reconstruct_delayed(rho, gamma), rho, atol=1e-15)

    def test_outer_complete_records_follow_the_exact_triple_correlation_formula(self):
        for q, rho in zip((1, -1), two_rebit_states()):
            for (flag, s, t), p in outer_tagged_records(rho).items():
                self.assertAlmostEqual(p, (1+flag*s*t*ALPHA**3*q)/8, places=14)

    def test_discarding_the_middle_flag_erases_the_outer_record_distinction(self):
        tables = [outer_tagged_records(rho) for rho in two_rebit_states()]
        for s, t in product((-1, 1), repeat=2):
            for records in tables:
                self.assertAlmostEqual(sum(records[flag, s, t] for flag in (-1, 1)), .25, places=14)
        distance = sum(abs(tables[0][key]-tables[1][key]) for key in tables[0])/2
        self.assertAlmostEqual(distance, ALPHA**3, places=14)

    def test_middle_postselection_alone_has_probability_half(self):
        for rho in two_rebit_states():
            records = outer_tagged_records(rho)
            self.assertAlmostEqual(sum(p for (flag, _, _), p in records.items() if flag == 1), .5, places=14)

    def test_ignored_middle_channel_cannot_change_outer_marginals_in_either_scalar_field(self):
        rng = np.random.default_rng(263)
        for complex_model in (False, True):
            sources = []
            for _ in range(2):
                matrix = rng.normal(size=(4, 4))
                if complex_model:
                    matrix = matrix+1j*rng.normal(size=(4, 4))
                density = matrix @ matrix.conj().T
                sources.append(density/np.trace(density))
            initial = np.kron(*sources)
            matrix = rng.normal(size=(16, 4))
            if complex_model:
                matrix = matrix+1j*rng.normal(size=(16, 4))
            isometry, _ = np.linalg.qr(matrix)
            kraus = [isometry[j:j+4] for j in range(0, 16, 4)]
            np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), np.eye(4), atol=8e-16)
            output = np.zeros((16, 16), dtype=complex)
            for k in kraus:
                extended = embed_operator(k, (1, 2), 4)
                output += extended @ initial @ extended.conj().T
            expected = np.kron(keep_systems(sources[0], (0,), 2), keep_systems(sources[1], (1,), 2))
            np.testing.assert_allclose(keep_systems(output, (0, 3), 4), expected, atol=4e-16)

    def test_record_equivalence_does_not_assert_equal_uncorrected_physical_states(self):
        uncorrected = sum(conditional_reference(flag, ALPHA)/2 for flag in (-1, 1))
        corrected = conditional_reference(1, ALPHA)
        self.assertAlmostEqual(np.trace(uncorrected @ pauli_word("YIY")).real, 0)
        self.assertAlmostEqual(np.trace(corrected @ pauli_word("YIY")).real, ALPHA)

    def test_zero_alignment_cannot_supply_a_complete_calibration(self):
        with self.assertRaises(ValueError):
            reconstruct_delayed(np.eye(8)/8, 0.)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DelayedReferenceRecordsTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 63, "scope": "Pauli tomography of an independent real target using a two-source reference preparation stage",
              "physical_feedforward_before_target_measurement_required": False,
              "middle_flag_must_be_available_for_record_analysis": True,
              "target_source_counted_as_an_additional_input": True,
              "three_target_homogeneous_coordinates_recovered": 36,
              "outer_complete_record_tv": ALPHA**3, "outer_record_tv_if_flag_discarded": 0,
              "postselection_success_probability_if_used": 0.5,
              "general_measurement_programs_classically_correctable_claimed": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("delayed_reference_records_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
