"""Round 164: native implementation and calibration budget for softer recovery.

The worst gate count remains 413 after the original join. Two reports use
two additional local rotations, so the expected count does increase. Angle
tolerance below concerns only the two reflection-sandwich Ry gates per run.
"""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import trace_distance
from certified_intervals import Interval as I, SCALE
from compiled_joining_frontier import compiled_sharing
from compiled_subject_joining import circuit_matrix, controlled_z, intake_circuit, ry
from compiled_three_subject_joining import inverse, gate_counts
from external_correlation_recovery_bound import channel
from global_history_recovery import LOWER as OLD_LOWER
from independent_source_alignment import keep_systems
from joint_relation_records import interval_fields
from minimal_coherent_joining_history import coherent_decoder_circuit
from noisy_classical_history import physical_read_branch, majority_error
from noisy_coherent_history_recovery import confusion_matrix, decode_label
from soft_history_decoder import PARAMETER, rational_circle, decoder, recovery_kraus, BIT_ERROR
from soft_decoder_certificate import UPPER
from walsh_coherent_history_compiler import walsh_physical_isometry, walsh_environment_circuit


MODIFIED = (0, 2, 3, 5)


def soft_decoder_circuit(report, parameter=PARAMETER, angle_errors=(0., 0.)):
    rational_circle(parameter)
    if type(report) is not int or not 0 <= report < 6:
        raise ValueError('Use a report in 0..5.')
    if len(angle_errors) != 2:
        raise ValueError('Specify the two reflection-sandwich rotation errors.')
    theta = 2*np.arctan(float(parameter))
    gates = controlled_z(0, 1, 2) if report < 3 else ()
    if report in MODIFIED:
        target = 2+report//3
        phi = theta if report % 3 == 0 else np.pi/2-theta
        before, after = -phi+angle_errors[0], phi+angle_errors[1]
        if before != 0:
            gates += (ry(target, before),)
        gates += controlled_z(0, target, 1)
        if after != 0:
            gates += (ry(target, after),)
    return gates+inverse(intake_circuit(2))


def actual_recovery(source, external_dimension=1, apply_decoder=True):
    if (type(external_dimension) is not int or external_dimension < 1
            or external_dimension & (external_dimension-1)):
        raise ValueError('Use a positive power-of-two external dimension.')
    if source.shape != (16*external_dimension, 16*external_dimension):
        raise ValueError('Wrong source dimension.')
    extra = external_dimension.bit_length()-1
    v = np.kron(walsh_physical_isometry(), np.eye(external_dimension))
    whole = v @ source @ v.T
    result = np.zeros_like(source)
    for bits in product((0, 1), repeat=3):
        branch = whole
        for site, bit in zip((4, 5, 6), bits):
            branch = physical_read_branch(branch, site, bit, .5)
        active = keep_systems(branch, (0, 1, 2, 3)+tuple(range(7, 7+extra)), 7+extra)
        if apply_decoder:
            label = decode_label(4*bits[0]+2*bits[1]+bits[2])
            repair = np.kron(circuit_matrix(soft_decoder_circuit(label), 4), np.eye(external_dimension))
            active = repair @ active @ repair.T
        result += active
    return result


def report_probabilities(error):
    matrix = confusion_matrix(error)
    return tuple(sum(matrix[r][g] for r in range(6))/6 for g in range(6))


def angle_certificate(epsilon=F(1, 1000)):
    if type(epsilon) not in (int, F) or epsilon < 0:
        raise ValueError('Use a nonnegative rational per-sandwich-Ry angle tolerance.')
    probabilities = report_probabilities(majority_error(1))
    affected = sum(probabilities[g] for g in MODIFIED)
    bound = I.exact(UPPER)+epsilon*affected
    return {'per_reflection_sandwich_Ry_tolerance_radians_exact': str(F(epsilon)),
            'affected_report_probability': interval_fields(affected),
            'external_error_upper_with_angle_tolerance': interval_fields(bound),
            'remaining_improvement_over_old_legal_lower': interval_fields(OLD_LOWER-bound),
            'still_strictly_beats_ideal_old_worst_error': F(bound.hi, SCALE) < OLD_LOWER,
            'covers_arbitrary_errors_of_all_native_gates': False,
            'other_gates_readout_contrast_and_coherent_storage_remain_ideal': True}


def resource_certificate():
    new = [gate_counts(soft_decoder_circuit(g)) for g in range(6)]
    old = [gate_counts(coherent_decoder_circuit(g+1)) for g in range(6)]
    environment = gate_counts(walsh_environment_circuit())
    probabilities = report_probabilities(majority_error(1))
    overhead = environment['all']+12
    expected_new = overhead+sum(p*c['all'] for p, c in zip(probabilities, new))
    expected_old = overhead+sum(p*c['all'] for p, c in zip(probabilities, old))
    extra = sum(p*(n['all']-o['all']) for p, n, o in zip(probabilities, new, old))
    return {'round': 164, 'decoder_gates_by_report': new, 'old_decoder_gates_by_report': old,
            'readout_count': 3, 'fresh_pointer_initializations': 3,
            'coherent_history_rebits': 1, 'raw_record_bits': 3, 'decoded_message_bits': 3,
            'new_environment_gates': environment['all'], 'new_yx_worst': environment['yx']+3+max(c['yx'] for c in new),
            'new_ry_worst': environment['ry']+9+max(c['ry'] for c in new),
            'new_native_gates_worst': overhead+max(c['all'] for c in new),
            'total_gates_including_initial_join_worst': overhead+max(c['all'] for c in new)+166,
            'expected_new_native_gates': interval_fields(expected_new),
            'expected_old_native_gates': interval_fields(expected_old),
            'expected_additional_native_gates': interval_fields(extra),
            'angle_tolerance_certificate': angle_certificate(),
            'old_common_interface_preserved_until_recovery': True,
            'all_eight_raw_outcomes_accepted_and_retained': True,
            'source_purification_and_noise_baths_closed_but_retained': True,
            'no_extra_pure_helpers_for_the_modified_decoders': True,
            'source_preparation_transport_calibration_and_classical_compute_costs_included': False,
            'new_global_error_bound_exact': str(UPPER)}


class CompiledSoftHistoryRecoveryTests(unittest.TestCase):
    def test_all_six_native_decoders_equal_the_exact_rational_maps(self):
        for parameter in (F(0), PARAMETER):
            for g in range(6):
                np.testing.assert_allclose(circuit_matrix(soft_decoder_circuit(g, parameter), 4), decoder(g, parameter), atol=4e-15)

    def test_same_worst_gate_count_and_honest_expected_increase(self):
        report = resource_certificate()
        self.assertEqual([r['all'] for r in report['decoder_gates_by_report']], [33, 21, 33, 23, 11, 23])
        self.assertEqual([r['all'] for r in report['old_decoder_gates_by_report']], [31, 21, 33, 21, 11, 23])
        self.assertEqual((report['new_yx_worst'], report['new_ry_worst'], report['new_native_gates_worst']), (293, 120, 413))
        self.assertEqual(report['total_gates_including_initial_join_worst'], 579)
        self.assertGreater(report['expected_additional_native_gates']['diagnostic'], .7)
        for e in (F(0), F(1, 4), F(1, 2)):
            p = report_probabilities(e)
            self.assertEqual(2*(p[0]+p[3]), (2+e-e*e)/3)
            self.assertEqual(sum(p[g] for g in MODIFIED), (4+e)/6)

    def test_actual_three_pointer_tree_matches_full_external_channel(self):
        rng = np.random.default_rng(164)
        raw = rng.normal(size=(32, 3))+1j*rng.normal(size=(32, 3))
        source = raw @ raw.conj().T; source /= np.trace(source)
        expected = channel(tuple(np.kron(k, np.eye(2)) for k in recovery_kraus()), source)
        np.testing.assert_allclose(actual_recovery(source, 2), expected, atol=4e-15)

    def test_nonselective_reads_preserve_the_old_common_interface(self):
        rng = np.random.default_rng(264)
        raw = rng.normal(size=(16, 4))+1j*rng.normal(size=(16, 4))
        source = raw @ raw.conj().T; source /= np.trace(source)
        stored = actual_recovery(source, apply_decoder=False)
        actual = keep_systems(stored, (1, 2, 3), 4)
        expected = channel(compiled_sharing(F(1, 2))[1], source)
        np.testing.assert_allclose(actual, expected, atol=4e-15)

    def test_two_rotation_angle_errors_obey_the_operator_bound(self):
        eps = .001
        for g in range(6):
            ideal = decoder(g)
            actual = circuit_matrix(soft_decoder_circuit(g, angle_errors=(eps, -eps)), 4)
            self.assertLessEqual(np.linalg.norm(actual-ideal, 2), eps+5e-15)

    def test_angle_error_channel_bound_includes_external_correlations(self):
        from minimal_coherent_joining_history import conditional_unitary
        rng = np.random.default_rng(364)
        vector = rng.normal(size=32)+1j*rng.normal(size=32); vector /= np.linalg.norm(vector)
        source = np.outer(vector, vector.conj()); matrix = confusion_matrix(BIT_ERROR); eps = .001
        perturbed = []
        for r in range(6):
            for g in range(6):
                repair = circuit_matrix(soft_decoder_circuit(g, angle_errors=((-1)**g*eps, (-1)**(g//2)*eps)), 4)
                perturbed.append(np.kron(np.sqrt(matrix[r][g]/6)*repair @ conditional_unitary(r+1), np.eye(2)))
        ideal = tuple(np.kron(k, np.eye(2)) for k in recovery_kraus())
        distance = trace_distance(channel(perturbed, source), channel(ideal, source))
        self.assertLessEqual(distance, eps*sum(report_probabilities(BIT_ERROR)[g] for g in MODIFIED)+2e-14)
        result = angle_certificate()
        self.assertTrue(result['still_strictly_beats_ideal_old_worst_error'])
        self.assertFalse(result['covers_arbitrary_errors_of_all_native_gates'])

    def test_resource_and_input_guards_do_not_hide_postselection_or_helpers(self):
        report = resource_certificate()
        self.assertTrue(report['all_eight_raw_outcomes_accepted_and_retained'])
        self.assertTrue(report['no_extra_pure_helpers_for_the_modified_decoders'])
        for value in (-F(1, 10), .001):
            with self.assertRaises(ValueError): angle_certificate(value)
        with self.assertRaises(ValueError): actual_recovery(np.eye(16)/16, 3)
        with self.assertRaises(ValueError): soft_decoder_circuit(True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CompiledSoftHistoryRecoveryTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = resource_certificate()
    report['automated_checks'] = {'run': checks.testsRun, 'failures': len(checks.failures), 'errors': len(checks.errors)}
    if args.write_results:
        Path(__file__).with_name('compiled_soft_history_recovery_results.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
