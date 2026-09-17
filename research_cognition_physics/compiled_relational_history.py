"""Round 170: compile the coherent relation correction in 18 original gates.

Two native conjugation chains turn a controlled Ry into the required real
Pauli rotations. The borrowed wires are arbitrary; no new pure helper or
history read is introduced. Four variable pulse errors are audited separately.
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
from compiled_soft_history_recovery import report_probabilities
from compiled_subject_joining import (Gate, circuit_matrix, controlled_ry,
    intake_circuit, ry, yx)
from compiled_three_subject_joining import inverse, gate_counts
from external_correlation_recovery_bound import channel
from four_angle_history_decoder import decoder_circuit as previous_decoder_circuit
from four_angle_minimax_bound import LOWER as PREVIOUS_FAMILY_LOWER
from independent_source_alignment import keep_systems
from joint_relation_records import interval_fields
from minimal_coherent_joining_history import conditional_unitary
from noisy_classical_history import physical_read_branch, majority_error
from noisy_coherent_history_recovery import decode_label, confusion_matrix
from relational_history_decoder import PARAMETER, correction, decoder, recovery_kraus, BIT_ERROR
from relational_history_certificate import UPPER
from soft_history_decoder import rational_circle
from walsh_coherent_history_compiler import walsh_physical_isometry, walsh_environment_circuit


def correction_circuit(parameter=PARAMETER, pulse_errors=(0.,0.,0.,0.)):
    rational_circle(parameter)
    if len(pulse_errors)!=4 or not all(np.isfinite(x) for x in pulse_errors):
        raise ValueError('Use four finite errors for the variable YX/Ry pulses.')
    if parameter==0 and all(x==0 for x in pulse_errors): return ()
    angle=2*np.arctan(float(parameter))
    gchain=(yx(2,1,np.pi/2),yx(3,2,np.pi/2),ry(2,np.pi/2))
    hchain=(yx(3,1,np.pi/2),yx(2,3,np.pi/2))
    def controlled(angle, errors):
        gates=list(controlled_ry(0,1,angle))
        for index,error in zip((1,3),errors):
            gate=gates[index];gates[index]=Gate(gate.kind,gate.sites,gate.angle+error)
        return tuple(gates)
    return (inverse(gchain)+controlled(-2*angle,pulse_errors[:2])+gchain
            +inverse(hchain)+controlled(2*angle,pulse_errors[2:])+hchain)


def decoder_circuit(report, parameter=PARAMETER, pulse_errors=(0.,0.,0.,0.)):
    previous=previous_decoder_circuit(report)
    tail=inverse(intake_circuit(2))
    if previous[-len(tail):]!=tail:
        raise ValueError('The previous decoder no longer ends with the required intake inverse.')
    return previous[:-len(tail)]+correction_circuit(parameter,pulse_errors)+tail


def actual_recovery(source, external_dimension=1, apply_decoder=True):
    if (type(external_dimension) is not int or external_dimension<1
            or external_dimension & (external_dimension-1)):
        raise ValueError('Use a positive power-of-two external dimension.')
    if source.shape!=(16*external_dimension,16*external_dimension):
        raise ValueError('Wrong source dimension.')
    extra=external_dimension.bit_length()-1
    v=np.kron(walsh_physical_isometry(),np.eye(external_dimension));whole=v @ source @ v.T
    result=np.zeros_like(source)
    for bits in product((0,1),repeat=3):
        branch=whole
        for site,bit in zip((4,5,6),bits):
            branch=physical_read_branch(branch,site,bit,.5)
        active=keep_systems(branch,(0,1,2,3)+tuple(range(7,7+extra)),7+extra)
        if apply_decoder:
            label=decode_label(4*bits[0]+2*bits[1]+bits[2])
            repair=np.kron(circuit_matrix(decoder_circuit(label),4),np.eye(external_dimension))
            active=repair @ active @ repair.T
        result+=active
    return result


def pulse_certificate(tolerance=F(1,20000)):
    if type(tolerance) not in (int,F) or tolerance<0:
        raise ValueError('Use a nonnegative rational pulse-angle tolerance.')
    # Each native rotation has generator norm 1/2. Four independent errors
    # add at most 2*tolerance in operator norm, hence also in external distance.
    upper=I.exact(UPPER+2*tolerance)
    return {'per_variable_native_pulse_error_radians_exact':str(F(tolerance)),
            'variable_pulses_affected':4,'variable_yx_pulses':2,'variable_ry_pulses':2,
            'external_error_upper':interval_fields(upper),
            'remaining_improvement_over_ideal_four_angle_family':interval_fields(PREVIOUS_FAMILY_LOWER-upper),
            'still_beats_every_ideal_four_angle_scheme':F(upper.hi,SCALE)<PREVIOUS_FAMILY_LOWER,
            'all_other_gates_readout_contrast_and_coherent_storage_ideal':True,
            'all_native_gate_noise_covered':False}


def resource_certificate():
    conditional=[gate_counts(decoder_circuit(g)) for g in range(6)]
    previous=[gate_counts(previous_decoder_circuit(g)) for g in range(6)]
    common=gate_counts(correction_circuit());environment=gate_counts(walsh_environment_circuit())
    p=report_probabilities(majority_error(1));overhead=environment['all']+12
    return {'round':170,'correction_native_counts':common,
            'conditional_native_counts':conditional,
            'new_native_gates_worst':overhead+max(c['all'] for c in conditional),
            'new_yx_worst':environment['yx']+3+max(c['yx'] for c in conditional),
            'new_ry_worst':environment['ry']+9+max(c['ry'] for c in conditional),
            'total_including_first_join_worst':166+overhead+max(c['all'] for c in conditional),
            'expected_new_native_gates':interval_fields(overhead+sum(prob*c['all'] for prob,c in zip(p,conditional))),
            'expected_previous_native_gates':interval_fields(overhead+sum(prob*c['all'] for prob,c in zip(p,previous))),
            'additional_gates_on_every_report':common['all'],
            'noisy_reads':3,'fresh_pointer_initializations':3,'coherent_history_rebits':1,
            'raw_record_bits':3,'decoded_message_bits':3,'new_pure_helpers':0,
            'source_purification_all_pointers_and_baths_retained':True,
            'all_eight_raw_records_accepted':True,'common_interface_before_recovery_unchanged':True,
            'compilation_gate_minimality_proved':False,
            'transport_classical_compute_and_physical_calibration_costs_included':False,
            'limited_pulse_error_certificate':pulse_certificate()}


class CompiledRelationalHistoryTests(unittest.TestCase):
    def test_correction_and_all_six_native_decoders_match_full_operators(self):
        full=np.zeros((16,16));full[:8,:8]=np.eye(8);full[8:,8:]=np.array(correction(),float)
        np.testing.assert_allclose(circuit_matrix(correction_circuit(),4),full,atol=3e-15)
        for t in (F(0),PARAMETER):
            for g in range(6):
                np.testing.assert_allclose(circuit_matrix(decoder_circuit(g,t),4),decoder(g,t),atol=5e-15)

    def test_complete_worst_and_expected_resource_costs(self):
        result=resource_certificate()
        self.assertEqual(result['correction_native_counts'],{'yx':10,'ry':8,'all':18})
        self.assertEqual([c['all'] for c in result['conditional_native_counts']],[51,39,51,41,29,41])
        self.assertEqual((result['new_yx_worst'],result['new_ry_worst'],result['new_native_gates_worst']),(303,128,431))
        self.assertEqual(result['total_including_first_join_worst'],597)
        self.assertAlmostEqual(result['expected_new_native_gates']['diagnostic']-result['expected_previous_native_gates']['diagnostic'],18.,places=12)

    def test_actual_pointer_tree_with_complex_external_system(self):
        rng=np.random.default_rng(170);raw=rng.normal(size=(32,3))+1j*rng.normal(size=(32,3))
        source=raw @ raw.conj().T;source/=np.trace(source)
        expected=channel(tuple(np.kron(k,np.eye(2)) for k in recovery_kraus()),source)
        np.testing.assert_allclose(actual_recovery(source,2),expected,atol=5e-15)

    def test_original_common_interface_still_present_before_recovery(self):
        rng=np.random.default_rng(270);raw=rng.normal(size=(16,3))+1j*rng.normal(size=(16,3))
        source=raw @ raw.conj().T;source/=np.trace(source)
        stored=actual_recovery(source,apply_decoder=False)
        np.testing.assert_allclose(keep_systems(stored,(1,2,3),4),channel(compiled_sharing(F(1,2))[1],source),atol=5e-15)

    def test_four_native_pulse_errors_bound_full_external_output(self):
        eps=5e-5;errors=(eps,-eps,eps,eps)
        ideal=circuit_matrix(correction_circuit(),4)
        noisy=circuit_matrix(correction_circuit(pulse_errors=errors),4)
        self.assertLessEqual(np.linalg.norm(noisy-ideal,2),2*eps+5e-15)
        rng=np.random.default_rng(370);v=rng.normal(size=32)+1j*rng.normal(size=32);v/=np.linalg.norm(v)
        source=np.outer(v,v.conj());p=confusion_matrix(BIT_ERROR)
        perturbed=tuple(np.kron(np.sqrt(p[r][g]/6)*circuit_matrix(decoder_circuit(g,pulse_errors=errors),4) @ conditional_unitary(r+1),np.eye(2)) for r in range(6) for g in range(6))
        perfect=tuple(np.kron(k,np.eye(2)) for k in recovery_kraus())
        self.assertLessEqual(trace_distance(channel(perturbed,source),channel(perfect,source)),2*eps+2e-14)

    def test_limited_calibration_margin_is_strict(self):
        result=pulse_certificate()
        self.assertTrue(result['still_beats_every_ideal_four_angle_scheme'])
        self.assertFalse(result['all_native_gate_noise_covered'])
        self.assertFalse(pulse_certificate(F(1,1000))['still_beats_every_ideal_four_angle_scheme'])

    def test_no_extra_helpers_or_postselection_and_guards(self):
        result=resource_certificate()
        self.assertEqual(result['new_pure_helpers'],0)
        self.assertTrue(result['all_eight_raw_records_accepted'])
        with self.assertRaises(ValueError): decoder_circuit(True)
        with self.assertRaises(ValueError): actual_recovery(np.eye(16)/16,3)
        with self.assertRaises(ValueError): correction_circuit(pulse_errors=(0,))
        with self.assertRaises(ValueError): pulse_certificate(-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CompiledRelationalHistoryTests))
    if not checks.wasSuccessful():raise SystemExit(1)
    report=resource_certificate()
    report['automated_checks']={'run':checks.testsRun,'failures':len(checks.failures),'errors':len(checks.errors)}
    if args.write_results:
        Path(__file__).with_name('compiled_relational_history_results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
