"""Round 150: four history labels and a sharp promised error for a real decoder.

Anticommuting repair pairs admit an orthogonal midpoint. The exact worst
external error of this construction follows from a two-dimensional block;
optimality among all promised-input recovery channels is not established.
"""

import argparse
import json
from pathlib import Path
import unittest

import numpy as np

from certified_intervals import Interval as I
from classical_joining_history import record_unitary, independent_encoding
from coarsened_history_bound import pure_extension_metrics, negative_witness, representative_kraus
from compiled_joining_frontier import compiled_sharing
from compiled_subject_joining import circuit_matrix, controlled_z, intake_circuit, ry
from compiled_three_subject_joining import inverse, gate_counts
from external_correlation_recovery_bound import channel
from flagged_subject_joining import joining_kraus, random_state
from independent_source_alignment import keep_systems
from joint_relation_records import interval_fields
from minimal_coherent_joining_history import conditional_unitary
from fractions import Fraction as F


GROUPS = ((0,2),(3,5),(1,),(4,))


def centered_unitary(group):
    if type(group) is not int or group not in range(4): raise ValueError("Use a group 0..3.")
    labels = GROUPS[group]
    return sum(record_unitary(r+1) for r in labels)/np.sqrt(len(labels))


def centered_decoder(group):
    ks = joining_kraus((2,2))
    return np.vstack((ks[(0,0)], centered_unitary(group) @ ks[(0,1)])).T


def centered_kraus():
    return tuple(centered_decoder(g) @ conditional_unitary(r+1)/np.sqrt(6)
                 for g, labels in enumerate(GROUPS) for r in labels)


def exact_error_interval():
    return (I.exact(1)+(I.exact(25)-16*I.exact(2).sqrt()).sqrt())/12


def centered_decoder_circuit(group):
    centered_unitary(group)
    gates = ()
    if group in (0,2): gates += controlled_z(0,1,2)
    if group in (0,1):
        target = 2+group
        gates += (ry(target,-np.pi/4),)+controlled_z(0,target,1)+(ry(target,np.pi/4),)
    return gates+inverse(intake_circuit(2))


def analytic_kraus():
    kp, km = (joining_kraus((2,2))[key] for key in ((0,0),(0,1)))
    pplus, pminus = kp.T @ kp, km.T @ km
    real_y = np.array([[0.,-1.],[1.,0.]])
    ja = km.T @ np.kron(np.eye(2),np.kron(real_y,np.eye(2))) @ km
    jb = km.T @ np.kron(np.eye(4),real_y) @ km
    ell = pplus+pminus/np.sqrt(2)
    return (np.eye(16)/np.sqrt(3), np.sqrt(2/3)*ell, ja/np.sqrt(6), jb/np.sqrt(6))


def purification(omega):
    eig, vectors = np.linalg.eigh(omega)
    return vectors * np.sqrt(np.maximum(eig,0))


class CenteredHistoryRecoveryTests(unittest.TestCase):
    def test_only_selected_distinct_pairs_anticommute_and_midpoints_are_real_orthogonal(self):
        unitaries = [record_unitary(r) for r in range(1,7)]
        pairs = [(r,s) for r in range(6) for s in range(r+1,6)
                 if np.max(abs(unitaries[r] @ unitaries[s]+unitaries[s] @ unitaries[r]))<1e-14]
        self.assertEqual(pairs,[(0,2),(3,5)])
        for g in range(4): np.testing.assert_allclose(centered_unitary(g).T @ centered_unitary(g),np.eye(8),atol=3e-16)

    def test_midpoint_channel_matches_the_four_term_analytic_channel(self):
        rng = np.random.default_rng(150)
        raw=rng.normal(size=(16,4))+1j*rng.normal(size=(16,4)); state=raw @ raw.conj().T
        np.testing.assert_allclose(channel(centered_kraus(),state),channel(analytic_kraus(),state),atol=6e-15)
        np.testing.assert_allclose(sum(k.T @ k for k in centered_kraus()),np.eye(16),atol=5e-16)

    def test_every_sampled_promised_purification_attains_the_same_exact_error(self):
        expected = sum(exact_error_interval().floats())/2
        rng=np.random.default_rng(250)
        states=[np.eye(16)/16]
        states += [independent_encoding((random_state(rng,2),random_state(rng,2))) for _ in range(12)]
        for state in states:
            c=purification(state)
            self.assertAlmostEqual(pure_extension_metrics(centered_kraus(),c)["trace_error"],expected)
            # A complex isometry on the inaccessible reference cannot change the error.
            raw=rng.normal(size=(20,16))+1j*rng.normal(size=(20,16)); v=np.linalg.qr(raw)[0]
            self.assertAlmostEqual(pure_extension_metrics(centered_kraus(),c @ v.T)["trace_error"],expected)

    def test_two_dimensional_signed_block_gives_the_interval_expression(self):
        off=-(2-np.sqrt(2))/6
        block=np.array([[0.,off],[off,-1/6]])
        error=(np.sum(abs(np.linalg.eigvalsh(block)))+1/6)/2
        self.assertAlmostEqual(error,(1+np.sqrt(25-16*np.sqrt(2)))/12)
        lo,hi=exact_error_interval().floats()
        self.assertLess(abs((lo+hi)/2-error),3e-16)

    def test_promised_gain_does_not_contradict_the_unrestricted_lower_bound(self):
        self.assertAlmostEqual(pure_extension_metrics(centered_kraus(),negative_witness())["trace_error"],1/3)
        c=np.eye(16)/4
        self.assertLess(pure_extension_metrics(centered_kraus(),c)["trace_error"],
                        pure_extension_metrics(representative_kraus(GROUPS),c)["trace_error"])

    def test_all_conditional_decoders_compile_into_the_original_real_gates(self):
        for g in range(4):
            np.testing.assert_allclose(circuit_matrix(centered_decoder_circuit(g),4),centered_decoder(g),atol=4e-15)
        self.assertEqual(gate_counts(centered_decoder_circuit(0)),{"yx":16,"ry":17,"all":33})

    def test_storing_coarse_labels_preserves_the_old_common_interface(self):
        rng=np.random.default_rng(350); raw=rng.normal(size=(16,5)); state=raw @ raw.T; state/=np.trace(state)
        retained=channel(tuple(conditional_unitary(r)/np.sqrt(6) for r in range(1,7)),state)
        old=channel(compiled_sharing(F(1,2))[1],state)
        np.testing.assert_allclose(keep_systems(retained,(1,2,3),4),old,atol=3e-15)
        for group in GROUPS:
            gram=sum(conditional_unitary(r+1).T @ conditional_unitary(r+1)/6 for r in group)
            np.testing.assert_allclose(gram,len(group)*np.eye(16)/6,atol=3e-16)

    def test_marginals_outside_the_promise_need_not_have_the_promised_constant(self):
        km=joining_kraus((2,2))[(0,1)]
        c=km.T[:,0:1]
        self.assertNotAlmostEqual(pure_extension_metrics(centered_kraus(),c)["trace_error"],
                                 sum(exact_error_interval().floats())/2)
        for bad in (-1,4,True):
            with self.assertRaises(ValueError): centered_decoder(bad)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CenteredHistoryRecoveryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":150,"coarse_groups_zero_based":GROUPS,"fixed_message_bits":2,
        "one_accessible_coherent_history_rebit":True,
        "exact_promised_external_error_formula":"(1+sqrt(25-16sqrt(2)))/12",
        "exact_promised_external_error":interval_fields(exact_error_interval()),
        "promise_for_theorem":"real symmetric input marginal, equal orientation weights; arbitrary external extension",
        "sharp_for_this_decoder":True,"optimal_among_all_promised_decoders":False,
        "decoder_gate_counts":[gate_counts(centered_decoder_circuit(g)) for g in range(4)],
        "old_common_marginal_preserved_until_recovery":True,
        "closed_residual_history_is_retained_not_deleted":True,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("centered_history_recovery_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
