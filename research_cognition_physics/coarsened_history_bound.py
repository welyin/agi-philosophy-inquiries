"""Round 149: exact universal error for coarsening the fixed six-label history.

This minimax theorem ranges over all input states and external extensions.
Its lower-bound witness is outside the independent-encoding promise. It is
therefore a benchmark, not a lower bound for the original promised task.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import numpy as np

from classical_joining_history import record_unitary
from flagged_subject_joining import joining_kraus
from minimal_coherent_joining_history import conditional_unitary


def partitions(items=tuple(range(6))):
    if not items:
        yield ()
        return
    first, *rest = items
    for partition in partitions(tuple(rest)):
        yield ((first,),) + partition
        for index in range(len(partition)):
            yield partition[:index] + ((first,) + partition[index],) + partition[index+1:]


def validate_partition(groups):
    flat = [r for group in groups for r in group]
    if (not groups or any(not group for group in groups)
            or any(type(r) is not int for r in flat) or sorted(flat) != list(range(6))):
        raise ValueError("Partition the six integer labels 0..5 into nonempty groups.")


def representative_kraus(groups):
    validate_partition(groups)
    return tuple(conditional_unitary(group[0]+1).T @ conditional_unitary(r+1)/np.sqrt(6)
                 for group in groups for r in group)


def pure_extension_metrics(kraus, coefficient):
    """Coefficient matrix C represents vec(C) on input x inaccessible reference."""
    columns = np.column_stack([(k @ coefficient).reshape(-1) for k in kraus]
                              + [coefficient.reshape(-1)])
    _, reduced = np.linalg.qr(columns, mode="reduced")
    signs = np.array([1] * len(kraus) + [-1])
    difference = (reduced * signs) @ reduced.conj().T
    overlap = sum(abs(np.vdot(coefficient, k @ coefficient))**2 for k in kraus)
    return {"overlap": float(overlap),
            "trace_error": float(np.sum(np.abs(np.linalg.eigvalsh(difference)))/2)}


def negative_witness():
    return joining_kraus((2,2))[(0,1)].T/np.sqrt(8)


def decoder_target_overlap(groups, decoders):
    """Independent numerical audit of the Bessel lower-bound calculation."""
    validate_partition(groups)
    target = joining_kraus((2,2))[(0,1)].T
    injection = np.vstack((np.zeros((8,8)), np.eye(8)))
    return float(sum(abs(np.trace(target.T @ decoder @ injection @ record_unitary(r+1)))**2
                     for group, decoder_ks in zip(groups, decoders)
                     for decoder in decoder_ks for r in group)/(6*64))


class CoarsenedHistoryBoundTests(unittest.TestCase):
    def test_six_negative_repairs_are_hilbert_schmidt_orthogonal(self):
        unitaries = [record_unitary(r) for r in range(1,7)]
        gram = [[np.trace(a.T @ b) for b in unitaries] for a in unitaries]
        np.testing.assert_allclose(gram, 8*np.eye(6), atol=1e-15)

    def test_all_203_partitions_are_unique_and_cover_every_label(self):
        values = list(partitions())
        self.assertEqual(len(values),203)
        self.assertEqual(len(set(values)),203)
        for value in values: validate_partition(value)

    def test_representative_decoders_are_trace_preserving(self):
        for groups in partitions():
            ks = representative_kraus(groups)
            np.testing.assert_allclose(sum(k.T @ k for k in ks),np.eye(16),atol=8e-16)

    def test_every_partition_attains_the_universal_lower_witness(self):
        for groups in partitions():
            metric = pure_extension_metrics(representative_kraus(groups),negative_witness())
            self.assertAlmostEqual(metric["overlap"],len(groups)/6)
            self.assertAlmostEqual(metric["trace_error"],1-len(groups)/6)

    def test_arbitrary_complex_external_states_obey_the_constructive_upper_bound(self):
        rng = np.random.default_rng(149)
        for count in range(1,7):
            groups = tuple((r,) for r in range(count-1)) + (tuple(range(count-1,6)),)
            c = rng.normal(size=(16,5))+1j*rng.normal(size=(16,5)); c /= np.linalg.norm(c)
            self.assertLessEqual(pure_extension_metrics(representative_kraus(groups),c)["trace_error"],
                                 1-count/6+2e-14)

    def test_bessel_bound_allows_general_complex_multi_kraus_decoders(self):
        rng = np.random.default_rng(249)
        for groups in (((0,1,2),(3,4,5)),((0,2),(3,5),(1,),(4,))):
            decoders = []
            for _ in groups:
                raw = rng.normal(size=(48,16))+1j*rng.normal(size=(48,16))
                q = np.linalg.qr(raw)[0]
                decoders.append(tuple(q[i:i+16] for i in range(0,48,16)))
            self.assertLessEqual(decoder_target_overlap(groups,decoders),len(groups)/6+1e-14)

    def test_lower_witness_is_explicitly_outside_the_original_source_promise(self):
        c = negative_witness(); omega = c @ c.T
        kp = joining_kraus((2,2))[(0,0)]
        self.assertAlmostEqual(np.trace(omega),1.)
        self.assertAlmostEqual(np.trace(kp @ omega @ kp.T),0.)
        self.assertNotAlmostEqual(np.trace(kp @ omega @ kp.T),.5)

    def test_partition_guards_reject_omissions_repeats_and_empty_groups(self):
        for bad in ((),((),tuple(range(6))),((0,1,2),(3,4,4)),((0,1,2),(3,4)),((False,1,2,3,4,5),)):
            with self.assertRaises(ValueError): representative_kraus(bad)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CoarsenedHistoryBoundTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":149,"fixed_six_history_partition_count":203,
        "unrestricted_input_and_external_minimax_error":"1-K/6",
        "frontier":[{"labels":k,"exact_error":str(F(6-k,6))} for k in range(1,7)],
        "all_conditional_CPTP_decoders_covered_by_lower_bound":True,
        "lower_witness_satisfies_independent_encoding_promise":False,
        "lower_bound_applies_to_arbitrary_new_environment_instruments":False,
        "fixed_interface_and_one_coherent_history_bit":True,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("coarsened_history_bound_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
