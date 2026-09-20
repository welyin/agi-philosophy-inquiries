"""Round 141: sharp external-correlation error with only classical history.

After the fixed balanced join the decoder has dimension 8, the input 16.
Every fine classical-feedback Kraus operator factors through dimension 8.
The rank/trace inequality proves the lower bound; numerical probes only audit
implementations. Independent local auxiliaries are allowed, shared entangled
resources and quantum transfers from the environment are not included.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import trace_distance
from classical_joining_history import history_kraus, recovery_isometry
from encoded_composition_audit import independent_encoding
from flagged_subject_joining import joining_kraus


def dimension_bound(input_dimension, quantum_dimension):
    for value in (input_dimension, quantum_dimension):
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError("Positive integer dimensions are required.")
    overlap = min(F(1), F(quantum_dimension, input_dimension))
    return {"max_entangled_overlap_upper": overlap,
            "worst_external_trace_error_lower": 1-overlap}


def channel(kraus, state):
    return sum(k @ state @ k.conj().T for k in kraus)


def entangled_overlap(kraus):
    d = kraus[0].shape[1]
    if any(k.shape != (d, d) for k in kraus):
        raise ValueError("A channel on the original space is required.")
    return float(sum(abs(np.trace(k))**2 for k in kraus)/d**2)


def max_entangled_metrics(kraus):
    """Exact finite matrix evaluation using a small support; no large Choi array.

    QR reduces X diag(1,...,1,-1) X^dagger without changing nonzero eigenvalues.
    Returned fidelity is the overlap, not its square root.
    """
    d = kraus[0].shape[1]
    columns = [k.ravel()/np.sqrt(d) for k in kraus]
    columns.append(np.eye(d).ravel()/np.sqrt(d))
    _, upper = np.linalg.qr(np.column_stack(columns), mode="reduced")
    signs = np.ones(len(columns)); signs[-1] = -1
    difference = (upper*signs) @ upper.conj().T
    return {"overlap": entangled_overlap(kraus),
            "trace_error": float(np.abs(np.linalg.eigvalsh(difference)).sum()/2)}


def corrected_history():
    return tuple(recovery_isometry(i) @ k for i, k in enumerate(history_kraus()))


def polar_feedback(kraus):
    """Legal isometric feedback for each given rectangular Kraus operator."""
    result = []
    for k in kraus:
        u, _, vh = np.linalg.svd(k, full_matrices=False)
        decoder = vh.conj().T @ u.conj().T
        result.append(decoder @ k)
    return tuple(result)


def mixed_environment_basis(seed, complex_basis=False):
    rng = np.random.default_rng(seed)
    original = history_kraus()
    raw = rng.normal(size=(7, 7))
    if complex_basis:
        raw = raw + 1j*rng.normal(size=(7, 7))
    basis, _ = np.linalg.qr(raw)
    return tuple(sum(basis[i, j]*original[j] for j in range(7)) for i in range(7))


class ExternalCorrelationRecoveryBoundTests(unittest.TestCase):
    def test_rank_trace_inequality_for_complex_non_normal_matrices(self):
        rng = np.random.default_rng(141)
        for rank in (1, 3, 8):
            a = rng.normal(size=(16, rank))+1j*rng.normal(size=(16, rank))
            b = rng.normal(size=(rank, 16))+1j*rng.normal(size=(rank, 16))
            matrix = a @ b
            self.assertLessEqual(abs(np.trace(matrix))**2,
                                 rank*np.linalg.norm(matrix)**2+1e-9)
        for p in (k.T @ k for k in joining_kraus((2, 2)).values()):
            self.assertAlmostEqual(abs(np.trace(p))**2, 8*np.linalg.norm(p)**2)

    def test_actual_classical_correction_is_trace_preserving_with_rank_eight_terms(self):
        ks = corrected_history()
        np.testing.assert_allclose(sum(k.T @ k for k in ks), np.eye(16), atol=2e-15)
        self.assertTrue(all(np.linalg.matrix_rank(k) == 8 for k in ks))
        self.assertAlmostEqual(entangled_overlap(ks), .5)

    def test_alternative_real_and_complex_environment_bases_obey_the_same_bound(self):
        for complex_basis in (False, True):
            for seed in (141, 241, 341):
                ks = polar_feedback(mixed_environment_basis(seed, complex_basis))
                np.testing.assert_allclose(sum(k.conj().T @ k for k in ks), np.eye(16), atol=3e-15)
                self.assertLessEqual(entangled_overlap(ks), .5+2e-15)
                self.assertTrue(all(np.linalg.matrix_rank(k) <= 8 for k in ks))

    def test_purification_witness_has_a_valid_independent_encoded_marginal(self):
        source = independent_encoding((np.eye(2)/2, np.eye(2)/2))
        np.testing.assert_array_equal(source, np.eye(16)/16)
        phi = np.eye(16).ravel()/4
        whole = np.outer(phi, phi)
        np.testing.assert_array_equal(np.trace(whole.reshape(16,16,16,16),axis1=1,axis2=3), source)

    def test_small_support_metrics_match_full_external_state_and_attain_the_bound(self):
        ks = corrected_history()
        phi = np.eye(16).ravel()/4
        original = np.outer(phi, phi)
        output = channel(tuple(np.kron(k,np.eye(16)) for k in ks), original)
        metrics = max_entangled_metrics(ks)
        self.assertAlmostEqual(metrics["overlap"], float(phi @ output @ phi))
        self.assertAlmostEqual(metrics["trace_error"], trace_distance(output, original))
        self.assertAlmostEqual(metrics["trace_error"], .5)

    def test_pinching_equals_a_real_random_reflection_on_arbitrary_external_states(self):
        ps = tuple(k.T @ k for k in joining_kraus((2,2)).values())
        reflection = ps[0]-ps[1]
        np.testing.assert_allclose(reflection @ reflection, np.eye(16), atol=1e-15)
        rng = np.random.default_rng(441)
        for _ in range(4):
            raw = rng.normal(size=(48,5))+1j*rng.normal(size=(48,5))
            source = raw @ raw.conj().T; source /= np.trace(source)
            extended = np.kron(reflection,np.eye(3))
            actual = channel(tuple(np.kron(k,np.eye(3)) for k in corrected_history()), source)
            np.testing.assert_allclose(actual, (source+extended @ source @ extended)/2, atol=1e-16)
            self.assertLessEqual(trace_distance(actual,source), .5+1e-14)

    def test_fresh_decoder_isometry_does_not_increase_kraus_rank(self):
        rng = np.random.default_rng(541)
        larger, _ = np.linalg.qr(rng.normal(size=(32,8)))
        for k in history_kraus():
            self.assertEqual(np.linalg.matrix_rank(larger @ k), 8)
        self.assertEqual(dimension_bound(16,8)["worst_external_trace_error_lower"], F(1,2))

    def test_dimensions_and_guards(self):
        self.assertEqual(dimension_bound(64,16)["max_entangled_overlap_upper"], F(1,4))
        self.assertEqual(dimension_bound(16,32)["worst_external_trace_error_lower"], 0)
        for dims in ((0,8),(16,-1),(16,1.5),(True,8)):
            with self.assertRaises(ValueError): dimension_bound(*dims)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ExternalCorrelationRecoveryBoundTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":141,"input_dimension":16,"accessible_quantum_dimension":8,
        "optimized_over":"all environment measurements and conditional CPTP decoders",
        "optimal_worst_external_trace_error_exact":"1/2",
        "max_entangled_overlap_upper_exact":"1/2",
        "achieved_by":"round 138 classical recovery, a two-sector pinching",
        "attaining_metrics":max_entangled_metrics(corrected_history()),
        "independent_local_decoder_auxiliaries_allowed":True,
        "unbounded_classical_message_alphabet_allowed":True,
        "pre_shared_entangled_helper_resources_included":False,
        "environment_to_decoder_quantum_transfer_allowed":False,
        "source_purification_accessed_by_decoder":False,
        "alternative_bases_are_numerical_audits_not_global_optimization":True,
        "whole_apparatus_retained":True,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("external_correlation_recovery_bound_results.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__ == "__main__": main()
