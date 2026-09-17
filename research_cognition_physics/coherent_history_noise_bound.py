"""Round 143: sharp recovery error for inaccessible flag-dephasing baths.

The new storage coherence lambdas are declared noise parameters, not the old
readout contrast eta. All joining history remains accessible coherently; only
the new baths are closed. All CPTP recovery channels obey the overlap bound.
No redundancy was encoded before the noise in this task.
"""

import argparse
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import trace_distance
from compiled_joining_frontier import compiled_sharing
from compiled_subject_joining import circuit_matrix, controlled_ry
from encoded_composition_audit import independent_encoding
from external_correlation_recovery_bound import channel, entangled_overlap, max_entangled_metrics
from flagged_subject_joining import intake_matrix, random_state
from multisubject_history_capacity import orientation_projectors, sign_reflections


def coherence_values(values):
    result=[]
    for value in values:
        if isinstance(value,bool): raise ValueError("A coherence in [0,1] is required.")
        try: value=F(value)
        except (TypeError,ValueError,OverflowError):
            raise ValueError("Finite rational coherence values are required.") from None
        if not 0 <= value <= 1: raise ValueError("A coherence in [0,1] is required.")
        result.append(value)
    return tuple(result)


def phase_probabilities(values):
    lambdas=coherence_values(values)
    result=[]
    for bits in product((0,1),repeat=len(lambdas)):
        p=F(1)
        for bit,value in zip(bits,lambdas): p *= (1+(-1)**bit*value)/2
        result.append(p)
    return tuple(result)


def optimal_external_error(values):
    return 1-phase_probabilities(values)[0]


def phase_kraus(values):
    lambdas=coherence_values(values)
    reflections=sign_reflections(orientation_projectors(len(lambdas)+1))
    return tuple(np.sqrt(float(p))*q for p,q in zip(phase_probabilities(lambdas),reflections))


def physical_two_subject_bath(value):
    """Original Ry/YX bath coupling on h after the actual seven-wire join.

    Returns the two subnormalized isometries conditioned on a bath basis.
    These labels are calculation indices and are not available to recovery.
    """
    value=float(coherence_values((value,))[0])
    v=compiled_sharing(F(1,2))[0]
    initial=np.zeros((256,16)); initial[::2]=v
    coupled=circuit_matrix(controlled_ry(0,7,2*math.acos(value)),8,initial)
    return tuple(coupled.reshape(128,2,16)[:,b,:] for b in (0,1))


class CoherentHistoryNoiseBoundTests(unittest.TestCase):
    def test_independent_phase_weights_are_exact_and_identity_is_most_likely(self):
        for values in ((),(F(0),),(F(1),),(F(9,10),),(F(1,2),F(4,5))):
            weights=phase_probabilities(values)
            self.assertEqual(sum(weights),1)
            self.assertEqual(max(weights),weights[0])
            self.assertEqual(optimal_external_error(values),1-weights[0])
        self.assertEqual(optimal_external_error((F(9,10),)*2),F(39,400))

    def test_native_bath_coupling_after_actual_balanced_join_matches_declared_noise(self):
        v=compiled_sharing(F(1,2))[0]
        rng=np.random.default_rng(143)
        raw=rng.normal(size=(16,3))+1j*rng.normal(size=(16,3))
        source=raw @ raw.conj().T; source/=np.trace(source)
        for value in (F(0),F(1,3),F(9,10),F(1)):
            branches=physical_two_subject_bath(value)
            decoded=tuple(v.T @ k for k in branches)
            for k,b in zip(decoded,branches): np.testing.assert_allclose(v @ k,b,atol=3e-14)
            np.testing.assert_allclose(channel(decoded,source),channel(phase_kraus((value,)),source),atol=7e-15)
            np.testing.assert_allclose(sum(k.T @ k for k in branches),np.eye(16),atol=2e-14)

    def test_distinct_phase_errors_are_hilbert_schmidt_orthogonal(self):
        for n in (2,3):
            qs=sign_reflections(orientation_projectors(n)); d=4**n
            gram=np.array([[np.trace(a.T @ b) for b in qs] for a in qs])
            np.testing.assert_allclose(gram,d*np.eye(len(qs)),atol=5e-14)

    def test_arbitrary_complex_cptp_recovery_obeys_the_bessel_overlap_bound(self):
        rng=np.random.default_rng(243); d=16
        for _ in range(4):
            raw=rng.normal(size=(2*d,d))+1j*rng.normal(size=(2*d,d))
            isometry,_=np.linalg.qr(raw)
            recovery=(isometry[:d],isometry[d:])
            np.testing.assert_allclose(sum(k.conj().T @ k for k in recovery),np.eye(d),atol=1e-15)
            for value in (F(0),F(1,2),F(9,10)):
                combined=tuple(r @ k for r in recovery for k in phase_kraus((value,)))
                self.assertLessEqual(entangled_overlap(combined),float((1+value)/2)+1e-14)

    def test_do_nothing_after_coherent_inverse_attains_the_bound_on_external_witness(self):
        for values in ((F(0),),(F(1),),(F(9,10),),(F(1,2),F(4,5))):
            metrics=max_entangled_metrics(phase_kraus(values))
            error=optimal_external_error(values)
            self.assertAlmostEqual(metrics["overlap"],float(1-error))
            self.assertAlmostEqual(metrics["trace_error"],float(error))

    def test_all_promised_marginals_are_unchanged_despite_external_error(self):
        rng=np.random.default_rng(343)
        for n in (2,3):
            source=independent_encoding([random_state(rng,2) for _ in range(n)])
            for value in (F(0),F(1,2),F(9,10)):
                np.testing.assert_allclose(channel(phase_kraus((value,)*(n-1)),source),source,atol=2e-16)

    def test_coherence_multiplies_by_hamming_distance_in_actual_flagged_layout(self):
        n=3; w=intake_matrix((2,)*n); q=2**(n+1)
        rng=np.random.default_rng(443)
        raw=rng.normal(size=(64,2)); source=raw @ raw.T; source/=np.trace(source)
        value=F(3,5)
        before=w @ source @ w.T
        after=w @ channel(phase_kraus((value,value)),source) @ w.T
        for a,b in product(range(4),repeat=2):
            block=np.s_[a*q:(a+1)*q,b*q:(b+1)*q]
            np.testing.assert_allclose(after[block],float(value**((a ^ b).bit_count()))*before[block],atol=1e-16)

    def test_uniform_external_bound_and_noise_parameter_guards(self):
        rng=np.random.default_rng(543)
        vector=rng.normal(size=48)+1j*rng.normal(size=48); vector/=np.linalg.norm(vector)
        source=np.outer(vector,vector.conj())
        value=F(4,5)
        output=channel(tuple(np.kron(k,np.eye(3)) for k in phase_kraus((value,))),source)
        self.assertLessEqual(trace_distance(output,source),float(optimal_external_error((value,)))+1e-14)
        for value in (-1,2,float("inf"),float("nan"),True):
            with self.assertRaises(ValueError): coherence_values((value,))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CoherentHistoryNoiseBoundTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    rows=[]
    for n in (2,3,5):
        for value in (F(0),F(1,2),F(9,10),F(99,100),F(1)):
            error=optimal_external_error((value,)*(n-1))
            rows.append({"subjects":n,"per_flag_coherence":str(value),
                         "optimal_external_trace_error_exact":str(error),"error_decimal":float(error)})
    report={"round":143,"optimal_error_formula":"1 - product_j ((1+lambda_j)/2)",
        "bound_optimizes_all_CPTP_recovery_after_noise":True,
        "full_old_joining_history_accessible_coherently":True,
        "new_dephasing_bath_accessible":False,
        "whole_including_new_bath_retained":True,
        "noise_is_the_original_readout_contrast_eta":False,
        "redundant_encoding_before_noise_included":False,
        "two_subject_bath_coupling_original_gates":{"pure_bath_rebits":1,"yx":1,"ry":3},
        "promised_independent_input_marginals_unchanged":True,"rows":rows,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("coherent_history_noise_bound_results.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__ == "__main__": main()
