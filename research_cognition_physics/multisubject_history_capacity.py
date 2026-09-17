"""Round 142: many-subject external recovery and coherent memory capacity.

Full classical history of the round 135 channel gives orientation pinching.
The separate memory tradeoff optimizes a general classical/quantum encoding;
its attaining flagged representation is not claimed to retain the fixed
optimal common-interface marginal while the memory is stored.
"""

import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import trace_distance
from encoded_composition_audit import independent_encoding
from external_correlation_recovery_bound import channel, dimension_bound, max_entangled_metrics
from flagged_subject_joining import joining_kraus, orientation_flags, random_state
from multisubject_joining_optimum import multisubject_kraus, orientation_plan, validate_count


def capacity(count, coherent_flag_bits=0):
    n=validate_count(count)
    s=coherent_flag_bits
    if not isinstance(s,int) or isinstance(s,bool) or not 0 <= s <= n-1:
        raise ValueError("Keep an integer number of flags between zero and n-1.")
    d=4**n; q=2**(n+1); k=2**s
    return {"subjects":n,"input_dimension":d,"common_dimension":q,
            "coherent_flag_bits":s,"quantum_memory_dimension":q*k,
            "classical_orientation_bits":n-1-s,
            "sharp_external_trace_error":dimension_bound(d,q*k)["worst_external_trace_error_lower"]}


def orientation_projectors(count):
    validate_count(count)
    return tuple(k.T @ k for k in joining_kraus((2,)*count).values())


def grouped_projectors(count, coherent_flag_bits=0):
    capacity(count,coherent_flag_bits)
    projectors=orientation_projectors(count)
    group_size=2**coherent_flag_bits
    return tuple(sum(projectors[j:j+group_size]) for j in range(0,len(projectors),group_size))


def sign_reflections(projectors):
    """Walsh signs: a real unitary mixture giving the same block pinching."""
    m=len(projectors)
    if m < 1 or m & (m-1): raise ValueError("Use a power-of-two number of sectors.")
    return tuple(sum((-1)**((s & b).bit_count())*p for b,p in enumerate(projectors))
                 for s in range(m))


def corrected_optimal_history(count):
    """Recover after every Kraus record of the actual round 135 channel."""
    ks=multisubject_kraus(count)
    result=[]; index=0
    for flag in orientation_flags(count):
        _,wrong=orientation_plan(flag)
        weight=F(1,3**len(wrong))
        for _ in product(range(3),repeat=len(wrong)):
            k=ks[index]; index+=1
            decoder=k.T/np.sqrt(float(weight))
            result.append((flag,weight,k,decoder))
    if index != len(ks): raise ArithmeticError("History ordering mismatch.")
    return tuple(result)


class MultisubjectHistoryCapacityTests(unittest.TestCase):
    def test_orientation_sectors_are_complete_orthogonal_equal_rank_projections(self):
        for n in (1,2,3):
            ps=orientation_projectors(n); d=4**n; q=2**(n+1)
            np.testing.assert_allclose(sum(ps),np.eye(d),atol=1e-15)
            for i,p in enumerate(ps):
                self.assertAlmostEqual(np.trace(p),q)
                for j,r in enumerate(ps):
                    np.testing.assert_allclose(p @ r,p if i==j else np.zeros((d,d)),atol=1e-15)

    def test_actual_optimal_join_history_recovers_exactly_the_sector_pinching(self):
        rng=np.random.default_rng(142)
        for n in (2,3):
            d=4**n; q=2**(n+1)
            raw=rng.normal(size=(d,4))+1j*rng.normal(size=(d,4))
            source=raw @ raw.conj().T; source/=np.trace(source)
            combined=[]
            for flag,weight,k,decoder in corrected_optimal_history(n):
                np.testing.assert_allclose(decoder.T @ decoder,np.eye(q),atol=1e-15)
                expected=joining_kraus((2,)*n)[flag]
                np.testing.assert_allclose(decoder @ k,np.sqrt(float(weight))*expected.T @ expected,atol=1e-15)
                combined.append(decoder @ k)
            np.testing.assert_allclose(channel(combined,source),channel(orientation_projectors(n),source),atol=1e-16)

    def test_every_memory_split_preserves_promised_unknown_marginals(self):
        rng=np.random.default_rng(242)
        for n in (2,3):
            source=independent_encoding([random_state(rng,2) for _ in range(n)])
            for s in range(n):
                np.testing.assert_allclose(channel(grouped_projectors(n,s),source),source,atol=1e-16)

    def test_real_walsh_mixture_proves_a_uniform_external_error_upper_bound(self):
        rng=np.random.default_rng(342)
        n=3; d=4**n
        raw=rng.normal(size=(d,5)); source=raw @ raw.T; source/=np.trace(source)
        for s in range(n):
            ps=grouped_projectors(n,s); reflections=sign_reflections(ps)
            np.testing.assert_allclose(reflections[0],np.eye(d),atol=1e-15)
            for q in reflections: np.testing.assert_allclose(q @ q,np.eye(d),atol=2e-15)
            actual=channel(ps,source)
            expected=sum(q @ source @ q.T for q in reflections)/len(reflections)
            np.testing.assert_allclose(actual,expected,atol=1e-16)
            self.assertLessEqual(trace_distance(actual,source),float(capacity(n,s)["sharp_external_trace_error"])+1e-14)

    def test_max_entangled_witness_attains_all_small_memory_tradeoffs(self):
        for n in (1,2,3):
            for s in range(n):
                metrics=max_entangled_metrics(grouped_projectors(n,s))
                error=capacity(n,s)["sharp_external_trace_error"]
                self.assertAlmostEqual(metrics["trace_error"],float(error))
                self.assertAlmostEqual(metrics["overlap"],float(1-error))

    def test_pure_external_extension_is_bounded_for_partial_coherent_storage(self):
        rng=np.random.default_rng(442)
        raw=rng.normal(size=(48,2))+1j*rng.normal(size=(48,2))
        source=raw @ raw.conj().T; source/=np.trace(source)
        for s in (0,1):
            ks=tuple(np.kron(p,np.eye(3)) for p in grouped_projectors(2,s))
            self.assertLessEqual(trace_distance(channel(ks,source),source),float(capacity(2,s)["sharp_external_trace_error"])+1e-14)

    def test_exact_dimension_ledger_scales_without_constructing_large_arrays(self):
        expected=(F(0),F(1,2),F(3,4),F(7,8),F(15,16))
        for n,error in enumerate(expected,1):
            self.assertEqual(capacity(n)["sharp_external_trace_error"],error)
        for n in range(2,21):
            self.assertEqual(capacity(n,n-1)["sharp_external_trace_error"],0)
            self.assertEqual(capacity(n,n-2)["sharp_external_trace_error"],F(1,2))

    def test_guards_and_difference_from_common_output_error(self):
        from multisubject_joining_optimum import optimal_error
        self.assertEqual(optimal_error(3),F(1,4))
        self.assertEqual(capacity(3)["sharp_external_trace_error"],F(3,4))
        for args in ((0,0),(3,3),(3,-1),(3,1.5),(3,True)):
            with self.assertRaises(ValueError): capacity(*args)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MultisubjectHistoryCapacityTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    rows=[]
    for n in range(1,7):
        for s in range(n):
            row=capacity(n,s); row["sharp_external_trace_error"]=str(row["sharp_external_trace_error"])
            rows.append(row)
    report={"round":142,"classical_history_optimal_external_error_formula":"1 - 2^(1-n)",
        "actual_round_135_full_classical_history_attains_bound":True,
        "general_memory_tradeoff_formula":"1 - 2^(s-(n-1))",
        "general_memory_tradeoff_scope":"arbitrary classical/quantum encoding with quantum dimension 2^(n+1+s)",
        "minimum_coherent_flag_bits_for_exact_universal_recovery":"n-1",
        "flagged_tradeoff_preserves_the_fixed_optimal_common_marginal":False,
        "simultaneous_optimal_common_interface_and_minimal_quantum_history_solved":False,
        "unrestricted_ideal_control_assumed_for_general_n":True,
        "native_gate_compilation_for_general_n_claimed":False,
        "shared_entangled_auxiliaries_assumed_free":False,"rows":rows,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("multisubject_history_capacity_results.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__ == "__main__": main()
