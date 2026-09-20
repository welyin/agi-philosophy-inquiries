"""Round 145: six classical labels are necessary with one coherent history bit.

The one-dimensional positive-sector Kraus support forces each reversible
16-to-16 branch to induce a unitary on the negative sector. Its reduced Choi
rank is six. This strengthens the previous rank-seven divided-by-two bound.
"""

import argparse
from fractions import Fraction as F
from itertools import combinations
import json
import math
from pathlib import Path
import unittest

import numpy as np

from classical_joining_history import history_kraus, record_unitary
from external_correlation_recovery_bound import channel
from flagged_subject_joining import joining_kraus
from minimal_coherent_joining_history import conditional_unitary


PREFIX_CODES=("00","01","100","101","110","111")


def negative_unitaries():
    return tuple(record_unitary(r) for r in range(1,7))


def integer_choi_numerator():
    columns=np.column_stack([u.astype(np.int64).ravel() for u in negative_unitaries()])
    return columns @ columns.T


def rotated_unitaries(angle):
    result=list(negative_unitaries())
    c,s=math.cos(angle),math.sin(angle)
    for i,j in ((0,2),(3,5)):
        a,b=result[i],result[j]
        result[i]=c*a+s*b; result[j]=-s*a+c*b
    return tuple(result)


def reversible_branch(unitary,quantum_rotation):
    base=joining_kraus((2,2))
    return np.kron(quantum_rotation,np.eye(8)) @ np.vstack((base[(0,0)],unitary @ base[(0,1)]))


class CoherentHistoryMessageMinimumTests(unittest.TestCase):
    def test_negative_unitaries_have_exact_integer_orthogonal_gram(self):
        us=[u.astype(np.int64) for u in negative_unitaries()]
        gram=np.array([[np.sum(a*b) for b in us] for a in us])
        np.testing.assert_array_equal(gram,8*np.eye(6,dtype=np.int64))
        for u in us: np.testing.assert_array_equal(u.T @ u,np.eye(8,dtype=np.int64))

    def test_normalized_negative_choi_has_exact_rank_six_and_spectral_cap(self):
        numerator=integer_choi_numerator()
        np.testing.assert_array_equal(numerator @ numerator,8*numerator)
        self.assertEqual(np.trace(numerator),48)
        self.assertEqual(np.trace(numerator)//8,6)
        # Choi state is numerator / 48, so each nonzero eigenvalue is 1/6.
        self.assertEqual(F(8,48),F(1,6))

    def test_any_old_kraus_support_element_is_scalar_on_positive_sector(self):
        rng=np.random.default_rng(145); basis=history_kraus()
        positive=joining_kraus((2,2))[(0,0)]
        for _ in range(6):
            coefficients=rng.normal(size=7)+1j*rng.normal(size=7)
            k=sum(a*b for a,b in zip(coefficients,basis))
            np.testing.assert_allclose(k @ positive.T,coefficients[0]*np.eye(8),atol=1e-15)

    def test_general_complex_quantum_bit_rotation_preserves_branch_factorization(self):
        rng=np.random.default_rng(245)
        base=joining_kraus((2,2)); kp,km=base[(0,0)],base[(0,1)]
        for u in rotated_unitaries(.37):
            raw=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); rotation,_=np.linalg.qr(raw)
            v=reversible_branch(u,rotation)
            np.testing.assert_allclose(v.conj().T @ v,np.eye(16),atol=1e-15)
            np.testing.assert_allclose(v @ kp.T,np.kron(rotation[:,0:1],np.eye(8)),atol=1e-15)
            np.testing.assert_allclose(v @ km.T,np.kron(rotation[:,1:2],u),atol=1e-15)

    def test_continuous_alternative_six_unitary_ensembles_realize_same_channel(self):
        rng=np.random.default_rng(345)
        raw=rng.normal(size=(8,3))+1j*rng.normal(size=(8,3))
        source=raw @ raw.conj().T; source/=np.trace(source)
        expected=channel(negative_unitaries(),source)/6
        for angle in (.13,.47,1.19):
            us=rotated_unitaries(angle)
            for u in us: np.testing.assert_allclose(u.T @ u,np.eye(8),atol=1e-15)
            np.testing.assert_allclose(channel(us,source)/6,expected,atol=1e-16)

    def test_a_single_label_cannot_exceed_one_sixth_probability(self):
        choi=integer_choi_numerator()/48
        u=rotated_unitaries(.27)[0].ravel()/np.sqrt(8)
        remainder=choi-np.outer(u,u)/6
        self.assertGreaterEqual(np.linalg.eigvalsh(remainder).min(),-1e-16)
        forbidden=choi-(F(1,6)+F(1,1000))*np.outer(u,u)
        self.assertAlmostEqual(float(u @ np.asarray(forbidden,dtype=float) @ u),-.001)

    def test_actual_six_reversible_branches_attain_the_count(self):
        base=joining_kraus((2,2))
        rng=np.random.default_rng(445); raw=rng.normal(size=(16,4)); source=raw @ raw.T; source/=np.trace(source)
        recovered_common=np.zeros((8,8))
        for r in range(1,7):
            v=conditional_unitary(r)
            recovered_common+=sum(k @ source @ k.T for k in (v[:8],v[8:]))/6
        np.testing.assert_allclose(recovered_common,channel(history_kraus(),source),atol=2e-16)
        np.testing.assert_allclose(base[(0,0)] @ base[(0,1)].T,np.zeros((8,8)),atol=1e-16)

    def test_optimal_prefix_code_and_fixed_length_count(self):
        for a,b in combinations(PREFIX_CODES,2):
            self.assertFalse(a.startswith(b) or b.startswith(a))
        self.assertEqual(sum(F(1,2**len(code)) for code in PREFIX_CODES),1)
        self.assertEqual(sum(F(len(code),6) for code in PREFIX_CODES),F(8,3))
        self.assertTrue(2**2<6<=2**3)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CoherentHistoryMessageMinimumTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":145,"fixed_balanced_common_channel":True,
        "accessible_coherent_history_dimension":2,"universal_external_recovery_required":True,
        "minimum_nonzero_classical_labels":6,"minimum_fixed_length_message_bits":3,
        "any_single_label_probability_upper_exact":"1/6",
        "minimum_classical_record_entropy_bits":"log2(6)",
        "minimum_classical_record_entropy_decimal":math.log2(6),
        "minimum_single_trial_binary_prefix_expected_bits_exact":"8/3",
        "prefix_codes":list(PREFIX_CODES),
        "classical_alphabet_minimum_is_a_physical_readout_count":False,
        "lower_bound_allows_complex_operations":True,
        "negative_normalized_choi_nonzero_eigenvalues_exact":["1/6"]*6,
        "more_coherent_memory_or_approximate_recovery_included":False,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("coherent_history_message_minimum_results.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=="__main__": main()
