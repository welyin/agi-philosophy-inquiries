"""Round 139: orientation plus repair-side records still cannot improve 11/36.

The fixed first-stage Pauli outcome is unavailable. The bound includes every
record-controlled CPTP second stage and any fixed independent auxiliaries.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import unittest

import numpy as np

from approximate_joining_optimality import six_probes
from approximate_subject_joining import apply_kraus, deterministic_joining_kraus, trace_distance
from classical_joining_history import AXES, history_kraus
from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from encoded_composition_audit import independent_encoding, product_state
from flagged_subject_joining import random_state
from joining_history_advantage import balanced_pair, common_ab_marginal
from multisubject_joining_optimum import (
    local_spectral_certificate, moment_matrices, repaired_state, reorder_qubits, y_projector,
)
from quantum_interface_audit import PAULI_Y


def coarse_kraus(h, selector):
    if h not in (0,1) or selector not in (0,1): raise ValueError("Use binary records.")
    ks=history_kraus()
    return (ks[0]/np.sqrt(2),) if h==0 else ks[1+3*selector:4+3*selector]


def conditional_pair(states,h,selector):
    coarse_kraus(h,selector)
    targets=list(states[:2])
    if h: targets[selector]=repaired_state(targets[selector])
    return np.kron(*targets)


def conditional_certificate(h):
    if h not in (0,1): raise ValueError("Use a binary orientation record.")
    local_spectral_certificate()
    a,b=F(1,2),F(1,3)
    a_gamma,b_gamma=a/3+F(1,6),b/3+F(1,6)
    aa=a*(a_gamma if h else a); bb=b*(b_gamma if h else b)
    same=max(aa*a,bb*b); opposed=max(aa*b,bb*a)
    return {"h":h,"matched_ab_norm_exact":str(aa),"mismatched_ab_norm_exact":str(bb),
            "aligned_block_norm_exact":str(same),"opposed_block_norm_exact":str(opposed),
            "conditional_support_upper_bound_exact":str(4*(same+opposed)),
            "one_coarse_record_probability_exact":"1/4"}


@lru_cache(maxsize=4)
def conditional_score_matrix(h,selector):
    a,b=(m/6 for m in moment_matrices())
    ag,bg=a/3+np.eye(4)/6,b/3+np.eye(4)/6
    pair_a=[a,a]; pair_b=[b,b]
    if h: pair_a[selector]=ag; pair_b[selector]=bg
    ca,cb=np.kron(*pair_a),np.kron(*pair_b)
    grouped=np.zeros((512,512),dtype=complex)
    for first,second,out in product((0,1),repeat=3):
        refs=product_state((y_projector(first).T,y_projector(second).T,y_projector(out)))
        targets=np.kron(ca if first==out else cb,a if second==out else b)
        grouped+=np.kron(refs,targets)/4
    result=reorder_qubits(grouped,(0,1,3,5,7,2,4,6,8))
    if np.max(abs(result.imag))>1e-14: raise ArithmeticError("Expected a real score.")
    return result.real


def conditional_dual(h):
    cert=conditional_certificate(h)
    same=float(F(cert["aligned_block_norm_exact"]))/4
    opposed=float(F(cert["opposed_block_norm_exact"]))/4
    aligned=(np.eye(4)+np.kron(PAULI_Y,PAULI_Y).real)/2
    return np.kron(same*aligned+opposed*(np.eye(4)-aligned),np.eye(8))


def coarse_join(states):
    ks=deterministic_joining_kraus(4,2)
    return sum(apply_kraus(ks,independent_encoding((conditional_pair(states,h,c),states[2])))/4
               for h,c in product((0,1),repeat=2))


class CoarseJoiningHistoryBoundTests(unittest.TestCase):
    def test_actual_environment_grouping_and_conditional_states(self):
        rng=np.random.default_rng(139)
        for _ in range(5):
            states=[random_state(rng,2) for _ in range(2)]
            source=independent_encoding(states)
            for h,c in product((0,1),repeat=2):
                branch=apply_kraus(coarse_kraus(h,c),source)
                self.assertAlmostEqual(np.trace(branch),.25)
                np.testing.assert_allclose(branch,encode_state(conditional_pair(states,h,c))/4,atol=3e-16)

    def test_exact_tensor_bounds_sum_to_the_same_no_history_bound(self):
        first,second=conditional_certificate(0),conditional_certificate(1)
        self.assertEqual(first["conditional_support_upper_bound_exact"],"5/6")
        self.assertEqual(second["conditional_support_upper_bound_exact"],"5/9")
        total=(F(first["conditional_support_upper_bound_exact"])+F(second["conditional_support_upper_bound_exact"]))/2
        self.assertEqual(total,F(25,36))

    def test_actual_preparations_match_the_nontrivial_record_score_matrix(self):
        direct=np.zeros((512,512))
        for states in product(six_probes(),repeat=3):
            source=independent_encoding((conditional_pair(states,1,0),states[2]))
            direct+=np.kron(source.T,real_lift(product_state(states)))/216
        np.testing.assert_allclose(direct,conditional_score_matrix(1,0),atol=4e-17)

    def test_duals_bound_every_record_dependent_decoder_and_are_saturated(self):
        ks=deterministic_joining_kraus(4,2)
        for h,c in product((0,1),repeat=2):
            r,y=conditional_score_matrix(h,c),conditional_dual(h)
            slack=np.kron(y,np.eye(16))-r
            self.assertGreaterEqual(np.linalg.eigvalsh(slack).min(),-3e-16)
            for k in ks: np.testing.assert_allclose(slack@k.T.ravel(),0.,atol=2e-16)
            self.assertAlmostEqual(np.trace(y),float(F(conditional_certificate(h)["conditional_support_upper_bound_exact"])))

    def test_all_216_product_probe_errors_attain_eleven_over_thirty_six(self):
        for states in product(six_probes(),repeat=3):
            actual=coarse_join(states)
            self.assertAlmostEqual(trace_distance(actual,encode_state(product_state(states))),11/36)

    def test_old_joint_interface_is_preserved_after_ignoring_or_using_these_records(self):
        rng=np.random.default_rng(239)
        for _ in range(8):
            states=[random_state(rng,2) for _ in range(3)]
            output=coarse_join(states)
            expected=encode_state(balanced_pair(*states[:2]))
            np.testing.assert_allclose(common_ab_marginal(output),expected,atol=5e-16)

    def test_coarse_records_cannot_restore_any_old_antipodal_axis_perfectly(self):
        for site,axis in product(range(2),range(3)):
            pair=[]
            for sign in (-1,1):
                states=[np.eye(2)/2,np.eye(2)/2]
                states[site]=(np.eye(2)+sign*AXES[axis])/2
                pair.append(independent_encoding(states))
            self.assertAlmostEqual(trace_distance(*pair),1.)
            value=sum(trace_distance(apply_kraus(coarse_kraus(h,c),pair[0]),
                                      apply_kraus(coarse_kraus(h,c),pair[1])) for h,c in product((0,1),repeat=2))
            self.assertAlmostEqual(value,5/6)

    def test_record_ignorance_is_a_valid_coarsening_not_physical_erasure(self):
        rng=np.random.default_rng(339)
        source=independent_encoding([random_state(rng,2) for _ in range(2)])
        output=sum(apply_kraus(coarse_kraus(h,c),source) for h,c in product((0,1),repeat=2))
        np.testing.assert_allclose(output,apply_kraus(history_kraus(),source),atol=4e-16)
        with self.assertRaises(ValueError): coarse_kraus(2,0)
        with self.assertRaises(ValueError): conditional_certificate(-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CoarseJoiningHistoryBoundTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":139,"available_records":["orientation h","repair-side selector c"],
        "unavailable_records":["Pauli outcome E1,E2"],
        "conditional_certificates":[conditional_certificate(h) for h in (0,1)],
        "all_record_controlled_second_stage_support_bound_exact":"25/36",
        "sharp_worst_joint_error_exact":"11/36",
        "old_ab_preservation_required_for_lower_bound":False,
        "all_coarsenings_of_h_and_c_share_same_optimum":True,
        "claim_about_every_possible_two_bit_summary":False,
        "full_classical_history_can_do_better":True,
        "environment_physically_deleted":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("coarse_joining_history_bound_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
