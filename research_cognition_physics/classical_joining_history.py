"""Round 138: classical Pauli history restores the promised independent inputs.

Recovery is an orientation pinching on unrestricted inputs, not identity on all
external purifications. Seven canonical records suffice and cannot be merged
losslessly when postprocessing this specific fine history.
"""

import argparse
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import apply_kraus, trace_distance
from common_orientation_structure import encode_state
from compiled_joining_frontier import compiled_sharing
from compiled_subject_joining import circuit_matrix, intake_circuit, ry, yx
from compiled_three_subject_joining import inverse, gate_counts
from encoded_composition_audit import independent_encoding, product_state
from flagged_subject_joining import joining_kraus, random_state
from joining_history_advantage import common_ab_marginal, balanced_pair
from multisubject_joining_optimum import multisubject_kraus, reorder_qubits
from quantum_interface_audit import PAULI_X, PAULI_Y, PAULI_Z


AXES = (PAULI_X, PAULI_Y, PAULI_Z)
NAMES = ("none", "A_X", "A_Y", "A_Z", "B_X", "B_Y", "B_Z")
WEIGHTS = (F(1,2),)+(F(1,12),)*6
PREFIX_CODES = ("0", "100", "101", "1100", "1101", "1110", "1111")


def physical_record_class(bits):
    h, e1, e2, selector = bits
    if any(bit not in (0,1) for bit in bits):
        raise ValueError("Use four binary environment labels.")
    if h == 0:
        if (e1,e2) != (0,0):
            raise ValueError("This fine label is outside the promised history support.")
        return 0
    axis = {(0,0):1, (1,0):2, (0,1):0}.get((e1,e2))
    if axis is None:
        raise ValueError("This fine label is outside the promised history support.")
    return 1+3*selector+axis


def record_unitary(record):
    if record not in range(7):
        raise ValueError("Canonical record must be in 0..6.")
    if record == 0:
        return np.eye(8)
    site, axis = divmod(record-1, 3)
    pauli = (PAULI_Z.real, np.eye(2), PAULI_X.real)[axis]
    if site == 0:
        return np.kron(PAULI_Z.real, np.kron(pauli, np.eye(2)))
    return np.kron(np.eye(4), pauli)


def history_kraus():
    base = joining_kraus((2,2))
    return (base[(0,0)],)+tuple(record_unitary(r)@base[(0,1)]/np.sqrt(6) for r in range(1,7))


def recovery_isometry(record):
    base = joining_kraus((2,2))[(0,int(record != 0))]
    return base.T@record_unitary(record).T


def local_pauli_gates(site, axis, helper):
    if axis == "X":
        return (yx(helper, site, np.pi), ry(helper, -np.pi))
    if axis == "Z":
        return (ry(site, np.pi/2), yx(helper, site, np.pi),
                ry(site, -np.pi/2), ry(helper, -np.pi))
    if axis == "I":
        return ()
    raise ValueError("Use real I,X,Z.")


def recovery_circuit(record):
    record_unitary(record)
    gates = ()
    if record:
        site, axis = divmod(record-1, 3)
        gates += (ry(0, np.pi),)  # Fresh reference flag |0> -> |1>.
        if site == 0:
            gates += local_pauli_gates(1, "Z", 2)
        gates += local_pauli_gates(2+site, ("Z","I","X")[axis], 1)
    return gates + inverse(intake_circuit(2))  # Unknown A is borrowed and returned.


def recovered_input(source):
    return sum(recovery_isometry(r)@k@source@k.T@recovery_isometry(r).T
               for r,k in enumerate(history_kraus()))


def reflected(state, axis):
    """Description of a conditional state, not a freely available CP reflection."""
    return state-np.trace(state@AXES[axis]).real*AXES[axis]


def conditional_state(states, record):
    first, second = states
    if record:
        site, axis = divmod(record-1,3)
        if site == 0:
            first = reflected(first, axis)
        else:
            second = reflected(second, axis)
    return encode_state(np.kron(first,second))


def coarse_distance(groups, source_plus, source_minus):
    ks = history_kraus()
    return sum(trace_distance(apply_kraus([ks[r] for r in group],source_plus),
                              apply_kraus([ks[r] for r in group],source_minus)) for group in groups)


class ClassicalJoiningHistoryTests(unittest.TestCase):
    def test_seven_record_instrument_is_the_actual_dephased_compiled_environment(self):
        physical = compiled_sharing(F(1,2))[1]
        rng = np.random.default_rng(138)
        raw = rng.normal(size=(16,16)); source=raw@raw.T; source/=np.trace(source)
        grouped = [np.zeros((8,8)) for _ in range(7)]
        for i,k in enumerate(physical):
            if np.linalg.norm(k)>1e-12:
                bits=tuple((i>>(3-j))&1 for j in range(4))
                grouped[physical_record_class(bits)] += k@source@k.T
        ks=history_kraus()
        np.testing.assert_allclose(grouped,[k@source@k.T for k in ks],atol=2e-15)
        np.testing.assert_allclose(sum(k.T@k for k in ks),np.eye(16),atol=1e-15)

    def test_record_probabilities_are_state_independent_and_reflection_labels_are_correct(self):
        rng=np.random.default_rng(238)
        for _ in range(8):
            states=[random_state(rng,2) for _ in range(2)]
            source=independent_encoding(states)
            for r,k in enumerate(history_kraus()):
                branch=k@source@k.T
                self.assertAlmostEqual(np.trace(branch),float(WEIGHTS[r]))
                np.testing.assert_allclose(branch/float(WEIGHTS[r]),conditional_state(states,r),atol=7e-16)

    def test_recovery_is_exact_pinching_and_identity_on_promised_inputs(self):
        base=joining_kraus((2,2)); ps=[k.T@k for k in base.values()]
        rng=np.random.default_rng(338)
        raw=rng.normal(size=(16,16)); unrestricted=raw@raw.T; unrestricted/=np.trace(unrestricted)
        np.testing.assert_allclose(recovered_input(unrestricted),sum(p@unrestricted@p for p in ps),atol=5e-16)
        for _ in range(8):
            source=independent_encoding([random_state(rng,2) for _ in range(2)])
            np.testing.assert_allclose(recovered_input(source),source,atol=5e-16)

    def test_each_conditional_recovery_is_compiled_into_original_real_gates(self):
        initial=np.eye(16)[:,:8]
        for r in range(7):
            actual=circuit_matrix(recovery_circuit(r),4,initial)
            expected=recovery_isometry(r)
            np.testing.assert_allclose(actual,expected,atol=2e-15)
            np.testing.assert_allclose(expected.T@expected,np.eye(8),atol=5e-16)
        self.assertEqual(max(gate_counts(recovery_circuit(r))["yx"] for r in range(7)),8)
        self.assertEqual(max(gate_counts(recovery_circuit(r))["ry"] for r in range(7)),12)

    def test_classical_recovery_then_three_party_join_matches_the_optimal_output(self):
        rng=np.random.default_rng(438)
        for _ in range(6):
            states=[random_state(rng,2) for _ in range(3)]
            restored=recovered_input(independent_encoding(states[:2]))
            source=reorder_qubits(np.kron(restored,encode_state(states[2])),(0,1,4,2,3,5))
            output=apply_kraus(multisubject_kraus(3),source)
            expected=apply_kraus(multisubject_kraus(3),independent_encoding(states))
            np.testing.assert_allclose(output,expected,atol=5e-16)
            np.testing.assert_allclose(common_ab_marginal(output),encode_state(balanced_pair(*states[:2])),atol=6e-16)

    def test_every_pairwise_record_merge_loses_an_old_perfect_distinction(self):
        ks=history_kraus()
        for left,right in combinations(range(7),2):
            chosen=right if left==0 else left
            site,axis=divmod(chosen-1,3)
            states=[np.eye(2)/2,np.eye(2)/2]
            pair=[]
            for sign in (-1,1):
                trial=list(states); trial[site]=(np.eye(2)+sign*AXES[axis])/2
                pair.append(independent_encoding(trial))
            self.assertAlmostEqual(trace_distance(*pair),1.)
            groups=[(left,right)]+[(r,) for r in range(7) if r not in (left,right)]
            self.assertAlmostEqual(coarse_distance(groups,*pair),5/6)

    def test_classical_recovery_does_not_preserve_arbitrary_external_purification(self):
        vector=np.eye(16).ravel()/4
        original=np.outer(vector,vector)
        ps=[k.T@k for k in joining_kraus((2,2)).values()]
        output=sum(np.kron(p,np.eye(16))@original@np.kron(p,np.eye(16)) for p in ps)
        self.assertAlmostEqual(trace_distance(output,original),.5)
        np.testing.assert_allclose(np.trace(output.reshape(16,16,16,16),axis1=1,axis2=3),np.eye(16)/16,atol=2e-16)

    def test_exact_record_coding_and_input_guards(self):
        for a,b in combinations(PREFIX_CODES,2):
            self.assertFalse(a.startswith(b) or b.startswith(a))
        self.assertEqual(sum(p*len(code) for p,code in zip(WEIGHTS,PREFIX_CODES)),F(7,3))
        self.assertEqual(sum(WEIGHTS),F(1))
        self.assertTrue(2**2 < len(NAMES) <= 2**3)
        with self.assertRaises(ValueError): physical_record_class((1,1,1,0))
        with self.assertRaises(ValueError): record_unitary(7)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ClassicalJoiningHistoryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":138,"records":list(NAMES),"probabilities_exact":[str(p) for p in WEIGHTS],
        "record_probabilities_depend_on_unknown_input":False,
        "ideal_classical_history_recovers_promised_independent_input":True,
        "recovery_on_unrestricted_inputs":"P_plus omega P_plus + P_minus omega P_minus",
        "recovery_of_every_record_conditioned_input_as_unconditional_source_claimed":False,
        "arbitrary_external_purification_preserved":False,
        "purification_counterexample_trace_error_exact":"1/2",
        "minimum_classical_messages_for_lossless_postprocessing_this_history":7,
        "minimum_fixed_length_bits_in_that_class":3,
        "explicit_prefix_code_expected_bits_exact":"7/3",
        "native_recovery_fresh_reference_initializations":1,
        "native_recovery_yx_upper_bound":8,"native_recovery_ry_upper_bound":12,
        "three_party_optimal_error_after_ideal_classical_recovery_exact":"1/4",
        "finite_noisy_history_readout_completed_in_this_round":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("classical_joining_history_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
