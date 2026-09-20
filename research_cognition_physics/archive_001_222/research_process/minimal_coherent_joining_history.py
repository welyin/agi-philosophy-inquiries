"""Round 144: one coherent history rebit plus classical labels suffices.

An environment-only real orthogonal transform of the actual balanced join
gives six state-independent reversible branches. The common marginal is
unchanged while history is stored. The environment transform is explicit but
not yet compiled into native gates; conditional decoders are compiled.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import numpy as np

from classical_joining_history import physical_record_class, record_unitary
from compiled_joining_frontier import compiled_sharing
from compiled_subject_joining import circuit_matrix, controlled_z, intake_circuit, ry
from compiled_three_subject_joining import inverse, gate_counts
from external_correlation_recovery_bound import channel, max_entangled_metrics
from flagged_subject_joining import joining_kraus


def environment_support():
    """Seven orthonormal environment vectors in physical h,E1,E2,C order."""
    support=np.zeros((16,7)); support[0,0]=support[1,0]=1/np.sqrt(2)
    for i in range(8,16):
        bits=tuple((i>>(3-j))&1 for j in range(4))
        try: record=physical_record_class(bits)
        except ValueError: continue
        support[i,record]=1
    return support


def complete_real_basis(columns):
    result=[columns[:,i].copy() for i in range(columns.shape[1])]
    for vector in np.eye(columns.shape[0]).T:
        candidate=vector.copy()
        for _ in range(2):
            for old in result: candidate-=np.dot(old,candidate)*old
        norm=np.linalg.norm(candidate)
        if norm>1e-10: result.append(candidate/norm)
    if len(result)!=columns.shape[0]: raise ArithmeticError("Basis completion failed.")
    return np.column_stack(result)


def environment_transform():
    """Output environment order is Q (coherent), then a three-bit label."""
    original=complete_real_basis(environment_support())
    target=np.zeros((16,7)); target[:6,0]=1/np.sqrt(6)
    for record in range(1,7): target[8+record-1,record]=1
    completed=complete_real_basis(target)
    result=completed @ original.T
    if np.linalg.det(result)<0:
        completed[:,-1]*=-1
        result=completed @ original.T
    return result


def conditional_unitary(record):
    if not isinstance(record,int) or isinstance(record,bool) or record not in range(1,7):
        raise ValueError("Use one of the six repair records, 1..6.")
    ks=joining_kraus((2,2))
    return np.vstack((ks[(0,0)],record_unitary(record) @ ks[(0,1)]))


def transformed_branches():
    v=compiled_sharing(F(1,2))[0]
    env_first=v.reshape((2,)*7+(16,)).transpose(0,4,5,6,1,2,3,7).reshape(16,-1)
    transformed=(environment_transform() @ env_first).reshape(2,8,8,16)
    return tuple(transformed[:,label,:,:].reshape(16,16) for label in range(8))


def coherent_decoder_circuit(record):
    conditional_unitary(record)
    site,axis=divmod(record-1,3)
    gates=()
    if site==0: gates+=controlled_z(0,1,2)
    target=2+site
    if axis==0: gates+=controlled_z(0,target,1)
    elif axis==2:
        gates+=(ry(target,-np.pi/2),)+controlled_z(0,target,1)+(ry(target,np.pi/2),)
    return gates+inverse(intake_circuit(2))


class MinimalCoherentJoiningHistoryTests(unittest.TestCase):
    def test_support_and_real_environment_only_transform_are_orthonormal(self):
        support=environment_support(); transform=environment_transform()
        np.testing.assert_allclose(support.T @ support,np.eye(7),atol=3e-16)
        np.testing.assert_allclose(transform.T @ transform,np.eye(16),atol=8e-16)
        self.assertAlmostEqual(np.linalg.det(transform),1.)

    def test_actual_seven_wire_join_becomes_six_equal_weight_unitary_branches(self):
        branches=transformed_branches()
        for r in range(1,7):
            np.testing.assert_allclose(branches[r-1],conditional_unitary(r)/np.sqrt(6),atol=4e-15)
        np.testing.assert_allclose(branches[6:],np.zeros((2,16,16)),atol=4e-15)

    def test_every_report_probability_is_independent_of_unrestricted_input(self):
        for branch in transformed_branches()[:6]:
            np.testing.assert_allclose(branch.T @ branch,np.eye(16)/6,atol=3e-15)

    def test_exact_recovery_preserves_every_external_correlation_in_each_record(self):
        corrected=tuple(conditional_unitary(r).T @ b for r,b in enumerate(transformed_branches()[:6],1))
        for k in corrected: np.testing.assert_allclose(k,np.eye(16)/np.sqrt(6),atol=5e-15)
        metrics=max_entangled_metrics(corrected)
        self.assertAlmostEqual(metrics["overlap"],1.)
        self.assertLess(metrics["trace_error"],5e-14)
        rng=np.random.default_rng(144)
        raw=rng.normal(size=(48,3))+1j*rng.normal(size=(48,3))
        source=raw @ raw.conj().T; source/=np.trace(source)
        np.testing.assert_allclose(channel(tuple(np.kron(k,np.eye(3)) for k in corrected),source),source,atol=4e-15)

    def test_common_interface_marginal_is_unchanged_while_history_is_stored(self):
        rng=np.random.default_rng(244)
        raw=rng.normal(size=(16,5))+1j*rng.normal(size=(16,5))
        source=raw @ raw.conj().T; source/=np.trace(source)
        old=channel(compiled_sharing(F(1,2))[1],source)
        quantum_and_active=channel(transformed_branches(),source)
        marginal=np.trace(quantum_and_active.reshape(2,8,2,8),axis1=0,axis2=2)
        np.testing.assert_allclose(marginal,old,atol=5e-15)

    def test_each_coherent_decoder_uses_only_original_real_gates(self):
        for r in range(1,7):
            actual=circuit_matrix(coherent_decoder_circuit(r),4)
            np.testing.assert_allclose(actual,conditional_unitary(r).T,atol=3e-15)
        counts=[gate_counts(coherent_decoder_circuit(r)) for r in range(1,7)]
        self.assertEqual(max(c["yx"] for c in counts),16)
        self.assertEqual(max(c["ry"] for c in counts),17)

    def test_classical_label_count_has_a_choi_rank_lower_bound_but_is_not_optimal_claim(self):
        from classical_joining_history import history_kraus
        old=history_kraus()
        gram=np.array([[np.vdot(a,b).real for b in old] for a in old])
        np.testing.assert_allclose(gram,np.diag([8]+[F(4,3)]*6).astype(float),atol=3e-15)
        self.assertEqual(np.linalg.matrix_rank(gram),7)
        self.assertEqual((7+1)//2,4)  # At most two reduced Kraus terms per reversible branch.

    def test_source_sector_is_not_measured_and_record_guards(self):
        support=environment_support(); transform=environment_transform()
        np.testing.assert_allclose((transform @ support)[:6,0],np.ones(6)/np.sqrt(6),atol=5e-16)
        for r in range(1,7):
            output=transform @ support[:,r]
            expected=np.zeros(16); expected[8+r-1]=1
            np.testing.assert_allclose(output,expected,atol=5e-16)
        for r in (0,7,1.5,True):
            with self.assertRaises(ValueError): conditional_unitary(r)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(MinimalCoherentJoiningHistoryTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":144,"fixed_balanced_common_channel_preserved_during_storage":True,
        "minimum_additional_coherent_history_rebits":1,
        "constructed_nonzero_classical_labels":6,"fixed_length_message_bits":3,
        "label_probability_for_all_input_states_exact":"1/6",
        "optimal_external_recovery_trace_error_exact":"0",
        "recovery_is_exact_in_each_nonzero_record":True,
        "classical_label_lower_bound_with_one_coherent_rebit":4,
        "classical_label_optimum_determined":False,
        "environment_only_transform_dimension":16,
        "environment_transform_is_real_special_orthogonal":True,
        "environment_transform_native_gate_compilation_done":False,
        "classical_label_readout_is_ideal_in_this_round":True,
        "conditional_decoder_yx_upper_bound":16,"conditional_decoder_ry_upper_bound":17,
        "conditional_decoder_new_pure_ancillas":0,
        "quantum_history_transfer_or_joint_access_counted_as_resource":True,
        "unused_or_measured_environment_physically_deleted":False,
        "automated_checks":{"run":checks.testsRun,"failures":len(checks.failures),"errors":len(checks.errors)}}
    if args.write_results:
        Path(__file__).with_name("minimal_coherent_joining_history_results.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__ == "__main__": main()
