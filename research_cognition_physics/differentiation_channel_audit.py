"""Round 78: audit a user-proposed whole-to-specialized-subsystem analogy.

Taking a state's real part is positive but fails on an entangled spectator.
Erasing Y while equally retaining X,Z by v is a physical channel iff v<=1/2.
Partitioning a shared whole, discarding a reference, and erasing Y are distinct.
"""

import argparse
import json
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bell_network_statistics import bell_projectors
from independent_source_alignment import keep_systems,pair_reference
from quantum_interface_audit import IDENTITY,PAULI_X,PAULI_Y,PAULI_Z,ALPHA
from role_symmetry_and_swap import pauli_word


def erase_y_map(matrix,retention=1.):
    """Complex-linear extension; on Hermitian inputs transpose equals conjugate."""
    if not 0<=retention<=1: raise ValueError("Use X/Z retention in [0,1].")
    return retention*(matrix+matrix.T)/2+(1-retention)*np.trace(matrix)*IDENTITY/2


def apply_to_second(density,retention=1.):
    tensor=density.reshape(2,2,2,2)
    output=np.zeros_like(tensor,dtype=complex)
    for i,j in product(range(2),repeat=2): output[i,:,j,:]=erase_y_map(tensor[i,:,j,:],retention)
    return output.reshape(4,4)


def bell_extension(retention=1.):
    return apply_to_second(bell_projectors()[0],retention)


def physical_real_kraus(retention):
    if not 0<=retention<=.5: raise ValueError("A physical complete Y erasure requires retention <= 1/2.")
    matrices=(IDENTITY,PAULI_X,-1j*PAULI_Y,PAULI_Z)
    weights=((1+2*retention)/4,1/4,(1-2*retention)/4,1/4)
    return tuple(np.sqrt(p)*matrix for p,matrix in zip(weights,matrices))


def distributed_whole(target_correlation,reference_correlation=1.):
    # Coordinate regrouping into A target,reference; C target,reference.
    return keep_systems(np.kron(pair_reference(target_correlation),pair_reference(reference_correlation)),(0,2,1,3),4)


def local_reference_records(target_correlation,reference_correlation=1.,visibility=ALPHA):
    whole=distributed_whole(target_correlation,reference_correlation)
    lifted_y=pauli_word("YY")
    return {(a,c):float(np.trace(whole@np.kron((np.eye(4)+a*visibility*lifted_y)/2,
                                             (np.eye(4)+c*visibility*lifted_y)/2)).real)
            for a,c in product((-1,1),repeat=2)}


class DifferentiationChannelAuditTests(unittest.TestCase):
    def test_real_part_rule_preserves_valid_isolated_states_and_xz_coordinates(self):
        rng=np.random.default_rng(78)
        for _ in range(15):
            raw=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
            state=raw@raw.conj().T
            state/=np.trace(state)
            output=erase_y_map(state)
            np.testing.assert_allclose(output,state.real,atol=0)
            self.assertGreaterEqual(np.linalg.eigvalsh(output).min(),-1e-16)
            for axis in (PAULI_X,PAULI_Z):
                self.assertAlmostEqual(np.trace(output@axis),np.trace(state@axis),places=15)
            self.assertAlmostEqual(np.trace(output@PAULI_Y),0.,places=15)

    def test_the_same_rule_on_half_a_real_bell_state_produces_a_negative_probability(self):
        antisymmetric=np.array([0,1,-1,0])/np.sqrt(2)
        output=bell_extension()
        self.assertAlmostEqual(antisymmetric@output@antisymmetric,-.25,places=15)
        self.assertAlmostEqual(np.trace(output).real,1.,places=15)

    def test_exact_bell_basis_eigenvalues_give_the_complete_positivity_threshold(self):
        for retention in (0.,.25,.5,.75,1.):
            actual=np.linalg.eigvalsh(bell_extension(retention))
            expected=sorted(((1+2*retention)/4,.25,.25,(1-2*retention)/4))
            np.testing.assert_allclose(actual,expected,atol=3e-16)

    def test_explicit_real_kraus_construction_proves_sufficiency_at_and_below_half(self):
        for retention in (0.,.2,.5):
            kraus=physical_real_kraus(retention)
            np.testing.assert_allclose(sum(k.conj().T@k for k in kraus),IDENTITY,atol=3e-16)
            for k in kraus: np.testing.assert_allclose(k.imag,0.,atol=0)
            # All matrix units, including non-Hermitian ones, verify the map.
            for i,j in product(range(2),repeat=2):
                unit=np.zeros((2,2),dtype=complex); unit[i,j]=1
                actual=sum(k@unit@k.conj().T for k in kraus)
                np.testing.assert_allclose(actual,erase_y_map(unit,retention),atol=2e-16)

    def test_valid_differentiation_does_not_become_reversible_on_new_y_states(self):
        plus,minus=(IDENTITY+PAULI_Y)/2,(IDENTITY-PAULI_Y)/2
        np.testing.assert_array_equal(erase_y_map(plus,.5),erase_y_map(minus,.5))
        with self.assertRaises(ValueError): physical_real_kraus(.6)

    def test_shared_reference_preserves_joint_distinction_under_local_partition(self):
        plus=local_reference_records(1.)
        minus=local_reference_records(-1.)
        distance=sum(abs(plus[key]-minus[key]) for key in plus)/2
        self.assertAlmostEqual(distance,ALPHA**2,places=15)
        for a,c in product((-1,1),repeat=2): self.assertAlmostEqual(plus[a,c],(1+a*c*ALPHA**2)/4,places=15)

    def test_independent_references_remove_the_local_reading_capability(self):
        plus=local_reference_records(1.,0.)
        minus=local_reference_records(-1.,0.)
        self.assertEqual(plus,minus)
        for p in plus.values(): self.assertEqual(p,.25)

    def test_discarding_references_does_not_erase_the_targets_global_yy_coordinate(self):
        for correlation in (-1.,1.):
            whole=distributed_whole(correlation)
            target=keep_systems(whole,(0,2),4)
            np.testing.assert_allclose(target,pair_reference(correlation),atol=0)
            self.assertAlmostEqual(np.trace(target@pauli_word("YY")).real,correlation,places=15)
            np.testing.assert_allclose(keep_systems(whole,(0,1),4),np.eye(4)/4,atol=0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DifferentiationChannelAuditTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={"round":78,"motivation":"User analogy: an undifferentiated whole, specialized real subsystems, and communication as collaboration",
            "naive_real_part_map_is_a_universal_physical_channel":False,
            "bell_test_negative_probability_exact":"-1/4",
            "equal_xz_retention_with_complete_y_erasure_cp_range":"0 <= v <= 1/2",
            "cp_sufficiency":"Explicit real Kraus matrices sqrt(p) times I,X,-iY,Z",
            "partition_alone_forces_loss_of_joint_reading_capability":False,
            "discard_reference_erases_global_target_yy":False,
            "shared_reference_local_record_TV_diagnostic":ALPHA**2,
            "proposed_next_object":"Whole preparation + chosen subsystem partition + permitted local instruments/resources + message graph; add a physical channel only when actual state change is intended",
            "complex_structure_is_derived_as_the_unique_undifferentiated_whole":False,
            "quantum_theory_derived_from_cognition":False,
            "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("differentiation_channel_audit_results.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__ == "__main__": main()
