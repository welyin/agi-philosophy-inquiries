"""Round 88: distinguish real-state preservation from real-shadow closure."""

import argparse
import json
import unittest
from pathlib import Path

import numpy as np

from complex_whole_real_interfaces import random_state, real_instrument
from partition_interface_closure import sector_bases
from quantum_interface_audit import IDENTITY, PAULI_X, PAULI_Y, PAULI_Z


def action(kraus,matrix):
    return sum(k@matrix@k.conj().T for k in kraus)


def hermitian_sectors(d):
    symmetric,skew=sector_bases(d)
    return symmetric,tuple(1j*a for a in skew)


def leakage_blocks(kraus):
    d=kraus[0].shape[1]
    symmetric,imaginary=hermitian_sectors(d)
    # Q Phi P measures creation of imaginary outputs from real inputs.
    # P Phi Q measures visible feedback from imaginary inputs.
    return (np.array([action(kraus,s).imag for s in symmetric]),
            np.array([action(kraus,a).real for a in imaginary]))


def choi(kraus):
    return sum(np.outer(k.reshape(-1,order="F"),k.reshape(-1,order="F").conj())
               for k in kraus)


def real_choi_kraus(matrix,d):
    if not np.allclose(matrix.imag,0,atol=1e-12,rtol=0):
        raise ValueError("The Choi matrix must be real.")
    values,vectors=np.linalg.eigh(matrix.real)
    if values.min() < -1e-12:
        raise ValueError("The Choi matrix must be positive.")
    return tuple(np.sqrt(v)*vectors[:,i].reshape(d,d,order="F")
                 for i,v in enumerate(values) if v>1e-13)


def unitary(hamiltonian,t):
    values,vectors=np.linalg.eigh(hamiltonian)
    return (vectors*np.exp(-1j*t*values))@vectors.conj().T


def hamiltonian_constraint_dimension(d):
    symmetric,imaginary=hermitian_sectors(d)
    basis=symmetric+imaginary
    constraints=np.array([
        np.concatenate([(-1j*(h@s-s@h)).imag.reshape(-1) for s in symmetric])
        for h in basis]).T
    return len(basis)-np.linalg.matrix_rank(constraints)


class RealInterfaceStabilityTests(unittest.TestCase):
    def test_real_outputs_can_depend_on_previously_hidden_imaginary_input(self):
        yp=np.array([1,1j])/np.sqrt(2)
        ym=yp.conj()
        kraus=(np.outer([1.,0.],yp.conj()),np.outer([0.,1.],ym.conj()))
        creation,feedback=leakage_blocks(kraus)
        np.testing.assert_allclose(creation,0,atol=1e-16)
        self.assertGreater(np.linalg.norm(feedback),1)
        np.testing.assert_allclose(action(kraus,PAULI_Y),PAULI_Z,atol=3e-16)
        np.testing.assert_allclose(sum(k.conj().T@k for k in kraus),IDENTITY,atol=3e-16)

    def test_shadow_closure_alone_can_generate_imaginary_states(self):
        yp=np.array([1,1j])/np.sqrt(2)
        kraus=tuple(np.outer(yp,e) for e in np.eye(2))
        creation,feedback=leakage_blocks(kraus)
        self.assertGreater(np.linalg.norm(creation),0)
        np.testing.assert_allclose(feedback,0,atol=0)
        np.testing.assert_allclose(action(kraus,IDENTITY/2),(IDENTITY+PAULI_Y)/2,atol=2e-16)

    def test_both_zero_blocks_give_real_choi_and_real_kraus_reconstruction(self):
        rng=np.random.default_rng(88)
        for d in (2,3,4):
            raw=real_instrument(rng,d)
            # Complex Kraus labels may represent the same real operation.
            mixed=((raw[0]+1j*raw[1])/np.sqrt(2),
                   (1j*raw[0]+raw[1])/np.sqrt(2))
            creation,feedback=leakage_blocks(mixed)
            np.testing.assert_allclose(creation,0,atol=2e-16)
            np.testing.assert_allclose(feedback,0,atol=2e-16)
            matrix=choi(mixed)
            np.testing.assert_allclose(matrix.imag,0,atol=3e-16)
            recovered=real_choi_kraus(matrix,d)
            for e in np.eye(d*d).reshape(d*d,d,d):
                np.testing.assert_allclose(action(recovered,e),action(raw,e),atol=3e-15)

    def test_unrecorded_average_can_be_real_while_each_recorded_branch_is_not(self):
        plus=np.diag([1,1j])/np.sqrt(2)
        minus=plus.conj()
        for k in (plus,minus):
            self.assertGreater(np.linalg.norm(leakage_blocks((k,))[0]),0)
            self.assertGreater(np.linalg.norm(leakage_blocks((k,))[1]),0)
        for block in leakage_blocks((plus,minus)):
            np.testing.assert_allclose(block,0,atol=0)
        rho=(IDENTITY+PAULI_X)/2
        np.testing.assert_allclose(action((plus,minus),rho),IDENTITY/2,atol=2e-16)

    def test_unitary_real_preservation_is_unchanged_by_global_phase(self):
        rng=np.random.default_rng(188)
        for d in (2,3,4):
            orthogonal,_=np.linalg.qr(rng.normal(size=(d,d)))
            u=np.exp(.731j)*orthogonal
            for block in leakage_blocks((u,)):
                np.testing.assert_allclose(block,0,atol=3e-16)
            np.testing.assert_allclose(u.conj().T@u.conj(),
                                       np.exp(-2*.731j)*np.eye(d),atol=7e-16)

    def test_complete_generator_kernel_is_scalar_plus_imaginary_skew(self):
        rng=np.random.default_rng(288)
        for d in (2,3,4):
            self.assertEqual(hamiltonian_constraint_dimension(d),1+d*(d-1)//2)
            raw=rng.normal(size=(d,d))
            h=.73*np.eye(d)+1j*(raw-raw.T)
            for t in (.1,.7,1.9):
                u=unitary(h,t)
                for block in leakage_blocks((u,)):
                    np.testing.assert_allclose(block,0,atol=8e-15)

    def test_a_real_hamiltonian_need_not_preserve_real_states_forward_in_time(self):
        u=unitary(PAULI_Z/2,np.pi/2)
        np.testing.assert_allclose(action((u,),(IDENTITY+PAULI_X)/2),
                                   (IDENTITY+PAULI_Y)/2,atol=2e-16)
        self.assertGreater(np.linalg.norm(leakage_blocks((u,))[0]),1)

    def test_real_branch_covariance_survives_a_correlated_spectator(self):
        rng=np.random.default_rng(388)
        state=random_state(rng,6)
        for k in real_instrument(rng,3):
            extended=np.kron(k,np.eye(2))
            actual=action((extended,),state)
            np.testing.assert_allclose(actual.conj(),action((extended,),state.conj()),atol=0)
            np.testing.assert_allclose(actual.real,action((extended,),state.real),atol=2e-17)
        # Isolated shadow closure alone fails after adding one untouched qubit.
        yp=np.array([1,1j])/np.sqrt(2)
        prepare=tuple(np.outer(yp,e) for e in np.eye(2))
        extended=tuple(np.kron(k,IDENTITY) for k in prepare)
        inputs=[np.kron(IDENTITY/2,(IDENTITY+s*PAULI_Y)/2) for s in (-1,1)]
        np.testing.assert_array_equal(inputs[0].real,inputs[1].real)
        outputs=[action(extended,rho).real for rho in inputs]
        self.assertGreater(np.linalg.norm(outputs[1]-outputs[0]),.5)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RealInterfaceStabilityTests))
    if not checks.wasSuccessful(): raise SystemExit(1)
    report={
        "round":88,
        "state_preservation_condition":"Q Phi P = 0",
        "predictive_shadow_closure_condition":"P Phi Q = 0",
        "conditions_are_distinct_for_general_complex_cp_maps":True,
        "both_conditions_equivalent_to_conjugation_covariance_and_real_choi":True,
        "covariant_cp_maps_admit_a_real_kraus_representation":True,
        "shadow_closure_with_one_untouched_qubit_forces_both_blocks_zero":True,
        "each_observed_instrument_branch_must_be_checked":True,
        "unitary_class":"exp(i theta) O, O real orthogonal",
        "continuous_hamiltonian_class":"h I + i A, A real antisymmetric",
        "hamiltonian_parameter_dimension":"1+d(d-1)/2",
        "real_hamiltonian_alone_suffices":False,
        "permission_stability_is_derived_from_cognition":False,
        "quantum_theory_derived_from_cognition":False,
        "automated_checks":{"run":checks.testsRun,"failures":0,"errors":0}}
    if args.write_results:
        Path(__file__).with_name("real_interface_stability_results.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
