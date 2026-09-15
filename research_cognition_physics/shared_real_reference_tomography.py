"""Round 60: preshared real references enable separated tomography of real states.

This consumes a declared joint preparation resource. It does not change the
intrinsic local state space or establish unassisted local tomography.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator, real_words
from bipartite_composition import interaction
from complex_bipartite_closure import local_effect
from local_record_axiom import tensor_power
from noisy_imaginarity_distillation import resource_state
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z, two_rebit_states
from role_symmetry_and_swap import pauli_word


def shared_reference(systems):
    if systems < 1:
        raise ValueError("At least one reference system is required.")
    return (tensor_power(resource_state(1), systems)+tensor_power(resource_state(-1), systems))/2


def prepare_shared_reference(systems):
    """One maximally mixed seed, real |0> resets, and n-1 old pair flows."""
    if systems < 1:
        raise ValueError("At least one reference system is required.")
    rho = np.kron(IDENTITY/2, tensor_power(np.diag([1., 0]), systems-1))
    for site in range(1, systems):
        gate = embed_operator(interaction(-math.pi/2, "YX"), (0, site), systems)
        rho = gate @ rho @ gate.conj().T
    return rho


def assisted_effect(word, outcomes, beta=ALPHA):
    """A real product effect across enlarged parties (A_i,R_i)."""
    n = len(word)
    effect = np.eye(4**n, dtype=complex)
    active = [i for i, axis in enumerate(word) if axis != "I"]
    if len(outcomes) != len(active):
        raise ValueError("One outcome for each nontrivial axis is required.")
    axes = {"X": PAULI_X, "Y": PAULI_Y, "Z": PAULI_Z}
    for site, sign in zip(active, outcomes):
        if word[site] == "Y":
            observable = embed_operator(np.kron(PAULI_Y, PAULI_Y), (site, n+site), 2*n)
        else:
            observable = embed_operator(axes[word[site]], (site,), 2*n)
        effect = effect @ ((np.eye(4**n)+sign*beta*observable)/2)
    return effect


def assisted_records(density, word, reference=None, beta=ALPHA):
    n = len(word)
    if reference is None:
        reference = shared_reference(n)
    joint = np.kron(density, reference)
    weight = sum(axis != "I" for axis in word)
    return {outcomes: float(np.trace(joint @ assisted_effect(word, outcomes, beta)).real)
            for outcomes in product((-1, 1), repeat=weight)}


def reconstruct_assisted(density):
    n = round(math.log2(len(density)))
    reconstruction = np.zeros_like(density, dtype=complex)
    for word in real_words(n):
        records = assisted_records(density, word)
        weight = sum(axis != "I" for axis in word)
        moment = sum(math.prod(outcomes)*probability for outcomes, probability in records.items())/ALPHA**weight
        reconstruction += moment*pauli_word(word)/(2**n)
    return reconstruction


class SharedRealReferenceTomographyTests(unittest.TestCase):
    def test_old_star_circuit_prepares_the_entire_reference_family(self):
        for n in range(1, 5):
            reference = shared_reference(n)
            np.testing.assert_allclose(prepare_shared_reference(n), reference, atol=3e-16)
            np.testing.assert_array_equal(reference.imag, np.zeros_like(reference.imag))
            self.assertGreaterEqual(np.linalg.eigvalsh(reference).min(), -2e-16)
            self.assertAlmostEqual(np.trace(reference).real, 1, places=14)

    def test_reference_even_y_moments_are_one_and_odd_moments_zero(self):
        for n in range(1, 5):
            reference = shared_reference(n)
            for letters in product("IY", repeat=n):
                word = "".join(letters)
                self.assertEqual(np.trace(reference @ pauli_word(word)).real, int(word.count("Y") % 2 == 0))

    def test_no_single_reference_contains_a_local_imaginary_direction(self):
        for n in range(1, 5):
            reference = shared_reference(n)
            for site in range(n):
                for axis in (PAULI_X, PAULI_Y, PAULI_Z):
                    self.assertEqual(np.trace(reference @ embed_operator(axis, (site,), n)).real, 0)
        local = (np.eye(4)+ALPHA*np.kron(PAULI_Y, PAULI_Y))/2
        effective = np.einsum("abcd,db->ac", local.reshape(2, 2, 2, 2), IDENTITY/2)
        np.testing.assert_array_equal(effective, IDENTITY/2)

    def test_actual_old_pair_gate_and_pointer_readout_implement_the_lifted_y_effect(self):
        gate = interaction(math.pi/2, "YX")
        for sign in (-1, 1):
            actual = gate.conj().T @ np.kron(IDENTITY, local_effect("Z", sign)) @ gate
            expected = (np.eye(4)+sign*ALPHA*np.kron(PAULI_Y, PAULI_Y))/2
            np.testing.assert_allclose(actual, expected, atol=3e-16)

    def test_four_system_circuit_activates_the_old_hidden_pair(self):
        plus, minus = two_rebit_states()
        for q, rho in ((1, plus), (-1, minus)):
            initial = np.kron(rho, shared_reference(2))
            gate_a = embed_operator(interaction(math.pi/2), (0, 2), 4)
            gate_b = embed_operator(interaction(math.pi/2), (1, 3), 4)
            output = gate_b @ gate_a @ initial @ gate_a.conj().T @ gate_b.conj().T
            actual = assisted_records(rho, "YY")
            for s, t in product((-1, 1), repeat=2):
                effect = embed_operator(np.kron(local_effect("Z", s), local_effect("Z", t)), (2, 3), 4)
                p = np.trace(output @ effect).real
                self.assertAlmostEqual(p, (1+s*t*ALPHA**2*q)/4, places=14)
                self.assertAlmostEqual(p, actual[s, t], places=14)

    def test_weak_shared_yy_correlation_has_exactly_linear_activation_strength(self):
        plus, minus = two_rebit_states()
        for c in (0., .2, .8, 1.):
            reference = (np.eye(4)+c*pauli_word("YY"))/4
            left, right = assisted_records(plus, "YY", reference), assisted_records(minus, "YY", reference)
            self.assertAlmostEqual(sum(abs(left[key]-right[key]) for key in left)/2, c*ALPHA**2, places=14)

    def test_all_ten_and_thirty_six_real_coordinates_reconstruct_from_actual_records(self):
        rng = np.random.default_rng(60)
        for n in (2, 3):
            matrix = rng.normal(size=(2**n, 2**n))
            rho = matrix @ matrix.T
            rho /= np.trace(rho)
            np.testing.assert_allclose(reconstruct_assisted(rho), rho, atol=6e-16)

    def test_enlarged_local_effects_are_real_positive_and_form_complete_measurements(self):
        for word in ("YY", "YX", "XZ"):
            effects = [assisted_effect(word, outcomes) for outcomes in product((-1, 1), repeat=2)]
            np.testing.assert_allclose(sum(effects), np.eye(16), atol=2e-16)
            for effect in effects:
                np.testing.assert_array_equal(effect.imag, np.zeros((16, 16)))
                self.assertGreaterEqual(np.linalg.eigvalsh(effect).min(), -1e-16)

    def test_shared_reference_does_not_create_access_to_odd_y_complex_state_differences(self):
        psi = np.array([1., 0, 0, 1j])/math.sqrt(2)
        rho = np.outer(psi, psi.conj())
        for word in ("XY", "YZ", "YY", "XX"):
            a, b = assisted_records(rho, word), assisted_records(rho.conj(), word)
            np.testing.assert_allclose(list(a.values()), list(b.values()), atol=2e-16)

    def test_remote_setting_changes_do_not_change_unconditioned_local_records(self):
        rho = (np.eye(4)+.3*pauli_word("YY")+.2*pauli_word("XI"))/4
        for a in "XYZ":
            marginals = [sum(p for (s, t), p in assisted_records(rho, a+b).items() if s == 1) for b in "XYZ"]
            np.testing.assert_allclose(marginals, marginals[0], atol=3e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SharedRealReferenceTomographyTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"round": 60, "reference": "tau_n=(P_plus_Y^tensor_n+P_minus_Y^tensor_n)/2, a real density matrix",
              "preparation": "I/2 seed, n-1 independent real |0> resets, n-1 original YX(-pi/2) pair gates before distribution",
              "isolated_reference_marginals": "I/2", "new_local_y_readout_required": False,
              "assisted_pair_record_tv": ALPHA**2, "unassisted_local_tomography_restored": False,
              "assisted_real_state_dimensions": {str(n): (4**n+2**n)//2 for n in range(1, 5)},
              "odd_y_complex_coordinates_accessible": False,
              "reference_independent_of_unknown_target_required": True,
              "reference_product_across_parties": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("shared_real_reference_tomography_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
