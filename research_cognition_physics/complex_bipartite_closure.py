"""Conditional completion of the two-system complex state and control spaces.

The Hermitian tensor framework and the added Y-sensitive readout remain
assumptions. Local complex control comes from rounds 54-56, not a new primitive.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from bilocal_record_tomography import embed_operator
from bipartite_composition import FULL_BASIS, interaction, visible_summary
from complex_control_from_reference import real_lift, trace_reference
from imaginary_readout_closure import y_effect
from noisy_imaginarity_distillation import resource_state
from quantum_interface_audit import ALPHA, IDENTITY, PAULI_X, PAULI_Y, PAULI_Z, effective_readout
from role_symmetry_and_swap import pauli_word


AXES = ("X", "Y", "Z")
WORDS = tuple(a+b for a, b in product("IXYZ", repeat=2))
NONIDENTITY_WORDS = tuple(word for word in WORDS if word != "II")


def local_effect(axis, sign, new_visibility=.2):
    if axis == "Y":
        return y_effect(sign, new_visibility)
    return effective_readout(math.pi/2 if axis == "X" else 0, (sign+1)//2, 1.)


def local_visibility(axis, new_visibility):
    return new_visibility if axis == "Y" else ALPHA


def product_record_table(density, new_visibility=.2):
    if not 0 < new_visibility <= 1:
        raise ValueError("Nonzero new Y visibility is needed for complete tomography.")
    return {(a, b, s, t): float(np.trace(density @ np.kron(local_effect(a, s, new_visibility), local_effect(b, t, new_visibility))).real)
            for a, b, s, t in product(AXES, AXES, (-1, 1), (-1, 1))}


def reconstruct_product_records(records, new_visibility=.2):
    moments = {"II": 1.}
    for a, b in product(AXES, repeat=2):
        moments[a+b] = sum(s*t*records[a, b, s, t] for s, t in product((-1, 1), repeat=2))/(local_visibility(a, new_visibility)*local_visibility(b, new_visibility))
    for a in AXES:
        moments[a+"I"] = sum(s*records[a, "Z", s, t] for s, t in product((-1, 1), repeat=2))/local_visibility(a, new_visibility)
        moments["I"+a] = sum(t*records["Z", a, s, t] for s, t in product((-1, 1), repeat=2))/local_visibility(a, new_visibility)
    return sum(value*pauli_word(word) for word, value in moments.items())/4


def schmidt_preparation(vector):
    vector = np.asarray(vector, dtype=complex)
    if vector.shape != (4,) or not np.isclose(np.linalg.norm(vector), 1.):
        raise ValueError("A normalized complex four-vector is required.")
    left, values, right_h = np.linalg.svd(vector.reshape(2, 2))
    angle = 2*math.atan2(values[1], values[0])
    right = right_h.T
    entangler = interaction(angle, "YX")
    ket = np.kron(left, right) @ entangler @ np.array([1., 0, 0, 0])
    return {"left": left, "right": right, "entangler_angle": angle, "prepared_vector": ket}


def preparation_with_shared_reference(vector, bias):
    recipe = schmidt_preparation(vector)
    initial = np.diag([1., 0, 0, 0])
    joint = np.kron(resource_state(bias), initial)
    entangler = np.kron(IDENTITY, interaction(recipe["entangler_angle"], "YX"))
    joint = entangler @ joint @ entangler.conj().T
    for gate, sites in ((recipe["left"], (0, 1)), (recipe["right"], (0, 2))):
        lifted = embed_operator(real_lift(gate), sites, 3)
        joint = lifted @ joint @ lifted.conj().T
    return joint


def su2_axis_rotations():
    # Conjugating Y by these gives the named axis, and similarly for X.
    rz = lambda angle: math.cos(angle/2)*IDENTITY - 1j*math.sin(angle/2)*PAULI_Z
    rx = lambda angle: math.cos(angle/2)*IDENTITY - 1j*math.sin(angle/2)*PAULI_X
    ry = lambda angle: math.cos(angle/2)*IDENTITY - 1j*math.sin(angle/2)*PAULI_Y
    from_y = {"X": rz(-math.pi/2), "Y": IDENTITY, "Z": rx(math.pi/2)}
    from_x = {"X": IDENTITY, "Y": rz(math.pi/2), "Z": ry(-math.pi/2)}
    return from_y, from_x


class ComplexBipartiteClosureTests(unittest.TestCase):
    def test_six_new_state_directions_complete_the_exact_sixteen_dimensional_basis(self):
        matrices = np.array([pauli_word(word) for word in WORDS])
        flat = matrices.reshape(16, 16)
        np.testing.assert_array_equal(flat.conj() @ flat.T, 4*np.eye(16))
        self.assertEqual([word for word in WORDS if word.count("Y") % 2], ["IY", "XY", "YI", "YX", "YZ", "ZY"])
        self.assertEqual(len(FULL_BASIS), 10)

    def test_old_interaction_and_local_rotations_provide_all_nine_couplings(self):
        first, second = su2_axis_rotations()
        for a, b in product(AXES, repeat=2):
            local = np.kron(first[a], second[b])
            np.testing.assert_allclose(local @ pauli_word("YX") @ local.conj().T, pauli_word(a+b), atol=4e-16)

    def test_fifteen_control_directions_have_exact_orthogonal_gram_certificate(self):
        matrices = np.array([pauli_word(word) for word in NONIDENTITY_WORDS])
        for matrix in matrices:
            self.assertEqual(np.trace(matrix), 0)
            np.testing.assert_array_equal(matrix, matrix.conj().T)
        flat = matrices.reshape(15, 16)
        np.testing.assert_array_equal(flat.conj() @ flat.T/4, np.eye(15))

    def test_schmidt_recipe_prepares_every_sample_pure_state_with_one_old_entangling_flow(self):
        rng = np.random.default_rng(57)
        vectors = [np.array([1, 0, 0, 1j])/math.sqrt(2), np.array([0, 1, 0, 0])]
        for _ in range(25):
            vector = rng.normal(size=4) + 1j*rng.normal(size=4)
            vectors.append(vector/np.linalg.norm(vector))
        for vector in vectors:
            recipe = schmidt_preparation(vector)
            np.testing.assert_allclose(recipe["prepared_vector"], vector, atol=8e-16)
            self.assertLessEqual(recipe["entangler_angle"], math.pi/2 + 1e-15)

    def test_shared_reference_physically_compiles_both_local_complex_rotations(self):
        rng = np.random.default_rng(157)
        for _ in range(8):
            vector = rng.normal(size=4) + 1j*rng.normal(size=4)
            vector /= np.linalg.norm(vector)
            target = np.outer(vector, vector.conj())
            for bias in (.2, .9, 1.):
                joint = preparation_with_shared_reference(vector, bias)
                expected = ((1+bias)/2*np.kron(resource_state(1.), target)
                            +(1-bias)/2*np.kron(resource_state(-1.), target.conj()))
                np.testing.assert_allclose(joint, expected, atol=6e-16)
                self.assertLessEqual(np.linalg.norm(trace_reference(joint)-target, ord="nuc")/2, (1-bias)/2 + 8e-16)

    def test_classical_mixing_of_pure_recipes_recovers_arbitrary_complex_positive_states(self):
        rng = np.random.default_rng(257)
        for _ in range(10):
            matrix = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
            density = matrix @ matrix.conj().T
            density /= np.trace(density)
            values, vectors = np.linalg.eigh(density)
            output = sum(weight*np.outer(schmidt_preparation(vector)["prepared_vector"], schmidt_preparation(vector)["prepared_vector"].conj())
                         for weight, vector in zip(values, vectors.T))
            np.testing.assert_allclose(output, density, atol=7e-16)

    def test_original_noisy_axes_and_weak_new_y_axis_give_complete_product_tomography(self):
        rng = np.random.default_rng(357)
        for visibility in (.05, .2, 1.):
            matrix = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
            density = matrix @ matrix.conj().T
            density /= np.trace(density)
            records = product_record_table(density, visibility)
            reconstructed = reconstruct_product_records(records, visibility)
            np.testing.assert_allclose(reconstructed, density, atol=2e-14)
            for a, b in product(AXES, repeat=2):
                self.assertAlmostEqual(sum(records[a, b, s, t] for s, t in product((-1, 1), repeat=2)), 1., places=14)

    def test_new_imaginary_bell_pair_is_invisible_to_old_real_effects_but_locally_distinct(self):
        vector = np.array([1, 0, 0, 1j])/math.sqrt(2)
        plus = np.outer(vector, vector.conj())
        minus = plus.conj()
        for matrix in FULL_BASIS:
            self.assertEqual(np.trace((plus-minus) @ matrix).real, 0.)
        records = [product_record_table(density, .2) for density in (plus, minus)]
        distance = sum(abs(records[0]["X", "Y", s, t]-records[1]["X", "Y", s, t]) for s, t in product((-1, 1), repeat=2))/2
        self.assertAlmostEqual(distance, .2*ALPHA, places=14)

    def test_product_measurements_are_nonsignaling_after_the_extension(self):
        vector = np.array([.2, .3j, .4, math.sqrt(.71)])
        density = np.outer(vector, vector.conj())
        records = product_record_table(density)
        for a, s in product(AXES, (-1, 1)):
            marginals = [sum(records[a, b, s, t] for t in (-1, 1)) for b in AXES]
            np.testing.assert_allclose(marginals, marginals[0], atol=2e-16)

    def test_zero_new_visibility_is_rejected_as_a_complete_local_calibration(self):
        with self.assertRaises(ValueError):
            product_record_table(np.eye(4)/4, 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ComplexBipartiteClosureTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"scope": "Conditional Hermitian tensor extension with the weak Y effect, independent resources, original real interactions, and the control/state closures proved in rounds 54-56",
              "complex_two_system_control_closure": "SU(4), up to physically irrelevant global phase",
              "control_generator_count": 15, "exact_normalized_gram_determinant": 1,
              "old_real_homogeneous_state_dimension": 10, "complex_homogeneous_state_dimension": 16,
              "complex_normalized_state_dimension": 15,
              "pure_preparation_recipe": "One original YX entangling flow plus two local complex rotations, compiled with one shared Y reference",
              "new_odd_y_state_words": [word for word in WORDS if word.count("Y") % 2],
              "finite_setting_local_tomography": "Nine pairs of X/Y/Z settings, original noisy X/Z and any nonzero new Y visibility",
              "new_readout_is_ideal_required_for_local_tomography": False,
              "closed_bipartite_state_set": "All complex Hermitian 4x4 density matrices within the assumed ambient framework",
              "time_hamiltonian_or_spacetime_derived": False,
              "new_y_readout_derived_from_cognition": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("complex_bipartite_closure_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
