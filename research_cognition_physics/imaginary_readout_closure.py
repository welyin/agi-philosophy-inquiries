"""Conditional states forced by an explicitly added Y-sensitive local effect.

Hermitian matrices, the tensor product and the trace probability rule are
assumptions of this extension audit, not derived from cognitive principles.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from quantum_interface_audit import (IDENTITY, PAULI_X, PAULI_Y, PAULI_Z,
                                     effective_readout, two_rebit_states)
from role_symmetry_and_swap import proper_rotation_to


def y_effect(sign, visibility):
    if sign not in (-1, 1) or not 0 <= visibility <= 1:
        raise ValueError("Use a sign +/-1 and visibility in [0,1].")
    return (IDENTITY + sign*visibility*PAULI_Y)/2


def isolate_y_readout(effect):
    """Real half-turn randomization and classical outcome postprocessing."""
    if np.linalg.eigvalsh(effect).min() < -1e-14 or np.linalg.eigvalsh(IDENTITY-effect).min() < -1e-14:
        raise ValueError("A valid binary effect is required.")
    e0 = np.trace(effect).real/2
    y = np.trace(effect @ PAULI_Y).real/2
    if abs(y) < 1e-15:
        raise ValueError("The effect has no Y-sensitive component.")
    if y < 0:
        effect, e0, y = IDENTITY-effect, 1-e0, -y
    half_turn = -1j*PAULI_Y
    averaged = (effect + half_turn.conj().T @ effect @ half_turn)/2
    scale = 1/(2*max(e0, 1-e0))
    offset = .5-scale*e0
    return {"effect": scale*averaged + offset*IDENTITY,
            "visibility": 2*scale*y, "accept_old_minus": offset, "accept_old_plus": offset+scale}


def conditional_remote(density, effect):
    """First subsystem has dimension 2. Return unnormalized remote state."""
    dimension = density.shape[0] // 2
    weighted = np.kron(effect, np.eye(dimension)) @ density
    return np.einsum("abad->bd", weighted.reshape(2, dimension, 2, dimension))


def real_encoding(vector):
    """|Omega> = |0> Re(psi) - |1> Im(psi), preparing psi on outcome +Y."""
    vector = np.asarray(vector, dtype=complex)
    if vector.ndim != 1 or not np.isclose(np.linalg.norm(vector), 1.):
        raise ValueError("The target vector must be normalized.")
    return np.r_[vector.real, -vector.imag]


def density_from_bloch(x, y, z):
    return (IDENTITY + x*PAULI_X + y*PAULI_Y + z*PAULI_Z)/2


def bloch(density):
    return np.array([np.trace(density @ p).real for p in (PAULI_X, PAULI_Y, PAULI_Z)])


def conjugation_mixture(density, visibility):
    return (1+visibility)*density/2 + (1-visibility)*density.conj()/2


def prepare_ellipsoid_point(point, visibility):
    """Construct a real bipartite mixed preparation for any attainable point."""
    x, y, z = point
    if not 0 <= visibility <= 1:
        raise ValueError("Use visibility in [0,1].")
    if visibility == 0:
        if abs(y) > 1e-14:
            raise ValueError("Zero Y visibility only attains the real disk.")
        target = density_from_bloch(x, 0, z)
    else:
        target = density_from_bloch(x, y/visibility, z)
    weights, vectors = np.linalg.eigh(target)
    if weights.min() < -1e-13:
        raise ValueError("Point lies outside the one-use ellipsoid.")
    return sum(max(0., weight)*np.outer(real_encoding(vector), real_encoding(vector))
               for weight, vector in zip(weights, vectors.T))


def pure_state_error(vector, visibility):
    overlap = abs(np.dot(vector, vector))**2
    return (1-visibility)*math.sqrt(max(0., 1-overlap))/2


class ImaginaryReadoutClosureTests(unittest.TestCase):
    def test_any_y_sensitive_effect_yields_an_unbiased_y_readout(self):
        for e0, y in product((.3, .5, .7), (-.1, .1)):
            effect = e0*IDENTITY + .07*PAULI_X + y*PAULI_Y + .04*PAULI_Z
            report = isolate_y_readout(effect)
            np.testing.assert_allclose(report["effect"], y_effect(1, report["visibility"]), atol=2e-16)
            self.assertGreater(report["visibility"], 0)
            self.assertGreaterEqual(report["accept_old_minus"], -1e-16)
            self.assertLessEqual(report["accept_old_plus"], 1+1e-16)
        with self.assertRaises(ValueError):
            isolate_y_readout(effective_readout(.2, 1, .5))

    def test_y_readout_is_a_fair_coin_on_every_old_real_preparation(self):
        for angle, radius, nu, sign in product((0., .4, 1.7), (0., .3, 1.), (0., .2, 1.), (-1, 1)):
            density = density_from_bloch(radius*math.sin(angle), 0, radius*math.cos(angle))
            self.assertAlmostEqual(np.trace(density @ y_effect(sign, nu)).real, .5, places=14)

    def test_both_parties_need_a_y_component_to_read_the_hidden_pair(self):
        plus, minus = two_rebit_states()
        for a, b, s, t in product((0., .2, 1.), (0., .7, 1.), (-1, 1), (-1, 1)):
            effect = np.kron(y_effect(s, a), y_effect(t, b))
            self.assertAlmostEqual(np.trace((plus-minus) @ effect).real, s*t*a*b/2, places=14)
            self.assertAlmostEqual(np.trace(plus @ effect).real, (1+s*t*a*b)/4, places=14)
        real_effect = effective_readout(.3, 1, .8)
        self.assertEqual(np.trace((plus-minus) @ np.kron(real_effect, y_effect(1, 1.))).real, 0.)

    def test_bell_conditional_state_is_transposed_effect_and_average_is_unchanged(self):
        vector = np.array([1., 0, 0, 1.])/math.sqrt(2)
        density = np.outer(vector, vector)
        for nu in (0., .2, 1.):
            outputs = [conditional_remote(density, y_effect(s, nu)) for s in (-1, 1)]
            for s, output in zip((-1, 1), outputs):
                np.testing.assert_allclose(output, y_effect(s, nu).T/2, atol=2e-16)
                np.testing.assert_allclose(2*output, density_from_bloch(0, -s*nu, 0), atol=3e-16)
            np.testing.assert_allclose(sum(outputs), IDENTITY/2, atol=2e-16)

    def test_every_complex_pure_vector_has_a_real_conditional_preparation(self):
        rng = np.random.default_rng(54)
        for dimension in (2, 3, 4, 8):
            for _ in range(8):
                vector = rng.normal(size=dimension) + 1j*rng.normal(size=dimension)
                vector /= np.linalg.norm(vector)
                encoded = real_encoding(vector)
                self.assertAlmostEqual(np.linalg.norm(encoded), 1., places=14)
                output = conditional_remote(np.outer(encoded, encoded), y_effect(1, 1.))
                np.testing.assert_allclose(output, np.outer(vector, vector.conj())/2, atol=2e-16)

    def test_single_target_real_encoding_can_be_prepared_by_the_existing_so4_closure(self):
        vector = np.array([.3 + .4j, math.sqrt(.75)])
        encoded = real_encoding(vector)
        gate = proper_rotation_to(encoded)
        np.testing.assert_allclose(gate[:, 0], encoded, atol=3e-16)
        self.assertAlmostEqual(np.linalg.det(gate), 1., places=14)

    def test_noisy_effect_gives_the_exact_conjugate_mixture_in_any_dimension(self):
        rng = np.random.default_rng(154)
        for dimension, nu in product((2, 5), (.1, .5, 1.)):
            vector = rng.normal(size=dimension) + 1j*rng.normal(size=dimension)
            vector /= np.linalg.norm(vector)
            encoded = real_encoding(vector)
            output = conditional_remote(np.outer(encoded, encoded), y_effect(1, nu))
            np.testing.assert_allclose(2*output, conjugation_mixture(np.outer(vector, vector.conj()), nu), atol=2e-16)
            self.assertAlmostEqual(np.trace(output).real, .5, places=14)

    def test_arbitrary_real_joint_states_obey_the_one_use_ellipsoid_bound(self):
        rng = np.random.default_rng(254)
        for nu in (.05, .4, 1.):
            for _ in range(25):
                matrix = rng.normal(size=(4, 4))
                density = matrix @ matrix.T
                density /= np.trace(density)
                output = 2*conditional_remote(density, y_effect(1, nu))
                x, y, z = bloch(output)
                self.assertLessEqual(x*x + z*z + (y/nu)**2, 1 + 5e-16)

    def test_all_ellipsoid_points_have_constructive_real_preparations(self):
        rng = np.random.default_rng(354)
        for nu in (0., .2, .8, 1.):
            for radius in (0., .3, 1.):
                direction = rng.normal(size=3)
                direction *= radius/np.linalg.norm(direction)
                direction[1] *= nu
                preparation = prepare_ellipsoid_point(direction, nu)
                self.assertGreaterEqual(np.linalg.eigvalsh(preparation).min(), -2e-16)
                output = 2*conditional_remote(preparation, y_effect(1, nu))
                np.testing.assert_allclose(bloch(output), direction, atol=1e-15)
        with self.assertRaises(ValueError):
            prepare_ellipsoid_point((0, .6, 0), .5)

    def test_pure_target_error_formula_and_worst_case_are_tight(self):
        rng = np.random.default_rng(454)
        vectors = [np.array([1, 1j])/math.sqrt(2)]
        for _ in range(10):
            vector = rng.normal(size=4) + 1j*rng.normal(size=4)
            vectors.append(vector/np.linalg.norm(vector))
        for vector, nu in product(vectors, (.2, .7, 1.)):
            target = np.outer(vector, vector.conj())
            actual = np.linalg.norm(conjugation_mixture(target, nu)-target, ord="nuc")/2
            self.assertAlmostEqual(actual, pure_state_error(vector, nu), places=14)
            self.assertLessEqual(actual, (1-nu)/2 + 2e-16)
        self.assertAlmostEqual(pure_state_error(vectors[0], .2), .4, places=14)

    def test_real_local_output_cannot_match_the_new_conditional_joint_statistics(self):
        conditional = density_from_bloch(0, .3, 0)
        imaginary_read = y_effect(1, .2)
        for angle in (0., .7, 1.5):
            real_state = density_from_bloch(math.sin(angle), 0, math.cos(angle))
            self.assertAlmostEqual(np.trace(real_state @ imaginary_read).real, .5, places=14)
        self.assertAlmostEqual(np.trace(conditional @ imaginary_read).real, .53, places=14)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ImaginaryReadoutClosureTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {"scope": "Added E_s=(I+s*nu*Y)/2 in a Hermitian trace-probability tensor framework; full old real joint preparations retained",
              "both_sides_need_y_sensitive_effects_for_product_discrimination": True,
              "hidden_pair_product_record_TV": "nu_A*nu_B",
              "conditional_preparation_probability": "1/2",
              "exact_one_use_remote_bloch_set": "x^2+z^2+(y/nu)^2 <= 1; nu=0 gives the old real disk",
              "ideal_y_readout_forces_all_complex_single_pure_states": True,
              "arbitrary_target_real_encoding": "|0> Re(psi) - |1> Im(psi)",
              "noisy_conditional_output": "(1+nu)/2 * |psi><psi| + (1-nu)/2 * |psi*><psi*|",
              "worst_case_trace_distance_to_pure_target": "(1-nu)/2",
              "one_noisy_use_forces_the_full_bloch_ball": False,
              "new_y_readout_derived_from_cognition": False,
              "quantum_theory_derived_from_cognition": False,
              "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0}}
    if args.write_results:
        Path(__file__).with_name("imaginary_readout_closure_results.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
