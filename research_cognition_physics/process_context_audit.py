"""Round 174: local prediction equality does not imply equality of dilations.

Two pure-environment real dilations agree on every single-rebit state but are
perfectly distinguishable using a real Bell input. This audits an assumption;
it does not refute purification in real quantum theory.
"""

import argparse
import json
from pathlib import Path
import unittest

import numpy as np


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])
J = np.array([[0., -1.], [1., 0.]])
PLUS = (I, J)
MINUS = (X, Z)


def channel(pair, matrix):
    return sum(k @ matrix @ k.T for k in pair) / 2


def isometry(pair):
    # Output index (system, environment); all Kraus weights are 1/2.
    return np.stack(pair, axis=1).reshape(4, 2) / np.sqrt(2)


def joint_output(pair, density):
    return sum(np.kron(k, I) @ density @ np.kron(k.T, I) for k in pair) / 2


def witness_states():
    bell = np.outer([1., 0., 0., 1.], [1., 0., 0., 1.]) / 2
    return joint_output(PLUS, bell), joint_output(MINUS, bell)


class ProcessContextAuditTests(unittest.TestCase):
    def test_channels_agree_on_full_real_symmetric_basis(self):
        for basis in (I, X, Z):
            expected = np.trace(basis) * I / 2
            np.testing.assert_array_equal(channel(PLUS, basis), expected)
            np.testing.assert_array_equal(channel(MINUS, basis), expected)

    def test_antisymmetric_actions_and_kappa_are_opposite(self):
        np.testing.assert_array_equal(channel(PLUS, J), J)
        np.testing.assert_array_equal(channel(MINUS, J), -J)
        self.assertEqual(sum(np.linalg.det(k) for k in PLUS) / 2, 1)
        self.assertEqual(sum(np.linalg.det(k) for k in MINUS) / 2, -1)

    def test_isometries_are_pure_seed_extensions_of_the_channels(self):
        rho = np.array([[.7, .1], [.1, .3]])
        for pair in (PLUS, MINUS):
            v = isometry(pair)
            np.testing.assert_allclose(v.T @ v, I, atol=3e-16)
            whole = v @ rho @ v.T
            output = np.trace(whole.reshape(2, 2, 2, 2), axis1=1, axis2=3)
            np.testing.assert_allclose(output, channel(pair, rho), atol=3e-16)

    def test_bell_context_gives_orthogonal_outputs_with_identical_local_records(self):
        first, second = witness_states()
        np.testing.assert_array_equal(first, (np.eye(4) + np.kron(J, J)) / 4)
        np.testing.assert_array_equal(second, (np.eye(4) - np.kron(J, J)) / 4)
        np.testing.assert_array_equal(first @ second, np.zeros((4, 4)))
        self.assertEqual(np.abs(np.linalg.eigvalsh(first - second)).sum() / 2, 1)
        for a in (I, X, Z):
            for b in (I, X, Z):
                self.assertEqual(np.trace(np.kron(a, b) @ (first - second)), 0)
        effect = (np.eye(4) + np.kron(J, J)) / 2
        np.testing.assert_array_equal(effect @ effect, effect)
        self.assertEqual(np.trace(effect @ first), 1)
        self.assertEqual(np.trace(effect @ second), 0)

    def test_no_environment_rotation_relates_the_two_dilations(self):
        # Vanishing cross Gram proves zero overlap for EVERY environment matrix.
        cross = np.array([[np.trace(a.T @ b) for b in MINUS] for a in PLUS])
        np.testing.assert_array_equal(cross, np.zeros((2, 2)))
        vp, vm = isometry(PLUS), isometry(MINUS)
        for theta in (0., .3, 1.2):
            rotation = np.cos(theta) * I + np.sin(theta) * J
            for env in (rotation, rotation @ Z):
                self.assertAlmostEqual(np.linalg.norm(vm - np.kron(I, env) @ vp), 2)

    def test_true_environment_gauge_change_preserves_all_matrix_directions(self):
        env = (3 * I + 4 * J) / 5
        changed = tuple(sum(env[e, f] * PLUS[f] for f in range(2)) for e in range(2))
        np.testing.assert_allclose(isometry(changed), np.kron(I, env) @ isometry(PLUS), atol=2e-16)
        for basis in (I, X, Z, J):
            np.testing.assert_allclose(channel(changed, basis), channel(PLUS, basis), atol=2e-16)

    def test_complex_local_y_state_already_distinguishes_the_channels(self):
        y = 1j * J
        rho = (I + y) / 2
        first, second = channel(PLUS, rho), channel(MINUS, rho)
        np.testing.assert_array_equal(first, rho)
        np.testing.assert_array_equal(second, (I - y) / 2)
        self.assertEqual(np.abs(np.linalg.eigvalsh(first - second)).sum() / 2, 1)

    def test_own_real_bell_probe_calibrates_missing_parameter_for_general_channels(self):
        rng = np.random.default_rng(174)
        bell = np.outer([1., 0., 0., 1.], [1., 0., 0., 1.]) / 2
        effect = (np.eye(4) + np.kron(J, J)) / 2
        for _ in range(12):
            dilation, _ = np.linalg.qr(rng.normal(size=(8, 2)))
            kraus = tuple(dilation.reshape(2, 4, 2)[:, j, :] for j in range(4))
            apply = lambda m: sum(k @ m @ k.T for k in kraus)
            kappa = sum(np.linalg.det(k) for k in kraus)
            joint = sum(np.kron(k, I) @ bell @ np.kron(k.T, I) for k in kraus)
            calibrated = 2 * np.trace(effect @ joint) - 1
            self.assertAlmostEqual(calibrated, kappa, places=13)
            # These four matrix directions determine action with any partner.
            for i, j in np.ndindex(2, 2):
                matrix_unit = np.zeros((2, 2))
                matrix_unit[i, j] = 1
                reconstructed = sum(np.trace(s.T @ matrix_unit) * apply(s) / 2
                                    for s in (I, X, Z))
                reconstructed += np.trace(J.T @ matrix_unit) * calibrated * J / 2
                np.testing.assert_allclose(reconstructed, apply(matrix_unit), atol=8e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ProcessContextAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 174,
        "plus_unnormalized_kraus": ["I", "J"],
        "minus_unnormalized_kraus": ["X", "Z"],
        "each_kraus_weight": "1/2",
        "all_local_real_state_outputs": "I/2",
        "real_bell_context_trace_distance_exact": "1",
        "environment_orthogonal_alignment_frobenius_distance_exact": "2",
        "kappa_plus": 1,
        "kappa_minus": -1,
        "real_state_purification_uniqueness_is_refuted": False,
        "local_only_process_description_is_composition_complete": False,
        "complete_contextual_equivalence_is_well_defined_in_real_theory": True,
        "own_real_bell_probe_calibrates_kappa_without_actual_future_partner": True,
        "calibration_probability": "p_plus=(1+kappa)/2",
        "cognition_forces_local_tomography": False,
        "global_ideal_effect_is_claimed_free_in_original_noisy_interface": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("process_context_audit_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
