"""Round 130: exact intake on the original real preparation plane.

Retain the old host reference, route the newcomer target into the active
interface, and keep the unused private reference in the whole. This requires
no ideal orientation measurement. It does not preserve arbitrary hidden Y
content when the newcomer is outside the real preparation promise.
"""

import argparse
import json
import math
import unittest
from pathlib import Path

import numpy as np

from approximate_subject_joining import join_output
from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from encoded_composition_audit import independent_encoding
from flagged_subject_joining import random_state
from operational_effect_closure import compiled_copy_gate
from quantum_interface_audit import (
    ALPHA, PAULI_X, PAULI_Y, PAULI_Z, apply_branch, branch_kraus,
    effective_readout,
)


def route_reference_last(host_dimension, newcomer_dimension):
    """Register order RA,RB,A,B -> RA,A,B,RB; a relabeling, not a free SWAP gate."""
    size = 4 * host_dimension * newcomer_dimension
    order = np.arange(size).reshape(2, 2, host_dimension, newcomer_dimension)
    order = order.transpose(0, 2, 3, 1).ravel()
    return np.eye(size)[order]


def routed_whole(host, newcomer):
    route = route_reference_last(len(host), len(newcomer))
    return route @ independent_encoding((host, newcomer)) @ route.T


def retained_reference_join(host, newcomer):
    whole = routed_whole(host, newcomer)
    size = len(whole) // 2
    return np.einsum("abcb->ac", whole.reshape(size, 2, size, 2))


def old_local_update(encoded, host_kraus, newcomer_kraus):
    return sum(
        lift @ encoded @ lift.T
        for a in host_kraus for b in newcomer_kraus
        for lift in (real_lift(np.kron(a, b)),)
    )


class RealInterfaceJoiningTests(unittest.TestCase):
    def test_real_newcomer_joins_arbitrary_unknown_host_exactly(self):
        rng = np.random.default_rng(130)
        for da, db in ((2, 2), (3, 2), (2, 3)):
            host = random_state(rng, da)
            raw = rng.normal(size=(db, db))
            newcomer = raw @ raw.T
            newcomer /= np.trace(newcomer)
            np.testing.assert_allclose(retained_reference_join(host, newcomer),
                                       encode_state(np.kron(host, newcomer)), atol=1e-16)

    def test_unused_reference_stays_in_whole_as_a_product_for_promised_inputs(self):
        host = (np.eye(2) + .7 * PAULI_Y) / 2
        new = (np.eye(2) + .4 * PAULI_X + .5 * PAULI_Z) / 2
        whole = routed_whole(host, new)
        np.testing.assert_allclose(whole, np.kron(encode_state(np.kron(host, new)), np.eye(2) / 2), atol=0.)
        route = route_reference_last(2, 2)
        np.testing.assert_array_equal(route.T @ whole @ route, independent_encoding((host, new)))

    def test_actual_old_noisy_instruments_preserve_entire_selected_updates(self):
        host = (np.eye(2) + .3 * PAULI_X + .6 * PAULI_Z) / 2
        new = (np.eye(2) - .5 * PAULI_X + .2 * PAULI_Z) / 2
        for a in (0, 1):
            for b in (0, 1):
                ka = branch_kraus(.31, a, .5)
                kb = branch_kraus(-.27 + .2 * a, b, .7)
                actual = old_local_update(retained_reference_join(host, new), ka, kb)
                pa = apply_branch(host, .31, a, .5)
                pb = apply_branch(new, -.27 + .2 * a, b, .7)
                np.testing.assert_allclose(actual, encode_state(np.kron(pa, pb)), atol=2e-16)

    def test_all_adaptive_old_real_histories_ignore_unavailable_y_component(self):
        original = (np.eye(2) + .3 * PAULI_X + .6 * PAULI_Y + .4 * PAULI_Z) / 2
        shadow = original.real.astype(complex)
        pairs = [(original, shadow, ())]
        for step in range(4):
            updated = []
            for first, second, history in pairs:
                angle = .2 * step + (.3 if history and history[-1] else -.4)
                for result in (0, 1):
                    a = apply_branch(first, angle, result, .5)
                    b = apply_branch(second, angle, result, .5)
                    self.assertAlmostEqual(np.trace(a).real, np.trace(b).real)
                    np.testing.assert_allclose(a.real, b.real, atol=2e-16)
                    updated.append((a, b, history + (result,)))
            pairs = updated

    def test_complex_newcomer_loses_y_only_in_active_interface_not_whole(self):
        host = np.diag([1., 0.])
        outputs = [retained_reference_join(host, (np.eye(2) + s * PAULI_Y) / 2) for s in (-1, 1)]
        np.testing.assert_array_equal(*outputs)
        first, second = [routed_whole(host, (np.eye(2) + s * PAULI_Y) / 2) for s in (-1, 1)]
        self.assertAlmostEqual(np.abs(np.linalg.eigvalsh(first - second)).sum() / 2, 1.)

    def test_generic_approximation_introduces_visible_avoidable_old_readout_error(self):
        host, new = np.eye(2) / 2, np.diag([1., 0.])
        effect = real_lift(np.kron(np.eye(2), effective_readout(0., 1, .5)))
        exact = retained_reference_join(host, new)
        approximate = join_output(host, new)
        gap = np.trace((exact - approximate) @ effect).real
        self.assertAlmostEqual(gap, ALPHA / 12)

    def test_old_compiled_real_gate_creates_the_exact_shared_bell_target(self):
        host = np.diag([1., 0.])  # Pointer A
        new = (np.eye(2) + PAULI_X) / 2  # Control B
        gate = real_lift(compiled_copy_gate())
        output = gate @ retained_reference_join(host, new) @ gate.T
        bell = np.outer([1., 0., 0., 1.], [1., 0., 0., 1.]) / 2
        np.testing.assert_allclose(output, encode_state(bell), atol=3e-16)
        for setting in (0., math.pi / 2):
            correlation = 0.
            for a in (0, 1):
                for b in (0, 1):
                    effect = real_lift(np.kron(effective_readout(setting, a, .5),
                                              effective_readout(setting, b, .5)))
                    correlation += (2*a-1)*(2*b-1)*np.trace(output @ effect).real
            self.assertAlmostEqual(correlation, (ALPHA / 2) ** 2)

    def test_active_output_is_exactly_the_real_shadow_for_any_newcomer(self):
        rng = np.random.default_rng(230)
        for _ in range(6):
            host, new = random_state(rng, 2), random_state(rng, 2)
            np.testing.assert_allclose(retained_reference_join(host, new),
                                       encode_state(np.kron(host, new.real)), atol=2e-16)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RealInterfaceJoiningTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 130,
        "newcomer_promise": "Unknown density matrix real in the existing X-Z preparation basis",
        "host_may_be_arbitrary_unknown_complex_logical_state": True,
        "exact_common_output": "E(rho_A tensor sigma_B)",
        "orientation_filter_or_new_readout_needed": False,
        "unused_private_reference_kept_in_whole": True,
        "register_relabeling_counted_as_a_free_physical_swap": False,
        "old_real_instruments_and_finite_adaptive_histories_preserved": True,
        "arbitrary_newcomer_active_output": "E(rho_A tensor Re(sigma_B))",
        "unavailable_y_content_claimed_preserved_in_active_output": False,
        "avoidable_error_of_round127_for_eta_half_diagnostic": ALPHA / 12,
        "old_compiled_joint_gate_produces_bell_target_exactly": True,
        "cognitive_origin_of_joint_gate_or_reference_inferred": False,
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("real_interface_joining_results.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
