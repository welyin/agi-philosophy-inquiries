"""Construction route C0: one declared quantum process kernel and its controls.

The matrix state space, Born rule, tensor composition, gates and graph are
working inputs. This integrates inherited interfaces; it does not derive
quantum necessity, physical spacetime, gravity or complete SoCA cognition.
"""

import argparse
import copy
import json
import math
import sys
import unittest
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "research_cognition_physics"))

from common_orientation_structure import encode_state
from complex_control_from_reference import real_lift
from independent_source_alignment import keep_systems
from bilocal_record_tomography import embed_operator
from quantum_interface_audit import IDENTITY, PAULI_X, PAULI_Y, PAULI_Z


HADAMARD = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
CNOT = np.array([[1., 0., 0., 0.], [0., 1., 0., 0.],
                 [0., 0., 0., 1.], [0., 0., 1., 0.]])
SWAP = np.array([[1., 0., 0., 0.], [0., 0., 1., 0.],
                 [0., 1., 0., 0.], [0., 0., 0., 1.]])
EDGES = {(0, 1), (1, 2), (1, 3)}
NAMES = ("A", "B", "C", "M")


def zero_state():
    density = np.zeros((16, 16), dtype=complex)
    density[0, 0] = 1
    return density


def binary_effect(outcome):
    return (IDENTITY + (1 if outcome == 0 else -1) * PAULI_Z) / 2


def checked_probability(density, effect):
    value = np.trace(effect @ density)
    if abs(value.imag) > 1e-12 or value.real < -1e-12 or value.real > 1 + 1e-12:
        raise ValueError("Invalid Born probability.")
    return float(np.clip(value.real, 0, 1))


@dataclass
class Whole:
    density: np.ndarray = field(default_factory=zero_state)
    belief: np.ndarray | None = None
    weight: float = 1.0
    history: list = field(default_factory=list)
    pointer_used: bool = False
    resources: dict = field(default_factory=lambda: {
        "logical_quantum_slots": 4, "initial_blank_record_slots": 1,
        "unitary_calls": 0, "two_slot_calls": 0, "pointer_reads": 0,
        "formal_record_bits": 0, "feedback_message_bits": 0, "resets": 0,
    })

    def __post_init__(self):
        self.density = np.array(self.density, dtype=complex, copy=True)
        if self.density.shape != (16, 16):
            raise ValueError("C0 uses four fixed logical slots.")
        if not np.allclose(self.density, self.density.conj().T, atol=1e-12):
            raise ValueError("State must be Hermitian.")
        if np.linalg.eigvalsh(self.density).min() < -1e-12 or abs(np.trace(self.density) - 1) > 1e-12:
            raise ValueError("State must be positive and normalized.")
        self.belief = np.array(self.density if self.belief is None else self.belief, dtype=complex, copy=True)
        if self.belief.shape != (16, 16) or not np.allclose(self.belief, self.belief.conj().T, rtol=0, atol=1e-12):
            raise ValueError("Prediction model must be a Hermitian state on the same slots.")
        if np.linalg.eigvalsh(self.belief).min() < -1e-12 or abs(np.trace(self.belief) - 1) > 1e-12:
            raise ValueError("Prediction model must be positive and normalized.")
        self.encoded = encode_state(self.density)

    def effect(self, observable, targets):
        return embed_operator(observable, tuple(targets), 4)

    def probability(self, slot, outcome=0, predicted=False):
        density = self.belief if predicted else self.density
        return checked_probability(density, self.effect(binary_effect(outcome), (slot,)))

    def expectation(self, observable, targets):
        return float(np.trace(self.effect(observable, targets) @ self.density).real)

    def apply(self, gate, targets, label):
        targets = tuple(targets)
        if len(targets) not in (1, 2) or len(set(targets)) != len(targets):
            raise ValueError("Use one slot or one permitted pair.")
        if any(slot < 0 or slot >= 4 for slot in targets):
            raise ValueError("Unknown slot.")
        if len(targets) == 2 and tuple(sorted(targets)) not in EDGES:
            raise ValueError("This pair is not an edge of the declared C0 graph.")
        gate = np.asarray(gate)
        dimension = 2 ** len(targets)
        if gate.shape != (dimension, dimension) or not np.allclose(gate.conj().T @ gate, np.eye(dimension), rtol=0, atol=1e-12):
            raise ValueError("A unitary gate of matching size is required.")
        expanded = self.effect(gate, targets)
        self.density = expanded @ self.density @ expanded.conj().T
        self.belief = expanded @ self.belief @ expanded.conj().T
        real_gate = real_lift(expanded)
        self.encoded = real_gate @ self.encoded @ real_gate.T
        self.resources["unitary_calls"] += 1
        self.resources["two_slot_calls"] += int(len(targets) == 2)
        self.history.append({"kind": "unitary", "label": label, "slots": [NAMES[slot] for slot in targets]})

    def record_B(self):
        if self.pointer_used or self.probability(3, 0) < 1 - 1e-12:
            raise ValueError("The blank record slot is unavailable; no free reset is allowed.")
        self.apply(CNOT, (1, 3), "write B into M")
        self.pointer_used = True
        branches = []
        for outcome in (0, 1):
            projector = self.effect(binary_effect(outcome), (3,))
            probability = checked_probability(self.density, projector)
            prediction = checked_probability(self.belief, projector)
            if probability == 0:
                continue
            if prediction == 0:
                raise ValueError("The predictor assigned zero probability to an actual branch.")
            branch = copy.deepcopy(self)
            branch.density = projector @ self.density @ projector / probability
            branch.belief = projector @ self.belief @ projector / prediction
            real_projector = real_lift(projector)
            branch.encoded = real_projector @ self.encoded @ real_projector.T / probability
            branch.weight *= probability
            branch.resources["pointer_reads"] += 1
            branch.resources["formal_record_bits"] += 1
            branch.history.append({"kind": "record", "slot": "M", "outcome": outcome,
                                   "predicted_probability": prediction, "branch_probability": probability})
            branches.append(branch)
        return branches

    def feedback_A_to_zero(self):
        if not any(event["kind"] == "record" for event in self.history):
            raise ValueError("Feedback requires a retained observation.")
        keep_success = self.probability(0, 0, predicted=True)
        flip_success = 1 - keep_success
        chosen = "flip" if flip_success > keep_success else "keep"
        self.resources["feedback_message_bits"] += 1
        self.history.append({"kind": "decision", "target": "A=0", "keep_prediction": keep_success,
                             "flip_prediction": flip_success, "action": chosen})
        if chosen == "flip":
            self.apply(PAULI_X, (0,), "feedback flip A")

    def validate(self):
        if np.linalg.eigvalsh(self.density).min() < -1e-12:
            raise ValueError("Negative state eigenvalue.")
        if abs(np.trace(self.density) - 1) > 1e-12:
            raise ValueError("State normalization failed.")
        if not np.allclose(self.encoded, encode_state(self.density), atol=1e-12):
            raise ValueError("The shared-J real encoding no longer matches the complex whole.")


def bell_world():
    world = Whole()
    world.apply(HADAMARD, (0,), "prepare A plus")
    world.apply(CNOT, (0, 1), "entangle A B")
    return world


def interference(theta):
    world = Whole()
    world.apply(HADAMARD, (0,), "first interference rotation")
    world.apply(np.diag([np.exp(-0.5j * theta), np.exp(0.5j * theta)]), (0,), "phase")
    world.apply(HADAMARD, (0,), "recombine")
    world.validate()
    return world.probability(0, 1)


def signal_records(bit):
    world = Whole()
    if bit:
        world.apply(PAULI_X, (0,), "prepare input bit")
    records = [world.probability(2, 1)]
    world.apply(SWAP, (0, 1), "propagation tick 1")
    records.append(world.probability(2, 1))
    world.apply(SWAP, (1, 2), "propagation tick 2")
    records.append(world.probability(2, 1))
    world.validate()
    return records


def run_suite():
    bell = bell_world()
    diagonal_plus = (PAULI_Z + PAULI_X) / np.sqrt(2)
    diagonal_minus = (PAULI_Z - PAULI_X) / np.sqrt(2)
    chsh = (bell.expectation(np.kron(PAULI_Z, diagonal_plus), (0, 1))
            + bell.expectation(np.kron(PAULI_Z, diagonal_minus), (0, 1))
            + bell.expectation(np.kron(PAULI_X, diagonal_plus), (0, 1))
            - bell.expectation(np.kron(PAULI_X, diagonal_minus), (0, 1)))
    branches = bell.record_B()
    feedback = []
    for branch in branches:
        branch.feedback_A_to_zero()
        branch.validate()
        feedback.append({"weight": branch.weight, "A_zero_probability": branch.probability(0),
                         "history": branch.history, "resources": branch.resources})
    return {
        "baseline": "C0", "route_round": "00", "kind": "conditional integration baseline, not a new physical derivation",
        "slots": list(NAMES), "edges": [[NAMES[first], NAMES[second]] for first, second in sorted(EDGES)],
        "interference": [{"phase": theta, "A_one_probability": interference(theta)} for theta in (0., np.pi / 2, np.pi)],
        "bell_CHSH": chsh,
        "feedback_branches": feedback,
        "all_feedback_branches_retained": True,
        "signal": {"input_zero_C_one": signal_records(0), "input_one_C_one": signal_records(1),
                   "earliest_C_distinguishing_tick": 2,
                   "tick_interpretation": "declared gate-layer index, not derived physical time"},
        "real_equivalent_whole_dimension": 32,
        "real_encoding_uses_one_global_reference_not_independent_local_encodings": True,
        "classical_predictor_uses_known_preparation_model_not_unknown_state_self_readout": True,
        "physical_energy_clock_calibration_and_reset_recycling_completed": False,
        "full_SoCA_quantum_necessity_or_gravity_derived": False,
        "samples_collected": 0,
    }


class ConstructionBaselineTests(unittest.TestCase):
    def test_interference_fringe_uses_the_same_whole_state_kernel(self):
        for theta in (0., 0.37, np.pi / 2, np.pi):
            self.assertAlmostEqual(interference(theta), np.sin(theta / 2) ** 2)

    def test_bell_correlations_and_entanglement_are_consistent(self):
        world = bell_world()
        self.assertAlmostEqual(world.expectation(np.kron(PAULI_X, PAULI_X), (0, 1)), 1)
        self.assertAlmostEqual(world.expectation(np.kron(PAULI_Z, PAULI_Z), (0, 1)), 1)
        pair = keep_systems(world.density, (0, 1), 4)
        partial_transpose = pair.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
        self.assertAlmostEqual(np.linalg.eigvalsh(partial_transpose).min(), -0.5)
        world.validate()

    def test_prediction_observation_feedback_closes_on_both_records(self):
        branches = bell_world().record_B()
        self.assertEqual(len(branches), 2)
        self.assertAlmostEqual(sum(branch.weight for branch in branches), 1)
        actions = []
        for branch in branches:
            branch.feedback_A_to_zero()
            branch.validate()
            self.assertAlmostEqual(branch.weight, 0.5)
            self.assertAlmostEqual(branch.probability(0), 1)
            self.assertEqual(branch.resources["pointer_reads"], 1)
            self.assertEqual(branch.resources["feedback_message_bits"], 1)
            actions.extend(event["action"] for event in branch.history if event["kind"] == "decision")
        self.assertEqual(actions, ["keep", "flip"])

    def test_nonzero_rare_record_is_not_discarded(self):
        density = zero_state()
        density[0, 0] = 1 - 1e-15
        density[4, 4] = 1e-15
        branches = Whole(density).record_B()
        self.assertEqual(len(branches), 2)
        self.assertGreater(branches[1].weight, 0)
        self.assertAlmostEqual(branches[1].weight / 1e-15, 1)
        for branch in branches:
            branch.validate()

    def test_unread_record_can_be_reversed_but_public_reading_changes_the_process(self):
        coherent = Whole()
        coherent.apply(HADAMARD, (1,), "prepare B plus")
        coherent.apply(CNOT, (1, 3), "couple B M")
        self.assertAlmostEqual(coherent.expectation(PAULI_X, (1,)), 0)
        self.assertEqual(coherent.density.shape, (16, 16))
        coherent.apply(CNOT, (1, 3), "undo coupling")
        self.assertAlmostEqual(coherent.expectation(PAULI_X, (1,)), 1)
        measured = Whole()
        measured.apply(HADAMARD, (1,), "prepare B plus")
        branches = measured.record_B()
        average_x = 0
        for branch in branches:
            branch.apply(CNOT, (1, 3), "undo after public record")
            average_x += branch.weight * branch.expectation(PAULI_X, (1,))
            branch.validate()
        self.assertAlmostEqual(average_x, 0)

    def test_declared_schedule_has_no_earlier_A_to_C_signal(self):
        np.testing.assert_allclose(signal_records(0), [0, 0, 0], atol=1e-12)
        np.testing.assert_allclose(signal_records(1), [0, 0, 1], atol=1e-12)

    def test_no_undeclared_edge_or_free_pointer_reset_is_allowed(self):
        with self.assertRaises(ValueError):
            Whole().apply(CNOT, (0, 2), "illegal shortcut")
        branch = bell_world().record_B()[0]
        with self.assertRaises(ValueError):
            branch.record_B()
        self.assertEqual(branch.resources["resets"], 0)

    def test_random_preparation_mixture_is_preserved_by_the_same_operation(self):
        first = Whole()
        second = Whole()
        second.apply(PAULI_X, (1,), "prepare alternative")
        mixed = Whole(0.3 * first.density + 0.7 * second.density)
        for world in (first, second, mixed):
            world.apply(CNOT, (1, 2), "shared action")
        np.testing.assert_allclose(mixed.density, 0.3 * first.density + 0.7 * second.density, atol=1e-12)
        mixed.validate()

    def test_real_interface_agrees_without_physically_erasing_the_whole(self):
        world = Whole()
        world.apply(HADAMARD, (0,), "prepare plus")
        world.apply(np.diag([1., 1j]), (0,), "expose complex direction")
        self.assertGreater(np.linalg.norm(world.density.imag), 0.5)
        for axis in (PAULI_X, PAULI_Z):
            effect = world.effect((IDENTITY + axis) / 2, (0,))
            self.assertAlmostEqual(np.trace(effect @ world.density).real,
                                   np.trace(effect @ world.density.real).real)
        self.assertAlmostEqual(world.expectation(PAULI_Y, (0,)), 1)
        world.validate()

    def test_invalid_states_and_nonunitary_commands_are_rejected(self):
        invalid = zero_state()
        invalid[0, 0] = 1.1
        invalid[1, 1] = -0.1
        with self.assertRaises(ValueError):
            Whole(invalid)
        with self.assertRaises(ValueError):
            Whole().apply(binary_effect(0), (0,), "not a unitary")

    def test_report_keeps_assumptions_and_resources_explicit(self):
        report = run_suite()
        self.assertAlmostEqual(report["bell_CHSH"], 2 * np.sqrt(2))
        self.assertFalse(report["full_SoCA_quantum_necessity_or_gravity_derived"])
        self.assertTrue(report["all_feedback_branches_retained"])
        self.assertEqual(report["signal"]["earliest_C_distinguishing_tick"], 2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ConstructionBaselineTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    report = run_suite()
    report["automated_checks"] = {"run": tests.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("baseline_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"baseline": "C0", "tests": tests.testsRun, "CHSH": report["bell_CHSH"],
                      "feedback_success": sum(branch["weight"] * branch["A_zero_probability"] for branch in report["feedback_branches"]),
                      "A_to_C_ticks": 2}, indent=2))


if __name__ == "__main__":
    main()