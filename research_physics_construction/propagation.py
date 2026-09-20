"""Construction 01 / C1: local exchange propagation and delayed feedback.

Extend the C0 whole without changing its source. The continuous Hamiltonian,
ideal instantaneous control/readout and per-edge classical delay are working
inputs, not derived clocks or relativistic propagation laws. Exact finite
matrix evolution is a reference for a separately costed split-step algorithm.
"""

import argparse
import copy
import json
import math
import unittest
from pathlib import Path

import numpy as np

from baseline import (Whole, HADAMARD, CNOT, IDENTITY, PAULI_X, PAULI_Y,
                      EDGES, NAMES, binary_effect, checked_probability, real_lift)


EXCHANGE = (np.kron(PAULI_X, PAULI_X) + np.kron(PAULI_Y, PAULI_Y)) / 2
DEFAULT_HOP_DELAY = math.pi / (4 * math.sqrt(2))


def unitary_from_hamiltonian(hamiltonian, duration):
    values, vectors = np.linalg.eigh(hamiltonian)
    return (vectors * np.exp(-1j * duration * values)) @ vectors.conj().T


def hamiltonian_parts(coupling=1.0):
    world = Whole()
    return tuple(coupling * world.effect(EXCHANGE, edge) for edge in ((0, 1), (1, 2)))


def transfer_probability(duration, coupling=1.0):
    return float(np.sin(coupling * duration / np.sqrt(2)) ** 4)


def threshold_arrival(threshold, coupling=1.0):
    if not 0 < threshold <= 1 or coupling <= 0:
        raise ValueError("Require 0 < probability threshold <= 1 and positive coupling.")
    return float(np.sqrt(2) * np.arcsin(threshold ** 0.25) / coupling)


def split_propagator(duration, steps, symmetric=False, coupling=1.0):
    if not isinstance(steps, int) or steps < 1:
        raise ValueError("Use a positive integer step count.")
    first, second = hamiltonian_parts(coupling)
    delta = duration / steps
    if symmetric:
        half = unitary_from_hamiltonian(first, delta / 2)
        single = half @ unitary_from_hamiltonian(second, delta) @ half
    else:
        single = unitary_from_hamiltonian(second, delta) @ unitary_from_hamiltonian(first, delta)
    return np.linalg.matrix_power(single, steps)


class PropagatingWhole(Whole):
    def __init__(self, *args, coupling=1.0, hop_delay=DEFAULT_HOP_DELAY, **kwargs):
        super().__init__(*args, **kwargs)
        if not np.isfinite(coupling) or coupling <= 0 or not np.isfinite(hop_delay) or hop_delay <= 0:
            raise ValueError("Coupling and classical hop delay must be finite and positive.")
        self.coupling = float(coupling)
        self.hop_delay = float(hop_delay)
        self.time = 0.0
        self.hamiltonian = sum(hamiltonian_parts(self.coupling))
        self.remote_belief = self.belief.copy()
        self.pending_message = None
        self.delivered_records = []
        self.feedback_done = False
        self.resources.update({"propagation_intervals": 0, "evolution_parameter": 0.0,
                               "active_edge_parameter_time": 0.0,
                               "classical_bit_hops": 0, "message_in_flight": 0})

    def apply(self, gate, targets, label):
        super().apply(gate, targets, label)
        expanded = self.effect(gate, targets)
        self.remote_belief = expanded @ self.remote_belief @ expanded.conj().T
        self.history[-1]["time"] = self.time

    def _drift(self, duration):
        if duration <= 0:
            return
        gate = unitary_from_hamiltonian(self.hamiltonian, duration)
        self.density = gate @ self.density @ gate.conj().T
        self.belief = gate @ self.belief @ gate.conj().T
        self.remote_belief = gate @ self.remote_belief @ gate.conj().T
        lifted = real_lift(gate)
        self.encoded = lifted @ self.encoded @ lifted.T
        start = self.time
        self.time += duration
        self.resources["propagation_intervals"] += 1
        self.resources["evolution_parameter"] += duration
        self.resources["active_edge_parameter_time"] += 2 * duration
        self.history.append({"kind": "propagation", "start": start, "end": self.time,
                             "active_edges": [["A", "B"], ["B", "C"]],
                             "coupling": self.coupling, "solver": "exact finite matrix reference"})

    def advance(self, duration):
        if not np.isfinite(duration) or duration < 0:
            raise ValueError("Forward elapsed parameter must be finite and nonnegative.")
        destination = self.time + duration
        message = self.pending_message
        if message is not None:
            while message["completed_hops"] < 2:
                hop = message["completed_hops"] + 1
                due = message["sent_at"] + hop * self.hop_delay
                if due > destination:
                    break
                self._drift(due - self.time)
                source, target = message["path"][hop - 1:hop + 1]
                self.history.append({"kind": "message_hop", "time": due,
                                     "source": NAMES[source], "target": NAMES[target],
                                     "record_id": message["record_id"]})
                self.resources["classical_bit_hops"] += 1
                message["completed_hops"] = hop
                if hop == 2:
                    self.delivered_records.append({"outcome": message["outcome"], "received_at": due,
                                                   "record_id": message["record_id"]})
                    self.remote_belief = self.belief.copy()
                    self.resources["message_in_flight"] = 0
        self._drift(destination - self.time)

    def record_B(self):
        if self.pending_message is not None:
            raise ValueError("C1 contains only one allocated record and message slot.")
        branches = super().record_B()
        remote_prior = sum(self.effect(binary_effect(outcome), (3,)) @ self.remote_belief
                           @ self.effect(binary_effect(outcome), (3,)) for outcome in (0, 1))
        path = (3, 1, 0)
        if any(tuple(sorted(edge)) not in EDGES for edge in zip(path, path[1:])):
            raise ValueError("The declared classical path must follow existing edges.")
        for branch in branches:
            record = branch.history[-1]
            record["time"] = self.time
            record["visibility"] = "local M record; not yet delivered to A"
            branch.remote_belief = remote_prior.copy()
            branch.pending_message = {"record_id": 0, "outcome": record["outcome"], "path": path,
                                      "sent_at": self.time, "completed_hops": 0}
            branch.resources["message_in_flight"] = 1
            branch.history.append({"kind": "message_sent", "time": self.time,
                                   "record_id": 0, "path": [NAMES[slot] for slot in path],
                                   "earliest_delivery": self.time + 2 * self.hop_delay})
        return branches

    def controller_probability(self):
        return checked_probability(self.remote_belief, self.effect(binary_effect(0), (0,)))

    def feedback_A_to_zero(self, stale=False):
        if not self.delivered_records or self.feedback_done:
            raise ValueError("Use a newly delivered record; early or repeated feedback is forbidden.")
        keep_probability = self.controller_probability()
        if stale:
            chosen = "flip" if self.delivered_records[-1]["outcome"] == 1 else "keep"
        else:
            chosen = "flip" if 1 - keep_probability > keep_probability else "keep"
        self.resources["feedback_message_bits"] += 1
        self.feedback_done = True
        self.history.append({"kind": "decision", "time": self.time, "action": chosen,
                             "keep_prediction": keep_probability, "flip_prediction": 1 - keep_probability,
                             "policy": "stale zero-delay action" if stale else "propagated conditional belief",
                             "record_id": self.delivered_records[-1]["record_id"]})
        if chosen == "flip":
            self.apply(PAULI_X, (0,), "delivered-record feedback")


def prepared_bell(hop_delay=DEFAULT_HOP_DELAY):
    world = PropagatingWhole(hop_delay=hop_delay)
    world.apply(HADAMARD, (0,), "known A plus preparation")
    world.apply(CNOT, (0, 1), "known Bell preparation")
    return world


def feedback_report(hop_delay):
    records = []
    for branch in prepared_bell(hop_delay).record_B():
        branch.advance(2 * hop_delay)
        stale = copy.deepcopy(branch)
        conditional_best = max(branch.controller_probability(), 1 - branch.controller_probability())
        branch.feedback_A_to_zero()
        stale.feedback_A_to_zero(stale=True)
        branch.validate()
        stale.validate()
        records.append({"outcome": branch.delivered_records[-1]["outcome"], "weight": branch.weight,
                        "calibrated_success": branch.probability(0), "stale_success": stale.probability(0),
                        "conditional_best_success": conditional_best,
                        "events": branch.history, "resources": branch.resources})
    transferred = transfer_probability(2 * hop_delay)
    return {"hop_delay": hop_delay, "arrival_parameter": 2 * hop_delay, "branches": records,
            "calibrated_average_success": sum(item["weight"] * item["calibrated_success"] for item in records),
            "stale_average_success": sum(item["weight"] * item["stale_success"] for item in records),
            "best_no_message_action_at_same_parameter": (1 + transferred) / 2}


def interference_report(phase):
    world = PropagatingWhole()
    world.apply(PAULI_X, (1,), "known central excitation")
    half_duration = math.pi / (2 * math.sqrt(2))
    world.advance(half_duration)
    world.apply(np.diag([1, np.exp(1j * phase)]), (2,), "local phase on C")
    world.advance(half_duration)
    world.validate()
    return {"phase": phase, "B_one_probability": world.probability(1, 1),
            "elapsed_parameter": world.time, "resources": world.resources,
            "events": world.history}


def feedback_approximation_report(hop_delay, steps=32):
    duration = 2 * hop_delay
    approximated_gate = split_propagator(duration, steps)
    first, second = hamiltonian_parts()
    bound = duration ** 2 * float(np.linalg.norm(first @ second - second @ first, ord=2)) / (2 * steps)
    rows = []
    for branch in prepared_bell(hop_delay).record_B():
        approximated_belief = approximated_gate @ branch.belief @ approximated_gate.conj().T
        approximated_keep = checked_probability(approximated_belief, branch.effect(binary_effect(0), (0,)))
        branch.advance(duration)
        exact_keep = branch.controller_probability()
        rows.append({"outcome": branch.delivered_records[-1]["outcome"],
                     "exact_keep_prediction": exact_keep, "approximate_keep_prediction": approximated_keep,
                     "decision_margin": abs(exact_keep - 0.5),
                     "same_action": (exact_keep < 0.5) == (approximated_keep < 0.5),
                     "action_equality_certified_by_bound": bound < abs(exact_keep - 0.5)})
    return {"arrival_parameter": duration, "steps": steps, "bond_pulses": 2 * steps,
            "analytic_operator_bound": bound, "branches": rows,
            "scope": "fixed initial record instrument and no additional quantum intervention during delivery"}


class PropagationTests(unittest.TestCase):
    def test_local_hamiltonian_conserves_excitation_and_excludes_the_record_slot(self):
        world = PropagatingWhole()
        excitation = sum(world.effect(binary_effect(1), (slot,)) for slot in (0, 1, 2))
        np.testing.assert_array_equal(world.hamiltonian, world.hamiltonian.conj().T)
        np.testing.assert_array_equal(world.hamiltonian @ excitation, excitation @ world.hamiltonian)
        pointer = world.effect(binary_effect(1), (3,))
        np.testing.assert_array_equal(world.hamiltonian @ pointer, pointer @ world.hamiltonian)

    def test_exact_transfer_curve_and_return_match_the_analytic_chain(self):
        for duration in (0., 0.1, 0.7, math.pi / math.sqrt(2), 2 * math.pi / math.sqrt(2)):
            world = PropagatingWhole()
            world.apply(PAULI_X, (0,), "input one")
            energy = np.trace(world.hamiltonian @ world.density).real
            world.advance(duration)
            self.assertAlmostEqual(world.probability(2, 1), transfer_probability(duration))
            self.assertAlmostEqual(np.trace(world.hamiltonian @ world.density).real, energy)
            self.assertAlmostEqual(world.probability(3, 0), 1)
            world.validate()

    def test_nonzero_tail_and_threshold_time_are_not_a_strict_light_cone(self):
        for duration in (0.01, 0.02, 0.04):
            probability = transfer_probability(duration)
            self.assertGreater(probability, 0)
            self.assertAlmostEqual(probability / duration ** 4, 0.25, delta=0.0002)
        for threshold in (1e-8, 1e-4, 0.01, 0.5, 1):
            self.assertAlmostEqual(transfer_probability(threshold_arrival(threshold)), threshold)
        with self.assertRaises(ValueError):
            threshold_arrival(0)

    def test_phase_interference_and_propagation_share_the_same_hamiltonian(self):
        for phase in (0., math.pi / 2, math.pi):
            report = interference_report(phase)
            self.assertAlmostEqual(report["B_one_probability"], (1 + math.cos(phase)) / 2)
            self.assertEqual(report["resources"]["unitary_calls"], 2)
            self.assertEqual(report["resources"]["propagation_intervals"], 2)

    def test_lie_split_has_the_commutator_bound_and_step_halving_order(self):
        first, second = hamiltonian_parts()
        duration = 0.8
        exact = unitary_from_hamiltonian(first + second, duration)
        commutator_norm = np.linalg.norm(first @ second - second @ first, ord=2)
        errors = []
        for steps in (4, 8, 16, 32):
            error = np.linalg.norm(split_propagator(duration, steps) - exact, ord=2)
            self.assertLessEqual(error, duration ** 2 * commutator_norm / (2 * steps) + 1e-12)
            errors.append(error)
        self.assertTrue(all(1.8 < coarse / fine < 2.2 for coarse, fine in zip(errors, errors[1:])))

    def test_symmetric_split_has_observed_second_order_not_a_claimed_global_certificate(self):
        duration = 0.8
        exact = unitary_from_hamiltonian(sum(hamiltonian_parts()), duration)
        errors = [np.linalg.norm(split_propagator(duration, steps, symmetric=True) - exact, ord=2)
                  for steps in (4, 8, 16, 32)]
        self.assertTrue(all(3.5 < coarse / fine < 4.5 for coarse, fine in zip(errors, errors[1:])))
        with self.assertRaises(ValueError):
            split_propagator(duration, 0)

    def test_controller_cannot_consume_a_record_before_two_edge_deliveries(self):
        branches = prepared_bell().record_B()
        self.assertAlmostEqual(sum(branch.weight for branch in branches), 1)
        for branch in branches:
            with self.assertRaises(ValueError):
                branch.feedback_A_to_zero()
            branch.advance(branch.hop_delay)
            self.assertEqual(branch.resources["classical_bit_hops"], 1)
            self.assertEqual(branch.delivered_records, [])
            with self.assertRaises(ValueError):
                branch.feedback_A_to_zero()
            branch.advance(branch.hop_delay)
            self.assertEqual(branch.resources["classical_bit_hops"], 2)
            self.assertEqual(len(branch.delivered_records), 1)
            branch.feedback_A_to_zero()
            with self.assertRaises(ValueError):
                branch.feedback_A_to_zero()

    def test_remote_prior_is_identical_across_unreceived_branches(self):
        branches = prepared_bell().record_B()
        for branch in branches:
            branch.advance(branch.hop_delay / 2)
        np.testing.assert_allclose(branches[0].remote_belief, branches[1].remote_belief, rtol=0, atol=1e-12)
        self.assertAlmostEqual(branches[0].controller_probability(), branches[1].controller_probability())
        self.assertGreater(abs(branches[0].probability(0) - branches[1].probability(0)), 0.9)

    def test_no_early_delivery_at_a_nearby_time_and_no_duplicate_delivery(self):
        branch = prepared_bell().record_B()[1]
        arrival = 2 * branch.hop_delay
        branch.advance(arrival - 1e-7)
        self.assertEqual(branch.delivered_records, [])
        with self.assertRaises(ValueError):
            branch.feedback_A_to_zero()
        branch.advance(arrival - branch.time)
        self.assertEqual(len(branch.delivered_records), 1)
        branch.advance(0.2)
        self.assertEqual(len(branch.delivered_records), 1)
        self.assertEqual(branch.resources["classical_bit_hops"], 2)
        branch.validate()

    def test_split_prediction_bound_is_smaller_than_both_feedback_decision_margins(self):
        for hop_delay in (DEFAULT_HOP_DELAY, math.pi / (3 * math.sqrt(2))):
            report = feedback_approximation_report(hop_delay)
            for branch in report["branches"]:
                self.assertTrue(branch["same_action"])
                self.assertTrue(branch["action_equality_certified_by_bound"])
                self.assertLessEqual(abs(branch["exact_keep_prediction"] - branch["approximate_keep_prediction"]),
                                     report["analytic_operator_bound"] + 1e-12)

    def test_delayed_feedback_matches_both_exact_rational_benchmarks(self):
        short = feedback_report(DEFAULT_HOP_DELAY)
        longer = feedback_report(math.pi / (3 * math.sqrt(2)))
        self.assertAlmostEqual(short["calibrated_average_success"], 7 / 8)
        self.assertAlmostEqual(short["stale_average_success"], 7 / 8)
        self.assertAlmostEqual(short["best_no_message_action_at_same_parameter"], 5 / 8)
        self.assertAlmostEqual(longer["calibrated_average_success"], 25 / 32)
        self.assertAlmostEqual(longer["stale_average_success"], 23 / 32)
        self.assertAlmostEqual(longer["best_no_message_action_at_same_parameter"], 25 / 32)

    def test_joint_records_remain_normalized_and_pointer_is_not_reset(self):
        report = feedback_report(DEFAULT_HOP_DELAY)
        for branch in report["branches"]:
            self.assertAlmostEqual(branch["weight"], 0.5)
            self.assertEqual(branch["resources"]["formal_record_bits"], 1)
            self.assertEqual(branch["resources"]["feedback_message_bits"], 1)
            self.assertEqual(branch["resources"]["resets"], 0)
            self.assertEqual(branch["resources"]["message_in_flight"], 0)
            sent = next(event for event in branch["events"] if event["kind"] == "message_sent")
            decision = next(event for event in branch["events"] if event["kind"] == "decision")
            self.assertGreaterEqual(decision["time"], sent["earliest_delivery"])

    def test_random_mixtures_evolve_affinely_and_shared_real_encoding_matches(self):
        first, second = PropagatingWhole(), PropagatingWhole()
        second.apply(PAULI_X, (0,), "alternative source")
        mixed = PropagatingWhole(0.3 * first.density + 0.7 * second.density)
        for world in (first, second, mixed):
            world.advance(0.37)
            world.validate()
        np.testing.assert_allclose(mixed.density, 0.3 * first.density + 0.7 * second.density, rtol=0, atol=1e-12)

    def test_time_and_coupling_inputs_are_validated(self):
        for kwargs in ({"coupling": 0}, {"hop_delay": -1}, {"hop_delay": float("nan")}):
            with self.assertRaises(ValueError):
                PropagatingWhole(**kwargs)
        with self.assertRaises(ValueError):
            PropagatingWhole().advance(-1)


def run_report():
    duration = 0.8
    first, second = hamiltonian_parts()
    exact = unitary_from_hamiltonian(first + second, duration)
    commutator_norm = float(np.linalg.norm(first @ second - second @ first, ord=2))
    split_rows = [{"steps": steps,
                   "lie_operator_error_diagnostic": float(np.linalg.norm(split_propagator(duration, steps) - exact, ord=2)),
                   "lie_analytic_error_bound": duration ** 2 * commutator_norm / (2 * steps),
                   "symmetric_operator_error_diagnostic": float(np.linalg.norm(split_propagator(duration, steps, True) - exact, ord=2)),
                   "lie_bond_pulses": 2 * steps, "symmetric_bond_pulses_uncompressed": 3 * steps,
                   "sequential_active_bond_duration_at_same_coupling": 2 * duration}
                  for steps in (4, 8, 16, 32)]
    return {
        "model_version": "C1", "route_round": "01", "parent": "C0 baseline.py unchanged",
        "research_scope": "sufficiency construction and consistency checks only; necessity questions recorded but not investigated",
        "inputs": {"Hamiltonian": "g/2*(XX_AB+YY_AB+XX_BC+YY_BC)", "g": 1,
                   "default_classical_hop_delay": DEFAULT_HOP_DELAY,
                   "classical_path": ["M", "B", "A"],
                   "control_and_projective_readout_durations": "ideal instantaneous interventions",
                   "time_units": "declared evolution parameter, not a reconstructed physical clock"},
        "transfer_formula": "P_C(1|A=1)=sin(g*t/sqrt(2))^4; vacuum input remains vacuum",
        "leading_tail": "g^4*t^4/4; no strict positive waiting time in the continuous model",
        "perfect_first_transfer_parameter": math.pi / math.sqrt(2),
        "threshold_arrivals": [{"threshold": threshold, "first_arrival_parameter": threshold_arrival(threshold)}
                               for threshold in (1e-8, 1e-4, 0.01, 0.5, 1)],
        "interference": [interference_report(phase) for phase in (0., math.pi / 2, math.pi)],
        "interference_formula": "central excitation, propagation pi/(2sqrt(2)), local phase at C, same propagation again: P_B=(1+cos(phi))/2",
        "split_comparison_parameter": duration, "commutator_operator_norm": commutator_norm,
        "split_role": "numerical prediction approximation; sequential physical pulses have different wall-time unless controls are rescaled",
        "split_comparison": split_rows,
        "feedback_cases": [feedback_report(DEFAULT_HOP_DELAY), feedback_report(math.pi / (3 * math.sqrt(2)))],
        "feedback_approximation": [feedback_approximation_report(DEFAULT_HOP_DELAY),
                       feedback_approximation_report(math.pi / (3 * math.sqrt(2)))],
        "unreceived_record_not_used_by_A_controller": True,
        "all_nonzero_measurement_branches_retained": True,
        "classical_latency_derived_from_quantum_Hamiltonian": False,
        "new_physical_light_cone_metric_or_gravity_derived": False,
        "full_control_detector_clock_energy_costs_closed": False,
        "measurements_or_feedback_conserve_drift_Hamiltonian_energy": "not assumed; only free drift conservation is checked",
        "samples_collected": 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PropagationTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    report = run_report()
    report["automated_checks"] = {"run": tests.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("propagation_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"model": "C1", "tests": tests.testsRun, "split_comparison": report["split_comparison"],
                      "delayed_successes": [{"calibrated": case["calibrated_average_success"], "stale": case["stale_average_success"]}
                                            for case in report["feedback_cases"]]}, indent=2))


if __name__ == "__main__":
    main()