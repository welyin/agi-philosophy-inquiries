"""Construction 02 / C2: finite relative clocks and propagation calibration.

Two independent qubit clocks extend the unchanged C1 data process. A marker
fans out to both endpoints, and measured records return to B with explicit
latency. The estimator accepts counts only; a declared phase window and the
first transfer lobe give a sufficient inversion and confidence guarantee.
No physical clock law, unique model, or necessity claim is derived here.
"""

import argparse
import math
import unittest
import json
from itertools import product
from pathlib import Path

import numpy as np

from baseline import (IDENTITY, PAULI_X, PAULI_Y, PAULI_Z, HADAMARD, CNOT,
                      binary_effect, checked_probability)
from common_orientation_structure import encode_state
from propagation import (PropagatingWhole, hamiltonian_parts, unitary_from_hamiltonian,
                         transfer_probability, real_lift)


OUTCOMES = tuple(product((0, 1), repeat=3))
DIRECTIONS = ("A_to_C", "C_to_A")
PHASE_WINDOW = (0.2, 1.3)


def clock_gate(rate, duration):
    return np.diag(np.exp(-0.5j * rate * duration * np.array([1., -1.])))


def clock_effect(basis, outcome):
    if basis not in ("X", "Y") or outcome not in (0, 1):
        raise ValueError("Clock readout requires X or Y and a binary outcome.")
    axis = PAULI_X if basis == "X" else PAULI_Y
    return (IDENTITY + (1 if outcome == 0 else -1) * axis) / 2


class ClockCalibrationTrial:
    def __init__(self, direction, coupling=1.0, hop_delay=0.4, clock_rates=(1.0, 1.2)):
        if direction not in DIRECTIONS:
            raise ValueError("Use A_to_C or C_to_A.")
        if len(clock_rates) != 2 or any(not np.isfinite(rate) or rate <= 0 for rate in clock_rates):
            raise ValueError("Two finite positive clock rates are required.")
        self.direction = direction
        self.source, self.receiver = (0, 2) if direction == "A_to_C" else (2, 0)
        self.clock_rates = tuple(float(rate) for rate in clock_rates)
        world = PropagatingWhole(coupling=coupling, hop_delay=hop_delay)
        world.apply(PAULI_X, (self.source,), "known calibration excitation")
        branches = world.record_B()
        if len(branches) != 1 or branches[0].pending_message["outcome"] != 0:
            raise ValueError("The declared blank B source must give a deterministic start marker.")
        self.data = branches[0]
        plus = (IDENTITY + PAULI_X) / 2
        self.clocks = np.kron(plus, plus)
        self.initial_joint = np.kron(self.data.density, self.clocks)
        self.initial_encoded = encode_state(self.initial_joint)
        self.complete = False
        self.readout_used = False
        self.audit_complete = False

    def arrive(self):
        if self.complete:
            raise ValueError("This start marker was already delivered.")
        duration = 2 * self.data.hop_delay
        clock_unitary = np.kron(*(clock_gate(rate, duration) for rate in self.clock_rates))
        self.clocks = clock_unitary @ self.clocks @ clock_unitary.conj().T
        self.data.advance(duration)
        self.complete = True
        self.joint = np.kron(self.data.density, self.clocks)
        data_unitary = unitary_from_hamiltonian(self.data.hamiltonian, duration)
        unitary = np.kron(data_unitary, clock_unitary)
        lifted = real_lift(unitary)
        self.encoded = lifted @ self.initial_encoded @ lifted.T
        self.data.validate()

    def projector(self, basis, outcome):
        data_projector = self.data.effect(binary_effect(outcome[2]), (self.receiver,))
        clock_projector = np.kron(clock_effect(basis, outcome[0]), clock_effect(basis, outcome[1]))
        return np.kron(data_projector, clock_projector)

    def probabilities(self, basis):
        if not self.complete:
            raise ValueError("Read the endpoint clocks only after the marker arrives.")
        probabilities = np.array([checked_probability(self.joint, self.projector(basis, outcome))
                                  for outcome in OUTCOMES])
        if abs(probabilities.sum() - 1) > 1e-12:
            raise ValueError("Incomplete joint record law.")
        return probabilities / probabilities.sum()

    def readout(self, basis):
        if self.readout_used:
            raise ValueError("The clocks are single-use; X and Y need independent fresh trials.")
        probabilities = self.probabilities(basis)
        records = []
        for outcome, probability in zip(OUTCOMES, probabilities):
            if probability == 0:
                continue
            projector = self.projector(basis, outcome)
            state = projector @ self.joint @ projector
            state /= np.trace(state).real
            records.append({"outcome": outcome, "probability": float(probability), "joint_state": state})
        self.readout_used = True
        return records

    def resource_record(self):
        return {"logical_quantum_slots": 6, "new_clock_slots": 2,
                "initial_pure_slot_preparations": 6, "clock_plus_preparations": 2,
                "start_record_slots": 1, "start_record_reads": 1,
                "endpoint_clock_reads": 2, "endpoint_data_reads": 1,
                "formal_outcome_bits": 4, "marker_bit_hops": 3,
                "audit_return_bit_hops": 3, "total_classical_bit_hops": 6,
                "classical_marker_fanouts": 1,
                "readout_parameter": 2 * self.data.hop_delay,
                "latest_audit_receipt_parameter": 4 * self.data.hop_delay,
                "active_data_edge_parameter_time_until_archive": 8 * self.data.hop_delay,
                "active_clock_parameter_time_sum_until_archive": 8 * self.data.hop_delay,
                "resets": 0, "whole_complex_dimension": 64, "whole_real_encoding_dimension": 128}

    def collect_after_audit(self, basis):
        branches = self.readout(basis)
        delay = 2 * self.data.hop_delay
        post_data_gate = unitary_from_hamiltonian(self.data.hamiltonian, delay)
        post_clock_gate = np.kron(*(clock_gate(rate, delay) for rate in self.clock_rates))
        post_gate = np.kron(post_data_gate, post_clock_gate)
        nonselective_encoded = np.zeros_like(self.encoded)
        for branch in branches:
            projector = real_lift(self.projector(basis, branch["outcome"]))
            nonselective_encoded += projector @ self.encoded @ projector.T
            branch["joint_state"] = post_gate @ branch["joint_state"] @ post_gate.conj().T
            branch["readout_at"] = 2 * self.data.hop_delay
            branch["all_records_at_B"] = 4 * self.data.hop_delay
        self.joint = sum(branch["probability"] * branch["joint_state"] for branch in branches)
        lifted = real_lift(post_gate)
        self.encoded = lifted @ nonselective_encoded @ lifted.T
        self.audit_complete = True
        self.audit_ready_parameter = 4 * self.data.hop_delay
        return branches

    def event_schedule(self):
        delay = self.data.hop_delay
        records = [{"kind": "synchronized_preparation", "time": 0.0,
                    "scope": "declared common initial phase and data launch"},
                   {"kind": "marker_at_B_and_fanout", "time": delay},
                   {"kind": "marker_at_A_and_C", "time": 2 * delay},
                   {"kind": "local_clock_and_data_readouts", "time": 2 * delay}]
        for site in ("A", "C"):
            labels = ["clock_" + site]
            if site == ("C" if self.receiver == 2 else "A"):
                labels.append("receiver_data")
            for index, label in enumerate(labels):
                records.append({"kind": "audit_record_at_B", "source": site, "label": label,
                                "sent_at": (2 + index) * delay, "time": (3 + index) * delay})
        return sorted(records, key=lambda record: record["time"])


def sample_budget(probability_error=0.01, risk=0.01):
    if not 0 < probability_error < 0.1 or not 0 < risk < 1:
        raise ValueError("Use a positive probability error below 0.1 and risk in (0,1).")
    per_group = math.ceil(math.log(20 / risk) / (2 * probability_error ** 2))
    return {"probability_error": probability_error, "risk": risk, "per_group": per_group,
            "groups": 4, "total_fresh_trials": 4 * per_group,
            "union_failure_bound": 20 * math.exp(-2 * per_group * probability_error ** 2)}


def frequencies_from_counts(groups):
    if set(groups) != {"X", "Y"}:
        raise ValueError("Both independent clock readout settings are required.")
    arrays = {basis: np.asarray(counts) for basis, counts in groups.items()}
    for counts in arrays.values():
        if counts.shape != (8,) or not np.all(np.isfinite(counts)) or np.any(counts < 0) or np.any(counts != np.floor(counts)):
            raise ValueError("Counts must be eight finite nonnegative integers.")
    totals = {basis: int(counts.sum()) for basis, counts in arrays.items()}
    if totals["X"] <= 0 or totals["X"] != totals["Y"]:
        raise ValueError("Each setting needs the same positive number of fresh trials.")
    clock_values = []
    for slot in (0, 1):
        for basis in ("X", "Y"):
            clock_values.append(float(sum(count for count, outcome in zip(arrays[basis], OUTCOMES)
                                          if outcome[slot] == 0) / totals[basis]))
    data_successes = sum(count for counts in arrays.values() for count, outcome in zip(counts, OUTCOMES)
                         if outcome[2] == 1)
    return np.array(clock_values + [float(data_successes / (2 * totals["X"]))])


def phase_interval(positive_x, positive_y, probability_error, phase_window=PHASE_WINDOW):
    lower_phase, upper_phase = phase_window
    if not 0 < lower_phase < upper_phase < math.pi / 2:
        raise ValueError("Declare a known working phase window inside the first quadrant.")
    cosine = 2 * positive_x - 1
    sine = 2 * positive_y - 1
    cosine_low = max(math.cos(upper_phase), cosine - 2 * probability_error)
    cosine_high = min(math.cos(lower_phase), cosine + 2 * probability_error)
    sine_low = max(math.sin(lower_phase), sine - 2 * probability_error)
    sine_high = min(math.sin(upper_phase), sine + 2 * probability_error)
    if (cosine_low > cosine_high or sine_low > sine_high
            or cosine_low ** 2 + sine_low ** 2 > 1 + 1e-12
            or cosine_high ** 2 + sine_high ** 2 < 1 - 1e-12):
        raise ValueError("Inconclusive calibration: clock record is incompatible with the declared window.")
    lower = max(lower_phase, math.atan2(sine_low, cosine_high))
    upper = min(upper_phase, math.atan2(sine_high, cosine_low))
    if lower > upper:
        raise ValueError("Inconclusive clock phase interval.")
    return [lower, upper]


def quotient_interval(numerator, denominator):
    if denominator[0] <= 0:
        raise ValueError("The denominator interval must stay positive.")
    return [numerator[0] / denominator[1], numerator[1] / denominator[0]]


def calibration_intervals(frequencies, probability_error=0.01, phase_window=PHASE_WINDOW):
    frequencies = np.asarray(frequencies, dtype=float)
    if (frequencies.shape != (5,) or not np.all(np.isfinite(frequencies))
            or np.any(frequencies < 0) or np.any(frequencies > 1)
            or not 0 < probability_error < 0.1):
        raise ValueError("Use five probabilities and a valid probability error.")
    phase_A = phase_interval(frequencies[0], frequencies[1], probability_error, phase_window)
    phase_C = phase_interval(frequencies[2], frequencies[3], probability_error, phase_window)
    transfer = [max(0., frequencies[4] - probability_error), min(1., frequencies[4] + probability_error)]
    coupling_times_duration = [float(np.sqrt(2) * np.arcsin(value ** 0.25)) for value in transfer]
    return {"phase_A": phase_A, "phase_C": phase_C, "transfer_probability": transfer,
            "g_times_arrival_duration": coupling_times_duration,
            "g_over_omega_A": quotient_interval(coupling_times_duration, phase_A),
            "g_over_omega_C": quotient_interval(coupling_times_duration, phase_C),
            "omega_C_over_omega_A": quotient_interval(phase_C, phase_A),
            "omega_A_times_hop_delay": [value / 2 for value in phase_A]}


def intervals_overlap(first, second):
    return max(first[0], second[0]) <= min(first[1], second[1])


def feedback_certificate(transfer_interval):
    if len(transfer_interval) != 2 or not 0 <= transfer_interval[0] <= transfer_interval[1] <= 1:
        raise ValueError("Use a valid transfer probability interval.")
    if transfer_interval[1] < 0.5:
        return "flip_after_record_one"
    if transfer_interval[0] > 0.5:
        return "keep_after_record_one"
    return "unresolved_at_this_budget"


def certified_feedback_trial(transfer_interval, coupling=1.0, hop_delay=0.4):
    certificate = feedback_certificate(transfer_interval)
    if certificate == "unresolved_at_this_budget":
        return {"status": certificate, "task_executed": False}
    world = PropagatingWhole(coupling=coupling, hop_delay=hop_delay)
    world.apply(HADAMARD, (0,), "fresh feedback A plus")
    world.apply(CNOT, (0, 1), "fresh feedback Bell source")
    records = []
    for branch in world.record_B():
        branch.advance(2 * hop_delay)
        if not branch.delivered_records:
            raise ValueError("The new task record has not arrived.")
        outcome = branch.delivered_records[-1]["outcome"]
        action = "flip" if outcome == 1 and certificate == "flip_after_record_one" else "keep"
        branch.history.append({"kind": "calibrated_interval_decision", "time": branch.time,
                               "action": action, "observed_record": outcome,
                               "transfer_interval": list(transfer_interval)})
        branch.feedback_done = True
        branch.resources["feedback_message_bits"] += 1
        if action == "flip":
            branch.apply(PAULI_X, (0,), "count-certified feedback action")
        branch.validate()
        records.append({"outcome": outcome, "weight": branch.weight, "action": action,
                        "decision_parameter": branch.time, "success_probability": branch.probability(0)})
    return {"status": certificate, "task_executed": True, "branches": records,
            "average_success_probability": sum(record["weight"] * record["success_probability"] for record in records),
            "policy_inputs": "delivered two-bit validity/action code and delivered task record only; no continuous parameter estimate transmitted",
            "stable_coupling_and_delay_between_calibration_and_task_assumed": True,
            "additional_task_resources": {"fresh_logical_slots": 4, "M_readouts": 1,
                          "task_record_bit_hops": 2, "calibration_action_code_bits": 2,
                          "action_code_bit_hops_B_to_A": 2, "formal_record_bits": 1,
                                          "feedback_message_bits": 1, "resets": 0}}


def simulate_calibration(coupling=1.0, hop_delay=0.4, clock_rates=(1.0, 1.2), seed=20260917):
    if (not np.isfinite(coupling) or not np.isfinite(hop_delay) or coupling <= 0 or hop_delay <= 0
            or len(clock_rates) != 2 or any(not np.isfinite(rate) or rate <= 0 for rate in clock_rates)):
        raise ValueError("The simulation requires finite positive rates and delay.")
    phases = [2 * hop_delay * rate for rate in clock_rates]
    if any(not PHASE_WINDOW[0] <= phase <= PHASE_WINDOW[1] for phase in phases):
        raise ValueError("Simulation outside the declared clock phase window; no calibration guarantee claimed.")
    if not 0 < 2 * coupling * hop_delay <= math.pi / math.sqrt(2):
        raise ValueError("Simulation outside the declared first transfer lobe; no calibration guarantee claimed.")
    budget = sample_budget()
    random = np.random.default_rng(seed)
    counts = {}
    probability_tables = {}
    summaries = {}
    event_examples = {}
    resource = None
    for direction in DIRECTIONS:
        counts[direction] = {}
        probability_tables[direction] = {}
        for basis in ("X", "Y"):
            trial = ClockCalibrationTrial(direction, coupling, hop_delay, clock_rates)
            trial.arrive()
            records = trial.collect_after_audit(basis)
            if not trial.audit_complete:
                raise ValueError("Cannot estimate from records that have not returned to B.")
            by_outcome = {record["outcome"]: record["probability"] for record in records}
            probabilities = np.array([by_outcome.get(outcome, 0.0) for outcome in OUTCOMES])
            counts[direction][basis] = random.multinomial(budget["per_group"], probabilities).tolist()
            probability_tables[direction][basis] = probabilities.tolist()
            event_examples[direction] = trial.event_schedule()
            resource = trial.resource_record()
        frequencies = frequencies_from_counts(counts[direction])
        intervals = calibration_intervals(frequencies, budget["probability_error"])
        summaries[direction] = {"frequencies": frequencies.tolist(), "intervals": intervals,
                                "feedback_action_certificate": feedback_certificate(intervals["transfer_probability"])}
    compared = ("g_over_omega_A", "omega_C_over_omega_A", "omega_A_times_hop_delay")
    compatibility = {name: intervals_overlap(summaries[DIRECTIONS[0]]["intervals"][name],
                                              summaries[DIRECTIONS[1]]["intervals"][name]) for name in compared}
    trials = budget["total_fresh_trials"]
    combined_transfer = [max(summaries[direction]["intervals"]["transfer_probability"][0] for direction in DIRECTIONS),
                         min(summaries[direction]["intervals"]["transfer_probability"][1] for direction in DIRECTIONS)]
    feedback = (certified_feedback_trial(combined_transfer, coupling, hop_delay)
                if combined_transfer[0] <= combined_transfer[1] else
                {"status": "incompatible_transfer_intervals", "task_executed": False})
    calibration_duration = trials * 4 * hop_delay
    if feedback["task_executed"]:
        feedback["action_code_sent_at_B"] = calibration_duration
        feedback["action_code_received_at_A"] = calibration_duration + 2 * hop_delay
        feedback["new_task_epoch"] = feedback["action_code_received_at_A"]
        for branch in feedback["branches"]:
            branch["global_decision_parameter"] = feedback["new_task_epoch"] + branch["decision_parameter"]
    return {"model_version": "C2", "route_round": "02", "scope": "sufficiency only",
            "simulation_parameters_not_given_to_estimator": {"coupling": coupling, "hop_delay": hop_delay,
                                                              "clock_rates": list(clock_rates)},
            "declared_inversion_domain": {"clock_phase_window": list(PHASE_WINDOW),
                                           "transfer_phase_g_times_duration": [0, math.pi / math.sqrt(2)]},
            "sampling": budget, "simulated_count_seed": seed, "outcome_order": [list(outcome) for outcome in OUTCOMES],
            "simulated_counts": counts, "Born_probability_tables": probability_tables,
            "count_only_estimates": summaries, "directional_interval_compatibility": compatibility,
            "resource_per_trial": resource,
            "batch_resource_counts": {"fresh_trial_kits": trials, "pure_slot_initializations": 6 * trials,
                                      "clock_reads": 2 * trials, "formal_outcome_bits_if_fully_archived": 4 * trials,
                                      "classical_bit_hops": 6 * trials,
                                      "histogram_storage_bits_sufficient_excluding_metadata": 32 * math.ceil(math.log2(budget["per_group"] + 1)),
                                      "serial_parameter_duration_excluding_preparation": calibration_duration},
            "timed_event_examples": event_examples,
            "subsequent_feedback": feedback,
            "calibration_finished_before_the_fresh_feedback_source_is_prepared": True,
            "clock_preparation_synchrony_rate_stability_and_ideal_readouts_are_inputs": True,
            "clocks_generate_the_control_schedule_autonomously": False,
            "all_nonzero_outcomes_in_probability_and_count_tables": True,
            "saved_histograms_are_not_full_per_trial_raw_archives": True,
            "experimental_samples_collected": 0,
            "simulated_fresh_trials": trials,
            "physical_SI_clock_metric_or_gravity_derived": False,
            "necessity_arguments_performed": False}


class ClockCalibrationTests(unittest.TestCase):
    def test_six_slot_state_is_the_same_C1_process_with_two_declared_clocks(self):
        trial = ClockCalibrationTrial("A_to_C")
        trial.arrive()
        self.assertEqual(trial.joint.shape, (64, 64))
        np.testing.assert_allclose(trial.encoded, encode_state(trial.joint), rtol=0, atol=2e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(trial.joint).min(), -1e-12)
        self.assertAlmostEqual(np.trace(trial.joint).real, 1)

    def test_factorized_evolution_matches_one_total_hamiltonian(self):
        trial = ClockCalibrationTrial("C_to_A")
        clock_hamiltonian = (np.kron(0.5 * trial.clock_rates[0] * PAULI_Z, IDENTITY)
                             + np.kron(IDENTITY, 0.5 * trial.clock_rates[1] * PAULI_Z))
        total = np.kron(sum(hamiltonian_parts()), np.eye(4)) + np.kron(np.eye(16), clock_hamiltonian)
        unitary = unitary_from_hamiltonian(total, 0.8)
        expected = unitary @ trial.initial_joint @ unitary.conj().T
        trial.arrive()
        np.testing.assert_allclose(trial.joint, expected, rtol=0, atol=2e-12)

    def test_joint_readout_table_matches_clock_phases_and_transfer_analytically(self):
        for direction in DIRECTIONS:
            for basis in ("X", "Y"):
                trial = ClockCalibrationTrial(direction)
                trial.arrive()
                function = math.cos if basis == "X" else math.sin
                positives = [(1 + function(2 * trial.data.hop_delay * rate)) / 2 for rate in trial.clock_rates]
                transfer = transfer_probability(2 * trial.data.hop_delay)
                expected = [np.prod([positives[0] if outcome[0] == 0 else 1 - positives[0],
                                     positives[1] if outcome[1] == 0 else 1 - positives[1],
                                     transfer if outcome[2] == 1 else 1 - transfer]) for outcome in OUTCOMES]
                np.testing.assert_allclose(trial.probabilities(basis), expected, rtol=0, atol=2e-12)

    def test_measurement_branches_are_normalized_and_clocks_not_reused(self):
        trial = ClockCalibrationTrial("A_to_C")
        with self.assertRaises(ValueError):
            trial.probabilities("X")
        trial.arrive()
        branches = trial.readout("Y")
        self.assertEqual(len(branches), 8)
        self.assertAlmostEqual(sum(branch["probability"] for branch in branches), 1)
        for branch in branches:
            self.assertAlmostEqual(np.trace(branch["joint_state"]).real, 1)
            self.assertGreaterEqual(np.linalg.eigvalsh(branch["joint_state"]).min(), -1e-12)
        with self.assertRaises(ValueError):
            trial.readout("X")

    def test_marker_and_audit_records_have_explicit_nonzero_delivery_times(self):
        for direction in DIRECTIONS:
            trial = ClockCalibrationTrial(direction)
            trial.arrive()
            events = trial.event_schedule()
            receipts = [event for event in events if event["kind"] == "audit_record_at_B"]
            self.assertEqual(len(receipts), 3)
            self.assertEqual(max(event["time"] for event in receipts), 1.6)
            self.assertTrue(all(event["time"] > event["sent_at"] >= 0.8 for event in receipts))
            resources = trial.resource_record()
            self.assertEqual(resources["total_classical_bit_hops"], 6)
            self.assertEqual(resources["logical_quantum_slots"], 6)
            self.assertEqual(resources["resets"], 0)

    def test_postmeasurement_quantum_branches_are_retained_during_audit_transport(self):
        trial = ClockCalibrationTrial("A_to_C")
        trial.arrive()
        records = trial.collect_after_audit("X")
        self.assertTrue(trial.audit_complete)
        self.assertEqual(trial.audit_ready_parameter, 1.6)
        np.testing.assert_allclose(trial.encoded, encode_state(trial.joint), rtol=0, atol=2e-12)
        for record in records:
            self.assertEqual(record["readout_at"], 0.8)
            self.assertEqual(record["all_records_at_B"], 1.6)
            self.assertAlmostEqual(np.trace(record["joint_state"]).real, 1)
            self.assertGreaterEqual(np.linalg.eigvalsh(record["joint_state"]).min(), -1e-12)

    def test_exact_probabilities_recover_declared_dimensionless_quantities(self):
        frequencies = [(1 + math.cos(0.8)) / 2, (1 + math.sin(0.8)) / 2,
                       (1 + math.cos(0.96)) / 2, (1 + math.sin(0.96)) / 2, transfer_probability(0.8)]
        intervals = calibration_intervals(frequencies, 0.001)
        expected = {"g_over_omega_A": 1., "g_over_omega_C": 1 / 1.2,
                    "omega_C_over_omega_A": 1.2, "omega_A_times_hop_delay": 0.4}
        for name, value in expected.items():
            self.assertLess(intervals[name][0], value)
            self.assertGreater(intervals[name][1], value)

    def test_component_error_corners_are_covered_by_nonlinear_interval_propagation(self):
        exact = np.array([(1 + math.cos(0.8)) / 2, (1 + math.sin(0.8)) / 2,
                          (1 + math.cos(0.96)) / 2, (1 + math.sin(0.96)) / 2, transfer_probability(0.8)])
        for signs in product((-1, 1), repeat=5):
            intervals = calibration_intervals(exact + 0.01 * np.array(signs))
            for name, value in (("g_over_omega_A", 1), ("omega_C_over_omega_A", 1.2), ("omega_A_times_hop_delay", 0.4)):
                self.assertLessEqual(intervals[name][0], value + 1e-12)
                self.assertGreaterEqual(intervals[name][1], value - 1e-12)

    def test_simulated_counts_retain_all_trials_and_give_compatible_endpoint_estimates(self):
        report = simulate_calibration()
        self.assertTrue(all(report["directional_interval_compatibility"].values()))
        for direction in DIRECTIONS:
            for counts in report["simulated_counts"][direction].values():
                self.assertEqual(sum(counts), report["sampling"]["per_group"])
            estimates = report["count_only_estimates"][direction]
            self.assertEqual(estimates["feedback_action_certificate"], "flip_after_record_one")
            for name, value in (("g_over_omega_A", 1), ("omega_C_over_omega_A", 1.2), ("omega_A_times_hop_delay", 0.4)):
                self.assertLessEqual(estimates["intervals"][name][0], value)
                self.assertGreaterEqual(estimates["intervals"][name][1], value)

    def test_union_bound_covers_ten_mean_estimates_without_independence_between_means(self):
        budget = sample_budget()
        self.assertLessEqual(budget["union_failure_bound"], 0.01)
        self.assertEqual(budget["total_fresh_trials"], 4 * budget["per_group"])
        tighter = (16 * math.exp(-2 * budget["per_group"] * 0.01 ** 2)
                   + 4 * math.exp(-4 * budget["per_group"] * 0.01 ** 2))
        self.assertLessEqual(tighter, budget["union_failure_bound"])

    def test_count_interval_drives_a_fresh_task_without_zero_delay_record_access(self):
        report = simulate_calibration()
        task = report["subsequent_feedback"]
        self.assertTrue(task["task_executed"])
        self.assertAlmostEqual(task["average_success_probability"], 1 - transfer_probability(0.8) / 2)
        self.assertEqual([branch["action"] for branch in task["branches"]], ["keep", "flip"])
        self.assertTrue(all(branch["decision_parameter"] >= 0.8 for branch in task["branches"]))
        self.assertGreater(task["action_code_received_at_A"], task["action_code_sent_at_B"])
        self.assertEqual(task["new_task_epoch"], task["action_code_received_at_A"])
        self.assertTrue(all(branch["global_decision_parameter"] > task["new_task_epoch"] for branch in task["branches"]))
        self.assertFalse(certified_feedback_trial([0.49, 0.51])["task_executed"])

    def test_unsupported_or_inconsistent_inputs_are_reported_not_silently_fitted(self):
        for groups in ({"X": [1] * 8}, {"X": [-1] * 8, "Y": [1] * 8},
                       {"X": [0] * 8, "Y": [0] * 8}, {"X": [1] * 8, "Y": [2] * 8}):
            with self.assertRaises(ValueError):
                frequencies_from_counts(groups)
        with self.assertRaises(ValueError):
            phase_interval(0.01, 0.01, 0.01)
        with self.assertRaises(ValueError):
            ClockCalibrationTrial("unsupported")
        self.assertEqual(feedback_certificate([0.49, 0.51]), "unresolved_at_this_budget")

    def test_simulation_rejects_examples_outside_the_predeclared_operating_contract(self):
        with self.assertRaises(ValueError):
            simulate_calibration(clock_rates=(10., 10.))
        with self.assertRaises(ValueError):
            simulate_calibration(coupling=10.)
        with self.assertRaises(ValueError):
            simulate_calibration(hop_delay=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ClockCalibrationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = simulate_calibration()
    report["automated_checks"] = {"run": checks.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("clock_calibration_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"C2_tests": checks.testsRun, "sampling": report["sampling"],
                      "estimates": report["count_only_estimates"],
                      "compatibility": report["directional_interval_compatibility"]}, indent=2))


if __name__ == "__main__":
    main()