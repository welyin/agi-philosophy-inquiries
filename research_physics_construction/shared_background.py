"""Construction 03 / C3: a shared static clock-and-link background.

Two declared probe types use different calibrated clock frequencies but the
same spatial coupling rule. Count-only inference reconstructs normalized
clock and bond parameters; endpoint transfer is a held-out consistency test.
The metric interpretation is a specified local-dispersion convention, not
a continuum or Einstein-equation derivation. Necessity is not investigated.
"""

import argparse
import json
import math
import unittest
from itertools import product
from pathlib import Path

import numpy as np

from baseline import binary_effect, checked_probability, HADAMARD, CNOT, PAULI_X
from clock_calibration import (ClockCalibrationTrial, clock_effect, phase_interval,
                               quotient_interval, feedback_certificate, PHASE_WINDOW, DIRECTIONS)
from common_orientation_structure import encode_state
from propagation import hamiltonian_parts, PropagatingWhole


PROBES = {"P": 1.0, "Q": 1.2}
RECORDS = tuple(product((0, 1), repeat=4))


def declared_background():
    positions = np.array([0., 1., 2.])
    lapse = 1 + 0.15 * positions
    midpoints = (positions[:-1] + positions[1:]) / 2
    mid_lapse = 1 + 0.15 * midpoints
    spatial = 1 + 0.2 * midpoints
    bonds = mid_lapse / (2 * spatial)
    return {"positions": positions.tolist(), "lapse": lapse.tolist(),
            "spatial_midpoints": spatial.tolist(), "bonds": bonds.tolist(),
            "hop_delay": 0.4, "coordinate_spacing": 1.0}


def asymmetric_probabilities(bond_times):
    left, right = np.asarray(bond_times, dtype=float)
    frequency = math.hypot(left, right)
    if frequency == 0:
        return np.zeros(3)
    ratio = left ** 2 / frequency ** 2
    middle_left = ratio * math.sin(frequency) ** 2
    middle_right = (1 - ratio) * math.sin(frequency) ** 2
    endpoint = ratio * (1 - ratio) * (1 - math.cos(frequency)) ** 2
    return np.array([middle_left, middle_right, endpoint])


def annotate_bond_history(world, bonds):
    for event in world.history:
        if event["kind"] == "propagation":
            event.pop("coupling", None)
            event["edge_couplings"] = {"A_B": float(bonds[0]), "B_C": float(bonds[1])}


class BackgroundTrial(ClockCalibrationTrial):
    def __init__(self, direction, probe, bonds, endpoint_lapses, hop_delay=0.4):
        if probe not in PROBES:
            raise ValueError("Unknown calibrated probe type.")
        if len(bonds) != 2 or len(endpoint_lapses) != 2:
            raise ValueError("Two bonds and two endpoint lapse factors are required.")
        values = np.array([*bonds, *endpoint_lapses, hop_delay], dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values <= 0):
            raise ValueError("Background parameters must be finite and positive.")
        rates = tuple(PROBES[probe] * value for value in endpoint_lapses)
        duration = 2 * hop_delay
        if any(not PHASE_WINDOW[0] <= rate * duration <= PHASE_WINDOW[1] for rate in rates):
            raise ValueError("Outside the predeclared clock phase window.")
        if not 0 < math.hypot(*bonds) * duration < math.pi / 2:
            raise ValueError("Use the first monotonic middle-population window.")
        super().__init__(direction, coupling=1., hop_delay=hop_delay, clock_rates=rates)
        first, second = hamiltonian_parts()
        self.data.hamiltonian = bonds[0] * first + bonds[1] * second
        self.probe = probe
        self.bonds = tuple(bonds)

    def arrive(self):
        super().arrive()
        annotate_bond_history(self.data, self.bonds)

    def projector(self, basis, outcome):
        middle = self.data.effect(binary_effect(outcome[2]), (1,))
        receiver = self.data.effect(binary_effect(outcome[3]), (self.receiver,))
        clocks = np.kron(clock_effect(basis, outcome[0]), clock_effect(basis, outcome[1]))
        return np.kron(middle @ receiver, clocks)

    def probabilities(self, basis):
        if not self.complete:
            raise ValueError("Read after the marker and local trigger are ready.")
        values = np.array([checked_probability(self.joint, self.projector(basis, outcome))
                           for outcome in RECORDS])
        if abs(values.sum() - 1) > 1e-12:
            raise ValueError("The joint record law is not normalized.")
        return values / values.sum()

    def readout(self, basis):
        if self.readout_used:
            raise ValueError("Fresh preparations are required for different clock bases.")
        records = []
        for outcome, probability in zip(RECORDS, self.probabilities(basis)):
            if probability == 0:
                continue
            projector = self.projector(basis, outcome)
            state = projector @ self.joint @ projector
            state /= np.trace(state).real
            records.append({"outcome": outcome, "probability": float(probability), "joint_state": state})
        self.readout_used = True
        return records

    def resource_record(self):
        resources = super().resource_record()
        resources.update({"central_data_reads": 1, "formal_outcome_bits": 5,
                          "local_delayed_trigger_at_B": 1,
                          "B_trigger_delay_is_calibrated_controller_input": True})
        return resources

    def event_schedule(self):
        records = super().event_schedule()
        records.append({"kind": "local_B_data_readout", "time": 2 * self.data.hop_delay,
                        "scope": "one-hop marker at B followed by declared local delay tau"})
        return sorted(records, key=lambda event: event["time"])


def means_from_counts(groups):
    if set(groups) != {"X", "Y"}:
        raise ValueError("Independent X/Y settings required.")
    arrays = {basis: np.asarray(counts) for basis, counts in groups.items()}
    for counts in arrays.values():
        if (counts.shape != (16,) or not np.all(np.isfinite(counts)) or np.any(counts < 0)
                or np.any(counts != np.floor(counts))):
            raise ValueError("Use sixteen nonnegative integer counts.")
    total = int(arrays["X"].sum())
    if total <= 0 or int(arrays["Y"].sum()) != total:
        raise ValueError("Both settings need equal positive sample counts.")
    means = []
    for slot in (0, 1):
        for basis in ("X", "Y"):
            means.append(sum(count for count, outcome in zip(arrays[basis], RECORDS) if outcome[slot] == 0) / total)
    for slot in (2, 3):
        means.append(sum(count for counts in arrays.values() for count, outcome in zip(counts, RECORDS)
                         if outcome[slot] == 1) / (2 * total))
    return np.asarray(means, dtype=float)


def inverse_bonds(middle_left, middle_right):
    total = middle_left + middle_right
    if not 0 < total < 1 or min(middle_left, middle_right) < 0:
        raise ValueError("Inconclusive middle-population inversion.")
    frequency = math.asin(math.sqrt(total))
    return frequency * np.sqrt(np.array([middle_left, middle_right]) / total)


def geometry_from_phases(phase_A, phase_C, bond_times, spacing=1.):
    if phase_A <= 0 or phase_C <= 0 or min(bond_times) <= 0 or spacing <= 0:
        raise ValueError("Positive scales are required.")
    ratio = phase_C / phase_A
    normalized_bonds = np.asarray(bond_times) / phase_A
    mid_lapse = np.array([0.75 + 0.25 * ratio, 0.25 + 0.75 * ratio])
    return {"lapse_C_over_A": float(ratio),
            "g_edges_over_lapse_A": normalized_bonds.tolist(),
            "normalized_lapse_midpoints": mid_lapse.tolist(),
            "spatial_midpoints": (mid_lapse / (2 * spacing * normalized_bonds)).tolist(),
            "candidate_metric_null_speeds_over_lapse_A": (2 * spacing * normalized_bonds).tolist(),
            "interpretation": "piecewise static metric convention with linear lapse interpolation, not derived continuum dynamics"}


def estimate_probe(left, right, bare_frequency, error=0.01):
    left, right = np.asarray(left), np.asarray(right)
    if (left.shape != (6,) or right.shape != (6,) or not np.isfinite(bare_frequency)
            or bare_frequency <= 0 or not 0 < error < 0.1
            or not np.all(np.isfinite(np.r_[left, right])) or np.any(np.r_[left, right] < 0)
            or np.any(np.r_[left, right] > 1)):
        raise ValueError("Use two six-mean records and declared positive calibration constants.")
    clocks = (left[:4] + right[:4]) / 2
    phase_A = np.array(phase_interval(clocks[0], clocks[1], error)) / bare_frequency
    phase_C = np.array(phase_interval(clocks[2], clocks[3], error)) / bare_frequency
    probabilities = [[max(0., mean - error), min(1., mean + error)] for mean in (left[4], right[4])]
    summed = np.sum(probabilities, axis=0)
    if not 0 < summed[0] <= summed[1] < 1:
        raise ValueError("Inconclusive calibration in the specified monotonic window.")
    frequency = np.arcsin(np.sqrt(summed))
    bond_intervals = [[float(frequency[0] * math.sqrt(probability[0] / summed[1])),
                       float(frequency[1] * math.sqrt(min(1., probability[1] / summed[0])))]
                      for probability in probabilities]
    lapse_ratio = quotient_interval(phase_C, phase_A)
    normalized = [quotient_interval(interval, phase_A) for interval in bond_intervals]
    lapse_mid = [[0.75 + 0.25 * lapse_ratio[0], 0.75 + 0.25 * lapse_ratio[1]],
                 [0.25 + 0.75 * lapse_ratio[0], 0.25 + 0.75 * lapse_ratio[1]]]
    spatial = [quotient_interval(lapse, [2 * value for value in bond]) for lapse, bond in zip(lapse_mid, normalized)]
    fraction = [probabilities[0][0] / summed[1], min(1., probabilities[0][1] / summed[0])]
    product_range = [min(value * (1 - value) for value in fraction),
                     0.25 if fraction[0] <= 0.5 <= fraction[1] else max(value * (1 - value) for value in fraction)]
    predicted_endpoint = [float(product_range[0] * (1 - math.cos(frequency[0])) ** 2),
                          float(product_range[1] * (1 - math.cos(frequency[1])) ** 2)]
    return {"lapse_A_times_T": phase_A.tolist(), "lapse_C_times_T": phase_C.tolist(),
            "bond_times_T": bond_intervals, "lapse_C_over_A": lapse_ratio,
            "g_edges_over_lapse_A": normalized, "spatial_midpoints": spatial,
            "predicted_endpoint_probability_not_used_in_fit": predicted_endpoint}


def metric_interpolation(x, lapse_ratio, spatial_midpoints):
    position = np.asarray(x, dtype=float)
    if (not np.all(np.isfinite(position)) or np.any(position < 0) or np.any(position > 2)
            or not np.isfinite(lapse_ratio) or lapse_ratio <= 0 or len(spatial_midpoints) != 2
            or not np.all(np.isfinite(spatial_midpoints)) or min(spatial_midpoints) <= 0):
        raise ValueError("Use positive sampled scales and coordinates in the declared patch [0,2].")
    lapse = 1 + (lapse_ratio - 1) * position / 2
    spatial = np.interp(position, [0.5, 1.5], spatial_midpoints)
    return lapse, spatial


def shared_count_witness(means):
    phase_A, phase_C, middle_left, middle_right = [], [], [], []
    for probe, bare in PROBES.items():
        for direction in DIRECTIONS:
            values = np.asarray(means[probe][direction])
            phase_A.append(math.atan2(2 * values[1] - 1, 2 * values[0] - 1) / bare)
            phase_C.append(math.atan2(2 * values[3] - 1, 2 * values[2] - 1) / bare)
        middle_left.append(means[probe]["A_to_C"][4])
        middle_right.append(means[probe]["C_to_A"][4])
    phase_A, phase_C = float(np.mean(phase_A)), float(np.mean(phase_C))
    bonds = inverse_bonds(float(np.mean(middle_left)), float(np.mean(middle_right)))
    probabilities = asymmetric_probabilities(bonds)
    residuals = []
    endpoint_residuals = []
    for probe, bare in PROBES.items():
        clock_means = [(1 + math.cos(bare * phase_A)) / 2, (1 + math.sin(bare * phase_A)) / 2,
                       (1 + math.cos(bare * phase_C)) / 2, (1 + math.sin(bare * phase_C)) / 2]
        for index, direction in enumerate(DIRECTIONS):
            prediction = np.array(clock_means + [probabilities[index], probabilities[2]])
            residuals.extend(np.abs(prediction - means[probe][direction]))
            endpoint_residuals.append(abs(prediction[-1] - means[probe][direction][-1]))
    return {"lapse_A_times_T": phase_A, "lapse_C_times_T": phase_C,
            "bond_times_T": bonds.tolist(), "geometry": geometry_from_phases(phase_A, phase_C, bonds),
            "predicted_endpoint_probability": float(probabilities[2]),
            "maximum_probability_residual": float(max(residuals)),
            "maximum_heldout_endpoint_residual": float(max(endpoint_residuals)),
            "fitting_inputs": "clock and middle-population counts only; endpoint counts excluded"}


def calibrated_background_feedback(bonds, hop_delay, endpoint_interval):
    certificate = feedback_certificate(endpoint_interval)
    if certificate == "unresolved_at_this_budget":
        return {"executed": False, "status": certificate}
    world = PropagatingWhole(hop_delay=hop_delay)
    first, second = hamiltonian_parts()
    world.hamiltonian = bonds[0] * first + bonds[1] * second
    world.apply(HADAMARD, (0,), "fresh Bell source A plus")
    world.apply(CNOT, (0, 1), "fresh AB Bell preparation")
    branches = []
    for branch in world.record_B():
        branch.advance(2 * hop_delay)
        annotate_bond_history(branch, bonds)
        if not branch.delivered_records:
            raise ValueError("Feedback requires the delivered new task record.")
        outcome = branch.delivered_records[-1]["outcome"]
        action = "flip" if outcome == 1 and certificate == "flip_after_record_one" else "keep"
        branch.resources["feedback_message_bits"] += 1
        branch.feedback_done = True
        branch.history.append({"kind": "interval_decision", "time": branch.time,
                               "record": outcome, "action": action, "calibration_code": certificate})
        if action == "flip":
            branch.apply(PAULI_X, (0,), "calibration-certified A flip")
        branch.validate()
        branches.append({"record": outcome, "action": action, "weight": branch.weight,
                         "decision_parameter_from_task_start": branch.time,
                         "success": branch.probability(0), "resources": branch.resources})
    predicted = ([1 - endpoint_interval[1] / 2, 1 - endpoint_interval[0] / 2]
                 if certificate == "flip_after_record_one" else
                 [(1 + endpoint_interval[0]) / 2, (1 + endpoint_interval[1]) / 2])
    return {"executed": True, "status": certificate, "branches": branches,
            "predicted_success_interval": predicted,
            "simulation_success": sum(branch["weight"] * branch["success"] for branch in branches),
            "decision_inputs": "two-bit validity/action code and delivered new task record only",
            "new_task_extra_resources": {"fresh_quantum_slots": 4, "task_record_reads": 1,
                                         "calibration_code_bit_hops": 2, "task_record_bit_hops": 2,
                                         "resets": 0}}


def simulate_shared_background(seed=20260918, error=0.01, risk=0.01):
    if not 0 < error < 0.1 or not 0 < risk < 1:
        raise ValueError("Use valid fixed precision and risk budgets.")
    parameters = declared_background()
    count_per_group = math.ceil(math.log(48 / risk) / (2 * error ** 2))
    random = np.random.default_rng(seed)
    counts, means, estimates = {}, {}, {}
    resources, schedules = None, {}
    for probe, bare in PROBES.items():
        counts[probe], means[probe] = {}, {}
        for direction in DIRECTIONS:
            counts[probe][direction] = {}
            for basis in ("X", "Y"):
                trial = BackgroundTrial(direction, probe, parameters["bonds"],
                                        [parameters["lapse"][0], parameters["lapse"][-1]], parameters["hop_delay"])
                trial.arrive()
                branches = trial.collect_after_audit(basis)
                if not trial.audit_complete:
                    raise ValueError("Counts unavailable before audit delivery.")
                by_outcome = {branch["outcome"]: branch["probability"] for branch in branches}
                law = np.array([by_outcome.get(outcome, 0.) for outcome in RECORDS])
                counts[probe][direction][basis] = random.multinomial(count_per_group, law / law.sum()).tolist()
                resources = trial.resource_record()
                schedules[direction] = trial.event_schedule()
            means[probe][direction] = means_from_counts(counts[probe][direction]).tolist()
        estimates[probe] = estimate_probe(means[probe][DIRECTIONS[0]], means[probe][DIRECTIONS[1]], bare, error)
    witness = shared_count_witness(means)
    witness["all_means_within_declared_error"] = witness["maximum_probability_residual"] <= error
    total = 8 * count_per_group
    endpoint_interval = [max(estimates[probe]["predicted_endpoint_probability_not_used_in_fit"][0] for probe in PROBES),
                         min(estimates[probe]["predicted_endpoint_probability_not_used_in_fit"][1] for probe in PROBES)]
    feedback = (calibrated_background_feedback(parameters["bonds"], parameters["hop_delay"], endpoint_interval)
                if endpoint_interval[0] <= endpoint_interval[1] else
                {"executed": False, "status": "incompatible_calibration_intervals"})
    if feedback["executed"]:
        feedback["calibration_archive_ready_at_B"] = 4 * parameters["hop_delay"] * total
        feedback["two_bit_code_ready_at_A_and_new_task_start"] = (
            feedback["calibration_archive_ready_at_B"] + 2 * parameters["hop_delay"])
        feedback["final_decision_parameter"] = feedback["two_bit_code_ready_at_A_and_new_task_start"] + 2 * parameters["hop_delay"]
    return {"model_version": "C3", "route_round": "03", "research_scope": "sufficiency only",
            "declared_background_not_given_to_estimator": parameters,
            "known_bare_clock_frequencies": PROBES,
            "shared_hopping_and_clock_coupling_are_model_inputs": True,
            "probe_types_differ_in_clock_transition_not_proved_relativistic_mass": True,
            "outcome_order": [list(outcome) for outcome in RECORDS],
            "endpoint_holdout_is_an_observable_not_an_independent_trial_split": True,
            "sample_budget": {"error_per_mean": error, "risk": risk, "means": 24,
                              "groups": 8, "trials_per_group": count_per_group, "total_fresh_trials": total,
                              "union_failure_bound": 48 * math.exp(-2 * count_per_group * error ** 2)},
            "seed": seed, "simulated_counts": counts, "count_means": means,
            "per_probe_intervals": estimates, "shared_parameter_witness": witness,
            "calibrated_fresh_feedback": feedback,
            "resource_per_trial": resources, "event_schedules": schedules,
            "batch_resources": {"fresh_trial_kits": total, "pure_slot_initializations": 6 * total,
                                "formal_outcome_bits_if_archived": 5 * total, "classical_bit_hops": 6 * total,
                                "histogram_storage_bits_excluding_metadata": 128 * math.ceil(math.log2(count_per_group + 1)),
                                "serial_parameter_duration_excluding_preparation": 4 * parameters["hop_delay"] * total},
            "same_spatial_samples_used_by_both_probes": True,
            "metric_convention": "ell_A=1; linear lapse interpolation; b_edge=ell_mid/(2*a*g_edge); spatial b linearly interpolated between bond centers with constant endpoint extension",
            "continuum_Dirac_dynamics_or_equivalence_principle_derived": False,
            "background_backreaction_Einstein_equations_or_metric_uniqueness_derived": False,
            "simulated_trials_not_laboratory_samples": total, "experimental_samples": 0,
            "necessity_arguments_performed": False}


class SharedBackgroundTests(unittest.TestCase):
    def test_candidate_metric_is_positive_on_patch_and_reproduces_both_bond_conventions(self):
        report = shared_count_witness({probe: {
            direction: [(1 + math.cos(bare * 0.8)) / 2, (1 + math.sin(bare * 0.8)) / 2,
                        (1 + math.cos(bare * 1.04)) / 2, (1 + math.sin(bare * 1.04)) / 2,
                        asymmetric_probabilities(0.8 * np.array(declared_background()["bonds"]))[index],
                        asymmetric_probabilities(0.8 * np.array(declared_background()["bonds"]))[2]]
            for index, direction in enumerate(DIRECTIONS)} for probe, bare in PROBES.items()})
        geometry = report["geometry"]
        lapse, spatial = metric_interpolation(np.linspace(0, 2, 21), geometry["lapse_C_over_A"], geometry["spatial_midpoints"])
        self.assertTrue(np.all(lapse > 0) and np.all(spatial > 0))
        midpoint_lapse, midpoint_spatial = metric_interpolation(np.array([0.5, 1.5]), geometry["lapse_C_over_A"], geometry["spatial_midpoints"])
        np.testing.assert_allclose(midpoint_lapse / midpoint_spatial,
                                   2 * np.array(geometry["g_edges_over_lapse_A"]), atol=1e-14)

    def test_speed_convention_matches_the_frozen_uniform_chain_dispersion(self):
        spacing = 1.
        for bond in declared_background()["bonds"]:
            offset = 1e-4
            carrier = math.pi / (2 * spacing)
            upper = 2 * bond * math.cos((carrier + offset) * spacing)
            lower = 2 * bond * math.cos((carrier - offset) * spacing)
            self.assertAlmostEqual(abs((upper - lower) / (2 * offset)), 2 * spacing * bond, delta=1e-8)

    def test_declared_background_has_shared_positive_nonuniform_coefficients(self):
        model = declared_background()
        self.assertNotEqual(model["bonds"][0], model["bonds"][1])
        mid_lapse = (np.array(model["lapse"][:-1]) + np.array(model["lapse"][1:])) / 2
        np.testing.assert_allclose(mid_lapse / np.array(model["spatial_midpoints"]), 2 * np.array(model["bonds"]), atol=1e-14)

    def test_joint_records_match_asymmetric_chain_and_clock_formulas(self):
        model = declared_background()
        expected = asymmetric_probabilities(0.8 * np.array(model["bonds"]))
        for index, direction in enumerate(DIRECTIONS):
            trial = BackgroundTrial(direction, "Q", model["bonds"], (1., 1.3))
            trial.arrive()
            table = trial.probabilities("Y")
            self.assertAlmostEqual(sum(prob for prob, outcome in zip(table, RECORDS) if outcome[2]), expected[index])
            self.assertAlmostEqual(sum(prob for prob, outcome in zip(table, RECORDS) if outcome[3]), expected[2])
            self.assertAlmostEqual(sum(prob for prob, outcome in zip(table, RECORDS) if outcome[1] == 0),
                                   (1 + math.sin(1.2 * 1.3 * 0.8)) / 2)
            self.assertAlmostEqual(sum(prob for prob, outcome in zip(table, RECORDS) if outcome[2] and outcome[3]), 0)

    def test_reverse_middle_records_recover_both_bonds_without_endpoint_fit(self):
        for bond_times in ((0.4, 0.2), (0.3, 0.6), (0.2, 0.2)):
            probabilities = asymmetric_probabilities(bond_times)
            np.testing.assert_allclose(inverse_bonds(*probabilities[:2]), bond_times, rtol=0, atol=1e-14)

    def test_clock_normalization_gives_the_same_geometry_for_both_probe_types(self):
        model = declared_background()
        time = 0.8
        for bare in PROBES.values():
            geometry = geometry_from_phases(bare * time / bare, bare * 1.3 * time / bare,
                                            time * np.array(model["bonds"]))
            self.assertAlmostEqual(geometry["lapse_C_over_A"], 1.3)
            np.testing.assert_allclose(geometry["spatial_midpoints"], [1.1, 1.3], atol=1e-14)

    def test_full_joint_state_and_post_readout_audit_still_match_shared_real_encoding(self):
        model = declared_background()
        trial = BackgroundTrial("A_to_C", "P", model["bonds"], (1., 1.3))
        trial.arrive()
        initial_count = len(trial.probabilities("X"))
        propagation_events = [event for event in trial.data.history if event["kind"] == "propagation"]
        self.assertTrue(all(event["edge_couplings"]["A_B"] == model["bonds"][0]
                    and "coupling" not in event for event in propagation_events))
        records = trial.collect_after_audit("X")
        self.assertEqual(initial_count, 16)
        self.assertAlmostEqual(sum(record["probability"] for record in records), 1)
        np.testing.assert_allclose(trial.encoded, encode_state(trial.joint), rtol=0, atol=3e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(trial.joint).min(), -1e-12)
        self.assertTrue(all(record["all_records_at_B"] == 1.6 for record in records))
        with self.assertRaises(ValueError):
            trial.readout("Y")

    def test_interval_propagation_covers_probability_error_box_corners(self):
        model = declared_background()
        middle = asymmetric_probabilities(0.8 * np.array(model["bonds"]))
        clocks = np.array([(1 + math.cos(0.8)) / 2, (1 + math.sin(0.8)) / 2,
                           (1 + math.cos(1.04)) / 2, (1 + math.sin(1.04)) / 2])
        for signs in product((-1, 1), repeat=6):
            observed = clocks + 0.01 * np.array(signs[:4])
            left = np.r_[observed, middle[0] + 0.01 * signs[4], middle[2]]
            right = np.r_[observed, middle[1] + 0.01 * signs[5], middle[2]]
            intervals = estimate_probe(left, right, 1.)
            for bounds, truth in zip(intervals["spatial_midpoints"], model["spatial_midpoints"]):
                self.assertLessEqual(bounds[0], truth + 1e-12)
                self.assertGreaterEqual(bounds[1], truth - 1e-12)
            bounds = intervals["predicted_endpoint_probability_not_used_in_fit"]
            self.assertLessEqual(bounds[0], middle[2] + 1e-12)
            self.assertGreaterEqual(bounds[1], middle[2] - 1e-12)

    def test_count_only_witness_explains_both_probe_types_and_heldout_endpoints(self):
        report = simulate_shared_background()
        self.assertTrue(report["shared_parameter_witness"]["all_means_within_declared_error"])
        self.assertLess(report["shared_parameter_witness"]["maximum_heldout_endpoint_residual"], 0.01)
        for probe in PROBES:
            intervals = report["per_probe_intervals"][probe]
            self.assertLess(intervals["lapse_C_over_A"][0], 1.3)
            self.assertGreater(intervals["lapse_C_over_A"][1], 1.3)
            for bounds, truth in zip(intervals["spatial_midpoints"], [1.1, 1.3]):
                self.assertLess(bounds[0], truth)
                self.assertGreater(bounds[1], truth)

    def test_reconstructed_common_background_is_executable_in_the_six_slot_engine(self):
        report = simulate_shared_background()
        witness = report["shared_parameter_witness"]
        geometry = witness["geometry"]
        for probe in PROBES:
            for index, direction in enumerate(DIRECTIONS):
                trial = BackgroundTrial(direction, probe, geometry["g_edges_over_lapse_A"],
                                        (1., geometry["lapse_C_over_A"]), witness["lapse_A_times_T"] / 2)
                trial.arrive()
                predicted = []
                for slot in (0, 1):
                    for basis in ("X", "Y"):
                        table = trial.probabilities(basis)
                        predicted.append(sum(probability for probability, record in zip(table, RECORDS) if record[slot] == 0))
                table = trial.probabilities("X")
                predicted.extend(sum(probability for probability, record in zip(table, RECORDS) if record[slot] == 1)
                                 for slot in (2, 3))
                self.assertLess(np.max(np.abs(np.array(predicted) - report["count_means"][probe][direction])), 0.01)

    def test_count_derived_prediction_guides_a_new_delayed_feedback_task(self):
        task = simulate_shared_background()["calibrated_fresh_feedback"]
        self.assertTrue(task["executed"])
        self.assertEqual([branch["action"] for branch in task["branches"]], ["keep", "flip"])
        self.assertLess(task["predicted_success_interval"][0], task["simulation_success"])
        self.assertGreater(task["predicted_success_interval"][1], task["simulation_success"])
        self.assertGreater(task["two_bit_code_ready_at_A_and_new_task_start"], task["calibration_archive_ready_at_B"])
        self.assertGreater(task["final_decision_parameter"], task["two_bit_code_ready_at_A_and_new_task_start"])
        self.assertTrue(all(branch["decision_parameter_from_task_start"] >= 0.8 for branch in task["branches"]))

    def test_sample_guarantee_counts_and_communication_are_explicit(self):
        report = simulate_shared_background()
        budget = report["sample_budget"]
        self.assertLessEqual(budget["union_failure_bound"], 0.01)
        self.assertEqual(budget["groups"], 8)
        for directions in report["simulated_counts"].values():
            for groups in directions.values():
                self.assertTrue(all(sum(counts) == budget["trials_per_group"] for counts in groups.values()))
        self.assertEqual(report["resource_per_trial"]["formal_outcome_bits"], 5)
        self.assertEqual(report["resource_per_trial"]["total_classical_bit_hops"], 6)
        self.assertEqual(report["resource_per_trial"]["resets"], 0)

    def test_input_contract_rejects_invalid_or_out_of_window_data(self):
        with self.assertRaises(ValueError):
            BackgroundTrial("A_to_C", "P", (10., 10.), (1., 1.3))
        with self.assertRaises(ValueError):
            BackgroundTrial("A_to_C", "P", (0.5, 0.5), (10., 1.))
        with self.assertRaises(ValueError):
            means_from_counts({"X": [0] * 16, "Y": [0] * 16})
        with self.assertRaises(ValueError):
            inverse_bonds(0.8, 0.8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    tests = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SharedBackgroundTests))
    if not tests.wasSuccessful():
        raise SystemExit(1)
    report = simulate_shared_background()
    report["automated_checks"] = {"run": tests.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("shared_background_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"tests": tests.testsRun, "budget": report["sample_budget"],
                      "estimates": report["per_probe_intervals"], "common_witness": report["shared_parameter_witness"]}, indent=2))


if __name__ == "__main__":
    main()