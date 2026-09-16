"""Round 133: coherent optimal intake compiled solely into old Ry and YX gates.

Input wires RA,RB,A,B become h,R,A,B; two pure environment rebits are added.
All flags and environment remain coherent. Only the three-wire common interface
is restricted. Angles and joint gate access are explicit ideal-control resources.
"""

import argparse
from dataclasses import dataclass
from functools import lru_cache
import json
import math
from pathlib import Path
import unittest

import numpy as np

from approximate_subject_joining import (
    apply_kraus, deterministic_joining_kraus, encoded_host_marginal,
    logical_newcomer, trace_distance,
)
from approximate_joining_optimality import choi_from_kraus
from bipartite_composition import interaction
from common_orientation_structure import encode_state
from encoded_composition_audit import independent_encoding
from flagged_subject_joining import intake_matrix, joining_kraus, random_state
from quantum_interface_audit import PAULI_X, PAULI_Y, PAULI_Z, rotation_unitary


@dataclass(frozen=True)
class Gate:
    kind: str
    sites: tuple
    angle: float

    def matrix(self):
        if self.kind == "Ry" and len(self.sites) == 1:
            return rotation_unitary(self.angle).real
        if self.kind == "YX" and len(self.sites) == 2:
            return interaction(self.angle, "YX").real
        raise ValueError("Only original Ry and oriented YX primitives are allowed.")


def ry(site, angle):
    return Gate("Ry", (site,), float(angle))


def yx(first, second, angle):
    return Gate("YX", (first, second), float(angle))


def apply_gate(gate, columns, systems):
    """Apply a physical gate to rows. Axis permutations are array indexing only."""
    sites = gate.sites
    if len(set(sites)) != len(sites) or any(s < 0 or s >= systems for s in sites):
        raise ValueError("Distinct physical wires in range are required.")
    rest = tuple(i for i in range(systems) if i not in sites)
    order = sites + rest + (systems,)
    tensor = columns.reshape((2,) * systems + (columns.shape[1],)).transpose(order)
    transformed = gate.matrix() @ tensor.reshape(2**len(sites), -1)
    return transformed.reshape(tensor.shape).transpose(np.argsort(order)).reshape(columns.shape)


def circuit_matrix(gates, systems, columns=None):
    result = np.eye(2**systems) if columns is None else columns.copy()
    for gate in gates:
        result = apply_gate(gate, result, systems)
    return result


def controlled_ry(control, target, angle):
    """P0 I + P1 Ry(angle): one YX and three original local rotations."""
    return (ry(control, math.pi/2), yx(target, control, -angle/2),
            ry(control, -math.pi/2), ry(target, angle/2))


def controlled_yx(control, first, second, angle):
    """Exact three-wire control, compiled into four two-wire YX rotations."""
    return (yx(first, control, -math.pi/2), yx(control, second, -angle/2),
            yx(first, control, math.pi/2), yx(first, second, angle/2))


def controlled_yz(control, first, second, angle):
    return ((ry(second, math.pi/2),) + controlled_yx(control, first, second, angle)
            + (ry(second, -math.pi/2),))


def controlled_z(control, target, helper, negative=False):
    """Borrowed helper is returned as an operator identity, even if entangled."""
    return (controlled_yz(control, helper, target, math.pi)
            + controlled_ry(control, helper, math.pi if negative else -math.pi))


def intake_circuit(helper=4):
    # W = diag(1,1,-1,1) U_XY(pi/2); its 4x4 determinant is -1.
    # A third, borrowed wire implements the conditional reflection with SO gates.
    return ((yx(1, 0, math.pi/2),) + controlled_z(0, 1, helper, negative=True))


def repair_circuit(target=3, first_environment=4, second_environment=5):
    middle = 2 * math.asin(1/math.sqrt(3))
    return (controlled_yx(0, first_environment, target, math.pi/4)
            + controlled_yz(0, second_environment, target, middle)
            + controlled_yx(0, first_environment, target, math.pi/4))


def endpoint_circuit(repair_host=False):
    gates = intake_circuit()
    if repair_host:
        gates += controlled_z(0, 1, 4)
    return gates + repair_circuit(2 if repair_host else 3)


def resource_counts(gates, extra_pure_rebits):
    return {
        "extra_pure_rebit_initializations": extra_pure_rebits,
        "joint_yx_rotations": sum(g.kind == "YX" for g in gates),
        "local_ry_rotations": sum(g.kind == "Ry" for g in gates),
        "gate_count": len(gates),
        "joining_readouts": 0,
        "joining_classical_outcome_messages": 0,
        "total_joining_rebits": 4 + extra_pure_rebits,
        "active_common_rebits": 3,
        "retained_environment_rebits_including_orientation": 1 + extra_pure_rebits,
        "gate_or_ancilla_minimality_claimed": False,
        "transport_topology_and_gate_noise_costs_included": False,
    }


def isometry_and_kraus(gates, systems):
    extra = systems - 4
    columns = np.eye(2**systems)[:, ::2**extra]
    isometry = circuit_matrix(gates, systems, columns)
    environment = (0,) + tuple(range(4, systems))
    order = environment + (1, 2, 3, systems)
    kraus = isometry.reshape((2,) * systems + (16,)).transpose(order).reshape(-1, 8, 16)
    return isometry, tuple(kraus)


@lru_cache(maxsize=2)
def compiled_endpoint(repair_host=False):
    return isometry_and_kraus(endpoint_circuit(repair_host), 6)


class CompiledSubjectJoiningTests(unittest.TestCase):
    def test_controlled_rotations_match_the_entire_operator(self):
        for angle in (-1.17, 0., .31, math.pi):
            expected = np.zeros((8, 8))
            expected[:4, :4] = np.eye(4)
            expected[4:, 4:] = interaction(angle, "YX").real
            np.testing.assert_allclose(circuit_matrix(controlled_yx(0, 1, 2, angle), 3),
                                       expected, atol=8e-16)
            expected[4:, 4:] = interaction(angle, "YZ").real
            np.testing.assert_allclose(circuit_matrix(controlled_yz(0, 1, 2, angle), 3),
                                       expected, atol=8e-16)
            expected2 = np.zeros((4, 4))
            expected2[:2, :2] = np.eye(2)
            expected2[2:, 2:] = rotation_unitary(angle).real
            np.testing.assert_allclose(circuit_matrix(controlled_ry(0, 1, angle), 2),
                                       expected2, atol=8e-16)

    def test_borrowed_helper_returns_exactly_without_a_state_promise(self):
        for negative in (False, True):
            actual = circuit_matrix(controlled_z(0, 1, 2, negative), 3)
            phase = np.diag([1., 1., -1. if negative else 1., 1. if negative else -1.])
            np.testing.assert_allclose(actual, np.kron(phase, np.eye(2)), atol=9e-16)
        actual = circuit_matrix(intake_circuit(2), 3)
        np.testing.assert_allclose(actual, np.kron(intake_matrix((1, 1)), np.eye(2)), atol=1e-15)
        self.assertAlmostEqual(np.linalg.det(intake_matrix((1, 1))), -1.)
        self.assertAlmostEqual(np.linalg.det(actual), 1.)

    def test_coherent_environment_has_exact_four_nonzero_kraus_labels(self):
        _, ks = compiled_endpoint()
        old = joining_kraus((2, 2))
        expected = [np.zeros((8, 16)) for _ in range(8)]
        expected[0] = old[(0, 0)]
        expected[4] = old[(0, 1)]/math.sqrt(3)
        expected[5] = np.kron(np.eye(4), PAULI_Z.real) @ old[(0, 1)]/math.sqrt(3)
        expected[6] = np.kron(np.eye(4), PAULI_X.real) @ old[(0, 1)]/math.sqrt(3)
        np.testing.assert_allclose(ks, expected, atol=1.2e-15)

    def test_full_channel_equals_round_127_including_unpromised_inputs(self):
        _, ks = compiled_endpoint()
        np.testing.assert_allclose(choi_from_kraus(ks),
                                   choi_from_kraus(deterministic_joining_kraus()), atol=2e-15)
        np.testing.assert_allclose(sum(k.T@k for k in ks), np.eye(16), atol=3e-15)

    def test_unknown_states_host_preservation_and_worst_error(self):
        rng = np.random.default_rng(133)
        _, ks = compiled_endpoint()
        for pauli in (PAULI_X, PAULI_Y, PAULI_Z):
            for sign in (-1, 1):
                host, new = random_state(rng, 2), (np.eye(2)+sign*pauli)/2
                output = apply_kraus(ks, independent_encoding((host, new)))
                np.testing.assert_allclose(output, encode_state(np.kron(host, logical_newcomer(new))), atol=1e-15)
                np.testing.assert_allclose(encoded_host_marginal(output, 2, 2), encode_state(host), atol=2e-15)
                self.assertAlmostEqual(trace_distance(output, encode_state(np.kron(host, new))), 1/6)

    def test_whole_circuit_and_external_purification_are_reversible(self):
        gates = endpoint_circuit()
        whole = circuit_matrix(gates, 6)
        np.testing.assert_allclose(whole.T@whole, np.eye(64), atol=4e-15)
        rng = np.random.default_rng(233)
        # Three columns are an untouched external system, with arbitrary correlations.
        original = rng.normal(size=(64, 3))
        original /= np.linalg.norm(original)
        changed = circuit_matrix(gates, 6, original)
        inverse = tuple(Gate(g.kind, g.sites, -g.angle) for g in reversed(gates))
        np.testing.assert_allclose(circuit_matrix(inverse, 6, changed), original, atol=5e-16)

    def test_orientation_branches_retain_their_distinct_conditional_errors(self):
        _, ks = compiled_endpoint()
        host = (np.eye(2)+.6*PAULI_Y)/2
        new = (np.eye(2)+PAULI_X)/2
        source, ideal = independent_encoding((host, new)), encode_state(np.kron(host, new))
        for indices, error in ((range(4), 0.), (range(4, 8), 1/3)):
            branch = apply_kraus([ks[i] for i in indices], source)
            self.assertAlmostEqual(np.trace(branch), .5)
            self.assertAlmostEqual(trace_distance(2*branch, ideal), error)

    def test_resources_are_actual_primitive_counts_and_no_measurements(self):
        gates = endpoint_circuit()
        counts = resource_counts(gates, 2)
        self.assertEqual(counts["joint_yx_rotations"], 18)
        self.assertEqual(counts["local_ry_rotations"], 7)
        self.assertEqual(counts["joining_readouts"], 0)
        self.assertEqual(counts["retained_environment_rebits_including_orientation"], 3)
        self.assertTrue(all(np.isrealobj(g.matrix()) for g in gates))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CompiledSubjectJoiningTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    report = {
        "round": 133,
        "input_order": ["RA", "RB", "A", "B", "E1_zero", "E2_zero"],
        "output_order": ["coherent_h", "R", "A", "B", "E1", "E2"],
        "primitives": ["Ry", "oriented YX"],
        "ideal_gate_angles_and_joint_access_assumed": True,
        "unknown_source_labels_or_purifications_accessed": False,
        "orientation_measurement_assumed": False,
        "entire_channel_matches_round_127": True,
        "old_host_preserved_on_promised_inputs": True,
        "worst_common_output_trace_error_exact": "1/6",
        "negative_orientation_conditional_error_exact": "1/3",
        "arbitrary_external_correlations_recoverable_by_whole_inverse": True,
        "resources": resource_counts(endpoint_circuit(), 2),
        "automated_checks": {"run": checks.testsRun, "failures": 0, "errors": 0},
    }
    if args.write_results:
        Path(__file__).with_name("compiled_subject_joining_results.json").write_text(
            json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
