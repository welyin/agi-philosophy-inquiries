"""Use an unwritten birth flag as a coherent controlled-controlled-NOT."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_consensus_dispatch import CNOT, equal
from quantum_loop_inheritance import I, X, Z, ZERO, controlled_flip
from quantum_recursive_ports import ACCEPT, ROUTER, V, apply_axes


RESULT = Path(__file__).with_name("seed_router_logical_computation_results.json")


def target_gate():
    matrix = np.zeros((8, 8))
    for a, b, c in product((0, 1), repeat=3):
        input_index = 4*a + 2*b + c
        output_index = 4*a + 2*b + (c ^ (a & b))
        matrix[output_index, input_index] = 1
    return matrix


def compute_gate():
    # Physical factors: flag, parent (S,R), spare a,b, independent target c.
    parent = V @ np.kron(ZERO, ZERO)
    equal(ACCEPT @ parent, parent)
    injection = np.kron(np.kron(ZERO, parent)[:, None], np.eye(8))
    invert_controls = np.kron(np.eye(8), np.kron(np.kron(X, X), I))
    router = np.kron(ROUTER, I)
    encoded = router @ invert_controls @ injection
    copied = apply_axes(encoded.reshape(2, 4, 2, 2, 2, 8), CNOT, (0, 4)).reshape(64, 8)
    decode = invert_controls @ router.T
    target = target_gate()
    equal(ROUTER.T, ROUTER)
    equal(decode @ copied, injection @ target)
    # A tagged external reference includes all eight unknown input columns.
    reference = np.diag(np.exp(1j * np.arange(8) / 3)) / np.sqrt(8)
    equal(decode @ copied @ reference, injection @ target @ reference)
    choi_fidelity = 0.0
    logical_branches = []
    for log in (0, 1):
        # Orthogonal, inaccessible copy of the intermediate flag. These are
        # its two Kraus branches after tracing out the copy.
        masked = copied.reshape(2, 32, 8).copy()
        masked[1 - log] = 0
        branch = decode @ masked.reshape(64, 8)
        logical = injection.T @ branch
        filter_diagonal = [int(((i >> 1) == 3) == log) for i in range(8)]
        equal(logical, target @ np.diag(filter_diagonal))
        equal(branch, injection @ logical)
        choi_fidelity += abs(np.trace(target.T @ logical) / 8)**2
        logical_branches.append(logical)
    equal(choi_fidelity, 5/8)
    witness = (np.eye(8)[:, 0] + np.eye(8)[:, 6]) / np.sqrt(2)
    ideal = target @ witness
    leakage_score = sum(abs(ideal @ op @ witness)**2 for op in logical_branches)
    equal(leakage_score, 0.5)


def encoded_phase():
    # Peel the two given CNOTs from V. R_y(pi/2) Z = H.
    extracted = controlled_flip("S") @ V @ controlled_flip("R") @ np.kron(Z, I)
    h = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    equal(extracted, np.kron(h, I))
    flip_real_imag = CNOT @ np.kron(I, h) @ CNOT @ np.kron(I, h)
    quarter_turn = X @ Z
    equal(flip_real_imag, np.block([[I, np.zeros((2, 2))], [np.zeros((2, 2)), quarter_turn]]))
    phase = np.diag([1, 1j])
    for state in (np.array([2-1j, 3j])/np.sqrt(14), np.array([1, 1])/np.sqrt(2)):
        def encode(vector):
            return np.column_stack((vector.real, vector.imag)).ravel()

        actual = flip_real_imag @ encode(state)
        equal(actual, encode(phase @ state))
        equal(np.square(actual).reshape(2, 2).sum(axis=1), np.abs(phase @ state)**2)


def real_control_ranks():
    dimension = 8
    identity = np.eye(dimension)
    uniform = np.ones(dimension) / np.sqrt(dimension)
    reflect_uniform = identity - 2 * np.outer(uniform, uniform)
    tangents = []
    for point in identity:
        reflect_point = identity - 2 * np.outer(point, point)
        rotation = reflect_point @ reflect_uniform
        equal(np.trace(rotation) - (dimension - 2), 4/dimension - 2)
        tangents.append(np.outer(point, uniform) - np.outer(uniform, point))
    brackets = [a @ b - b @ a for i, a in enumerate(tangents) for b in tangents[:i]]
    initial_rank = np.linalg.matrix_rank(np.reshape(tangents, (-1, dimension**2)), tol=1e-10)
    enlarged_rank = np.linalg.matrix_rank(np.reshape(tangents + brackets, (-1, dimension**2)), tol=1e-10)
    assert (initial_rank, enlarged_rank) == (7, 28)
    return [int(initial_rank), int(enlarged_rank)]


def run():
    compute_gate()
    encoded_phase()
    return {
        "unknown_input_qubits": 3,
        "basis_columns_with_relative_phases_verified": 8,
        "partner_dimension_in_test": 8,
        "birth_scatterings_per_toffoli": 2,
        "flag_and_parent_restored_on_entire_input_code": True,
        "ideal_choi_fidelity_without_flag_leak": 1,
        "choi_fidelity_if_flag_log_is_inaccessible": 0.625,
        "equal_predicate_sector_superposition_fidelity_after_flag_leak": 0.5,
        "exact_hadamard_on_native_target_with_spectator": True,
        "exact_encoded_quarter_phase_via_one_extra_real_bit": True,
        "three_qubit_real_skew_ranks_before_and_after_brackets": real_control_ranks(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    observed = run()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != observed:
            raise SystemExit("Recomputation differs from saved results.")
        print("Verified coherent birth flag, exact Toffoli, and encoded phase.")
    else:
        RESULT.write_text(json.dumps(observed, indent=2) + "\n", encoding="utf-8")
