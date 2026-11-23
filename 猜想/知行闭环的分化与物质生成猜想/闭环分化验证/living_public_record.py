"""Check independent local birth records and the quantum data that the records leave intact."""

import argparse
import json
from pathlib import Path

import numpy as np

from quantum_consensus_dispatch import equal, route
from quantum_living_tree_recovery import (
    INPUT_DIM, leaf_correction, sample_x_and_keep_others, target_encoding,
)


RESULT = Path(__file__).with_name("living_public_record_results.json")
I_DATA = np.eye(INPUT_DIM)
I_SITE = np.eye(2)
X = I_SITE[::-1]


def outer(left, right):
    return np.outer(left, right.conj())


def success(left, right):
    """Equal-prior optimal binary discrimination for density matrices."""
    return 0.5 + np.sum(np.abs(np.linalg.eigvalsh(left - right))) / 4


def copy_location_from_flags():
    """A and B independently copy their own spare flag to their own blank reader."""
    history, ancestor_code = route()
    images = []
    for location in range(2):
        # Axes: f_A, A/spare data, C, f_B, B/spare data, two evidence bits, input.
        tree = target_encoding(ancestor_code, location)
        images.append(tree.reshape(-1, INPUT_DIM))
    encoding = np.column_stack(images)
    equal(encoding.T @ encoding, np.eye(2 * INPUT_DIM))

    # The new reader bits (a,b) have values (1-f_A, f_B).
    # Perform the two controlled writes on physical rows, then compress back
    # into the direct sum of the known trees. Each Kraus operator also has a
    # different environment index: these maps correspond to tracing both readers.
    rows = encoding.reshape(2, 16, 4, 2, 16, 2, 2, 2 * INPUT_DIM)
    # Use the valid role assignment W at C, with unknown evidence at A and B.
    # A's evidence is on axis 5. Its X and Z are actual local operations on
    # the data hosted by the same participant who reads a definite flag.
    z = np.diag([1, -1])
    data_x = np.kron(np.eye(4), np.kron(X, I_SITE))
    data_z = np.kron(np.eye(4), np.kron(z, I_SITE))
    equal(np.flip(rows, axis=5).reshape(encoding.shape), encoding @ data_x)
    physical_z = np.diag(z).reshape(1, 1, 1, 1, 1, 2, 1, 1)
    equal((rows * physical_z).reshape(encoding.shape), encoding @ data_z)
    equal(data_x @ data_z, -data_z @ data_x)
    maps = {}
    for a in range(2):
        for b in range(2):
            physical_sector = rows[1 - a, :, :, b].reshape(-1, 2 * INPUT_DIM)
            maps[a, b] = physical_sector.T @ physical_sector
            expected = np.kron(outer(I_SITE[a], I_SITE[a]), I_DATA) if a == b else 0
            equal(maps[a, b], expected)
    # Actual corrected birth has four recorded outcomes (location, X result).
    # Compress its maps into the same logical direct sum rather than granting
    # the pre-birth data arbitrary state preparation on the location factor.
    birth_maps = []
    for location in range(2):
        correction = leaf_correction(location)
        for outcome in range(2):
            raw = sample_x_and_keep_others(history, 1, location, outcome)
            actual = raw * correction[..., None] if outcome else raw
            actual = actual.reshape(-1, INPUT_DIM)
            compressed = encoding.T @ actual
            equal(encoding @ compressed, actual)
            birth_maps.append(compressed)
    return maps, birth_maps


def check_reader_effects(maps, birth_maps):
    # Work with all four matrix units of a single-reader effect (linearity).
    # Either reader reports the same diagonal observable. Its probabilities
    # at the older three-input interface are independent of the input.
    for row in range(2):
        for col in range(2):
            from_a = sum(maps[row, b].T @ maps[col, b] for b in range(2))
            from_b = sum(maps[a, row].T @ maps[a, col] for a in range(2))
            expected = np.kron(outer(I_SITE[row], I_SITE[row]), I_DATA) if row == col else 0
            equal(from_a, expected)
            equal(from_b, expected)
            before_birth = sum(k.T @ from_a @ k for k in birth_maps)
            equal(before_birth, I_DATA / 2 if row == col else 0)


def check_public_channel(maps):
    # If B is inaccessible but A is retained, the remaining channel takes
    # each input to tree + A. The second copy alone traces out the cross term.
    one_missing = [
        np.vstack([maps[a, b] for a in range(2)]) for b in range(2)
    ]
    two_missing = list(maps.values())
    logical_dim = 2 * INPUT_DIM
    within_count = 0
    for i in range(logical_dim):
        for j in range(logical_dim):
            unit = np.zeros((logical_dim, logical_dim))
            unit[i, j] = 1
            within = i // INPUT_DIM == j // INPUT_DIM
            within_count += int(within)
            recovered = sum(op @ unit @ op.T for op in two_missing)
            equal(recovered, unit if within else 0)
            delivered = sum(op @ unit @ op.T for op in one_missing)
            reader_state = I_SITE[i // INPUT_DIM]
            equal(delivered, np.kron(outer(reader_state, reader_state), unit) if within else 0)

    # A diagonal observable commutes with each within-tree matrix unit only
    # if the corresponding diagonal entries agree. Commuting with all basis
    # projections already forced it to be diagonal.
    adjacent_differences = [
        np.eye(logical_dim)[i] - np.eye(logical_dim)[i + 1]
        for i in range(logical_dim - 1) if (i + 1) % INPUT_DIM
    ]
    center_dimension = logical_dim - np.linalg.matrix_rank(adjacent_differences)
    equal(center_dimension, 2)
    for location in range(2):
        # Input data is maximally entangled with an untouched 8-dim partner.
        amplitudes = np.kron(I_SITE[location, :, None], I_DATA / np.sqrt(INPUT_DIM))
        saved_overlap = sum(abs(np.vdot(amplitudes, op @ amplitudes)) ** 2 for op in two_missing)
        equal(saved_overlap, 1)
    return within_count, int(center_dimension)


def record_comparisons():
    ket0, ket1 = I_SITE
    plus = (ket0 + ket1) / np.sqrt(2)
    minus = (ket0 - ket1) / np.sqrt(2)
    # A diagnostic coherent input across the two alternatives requires
    # preparation/control beyond the permanent, public location in note 43.
    markers = [np.array([0.6, 0.8]), np.array([0.8, 0.6])]
    r0 = np.kron(ket0, ket0)
    r1 = np.kron(*markers)
    recorded = np.stack((r0, r1)) / np.sqrt(2)
    location = recorded @ recorded.conj().T
    visibility = abs(2 * location[0, 1])
    equal(visibility, 0.48)
    guesses = [success(outer(ket0, ket0), outer(m, m)) for m in markers]
    equal(guesses, [0.9, 0.8])

    # Orthogonal Bell markers dephase the location but store the result only
    # in a joint correlation. Neither heir can read the location alone.
    bell_plus = np.diag(plus)
    bell_minus = np.diag(minus)
    equal(np.vdot(bell_minus, bell_plus), 0)
    equal(success(bell_plus @ bell_plus.T, bell_minus @ bell_minus.T), 0.5)
    equal(success(bell_plus.T @ bell_plus, bell_minus.T @ bell_minus), 0.5)
    equal(success(outer(bell_plus.ravel(), bell_plus.ravel()),
                  outer(bell_minus.ravel(), bell_minus.ravel())), 1)

    # Ideal redundant records: an X-basis read on one public copy spoils its
    # later Z agreement with the other. Orthogonal records still untouched
    # elsewhere keep the tree block diagonal. No tree operation is used.
    for bit in range(2):
        written = outer(I_SITE[bit], I_SITE[bit])
        x_unread = (written + X @ written @ X) / 2
        equal(x_unread[bit, bit], 0.5)

    # Before any durable copy, a coherent write followed by its inverse
    # restores the interference. Every register including the blank is returned.
    cnot = np.eye(4)[[0, 1, 3, 2]]
    before = np.kron(plus, ket0)
    after = cnot @ before
    equal(after.reshape(2, 2) @ after.reshape(2, 2).T, I_SITE / 2)
    equal(cnot @ after, before)
    return {
        "independent_marker_overlaps": [0.6, 0.8],
        "single_reader_optimal_success": np.round(guesses, 12).tolist(),
        "joint_trace_visibility": round(float(visibility), 12),
        "bell_marker_visibility": 0.0,
        "bell_marker_single_reader_success": [0.5, 0.5],
        "bell_marker_joint_reader_success": 1.0,
        "copy_Z_agreement_after_nonselective_X_read": 0.5,
        "coherent_write_and_inverse_restores_both_registers": True,
    }


def verify():
    maps, birth_maps = copy_location_from_flags()
    check_reader_effects(maps, birth_maps)
    within, center_dimension = check_public_channel(maps)
    return {
        "physical_source": "note 43 location-specific four-active-leaf encoders",
        "physical_local_records": ["1 - A spare activity flag", "B spare activity flag"],
        "new_blank_reader_bits": 2,
        "logical_data_dimension_per_location": INPUT_DIM,
        "matrix_units_checked": (2 * INPUT_DIM) ** 2,
        "block_diagonal_units_fixed": within,
        "off_diagonal_units_killed_with_one_inaccessible_copy": (2 * INPUT_DIM) ** 2 - within,
        "dimension_of_center": center_dimension,
        "each_readers_four_effect_units_pull_back_to_the_location_diagonal": True,
        "old_input_probability_effect_for_each_postbirth_location": "identity_8 / 2",
        "within_location_maximally_entangled_data_fidelity": 1.0,
        "local_A_evidence_X_and_Z_anticommute_and_preserve_both_flag_sectors": True,
        "coherent_location_preparation_is_a_separate_diagnostic": True,
        "comparisons": record_comparisons(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    computed = verify()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != computed:
            raise SystemExit("Saved public record results differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(computed, indent=2) + "\n", encoding="utf-8")
