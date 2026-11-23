"""Recover an entangled root and two unknown evidence registers while the new tree stays active."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_consensus_dispatch import birth_map, equal, route
from quantum_recursive_ports import READY, V, Z_RETURN


RESULT = Path(__file__).with_name("quantum_living_tree_recovery_results.json")
SITES = ("A", "C", "B")
Z = np.diag([1, -1])
INPUT_DIM = 8  # root and two evidence qubits


def grow(incoming, where):
    """A or B splits; the other candidate and C continue their feedback."""
    return np.einsum(
        "ai,cj,bk,ijkx->acbx",
        birth_map(where == 0), V, birth_map(where == 1), incoming, optimize=True,
    )


def target_encoding(ancestor_code, where):
    ready_leaves = np.einsum("ia,jc,kb,acbx->ijkx", READY, READY, READY, ancestor_code)
    grown = grow(ready_leaves, where)
    target = np.zeros((32, 4, 32, 2, 2, INPUT_DIM))
    for u, v in product(range(2), repeat=2):
        target[:, :, :, u, v, 2 * u + v::4] = grown
    return target


def sample_x_and_keep_others(history, witness, where, outcome):
    """Lift note 41's coherent circuit to two unknown input wires by linearity.

    All of its message/control gates preserve the computational values of M.
    The stored three-plus amplitudes have coefficient 1/sqrt(8) in each M
    sector. For fixed u,v and a plus witness it is 1/sqrt(2) instead.
    """
    output = np.zeros((32, 4, 32, 2, 2, INPUT_DIM))
    other_sites = [site for site in range(3) if site != witness]
    for u, v in product(range(2), repeat=2):
        evidence = [0] * 3
        evidence[other_sites[0]], evidence[other_sites[1]] = u, v
        projected = np.zeros((4, 4, 4, 2))
        for w in range(2):
            evidence[witness] = w
            sector = (slice(None),) * 3 + tuple(evidence) + (where, where)
            projected += (2 / np.sqrt(2)) * (-1)**(outcome * w) * history[sector]
        output[:, :, :, u, v, 2 * u + v::4] = grow(projected, where)
    return output


def leaf_correction(where):
    """Only touch the four active leaves and two retained evidence registers."""
    ticked_z = V @ Z_RETURN @ V.T
    equal(ticked_z, np.kron(Z, Z))

    def endpoint(is_parent):
        cell = Z_RETURN if is_parent else ticked_z
        spare = Z_RETURN if is_parent else np.eye(4)
        # Physical register order: activity flag, old cell, candidate spare.
        return np.diag(np.kron(np.eye(2), np.kron(cell, spare)))

    return np.einsum(
        "a,c,b,u,v->acbuv", endpoint(where == 0), np.diag(ticked_z),
        endpoint(where == 1), np.diag(Z), np.diag(Z),
    )


def flag_and_code_checks(tensor, where):
    blocks = tensor.reshape(2, 16, 4, 2, 16, 2, 2, INPUT_DIM)
    for left, right in product(range(2), repeat=2):
        if (left, right) != (1 - where, where):
            equal(blocks[left, :, :, right], 0)


def verify():
    history, ancestor_code = route()
    logical_parity = np.kron(np.kron(Z, Z), Z)
    identity = np.eye(INPUT_DIM)
    for where in range(2):
        target = target_encoding(ancestor_code, where)
        flag_and_code_checks(target, where)
        encoder = target.reshape(-1, INPUT_DIM)
        equal(encoder.T @ encoder, identity)
        correcting_sign = leaf_correction(where)
        for witness in range(3):
            fidelities = np.zeros(2)
            for s in range(2):
                raw = sample_x_and_keep_others(history, witness, where, s)
                flag_and_code_checks(raw, where)
                uncorrected = raw.reshape(-1, INPUT_DIM)
                phase = (-1)**(where * s)
                error = logical_parity if s else identity
                equal(uncorrected, phase / 2 * (encoder @ error))
                equal(uncorrected.T @ uncorrected, identity / 4)
                corrected = raw * correcting_sign[..., None] if s else raw
                flag_and_code_checks(corrected, where)
                fixed = corrected.reshape(-1, INPUT_DIM)
                equal(fixed, phase / 2 * encoder)
                # Choi entanglement fidelity ignoring s, conditional on where.
                # Both outcomes together have probability 1/2.
                for index, matrix in enumerate((uncorrected, fixed)):
                    fidelities[index] += 2 * abs(np.trace(encoder.T @ matrix) / INPUT_DIM)**2
            equal(fidelities, [0.5, 1])
    return {
        "witness_sites_individually_checked": list(SITES),
        "unknown_qubits": ["root", "other evidence 1", "other evidence 2"],
        "input_code_dimension": INPUT_DIM,
        "location_and_X_outcome_blocks_per_witness_choice": 4,
        "input_columns_per_block": INPUT_DIM,
        "all_branch_probabilities_for_any_input": 0.25,
        "corrected_branch_isometry_factor": "(-1)^(location times X outcome)/2",
        "new_activity_flags_A_B_for_location_0_1": [[1, 0], [0, 1]],
        "active_cells_retained": 4,
        "correction_on_location_flags_and_unused_spare": "identity",
        "feedforward": "one X outcome bit to the four active leaves and the two retained evidence ports; gates also use the existing birth location",
        "ancillary_X_witness_consumed": 1,
        "reference_dimension_in_Choi_fidelity": INPUT_DIM,
        "Choi_fidelity_given_location_after_feedforward": 1.0,
        "Choi_fidelity_given_location_discarding_outcome_without_feedforward": 0.5,
        "unknown_evidence_states": "extension replacing two of note 41's plus evidence preparations; these inputs must be supplied",
        "all_full_physical_output_column_and_activity_flag_checks_passed": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = verify()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != result:
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
