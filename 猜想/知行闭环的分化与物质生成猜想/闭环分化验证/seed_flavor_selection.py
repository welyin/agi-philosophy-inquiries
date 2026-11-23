"""Test diagonal symmetry of the actual split and independent flavor multiplicity."""

import argparse
import json
from pathlib import Path

import numpy as np

from quantum_recursive_lineage import Y, Y_MINUS, Y_PLUS, tree_code
from quantum_recursive_ports import ENCODE, I, READY, V


RESULT = Path(__file__).with_name("seed_flavor_selection_results.json")
X = np.array([[0, 1], [1, 0]])
Z = np.diag([1, -1])


def equal(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)


def perpendicular(matrix, direction):
    return matrix - direction * (np.vdot(direction, matrix) / np.vdot(direction, direction))


def shared_rotation(split):
    # Central phase is immaterial for the induced isometric channel.
    drifts = []
    for generator in (X, Y, Z):
        output_action = np.kron(generator, I) + np.kron(I, generator)
        drift = output_action @ split - split @ generator
        drifts.append(perpendicular(drift, split).reshape(-1))
    rank = np.linalg.matrix_rank(np.stack(drifts, axis=1), tol=1e-11)
    assert rank == 3
    checks = {}
    for name, gate in (("I", I), ("X", X), ("Y", Y), ("Z", Z),
                       ("Hadamard", (X + Z) / np.sqrt(2))):
        left, right = np.kron(gate, gate) @ split, split @ gate
        joint_gate = split @ gate @ split.conj().T + np.eye(4) - split @ split.conj().T
        equal(joint_gate.conj().T @ joint_gate, np.eye(4))
        equal(joint_gate @ split, right)
        norm = np.linalg.norm(perpendicular(left, right))
        compatible = bool(norm < 1e-12)
        assert compatible == (name in ("I", "Z"))
        checks[name] = compatible
    # Product output rays and a non-product probe in the same physical code.
    for vec, other in ((Y_PLUS, Y_MINUS), (Y_MINUS, Y_PLUS)):
        equal(split @ vec, np.kron(vec, other))
    coeff = np.array([np.sqrt(0.3), np.exp(0.4j) * np.sqrt(0.7)])
    probe = coeff[0] * Y_PLUS + coeff[1] * Y_MINUS
    equal(np.linalg.svd((split @ probe).reshape(2, 2), compute_uv=False),
          sorted(abs(coeff), reverse=True))
    return int(rank), checks


def flavor_extension(split, families):
    # Take one independent gauge doublet with N equivalent family labels.
    # Passing the entire spectator through preserves the split on any loop input.
    identity = np.eye(families)
    gauge = np.kron(1j * (X + Z) / np.sqrt(2), identity)
    encode = np.kron(split, np.eye(2 * families))
    equal(np.kron(np.eye(4), gauge) @ encode, encode @ np.kron(I, gauge))
    # A spectator may even be correlated with the seed. Changing the spectrum
    # only induces a spectator unitary and cannot alter any loop-only statistic.
    joint = np.arange(1, 4 * families + 1) * np.exp(0.19j * np.arange(4 * families))
    joint /= np.linalg.norm(joint)
    root = joint.reshape(2, -1)
    expected = split @ (root @ root.conj().T) @ split.conj().T
    spectra = []
    for last_shift in (0, 1):
        diagonal = np.arange(1, families + 1)
        diagonal[-1] += last_shift
        mass = np.kron(I, np.diag(diagonal))
        equal(gauge @ mass, mass @ gauge)
        evolution = np.diag(np.exp(-0.37j * np.diag(mass)))
        output = (encode @ joint).reshape(4, -1) @ evolution.T
        equal(output @ output.conj().T, expected)
        spectra.append(diagonal.tolist())
    return {
        "families": families,
        "two_distinct_family_spectra_in_test_units": spectra,
        "conditional_one_loop_b_left_from_round1545": 7 - 4 * families / 3,
    }


def run():
    split = V @ READY
    equal(split, ENCODE)
    tree, _ = tree_code(((None, None), None), split)
    equal(tree.conj().T @ tree, I)
    rank, checks = shared_rotation(split)
    return {
        "three_actual_leaf_ports_encoded_root_rank": int(np.linalg.matrix_rank(tree)),
        "traceless_generators_tested_independently": rank,
        "dimension_of_projective_diagonal_symmetry_lie_algebra": 3 - rank,
        "sample_identical_parent_child_actions_up_to_phase": checks,
        "corresponding_joint_code_actions_all_match_the_parent": True,
        "unobserved_independent_flavor_extension": [flavor_extension(split, n) for n in (3, 6)],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        if result != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified split symmetry, three-leaf encoding and flavor extensions.")
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
