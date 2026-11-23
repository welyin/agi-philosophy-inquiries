"""Test the additive time generators on the actual recursive birth isometry."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
SOURCE = HERE.parents[2] / CFG["birth_implementation"]
sys.path.insert(0, str(SOURCE.parent))
import quantum_recursive_lineage as lineage


def image_rank(matrices):
    features = np.column_stack([m.reshape(-1) for m in matrices])
    return int(np.linalg.matrix_rank(np.vstack([features.real, features.imag]), tol=1e-10))


def run():
    identity = np.eye(2)
    x = np.array([[0, 1], [1, 0]])
    y = lineage.Y
    z = lineage.Z
    rotation = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
    cycle = lineage.controlled_flip("S") @ np.kron(rotation, identity) @ lineage.controlled_flip("R")
    birth = cycle @ lineage.K
    np.testing.assert_allclose(birth.conj().T @ birth, identity, atol=1e-13)
    for sign, state in [(1, lineage.Y_PLUS), (-1, lineage.Y_MINUS)]:
        daughter = lineage.Y_MINUS if sign == 1 else lineage.Y_PLUS
        np.testing.assert_allclose(birth @ state, np.kron(state, daughter), atol=1e-13)
    # Remove energy zeros. The kernel describes all three traceless local generators.
    basis = (x, y, z)
    defects = ([birth @ p for p in basis]
               + [-np.kron(p, identity) @ birth for p in basis]
               + [-np.kron(identity, p) @ birth for p in basis])
    rank = image_rank(defects)
    assert rank == 7
    root = CFG["root_gap"]
    examples = []
    for fraction in CFG["examples_of_oriented_left_fractions"]:
        h_left = root * fraction * y / 2
        h_right = -root * (1-fraction) * y / 2
        error = (np.kron(h_left, identity)+np.kron(identity, h_right)) @ birth - birth @ (root*y/2)
        assert np.linalg.norm(error) < 1e-13
        examples.append({"oriented_left": root*fraction, "oriented_right": root*(1-fraction),
                         "sum_absolute_gaps": root*(abs(fraction)+abs(1-fraction))})
    # Compare the global pure outputs, with an admissible joint readout on the code.
    t = CFG["ramsey_time_in_inverse_gap_units"]
    zero = np.array([1, 0])
    u = np.cos(root*t/2)*identity-1j*np.sin(root*t/2)*y
    after = np.kron(u,u) @ birth @ zero
    before = birth @ u @ zero
    trace_distance = np.sqrt(max(0, 1-abs(np.vdot(after,before))**2))
    np.testing.assert_allclose(trace_distance,abs(np.sin(root*t/2)),atol=1e-13)
    for state in (after,before):
        amplitudes=state.reshape(2,2)
        np.testing.assert_allclose(amplitudes @ amplitudes.conj().T, identity/2, atol=1e-13)
    depth = CFG["balanced_depth"]
    gap = root/(2**depth)
    assert gap*(2**depth)==root
    return dict(traceless_unknowns=9, constraint_rank=rank, solution_dimension=9-rank,
                examples=examples, identical_port_global_trace_distance=round(float(trace_distance),12),
                single_port_trace_distance=0, balanced_leaves=2**depth, balanced_leaf_gap=gap,
                max_resolved_positive_leaves=int(np.floor(root/CFG["minimum_resolved_gap"]+1e-12)))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=run()
    result=HERE/"results.json"
    if args.check:
        assert json.loads(result.read_bytes())==out
    else:
        with result.open("x", encoding="utf8") as f:
            json.dump(out,f,ensure_ascii=False,indent=2)
            f.write("\n")
    print("1819 PASS")
