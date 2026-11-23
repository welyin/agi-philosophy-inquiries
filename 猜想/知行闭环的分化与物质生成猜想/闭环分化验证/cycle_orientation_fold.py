"""Heralded exact orientation-sector reflection and reversible folding of cyclic planes."""

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from cycle_memory_readout import CYCLE, H, W


def close(left, right):
    np.testing.assert_allclose(left, right, rtol=0, atol=3e-13)


def controlled(unitary):
    return np.block([[np.eye(len(unitary)), np.zeros_like(unitary)],
                     [np.zeros_like(unitary), unitary]])


def code_data():
    cycle = W.T @ CYCLE @ W
    orientation = (cycle - cycle.T)/np.sqrt(3)
    # Interchange the two physical bits, swapping 01 and 10.
    flip_bits = np.eye(4)[:, [0, 2, 1, 3]]
    flip = W.T @ flip_bits @ W
    close(flip @ cycle @ flip.T, cycle.T)
    k = np.kron(cycle, cycle.T)
    projector = (np.eye(4) + k + k@k)/3
    close(projector, (np.eye(4) - np.kron(orientation, orientation))/2)
    return cycle, orientation, np.kron(np.eye(2), flip), k, projector


def walk_check(k, projector):
    data_dim = 8   # Main reflection control and two real planes.
    # Three address bits, one invalid flag, data. Walsh prepare, select, recombine.
    walsh = np.kron(np.kron(H, H), H)
    h_data = np.kron(np.kron(walsh, np.eye(2)), np.eye(data_dim))
    selection = np.zeros_like(h_data)
    for addr in range(8):
        if addr == 0:
            op = -np.eye(data_dim)
        elif addr in (1, 2):
            op = controlled(k)
        elif addr in (3, 4):
            op = controlled(k@k)
        else:
            op = np.eye(data_dim)
        for invalid in range(2):
            src = (2*addr+invalid)*data_dim
            dst = (2*addr+(invalid ^ (addr >= 5)))*data_dim
            selection[dst:dst+data_dim, src:src+data_dim] = op
    walk = h_data @ selection @ h_data
    close(walk.T @ walk, np.eye(len(walk)))
    columns = walk[:, :data_dim]
    back = columns @ columns.T
    good = np.zeros_like(walk)
    good[:data_dim, :data_dim] = np.eye(data_dim)
    bad = np.eye(len(walk))-good
    reflect = 2*projector - np.eye(4)
    close(columns[:data_dim], (3/8)*controlled(reflect))
    success = 9/64
    for outcome, sign in ((back, 1), (np.eye(len(walk))-back, -1)):
        # Every accepted retry implements the same operation on all input columns.
        close(good @ outcome @ bad @ columns, sign*(1-success)*good @ columns)
        factor = 1-success if sign == 1 else success
        close(bad @ outcome @ bad @ columns, factor*bad @ columns)
    retry = 2*Fraction(9, 64)*Fraction(55, 64)
    return {"first_success": "9/64", "success_per_retry": str(retry),
            "mean_good_tests": str(1 + Fraction(55,64)/retry),
            "mean_forward_or_inverse_walk_calls": str(1 + 2*Fraction(55,64)/retry),
            "unfinished_after_8_good_tests": round((1-success)*float(1-retry)**7, 12)}


def fold_check(cycle, flip, proj):
    eye = np.eye(4)
    complement = eye-proj
    close(flip @ proj @ flip.T, complement)
    folding = np.vstack((proj, flip @ complement))
    close(folding.T @ folding, eye)
    close(folding @ folding.T, np.kron(np.eye(2), proj))
    h = np.kron(H, eye)
    clean = np.vstack((eye, np.zeros_like(eye)))
    close(controlled(flip) @ h @ controlled(2*proj-eye) @ h @ clean, folding)
    b_gate = np.kron(np.eye(2), cycle)
    effective = np.block([[b_gate, np.zeros_like(eye)], [np.zeros_like(eye), b_gate.T]])
    close(effective @ folding, folding @ b_gate)
    a_gate = np.kron(cycle, np.eye(2))
    close(np.kron(np.eye(2), a_gate) @ folding, folding @ a_gate)
    close(proj @ a_gate, proj @ b_gate)
    # Choose two orthogonal states whose differences reside wholly in the flag.
    _, vectors = np.linalg.eigh(proj)
    common = vectors[:, -1]
    for sign in (-1, 1):
        psi = (common + sign*flip @ common)/np.sqrt(2)
        flagged = (folding @ psi).reshape(2, 4)
        close(flagged.T @ flagged, np.outer(common, common))
        close(np.sum(flagged**2, axis=1), [0.5, 0.5])
    units = 0
    for i in range(4):
        for j in range(4):
            unit = eye[:, [i]] @ eye[[j], :]
            close(folding.T @ (folding @ unit @ folding.T) @ folding, unit)
            units += 1
    return {"input_real_dimension": 4, "common_plane_dimension": round(np.trace(proj)),
            "matrix_units_returned": units, "single_pair_classical_flag_overlap": "1/2"}


def physical_dispatch(flip, proj):
    # Physical labels 01,11,10 carry cyclic coordinates 0,1,2.
    active = [1, 3, 2]
    addition = np.eye(16)
    orbit_columns = np.zeros((16, 3))
    for j, left in enumerate(active):
        for k, right in enumerate(active):
            source = 4*left + right
            total = (j+k) % 3
            destination = 4*left + active[total]
            addition[:, source] = np.eye(16)[:, destination]
            orbit_columns[source, total] = 1/np.sqrt(3)
    close(addition.T @ addition, np.eye(16))
    uniform = np.zeros((4, 1))
    uniform[active, 0] = 1/np.sqrt(3)
    close(addition @ orbit_columns, np.kron(uniform, np.eye(4)[:, active]))
    raw_code = np.kron(W, W)
    common = orbit_columns @ W[active]
    close(common @ common.T, raw_code @ proj @ raw_code.T)
    prepared = np.kron(H, H)[:, [0]]
    prepared[0, 0] = 0
    close(prepared.T @ prepared, [[3/4]])
    close(prepared/np.sqrt(3/4), uniform)
    folded = np.vstack((proj, flip @ (np.eye(4)-proj)))
    outgoing = np.kron(np.eye(2), addition @ raw_code) @ folded
    # Dimensions are flag × A's two physical bits × B's two physical bits.
    attached = np.kron(np.eye(2), np.kron(uniform, W))
    delivered = attached.T @ outgoing
    close(outgoing, attached @ delivered)
    close(delivered @ delivered.T, np.eye(4))
    decoder = folded.T @ np.kron(np.eye(2), raw_code.T @ addition.T) @ attached
    for i in range(4):
        for j in range(4):
            unit = np.eye(4)[:, [i]] @ np.eye(4)[[j], :]
            returned = decoder @ delivered @ unit @ delivered.T @ decoder.T
            close(returned, unit)
    return {"physical_permutation_dimension": 16, "inactive_label_pairs_fixed": 7,
            "delivered_registers": ["B two-bit cycle", "coherent binary flag"],
            "three_label_reference_preparation": "3/4", "matrix_units_recovered": 16}


def sector_counts(orientation):
    rows = []
    for count in range(2, 6):
        size = 2**count
        projectors = [np.eye(size)]
        for site in range(1, count):
            pair = np.array([[1.0]])
            for index in range(count):
                local = orientation if index in (0, site) else np.eye(2)
                pair = np.kron(pair, local)
            passed = (np.eye(size)-pair)/2
            projectors = [p @ choice for p in projectors for choice in (passed, np.eye(size)-passed)]
        close(sum(projectors), np.eye(size))
        ranks = [int(np.linalg.matrix_rank(p, tol=1e-10)) for p in projectors]
        assert ranks == [2] * (size//2)
        for p in projectors:
            close(p @ p, p)
        # Each recorded sector contributes squared overlap (rank / size)^2
        # to the returned maximally entangled reference test.
        reference_overlap = sum(Fraction(rank, size)**2 for rank in ranks)
        assert reference_overlap == Fraction(2, size)
        rows.append({"loops": count, "sectors": len(ranks), "dimension_per_sector": 2,
                     "best_classical_flag_reference_overlap": str(reference_overlap)})
    return rows


def run():
    cycle, orientation, flip, k, proj = code_data()
    return {"controlled_reflection": walk_check(k, proj),
            "coherent_fold": fold_check(cycle, flip, proj),
            "physical_dispatch": physical_dispatch(flip, proj),
            "multiple_loops": sector_counts(orientation)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("cycle_orientation_fold_results.json")
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
