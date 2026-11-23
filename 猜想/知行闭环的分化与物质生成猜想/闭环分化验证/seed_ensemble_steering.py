"""Steering the two return ports of the first physical birth in note 40."""

import argparse
import json
from pathlib import Path

import numpy as np

from public_choice_quantum_return import H, I, S, X, Y, Z, equal
from quantum_recursive_ports import seed, tick


RESULT = Path(__file__).with_name("seed_ensemble_steering_results.json")


def projector(v):
    return np.outer(v, v.conj())


def born_pair():
    state, _ = tick(seed(0, 0), 0, 7)
    # Active B, dormant C/D, both action ports blank, old return initially 0.
    pair = state[1, 0, 0, :2, :2, 0, 0, 0]
    equal(pair, I / np.sqrt(2))
    equal(np.linalg.norm(pair), np.linalg.norm(state[..., 0]))
    return pair


def remote_state(amplitudes, effect):
    # Rows are B, columns are A. Contract the A indices with its effect.
    return np.einsum("ba,dc,ca->bd", amplitudes, amplitudes.conj(), effect)


def basis_runs(pair):
    for axis, gate in ((Z, I), (X, H), (Y, H @ S.conj().T)):
        for bit in (0, 1):
            effect = projector(gate.conj().T[:, bit])
            equal(effect, (I + (-1)**bit * axis) / 2)
            remote = remote_state(pair, effect)
            equal(remote, (I + (-1)**bit * axis.T) / 4)
            equal(np.trace(remote), 0.5)
    target = (I + Y) / 2
    nearest = I / 2
    equal(np.linalg.norm(target - nearest, ord="nuc") / 2, 0.5)


def three_result_instrument(pair):
    keep, choose = np.sqrt([1/3, 2/3])
    rotation = np.array([[keep, -choose], [choose, keep]])
    # A controls an actual rotation of the new blank flag F only when A=0.
    gate = np.kron(np.diag([1, 0]), rotation) + np.kron(np.diag([0, 1]), I)
    equal(gate.T @ gate, np.eye(4))
    flagged = np.einsum("cfak,ba,k->bcf", gate.reshape(2, 2, 2, 2), pair, [1, 0])
    recorded = [flagged[:, :, 1] @ flagged[:, :, 1].conj().T]
    rotate_after_zero = flagged[:, :, 0] @ H.T
    recorded.extend(projector(v) for v in rotate_after_zero.T)
    rays = np.array([[1, 0], [0.5, np.sqrt(3)/2], [0.5, -np.sqrt(3)/2]])
    for conditional_state, ray in zip(recorded, rays):
        equal(conditional_state, projector(ray) / 3)
    equal(sum(recorded), I / 2)
    equal(np.abs(rays @ rays.T)**2, np.eye(3)*0.75 + 0.25)


def general_ensemble(parts):
    rho = sum(parts)
    eig, vec = np.linalg.eigh(rho)
    support = eig > 1e-12
    columns = vec[:, support]
    sqrt_rho = (columns * np.sqrt(eig[support])) @ columns.conj().T
    # Conditional construction of the purification from a maximally entangled
    # pair, once tunable weak measurements and two-sided feedforward are granted.
    d = len(rho)
    amplitudes = np.sqrt(np.where(support, eig, 0.0))
    kraus = [np.diag(np.roll(amplitudes, k)) for k in range(d)]
    equal(sum(a.T @ a for a in kraus), np.eye(d))
    for k, operator in enumerate(kraus):
        correction = np.roll(np.eye(d), -k, axis=0)
        branch = operator / np.sqrt(d)
        corrected = vec @ correction @ branch @ correction.T @ vec.conj().T
        equal(corrected, sqrt_rho / np.sqrt(d))
    inverse = (columns / np.sqrt(eig[support])) @ columns.conj().T
    complement = np.eye(len(rho)) - columns @ columns.conj().T
    effects = [(inverse @ part @ inverse).T for part in parts]
    effects[0] += complement.T  # Complete the POVM on an unused kernel.
    equal(sum(effects), np.eye(len(rho)))
    equal(sqrt_rho @ sqrt_rho.conj().T, rho)
    for effect, part in zip(effects, parts):
        assert min(np.linalg.eigvalsh(effect)) > -1e-12
        equal(remote_state(sqrt_rho, effect), part)
    return int(sum(support))


def run():
    pair = born_pair()
    basis_runs(pair)
    three_result_instrument(pair)
    # Complex qutrit eigenvectors, with off-diagonal ensemble terms and
    # mixed members; the same purification handles two partitions of each rho.
    d = 3
    fourier = np.exp(2j*np.pi*np.outer(np.arange(d), np.arange(d))/d) / np.sqrt(d)
    u, v, w = fourier.T
    states = [projector(u), projector(v), projector((u + 1j*v)/np.sqrt(2))]
    ranks = []
    for third_weight in (0.0, 0.2):
        parts = [(1 - third_weight) * a / 3 for a in states]
        parts.append(third_weight * projector(w))
        rank = general_ensemble(parts)
        equal(general_ensemble([parts[0] + parts[2], parts[1] + parts[3]]), rank)
        ranks.append(rank)
    assert ranks == [2, 3]
    return {
        "seed": "note 40 first birth with old return zero",
        "remote_marginal": "I/2",
        "binary_reads_checked": ["Z -> Z", "X -> X", "Y -> -Y"],
        "closest_real_state_trace_distance_to_pure_Y": 0.5,
        "trine_outcome_probabilities": [1/3] * 3,
        "trine_pairwise_squared_overlap": 0.25,
        "extra_flag_blank_qubits": 1,
        "squared_flag_rotation_keep_amplitude": 1/3,
        "qutrit_test_ranks": ranks,
        "deterministic_schmidt_preparation_checked": True,
        "ensemble_partitions_per_qutrit_state": 2,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    result = run()
    if arguments.check:
        if result != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
