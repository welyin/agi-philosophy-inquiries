"""Use a public bit to choose a probe basis and deliver the next common reading."""

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np

from quantum_recursive_ports import apply_axes


RESULT = Path(__file__).with_name("public_choice_quantum_return_results.json")
I = np.eye(2)
X = I[::-1]
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1])
H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
S = np.diag([1, 1j])
PLUS = H[:, 0]
CNOT = np.eye(4)[[0, 1, 3, 2]]
BASIS_READ = np.block([[H, np.zeros((2, 2))], [np.zeros((2, 2)), H @ S.conj().T]])


def equal(observed, expected):
    np.testing.assert_allclose(observed, expected, atol=1e-12, rtol=0)


def issue_reading(state, phase, control, readers):
    # W is axis 1. The basis-controlled instrument ends with W in the Z basis.
    state = apply_axes(state, np.diag([1, np.exp(1j * phase)]), (1,))
    state = apply_axes(state, BASIS_READ, (control, 1))
    first, second = readers
    state = apply_axes(state, CNOT, (1, first))
    return apply_axes(state, CNOT, (first, second))


def prob(bit, old, phase):
    return (1 + (-1)**bit * np.cos(phase - old * np.pi / 2)) / 2


def check_histories(old, phases):
    # W, one old record, and two separate pairs of initially blank new readers.
    state = np.zeros((2,) * 6, dtype=complex)
    state[old, :, 0, 0, 0, 0] = PLUS
    state = issue_reading(state, phases[0], 0, (2, 3))
    # Reuse this actual probe; the old reading stays in both parties' registers.
    state = apply_axes(state, CNOT, (2, 1))
    state = apply_axes(state, H, (1,))
    state = issue_reading(state, phases[1], 2, (4, 5))
    expected = np.zeros_like(state.real)
    for r, t in product(range(2), repeat=2):
        expected[old, t, r, r, t, t] = prob(r, old, phases[0]) * prob(t, r, phases[1])
    equal(np.abs(state)**2, expected)


def late_choice():
    # Calibration source hides its sign and basis; score matching-basis trials.
    strategy = (X + Y) / np.sqrt(2)
    effects = [(I + sign * strategy) / 2 for sign in (1, -1)]
    scores = []
    for basis, observable in enumerate((X, Y)):
        read = (H, H @ S.conj().T)[basis]
        equal(read.conj().T @ Z @ read, observable)
        for bit in (0, 1):
            rho = (I + (-1)**bit * observable) / 2
            equal((read @ rho @ read.conj().T)[bit, bit], 1)
            scores.append(np.trace(effects[bit] @ rho).real)
    bound = (1 + 1 / np.sqrt(2)) / 2
    equal(scores, [bound] * 4)
    for sx, sy in product((-1, 1), repeat=2):
        equal(np.linalg.eigvalsh(sx * X + sy * Y), [-np.sqrt(2), np.sqrt(2)])
    return round(float(bound), 12)


def run():
    for gate in (BASIS_READ, CNOT, H, S):
        equal(gate.conj().T @ gate, np.eye(len(gate)))
    for s in (0, 1):
        equal(np.linalg.matrix_power(Z, s) @ H[:, s], PLUS)
    for old in (0, 1):
        for phases in ((0, np.pi/2), (np.pi, -np.pi/2), (np.pi/3, np.pi/5)):
            check_histories(old, phases)
    # Same Z populations; the off-diagonal term changes the resulting record.
    for basis in (0, 1):
        for phase in (0, np.pi/3, np.pi/2, np.pi):
            probe = np.array([1, np.exp(1j * phase)]) / np.sqrt(2)
            rho = np.outer(probe, probe.conj())
            gate = (H, H @ S.conj().T)[basis]
            equal((gate @ rho @ gate.conj().T)[0, 0], prob(0, basis, phase))
            equal((gate @ np.diag(np.diag(rho)) @ gate.conj().T)[0, 0], 0.5)
    table = [[round(float(prob(0, old, phase)), 12) for old in (0, 1)]
             for phase in (0, np.pi/2, np.pi/3)]
    return {
        "basis_selected_by_public_bit_0_1": ["X", "Y"],
        "probability_next_bit_zero_for_phases_0_pi_over_2_pi_over_3": table,
        "two_round_schedules_verified_including_both_initial_sites": 6,
        "each_new_reading_has_two_matching_recipients": True,
        "original_site_and_first_transcript_preserved": True,
        "dephasing_probe_in_Z_before_reading_gives_probability_zero": 0.5,
        "late_basis_choice_success_with_quantum_probe": 1,
        "late_basis_choice_optimal_success_with_only_prior_classical_record": late_choice(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != result:
            raise SystemExit("Saved results differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
