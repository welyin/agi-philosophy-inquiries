"""Check closed unitary routing from one active loop to two over a clock window."""

import argparse
import json
from pathlib import Path

import numpy as np

from quantum_loop_inheritance import controlled_flip, Z


RESULT = Path(__file__).with_name("quantum_split_router_results.json")
TOL = 1e-12


def check_equal(actual, expected):
    np.testing.assert_allclose(actual, expected, atol=TOL, rtol=0)


def pointer_embedding(position, flag, data, length):
    """Embed columns of data in the internal (clock, activity flag) sector."""
    columns = np.zeros((length, 2, 16, data.shape[1]), dtype=data.dtype)
    columns[position, flag] = data
    return columns.reshape(32 * length, data.shape[1])


def clocked_tick(length, scatter, dynamics):
    tick = np.zeros((32 * length, 32 * length))
    for j in range(length):
        next_j = (j + 1) % length
        block = scatter @ dynamics if j == 0 else dynamics
        tick[32 * next_j:32 * (next_j + 1), 32 * j:32 * (j + 1)] = block
    return tick


def check_router(length=8):
    eye4 = np.eye(4)
    eye16 = np.eye(16)
    eye32 = np.eye(32)
    rotation = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
    loop = controlled_flip("S") @ np.kron(rotation, np.eye(2)) @ controlled_flip("R")
    parent_only = np.kron(loop, eye4)
    local_pair = np.kron(loop, loop)
    initial_dynamics = np.block([[parent_only, np.zeros_like(eye16)],
                                 [np.zeros_like(eye16), local_pair]])

    # The data order is S_A,R_A,S_B,R_B. Swap only the middle two factors.
    swap = eye16.reshape(2, 2, 2, 2, 16).swapaxes(1, 2).reshape(16, 16)
    pauli_y = np.array([[0, -1j], [1j, 0]])
    minus_yy = (eye4 - np.kron(pauli_y, pauli_y).real) / 2
    blank_b = np.diag([1, 0, 0, 0])
    source = np.kron(minus_yy, blank_b)
    arrival = swap @ source @ swap.T
    scatter = np.block([[eye16 - source, source @ swap.T],
                        [swap @ source, eye16 - arrival]])
    check_equal(source @ source, source)
    check_equal(scatter.T @ scatter, eye32)
    check_equal(scatter, scatter.T)
    assert np.linalg.matrix_rank(source) == 2

    # One and the same unitary advances the internal clock on every application.
    tick = clocked_tick(length, scatter, initial_dynamics)
    check_equal(tick.T @ tick, np.eye(32 * length))

    seeds = eye16[:, [0, 4]]  # |0,x,0,0> for x=0,1; test the whole logical span.
    before_route = parent_only @ seeds
    parent_outputs = before_route.reshape(4, 4, 2)[:, 0, :]
    expected_bells = np.array([[1, 0], [0, -1], [0, 1], [1, 0]]) / np.sqrt(2)
    check_equal(parent_outputs, expected_bells)
    check_equal(source @ before_route, before_route)
    split = swap @ before_route
    initial_code = pointer_embedding(0, 0, seeds, length)
    trajectory = tick @ initial_code
    code_orbit = [initial_code]
    pair_data = split.copy()
    for k in range(length):
        j = (k + 1) % length
        check_equal(trajectory, pointer_embedding(j, 1, pair_data, length))
        code_orbit.append(trajectory)
        if k + 1 < length:
            trajectory = tick @ trajectory
            pair_data = local_pair @ pair_data
    orbit_columns = np.column_stack(code_orbit)
    check_equal(orbit_columns.T @ orbit_columns, np.eye(2 * (length + 1)))
    check_equal(np.linalg.matrix_power(loop, 6), -eye4)
    check_equal(np.linalg.matrix_power(local_pair, 6) @ split, split)

    # Let the parent complete a full six-tick orbit before the routing opportunity.
    # Both logical columns acquire the same minus sign during these parent-only ticks.
    assert length >= 7
    delayed_run = pointer_embedding(length - 6, 0, seeds, length)
    parent_data = seeds.copy()
    for step in range(1, 7):
        delayed_run = tick @ delayed_run
        parent_data = parent_only @ parent_data
        check_equal(delayed_run, pointer_embedding((length - 6 + step) % length,
                                                   0, parent_data, length))
    check_equal(parent_data, -seeds)
    delayed_run = tick @ delayed_run
    child_data = -split
    check_equal(delayed_run, pointer_embedding(1, 1, child_data, length))
    for step in range(1, length):
        delayed_run = tick @ delayed_run
        child_data = local_pair @ child_data
        check_equal(delayed_run, pointer_embedding((1 + step) % length, 1, child_data, length))

    # A counterfactual phase change just after return; identical clock and router.
    flip_return_phase = np.kron(np.kron(np.eye(2), Z), eye4)
    rejected = flip_return_phase @ before_route
    post_scatter = scatter @ np.vstack((rejected, np.zeros_like(rejected)))
    check_equal(post_scatter, np.vstack((rejected, np.zeros_like(rejected))))
    sleeping = pointer_embedding(1, 0, rejected, length)
    for k in range(1, length):
        sleeping = tick @ sleeping
        rejected = parent_only @ rejected
        check_equal(sleeping, pointer_embedding((k + 1) % length, 0, rejected, length))

    # Test all twelve basis columns spanning "B occupied", including their coherences
    # by linearity; B is preserved in this entire subspace over the same window.
    occupied = eye16[:, [n for n in range(16) if n % 4 != 0]]
    occupied_run = pointer_embedding(0, 0, occupied, length)
    occupied_data = occupied.copy()
    for k in range(1, length + 1):
        occupied_run = tick @ occupied_run
        occupied_data = parent_only @ occupied_data
        check_equal(occupied_run, pointer_embedding(k % length, 0, occupied_data, length))

    # For lengths divisible by six the child pair returns to its arrival code
    # precisely when the router reopens. Check all ticks of two full periods.
    recurrence = {}
    for period_clock in (6, 12):
        repeating_tick = clocked_tick(period_clock, scatter, initial_dynamics)
        state = pointer_embedding(0, 0, seeds, period_clock)
        cycle_phase = (-1) ** (period_clock // 6)
        for period_index in range(2):
            phase = cycle_phase ** period_index
            for step in range(1, period_clock + 1):
                state = repeating_tick @ state
                local = np.linalg.matrix_power(local_pair, step - 1) @ split
                check_equal(state, pointer_embedding(step % period_clock, 1,
                                                    phase * local, period_clock))
            for step in range(1, period_clock + 1):
                state = repeating_tick @ state
                local = np.linalg.matrix_power(parent_only, step) @ seeds
                check_equal(state, pointer_embedding(step % period_clock, 0,
                                                    phase * local, period_clock))
        check_equal(repeating_tick @ state, pointer_embedding(1, 1, split, period_clock))
        recurrence[str(period_clock)] = {
            "projective_return_ticks": 2 * period_clock,
            "logical_code_phase_per_return": cycle_phase,
            "periods_verified": 2,
            "renewed_birth_ticks_verified": [1, 2 * period_clock + 1, 4 * period_clock + 1],
        }

    return {
        "clock_length": length,
        "fixed_step_dimension": tick.shape[0],
        "accepted_logical_code_dimension": 2,
        "returned_code_basis": ["Phi+", "-Psi-"],
        "source_scattering_rank": int(np.linalg.matrix_rank(source)),
        "closed_tick_and_scattering_unitary_within_1e-12": True,
        "local_pair_ticks_before_router_returns": length - 1,
        "orthogonal_code_orbit_columns": int(orbit_columns.shape[1]),
        "delayed_start_parent_only_ticks": 6,
        "delayed_start_birth_tick": 7,
        "rejected_phase_flipped_code_stays_dormant_until_next_window": True,
        "occupied_B_subspace_dimension_tested": occupied.shape[1],
        "pair_code_restored_after_six_local_ticks": True,
        "closed_split_reunion_recurrence": recurrence,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = check_router()
    if args.check:
        if json.loads(RESULT.read_text(encoding="utf-8")) != result:
            raise SystemExit("Saved answers differ from recomputation.")
        print("Verified", RESULT.name)
    else:
        RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
