"""Round 1272 finite thermal protocol. Default: verify saved JSON, no writes.

The 15360-dimensional carrier is represented by state columns and permutation
actions. No dense 15360 x 15360 operator or density matrix is constructed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DIMS = (2, 2, 2, 2, 4, 5, 2, 3, 2, 2, 2)
LABELS = ("A", "B", "S_A", "S_B", "Z_A", "Z_B", "K", "L", "M", "E", "J_out")
DIM = int(np.prod(DIMS))
COORD = np.indices(DIMS, dtype=np.int16).reshape(len(DIMS), DIM)
ENERGY = np.sum(COORD[:4], axis=0)
DELTA = 1.0e-9


def make_gate(name, theta, swap=None, flip=None, controls=()):
    active = np.ones(DIM, dtype=bool)
    for axis, values in controls:
        active &= np.isin(COORD[axis], values)
    changed = COORD.copy()
    if swap is not None:
        a, b = swap
        changed[a, active] = COORD[b, active]
        changed[b, active] = COORD[a, active]
    else:
        changed[flip, active] = 1 - COORD[flip, active]
    permutation = np.ravel_multi_index(tuple(changed), DIMS)
    assert np.array_equal(permutation[permutation], np.arange(DIM))
    assert np.array_equal(ENERGY[permutation], ENERGY)
    assert np.array_equal(active[permutation], active)
    return dict(name=name, theta=float(theta), active=active, permutation=permutation)


GATES = (
    make_gate("prepare_A", np.pi / 2, swap=(0, 2), controls=((4, (0,)),)),
    make_gate("prepare_B", np.pi / 2, swap=(1, 3), controls=((5, (0, 1)),)),
    make_gate("controlled_half_swap", np.pi / 4, swap=(0, 1), controls=((6, (1,)),)),
    make_gate("A_to_M", np.pi / 2, flip=8, controls=((7, (1, 2)), (0, (1,)))),
    make_gate("M_to_E", np.pi / 2, flip=9, controls=((7, (1, 2)), (8, (1,)))),
    make_gate("undo_M_to_E", np.pi / 2, flip=9, controls=((7, (2,)), (8, (1,)))),
    make_gate("undo_A_to_M", np.pi / 2, flip=8, controls=((7, (2,)), (0, (1,)))),
    make_gate("half_swap", np.pi / 4, swap=(0, 1)),
    make_gate("A_to_output", np.pi / 2, flip=10, controls=((0, (1,)),)),
)


def apply_gate(state, gate, error=0.0):
    """exp[-i (1+error) theta (I-P)], with inactive sectors exactly fixed."""
    out = state.copy()
    active = gate["active"]
    angle = (1.0 + error) * gate["theta"]
    out[active] = np.exp(-1j * angle) * (
        np.cos(angle) * state[active]
        + 1j * np.sin(angle) * state[gate["permutation"][active]]
    )
    return out


def ready(k, ell):
    state = np.zeros((DIM, 1), dtype=complex)
    for za in range(4):
        for zb in range(5):
            index = np.ravel_multi_index((0, 0, 1, 1, za, zb, k, ell, 0, 0, 0), DIMS)
            state[index, 0] = 1 / np.sqrt(20)
    return state


def run(state, errors=None):
    errors = np.zeros(9) if errors is None else errors
    snapshots = [state.copy()]
    for gate, error in zip(GATES, errors):
        state = apply_gate(state, gate, error)
        snapshots.append(state.copy())
    # Every pulse commutes with the physical free energy. Keep its actual
    # proper-clock phase, rather than resetting it between gate windows.
    free_phase = np.exp(-1j * (11 + ENERGY) * (17 / 8))[:, None]
    return free_phase * state, snapshots


def mean_axis(state, axis):
    return float(np.sum(COORD[axis, :, None] * np.abs(state) ** 2).real)


def reduced(state, keep):
    """Pure state on the carrier and a retained reference-column factor."""
    all_dims = DIMS + (state.shape[1],)
    rest = tuple(i for i in range(len(all_dims)) if i not in keep)
    matrix = state.reshape(all_dims).transpose(tuple(keep) + rest)
    rows = int(np.prod([all_dims[i] for i in keep]))
    matrix = matrix.reshape(rows, -1)
    return matrix @ matrix.conj().T


def trace_distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a - b))) / 2)


def entropy(rho):
    values = np.linalg.eigvalsh(rho)
    values = values[values > 1e-14]
    return float(-np.sum(values * np.log(values)))


def rat(value):
    value = F(value)
    return f"{value.numerator}/{value.denominator}"


def compute():
    gate_checks = []
    for gate in GATES:
        perm, active = gate["permutation"], gate["active"]
        gate_checks.append({
            "name": gate["name"],
            "involution_all_basis": bool(np.array_equal(perm[perm], np.arange(DIM))),
            "energy_commutator_zero_all_basis": bool(np.array_equal(ENERGY[perm], ENERGY)),
            "active_invariant_all_basis": bool(np.array_equal(active[perm], active)),
            "generator_norm_upper_pi_units": round(2 * gate["theta"] / np.pi, 12),
        })

    branches, columns, actual_columns = [], [], []
    errors = DELTA * np.array((1, -1, 1, -1, 1, -1, 1, -1, 1))
    support_leak = 0.0
    for k in range(2):
        for ell in range(3):
            initial = ready(k, ell)
            output, snapshots = run(initial)
            actual, actual_snapshots = run(initial, errors)
            expected = F(13, 40) if k == 0 or ell == 1 else F(2, 5)
            probability = mean_axis(output, 10)
            assert abs(probability - float(expected)) < 2e-14
            assert abs(mean_axis(snapshots[2], 0) - 1 / 4) < 2e-14
            assert abs(mean_axis(snapshots[2], 1) - 2 / 5) < 2e-14
            for item in snapshots + actual_snapshots:
                support_leak = max(support_leak, float(np.sum(np.abs(item[ENERGY != 2]) ** 2)))
                assert abs(np.vdot(item, item).real - 1) < 2e-14
            record_after_undo = mean_axis(snapshots[7], 8) + mean_axis(snapshots[7], 9)
            if ell in (0, 2):
                assert record_after_undo < 2e-28
            branches.append({
                "K": k, "L": ell,
                "prepared_qA": mean_axis(snapshots[2], 0),
                "prepared_qB": mean_axis(snapshots[2], 1),
                "first_exchange_delta_A": mean_axis(snapshots[3], 0) - 1 / 4,
                "first_exchange_delta_B": mean_axis(snapshots[3], 1) - 2 / 5,
                "final_probability": probability,
                "exact_final_probability": rat(expected),
                "actual_amplitude_test_probability": mean_axis(actual, 10),
                "M_plus_E_after_step7": record_after_undo,
                "full_pure_norm": float(np.vdot(output, output).real),
            })
            columns.append(output[:, 0])
            actual_columns.append(actual[:, 0])

    isometry = np.column_stack(columns)
    actual_isometry = np.column_stack(actual_columns)
    gram_error = float(np.max(np.abs(isometry.conj().T @ isometry - np.eye(6))))
    actual_gram_error = float(np.max(np.abs(actual_isometry.conj().T @ actual_isometry - np.eye(6))))
    operator_error = float(np.linalg.norm(actual_isometry - isometry, ord=2))
    assert gram_error < 3e-14 and actual_gram_error < 3e-14
    assert operator_error < 9 * np.pi * DELTA
    assert support_leak == 0.0

    _, prep_snapshots = run(ready(1, 0))
    rho_a = reduced(prep_snapshots[2], (0,))
    rho_b = reduced(prep_snapshots[2], (1,))
    rho_ab = reduced(prep_snapshots[2], (0, 1))
    preparation_product_error = float(np.max(np.abs(rho_ab - np.kron(rho_a, rho_b))))
    assert preparation_product_error < 3e-14
    assert abs(rho_a[1, 1].real - 1 / 4) < 2e-14
    assert abs(rho_b[1, 1].real - 2 / 5) < 2e-14

    bell_checks = []
    for ell in range(3):
        bell = np.column_stack((ready(0, ell)[:, 0], ready(1, ell)[:, 0])) / np.sqrt(2)
        output, _ = run(bell)
        initial_kr = reduced(bell, (6, 11))
        final_kr = reduced(output, (6, 11))
        r_error = float(np.max(np.abs(reduced(output, (11,)) - np.eye(2) / 2)))
        # Complete J_out instrument: retain both subnormalised conditional
        # states and sum them only for the unconditioned reference statement.
        conditional_r = np.zeros((2, 2), dtype=complex)
        branch_probabilities = []
        for bit in (0, 1):
            branch = output.copy()
            branch[COORD[10] != bit] = 0
            conditional_r += reduced(branch, (11,))
            branch_probabilities.append(float(np.vdot(branch, branch).real))
        assert r_error < 3e-14
        assert abs(sum(branch_probabilities) - 1) < 3e-14
        assert np.max(np.abs(conditional_r - np.eye(2) / 2)) < 3e-14
        bell_checks.append({
            "L": ell,
            "reference_only_error": r_error,
            "KR_trace_distance_from_input": trace_distance(final_kr, initial_kr),
            "output_probabilities": branch_probabilities,
            "reference_after_complete_instrument_error": float(np.max(np.abs(conditional_r - np.eye(2) / 2))),
        })

    start = F(1, 2) - F(1, 256)
    latest_end = F(1, 2) + 9 * F(1, 16) + F(1, 256)
    capacity = F(21, 400)
    sample_start = F(2) - F(1, 100000000)
    output_time = F(17, 8)
    q_integral = 140 * (F(1, 4) - 3 * F(1, 5) + 3 * F(1, 6) - F(1, 7))
    assert q_integral == 1
    assert start - capacity > 0 and sample_start - latest_end > 0
    assert output_time - latest_end - capacity > 0
    error_upper = (1 + 9 * F(22, 7)) * F(1, 10**9)
    coarse_error = F(3, 10**8)
    witness_lower = F(3, 40) - 2 * coarse_error
    assert error_upper < coarse_error and witness_lower > 0

    return {
        "status": "PASS_AUTHOR_PROTOCOL_WITH_SCOPE",
        "scope": "New finite comoving effective task; not a Jacobson identity or gravity-selection proof.",
        "dimensions": list(DIMS), "labels": list(LABELS), "total_dimension": DIM,
        "method": "state columns and permutation actions; no large dense matrix",
        "gate_checks": gate_checks, "branches": branches,
        "isometry_gram_error": gram_error,
        "actual_isometry_gram_error": actual_gram_error,
        "actual_test_operator_error": operator_error,
        "excitation_support_leak": support_leak,
        "preparation": {
            "target_AB_product_error": preparation_product_error,
            "beta_A_times_gap": float(np.log(3)),
            "beta_B_times_gap": float(np.log(3 / 2)),
            "target_entropy_A": entropy(rho_a), "target_entropy_B": entropy(rho_b),
            "full_preparation_is_pure_and_not_joint_Gibbs": True,
        },
        "bell_reference_checks": bell_checks,
        "exact_certificates": {
            "pulse_integral": rat(q_integral), "pulse_peak": "35/1",
            "earliest_start_minus_capacity": rat(start - capacity),
            "latest_completion": rat(latest_end),
            "geometry_sample_gap": rat(sample_start - latest_end),
            "output_gap_after_capacity": rat(output_time - latest_end - capacity),
            "total_half_channel_error_upper_using_pi_lt_22_over_7": rat(error_upper),
            "coarse_half_channel_error": rat(coarse_error),
            "record_effect_probability_gap_lower": rat(witness_lower),
            "initial_excitation": "2/1", "positive_baseline_units": "11/1",
        },
        "not_certified": [
            "all ideal TR operations", "intermediate external fine ports",
            "unchanged KR marginal", "external or irreversible record erasure",
            "Unruh temperature or Wald thermal identity", "new heat-driven geometry",
        ],
    }


def equal_saved(actual, saved, location="root"):
    if isinstance(actual, dict):
        assert set(actual) == set(saved), location
        for key in actual:
            equal_saved(actual[key], saved[key], f"{location}.{key}")
    elif isinstance(actual, list):
        assert len(actual) == len(saved), location
        for index, (a, b) in enumerate(zip(actual, saved)):
            equal_saved(a, b, f"{location}[{index}]")
    elif isinstance(actual, float):
        assert np.isclose(actual, saved, rtol=1e-11, atol=5e-14), location
    else:
        assert actual == saved, location


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="exclusively create first results JSON")
    args = parser.parse_args()
    result = compute()
    destination = HERE / "protocol_results.json"
    if args.write:
        with destination.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        equal_saved(result, json.loads(destination.read_text(encoding="utf-8")))
    print(json.dumps({"status": result["status"], "dimension": DIM,
                      "excitation_support_leak": result["excitation_support_leak"],
                      "isometry_gram_error": result["isometry_gram_error"],
                      "reference_cases": len(result["bell_reference_checks"]),
                      "mode": "first-save" if args.write else "read-only-replay"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
