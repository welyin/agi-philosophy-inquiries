"""Verify the physical exchange on the qubit actually passed into the original birth detector."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[2] / "猜想" / "知行闭环的分化与物质生成猜想" / "闭环分化验证"
sys.path.insert(0, str(MODEL))
from quantum_loop_inheritance import I, X, Z, ZERO, controlled_flip
from quantum_split_retry import Y, clock_blocks
from living_ledger_succession import handover_columns, quantum_handover


PAULI = np.array([X, Y, Z])
I4 = np.eye(4)
SWAP = I4.reshape(2, 2, 4).swapaxes(0, 1).reshape(4, 4)
BLANK = np.kron(I, ZERO[:, None])
Z_AXIS = np.array([0., 0., 1.])


def check_equal(left, right):
    np.testing.assert_allclose(left, right, rtol=0, atol=3e-12)


def sigma(coefficients):
    return np.tensordot(coefficients, PAULI, axes=1)


def birth_effects():
    rotate = np.array([[1., -1.], [1., 1.]]) / np.sqrt(2)
    loop = controlled_flip("S") @ np.kron(rotate, I) @ controlled_flip("R")
    early, _, idle = clock_blocks(loop)
    source = np.eye(64)[:, [0, 8]]
    effects = []
    for steps in range(6):
        path = early @ np.linalg.matrix_power(idle, steps) @ source
        success = path.reshape(4, 16, 2)[2]
        effect = success.T @ success
        effects.append(effect)
    for effect, signed_axis in zip(effects, (Z, -X, 0 * I, -Z, X, 0 * I)):
        check_equal(effect, (I + signed_axis) / 2)
    return effects


def source_operator(n, theta, gap):
    h = gap * (I - sigma(n)) / 2
    propagator = np.cos(theta) * I4 - 1j * np.sin(theta) * SWAP
    old, new = np.kron(h, I), np.kron(I, h)
    check_equal(propagator.conj().T @ (old + new) @ propagator, old + new)
    # Free h evolution commutes with both the exchange and the measured energy.
    return BLANK.T @ (propagator.conj().T @ new @ propagator - new) @ BLANK


def rounded(values):
    return np.asarray(values).round(12).tolist()


def calculate():
    inputs = json.loads((HERE / "inputs.json").read_bytes())
    gap = inputs["physical_inputs"]["gap"]
    check_equal(handover_columns(2), I)
    original_corrected_branches = quantum_handover()
    effects = birth_effects()
    hermitian_coordinates = lambda op: [np.trace(item @ op).real for item in [I, *PAULI]]
    rank = int(np.linalg.matrix_rank([hermitian_coordinates(op) for op in [I, *effects]]))
    assert rank == 3
    y_inputs = [(I + sign * Y) / 2 for sign in [1, -1]]
    probabilities = [[np.trace(rho @ effect).real for effect in effects] for rho in y_inputs]
    check_equal(probabilities, np.full((2, 6), .5))
    identities = 0
    for vector in inputs["probe_axes"]:
        n = np.array(vector) / np.linalg.norm(vector)
        for theta in np.linspace(0, np.pi / 2, 37):
            c, s = np.cos(theta), np.sin(theta)
            coeff = gap * (-s*s * n + c*s * np.cross(n, Z_AXIS)) / 2
            constant = gap * s*s * n[2] / 2
            whole = source_operator(n, theta, gap)
            visible = constant * I + coeff[0] * X + coeff[2] * Z
            check_equal(whole, visible + coeff[1] * Y)
            check_equal(np.linalg.norm(whole - visible, ord=2), abs(coeff[1]))
            gains = [np.trace(rho @ whole).real for rho in y_inputs]
            check_equal(abs(gains[0] - gains[1]) / 2, abs(coeff[1]))
            identities += 1
    example = {}
    xaxis = np.array([1., 0., 0.])
    for name, theta in [("initial", 0.), ("half_exchange", np.pi/4), ("finished", np.pi/2)]:
        op = source_operator(xaxis, theta, gap)
        example[name] = rounded([np.trace(rho @ op).real for rho in y_inputs])
    # Linear window weight in calibrated exchange angle. By parts the integral
    # equals endpoint gain minus the unweighted average of the gain.
    nodes, weights = np.polynomial.legendre.leggauss(20)
    weighted = source_operator(xaxis, np.pi/2, gap)
    for node, weight in zip(nodes, weights):
        weighted -= weight * source_operator(xaxis, (node + 1) * np.pi / 4, gap) / 2
    moment = [np.trace(rho @ weighted).real for rho in y_inputs]
    check_equal(moment, [gap/(2*np.pi), -gap/(2*np.pi)])
    # The positive single-run Z record and a distinct mixed-axis ensemble task.
    z_final = source_operator(Z_AXIS, np.pi / 2, gap)
    check_equal(z_final, gap * (I - effects[0]))
    diagonal_axis = np.array([1., 0., 1.]) / np.sqrt(2)
    test_states = [(I + sign * sigma(diagonal_axis)) / 2 for sign in (1, -1)]
    mixed = [op / 2 for effect in (effects[0], effects[4]) for op in (effect, I - effect)]
    optimal_success = sum(max(np.trace(rho @ op).real for rho in test_states) / 2 for op in mixed)
    check_equal(optimal_success, (1 + 1/np.sqrt(2))/2)
    rank_with_phase = int(np.linalg.matrix_rank([hermitian_coordinates(op) for op in [I, *effects, Y]]))
    p = inputs["thermal_test"]["old_excitation"]
    old_reference = np.diag([1-p, p])
    heat = np.trace(z_final @ old_reference).real
    entropy = -p*np.log(p) - (1-p)*np.log1p(-p)
    beta = np.log((1-p)/p) / gap
    excess = entropy - beta*heat
    check_equal(excess, -np.log1p(-p))
    return {"round": 1818, "original_successor_branches_verified": original_corrected_branches,
            "six_window_signed_axes": ["Z","-X","0","-Z","X","0"],
            "real_accessible_dimension": rank, "with_named_Y_reader_dimension": rank_with_phase,
            "birth_probabilities_for_plus_and_minus_Y": rounded(probabilities),
            "same_real_X_energy_in_gap_units": example,
            "linear_weight_flux_plus_minus_Y_in_gap_units": rounded(moment),
            "thermal_Z_example": {"gain_over_gap": rounded(heat/gap), "beta_times_gap": rounded(beta*gap),
                                  "receiver_entropy_change_nats": rounded(entropy),
                                  "entropy_minus_beta_Q_nats": rounded(excess)},
            "operator_identity_test_cases": identities,
            "sharp_Z_event_from_window_1": True,
            "random_XZ_record_for_diagonal_axis_max_success": rounded(optimal_success)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    path = HERE / "results.json"
    if args.check:
        if result != json.loads(path.read_bytes()):
            raise SystemExit("Stored results differ from the independent gate calculation")
        print("PASS 1818: all six original birth windows, four handover branches, 185 source tests")
    else:
        path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf8")
