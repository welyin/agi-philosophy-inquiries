"""Independent finite algebra and rational-budget audit for round 1259.

No author code imports.  This checks a finite code representation of the
continuous interval encoding; normal-state and mother-action claims belong
to the analytic proof. Default execution only compares saved results.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np


def maximum_error(value):
    return float(np.max(np.abs(value)))


def show(value):
    return float(f"{float(value):.13g}")


def ket(dimension, index):
    vector = np.zeros(dimension, complex)
    vector[index] = 1
    return vector


def projector(dimension, index):
    vector = ket(dimension, index)
    return np.outer(vector, vector.conj())


def group_code_and_source():
    # Axes: local K, original programme data, physical B flag, reference.
    nk, nd, nf, nr = 3, 2, 3, 6
    blank = 2
    angles = [math.pi / 2, math.pi / 3, 0.0]
    targets = [0, 1, None]
    generator = np.zeros((nk * nd * nf, nk * nd * nf), complex)
    amplitudes = np.zeros((nk, nf), float)
    for k in range(nk):
        rotation = np.zeros((nf, nf), complex)
        if targets[k] is None:
            amplitudes[k, blank] = 1
        else:
            a = targets[k]
            rotation[a, blank] = 1j * angles[k]
            rotation[blank, a] = -1j * angles[k]
            amplitudes[k, a] = math.sin(angles[k])
            amplitudes[k, blank] = math.cos(angles[k])
        generator += np.kron(projector(nk, k), np.kron(np.eye(nd), rotation))
    assert maximum_error(generator - generator.conj().T) < 1e-14
    eigenvalues, basis = np.linalg.eigh(generator)
    unitary = (basis * np.exp(-1j * eigenvalues)) @ basis.conj().T

    embedding = np.zeros((nk * nd * nf, nk * nd), complex)
    ideal = np.zeros_like(embedding)
    for k in range(nk):
        for d in range(nd):
            column = k * nd + d
            embedding[(k * nd + d) * nf + blank, column] = 1
            for a in range(nf):
                ideal[(k * nd + d) * nf + a, column] = amplitudes[k, a]
    actual = unitary @ embedding
    map_error = maximum_error(actual - ideal)
    isometry_error = maximum_error(actual.conj().T @ actual - np.eye(nk * nd))
    assert map_error < 2e-14 and isometry_error < 2e-14

    # Nonuniform Schmidt coefficients with complex phases: data/K entangle R.
    coefficients = np.array(
        [complex(j + 1, (-1) ** j * (j + 2)) for j in range(nk * nd)]
    )
    coefficients /= np.linalg.norm(coefficients)
    logical_reference = np.zeros(nk * nd * nr, complex)
    for j, coefficient in enumerate(coefficients):
        logical_reference[j * nr + j] = coefficient
    before = np.kron(embedding, np.eye(nr)) @ logical_reference
    after = np.kron(actual, np.eye(nr)) @ logical_reference
    expected = np.kron(ideal, np.eye(nr)) @ logical_reference
    reference_map_error = maximum_error(after - expected)
    assert reference_map_error < 2e-14

    def data_reference_density(state):
        matrix = state.reshape(nk, nd, nf, nr).transpose(1, 3, 0, 2)
        matrix = matrix.reshape(nd * nr, nk * nf)
        return matrix @ matrix.conj().T

    data_reference_error = maximum_error(
        data_reference_density(after) - data_reference_density(before)
    )
    assert data_reference_error < 2e-14

    physical_data_values = np.array([1.0, 2.0])
    shifts = np.array([4.0, 8.0, 0.0])
    physical_B = np.kron(
        np.eye(nk),
        np.diag([b + shift for b in physical_data_values for shift in shifts]),
    )
    flag_shift = np.kron(np.eye(nk * nd), np.diag(shifts))
    decoded_B = physical_B - flag_shift
    decoded_error = maximum_error(
        unitary.conj().T @ decoded_B @ unitary - decoded_B
    )
    assert decoded_error < 2e-13
    bare_operator_change = np.linalg.norm(
        unitary.conj().T @ physical_B @ unitary - physical_B, ord=2
    )
    assert bare_operator_change > 3
    B_lift = np.kron(physical_B, np.eye(nr))
    bare_mean_shift = float(np.vdot(after, B_lift @ after).real
                            - np.vdot(before, B_lift @ before).real)
    assert bare_mean_shift > 1

    tensor = after.reshape(nk, nd, nf, nr)
    actual_probabilities = np.sum(abs(tensor) ** 2, axis=(0, 1, 3))
    weights = abs(coefficients.reshape(nk, nd)) ** 2
    expected_probabilities = np.sum(
        weights[:, :, None] * amplitudes[:, None, :] ** 2, axis=(0, 1)
    )
    assert maximum_error(actual_probabilities - expected_probabilities) < 2e-14
    # Drop only the blank-ready precondition, using the same unitary/effect.
    nonblank = np.kron(ket(nk, 0), np.kron(ket(nd, 0), ket(nf, 0)))
    bad_output = (unitary @ nonblank).reshape(nk, nd, nf)
    bad_p_result0 = float(np.sum(abs(bad_output[:, :, 0]) ** 2))
    bad_p_blank = float(np.sum(abs(bad_output[:, :, blank]) ** 2))
    assert bad_p_result0 < 1e-28 and bad_p_blank > 1 - 2e-14

    # Abstract noncommuting spectator/source matrices test the exact algebra.
    # They are not a numerical construction of the Einstein stress tensor.
    other_dimension = 2
    identity = np.eye(nk * nd * nf)
    root = np.array([[0.94, 0.01], [0.01, 0.96]], complex)
    U_internal = np.diag([0.835, 0.830, 0.835])
    U_hold = np.kron(np.kron(U_internal, np.eye(nd * nf)), np.eye(2))
    mass_root = np.kron(identity, root)
    J_hold = mass_root - U_hold
    G_full = np.kron(generator, np.eye(other_dimension))
    U_full = np.kron(unitary, np.eye(other_dimension))
    epsilon = float(F(135005697, 40000000000000))
    pulse_peak = 280
    interaction = epsilon * pulse_peak * G_full
    J_during = J_hold - interaction
    complete_mass = J_during + U_hold + interaction
    source_identity_error = maximum_error(complete_mass - mass_root)
    assert source_identity_error < 2e-15
    assert float(np.linalg.eigvalsh(J_during).min()) > 0.01

    sigma_x = np.array([[0, 1], [1, 0]], complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], complex)
    sigma_z = np.diag([1, -1]).astype(complex)
    source_menu = [mass_root, J_hold]
    source_menu += [np.kron(identity, b) for b in (sigma_x, sigma_y, sigma_z)]
    source_error = max(
        maximum_error(U_full.conj().T @ b @ U_full - b) for b in source_menu
    )
    assert source_error < 2e-13
    omitted_J_residual = np.linalg.norm(
        J_hold + U_hold + interaction - mass_root, ord=2
    )
    interaction_norm = np.linalg.norm(interaction, ord=2)
    assert abs(omitted_J_residual - interaction_norm) < 1e-14
    assert omitted_J_residual > 0.001

    return {
        "status": "PASS",
        "dimensions_K_data_flag_reference": [nk, nd, nf, nr],
        "unitary_error": show(maximum_error(unitary.conj().T @ unitary - np.eye(18))),
        "whole_ready_code_map_error": show(map_error),
        "isometry_error": show(isometry_error),
        "entangled_reference_map_error": show(reference_map_error),
        "retained_data_reference_density_error": show(data_reference_error),
        "decoded_B_operator_error": show(decoded_error),
        "bare_B_operator_change": show(bare_operator_change),
        "bare_B_mean_shift": show(bare_mean_shift),
        "three_flag_probabilities": [show(x) for x in actual_probabilities],
        "blank_deleted_result0_probability": show(bad_p_result0),
        "blank_deleted_failure_probability": show(bad_p_blank),
        "complete_mass_identity_error": show(source_identity_error),
        "source_commutant_transport_error": show(source_error),
        "omit_J_compensation_mass_residual": show(omitted_J_residual),
        "limitations": [
            "Finite code algebra, not all normal L2(B) states or a spatial mother model.",
            "The other-sector matrices illustrate commutant identities, not Einstein-source data.",
            "Decoded original programme data are preserved; bare B is intentionally changed.",
            "Nonselective source preservation does not imply every conditional source state is unchanged.",
        ],
    }


def evaluate_polynomial(coefficients, x):
    return sum((coefficient * x ** power for power, coefficient in enumerate(coefficients)), F(0))


def differentiate(coefficients):
    return [F(i) * coefficients[i] for i in range(1, len(coefficients))]


def group_exact_resource_window():
    n, M, s = F(1), 513, 200000
    volume = F(1000, M ** 3)
    hbar = F(1)
    epsilon = hbar / (s ** 2 * n * volume)
    assert epsilon == F(135005697, 40000000000000)
    gate_time, read_time = F(1, 128), F(1, 128)
    processing = gate_time + read_time
    assert processing == F(1, 64)
    polynomial = [F(0), F(0), F(0), F(1), F(-3), F(3), F(-1)]
    integral = sum((value / (power + 1) for power, value in enumerate(polynomial)), F(0))
    assert 140 * integral == 1
    derivatives = polynomial
    endpoint_values = []
    for order in range(3):
        pair = [evaluate_polynomial(derivatives, F(0)),
                evaluate_polynomial(derivatives, F(1))]
        assert pair == [0, 0]
        endpoint_values.append([str(x) for x in pair])
        derivatives = differentiate(derivatives)
    peak_polynomial = evaluate_polynomial(polynomial, F(1, 2))
    assert peak_polynomial == F(1, 64)
    maximum_g = 140 * peak_polynomial / gate_time
    assert maximum_g == 280
    # The analytic input pi<22/7 gives ||G||<=pi/2<11/7.
    interaction_upper = epsilon * maximum_g * F(11, 7)
    assert interaction_upper == 440 * epsilon
    assert interaction_upper < F(3, 2000)
    J_input = F(1, 100)
    J_remainder_lower = J_input - interaction_upper
    assert J_remainder_lower > F(17, 2000)
    # Necessary action for a blank-to-orthogonal branch is pi/2;
    # the declared amplitude/time allowance is not automatically arbitrary.
    maximum_allowed_action = J_input * gate_time / epsilon
    assert maximum_allowed_action > F(11, 7)
    travel_upper = F(35, 342)
    access_upper = travel_upper + processing
    assert access_upper == F(1291, 10944)
    slack = F(1, 8) - access_upper
    assert slack == F(77, 10944) > 0
    trigger_min, trigger_max = F(65, 32), F(67, 32)
    end_upper = trigger_max + access_upper
    assert F(2) < trigger_min
    assert end_upper < F(71, 32) < F(9, 4)
    return {
        "status": "PASS",
        "hbar_fixed": str(hbar),
        "physical_similarity_s": s,
        "cell_volume": str(volume),
        "effective_epsilon": str(epsilon),
        "pulse_integral": str(140 * integral),
        "pulse_zero_extension_derivatives_0_to_2": endpoint_values,
        "gate_time": str(gate_time),
        "flag_access_time": str(read_time),
        "processing_time": str(processing),
        "physical_gate_time": str(s * gate_time),
        "pulse_maximum_g": str(maximum_g),
        "interaction_strict_upper": str(interaction_upper),
        "declared_comparison_upper": str(F(3, 2000)),
        "J_input_lower": str(J_input),
        "J_remaining_strict_lower": str(J_remainder_lower),
        "J_remaining_comparison": str(F(17, 2000)),
        "allowed_action_upper_budget_over_epsilon": str(maximum_allowed_action),
        "access_upper": str(access_upper),
        "access_slack_to_one_eighth": str(slack),
        "latest_end_strict_upper": str(end_upper),
        "limitations": [
            "The original finite-terminal geometry/transport allowance is inherited conditionally.",
            "No 513^3 grid, propagation channel, global clock or quantum geometry has been constructed.",
            "J_input is a declared quantum baseline-window condition, not the old classical margin.",
            "The held-block bound does not certify all local J_i under a sum of spatial roots.",
            "The physical similarity changes the complete material/resource model; hbar is not tuned.",
        ],
    }


def calculate():
    return {
        "round": 1259,
        "status": "PASS_INDEPENDENT_1259",
        "independent_groups": 2,
        "groups": {
            "finite_code_reference_and_source": group_code_and_source(),
            "exact_pulse_resource_and_window": group_exact_resource_window(),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = calculate()
    target = Path(__file__).with_name("independent_results.json")
    if args.write:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        action = "FIRST_WRITE"
    else:
        assert json.loads(target.read_text(encoding="utf-8")) == result
        action = "READ_ONLY_MATCH"
    print(json.dumps({"status": result["status"], "action": action,
                      "independent_groups": result["independent_groups"]}))


if __name__ == "__main__":
    main()
