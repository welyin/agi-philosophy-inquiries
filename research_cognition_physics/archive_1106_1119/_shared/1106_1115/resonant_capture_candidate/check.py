"""Read-only-by-default audit of one typed resonant-capture candidate.

This does not construct a normal wave packet, spatial capture interaction, or
new research round. --write-results explicitly saves the candidate audit.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.json"
TOL = 2e-11


def opnorm(matrix: np.ndarray) -> float:
    return float(np.linalg.norm(matrix, 2))


def require_small(value: float, label: str, tolerance: float = TOL) -> None:
    if not value <= tolerance:
        raise AssertionError(f"{label}: {value} > {tolerance}")


def mdot(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[0] - a[1] * b[1])


def exponential(h: np.ndarray, t: float) -> np.ndarray:
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * t * values)) @ vectors.conj().T


def tensor(*factors: np.ndarray) -> np.ndarray:
    answer = np.array([1.0 + 0.0j])
    for factor in factors:
        answer = np.kron(answer, factor)
    return answer


def receiver_image(w: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    # Physical factor order: history, logical material, D, M, Q_C, Q_r.
    kraus = w.reshape(64, 2, 2)
    return sum(k @ matrix @ k.conj().T for k in kraus)


def calculate() -> dict:
    omega = 1.0
    length = 1.0
    capture_delta = math.pi / 4
    completion_proper = math.pi
    window = 0.02
    resident_base = 10.0
    launch_base = 5.0

    uc = np.array([5 / 4, 3 / 4])
    ua = np.array([5 / 4, -3 / 4])
    q = [np.array([5 / 2, (n - 1) * 3 / 2]) for n in range(3)]
    momenta = np.array([
        (2 + c) * uc + (2 + a) * ua + q[n]
        for c in range(2) for a in range(2) for n in range(3)
    ])
    # a=|01,R1>, b=|10,R0>, in the original 12-dimensional material space.
    pair_embedding = np.zeros((12, 2), complex)
    pair_embedding[4, 0] = 1
    pair_embedding[6, 1] = 1
    k = momenta[4]
    pair_error = float(np.max(abs(momenta[6] - k)))
    require_small(pair_error, "equal incoming total four-momentum")
    invariant_mass = math.sqrt(mdot(k, k))
    ustar = k / invariant_mass
    vstar = ustar[1] / ustar[0]
    to_rest = ustar[0] * np.array([[1.0, -vstar], [-vstar, 1.0]])
    ucr = to_rest @ uc
    uar = to_rest @ ua
    qr = [to_rest @ p for p in q[:2]]
    arrival = 23 * length / 15
    starts = [0.0, 152 * length / 83, 152 * length / 83,
              152 * length / 175]
    velocities = [ucr[1] / ucr[0], uar[1] / uar[0],
                  qr[0][1] / qr[0][0], qr[1][1] / qr[1][0]]
    arrival_residual = max(abs(x + v * arrival - length)
                           for x, v in zip(starts, velocities))
    require_small(arrival_residual, "typed incoming trajectories meet")
    expected_uc = np.array([23, 15]) / (4 * math.sqrt(19))
    expected_ua = np.array([83, -45]) / (16 * math.sqrt(19))
    require_small(float(np.max(abs(ucr - expected_uc))), "carrier velocity")
    require_small(float(np.max(abs(uar - expected_ua))), "clock A velocity")
    source = np.array([0.0, 0.0])
    capture_begin = np.array([arrival, length])
    capture_end = capture_begin + np.array([capture_delta, 0.0])
    completed = capture_begin + np.array([completion_proper, 0.0])
    comparison_boost = 1.25 * np.array([[1.0, -0.6], [-0.6, 1.0]])
    transformed_events = [comparison_boost @ x for x in
                          [source, capture_begin, capture_end, completed]]
    proper_invariance_error = abs(
        mdot(transformed_events[3] - transformed_events[1],
             transformed_events[3] - transformed_events[1])
        - completion_proper ** 2)
    require_small(proper_invariance_error, "transported processing proper time")
    assert all(abs(v) < 1 for v in velocities)
    assert mdot(completed - source, completed - source) > 0

    i2 = np.eye(2, dtype=complex)
    i32 = np.eye(32, dtype=complex)
    z = np.array([1, 0], complex)
    one = np.array([0, 1], complex)
    p0, p1 = np.outer(z, z), np.outer(one, one)
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    b = -y
    btilde = np.zeros((12, 12), complex)
    for n in range(2):
        btilde[3 + n + 1, 6 + n] = 1j
        btilde[6 + n, 3 + n + 1] = -1j
    material_restriction_error = opnorm(
        pair_embedding.conj().T @ btilde @ pair_embedding - b)
    require_small(material_restriction_error, "1109 pair restriction")
    selected_p_errors = []
    for mu in range(2):
        selected_p_errors.append(opnorm(
            np.diag(momenta[:, mu]) @ pair_embedding -
            pair_embedding * k[mu]))
    require_small(max(selected_p_errors), "selected full material P")

    swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                     [0, 1, 0, 0], [0, 0, 0, 1]], complex)
    ap, am = (i2 + b) / 2, (i2 - b) / 2
    ur = np.kron(ap, i2) - 1j * np.kron(am, y)
    um = np.kron(p0, i2) - 1j * np.kron(p1, x)
    g1 = np.kron(np.eye(8), swap)
    g2 = np.kron(ur, np.eye(8))
    g3 = np.kron(np.kron(i2, um), np.eye(4))
    gates = [g1, g2, g3]
    gate_error = max(opnorm(g.conj().T @ g - i32) for g in gates)
    require_small(gate_error, "given gates unitary")

    # Hc and HFK are specified directly, BEFORE constructing Gamma.
    hc = 1.5 * omega * np.eye(4)
    hfk = 1.5 * omega * np.eye(128, dtype=complex)
    for j, g in enumerate(gates):
        coupling = omega * math.sqrt((j + 1) * (3 - j)) / 2
        hc[j + 1, j] = hc[j, j + 1] = coupling
        hfk[(j + 1) * 32:(j + 2) * 32, j * 32:(j + 1) * 32] = coupling * g
        hfk[j * 32:(j + 1) * 32, (j + 1) * 32:(j + 2) * 32] = coupling * g.conj().T
    hin_internal = np.kron(hc, i32)
    # The resident launch clock is running before arrival; prepare it so that
    # it reaches head 0 at the predeclared arrival event, independently of Q.
    head0 = np.array([1, 0, 0, 0], complex)
    launch_prepared = exponential(hc, -arrival) @ head0
    launch_arrival_error = float(np.linalg.norm(
        exponential(hc, arrival) @ launch_prepared - head0))
    require_small(launch_arrival_error, "predeployed running launch clock")
    gamma = np.zeros((128, 128), complex)
    products = [i32]
    for j in range(1, 4):
        products.append(gates[j - 1] @ products[-1])
    for j, product in enumerate(products):
        gamma[j * 32:(j + 1) * 32, j * 32:(j + 1) * 32] = product
    h_intertwiner = opnorm(hfk @ gamma - gamma @ hin_internal)
    require_small(h_intertwiner, "independently specified Hamiltonians intertwine")
    ready = np.zeros((128, 32), complex)
    ready[:32] = i32
    require_small(opnorm(gamma @ ready - ready), "ready subspace identity")
    full_endpoint_target = np.zeros_like(ready)
    full_endpoint_target[96:] = products[-1]
    full_endpoint_error = opnorm(
        exponential(hfk, completion_proper) @ ready + full_endpoint_target)
    require_small(full_endpoint_error, "all-data autonomous endpoint")
    eigval, eigvec = np.linalg.eigh(hfk)
    spectral_errors = []
    for level in range(4):
        block = eigvec[:, abs(eigval - level * omega) < 1e-8]
        compression = ready.conj().T @ block @ block.conj().T @ ready
        spectral_errors.append(opnorm(
            compression - math.comb(3, level) / 8 * i32))
    require_small(max(spectral_errors), "full prepared spectral measure")

    baseline = resident_base + launch_base + invariant_mass
    hin = baseline * np.eye(128) + hin_internal
    hout = baseline * np.eye(128) + hfk
    cap = gamma @ exponential(hin, capture_delta)
    capture_intertwiner = opnorm(hout @ cap - cap @ hin)
    require_small(capture_intertwiner, "positive-time capture energy")
    frame_p_errors = []
    # Common rest frame and its one specified finite boost.
    for frame_u in [np.array([1.0, 0.0]), comparison_boost @ np.array([1.0, 0.0])]:
        for mu in range(2):
            frame_p_errors.append(opnorm(
                frame_u[mu] * hout @ cap - cap @ (frame_u[mu] * hin)))
    require_small(max(frame_p_errors), "capture in two frames, complete P")
    end = exponential(hout, completion_proper - capture_delta) @ cap @ ready
    expected_end = -np.exp(-1j * baseline * completion_proper) * full_endpoint_target
    positive_window_endpoint_error = opnorm(end - expected_end)
    require_small(positive_window_endpoint_error, "capture plus autonomous interval")

    # Known clock B=-1; arbitrary fresh Q remains a two-column operator input.
    known_clock = np.array([1, 1j]) / math.sqrt(2)
    data_input = np.column_stack([
        tensor(known_clock, z, z, qstate, z) for qstate in [z, one]
    ])
    input_q = ready @ data_input
    effect_m = tensor(np.eye(4), i2, i2, p1, i2, i2)
    commutator_m = opnorm(hout @ effect_m - effect_m @ hout)
    assert commutator_m > 0.1
    times = [capture_delta, math.pi / 2, completion_proper]
    record_rows = []
    energy_errors = []
    reference_energy = input_q.conj().T @ hin @ input_q
    reference_energy2 = input_q.conj().T @ (hin @ hin) @ input_q
    for tau in times:
        w = exponential(hout, tau - capture_delta) @ cap @ input_q
        actual = w.conj().T @ effect_m @ w
        probability = math.sin(omega * tau / 2) ** 6
        residual = opnorm(actual - probability * i2)
        require_small(residual, "nonconstant internal record, all Q")
        energy_errors.append(opnorm(w.conj().T @ hout @ w - reference_energy))
        energy_errors.append(opnorm(w.conj().T @ (hout @ hout) @ w - reference_energy2))
        record_rows.append({"proper_time": tau, "M_probability": probability,
                            "all_Q_effect_residual": residual})
    require_small(max(energy_errors), "energy and energy-square despite changing record")
    wfinal = exponential(hout, completion_proper - capture_delta) @ cap @ input_q
    direct_received = np.column_stack([
        tensor(known_clock, one, one, z, qstate) for qstate in [z, one]
    ])
    direct_target = np.zeros_like(input_q)
    direct_target[96:] = direct_received
    # Hc transfer contributes -1; U_m contributes -i.
    direct_target *= 1j * np.exp(-1j * baseline * completion_proper)
    true_receiver_error = opnorm(wfinal - direct_target)
    require_small(true_receiver_error, "actual receiver factor, unknown Q/Z")

    window_rows = []
    basis_matrices = [np.outer(a, b.conj()) for a in [z, one] for b in [z, one]]
    for offset in [-window, window]:
        tau = completion_proper + offset
        w = exponential(hout, tau - capture_delta) @ cap @ input_q
        p_terminal = math.cos(omega * offset / 2) ** 6
        overlap = wfinal.conj().T @ w
        scalar = np.trace(overlap) / 2
        require_small(opnorm(overlap - scalar * i2), "complete window overlap scalar")
        require_small(abs(abs(scalar) ** 2 - p_terminal), "complete window fidelity")
        halfdiamond = math.sqrt(max(0.0, 1 - p_terminal))
        bound = math.sqrt(3) * omega * abs(offset) / 2
        assert halfdiamond <= bound + TOL
        return_probability = math.sin(omega * offset / 2) ** 6
        receiver_residual = 0.0
        for rho in basis_matrices:
            expected = (1 - return_probability) * rho + return_probability * np.trace(rho) * p0
            receiver_residual = max(receiver_residual,
                                    opnorm(receiver_image(w, rho) - expected))
        require_small(receiver_residual, "receiver reduced channel in finite window")
        window_rows.append({"offset": offset,
                            "complete_output_halfdiamond": halfdiamond,
                            "complete_output_bound": bound,
                            "M_missing_probability": 1 - p_terminal,
                            "receiver_halfdiamond": return_probability,
                            "receiver_channel_residual": receiver_residual})

    # This distinguishes algebraic compatibility from a physically supplied capture law.
    bare_ready_energy_defect = opnorm((hout - hin) @ input_q)
    require_small(abs(bare_ready_energy_defect - math.sqrt(1.5) * omega),
                  "bare same-basis loading is not an energy intertwiner")
    wcap = cap @ input_q
    barecap = exponential(hin, capture_delta) @ input_q
    cap_m = wcap.conj().T @ effect_m @ wcap
    bare_m = barecap.conj().T @ effect_m @ barecap
    partial_record_probability = math.sin(omega * capture_delta / 2) ** 6
    require_small(opnorm(cap_m - partial_record_probability * i2), "capture already includes prefix")
    require_small(opnorm(bare_m), "mere loading does not create the prefix record")
    assert partial_record_probability > 0

    return {
        "status": "PASS_TYPED_CANDIDATE_NOT_ACTUAL_CAPTURE",
        "formal_round": None,
        "new_adopted_axioms": 0,
        "scientific_calibration": 0,
        "numerical_check_groups": 7,
        "parameters": {"Omega": omega, "L": length, "capture_delta": capture_delta,
                       "completion_proper_time": completion_proper, "window": window,
                       "resident_baseline_mass": resident_base, "launch_baseline_mass": launch_base},
        "group1_typed_geometry": {
            "incoming_pair_total_P_original": k.tolist(), "invariant_mass": invariant_mass,
            "u_star_original": ustar.tolist(), "u_C_rest": ucr.tolist(), "u_A_rest": uar.tolist(),
            "initial_positions_C_A_R0_R1": starts, "arrival_time": arrival,
            "source_capturebegin_captureend_complete": [p.tolist() for p in [source, capture_begin, capture_end, completed]],
            "same_events_second_frame": [p.tolist() for p in transformed_events],
            "arrival_residual": arrival_residual, "proper_invariance_residual": proper_invariance_error,
            "all_external_legs_future_timelike": True,
            "exact_sharp_positions_and_momenta_not_claimed": True},
        "group2_material_sector": {
            "selected_original_basis_indices": [4, 6], "logical_B_restriction_residual": material_restriction_error,
            "full_P_restriction_residual": max(selected_p_errors),
            "full_1109_product_input_menu_not_claimed": True,
            "payload_is_original_C_degenerate_factor": True},
        "group3_independent_H": {
            "dimension": 128, "gate_unitarity_residual": gate_error,
            "H_defined_before_Gamma": True, "intertwiner_residual": h_intertwiner,
            "endpoint_all_data_residual": full_endpoint_error,
            "launch_spectrum": [0.0, omega, 2 * omega, 3 * omega],
            "launch_probabilities": [1 / 8, 3 / 8, 3 / 8, 1 / 8],
            "launch_mean": 1.5 * omega, "launch_variance": 0.75 * omega ** 2,
            "launch_preparation_independent_of_Q": True,
            "running_launch_ready_at_arrival_residual": launch_arrival_error,
            "complete_prepared_spectral_residual": max(spectral_errors),
            "positive_mass_range": [baseline, baseline + 3 * omega]},
        "group4_positive_capture": {
            "capture_energy_intertwiner_residual": capture_intertwiner,
            "two_frame_P_intertwiner_residual": max(frame_p_errors),
            "capture_plus_execution_endpoint_residual": positive_window_endpoint_error,
            "Vcap_is_explicit_material_input_not_derived_collision": True},
        "group5_clock_and_true_receiver": {
            "record_rows": record_rows, "record_commutator_norm": commutator_m,
            "energy_and_energy_square_residual": max(energy_errors),
            "source_to_resident_receiver_operator_residual": true_receiver_error,
            "all_unknown_Q_and_passive_reference_by_operator_identity": True,
            "no_freeflight_clock_difference_claim": True},
        "group6_finite_completion_window": {
            "rows": window_rows, "complete_processor_is_retained": True,
            "permanent_release_not_claimed": True},
        "group7_admission_boundary": {
            "bare_ready_capture_energy_defect": bare_ready_energy_defect,
            "positive_capture_already_writes_M_probability": partial_record_probability,
            "bare_load_M_probability": 0.0,
            "Gamma_ready_identity_does_not_make_capture_free": True,
            "finite_space_causal_capture_response_missing": True,
            "normal_wavepacket_or_alternative_actual_port_contract_missing": True,
            "full_translation_representation_not_established_by_label_transport": True,
            "continue_only_if_capture_response_independently_supplied_or_verified": True,
            "do_not_mark_1106_complete": True}
    }


def compare(saved, actual, path="root"):
    if isinstance(actual, dict):
        if not isinstance(saved, dict) or saved.keys() != actual.keys():
            raise AssertionError(f"{path}: object keys differ")
        for key in actual:
            compare(saved[key], actual[key], f"{path}.{key}")
    elif isinstance(actual, list):
        if not isinstance(saved, list) or len(saved) != len(actual):
            raise AssertionError(f"{path}: list shape differs")
        for i, (old, new) in enumerate(zip(saved, actual)):
            compare(old, new, f"{path}[{i}]")
    elif isinstance(actual, float):
        if not isinstance(saved, (float, int)) or not math.isclose(saved, actual, rel_tol=2e-10, abs_tol=2e-10):
            raise AssertionError(f"{path}: {saved!r} != {actual!r}")
    elif saved != actual:
        raise AssertionError(f"{path}: {saved!r} != {actual!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    actual = calculate()
    if args.write_results:
        RESULTS.write_text(json.dumps(actual, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print("WROTE results.json")
    else:
        compare(json.loads(RESULTS.read_text(encoding="utf-8")), actual)
        print("SAVED_RESULTS_MATCH")
    print("PASS: 7 numerical groups; typed capture candidate only")
    print("Formal round: none; scientific calibration: 0; 1106 completion: false")
    print("Actual capture region, source/control response and port realization remain unverified")


if __name__ == "__main__":
    main()
