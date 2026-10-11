"""Round 1106: finite effective rail/capture compilation, read-only by default.

Only 288-dimensional local dense matrices and an 82944-by-4 full-network
isometry are used. The rail and local capture laws are declared material
inputs; this is not a microscopic collision construction or full PT1 proof.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.json"
TOL = 4e-10


def norm(a):
    return float(np.linalg.norm(a, 2))


def small(a, label, tol=TOL):
    value = float(a)
    if value > tol:
        raise AssertionError(f"{label}: {value}")
    return value


def tensor(*args):
    result = np.array([1.0 + 0.0j])
    for arg in args:
        result = np.kron(result, arg)
    return result


def mdot(a, b):
    return float(a[0] * b[0] - a[1] * b[1])


def exp_generator(h, t):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * e)) @ v.conj().T


def exp_action(e, v, t, a):
    return v @ (np.exp(-1j * t * e)[:, None] * (v.conj().T @ a))


def network(w0, w1):
    # Input order: source port 0 with qubit 0/1, then source port D.
    return np.column_stack([np.kron(w0[:, q], w1[:, 0]) for q in (1, 2)]
                           + [np.kron(w0[:, 0], w1[:, q]) for q in (1, 2)])


def calculate():
    omega, length, spacing = 1.0, 1.0, 10.0
    delta, completion = math.pi / 4, math.pi
    m_resident, m_launch = 10.0, 5.0
    arrival = 23 * length / 15
    uc = np.array([5 / 4, 3 / 4])
    ua = np.array([5 / 4, -3 / 4])
    q = [np.array([5 / 2, (n - 1) * 3 / 2]) for n in range(3)]
    momenta_original = np.array([(2 + c) * uc + (2 + a) * ua + q[n]
                                for c in range(2) for a in range(2) for n in range(3)])
    k = momenta_original[4]
    small(np.max(abs(k - momenta_original[6])), "resonant pair")
    mass = math.sqrt(mdot(k, k))
    vstar = k[1] / k[0]
    rest = k[0] / mass * np.array([[1, -vstar], [-vstar, 1]])
    ucr, uar = rest @ uc, rest @ ua
    qr0, qr1 = rest @ q[0], rest @ q[1]
    p_a = [2 * ucr, 3 * uar, qr1]
    p_b = [3 * ucr, 2 * uar, qr0]
    starts_a = [0.0, 152 * length / 83, 152 * length / 175]
    starts_b = [0.0, 152 * length / 83, 152 * length / 83]
    actions, chis, arrival_errors = [], [], []
    boost = 1.25 * np.array([[1, -0.6], [-0.6, 1]])
    action_boost_errors = []
    for material, starts in [(p_a, starts_a), (p_b, starts_b)]:
        small(np.max(abs(np.sum(material, axis=0) - np.array([mass, 0]))), "rest total P")
        action, chi = 0.0, 0.0
        for p, initial_x in zip(material, starts):
            displacement = np.array([arrival, length - initial_x])
            arrival_errors.append(abs(initial_x + p[1] / p[0] * arrival - length))
            assert mdot(p, p) > 0
            assert mdot(displacement, displacement) > 0
            action += mdot(p, displacement)
            chi += p[1] * initial_x
            action_boost_errors.append(abs(mdot(boost @ p, boost @ displacement) - mdot(p, displacement)))
        actions.append(action)
        chis.append(chi)
    small(max(arrival_errors), "all incoming rails arrive")
    small(max(action_boost_errors), "whole-material action under boost")
    phase_identity = np.max(abs(np.asarray(actions) - (mass * arrival + np.asarray(chis))))
    small(phase_identity, "same energy and port-basis action")
    bin_phase = np.diag(np.exp(-1j * np.asarray(chis)))
    bout_phase = np.eye(2)
    physical_free = np.exp(-1j * mass * arrival) * np.eye(2)
    rail_free = np.diag(np.exp(-1j * np.asarray(actions)))
    same_free_identity = norm(bout_phase.conj().T @ physical_free @ bin_phase - rail_free)
    small(same_free_identity, "free rail compilation identity")
    event_rows = []
    for origin in [0.0, spacing]:
        events = [np.array([0.0, origin]), np.array([arrival, origin + length]),
                  np.array([arrival + delta, origin + length]),
                  np.array([arrival + completion, origin + length])]
        for first, second in zip(events, events[1:]):
            assert second[0] > first[0] and mdot(second - first, second - first) > 0
        transformed = [boost @ e for e in events]
        small(abs(mdot(transformed[-1] - transformed[1], transformed[-1] - transformed[1])
                  - completion ** 2), "same receiver proper duration")
        event_rows.append({"origin": origin, "events_o_a_b_f": [x.tolist() for x in events],
                           "boosted_same_events": [x.tolist() for x in transformed]})
    assert spacing > max(starts_a + starts_b) + length
    assert abs(event_rows[0]["boosted_same_events"][0][0] - event_rows[1]["boosted_same_events"][0][0]) > 0

    i2, i3, i72 = np.eye(2, dtype=complex), np.eye(3, dtype=complex), np.eye(72, dtype=complex)
    z, one = np.array([1, 0], complex), np.array([0, 1], complex)
    vac = np.array([1, 0, 0], complex)
    qbasis = [i3[:, i] for i in range(3)]
    y = np.array([[0, -1j], [1j, 0]], complex)
    x = np.array([[0, 1], [1, 0]], complex)
    b = -y
    p0, p1 = np.outer(z, z), np.outer(one, one)
    ap, am = (i2 + b) / 2, (i2 - b) / 2
    ur = np.kron(ap, i2) - 1j * np.kron(am, y)
    um = np.kron(p0, i2) - 1j * np.kron(p1, x)
    swap3 = np.zeros((9, 9), complex)
    for a in range(3):
        for bindex in range(3):
            swap3[bindex * 3 + a, a * 3 + bindex] = 1
    gates = [np.kron(np.eye(8), swap3), np.kron(ur, np.eye(18)),
             np.kron(np.kron(i2, um), np.eye(9))]
    gate_error = max(norm(g.conj().T @ g - i72) for g in gates)
    small(gate_error, "given full gates")
    hc = 1.5 * omega * np.eye(4)
    hfk = 1.5 * omega * np.eye(288, dtype=complex)
    for j, gate in enumerate(gates):
        coupling = omega * math.sqrt((j + 1) * (3 - j)) / 2
        hc[j + 1, j] = hc[j, j + 1] = coupling
        hfk[(j + 1) * 72:(j + 2) * 72, j * 72:(j + 1) * 72] = coupling * gate
        hfk[j * 72:(j + 1) * 72, (j + 1) * 72:(j + 2) * 72] = coupling * gate.conj().T
    hin_internal = np.kron(hc, i72)
    products = [i72]
    gamma = np.zeros((288, 288), complex)
    for j in range(1, 4):
        products.append(gates[j - 1] @ products[-1])
    for j in range(4):
        gamma[j * 72:(j + 1) * 72, j * 72:(j + 1) * 72] = products[j]
    small(norm(hfk @ gamma - gamma @ hin_internal), "independent H intertwiner")
    ready = np.zeros((288, 72), complex)
    ready[:72] = i72
    efk, vfk = np.linalg.eigh(hfk)
    endpoint_all = np.zeros_like(ready)
    endpoint_all[216:] = products[-1]
    all_endpoint_error = norm(exp_action(efk, vfk, completion, ready) + endpoint_all)
    small(all_endpoint_error, "all 72-dimensional data endpoint")
    spectral_errors = []
    for level in range(4):
        block = vfk[:, abs(efk - level * omega) < 1e-8]
        spectral_errors.append(norm(ready.conj().T @ block @ block.conj().T @ ready
                                    - math.comb(3, level) / 8 * i72))
    small(max(spectral_errors), "ready spectral measure")
    baseline = m_resident + m_launch + mass
    hin = baseline * np.eye(288) + hin_internal
    hout = baseline * np.eye(288) + hfk
    eout = efk + baseline
    cap = gamma @ exp_generator(hin, delta)
    cap_errors = []
    for u in [np.array([1.0, 0.0]), boost @ np.array([1.0, 0.0])]:
        for mu in range(2):
            cap_errors.append(norm(u[mu] * hout @ cap - cap @ (u[mu] * hin)))
    small(max(cap_errors), "all material P intertwiner")

    clock = np.array([1, 1j]) / math.sqrt(2)
    head0 = np.array([1, 0, 0, 0], complex)
    launch_prepared = exp_generator(hc, -arrival) @ head0
    logical_prepared = np.diag(np.exp(1j * np.asarray(actions))) @ clock
    initial = np.column_stack([tensor(launch_prepared, logical_prepared, z, z, qv, vac)
                               for qv in qbasis])
    free_before = np.exp(-1j * (m_resident + m_launch) * arrival) * tensor(
        exp_generator(hc, arrival), rail_free, i2, i2, i3, i3)
    at_a = free_before @ initial
    arrival_phase = np.exp(-1j * (m_resident + m_launch) * arrival)
    declared_ready = np.column_stack([tensor(head0, clock, z, z, qv, vac) for qv in qbasis])
    arrival_preparation_error = norm(at_a - arrival_phase * declared_ready)
    small(arrival_preparation_error, "running predeployment including all rail phases")
    wcap = cap @ at_a

    def local_w(tau):
        if tau < delta:
            raise AssertionError("unspecified capture internal path")
        return exp_action(eout, vfk, tau - delta, wcap)

    wt = local_w(completion)
    head3 = np.array([0, 0, 0, 1], complex)
    local_phase = 1j * arrival_phase * np.exp(-1j * baseline * completion)
    aux = local_phase * tensor(head3, clock, one, one, vac)
    target_local = np.column_stack([np.kron(aux, qv) for qv in qbasis])
    local_endpoint_error = norm(wt - target_local)
    small(local_endpoint_error, "true resident receiver qutrit")
    m_effect = tensor(np.eye(4), i2, i2, p1, i3, i3)
    record_rows = []
    energy_errors = []
    ein = initial.conj().T @ hin @ initial
    ein2 = initial.conj().T @ hin @ hin @ initial
    for tau in [delta, math.pi / 2, completion]:
        w = local_w(tau)
        expected = math.sin(omega * tau / 2) ** 6
        residual = norm(w.conj().T @ m_effect @ w - expected * i3)
        small(residual, "M response all vac/Q")
        energy_errors.extend([norm(w.conj().T @ hout @ w - ein),
                              norm(w.conj().T @ hout @ hout @ w - ein2)])
        record_rows.append({"tau": tau, "M_probability": expected, "all_input_effect_residual": residual})
    small(max(energy_errors), "complete energy and squared energy")

    actual_network = network(wt, wt)
    target_network = network(target_local, target_local)
    network_endpoint_error = norm(actual_network - target_network)
    small(network_endpoint_error, "same original CO code full output compilation")
    network_iso_error = norm(actual_network.conj().T @ actual_network - np.eye(4))
    small(network_iso_error, "full source code arbitrary state/reference")
    # Reorder physical output [aux0,Qr0,aux1,Qr1] -> [aux0,aux1,Qr0,Qr1].
    ordered = actual_network.reshape(96, 3, 96, 3, 4).transpose(0, 2, 1, 3, 4).reshape(96 * 96, 9, 4)
    js_out = np.zeros((9, 4), complex)
    js_out[3, 0] = js_out[6, 1] = js_out[1, 2] = js_out[2, 3] = 1
    factor_target = np.kron(aux, aux)[:, None, None] * js_out[None, :, :]
    auxiliary_factor_error = norm((ordered - factor_target).reshape(-1, 4))
    small(auxiliary_factor_error, "same full auxiliary output for both occupied rails")
    coherent_input = np.array([1, 0, 1, 0], complex) / math.sqrt(2)
    coherent_image = np.einsum("arc,c->ar", ordered, coherent_input)
    receiver_state = coherent_image.T @ coherent_image.conj()
    coherent_target = js_out @ coherent_input
    coherent_receiver_error = norm(receiver_state - np.outer(coherent_target, coherent_target.conj()))
    small(coherent_receiver_error, "no which-rail decoherence")

    window_rows = []
    for e0, e1 in [(0.02, 0.02), (-0.02, 0.02), (0.02, 0.01)]:
        we = network(local_w(completion + e0), local_w(completion + e1))
        overlap = actual_network.conj().T @ we
        scalar = np.trace(overlap) / 4
        small(norm(overlap - scalar * np.eye(4)), "all-position/reference overlap scalar")
        f = math.cos(omega * e0 / 2) ** 6 * math.cos(omega * e1 / 2) ** 6
        small(abs(abs(scalar) ** 2 - f), "two-complete-cell fidelity")
        exact = math.sqrt(max(0.0, 1 - f))
        upper = math.sqrt(3) * omega / 2 * math.sqrt(e0 * e0 + e1 * e1)
        assert exact <= upper + TOL
        window_rows.append({"e0": e0, "e1": e1, "complete_halfdiamond": exact,
                            "analytic_upper_bound": upper,
                            "overlap_scalar_residual": norm(overlap - scalar * np.eye(4))})

    # Whole-network first and second P moments without a global dense matrix.
    hin_initial = hin @ initial
    hout_final = hout @ wt
    def sum_action(w, hw):
        pairs = [(1, 0), (2, 0), (0, 1), (0, 2)]
        return np.column_stack([np.kron(hw[:, a], w[:, b]) + np.kron(w[:, a], hw[:, b])
                                for a, b in pairs])
    win = network(initial, initial)
    hin_win = sum_action(initial, hin_initial)
    hout_wout = sum_action(wt, hout_final)
    total_first_error = norm(actual_network.conj().T @ hout_wout - win.conj().T @ hin_win)
    total_second_error = norm(hout_wout.conj().T @ hout_wout - hin_win.conj().T @ hin_win)
    small(total_first_error, "full network energy")
    small(total_second_error, "full network energy square")

    return {
        "status": "PASS_CONDITIONAL_EFFECTIVE_CO_CLOCK_CAPTURE_BRIDGE",
        "round": 1106, "new_adopted_axioms": 0, "scientific_calibration": 0,
        "numerical_check_groups": 8,
        "scope": {"effective_rail_and_local_capture_laws_declared": True,
                  "complete_PT1_or_all_SR_not_claimed": True,
                  "no_arbitrary_arrival_filter_required": True,
                  "fresh_single_source_and_joint_source_code_distinguished": True,
                  "one_finite_exact_model_not_universal_finite_processor": True},
        "parameters": {"Omega": omega, "L": length, "D": spacing, "delta": delta,
                       "T": completion, "m_resident": m_resident, "m_launch": m_launch},
        "group1_material_and_two_port_events": {
            "selected_total_P_original": k.tolist(), "incoming_invariant_mass": mass,
            "initial_positions_a": starts_a, "initial_positions_b": starts_b,
            "arrival_error": max(arrival_errors), "events": event_rows,
            "boost_action_error": max(action_boost_errors),
            "transported_source_events_not_resliced": True},
        "group2_same_energy_and_port_phase": {
            "actions_a_b": actions, "input_spatial_phases_a_b": chis,
            "action_identity_residual": float(phase_identity),
            "physical_to_rail_free_identity_residual": same_free_identity,
            "logical_vac_and_Q_have_equal_mechanical_material": True},
        "group3_independently_defined_local_material": {
            "local_dimension": 288, "data_dimension": 72,
            "gate_unitarity_residual": gate_error,
            "full_72_data_endpoint_residual": all_endpoint_error,
            "ready_spectral_residual": max(spectral_errors),
            "positive_mass_range": [baseline, baseline + 3 * omega],
            "launch_spectrum": [0.0, omega, 2 * omega, 3 * omega],
            "launch_probabilities": [1 / 8, 3 / 8, 3 / 8, 1 / 8]},
        "group4_predeployment_positive_capture": {
            "running_launch_and_all_rail_preparation_residual": arrival_preparation_error,
            "complete_four_P_intertwiner_residual": max(cap_errors),
            "delta_positive": delta > 0, "declared_local_capture_not_free_label": True,
            "source_Q_not_used_to_prepare_resources": True},
        "group5_clock_record_and_real_receiver": {
            "clock_record_rows": record_rows,
            "local_true_receiver_operator_residual": local_endpoint_error,
            "local_energy_and_square_residual": max(energy_errors),
            "no_unaccounted_midrun_measurement": True},
        "group6_original_CO_complete_compilation": {
            "whole_isometry_shape": list(actual_network.shape),
            "full_endpoint_residual": network_endpoint_error,
            "full_isometry_residual": network_iso_error,
            "same_full_auxiliary_factor_residual": auxiliary_factor_error,
            "position_coherent_receiver_residual": coherent_receiver_error,
            "arbitrary_passive_reference_by_operator_identity": True},
        "group7_whole_process_finite_window": {"cases": window_rows,
            "all_two_cell_auxiliary_outputs_retained": True,
            "permanent_decoupling_not_claimed": True},
        "group8_joint_resources_and_limit": {
            "whole_energy_operator_residual": total_first_error,
            "whole_energy_square_operator_residual": total_second_error,
            "vac_rails_not_free_absent_material": True,
            "full_local_mother_space_capture_intertwines": True,
            "no_82944_square_matrix_constructed": True,
            "complete_stage_success_not_claimed": True}
    }


def compare(saved, actual, path="root"):
    if isinstance(actual, dict):
        if not isinstance(saved, dict) or saved.keys() != actual.keys():
            raise AssertionError(f"{path}: keys differ")
        for k in actual:
            compare(saved[k], actual[k], f"{path}.{k}")
    elif isinstance(actual, list):
        if not isinstance(saved, list) or len(saved) != len(actual):
            raise AssertionError(f"{path}: length differs")
        for i, (a, b) in enumerate(zip(saved, actual)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(actual, float):
        if not math.isclose(saved, actual, rel_tol=3e-9, abs_tol=3e-9):
            raise AssertionError(f"{path}: {saved} != {actual}")
    elif saved != actual:
        raise AssertionError(f"{path}: {saved} != {actual}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_results:
        RESULTS.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print("WROTE results.json")
    else:
        compare(json.loads(RESULTS.read_text(encoding="utf-8")), result)
        print("SAVED_RESULTS_MATCH")
    print("PASS: 8 numerical groups; local 288; whole isometry 82944 x 4")
    print("Conditional effective CO/clock/capture bridge; scientific calibration: 0")
    print("Complete PT1, all influences/matter and whole SR completion: not claimed")


if __name__ == "__main__":
    main()
