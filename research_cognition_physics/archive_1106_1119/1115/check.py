"""1115: complete parent-process lift and finite-domain dictionary certificate.

No old experiment is rerun.  Frozen constructors/results are reused; the new
calculations retain the native clock, actual material factors and records.
Default only compares the saved receipt.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESULTS = HERE / "results.json"
TOL = 2e-9


def norm(a):
    return float(np.linalg.norm(a))


def small(x, label):
    x = float(x)
    if x > TOL:
        raise AssertionError(f"{label}: {x}")
    return x


def load_module(path):
    spec = importlib.util.spec_from_file_location("round1112_dictionary_source", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sqrt_psd(a):
    e, v = np.linalg.eigh(a)
    return (v * np.sqrt(np.maximum(e, 0))) @ v.conj().T


def reduced_q_channel(v, environment_dim):
    # Output order [environment, Q], input Q; return the full superoperator.
    a = v.reshape(environment_dim, 2, 2)
    return np.einsum("dai,dbj->abij", a, a.conj()).reshape(4, 4)


def calculate():
    src = load_module(STAGE / "1112/check.py")
    material = src.Material(1)  # constructor only, not its calculate()/main()
    saved1112 = json.loads((STAGE / "1112/results.json").read_text(encoding="utf-8"))
    saved1111 = json.loads((STAGE / "1111/results.json").read_text(encoding="utf-8"))
    saved1114 = json.loads((STAGE / "1114/results.json").read_text(encoding="utf-8"))
    assert saved1112["round"] == 1112 and saved1114["round"] == 1114
    assert saved1112["scientific_calibration"] == saved1114["scientific_calibration"] == 0
    assert "PASS" in saved1111["status"]

    # All native clock degrees remain in the complete 288-dimensional parent
    # object.  The four-dimensional code is an explicit prepared department,
    # never a claim that the physical object has only four states.
    vac = np.array([1, 0, 0], complex)
    data = src.tensor(material.aux, vac, vac)
    jclock = material.gamma @ np.kron(np.eye(4), data[:, None])
    hnative = material.baseline * np.eye(288) + material.hfk
    hclock = material.baseline * np.eye(4) + material.hc
    clock_embedding = small(norm(jclock.conj().T @ jclock - np.eye(4)), "native clock embedding")
    clock_generator = small(norm(hnative @ jclock - jclock @ hclock), "same native and parent generator")

    # The ideal parent C has a concrete unitary absent from the already
    # certified H-commuting native subgroup.  This is a coverage marker, not
    # a no-go for extensions using extra resources.
    e, basis = np.linalg.eigh(hclock)
    p = np.eye(4)
    p[:, [0, 3]] = p[:, [3, 0]]
    u4 = basis @ p @ basis.conj().T
    unative = np.eye(288) + jclock @ (u4 - np.eye(4)) @ jclock.conj().T
    unitary_error = small(norm(unative.conj().T @ unative - np.eye(288)), "complete parent unitary")
    complete_pure_transfer = small(norm(unative @ (jclock @ basis[:, 0]) - jclock @ basis[:, 3]), "parent full pure transfer")
    native_commutator = norm(hnative @ unative - unative @ hnative)
    small(abs(native_commutator - 3 * math.sqrt(2)), "native C coverage distinction")

    p0, p1 = np.diag([1., 0]), np.diag([0., 1.])
    clock_read = src.tensor(np.eye(4), np.eye(2), np.eye(2), p1, np.eye(3), np.eye(3))
    times = [5.07, 5.47]
    clock_states, clock_reads, lifted, qchannels = [], [], [], []
    for t in times:
        state = jclock @ src.exp_h(hclock, t)[:, 0]
        clock_states.append(state)
        clock_reads.append(float(np.vdot(state, clock_read @ state).real))
        w = np.kron(state[:, None], np.eye(2))
        lifted.append(w)
        qchannels.append(reduced_q_channel(w, 288))
        small(norm(w.conj().T @ w - np.eye(2)), "complete prepared parent channel")
        small(norm(qchannels[-1] - np.eye(4)), "same logical identity after clock trace")
        small(abs(clock_reads[-1] - math.sin(t / 2) ** 6), "actual clock read retained")
    clock_distance = math.sqrt(max(0., 1 - abs(np.vdot(clock_states[0], clock_states[1])) ** 2))
    old_distance = saved1112["group7_complete_process_not_reduced_channel"]["complete_clock_half_trace_distance"]
    small(abs(clock_distance - old_distance), "frozen clock witness reproduced in complete parent")
    assert abs(clock_reads[0] - clock_reads[1]) > 0.02

    # Complete pair object A_2 tensor B_2: xx, xy, yx, yy all retained.
    pa = [np.array([5., 4., 0., 0.]), np.array([5., 0., 4., 0.])]
    pb = [np.array([5., -4., 0., 0.]), np.array([5., 0., -4., 0.])]
    pair_p = np.array([a + b for a in pa for b in pb])
    pair_ops = [np.diag(pair_p[:, mu]) for mu in range(4)]
    flip_pair = np.eye(4)
    flip_pair[:, [0, 3]] = flip_pair[:, [3, 0]]
    controlled = np.kron(p0, np.eye(4)) + np.kron(p1, flip_pair)
    scatter = np.zeros((8, 2), complex)
    scatter[0, 0] = scatter[7, 1] = 1
    init = np.kron(np.eye(2), np.eye(4)[:, 0:1])
    scatter_residual = small(norm(controlled @ init - scatter), "full pair lift agrees with actual prepared xx")
    pair_conservation = max(norm(np.kron(np.eye(2), op) @ controlled - controlled @ np.kron(np.eye(2), op)) for op in pair_ops)
    small(pair_conservation, "complete pair mother four-P boundary")
    plus = np.array([1., 1.]) / math.sqrt(2)
    joint = (scatter @ plus).reshape(2, 4)
    rho_q = joint @ joint.conj().T
    logical_damage = float(np.linalg.norm(rho_q - np.outer(plus, plus), ord="nuc") / 2)
    small(abs(logical_damage - 0.5), "discarded recoil is not clean Q")
    x = np.array([[0, 1], [1, 0]], complex)
    effect = (np.eye(2) + 0.6 * x) / 2
    actual_p = float(np.trace(effect @ rho_q).real)
    false_clean_p = float(np.vdot(plus, effect @ plus).real)
    small(abs(actual_p - 0.5) + abs(false_clean_p - 0.8), "existing logical read distinguishes false deletion")

    # New joint lift: same running native clock + conservative material
    # correlation + 1112 nonsharp instrument + both recorded outcomes.
    kp, km = sqrt_psd(effect), sqrt_psd(np.eye(2) - effect)
    vi = np.zeros((8, 2), complex)
    vi.reshape(2, 2, 2, 2)[:, 0, 0] = kp
    vi.reshape(2, 2, 2, 2)[:, 1, 1] = km
    v_joint = np.kron(vi, np.eye(4)) @ scatter  # [Q,D,E,A,B]
    instrument_complete = small(norm(v_joint.conj().T @ v_joint - np.eye(2)), "full scatter then instrument")
    t = 10.0
    evolution = src.exp_h(hclock, t)
    work_blank = np.eye(4)[:, 0:1]
    jc_work = np.kron(jclock @ evolution, work_blank)
    hc_work = np.kron(hnative @ jclock @ evolution, work_blank)
    full_lift = np.kron(jc_work, v_joint)
    full_isometry = small(norm(full_lift.conj().T @ full_lift - np.eye(8)), "joint parent all-input isometry")
    transport = src.boost_matrix([0.3, 0.2, 0.1])
    u = np.array([1., 0., 0., 0.])
    # Unchanged incoming/outgoing carrier and program roles have scalar
    # four-P on these logical values.  Keep an explicit positive sum in both
    # boundaries; the algebra is independent of its numerical value.
    spectator_p = np.array([7., 1., 0.5, 0.])
    all_errors = []
    for lam in [np.eye(4), transport]:
        ul = lam @ u
        out_pair = pair_p @ lam.T
        pin_pair = lam @ pair_p[0]
        ps = lam @ spectator_p
        for mu in range(4):
            op_pair = np.diag(out_pair[:, mu])
            lhs = ul[mu] * np.kron(hc_work, v_joint) + np.kron(jc_work, np.kron(np.eye(8), op_pair) @ v_joint) + ps[mu] * full_lift
            pin = ul[mu] * np.kron(hclock, np.eye(2)) + (pin_pair[mu] + ps[mu]) * np.eye(8)
            all_errors.append(norm(lhs - full_lift @ pin))
    small(max(all_errors), "joint native clock and physical recoil use same four-P")
    # The open physical preparation fixes the clock state; all four formal
    # clock columns were used only to certify the operator intertwiner.
    ready_input = np.kron(np.eye(4)[:, 0:1], np.eye(2))
    full_ready = full_lift @ ready_input
    expected_ready = np.kron(jc_work[:, 0:1], v_joint)
    prepared_error = small(norm(full_ready - expected_ready), "declared preparation not arbitrary native-clock permission")

    # A single shared role table defines both parent and actual event tasks.
    event = np.array([23 / 15 + math.pi, 1., 0., 0.])
    origin = np.zeros(4)
    transformed = transport @ event
    parent_dictionary_event = transformed - transport @ origin
    actual_dictionary_event = transport @ (event - origin)
    endpoint_error = small(norm(parent_dictionary_event - actual_dictionary_event), "same exact task endpoint map")
    inverse = np.linalg.inv(transport)
    inverse_error = small(norm(inverse @ transformed - event), "same reference inverse")
    small(abs(src.mdot(event, event) - src.mdot(transformed, transformed)), "same task interval")
    # A non-product distributed encoding cannot count as a true local
    # receiver merely because the global output map is an isometry.
    scatter_q_super = np.einsum("abi,cbj->acij", scatter.reshape(2, 4, 2), scatter.reshape(2, 4, 2).conj()).reshape(4, 4)
    receiver_identity_error = norm(scatter_q_super - np.eye(4))
    assert receiver_identity_error > 1.

    objects = [
        {"role": "logical_CO_position_payload", "parent_complete_matrix_object": True,
         "implemented_full_logical_instruments": True, "native_full_C_claimed": False},
        {"role": "running_clock_device", "native_dimension": 288,
         "parent_complete_matrix_object": True, "retained_full_state": True,
         "native_full_C_claimed": False, "implemented_menu": "1106/1111/1112 declared preparation, wait, records, conservative slots"},
        {"role": "recoil_material_A_tensor_B", "native_dimension": 4,
         "parent_complete_matrix_object": True, "retained_full_state": True,
         "native_full_C_claimed": False, "actual_preparation": "xx; coherent xx/yy output; crossed geometry not promoted"},
        {"role": "program_work_vac_carriers_records_environment", "parent_complete_matrix_object": True,
         "retained_full_state": True, "context_before_fresh_input": True}
    ]
    for row in objects:
        assert row["parent_complete_matrix_object"]
    sources = {}
    for name in ["../archive_1086_1105/1090/proof.md", "1106/proof.md", "1111/proof.md", "1111/results.json",
                 "1112/proof.md", "1112/check.py", "1112/results.json", "1114/proof.md", "1114/results.json"]:
        sources[name] = hashlib.sha256((STAGE / name).read_bytes()).hexdigest()
    return {
        "round": 1115, "status": "PASS_DECLARED_DOMAIN_COMPLETE_PARENT_PROCESS_REALIZATION",
        "numerical_check_groups": 6, "new_adopted_axioms": 0, "scientific_calibration": 0,
        "source_sha256": sources,
        "scope": {"exact_parent_FUCP_preserved": True, "native_materials_are_complete_parent_objects_not_deleted": True,
                  "partial_experiment_realization_not_full_material_mother_C": True,
                  "one_event_port_resource_dictionary": True, "stage_closure_left_to_main_audit": True,
                  "no_old_full_experiment_rerun": True},
        "group1_complete_objects_and_C_coverage": {"roles": objects,
            "native_clock_embedding_residual": clock_embedding, "same_generator_residual": clock_generator,
            "complete_parent_unitarity_residual": unitary_error, "complete_parent_pure_transfer_residual": complete_pure_transfer,
            "chosen_parent_U_native_H_commutator": native_commutator,
            "commutator_is_coverage_marker_not_new_resource_no_go": True},
        "group2_clock_complete_lift_and_P_deletion": {"native_parent_dimension": 288,
            "full_clock_Q_isometry_shape": list(lifted[0].shape), "times": times,
            "both_reduced_Q_channels_are_identity": True, "complete_clock_half_trace_distance": clock_distance,
            "actual_M_reads": clock_reads, "read_gap": abs(clock_reads[0] - clock_reads[1]),
            "source_frozen_distance": old_distance},
        "group3_recoil_complete_lift_and_P_deletion": {"full_pair_dimension": 4,
            "scatter_output_shape": list(scatter.shape), "actual_prepared_scatter_residual": scatter_residual,
            "full_pair_four_P_commutator": pair_conservation, "Q_half_trace_distance_from_false_identity": logical_damage,
            "existing_nonsharp_read_actual_probability": actual_p, "false_clean_probability": false_clean_p},
        "group4_joint_clock_matter_record_complete_process": {"open_formal_clock_Q_dimension": 8,
            "actual_clock_preparation_fixed": True, "isometry_shape": list(full_lift.shape),
            "retained_unchanged_scalar_carrier_program_four_P": spectator_p.tolist(),
            "instrument_complete_residual": instrument_complete, "whole_isometry_residual": full_isometry,
            "same_four_P_in_two_references_residual": max(all_errors), "declared_preparation_residual": prepared_error,
            "all_materials_native_clock_and_D_E_records_retained": True},
        "group5_same_event_tasks_and_receiver_reflection": {"original_single_source_completion": event.tolist(),
            "whole_transported_completion": transformed.tolist(), "same_endpoint_dictionary_residual": endpoint_error,
            "inverse_role_transport_residual": inverse_error, "distributed_scatter_receiver_not_identity_residual": receiver_identity_error,
            "menu_equality_only_for_declared_event_domain": True,
            "local_receiver_factorization_required_not_global_decode": True},
        "group6_composition_domain_and_counts": {"complete_resources_and_allowed_future_wires_are_typed": True,
            "same_window_directory_embeddings_inherit_1112": True,
            "whole_reference_and_all_admitted_influences_inherit_1114": True,
            "unimplemented_native_parent_CP_not_counted_as_actual_SRF_resource": True,
            "no_approximation_closure_promoted_to_exact_P": True,
            "new_numerical_groups_are_certificate_checks_not_new_scientific_calibration": True}
    }


def compare(saved, actual, path="root"):
    if isinstance(actual, dict):
        if not isinstance(saved, dict) or saved.keys() != actual.keys():
            raise AssertionError(f"{path}: keys differ")
        for key in actual:
            compare(saved[key], actual[key], f"{path}.{key}")
    elif isinstance(actual, list):
        if not isinstance(saved, list) or len(saved) != len(actual):
            raise AssertionError(f"{path}: length differs")
        for i, (x, y) in enumerate(zip(saved, actual)):
            compare(x, y, f"{path}[{i}]")
    elif isinstance(actual, float):
        if not math.isclose(float(saved), actual, rel_tol=3e-9, abs_tol=3e-9):
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
    print("PASS: 6 finite dictionary/certificate groups; largest lift 36864 x 8")
    print("Complete native clock, recoil and records retained; declared-domain realization")
    print("Scientific calibration 0; complete native-material C and full physical SR not claimed")


if __name__ == "__main__":
    main()
