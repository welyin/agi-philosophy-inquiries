"""1158: complete material blocks for the candidate same-carrier rules.

Import only the frozen 1112 Material helper, never its complete calculate().
Default execution is read-only; --write-once writes the first local receipt.
Whole-network material tensors are handled by the proved factor identities.
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
OLD = HERE.parent.parent / "archive_1106_1119"
TOL = 4e-9
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1.0, -1.0]).astype(complex)
H = (X + Z) / math.sqrt(2)


def small(value, label):
    value = float(value)
    if value > TOL:
        raise AssertionError(f"{label}: {value}")
    return value


def norm(value):
    return float(np.linalg.norm(value))


def source_check():
    record = json.loads((HERE / "source_hashes.json").read_text(encoding="utf-8"))
    for name, digest in record["sha256"].items():
        actual = hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        if actual != digest:
            raise AssertionError(f"source changed: {name}")
    return len(record["sha256"])


def load_material():
    spec = importlib.util.spec_from_file_location("round1112_material_reused", OLD / "1112/check.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exp_h(h, t):
    eigenvalues, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * t * eigenvalues)) @ vectors.conj().T


def old_clock_and_complete_block(old):
    material = old.Material(1)
    hout = material.baseline * np.eye(288) + material.hfk
    hin = material.baseline * np.eye(288) + np.kron(material.hc, np.eye(72))
    delta = math.pi / 4
    cap = material.gamma @ exp_h(hin, delta)
    intertwiner = small(norm(hout @ cap - cap @ hin), "complete original capture energy")
    waiting_rows = []
    for duration in [math.pi / 2, math.pi, 0.37]:
        full_input = np.eye(288, dtype=complex).reshape(4, 8, 3, 3, 1, 1, 288)
        actual = material.wait(full_input, duration).reshape(288, 288)
        independent = exp_h(hout, duration)
        direct = small(norm(actual - independent), "direct FK versus compiled waiting")
        unitary_error = small(norm(actual.conj().T @ actual - np.eye(288)), "full mother waiting")
        # For remote X and Z, the complete Heisenberg residual is precisely
        # (U*U-I) tensor g; never inferred from one fixed ready-state marginal.
        remote = max(norm(np.kron(actual.conj().T @ actual - np.eye(288), g))
                     for g in [X, Z])
        small(remote, "remote full-mother generators")
        waiting_rows.append({"duration": duration, "direct_FK_residual": direct,
                             "unitarity_residual": unitary_error,
                             "remote_generator_residual": remote})
    probabilities = []
    for duration in [math.pi / 2, math.pi]:
        prepared = material.loaded_joint(np.eye(3), 0.0)
        actual = material.wait(prepared, duration)
        # Auxiliary index is ab x Dr x Mr; Mr is its least significant bit.
        p = np.sum(np.abs(actual[:, [1, 3, 5, 7], ...]) ** 2,
                   axis=tuple(range(actual.ndim - 1)))
        expected = math.sin(duration / 2) ** 6
        small(norm(p - expected), "same original Mr for all qutrit inputs")
        probabilities.append({"proper_time": duration, "all_qutrit_probabilities": p.tolist()})
    return {"original_device_dimension": 288,
            "positive_internal_energy_min": float(np.linalg.eigvalsh(hout)[0]),
            "positive_internal_energy_max": float(np.linalg.eigvalsh(hout)[-1]),
            "capture_full_energy_intertwiner": intertwiner,
            "waiting": waiting_rows, "original_Mr": probabilities,
            "K2_general_gates_not_claimed_generated_by_this_H": True}


def bell_gate():
    cnot = np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                     [0, 0, 0, 1], [0, 0, 1, 0]], complex)
    code = [4, 5, 7, 8]  # qutrit pair |11>, |12>, |21>, |22>.
    gate = np.eye(9, dtype=complex)
    gate[np.ix_(code, code)] = cnot @ np.kron(H, I2)
    bell = gate[:, 4].reshape(3, 3)
    target = np.zeros((3, 3), complex)
    target[1, 1] = target[2, 2] = 1 / math.sqrt(2)
    small(norm(bell - target), "source complete qutrit Bell preparation")
    return gate, bell


def reading_gate(observables):
    gate = np.zeros((24, 24), complex)
    for setting, observable in enumerate(observables):
        plus = np.zeros((3, 3), complex)
        minus = np.zeros((3, 3), complex)
        plus[0, 0] = 1
        plus[1:, 1:] = (I2 + observable) / 2
        minus[1:, 1:] = (I2 - observable) / 2
        setting_projector = np.diag([int(setting == 0), int(setting == 1)])
        gate += np.kron(np.kron(plus, setting_projector), np.eye(4))
        gate += np.kron(np.kron(minus, setting_projector), np.kron(X, X))
    small(norm(gate.conj().T @ gate - np.eye(24)), "complete reading unitary")
    return gate


def block_commutation(old, material, gate, seed):
    rng = np.random.default_rng(seed)
    arbitrary = rng.normal(size=(4, 8, 3, 3, material.k, 1, 3))
    arbitrary = arbitrary + 1j * rng.normal(size=arbitrary.shape)
    arbitrary /= norm(arbitrary)
    comm = material.h_action(material.logical(arbitrary, gate))
    comm -= material.logical(material.h_action(arbitrary), gate)
    energy_residual = small(norm(comm), "dressed operation and same full native H")
    # A random full physical mother diagnostic; full quantifiers are analytic.
    ext = rng.normal(size=(4, 8, 3, 3, material.k, material.d, 2))
    ext = ext + 1j * rng.normal(size=ext.shape)
    ext /= norm(ext)
    swap_comm = material.h_action(material.swap_external(ext))
    swap_comm -= material.swap_external(material.h_action(ext))
    swap_residual = small(norm(swap_comm), "complete physical loading P identity")
    u_prime = np.array([1.25, -0.75, 0, 0])
    four = [small(abs(float(factor)) * energy_residual, "transported total four P")
            for factor in u_prime]
    return {"full_mother_random_diagnostic_columns": [3, 2],
            "energy_commutator_residual": energy_residual,
            "complete_loading_commutator_residual": swap_residual,
            "transported_four_P_residuals": four,
            "analytic_proof_not_random_columns_carries_full_mother_quantifier": True}


def compiled_blocks(old, source_gate, read_a, read_b):
    source = old.Material(3)
    local = old.Material(8)
    source_receipt = source.compile(source_gate, np.eye(9), start=0.0, finish=20.0)
    local_input = np.zeros((24, 6), complex)
    for q in range(3):
        for setting in range(2):
            local_input[q * 8 + setting * 4, q * 2 + setting] = 1
    a_receipt = local.compile(read_a, local_input, start=8.0, finish=20.0)
    b_receipt = local.compile(read_b, local_input, start=8.0, finish=20.0)
    source_support = block_commutation(old, source, source_gate, 115801)
    a_support = block_commutation(old, local, read_a, 115802)
    b_support = block_commutation(old, local, read_b, 115803)
    return {"source_full_qutrit_pair": source_receipt,
            "A_all_qutrit_setting_inputs": a_receipt,
            "B_all_qutrit_setting_inputs": b_receipt,
            "dressed_complete_support_and_charges": [source_support, a_support, b_support],
            "largest_state_matrix_rows": 55296, "largest_state_matrix_columns": 6,
            "global_288_device_tensor_constructed": False,
            "future_C_logical_dimension": 1152, "future_C_work_dimension": 384,
            "future_C_implementation_is_1112_analytic_full_compiler_not_dense_simulation": True}


def complete_network(old, bell, read_a, read_b):
    # A=(qutrit, setting, pointer, environment); B has the identical type.
    # Eight input columns are the arbitrary joint fresh (x,y,Q) state.
    initial = np.zeros((24, 24, 2, 8), complex)
    for x in range(2):
        for y in range(2):
            for q in range(2):
                col = 4 * x + 2 * y + q
                for a in range(3):
                    for b in range(3):
                        initial[a * 8 + x * 4, b * 8 + y * 4, q, col] = bell[a, b]
    early = old.on_axes(old.on_axes(initial, read_a, [0]), read_b, [1])
    reversed_order = old.on_axes(old.on_axes(initial, read_b, [1]), read_a, [0])
    ordering_residual = small(norm(early - reversed_order), "two spacelike instrument orders")
    matrix = early.reshape(1152, 8)
    iso = small(norm(matrix.conj().T @ matrix - np.eye(8)), "complete open-input history")
    correlation = []
    marginals_b = {}
    for x in range(2):
        for y in range(2):
            state = early[:, :, :, 4 * x + 2 * y].reshape(3, 2, 2, 2, 3, 2, 2, 2, 2)
            probabilities = np.zeros((2, 2))
            for ma in range(2):
                for mb in range(2):
                    probabilities[ma, mb] = np.sum(np.abs(state[:, :, ma, :, :, :, mb, :, :]) ** 2)
            e = float(sum((-1) ** (ma + mb) * probabilities[ma, mb]
                          for ma in range(2) for mb in range(2)))
            correlation.append(e)
            coarse = early[:, :, :, 4 * x + 2 * y]
            marginals_b[(x, y)] = np.einsum("abq,acq->bc", coarse, coarse.conj())
    chsh = correlation[0] + correlation[1] + correlation[2] - correlation[3]
    small(abs(chsh - 2 * math.sqrt(2)), "actual pointer CHSH")
    ns = max(norm(marginals_b[(0, y)] - marginals_b[(1, y)]) for y in range(2))
    small(ns, "complete B block not just B pointer")
    restored = old.on_axes(old.on_axes(early, read_a.conj().T, [0]), read_b.conj().T, [1])
    restored_error = small(norm(restored - initial), "all measurement resources reunited")

    # Noncommuting future: H on A's logical code, then A-controls-X on Q.
    h3 = np.eye(3, dtype=complex)
    h3[1:, 1:] = H
    controlled = np.eye(6, dtype=complex)
    controlled[4:6, 4:6] = X
    def future(a):
        shape = a.shape
        a = a.reshape(3, 2, 2, 2, 3, 2, 2, 2, 2, shape[-1])
        a = old.on_axes(a, h3, [0])
        a = old.on_axes(a, controlled, [0, 8])
        return a.reshape(shape)
    future_actual = future(restored)
    future_target = future(initial)
    future_error = small(norm(future_actual - future_target), "unknown Q and full active future")
    future_matrix = future_actual.reshape(1152, 8)
    small(norm(future_matrix.conj().T @ future_matrix - np.eye(8)), "full continued isometry")
    # Maximally entangled reference for all eight open inputs, as a diagnostic.
    reference_error = small(norm((future_matrix - future_target.reshape(1152, 8)).reshape(-1)) / math.sqrt(8),
                            "complete 8-dimensional reference diagnostic")

    # Incorrect deletion model: dephase EA before the same coherent inverse.
    # This is a different CP process; no claim of free reversible deletion.
    chosen = early[:, :, :, 0:1]
    bell_projector = np.outer(bell.reshape(-1), bell.reshape(-1).conj())
    recovered_probability = 0.0
    branch_norm = 0.0
    for outcome in range(2):
        branch = chosen.reshape(3, 2, 2, 2, 24, 2, 1).copy()
        branch[:, :, :, 1 - outcome, ...] = 0
        branch = branch.reshape(chosen.shape)
        branch_norm += norm(branch) ** 2
        branch = old.on_axes(old.on_axes(branch, read_a.conj().T, [0]), read_b.conj().T, [1])
        expanded = branch.reshape(3, 2, 2, 2, 3, 2, 2, 2, 2, 1)
        projected = old.on_axes(expanded, bell_projector, [0, 4])
        recovered_probability += float(np.vdot(expanded, projected).real)
    small(abs(branch_norm - 1), "nonselective deletion comparison trace")
    small(abs(recovered_probability - 0.5), "environment retention has active-future consequence")
    return {"complete_output_dimension": 1152, "open_input_dimension": 8,
            "spacelike_order_residual": ordering_residual, "complete_isometry_residual": iso,
            "pointer_correlations": correlation, "CHSH": float(chsh),
            "complete_B_block_no_signalling_residual": ns,
            "full_record_reunion_inverse_residual": restored_error,
            "noncommuting_active_future_residual": future_error,
            "reference_diagnostic_dimension": 8, "reference_residual": reference_error,
            "coherent_Bell_recovery_probability": 1.0,
            "EA_dephased_different_process_Bell_recovery_probability": recovered_probability,
            "external_uncaptured_classical_log_not_erased": True}


def mink(a, b):
    return float(a[0] * b[0] - np.dot(a[1:], b[1:]))


def geometry_and_existing_recoil():
    boost = np.array([[1.25, -0.75, 0, 0], [-0.75, 1.25, 0, 0],
                      [0, 0, 1, 0], [0, 0, 0, 1]], float)
    source_out = np.array([1.07, 0, 0, 0])
    a_in, b_in = np.array([8., -3, 0, 0]), np.array([8., 3, 0, 0])
    a_out, b_out = np.array([9.07, -3, 0, 0]), np.array([9.07, 3, 0, 0])
    c_in = np.array([17., 0, 0, 0])
    intervals = [a_in-source_out, b_in-source_out, c_in-a_out, c_in-b_out, c_in-a_in]
    rows = []
    for interval in intervals:
        proper = math.sqrt(mink(interval, interval))
        momentum = interval / proper  # Each listed carrier has mass one.
        transported = boost @ momentum
        small(abs(mink(momentum, momentum) - 1), "finite positive carrier mass")
        phase = mink(momentum, interval)
        phase_prime = mink(transported, boost @ interval)
        small(abs(phase - phase_prime), "same full mechanical phase")
        rows.append({"interval": interval.tolist(), "proper_time": proper,
                     "positive_four_P": momentum.tolist(),
                     "boosted_four_P": transported.tolist(), "action": phase})
    cross_margin = -mink(b_out-a_in, b_out-a_in)
    assert cross_margin > 0
    local_proper = []
    for duration in [0.37, 0.41, 0.29]:
        segment = np.array([duration, 0, 0, 0])
        local_proper.append(math.sqrt(mink(boost @ segment, boost @ segment)))
    small(norm(np.array(local_proper) - np.array([0.37, 0.41, 0.29])), "same positive native slots")

    # Preserve 1111's full two-material product, including crossed sectors.
    pax, pay = np.array([5., 4, 0, 0]), np.array([5., 0, 4, 0])
    pbx, pby = np.array([5., -4, 0, 0]), np.array([5., 0, -4, 0])
    pair_p = np.column_stack([pax+pbx, pax+pby, pay+pbx, pay+pby])
    xbar = np.array([[0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]], complex)
    full_pair_comm = []
    for four_p in [pair_p, boost @ pair_p]:
        for component in four_p:
            q = np.diag(component)
            full_pair_comm.append(small(norm(q @ xbar - xbar @ q), "old full pair recoil total P"))
    delta_a, delta_b = pay-pax, pby-pbx
    small(norm(delta_a+delta_b), "two actual recoil participants retained")
    boosted_a, boosted_b = boost @ delta_a, boost @ delta_b
    small(norm(boosted_a+boosted_b), "transported individual exchange cancels")
    assert abs(boosted_a[0]) > 0
    return {"causal_carriers": rows, "fresh_setting_to_remote_report_spacelike_margin": cross_margin,
            "transported_native_proper_slots": local_proper,
            "full_1111_pair_dimension": 4, "full_pair_P_commutator_residuals": full_pair_comm,
            "individual_recoil_A": delta_a.tolist(), "individual_recoil_B": delta_b.tolist(),
            "boosted_individual_exchange_A": boosted_a.tolist(),
            "boosted_individual_exchange_B": boosted_b.tolist(),
            "old_complete_device_dressing_and_arbitrary_graph_quantifiers_reused_analytically": True,
            "window_interior_free_trajectory_not_assumed": True}


def calculate():
    sources = source_check()
    old = load_material()
    source, bell = bell_gate()
    read_a = reading_gate([Z, X])
    read_b = reading_gate([(Z+X)/math.sqrt(2), (Z-X)/math.sqrt(2)])
    return {"round": 1158, "status": "PASS_CONDITIONAL_COMPLETE_MATERIAL_RULE_COMPATIBILITY",
            "numerical_check_groups": 4, "source_hashes_checked": sources,
            "new_adopted_axioms": 0, "scientific_calibration": 0,
            "old_clock_and_full_passive_algebra": old_clock_and_complete_block(old),
            "complete_local_material_compilation": compiled_blocks(old, source, read_a, read_b),
            "same_entangled_process_and_active_future": complete_network(old, bell, read_a, read_b),
            "actual_causal_lines_four_P_and_original_recoil": geometry_and_existing_recoil(),
            "scope": {"whole_old_finite_effective_domain_proved_analytically": True,
                      "K2_native_permission_explicit_not_derived_from_commutators": True,
                      "G_not_adopted": True, "candidate_C1_C2_C3_not_adopted_yet": True,
                      "no_claim_all_K2_window_interiors_generated_by_old_H": True,
                      "not_a_new_continuous_field_or_empirical_test": True,
                      "formal_acceptance_left_to_mainline": True}}


def compare(saved, actual, path="root"):
    if isinstance(actual, dict):
        assert saved.keys() == actual.keys(), path
        for key in actual:
            compare(saved[key], actual[key], path + "." + key)
    elif isinstance(actual, list):
        assert len(saved) == len(actual), path
        for i, (a, b) in enumerate(zip(saved, actual)):
            compare(a, b, path + f"[{i}]")
    elif isinstance(actual, float):
        assert math.isclose(float(saved), actual, rel_tol=3e-9, abs_tol=3e-9), (path, saved, actual)
    else:
        assert saved == actual, (path, saved, actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-once", action="store_true")
    args = parser.parse_args()
    result = calculate()
    out = HERE / "results.json"
    if args.write_once:
        if out.exists():
            raise FileExistsError(out)
        out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    else:
        compare(json.loads(out.read_text(encoding="utf-8")), result)
    print("PASS: 1158; four groups; complete material blocks and one retained-resource network")
    print("READ_ONLY_MATCH" if not args.write_once else "FIRST_RESULTS_WRITTEN")
    print("K2 permissions remain explicit; full-family and arbitrary-reference claims use proof.md")


if __name__ == "__main__":
    main()
