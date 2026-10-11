"""1112: conservative logical compilation with complete loading and offloading.

The same 1106 Hamiltonian is retained.  Only its finite degenerate work
extension and declared local dressed operations are used.  Default is
read-only comparison; --write-results writes the deterministic receipt.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.json"
TOL = 2e-9


def residual(a):
    return float(np.linalg.norm(a))


def small(value, label):
    value = float(value)
    if value > TOL:
        raise AssertionError(f"{label}: {value}")
    return value


def tensor(*items):
    out = np.array([1.0 + 0j])
    for item in items:
        out = np.kron(out, item)
    return out


def exp_h(h, t):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * t * e)) @ v.conj().T


def on_axes(a, op, axes):
    others = [i for i in range(a.ndim) if i not in axes]
    order = list(axes) + others
    shape = [a.shape[i] for i in order]
    moved = np.transpose(a, order).reshape(op.shape[1], -1)
    out = (op @ moved).reshape(shape)
    return np.transpose(out, np.argsort(order))


class Material:
    """Physical axes: head, ab/D/M, QC, Qr, W, external memory, input."""

    def __init__(self, k):
        self.k, self.d = k, 3 * k
        i2 = np.eye(2, dtype=complex)
        y = np.array([[0, -1j], [1j, 0]], complex)
        x = np.array([[0, 1], [1, 0]], complex)
        p0, p1 = np.diag([1, 0]), np.diag([0, 1])
        b = -y
        ur = np.kron((i2 + b) / 2, i2) - 1j * np.kron((i2 - b) / 2, y)
        um = np.kron(p0, i2) - 1j * np.kron(p1, x)
        sw = np.zeros((9, 9), complex)
        for q in range(3):
            for r in range(3):
                sw[3 * r + q, 3 * q + r] = 1
        self.gates = [np.kron(np.eye(8), sw), np.kron(ur, np.eye(18)),
                      np.kron(np.kron(i2, um), np.eye(9))]
        self.products = [np.eye(72, dtype=complex)]
        self.hc = 1.5 * np.eye(4)
        self.hfk = 1.5 * np.eye(288, dtype=complex)
        for j, gate in enumerate(self.gates):
            self.products.append(gate @ self.products[-1])
            coupling = math.sqrt((j + 1) * (3 - j)) / 2
            self.hc[j + 1, j] = self.hc[j, j + 1] = coupling
            self.hfk[(j + 1) * 72:(j + 2) * 72, j * 72:(j + 1) * 72] = coupling * gate
            self.hfk[j * 72:(j + 1) * 72, (j + 1) * 72:(j + 2) * 72] = coupling * gate.conj().T
        self.gamma = np.zeros((288, 288), complex)
        for j, v in enumerate(self.products):
            self.gamma[j * 72:(j + 1) * 72, j * 72:(j + 1) * 72] = v
        self.baseline = 15.0 + 2 * math.sqrt(19) + 2.0  # finite work material rest mass
        self.aux = tensor(np.array([1, 1j]) / math.sqrt(2), [1, 0], [1, 0])

    def gamma_apply(self, a, adjoint=False):
        shape = a.shape
        q = a.reshape(4, 72, -1)
        out = np.empty_like(q)
        for j, v in enumerate(self.products):
            out[j] = (v.conj().T if adjoint else v) @ q[j]
        return out.reshape(shape)

    def h_action(self, a):
        # Direct independently specified FK blocks, not a definition by conjugation.
        shape = a.shape
        q = a.reshape(4, 72, -1)
        out = (self.baseline + 1.5) * q.copy()
        for j, gate in enumerate(self.gates):
            coupling = math.sqrt((j + 1) * (3 - j)) / 2
            out[j + 1] += coupling * (gate @ q[j])
            out[j] += coupling * (gate.conj().T @ q[j + 1])
        return out.reshape(shape)

    def wait(self, a, duration):
        q = self.gamma_apply(a, True)
        q = on_axes(q, exp_h(self.hc, duration), [0])
        return np.exp(-1j * self.baseline * duration) * self.gamma_apply(q)

    def logical(self, a, x):
        return self.gamma_apply(on_axes(self.gamma_apply(a, True), x, [2, 4]))

    def swap_external(self, a):
        if a.shape[5] != self.d:
            raise AssertionError("a true external full logical slot is required")
        q = self.gamma_apply(a, True)
        q = q.reshape(4, 8, 3, 3, self.k, 3, self.k, a.shape[-1])
        q = np.transpose(q, (0, 1, 5, 3, 6, 2, 4, 7)).copy()
        return self.gamma_apply(q.reshape(a.shape))

    def encode_joint(self, joint, t):
        # joint has [logical, external, arbitrary input columns].
        phi = np.exp(-1j * self.baseline * t) * exp_h(self.hc, t)[:, 0]
        j = joint.reshape(3, self.k, joint.shape[1], joint.shape[2])
        out = np.einsum("h,a,cwen,r->hacrwen", phi, self.aux, j, [1, 0, 0])
        return self.gamma_apply(out)

    def blank_joint(self, external, t):
        j = np.zeros((self.d, external.shape[0], external.shape[1]), complex)
        j[0] = external
        return self.encode_joint(j, t)

    def loaded_joint(self, logical, t):
        j = np.zeros((self.d, self.d, logical.shape[1]), complex)
        j[:, 0] = logical
        return self.encode_joint(j, t)

    def compile(self, x, vin, start=4.0, durations=(0.37, 0.41, 0.29), finish=10.0):
        dl, dx, du = durations
        input_full = self.blank_joint(vin, start)
        loaded = self.swap_external(self.wait(input_full, dl))
        load_err = residual(loaded - self.loaded_joint(vin, start + dl))
        # The incoming true carrier is exactly blank and uncorrelated; it is
        # retained as a known spectator, not traced from an entangled state.
        device = loaded[:, :, :, :, :, 0:1, :].copy()
        leakage = residual(loaded[:, :, :, :, :, 1:, :])
        operated = self.logical(self.wait(device, dx), x)
        fresh_out = np.zeros((4, 8, 3, 3, self.k, self.d, vin.shape[1]), complex)
        fresh_out[:, :, :, :, :, 0:1, :] = operated
        released = self.swap_external(self.wait(fresh_out, du))
        end = start + dl + dx + du
        target = self.blank_joint(x @ vin, end)
        endpoint_err = residual(released - target)
        if finish < end:
            raise AssertionError("the common completion precedes the declared slots")
        after = self.wait(released, finish - end)
        future_err = residual(after - self.blank_joint(x @ vin, finish))
        for val, label in [(load_err, "all-input load"), (leakage, "incoming carrier product blank"),
                           (endpoint_err, "actual external output"), (future_err, "continuing H after offload")]:
            small(val, label)
        return {"load_residual": load_err, "incoming_carrier_leakage": leakage,
                "external_complete_output_residual": endpoint_err,
                "continuing_device_future_residual": future_err,
                "largest_isometry_shape": [288 * self.k * self.d, vin.shape[1]],
                "physical_offload_time": end, "common_completion": finish}


def complete_columns(v):
    n, r = v.shape
    small(residual(v.conj().T @ v - np.eye(r)), "isometry supplied for extension")
    columns = [v[:, i].copy() for i in range(r)]
    for seed in np.eye(n, dtype=complex).T:
        z = seed.copy()
        for _ in range(2):
            for c in columns:
                z -= c * np.vdot(c, z)
        size = np.linalg.norm(z)
        if size > 1e-10:
            columns.append(z / size)
        if len(columns) == n:
            break
    out = np.column_stack(columns)
    small(residual(out.conj().T @ out - np.eye(n)), "full unitary extension")
    return out


def sqrt_positive(a):
    e, v = np.linalg.eigh(a)
    assert min(e) > -1e-12
    return (v * np.sqrt(np.maximum(0, e))) @ v.conj().T


def one_occupied(count, payload_dim):
    k = payload_dim + 1
    out = np.zeros((k ** count, count * payload_dim), complex)
    for port in range(count):
        for value in range(payload_dim):
            digits = [0] * count
            digits[port] = value + 1
            row = sum(a * k ** (count - i - 1) for i, a in enumerate(digits))
            out[row, port * payload_dim + value] = 1
    return out


def instrument(kraus):
    # Keep both a readable outcome and a distinct dilation environment.
    nout, nin = kraus[0].shape
    out = np.zeros((nout * 4, nin), complex)
    out.reshape(nout, 2, 2, nin)[:, 0, 0] = kraus[0]
    out.reshape(nout, 2, 2, nin)[:, 1, 1] = kraus[1]
    return out


def boost_matrix(v):
    v = np.asarray(v, float)
    b2 = float(v @ v)
    gamma = 1 / math.sqrt(1 - b2)
    out = np.eye(4)
    out[0, 0] = gamma
    out[0, 1:] = out[1:, 0] = -gamma * v
    out[1:, 1:] += (gamma - 1) * np.outer(v, v) / b2
    return out


def mdot(a, b):
    return float(a[0] * b[0] - a[1:] @ b[1:])


def calculate():
    base = Material(1)
    ih = np.kron(base.hc, np.eye(72))
    inter = small(residual(base.hfk @ base.gamma - base.gamma @ ih), "same fixed 1106 H")
    full = np.eye(288, dtype=complex).reshape(4, 8, 3, 3, 1, 1, 288)
    q = np.arange(3)
    fourier = np.exp(2j * np.pi * np.outer(q, q) / 3) / math.sqrt(3)
    phase = np.diag(np.exp(1j * np.array([0.13, 0.41, -0.29])))
    ox = base.logical(full, fourier)
    algebra = small(residual(base.logical(base.logical(full, phase), fourier)
                             - base.logical(full, fourier @ phase)), "finite algebra product")
    comm = small(residual(base.h_action(ox) - base.logical(base.h_action(full), fourier)), "full-mother H commutator")
    inverse = small(residual(base.logical(ox, fourier.conj().T) - full), "full inverse")
    time_law = small(residual(base.logical(base.wait(base.logical(base.wait(full, 0.17), phase), 0.23), fourier)
                              - base.logical(base.wait(full, 0.40), fourier @ phase)), "word includes all waits")

    # The complete external 3-dimensional slot, not just two selected states.
    ext_full = np.eye(864, dtype=complex).reshape(4, 8, 3, 3, 1, 3, 864)
    sx = base.swap_external(ext_full)
    swap_comm = small(residual(base.h_action(sx) - base.swap_external(base.h_action(ext_full))), "full load/offload commutator")
    u = np.array([5 / 4, -3 / 4, 0, 0])
    pcarrier = np.array([5., 3., 0., 0.])
    four_p_errors = []
    for mu in range(4):
        pa = u[mu] * base.h_action(ext_full) + pcarrier[mu] * ext_full
        ps = u[mu] * base.h_action(sx) + pcarrier[mu] * sx
        four_p_errors.append(residual(ps - base.swap_external(pa)))
    small(max(four_p_errors), "all four components with moving true carrier")
    bridge = base.compile(fourier, np.eye(3), start=4.0)
    del ext_full, sx, full, ox

    # A genuinely quantum, non-SWAP instrument on the occupied source qubit.
    pauli_x = np.array([[0, 1], [1, 0]], complex)
    kp2 = sqrt_positive((np.eye(2) + 0.6 * pauli_x) / 2)
    km2 = sqrt_positive((np.eye(2) - 0.6 * pauli_x) / 2)
    kp, km = np.zeros((3, 3), complex), np.zeros((3, 3), complex)
    kp[0, 0] = 1
    kp[1:, 1:], km[1:, 1:] = kp2, km2
    vinst = instrument([kp, km])
    ui = complete_columns(vinst)
    quantum = Material(4).compile(ui, np.eye(12)[:, :3])
    plus = np.array([0, 1, 1], complex) / math.sqrt(2)
    minus = np.array([0, 1, -1], complex) / math.sqrt(2)
    probabilities = [float(np.linalg.norm(kp @ a) ** 2) for a in [plus, minus]]
    small(abs(probabilities[0] - 0.8) + abs(probabilities[1] - 0.2), "nontrivial clock-independent instrument")
    choi = (np.kron(vinst[:, 1], [1, 0]) + np.kron(vinst[:, 2], [0, 1])) / math.sqrt(2)
    small(abs(np.linalg.norm(choi) - 1), "retained arbitrary-reference witness")

    # Original public probe and tent, with coherent positions and full records.
    positions = [-0.5, 0.75, 1.25]
    public_k = [np.zeros((6, 6), complex), np.zeros((6, 6), complex)]
    tent_k = [np.zeros((6, 6), complex), np.zeros((6, 6), complex)]
    effect_eigenvalues = []
    for i, pos in enumerate(positions):
        effect = np.eye(2) / 2 + pos / math.sqrt(1 + pos * pos) * pauli_x / 4
        effect_eigenvalues.append(np.linalg.eigvalsh(effect).tolist())
        public_k[0][2 * i:2 * i + 2, 2 * i:2 * i + 2] = sqrt_positive(effect)
        public_k[1][2 * i:2 * i + 2, 2 * i:2 * i + 2] = sqrt_positive(np.eye(2) - effect)
        tau = max(0., 1 - abs(pos))
        tent_k[0][2 * i:2 * i + 2, 2 * i:2 * i + 2] = math.sqrt(tau) * np.eye(2)
        tent_k[1][2 * i:2 * i + 2, 2 * i:2 * i + 2] = math.sqrt(1 - tau) * np.eye(2)
    vp, vt = instrument(public_k), instrument(tent_k)
    small(residual(vp.conj().T @ vp - np.eye(6)), "public complete instrument")
    small(residual(vt.conj().T @ vt - np.eye(6)), "tent complete instrument")
    js = one_occupied(3, 2)
    port_public = np.kron(js, np.eye(4)) @ vp
    small(residual(port_public.conj().T @ port_public - np.eye(6)), "actual public-probe output ports")
    # The larger hardware already exists; subcatalogue is a support restriction.
    jsmall = np.eye(6)[:, :4]
    smaller = instrument([k[:4, :4] for k in public_k])
    natural = small(residual(vp @ jsmall - np.kron(jsmall, np.eye(4)) @ smaller), "same-context complete naturality")

    # CO meeting: keep participant X's unknown qubit, both positions, and add
    # the actual midpoint marker.  All three groups are real output carriers.
    sxpos, typos, zpos = [-1., 1.], [-1., 1.], [-1., 0., 1.]
    vm = np.zeros((24, 8), complex)
    for ix in range(2):
        for qbit in range(2):
            for iy in range(2):
                col = (2 * ix + qbit) * 2 + iy
                iz = zpos.index((sxpos[ix] + typos[iy]) / 2)
                row = ((2 * ix + qbit) * 2 + iy) * 3 + iz
                vm[row, col] = 1
    jx, jy, jz = one_occupied(2, 2), one_occupied(2, 1), one_occupied(3, 1)
    jm = tensor(jx, jy, jz)
    actual_meeting = jm @ vm
    meeting_iso = small(residual(actual_meeting.conj().T @ actual_meeting - np.eye(8)), "all unknown meeting input")
    midpoint_error, participant_error = 0., 0.
    for ix in range(2):
        for qbit in range(2):
            for iy in range(2):
                col = (2 * ix + qbit) * 2 + iy
                output = actual_meeting[:, col].reshape(9, 4, 8)
                px, py, pz = np.sum(abs(output) ** 2, axis=(1, 2)), np.sum(abs(output) ** 2, axis=(0, 2)), np.sum(abs(output) ** 2, axis=(0, 1))
                expected_x, expected_y = jx[:, 2 * ix + qbit], jy[:, iy]
                participant_error = max(participant_error, residual(px - abs(expected_x) ** 2), residual(py - abs(expected_y) ** 2))
                zmean = sum(zpos[j] * float(abs(jz[:, j]) ** 2 @ pz) for j in range(3))
                midpoint_error = max(midpoint_error, abs(zmean - (sxpos[ix] + typos[iy]) / 2))
    small(participant_error, "original participants really restored")
    small(midpoint_error, "actual new meeting port")
    meeting = Material(8).compile(complete_columns(vm), np.eye(24)[:, :8])
    coherent = np.zeros(8, complex)
    coherent[0], coherent[7] = 1 / math.sqrt(2), 1j / math.sqrt(2)
    small(abs(np.linalg.norm(actual_meeting @ coherent) - 1), "coherent two-position meeting")

    # All comparisons use the same deadline; moving whole event slices does
    # not mean selecting a fresh t'=constant surface.
    eta = np.diag([1., -1., -1., -1.])
    lorentz = boost_matrix([0.3, 0.2, 0.1])
    small(residual(lorentz.T @ eta @ lorentz - eta), "reference matrix")
    start, end, completion = 4., 5.07, 10.
    edges, invariant_errors, momentum_residuals = [], [], []
    mass = 3.
    for xcoord in sorted(set(sxpos + typos + positions)):
        edges.append((np.array([0., xcoord, 0, 0]), np.array([start, 0, 0, 0])))
    for xcoord in sorted(set(sxpos + typos + zpos + positions)):
        edges.append((np.array([end, 0, 0, 0]), np.array([completion, xcoord, 0, 0])))
    minimum_interval = float("inf")
    total_action = 0.
    for a, b in edges:
        dx = b - a
        interval = mdot(dx, dx)
        assert dx[0] > 0 and interval > 0
        minimum_interval = min(minimum_interval, interval)
        p = mass * dx / math.sqrt(interval)
        total_action += mdot(p, dx)
        invariant_errors.append(abs(mdot(lorentz @ p, lorentz @ dx) - mdot(p, dx)))
        momentum_residuals.append(abs(mdot(lorentz @ p, lorentz @ p) - mass * mass))
    # Local positive windows remain on the same processor worldline.
    slots = [start, start + 0.37, start + 0.37 + 0.41, end]
    slot_errors = [abs(mdot(lorentz @ np.array([b - a, 0, 0, 0]), lorentz @ np.array([b - a, 0, 0, 0])) - (b - a) ** 2)
                   for a, b in zip(slots, slots[1:])]
    small(max(invariant_errors + momentum_residuals + slot_errors), "same rail/P/clock transport")
    source_times = [(lorentz @ np.array([0., x, 0, 0]))[0] for x in sxpos]
    assert abs(source_times[0] - source_times[1]) > 0.1

    # Identity logical data do not erase a different complete clock history.
    t_a, t_b = 5.07, 5.47
    phi_a, phi_b = exp_h(base.hc, t_a)[:, 0], exp_h(base.hc, t_b)[:, 0]
    distinguishability = math.sqrt(max(0., 1 - abs(np.vdot(phi_a, phi_b)) ** 2))
    expected_dist = math.sqrt(1 - math.cos((t_b - t_a) / 2) ** 6)
    small(abs(distinguishability - expected_dist), "complete history cannot be dropped")
    assert distinguishability > 0.3
    sources = {}
    for name in ["1106/check.py", "1106/proof.md", "1111/proof.md"]:
        sources[name] = hashlib.sha256((HERE.parent / name).read_bytes()).hexdigest()
    return {
        "status": "PASS_CONDITIONAL_CONSERVATIVE_LOGICAL_PORT_CATEGORY", "round": 1112,
        "new_adopted_axioms": 0, "scientific_calibration": 0, "numerical_check_groups": 8,
        "sources_sha256": sources,
        "scope": {"fixed_1106_H_retained": True, "declared_positive_local_dressed_slots_are_material_permission": True,
                  "exact_logical_category_not_complete_nonscalar_mother_C": True,
                  "all_complete_records_resources_and_continuing_device_retained": True,
                  "original_ports_not_only_central_log": True, "all_SR_or_reference_family_closure_not_claimed": True},
        "group1_same_H_and_algebra": {"intertwiner_residual": inter, "composition_residual": algebra,
                                     "full_mother_commutator_residual": comm, "inverse_residual": inverse,
                                     "positive_slot_word_residual": time_law},
        "group2_true_external_complete_type_bridge": {"full_swap_dimension": 864,
            "all_input_swap_commutator_residual": swap_comm, "four_P_residual": max(four_p_errors), **bridge},
        "group3_non_SWAP_quantum_instrument": {"complete_stinespring_dimension": 12,
            "plus_effect_probabilities_on_X_eigenstates": probabilities,
            "readable_outcome_and_environment_both_retained": True, **quantum},
        "group4_original_public_probe_tent_and_naturality": {"positions": positions,
            "public_effect_eigenvalues": effect_eigenvalues, "public_port_isometry_shape": list(port_public.shape),
            "tent_complete_isometry_residual": residual(vt.conj().T @ vt - np.eye(6)),
            "public_complete_isometry_residual": residual(vp.conj().T @ vp - np.eye(6)),
            "same_window_directory_naturality_residual": natural},
        "group5_actual_CO_meeting_full_input": {"input_dimension": 8, "output_physical_ports_shape": list(actual_meeting.shape),
            "original_participant_port_residual": participant_error, "new_midpoint_port_residual": midpoint_error,
            "complete_isometry_residual": meeting_iso, "unknown_qubit_and_position_coherence_retained": True, **meeting},
        "group6_fixed_window_whole_reference_transport": {"source_time": 0., "hub_slot_boundaries": slots,
            "original_target_completion": completion, "minimum_rail_interval_squared": minimum_interval,
            "mechanical_action_sum_all_rails": total_action, "action_transport_residual": max(invariant_errors),
            "positive_mass_shell_residual": max(momentum_residuals), "proper_slot_residual": max(slot_errors),
            "transformed_original_source_times": source_times},
        "group7_complete_process_not_reduced_channel": {"different_completion_times": [t_a, t_b],
            "complete_clock_half_trace_distance": distinguishability,
            "same_logical_channel_is_not_complete_process_identity": True,
            "only_predeclared_logical_slots_and_record_instruments_compiled": True},
        "group8_permissions_and_remaining_coverage": {"local_geometry_not_inferred_from_commutator": True,
            "same_materials_for_vac_and_occupied_rails": True,
            "full_mother_operation_extension_defined_but_no_actual_arbitrary_arrival_filter": True,
            "finite_comparison_window_fixed_before_unknown_input": True,
            "no_universal_processor_or_normal_wavepacket_requirement": True,
            "no_new_scientific_calibration": True}
    }


def compare(saved, actual, path="root"):
    if isinstance(actual, dict):
        if not isinstance(saved, dict) or saved.keys() != actual.keys():
            raise AssertionError(f"{path}: keys differ")
        for key in actual:
            compare(saved[key], actual[key], f"{path}.{key}")
    elif isinstance(actual, list):
        if not isinstance(saved, list) or len(saved) != len(actual):
            raise AssertionError(f"{path}: lengths differ")
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
    print("PASS: 8 numerical groups; largest full isometry 55296 x 8")
    print("Same continuing H, complete conservative loading/offloading and original CO ports")
    print("Declared effective category; scientific calibration 0; whole SR completion not claimed")


if __name__ == "__main__":
    main()
