"""Round 1118: finite rest types, the same transported worldline, finite recoil.

Run without flags to recompute and compare the saved certificate.
--write-results writes only this round's results.json.
No older experiment is imported or executed. All numerical units have c = 1.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
SOURCES = (
    "archive_1106_/1111/proof.md",
    "archive_1106_/1111/results.json",
    "archive_1106_/1114/proof.md",
    "archive_1106_/1108/proof.md",
    "archive_1086_1105/1105/proof.md",
    "archive_1086_1105/_admission/dynamic_contract_audit.md",
)
TOL = 2e-11
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
PAULI = (X, Y, Z)
ETA = np.diag([1., -1., -1., -1.])


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def small(value):
    assert value < TOL, value
    return float(value)


def boost(v):
    """Plus-sign convention: the rest worldline is mapped to slope +v."""
    v = np.asarray(v, float)
    b2 = float(v @ v)
    assert b2 < 1
    gamma = 1 / np.sqrt(1 - b2)
    result = np.eye(4)
    result[0, 0] = gamma
    result[0, 1:] = gamma * v
    result[1:, 0] = gamma * v
    if b2:
        result[1:, 1:] += (gamma - 1) * np.outer(v, v) / b2
    small(np.linalg.norm(result.T @ ETA @ result - ETA))
    return result


def rotation(axis, theta):
    axis = np.asarray(axis, float)
    axis /= np.linalg.norm(axis)
    nx, ny, nz = axis
    cross = np.array([[0, -nz, ny], [nz, 0, -nx], [-ny, nx, 0]])
    r = np.cos(theta) * np.eye(3)
    r += (1 - np.cos(theta)) * np.outer(axis, axis) + np.sin(theta) * cross
    s = sum(axis[i] * PAULI[i] for i in range(3))
    u = np.cos(theta / 2) * I2 - 1j * np.sin(theta / 2) * s
    return r, u


def group1():
    # Internal clock and nontrivial spin coexist; the theorem needs no spin-0.
    m_internal = 3 * I2 + .2 * Z + .3 * X
    h_rest = np.kron(m_internal, I2)
    residuals = []
    vector_residuals = []
    for axis, theta in [([1, 0, 0], .41), ([0, 1, 0], .72), ([1, 2, 3], .63)]:
        r, u = rotation(axis, theta)
        full_u = np.kron(I2, u)
        residuals.append(opnorm(full_u.conj().T @ h_rest @ full_u - h_rest))
        # Deliberately omit commuting translations: Pauli vectors then give
        # a finite-dimensional rotational vector that is not forced to zero.
        for i in range(3):
            vector_residuals.append(opnorm(
                u.conj().T @ PAULI[i] @ u - sum(r[i, j] * PAULI[j] for j in range(3))))
    commutator = opnorm(X @ Y - Y @ X)
    finite_catalog = np.vstack((np.eye(3), -np.eye(3), np.zeros((1, 3))))
    r, _ = rotation([0, 0, 1], np.pi / 7)
    missed = min(np.linalg.norm(r @ np.array([1., 0., 0.]) - p) for p in finite_catalog)
    assert commutator > 1 and missed > .1
    return {
        "status": "PASS",
        "rest_type_dimension": 4,
        "rest_mass_eigenvalues": np.linalg.eigvalsh(h_rest).tolist(),
        "nontrivial_spin_allowed": True,
        "rotated_complete_H_residual": small(max(residuals)),
        "noncommuting_Pauli_vector_covariance_residual": small(max(vector_residuals)),
        "dropped_commutativity_counterexample_commutator_norm": commutator,
        "discrete_catalog_missing_rotated_direction_distance": float(missed),
        "finite_samples_do_not_prove_full_rotation_quantifier": True,
    }


def group2():
    m = np.kron(3 * I2 + .2 * Z + .3 * X, I2)
    zero = np.zeros_like(m)
    rest_p = [m, zero, zero, zero]
    rows = []
    for v in ([.6, 0., 0.], [.3, .2, .1]):
        a = boost(v)
        p = [sum(a[mu, nu] * rest_p[nu] for nu in range(4)) for mu in range(4)]
        event = a @ np.array([1.2, 0., 0., 0.])
        actual_v = event[1:] / event[0]
        shell = p[0] @ p[0] - sum(pi @ pi for pi in p[1:]) - m @ m
        rows.append({
            "declared_velocity": v,
            "transported_endpoint": event.tolist(),
            "same_worldline_velocity_residual": small(np.linalg.norm(actual_v - v)),
            "operator_p_minus_Ev_residual": small(max(opnorm(p[i + 1] - actual_v[i] * p[0]) for i in range(3))),
            "full_internal_mass_shell_residual": small(opnorm(shell)),
            "proper_interval_residual": small(abs(event @ ETA @ event - 1.2**2)),
        })
    return {
        "status": "PASS",
        "rows": rows,
        "positive_mass_coherence_operator_identity_only": True,
        "unknown_mass_preparation_or_acceleration_not_certified": True,
        "reference_change_is_not_scattering": True,
    }


def group3(source):
    params = source["parameters"]
    old = source["group4_material_recoil_and_actual_geometry"]
    m, energy, momentum = (float(params[k]) for k in ("m", "E", "k"))
    assert (m, energy, momentum) == (3., 5., 4.)
    assert float(params["collision_duration"]) > 0
    final_event = np.asarray(params["scatter_complete"], float)
    start_event = np.asarray(params["scatter_start"], float)
    input_p = [
        np.array([energy, momentum, 0., 0.]),
        np.array([energy, -momentum, 0., 0.]),
    ]
    output_p = [
        np.array([energy, 0., momentum, 0.]),
        np.array([energy, 0., -momentum, 0.]),
    ]
    impulse = [output_p[i] - input_p[i] for i in range(2)]
    dv = [output_p[i][1:] / energy - input_p[i][1:] / energy for i in range(2)]
    gamma = energy / m
    source_endpoint_residual = 0.
    source_velocity_residual = 0.
    for row in old["rows"]:
        p = np.asarray(row["P"], float)
        delta = np.asarray(row["endpoint"], float) - final_event
        sigma = float(row["proper_flight"])
        source_endpoint_residual = max(source_endpoint_residual, np.linalg.norm(delta - sigma * p / m))
        source_velocity_residual = max(source_velocity_residual, np.linalg.norm(delta[1:] / delta[0] - p[1:] / p[0]))
    a = boost([.6, 0., 0.])
    pin_b, pout_b = ([a @ p for p in pair] for pair in (input_p, output_p))
    boosted_impulse = [pout_b[i] - pin_b[i] for i in range(2)]
    rows = []
    for i in range(2):
        delta_event = a @ (float(params["flight"]) * output_p[i] / energy)
        actual_v = delta_event[1:] / delta_event[0]
        rows.append({
            "body": i,
            "input_four_momentum": pin_b[i].tolist(),
            "output_four_momentum": pout_b[i].tolist(),
            "impulse": boosted_impulse[i].tolist(),
            "output_event": (a @ final_event + delta_event).tolist(),
            "actual_velocity_matches_P_residual": small(np.linalg.norm(actual_v - pout_b[i][1:] / pout_b[i][0])),
        })
    pair_residual = np.linalg.norm(sum(impulse))
    boosted_pair_residual = np.linalg.norm(sum(boosted_impulse))
    # Deliberately retaining only A loses its compensating partner.
    omitted_partner_defect = np.linalg.norm(impulse[0])
    assert omitted_partner_defect > 1
    return {
        "status": "PASS",
        "source_round": 1111,
        "source_scattering_duration": float(params["collision_duration"]),
        "source_positive_window_check": bool(final_event[0] > start_event[0]),
        "m": m, "E": energy, "gamma": gamma,
        "actual_unboosted_impulses": [p.tolist() for p in impulse],
        "actual_velocity_changes": [v.tolist() for v in dv],
        "finite_turn_p_response_residual": small(max(np.linalg.norm(impulse[i][1:] - gamma * m * dv[i]) for i in range(2))),
        "source_endpoint_identity_residual": small(source_endpoint_residual),
        "source_velocity_identity_residual": small(source_velocity_residual),
        "whole_pair_conservation_residual": small(pair_residual),
        "boosted_pair_conservation_residual": small(boosted_pair_residual),
        "boosted_body_energy_transfers": [float(p[0]) for p in boosted_impulse],
        "boosted_rows": rows,
        "omitting_recoil_partner_defect": float(omitted_partner_defect),
        "unchanged_output_port_defect": float(old["wrong_unchanged_rail_endpoint_distance"]),
        "new_detector_or_bare_Q_readout_not_claimed": True,
        "old_device_experiment_not_rerun": True,
    }


def group4():
    m, bmax = 3., .1
    rows = []
    for b in (.01, .03, .05, .1):
        gamma = 1 / np.sqrt(1 - b*b)
        p, newton_p = gamma * m * b, m * b
        kinetic, newton_kinetic = m * (gamma - 1), .5 * m * b*b
        p_error = p - newton_p
        p_bound = m * b**3 / (np.sqrt(1-bmax*bmax) * (1+np.sqrt(1-bmax*bmax)))
        energy_error = kinetic - newton_kinetic
        energy_bound = 3*m*b**4/(8*(1-bmax*bmax)**2.5)
        assert -TOL <= p_error <= p_bound + TOL
        assert -TOL <= energy_error <= energy_bound + TOL
        rows.append({
            "speed": b,
            "momentum": float(p),
            "newton_momentum_error": float(p_error),
            "uniform_momentum_error_bound": float(p_bound),
            "kinetic_energy": float(kinetic),
            "newton_kinetic_error": float(energy_error),
            "uniform_kinetic_error_bound": float(energy_bound),
        })
    # Mathematical response tensor. No arbitrary small-impulse hardware is claimed.
    v = np.array([.06, .08, 0.])
    gamma = 1 / np.sqrt(1-v@v)
    jac = m * (gamma*np.eye(3) + gamma**3*np.outer(v, v))
    expected = np.array([m*gamma, m*gamma, m*gamma**3])
    return {
        "status": "PASS",
        "rows": rows,
        "response_tensor_eigenvalues": np.linalg.eigvalsh(jac).tolist(),
        "response_tensor_spectrum_residual": small(np.linalg.norm(np.linalg.eigvalsh(jac) - expected)),
        "uniform_speed_domain": [0., bmax],
        "arbitrary_impulse_actuator_not_built": True,
    }


def group5():
    # The very same pair has total rest energy 10, not the sum of bare masses 6.
    m, energy, k = 3., 5., 4.
    initial_pair = np.array([2*energy, 0., 0., 0.])
    a = boost([.6, 0., 0.])
    moved_pair = a @ initial_pair
    composite_mass = np.sqrt(moved_pair @ ETA @ moved_pair)
    h0 = .3*np.kron(Z, I2) + .2*np.kron(I2, Z)
    interaction = .15*np.kron(X, X)
    complete_h = 6*np.eye(4) + h0 + interaction
    masses = np.linalg.eigvalsh(complete_h)
    assert min(masses) > 0
    gamma = 1.25
    full_p = [gamma*complete_h, .6*gamma*complete_h, np.zeros((4, 4)), np.zeros((4, 4))]
    shell = full_p[0] @ full_p[0] - full_p[1] @ full_p[1] - complete_h @ complete_h
    scalar_binding_shift = -.1
    shifted_h = complete_h + scalar_binding_shift*np.eye(4)
    shifted_mass_difference = np.linalg.eigvalsh(shifted_h) - masses
    return {
        "status": "PASS",
        "actual_pair_composite_rest_mass": float(composite_mass),
        "bare_mass_sum": 2*m,
        "actual_pair_internal_relative_kinetic_energy": 2*(energy-m),
        "composite_mass_vs_total_energy_residual": small(abs(composite_mass - 2*energy)),
        "toy_complete_rest_spectrum": masses.tolist(),
        "complete_interaction_included_shell_residual": small(opnorm(shell)),
        "declared_binding_shift_only": scalar_binding_shift,
        "binding_shift_mass_difference": shifted_mass_difference.tolist(),
        "binding_shift_residual": small(np.linalg.norm(shifted_mass_difference - scalar_binding_shift)),
        "binding_dynamics_or_stability_not_derived": True,
        "null_sector_has_no_positive_rest_worldline_type": True,
    }


def calculate():
    hashes = {name: hashlib.sha256((RESEARCH/name).read_bytes()).hexdigest() for name in SOURCES}
    source = json.loads((RESEARCH/SOURCES[1]).read_text(encoding="utf-8"))
    return {
        "round": 1118,
        "status": "PASS_CONDITIONAL_FINITE_REST_TYPE_INERTIA_BRIDGE",
        "new_adopted_axioms": 0,
        "scientific_calibrations": 0,
        "check_groups": 5,
        "units": "c=1; four-vector ordering (E,p_x,p_y,p_z)",
        "scope": {
            "finite_internal_rest_type_not_full_quantum_position": True,
            "full_SO3_and_commuting_translation_generators_are_department_conditions": True,
            "same_transport_of_four_P_and_original_CO_worldline_required": True,
            "massive_declared_material_sector_only": True,
            "existing_1111_effective_scattering_permission_reused": True,
            "not_full_joint_FUCP_CO_SR_or_all_matter_completion": True,
            "no_new_axiom_adopted": True,
        },
        "source_sha256": hashes,
        "groups": {
            "1_finite_rest_rotations_and_missing_conditions": group1(),
            "2_same_four_momentum_worldline_transport": group2(),
            "3_actual_conservative_turn_and_recoil": group3(source),
            "4_finite_low_velocity_error_bounds": group4(),
            "5_complete_composite_energy_and_mass": group5(),
        },
    }


def compare(a, b, path="$"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b), path
        for key in a:
            compare(a[key], b[key], path+"."+key)
    elif isinstance(a, (list, tuple)):
        assert isinstance(b, (list, tuple)) and len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, f"{path}[{i}]")
    elif isinstance(a, bool) or isinstance(a, str) or a is None:
        assert type(a) is type(b) and a == b, (path, a, b)
    elif isinstance(a, (int, float)):
        assert isinstance(b, (int, float)) and not isinstance(b, bool), path
        assert np.isfinite(a) and np.isfinite(b), path
        assert abs(a-b) <= 2e-11*(1+max(abs(a), abs(b))), (path, a, b)
    else:
        raise TypeError((path, type(a)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    calculated = calculate()
    target = HERE/"results.json"
    if args.write_results:
        target.write_text(json.dumps(calculated, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        mode = "saved"
    else:
        compare(calculated, json.loads(target.read_text(encoding="utf-8")))
        mode = "saved results matched"
    for name, group in calculated["groups"].items():
        print(f"PASS {name}: {group['status']}")
    print(f"PASS 5 groups; {mode}; no older device rerun; conditional department only.")


if __name__ == "__main__":
    main()