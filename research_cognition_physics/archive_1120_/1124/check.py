"""1124: continuous cap geometry and finite relational scattering checks.

Default recomputes and compares results.json without writing. --write-results
creates the initial artifact. C_N and three-mode examples are algebra probes,
not finite approximations claiming a continuum relativistic controller.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
TOL = 3e-12
ALPHA = math.pi / 6


def norm(matrix):
    return float(np.linalg.norm(matrix))


def cap_average(theta, alpha):
    sine = math.sin(theta)
    if sine < math.cos(alpha):
        return 0.0
    return math.acos(min(1.0, max(-1.0, math.cos(alpha) / sine))) / math.pi


def geometry_group():
    # Independent direct orbit counting, including rotations of the reference.
    count = 65536
    phi = 2 * math.pi * (np.arange(count) + .5) / count
    rows = []
    largest = 0.0
    for theta in [0.0, .7, math.pi / 3, 1.17, 1.4, math.pi / 2, 1.9, 2.3]:
        expected = cap_average(theta, ALPHA)
        values = []
        for reference in [0.0, .23, 1.91]:
            observed = float(np.mean(math.sin(theta) * np.cos(phi - reference)
                                     >= math.cos(ALPHA)))
            values.append(observed)
            largest = max(largest, abs(observed - expected))
        rows.append({"theta": theta, "formula": expected,
                     "direct_orbit_fractions": values})
    assert largest <= 2 / count
    # The SO(3) cap-area fraction differs from the SO(2)_x orbit supremum.
    nodes, weights = leggauss(512)
    u = math.sin(ALPHA) * nodes  # u = cos(theta), on the nonzero interval
    values = np.arccos(np.clip(math.cos(ALPHA) / np.sqrt(1 - u * u), -1, 1)) / math.pi
    average_sphere = float(.5 * math.sin(ALPHA) * np.dot(weights, values))
    area_fraction = (1 - math.cos(ALPHA)) / 2
    assert abs(average_sphere - area_fraction) < 1e-8
    analytic_sup = ALPHA / math.pi
    assert abs(cap_average(math.pi / 2, ALPHA) - analytic_sup) < TOL
    return {"alpha": ALPHA, "analytic_orbit_supremum": analytic_sup,
            "no_resource_task_error_floor": 1 - analytic_sup,
            "orbit_grid_points": count, "max_direct_orbit_error": largest,
            "sphere_average_quadrature": average_sphere,
            "sphere_cap_area_fraction": area_fraction,
            "sphere_average_residual": abs(average_sphere - area_fraction),
            "samples": rows,
            "scope": "Sampling diagnoses formula (6); the continuum supremum is analytic."}


def perm_matrix(permutation):
    out = np.zeros((len(permutation), len(permutation)), dtype=complex)
    out[permutation, np.arange(len(permutation))] = 1
    return out


def permute_density(rho, permutation):
    result = np.empty_like(rho)
    result[np.ix_(permutation, permutation)] = rho
    return result


def finite_reference_group():
    n = 6
    modes = n + 1  # x is 0; y_j is j+1
    dim_active = 2 * n * modes  # source Q, active controller C, material direction
    dim = dim_active * n  # additional inaccessible reference R

    def active_index(q, c, m):
        return (q * n + c) * modes + m

    def index(q, c, m, r):
        return active_index(q, c, m) * n + r

    def rotated_m(m, shift):
        return 0 if m == 0 else 1 + (m - 1 + shift) % n

    active_rotation = np.empty(dim_active, dtype=int)
    executor = np.empty(dim_active, dtype=int)
    full_executor = np.empty(dim, dtype=int)
    reference_executor = np.empty(dim, dtype=int)
    active_permutations, global_permutations = [], []
    for shift in range(n):
        pa, pg = np.empty(dim, dtype=int), np.empty(dim, dtype=int)
        for q in range(2):
            for c in range(n):
                for m in range(modes):
                    mm = c + 1 if m == 0 else (0 if m == c + 1 else m)
                    mm = mm if q == 1 else m
                    ai = active_index(q, c, m)
                    executor[ai] = active_index(q, c, mm)
                    active_rotation[ai] = active_index(q, (c + 1) % n, rotated_m(m, 1))
                    for r in range(n):
                        ii = index(q, c, m, r)
                        full_executor[ii] = index(q, c, mm, r)
                        rm = r + 1 if m == 0 else (0 if m == r + 1 else m)
                        reference_executor[ii] = index(q, c, rm if q else m, r)
                        pa[ii] = index(q, (c + shift) % n, rotated_m(m, shift), r)
                        pg[ii] = index(q, (c + shift) % n, rotated_m(m, shift), (r + shift) % n)
        active_permutations.append(pa)
        global_permutations.append(pg)
    u, rotation = perm_matrix(executor), perm_matrix(active_rotation)
    unitarity = norm(u.conj().T @ u - np.eye(dim_active))
    covariance = norm(u @ rotation - rotation @ u)
    assert max(unitarity, covariance) < TOL

    e = np.zeros(dim)
    for q in range(2):
        for c in range(n):
            for r in range(n):
                e[index(q, c, r + 1, r)] = 1.0
    average_effect = sum(e[pa] for pa in active_permutations) / n
    global_effect_residual = max(float(np.max(np.abs(e[pg] - e))) for pg in global_permutations)
    assert abs(float(np.max(average_effect)) - 1 / n) < TOL
    assert global_effect_residual < TOL

    # Complete isometry on Q and every active program value, not only basis probabilities.
    jin = np.zeros((dim_active, 2 * n), dtype=complex)
    jout = np.zeros_like(jin)
    for q in range(2):
        for c in range(n):
            jin[active_index(q, c, 0), q * n + c] = 1
            jout[active_index(q, c, c + 1 if q else 0), q * n + c] = 1
    isometry = norm(jout.conj().T @ jout - np.eye(2 * n))
    complete_gate = norm(u @ jin - jout)
    rng = np.random.default_rng(1124)
    amplitudes = rng.normal(size=(2 * n, n * 3)) + 1j * rng.normal(size=(2 * n, n * 3))
    amplitudes /= np.linalg.norm(amplitudes)
    passive_reference_residual = norm(u @ jin @ amplitudes - jout @ amplitudes)
    assert max(isometry, complete_gate, passive_reference_residual) < TOL

    def initial_density(q_state, distribution):
        rho = np.zeros((dim, dim), dtype=complex)
        for c in range(n):
            for r in range(n):
                for q in range(2):
                    for qp in range(2):
                        rho[index(q, c, 0, r), index(qp, c, 0, r)] = q_state[q, qp] * distribution[c, r]
        return rho

    q1 = np.diag([0., 1.])
    plus = np.ones((2, 2)) / 2
    aligned = np.eye(n) / n
    product = np.ones((n, n)) / n ** 2
    rows = []
    largest_global = 0.0
    largest_recovery = 0.0
    largest_saturation = 0.0
    largest_matrix_unit_recovery = 0.0
    for eta in [0., .25, .7, 1.]:
        distribution = (1 - eta) * product + eta * aligned
        relative_average = np.ones((n, 1)) @ np.sum(distribution, axis=0, keepdims=True) / n
        a = float(np.sum(np.abs(distribution - relative_average)) / 2)
        rho = initial_density(q1, distribution)
        output = permute_density(rho, full_executor)
        probability = float(np.dot(e, np.diag(output).real))
        averaged_output = sum(permute_density(output, pa) for pa in active_permutations) / n
        probability_average = float(np.dot(e, np.diag(averaged_output).real))
        global_output = sum(permute_density(output, pg) for pg in global_permutations) / n
        largest_global = max(largest_global, norm(output - global_output))
        largest_saturation = max(largest_saturation, abs(probability - (1 / n + a)))
        assert abs(probability_average - 1 / n) < TOL
        # Active marginal can be invariant while joint reference correlations carry direction.
        assert np.max(np.abs(np.sum(distribution, axis=1) - 1 / n)) < TOL
        coherent_input = initial_density(plus, distribution)
        coherent_output = permute_density(coherent_input, full_executor)
        recovered = permute_density(coherent_output, full_executor)  # U is self-adjoint
        largest_recovery = max(largest_recovery, norm(recovered - coherent_input))
        q_marginal = np.einsum("aibi->ab", coherent_output.reshape(2, dim // 2, 2, dim // 2))
        assert abs(q_marginal[0, 1]) < TOL  # path record, not loss of complete Q
        for q in range(2):
            for qp in range(2):
                unit = np.zeros((2, 2), dtype=complex)
                unit[q, qp] = 1
                original = initial_density(unit, distribution)
                recovered_unit = permute_density(permute_density(original, full_executor), full_executor)
                largest_matrix_unit_recovery = max(largest_matrix_unit_recovery, norm(recovered_unit - original))
        rows.append({"eta": eta, "relative_trace_distance_a": a,
                     "predicted_a": eta * (1 - 1 / n), "success": probability,
                     "success_after_active_twirl": probability_average,
                     "success_after_global_twirl": float(np.dot(e, np.diag(global_output).real)),
                     "source_Q_offdiagonal_after_tracing_path": float(abs(q_marginal[0, 1]))})
    assert max(largest_global, largest_recovery, largest_saturation, largest_matrix_unit_recovery) < TOL

    # A participating reference gives an explicitly different process: a=0, p=1.
    # Its kappa cannot be set to zero. This only certifies a lower bound, not the diamond norm.
    invariant_input = initial_density(q1, product)
    used_reference_output = permute_density(invariant_input, reference_executor)
    used_reference_average = sum(permute_density(used_reference_output, pa) for pa in active_permutations) / n
    kappa_witness = float(np.sum(np.abs(np.diag(used_reference_output - used_reference_average))) / 2)
    reference_success = float(np.dot(e, np.diag(used_reference_output).real))
    assert abs(reference_success - 1) < TOL
    assert abs(kappa_witness - (1 - 1 / n)) < TOL
    return {"cyclic_group_order": n, "scope": "Finite C_N algebra probe, not continuous SO(2) field dynamics",
            "unitarity_residual": unitarity, "active_covariance_residual": covariance,
            "global_effect_invariance_residual": global_effect_residual,
            "active_averaged_effect_norm": float(np.max(average_effect)),
            "complete_isometry_residual": isometry, "complete_control_identity_residual": complete_gate,
            "arbitrary_reference_sample_residual": passive_reference_residual,
            "max_global_state_invariance_residual": largest_global,
            "max_complete_Q_recovery_residual": largest_recovery,
            "max_matrix_unit_recovery_residual": largest_matrix_unit_recovery,
            "max_a_bound_saturation_residual": largest_saturation,
            "correlation_samples": rows,
            "participating_reference_negative_control": {
                "initial_a": 0.0, "success": reference_success,
                "certified_kappa_lower_bound": kappa_witness,
                "not_claimed": "No exact diamond-norm value is inferred from this input witness."}}


def s_wave_group():
    qx = qy = (1 - math.cos(ALPHA)) / 2
    angular_vector = np.sqrt([qx, qy, 1 - qx - qy])
    projector = np.outer(angular_vector, angular_vector)
    energies = np.array([4.7, 5.0, 5.3])
    phases = [.37, math.pi / 2, 2.1]
    delta = math.pi / 8
    full = np.zeros((18, 18), dtype=complex)
    h = np.diag(np.repeat(2 * energies, 6))
    rows = []
    largest_unitarity = 0.0
    largest_formula = 0.0
    for j, (energy, phase) in enumerate(zip(energies, phases)):
        scattering = np.eye(3) + (np.exp(2j * phase) - 1) * projector
        largest_unitarity = max(largest_unitarity, norm(scattering.conj().T @ scattering - np.eye(3)))
        probabilities = np.abs(scattering[:, 0]) ** 2
        predicted = 4 * math.sin(phase) ** 2 * qx * qy
        largest_formula = max(largest_formula, abs(float(probabilities[1]) - predicted))
        assert abs(float(sum(probabilities)) - 1) < TOL
        block = np.zeros((6, 6), dtype=complex)
        block[:3, :3] = np.eye(3)
        block[3:, 3:] = scattering
        full[6*j:6*j+6, 6*j:6*j+6] = np.exp(-2j * delta * energy) * block
        rows.append({"one_body_energy": float(energy), "phase_shift": phase,
                     "input_x_output_probabilities_x_y_remainder": [float(x) for x in probabilities],
                     "expected_y_probability": predicted})
    assert max(largest_unitarity, largest_formula) < TOL
    weights = np.array([1, 1j, -.7 + .2j], dtype=complex)
    weights /= np.linalg.norm(weights)
    j_in = np.zeros((18, 2), dtype=complex)
    for energy_index, amplitude in enumerate(weights):
        for source in range(2):
            j_in[energy_index * 6 + source * 3, source] = amplitude
    j_out = full @ j_in
    complete_isometry = norm(j_out.conj().T @ j_out - np.eye(2))
    conservation = norm(full @ h - h @ full)
    assert max(complete_isometry, conservation) < TOL
    q1_vector = j_out[:, 1]
    py = float(sum(abs(q1_vector[6*j + 4]) ** 2 for j in range(3)))
    bound = 4 * qx * qy
    assert py <= bound + TOL
    coherent_density = np.outer(q1_vector, q1_vector.conj())
    dephased = np.zeros_like(coherent_density)
    for j in range(3):
        dephased[6*j:6*j+6, 6*j:6*j+6] = coherent_density[6*j:6*j+6, 6*j:6*j+6]
    discarded_coherence = norm(coherent_density - dephased)
    assert discarded_coherence > .5
    return {"angular_fractions": {"qx": qx, "qy": qy, "remainder": 1 - qx - qy},
            "rank_one_projector_residual": norm(projector @ projector - projector),
            "max_unitarity_residual": largest_unitarity,
            "max_probability_formula_residual": largest_formula,
            "per_energy_samples": rows,
            "all_energy_complete_isometry_residual": complete_isometry,
            "energy_commutator_residual": conservation,
            "coherent_energy_packet_y_probability": py,
            "universal_y_probability_bound": bound,
            "task_error_floor": 1 - bound,
            "higher_mode_remainder_needed_for_99_percent": math.sqrt(.99) - 2 * math.sqrt(qx * qy),
            "norm_of_energy_coherence_lost_if_dephased": discarded_coherence,
            "scope": "Three angular sectors and three energy labels only calibrate full-output algebra."}


def compare(left, right, path="root"):
    if isinstance(left, dict):
        assert left.keys() == right.keys(), path
        for key in left:
            compare(left[key], right[key], path + "." + key)
    elif isinstance(left, list):
        assert len(left) == len(right), path
        for i, value in enumerate(left):
            compare(value, right[i], path + f"[{i}]")
    elif isinstance(left, (int, float)) and not isinstance(left, bool):
        assert abs(left - right) <= 1e-10 * max(1, abs(left), abs(right)), path
    else:
        assert left == right, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = {"round": 1124, "groups_passed": 3,
              "continuous_cap": geometry_group(),
              "finite_relational_reference": finite_reference_group(),
              "complete_s_wave": s_wave_group()}
    path = HERE / "results.json"
    if args.write_results:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        status = "written"
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
        status = "matched read-only"
    print(json.dumps({"round": 1124, "groups_passed": 3, "saved_results": status,
                      "continuous_no_resource_error_floor": 1 - ALPHA / math.pi,
                      "finite_reference_saturation_residual": result["finite_relational_reference"]["max_a_bound_saturation_residual"]},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
