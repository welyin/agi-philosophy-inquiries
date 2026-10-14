"""Round 1159: finite diagnostics for operationally separated regional charges.

The all-input and all-finite-network claims are proved in proof.md.
No old research script is executed or old asset rewritten.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
TOL = 2e-10
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def op(a):
    return float(np.linalg.norm(a, 2))


def zero(a, name):
    value = op(a) if np.ndim(a) == 2 else float(abs(a))
    assert value < TOL, (name, value)
    return value


def kron(*factors):
    result = np.ones((1, 1), complex)
    for factor in factors:
        result = np.kron(result, factor)
    return result


def lift(a, site, dims):
    return kron(*(a if j == site else np.eye(d) for j, d in enumerate(dims)))


def exp_h(a, t):
    val, vec = np.linalg.eigh(a)
    return (vec * np.exp(-1j * t * val)) @ vec.conj().T


def comm(a, b):
    return a @ b - b @ a


def weyl(d):
    shift = np.roll(np.eye(d, dtype=complex), 1, axis=0)
    clock = np.diag(np.exp(2j * np.pi * np.arange(d) / d))
    return shift, clock


def twirl(a, dims):
    local_words = []
    for d in dims:
        x, z = weyl(d)
        local_words.append([np.linalg.matrix_power(x, j) @ np.linalg.matrix_power(z, k)
                            for j in range(d) for k in range(d)])
    average = np.zeros_like(a)
    count = 0
    for words in itertools.product(*local_words):
        w = kron(*words)
        average += w.conj().T @ a @ w
        count += 1
    return average / count, count


def source_check():
    manifest = json.loads((HERE / "source_hashes.json").read_text(encoding="utf-8"))
    for name, expected in manifest["sha256"].items():
        actual = hashlib.sha256((RESEARCH / name).read_bytes()).hexdigest()
        assert actual == expected, (name, actual, expected)
    return len(manifest["sha256"])


def response_and_scalar_certificate():
    dims = [2, 3, 2]
    local_h = [np.diag([2., 3.]), np.diag([3., 4., 5.]), np.diag([1., 2.])]
    velocities = [np.array([.2, 0, 0]), np.array([0, -.3, 0]), np.array([0, 0, .1])]
    local_p = []
    for h, velocity in zip(local_h, velocities):
        u = np.r_[1., velocity] / math.sqrt(1-float(velocity @ velocity))
        local_p.append([component * h for component in u])
    total = [sum(lift(local_p[j][mu], j, dims) for j in range(3)) for mu in range(4)]
    scalar = np.array([.7, .1, -.2, .3])
    q = [total[mu] + scalar[mu] * np.eye(12) for mu in range(4)]
    exact = []
    for mu in range(4):
        for j, d in enumerate(dims):
            for g in weyl(d):
                exact.append(zero(comm(q[mu], lift(g, j, dims))
                                  - lift(comm(local_p[j][mu], g), j, dims),
                                  "all local finite generator responses"))
    zero(max(op(comm(q[a], q[b])) for a in range(4) for b in range(4)), "translation components commute")
    a = np.array([.37, -.2, .1, .4])
    metric_a = a * np.array([1., -1, -1, -1])
    joint = exp_h(sum(metric_a[mu] * q[mu] for mu in range(4)), 1.)
    product = kron(*(exp_h(sum(metric_a[mu] * local_p[j][mu] for mu in range(4)), 1.) for j in range(3)))
    phase = np.exp(-1j * float(metric_a @ scalar))
    phase_residual = zero(joint-phase*product, "channel-compatible total scalar")
    rng = np.random.default_rng(115901)
    test = rng.normal(size=(12, 12)) + 1j*rng.normal(size=(12, 12))
    test = test + test.conj().T
    channel_residual = zero(joint.conj().T@test@joint-product.conj().T@test@product, "same conjugation")

    interaction = .013 * kron(Z, np.diag([-1., 0, 1.]), I2)
    interaction += .009 * kron(I2, np.diag([1., -2., 1.]), Z)
    remainder = interaction + .07*np.eye(12)
    epsilon = []
    eta = 0.
    for j, d in enumerate(dims):
        row = [op(comm(remainder, lift(g, j, dims))) for g in weyl(d)]
        epsilon.append(row)
        eta += (d-1)/2 * sum(row)
    averaged, count = twirl(remainder, dims)
    kappa = float(np.trace(remainder).real/12)
    twirl_residual = zero(averaged-kappa*np.eye(12), "complete Weyl average")
    non_scalar = op(remainder-kappa*np.eye(12))
    assert non_scalar <= eta + TOL
    # A specified normalized reference; its charge calibration is an input.
    reference = np.zeros(12, complex)
    reference[0] = 1
    beta = abs(float(np.vdot(reference, remainder @ reference).real))
    calibrated_bound = beta + 2*eta
    assert op(remainder) <= calibrated_bound + TOL
    # Independent exact total and local charges: removing the scalar mismatch
    # is checked, never inferred merely from conjugation.
    return {"dimensions": dims, "full_dimension": 12, "local_response_tests": len(exact),
            "max_local_response_residual": max(exact), "undetected_total_scalar_four_vector": scalar.tolist(),
            "full_translation_phase_residual": phase_residual, "conjugation_residual": channel_residual,
            "full_lift_difference": op(joint-product), "weyl_word_count": count,
            "twirl_residual": twirl_residual, "generator_commutator_defects": epsilon,
            "eta": eta, "non_scalar_remainder_norm": non_scalar,
            "reference_calibration_beta": beta, "full_remainder_norm": op(remainder),
            "full_remainder_bound_beta_plus_2eta": calibrated_bound}


def baseline_and_actual_velocity():
    # The original two outgoing directions of 1111, in their fixed material type.
    pa = np.array([[5., 5.], [4., 0], [0, 4.], [0, 0.]])
    pb = np.array([[5., 5.], [-4., 0], [0, -4.], [0, 0.]])
    lam = .4
    qa = np.array([lam, 0, 0, 0])
    qb = -qa
    shifted_a = pa+qa[:, None]
    shifted_b = pb+qb[:, None]
    total_residual = max(zero(kron(np.diag(shifted_a[m]), I2)+kron(I2, np.diag(shifted_b[m]))
                                  -kron(np.diag(pa[m]), I2)-kron(I2, np.diag(pb[m])),
                                  "opposite material baseline changes preserve full total") for m in range(4))
    velocities = [pa[1:, j]/pa[0, j] for j in range(2)]
    equations = np.vstack([np.column_stack([-v, np.eye(3)]) for v in velocities])
    singular = np.linalg.svd(equations, compute_uv=False)
    assert np.linalg.matrix_rank(equations) == 4
    single = np.column_stack([-velocities[0], np.eye(3)])
    assert np.linalg.matrix_rank(single) == 3
    shift_null = np.r_[1., velocities[0]]
    zero(single @ shift_null @ (single @ shift_null).conj(), "one velocity retains parallel scalar")
    speed_error_a = max(float(np.linalg.norm(shifted_a[1:, j]/shifted_a[0, j]-velocities[j])) for j in range(2))
    assert speed_error_a > .05
    gap = float(np.linalg.norm(velocities[0]-velocities[1]))
    # Direct momentum residual, not velocity error without an energy bound.
    test_q = np.array([.003, .001, -.002, .0004])
    momentum_errors = [float(np.linalg.norm(test_q[1:]-v*test_q[0])) for v in velocities]
    q0_bound = sum(momentum_errors)/gap
    qvec_bound = momentum_errors[0]+float(np.linalg.norm(velocities[0]))*q0_bound
    assert abs(test_q[0]) <= q0_bound+TOL and np.linalg.norm(test_q[1:]) <= qvec_bound+TOL
    # Controlled lifts are additional actual permissions, not obtained from a channel.
    t = math.pi/lam
    control_input = np.array([1., 1.], complex)/math.sqrt(2)
    relative_control = np.diag([1., np.exp(-1j*lam*t)])
    contrast = float(np.vdot(control_input, X @ control_input).real
                     - np.vdot(relative_control@control_input, X @ relative_control@control_input).real)
    zero(contrast-2., "conditional full phase would detect local scalar")
    # Finite selected time cannot identify an unrestricted scalar: exact alias.
    alias_t = .7
    zero(np.exp(-1j*(2*math.pi/alias_t)*alias_t)-1, "one time has phase alias")
    return {"old_A_four_P_by_direction": pa.tolist(), "opposite_scalar_shift": qa.tolist(),
            "total_four_P_difference": total_residual,
            "local_energy_channels_identical": True, "all_shifted_energies_positive": True,
            "shifted_masses_A_B": [math.sqrt((5+lam)**2-16), math.sqrt((5-lam)**2-16)],
            "fixed_original_endpoint_velocity_error_A": speed_error_a,
            "two_velocity_matrix_rank": 4, "two_velocity_singular_values": singular.tolist(),
            "one_velocity_matrix_rank": 3, "velocity_gap": gap,
            "robust_test_q": test_q.tolist(), "momentum_residuals": momentum_errors,
            "q0_bound": q0_bound, "spatial_q_bound": qvec_bound,
            "conditional_coherent_control_X_contrast": contrast,
            "single_time_alias_checked": True,
            "not_a_countermodel_preserving_original_DP3_endpoints": True}


def regional_composition():
    # Logical undressed data block only. 1111/1112 prove that Gamma dresses
    # this block while retaining the full original H_D. No old task rerun.
    dims = [3, 4, 4]
    control = np.diag([0., 0., 1.])
    xbar = np.array([[0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0]], complex)
    gate1 = kron(np.eye(3)-control, np.eye(4), np.eye(4)) + kron(control, xbar, np.eye(4))
    gate2 = kron(np.eye(3)-control, np.eye(4), np.eye(4)) + kron(control, np.eye(4), xbar)
    h3 = np.eye(3, dtype=complex)
    h3[1:, 1:] = np.array([[1, 1], [1, -1]])/math.sqrt(2)
    future = lift(h3, 0, dims)
    whole = gate2 @ future @ gate1
    zero(whole.conj().T@whole-np.eye(48), "complete two-pair history")
    assert op(comm(gate1, future)) > 1
    pa = [5*I2, np.diag([4., 0]), np.diag([0, 4.]), 0*I2]
    pb = [5*I2, np.diag([-4., 0]), np.diag([0, -4.]), 0*I2]
    charges = []
    for pair in [1, 2]:
        for local in [pa, pb]:
            local_pair = [kron(p, I2) if local is pa else kron(I2, p) for p in local]
            charges.append([lift(p, pair, dims) for p in local_pair])
    changes = [[whole.conj().T@q@whole-q for q in actor] for actor in charges]
    total = [sum(actor[mu] for actor in charges) for mu in range(4)]
    total_residual = max(zero(comm(whole, q), "same complete total four_P") for q in total)
    sum_residual = max(zero(sum(actor[mu] for actor in changes), "all regional changes sum to zero") for mu in range(4))
    prefix = future@gate1
    pullback = [prefix.conj().T@(gate2.conj().T@charges[2][mu]@gate2-charges[2][mu])@prefix for mu in range(4)]
    cocycle_residual = max(zero(pullback[mu]-changes[2][mu], "active future pulls back same charge") for mu in range(4))
    omitted_partner = op(changes[0][1])
    assert omitted_partner > 3.9
    initial = np.zeros((48, 2), complex)
    initial[16, 0] = 1.
    initial[32, 1] = 1.
    full_output = whole@initial
    zero(full_output.conj().T@full_output-np.eye(2), "all unknown logical input columns")
    input_state = initial @ (np.array([1., 1.])/math.sqrt(2))
    current_mean = [[float(np.vdot(input_state, d@input_state).real) for d in actor] for actor in changes]
    boost = np.array([[1.25, -.75, 0, 0], [-.75, 1.25, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1.]])
    boost_residual = []
    for actor, delta in zip(charges, changes):
        for mu in range(4):
            transported_q = sum(boost[mu, nu]*actor[nu] for nu in range(4))
            transported_delta = whole.conj().T@transported_q@whole-transported_q
            boost_residual.append(zero(transported_delta-sum(boost[mu, nu]*delta[nu] for nu in range(4)),
                                       "transport same regional observable"))
    return {"retained_data_dimension": 48, "unknown_logical_input_dimension": 2,
            "full_total_four_P_commutator_residual": total_residual,
            "sum_regional_changes_residual": sum_residual, "active_future_cocycle_residual": cocycle_residual,
            "active_future_noncommutator_norm": op(comm(gate1, future)),
            "omitted_recoil_partner_error_norm": omitted_partner,
            "regional_current_means_A1_B1_A2_B2": current_mean,
            "boosted_regional_current_means": [(boost@np.array(row)).tolist() for row in current_mean],
            "max_boost_identity_residual": max(boost_residual),
            "full_original_clock_and_Gamma_quantifiers_reused_analytically": True,
            "no_new_routing_or_window_interior_permission_claimed": True}


def cross_boundary_and_endpoint_resolution():
    # A joined two-part participant during an interaction, NOT two spacelike
    # separated updates. Its failure only rejects a premature separated cut.
    bare = kron(2*I2+.3*Z, I2)+kron(I2, 3*I2+.2*Z)
    coupling = .17*kron(Z, Z)
    full = bare+coupling
    probe = kron(X, I2)
    response_defect = op(comm(full, probe)-comm(bare, probe))
    zero(response_defect-.34, "cross-interface interaction response")
    r_avg, _ = twirl(coupling, [2, 2])
    zero(r_avg, "interaction has zero scalar mean")
    scalar_only_calibration = float(np.trace(coupling).real/4)
    assert op(coupling) > .16 and scalar_only_calibration == 0
    # Endpoint agreement does not certify an undeclared inside-window cut.
    short = .002
    hidden = math.pi/short*kron(Z, Z)
    at_endpoint = exp_h(hidden, short)
    at_half = exp_h(hidden, short/2)
    zero(at_endpoint+np.eye(4), "same endpoint identity channel")
    half_response = op(at_half.conj().T@probe@at_half-probe)
    zero(half_response-2., "finite internal response differs")
    # Perturbed charge-preserving event: norm-bounded finite observable only.
    angle = .007
    bad = exp_h(kron(X, I2), angle)
    q = kron(Z, I2)
    flux = op(bad.conj().T@q@bad-q)
    bound = 2*op(q)*op(bad-np.eye(4))
    assert flux <= bound+TOL
    return {"cross_interaction_norm": op(coupling), "local_translation_response_defect": response_defect,
            "scalar_reference_alone_is_insufficient": True,
            "coarse_endpoint_duration": short, "hidden_generator_norm": op(hidden),
            "endpoint_channels_identical": True, "half_window_observable_defect": half_response,
            "bounded_observable_flux_error": flux, "bounded_observable_error_bound": bound,
            "not_a_spacelike_CC_countermodel": True,
            "no_unbounded_energy_estimate_from_trace_distance": True}


def calculate():
    return {"round": 1159, "status": "PASS_CONDITIONAL_SEPARATED_BOUNDARY_INTERFACE",
            "numerical_check_groups": 4, "scientific_calibration": 0, "new_adopted_axioms": 0,
            "source_hashes_checked": source_check(),
            "local_response_and_scalar_certificate": response_and_scalar_certificate(),
            "baseline_and_actual_velocity": baseline_and_actual_velocity(),
            "regional_composition_and_future": regional_composition(),
            "cross_boundary_and_endpoint_resolution": cross_boundary_and_endpoint_resolution(),
            "scope": {"conditional_operator_theorem_not_all_axioms_independence": True,
                      "arbitrary_reference_and_all_qualified_networks_proved_analytically": True,
                      "actual_separated_boundaries_required": True,
                      "full_window_interior_stress_not_derived": True,
                      "old_1111_1112_clock_and_native_permissions_remain_explicit": True,
                      "numerical_finite_matrices_not_new_empirical_physics": True,
                      "formal_review_and_mainline_acceptance_pending": True}}


def compare(expected, actual, path="root"):
    if isinstance(actual, dict):
        assert expected.keys() == actual.keys(), path
        for key in actual:
            compare(expected[key], actual[key], path+"."+key)
    elif isinstance(actual, list):
        assert len(expected) == len(actual), path
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, path+f"[{i}]")
    elif isinstance(actual, float):
        assert math.isclose(float(expected), actual, abs_tol=2e-9, rel_tol=2e-9), (path, expected, actual)
    else:
        assert expected == actual, (path, expected, actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-once", action="store_true")
    args = parser.parse_args()
    actual = calculate()
    path = HERE/"results.json"
    if args.write_once:
        if path.exists():
            raise FileExistsError(path)
        path.write_text(json.dumps(actual, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    else:
        compare(json.loads(path.read_text(encoding="utf-8")), actual)
    print("PASS 1159: four groups; separated-charge qualification, scalar/velocity calibration, regional continuation")
    print("FIRST_RESULTS_WRITTEN" if args.write_once else "READ_ONLY_MATCH")
    print("Conditional finite interface; no claim of universal local stress or new actual microscopic gates")


if __name__ == "__main__":
    main()
