"""Round 1117: finite checks of a position-readout bridge, not a detector model.

Default: recompute and compare saved results without writing.
--write-results: write the explicitly requested result artifact.
All units in the checks are c=hbar=1. General statements are proved in proof.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESEARCH = STAGE.parent
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
MASSES = np.array([2.0, 3.0])
ELL = 0.25


def norm(a):
    return float(np.linalg.norm(a, ord=2))


def energy(p, masses=MASSES):
    return np.sqrt(np.dot(p, p) + np.asarray(masses) ** 2)


def unitary(p, t):
    return np.diag(np.exp(-1j * t * energy(p)))


def covariant_readout_unitary(p):
    z = ELL * float(p[0])
    return np.cos(z) * I2 - 1j * np.sin(z) * SX


def tensor_packet(n, a=0.2, center=(0.0, 0.0, 0.0)):
    """Normal H^1 packet: phi=a^-3/2 product cos(pi q_i/(2a)), |q_i|<=a."""
    nodes, weights = leggauss(n)
    q1 = a * nodes
    w1 = weights * np.cos(np.pi * nodes / 2) ** 2
    grids = np.meshgrid(q1, q1, q1, indexing="ij")
    wg = np.meshgrid(w1, w1, w1, indexing="ij")
    p = np.stack(grids, axis=-1).reshape(-1, 3) + np.array(center)
    w = (wg[0] * wg[1] * wg[2]).reshape(-1)
    raw_norm = float(w.sum())
    w /= raw_norm
    return p, w, raw_norm


def packet_response(n, times, a=0.2):
    p, w, raw_norm = tensor_packet(n, a)
    es = np.sqrt(np.sum(p * p, axis=1)[:, None] + MASSES[None, :] ** 2)
    gap = es[:, 1] - es[:, 0]
    velocity = 0.5 * np.sum(p[:, 0, None] / es, axis=1)
    out = []
    ket = np.array([1.0, 1j]) / np.sqrt(2)
    for t in times:
        state = np.exp(-1j * t * es) * ket[None, :]
        # Independent 2x2 expectation, not the sin expression used for comparison.
        sx_mean = np.einsum("ni,ij,nj->n", state.conj(), SX, state).real
        full_mean = float(np.dot(w, t * velocity + ELL * sx_mean))
        sin_mean = float(ELL * np.dot(w, np.sin(gap * t)))
        out.append((full_mean, sin_mean))
    return p, w, gap, raw_norm, out


def calculate():
    groups = []
    # 1. A genuine translated PVM, or its centered probability smearing.
    p = np.array([0.37, -0.22, 0.11])
    W = covariant_readout_unitary(p)
    dW = -1j * ELL * SX @ W
    connection = 1j * W.conj().T @ dW
    A = ELL * SX
    H = np.diag(energy(p))
    Adiag = np.diag([0.2, -0.15])
    phase = np.exp(-1j * np.dot(p, [0.4, -0.1, 0.2])) * I2
    residual = max(norm(W.conj().T @ W - I2), norm(connection - A),
                   norm(W @ phase - phase @ W))
    assert residual < 1e-12
    assert norm(H @ Adiag - Adiag @ H) == 0.0
    mixed = norm(H @ A - A @ H)
    assert mixed > 0.2
    groups.append({
        "id": "translation_covariance_does_not_remove_mass_mixing",
        "status": "PASS",
        "unitary_connection_and_translation_residual": residual,
        "mixed_connection_commutator_norm": mixed,
        "mass_preserving_connection_commutator_norm": 0.0,
        "centered_smearing_width": 0.4,
        "smearing_preserves_first_moment": True,
        "causal_local_detector_implementation_claimed": False,
    })

    # 2. Derivative identity tested as a complete two-column operator identity.
    eps = 2e-5
    errs = []
    finite_time = []
    for point in ([0.17, -0.28, 0.09], [0.6, 0.2, -0.1]):
        point = np.array(point)
        for t in (0.3, 1.7):
            U = unitary(point, t)
            for j in range(3):
                h = np.zeros(3)
                h[j] = eps
                derivative_U = (unitary(point + h, t) - unitary(point - h, t)) / (2 * eps)
                Aj = A if j == 0 else np.zeros((2, 2), complex)
                v0 = np.diag(point[j] / energy(point))
                lhs_without_derivative = 1j * U.conj().T @ derivative_U + U.conj().T @ Aj @ U - Aj
                rhs_without_derivative = t * v0 + U.conj().T @ Aj @ U - Aj
                errs.append(norm(lhs_without_derivative - rhs_without_derivative))
                finite_time.append(norm(U.conj().T @ Aj @ U - Aj))
    assert max(errs) < 2e-8
    assert max(finite_time) <= 2 * ELL + 1e-13
    groups.append({
        "id": "complete_first_moment_free_displacement_identity",
        "status": "PASS",
        "central_difference_step": eps,
        "max_independent_derivative_residual": max(errs),
        "max_connection_displacement_in_samples": max(finite_time),
        "analytic_connection_displacement_bound": 2 * ELL,
        "arbitrary_passive_reference_via_operator_identity": True,
    })

    # 3. Normal finite-energy, finite-position-second-moment packet.
    times = [0.0, np.pi / 2, np.pi, 10.0, 50.0]
    p24, w24, gap24, raw24, out24 = packet_response(24, times)
    p, w, gap, raw, out = packet_response(40, times)
    convergence = max(abs(a[0] - b[0]) for a, b in zip(out24, out))
    identity_resid = max(abs(a - b) for a, b in out)
    a = 0.2
    gap_min = (MASSES[1] ** 2 - MASSES[0] ** 2) / (
        np.sqrt(MASSES[1] ** 2 + 3 * a * a)
        + np.sqrt(MASSES[0] ** 2 + 3 * a * a))
    assert convergence < 2e-12
    assert identity_resid < 2e-13
    assert abs(raw - 1) < 1e-12 and abs(raw24 - 1) < 1e-12
    assert out[1][0] > 0.2499
    rows = []
    for t, (value, _) in zip(times, out):
        center = ELL * np.sin(t)
        bound = ELL * abs(t) * (1 - gap_min)
        assert abs(value - center) <= bound + 2e-13
        rows.append({"time": t, "measured_first_moment": value,
                     "center_fiber_formula": float(center),
                     "normal_packet_to_center_bound": float(bound)})
    groups.append({
        "id": "normal_packet_zero_momentum_centroid_oscillation",
        "status": "PASS",
        "support_half_width_each_axis": a,
        "raw_packet_norm": raw,
        "quadrature_orders": [24, 40],
        "quadrature_convergence_difference": convergence,
        "matrix_expectation_vs_sine_identity_residual": identity_resid,
        "mean_momentum_norm": float(np.linalg.norm(w @ p)),
        "finite_canonical_position_variance_each_axis": float(np.pi ** 2 / (4 * a * a)),
        "exact_gap_lower_bound_on_support": float(gap_min),
        "rows": rows,
        "normal_packet_is_not_a_sharp_momentum_eigenstate": True,
        "centered_readout_noise_changes_resolution_not_mean": True,
    })

    # 4. Bounded readout offsets need not invalidate a finite averaged inertia test.
    coarse = []
    for t, (value, _) in zip(times[1:], out[1:]):
        drift = value / t
        bound = 2 * ELL / t
        assert abs(drift) <= bound + 1e-13
        coarse.append({"time": t, "readout_average_speed": drift,
                       "arbitrary_state_connection_bound": bound})
    T = 50.0
    position_bias_each_end = 0.005
    mean_momentum_radius = 0.004
    expected_total = (2 * ELL + 2 * position_bias_each_end) / T + mean_momentum_radius / 2
    window_max = 60.0
    assert T <= window_max
    assert expected_total < 0.013
    groups.append({
        "id": "finite_window_weaker_operational_branch",
        "status": "PASS",
        "rows": coarse,
        "conditional_example": {
            "declared_free_window_max": window_max, "chosen_time": T,
            "assumed_endpoint_mean_bias": position_bias_each_end,
            "assumed_mean_momentum_radius": mean_momentum_radius,
            "fixed_mass": 2.0, "total_speed_error_bound": expected_total,
            "accuracy_budget": 0.013,
            "is_a_detector_or_preparation_certificate": False,
        },
        "no_infinite_time_limit_required": True,
    })

    # 5. Same M enters exact spectral v, p=gamma M v and energy; means need care.
    errs = []
    derivative_errors = []
    for mass in MASSES:
        for point in ([0.0, 0.0, 0.0], [0.4, -0.2, 0.1], [1.0, 0.6, -0.3]):
            point = np.array(point)
            E = np.sqrt(mass * mass + point @ point)
            v = point / E
            gamma = 1.0 / np.sqrt(1 - v @ v)
            errs.extend([abs(E - gamma * mass), float(np.linalg.norm(point - gamma * mass * v))])
            J = np.eye(3) / E - np.outer(point, point) / E ** 3
            epsj = 1e-5
            Jfd = np.empty((3, 3))
            for j in range(3):
                d = np.eye(3)[j] * epsj
                vp = (point + d) / np.sqrt(mass * mass + np.dot(point + d, point + d))
                vm = (point - d) / np.sqrt(mass * mass + np.dot(point - d, point - d))
                Jfd[:, j] = (vp - vm) / (2 * epsj)
            derivative_errors.append(norm(Jfd - J))
            assert np.linalg.eigvalsh(J).min() > 0
            assert norm(J) <= 1.0 / mass + 1e-13
    small_p = 0.05
    mixture_v = float(np.mean(small_p / np.sqrt(MASSES ** 2 + small_p ** 2)))
    average_mass_v = float(small_p / np.sqrt(np.mean(MASSES) ** 2 + small_p ** 2))
    assert max(errs) < 2e-14
    assert max(derivative_errors) < 2e-10
    assert abs(mixture_v - average_mass_v) > 0.0008
    groups.append({
        "id": "same_mass_in_velocity_and_low_speed_response",
        "status": "PASS",
        "max_exact_spectral_relation_residual": max(errs),
        "max_response_jacobian_finite_difference_residual": max(derivative_errors),
        "equal_weight_mass_mixture_speed": mixture_v,
        "speed_using_mean_mass_instead": average_mass_v,
        "incorrect_mean_mass_speed_gap": mixture_v - average_mass_v,
        "ratio_of_means_not_asserted": True,
    })

    # 6. Finite impulse response: algebra is conditional on an actual conservative kick.
    mass = 2.0
    kick = np.array([0.02, 0.0, 0.0])
    p, w, _ = tensor_packet(28, a=0.002)
    before = p / np.sqrt(mass * mass + np.sum(p * p, axis=1))[:, None]
    p_after = p + kick
    after = p_after / np.sqrt(mass * mass + np.sum(p_after * p_after, axis=1))[:, None]
    response = w @ (after - before)
    central_response = kick / np.sqrt(mass * mass + kick @ kick)
    linear = kick / mass
    finite_impulse_bound = float(np.linalg.norm(kick) ** 3 / (2 * mass ** 3))
    packet_bound = float(2 * np.dot(w, np.linalg.norm(p, axis=1)) / mass)
    assert np.linalg.norm(central_response - linear) <= finite_impulse_bound + 1e-14
    assert np.linalg.norm(response - central_response) <= packet_bound + 1e-14
    groups.append({
        "id": "finite_conservative_impulse_inertia_protocol",
        "status": "PASS",
        "mass": mass, "declared_momentum_transfer": kick.tolist(),
        "packet_mean_velocity_change": response.tolist(),
        "central_exact_change": central_response.tolist(),
        "low_speed_linear_change": linear.tolist(),
        "central_finite_impulse_error": float(np.linalg.norm(central_response - linear)),
        "central_finite_impulse_bound": finite_impulse_bound,
        "normal_packet_error": float(np.linalg.norm(response - central_response)),
        "normal_packet_bound": packet_bound,
        "actual_impulse_source_and_recoil_not_constructed": True,
    })

    sources = [
        "archive_1106_/1108/proof.md",
        "archive_935_955/research_note_941.md",
        "archive_956_989/research_note_971.md",
        "archive_1009_1043/research_note_1038.md",
        "archive_1106_/1106/proof.md",
        "archive_1106_/1111/proof.md",
        "archive_1106_/1115/proof.md",
    ]
    hashes = {p: hashlib.sha256((RESEARCH / p).read_bytes()).hexdigest() for p in sources}
    return {
        "round": 1117,
        "status": "PASS",
        "numerical_check_groups": 6,
        "scientific_calibration_increment": 0,
        "new_adopted_axioms": 0,
        "units_in_checks": "c=hbar=1",
        "scope": {
            "conditional_operational_position_bridge": True,
            "PT_source_or_Lorentz_classification_reproved": False,
            "ideal_envelope_is_one_finite_F_object": False,
            "new_local_detector_or_full_conservative_readout_built": False,
            "Lorentz_covariant_sharp_position_PVM_required": False,
            "mass_preserving_readout_required_for_effective_inertia": False,
            "all_UV_or_infinite_free_window_required": False,
            "full_cognitive_axiom_countermodel_claimed": False,
        },
        "groups": groups,
        "source_sha256": hashes,
    }


def compare(a, b, trail=""):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys(), trail
        for k in a:
            compare(a[k], b[k], trail + "/" + k)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), trail
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, trail + "/" + str(i))
    elif isinstance(a, float):
        assert isinstance(b, (float, int)) and np.isclose(a, b, rtol=2e-8, atol=2e-11), (trail, a, b)
    else:
        assert a == b, (trail, a, b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-results", action="store_true")
    args = ap.parse_args()
    result = calculate()
    path = HERE / "results.json"
    if args.write_results:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print("WROTE_RESULTS")
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
        print("SAVED_RESULTS_MATCH")
    print("PASS: 6 finite position/inertia bridge groups")
    print("Normal packets; exact instantaneous and bounded finite-window branches separated")
    print("Scientific calibration 0; no new detector, universal localization or whole-theory countermodel")


if __name__ == "__main__":
    main()
