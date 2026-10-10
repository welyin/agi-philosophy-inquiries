"""1055: finite-mass recoil and the same optical task, without a motion grid.

Fraction identities certify the continuous bounds. Quadratures below check
Gaussian moments, the exact forced propagator and source/parity identities;
they do not simulate the continuous massive optical process.
Default is read-only. --write creates the first result exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / "results.json"
HISTORY = (
    "archive_370_428/research_note_382.md",
    "archive_370_428/research_note_383.md",
    "archive_370_428/research_note_384.md",
    "archive_370_428/research_note_386.md",
    "archive_370_428/research_note_412.md",
    "archive_370_428/research_note_425.md",
    "archive_467_530/research_note_488.md",
    "archive_467_530/research_note_522.md",
    "archive_467_530/research_note_523.md",
    "archive_935_955/research_note_941.md",
    "archive_956_989/research_note_978.md",
    "archive_956_989/research_note_979.md",
    "archive_1009_1043/research_note_1038.md",
    "archive_1046_/research_note_1046.md",
    "archive_1046_/1046/proof.md",
    "archive_1046_/1046/photon_direction_receiver_results.json",
    "archive_1046_/1046/research_round_1046_checks.json",
    "archive_1046_/research_note_1051.md",
    "archive_1046_/research_note_1054.md",
    "archive_1046_/1054/proof.md",
    "archive_1046_/1054/results.json",
    "archive_1046_/1054/research_round_1054_checks.json",
    "archive_1046_/_admission/displacement_after1054/selection.md",
    "archive_1046_/_admission/displacement_after1054/pre_review.md",
    "archive_1046_/_admission/displacement_after1054/mainline_admission.md",
)
PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                  [[1, 0], [0, -1]]], complex)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational_certificate():
    m, sigma, t, radius = F(10**20), F(1, 10**6), F(1, 100), F(1, 10**5)
    g1, g2, h1 = 1/sigma, 1/sigma**2, F(2, 7)
    a0 = F(5, 4)*g2+3*g1+2+4*t*g1+3*t+3*t*t+4*t
    a1 = F(5, 4)*g2+11*g1+19
    common_a = F(2*10**12)
    assert max(a0, a1) < common_a
    a0_dir = F(5, 4)*g2+7*g1+12
    assert a0_dir < common_a
    # sqrt(3)<2 and the natural free Gaussian broadening, without a reset.
    spread_upper = 2*sigma+2/(m*sigma)
    mean_radius_upper = F(3, 10**6)
    assert spread_upper < mean_radius_upper
    t1_rms_upper = radius+mean_radius_upper
    e0 = t*common_a/m
    entry = t*t*(common_a+3*common_a)/m
    kinetic_jac = 3*entry
    c = F(1401, 42875000)
    blur_jac = 6*c*mean_radius_upper
    jac_error = kinetic_jac+blur_jac
    old_kappa = F(589541, 39711168000000)
    old_ball_error = 6*c*radius
    lower_slope = old_kappa-old_ball_error-jac_error
    probability_error = 2*e0+8*t*t/F(49)*mean_radius_upper
    probability_lower = F(802069, 2206176000000)-probability_error
    direction_effect = 2*(common_a/m+h1*t1_rms_upper)
    direction_effect_simple = F(3, 400000)
    direction_gap = F(6889, 4000000)-2*direction_effect_simple
    assert jac_error == F(26247, 42875000000000) < F(1, 10**9)
    assert probability_error == F(11, 24500000000)
    assert lower_slope == F(7615215221, 620487000000000000) > F(1, 10**8)
    assert probability_lower > F(1, 3*10**6)
    assert direction_effect < direction_effect_simple
    assert direction_gap == F(6829, 4000000) > F(17, 10000)
    # Actual transport budget at the edge of the declared displacement ball.
    mu = m/2
    impulse_peak = F(3, 2)*mu*radius
    velocity_peak = impulse_peak/m  # one body's speed added by the force
    max_phase = F(3, 5)*mu*radius**2
    action_upper = 3*mu*radius
    values = {
        "mass_each": m, "sigma_each_coordinate": sigma, "transport_T": F(1),
        "position_read_t": t, "direction_read_t": F(1), "displacement_radius": radius,
        "A0_position_envelope": a0, "A1_envelope": a1,
        "A0_direction_envelope": a0_dir, "A_common_upper": common_a,
        "natural_spread_upper": spread_upper, "normal_mean_radius_upper": mean_radius_upper,
        "state_error_position_upper": e0, "C1_entry_upper": entry,
        "kinetic_jacobian_upper": kinetic_jac, "normal_average_jacobian_upper": blur_jac,
        "total_jacobian_error_upper": jac_error, "old_slope_lower": old_kappa,
        "old_ball_jacobian_variation": old_ball_error, "new_co_lipschitz_lower": lower_slope,
        "simple_co_lipschitz_lower": F(1, 10**8), "simple_inverse_lipschitz_upper": F(10**8),
        "probability_error_upper": probability_error,
        "probability_at_origin_lower": probability_lower,
        "direction_state_mass_error_upper": common_a/m,
        "direction_position_rms_upper": t1_rms_upper,
        "direction_effect_error_upper": direction_effect,
        "direction_effect_simple_upper": direction_effect_simple,
        "direction_antipodal_gap_lower": direction_gap,
        "peak_force_induced_momentum_upper": impulse_peak,
        "peak_individual_velocity_increment_upper": velocity_peak,
        "force_area_upper": action_upper, "controlled_phase_magnitude_upper": max_phase,
    }
    return {"exact": {k: str(v) for k, v in values.items()},
            "decimal": {k: float(v) for k, v in values.items()}}


def integrate_polynomial(coefficients):
    return sum((c/F(i+1) for i, c in enumerate(coefficients)), F(0))


def multiply(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def force_certificate():
    mu = F(10**20, 2)
    ds = [F(1, 10**6), F(-2, 10**6), F(1, 3*10**6)]
    force_area = sum(abs(x) for x in ds)*3*mu  # componentwise budget diagnostic
    phase = F(0)
    for d in ds:
        force = [6*mu*d, -12*mu*d]
        impulse = [F(0), 6*mu*d, -6*mu*d]
        assert integrate_polynomial(force) == 0
        assert integrate_polynomial(impulse)/mu == d
        phase += integrate_polynomial(multiply(impulse, impulse))/(2*mu)
    expected = F(3, 5)*mu*sum((d*d for d in ds), F(0))
    assert phase == expected
    # A separate moderate-mass calibration avoids meaningless large-phase float reduction.
    nodes, weights = np.polynomial.legendre.leggauss(12)
    s, weights = (nodes+1)/2, weights/2
    mu_c = 1.7
    d = np.array([.13, -.21, .07])
    impulses = 6*mu_c*s[:, None]*(1-s[:, None])*d
    ps = np.array([[.2, .5, -.1], [-.7, .3, .4], [0., 0., 0.]])
    phases_quad = np.sum(weights[:, None]*np.sum(
        (ps[None, :, :]+impulses[:, None, :])**2, axis=2)/(2*mu_c), axis=0)
    phases_formula = np.sum(ps*ps, axis=1)/(2*mu_c)+ps@d+3*mu_c*np.dot(d, d)/5
    amplitude_error = float(np.max(np.abs(np.exp(-1j*phases_quad)-np.exp(-1j*phases_formula))))
    assert amplitude_error < 2e-14
    # A following -d cancels displacement but retains the free propagation and phase.
    phase_c = 3*mu_c*np.dot(d, d)/5
    two_actual = np.exp(-2j*(np.sum(ps*ps, axis=1)/(2*mu_c)+phase_c))
    false_identity = float(np.max(np.abs(two_actual-1)))
    assert false_identity > .1
    return {"exact_task_displacement": [str(d) for d in ds],
            "exact_task_phase": str(phase), "exact_componentwise_force_area": str(force_area),
            "calibration_scope": "moderate mass characteristic momentum solution, not the massive optical evolution",
            "quadrature_vs_factorized_amplitude_max_error": amplitude_error,
            "opposite_force_composition_distance_from_identity_on_samples": false_identity,
            "coherent_control_phase_retained": True}


def gaussian_calibration():
    # Independent Gauss-Hermite integration in momentum, not a position-motion grid.
    x, weights = np.polynomial.hermite.hermgauss(12)
    weights = weights/np.sqrt(np.pi)
    z = np.sqrt(2)*x  # N(0,1)
    normal = float(weights.sum())
    moment2 = float(np.dot(weights, z*z))
    moment4 = float(np.dot(weights, z**4))
    p2_sigma_scaled = 3*moment2/4
    p4_sigma_scaled = (3*moment4+6*moment2**2)/16
    assert abs(normal-1) < 2e-15
    assert abs(p2_sigma_scaled-.75) < 3e-15
    assert abs(p4_sigma_scaled-15/16) < 3e-15
    sigma, mass = 1e-6, 1e20
    variance_after = sigma*sigma+1/(mass*mass*sigma*sigma)
    assert np.sqrt(3*variance_after) < 3e-6
    # A six-dimensional tensor quadrature independently checks the recoil K^2 identity.
    x, weights = np.polynomial.hermite.hermgauss(5)
    z, weights = np.sqrt(2)*x, weights/np.sqrt(np.pi)
    index = np.indices((5,)*6).reshape(6, -1).T
    zz = z[index]
    ww = np.prod(weights[index], axis=1)
    sigma_c, k = .7, np.array([.3, -.4, 1.2])
    ps, p = zz[:, :3]/(2*sigma_c), zz[:, 3:]/(2*sigma_c)
    kinetic = np.sum((ps-k)**2, axis=1)/4+np.sum((p+k/2)**2, axis=1)
    v = 1/(4*sigma_c*sigma_c)
    k2 = float(k@k)
    expected1 = 15*v/4+k2/2
    expected2 = 327*v*v/16+5*v*k2+k2*k2/4
    error = max(abs(float(ww@kinetic)-expected1), abs(float(ww@(kinetic**2))-expected2))
    assert error < 3e-13
    return {"normalization": normal, "sigma2_times_p2": p2_sigma_scaled,
            "sigma4_times_p4": p4_sigma_scaled,
            "task_relative_rms_after_transport": float(np.sqrt(3*variance_after)),
            "kinetic_K_moment_calibration_max_error": error,
            "kinetic_K_mean": expected1, "kinetic_K_second_moment": expected2,
            "moment_calibration_scope": "one fixed moderate-width Gaussian; task bounds use exact moments"}


def source_and_parity_calibration():
    n = 16
    nodes, weights = np.polynomial.legendre.leggauss(n)
    radial_k = 1.5+nodes/2
    radial_w = weights/2*(8/3)*np.sin(np.pi*(radial_k-1))**4
    radial_norm_error = abs(float(radial_w.sum())-1)
    assert radial_norm_error < 2e-14
    directions = np.vstack((np.eye(3), -np.eye(3)))
    ks = np.repeat(radial_k, 6)[:, None]*np.tile(directions, (n, 1))
    ww = np.repeat(radial_w/6, 6)
    dirs = np.tile(directions, (n, 1))
    projectors = np.eye(3)[None, :, :]-dirs[:, :, None]*dirs[:, None, :]
    h = np.sqrt(ww)[:, None, None]*projectors
    # Column j is the transverse wavepacket h_j. Sum over j/2 is normalized.
    h_flat = h.reshape(-1, 3)
    gram = h_flat.conj().T@h_flat
    gram_error = float(np.max(np.abs(gram-2*np.eye(3)/3)))
    assert gram_error < 2e-14
    vector_j = h.reshape(-1)/np.sqrt(2)
    assert abs(float(np.vdot(vector_j, vector_j).real)-1) < 2e-14
    # Source phases preserve normalization and correlate reference recoil with photon k.
    y, a = np.array([.2, -.1, .05]), np.array([.25, 0., 0.])
    source = np.exp(-1j*ks@(y+a))[:, None, None]*h/np.sqrt(2)
    assert abs(float(np.vdot(source, source).real)-1) < 2e-14
    p_cm, p_rel = np.array([.2, -.3, .4]), np.array([-.4, .1, .6])
    physical_x = p_cm/2+p_rel
    physical_y = p_cm/2-p_rel-ks
    total = physical_x+physical_y+ks
    momentum_error = float(np.max(np.abs(total-p_cm)))
    # Actual absorption transfers +k to X, preserving total momentum also at the vertex.
    total_after = physical_x+ks+physical_y
    momentum_error = max(momentum_error, float(np.max(np.abs(total_after-p_cm))))
    k_original = (np.sum(physical_x**2)+np.sum(physical_y**2, axis=1))/2
    k_reference = np.sum((p_cm-ks)**2, axis=1)/4+np.sum((p_rel+ks/2)**2, axis=1)
    kinetic_error = float(np.max(np.abs(k_original-k_reference)))
    assert max(momentum_error, kinetic_error) < 3e-14
    # Kernel-level parity; the full continuum identity is proved analytically.
    r = np.array([.17, -.11, .09])
    coupling = np.einsum("lij,jab->liab", projectors, PAULI)
    a_r = np.sqrt(ww)[:, None, None, None]*np.exp(1j*ks@r)[:, None, None, None]*coupling
    a_minus = np.sqrt(ww)[:, None, None, None]*np.exp(-1j*ks@r)[:, None, None, None]*coupling
    opposite = np.concatenate([np.arange(6*i+3, 6*i+6).tolist()+np.arange(6*i, 6*i+3).tolist()
                               for i in range(n)])
    parity_error = float(np.max(np.abs(-a_r[opposite]+a_minus)))
    assert parity_error < 2e-14
    # Spin/reference source preparation is an isometry; the old reference is not traced.
    jmap = np.kron(np.eye(2), source.reshape(-1, 1))
    iso_error = float(np.max(np.abs(jmap.conj().T@jmap-np.eye(2))))
    assert iso_error < 2e-14
    bell = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
    joint = np.kron(jmap, np.eye(2))@bell
    joint = joint.reshape(2, source.size, 2)
    rho_qr = np.einsum("akr,bks->arbs", joint, joint.conj()).reshape(4, 4)
    reference_error = float(np.max(np.abs(rho_qr-np.outer(bell, bell.conj()))))
    assert reference_error < 2e-14
    # Retaining L gives a pure coherent vector, discarding L gives rank three here.
    source_columns = source.reshape(-1, 3)
    eigenvalues = np.linalg.eigvalsh(source_columns.conj().T@source_columns)
    assert np.max(np.abs(eigenvalues-1/3)) < 2e-14
    return {"radial_normalization_error": radial_norm_error,
            "source_gram_error": gram_error, "source_isometry_error": iso_error,
            "total_momentum_intertwining_error": momentum_error,
            "reference_kinetic_identity_error": kinetic_error,
            "kernel_parity_error": parity_error,
            "passive_reference_reduction_error": reference_error,
            "source_reduced_nonzero_eigenvalues": eigenvalues.tolist(),
            "quadrature_scope": "angular 2-design for source Gram/parity only; no optical dynamics discretization"}


def run():
    return {
        "round": 1055, "date": "2026-10-09", "science_checks_passed": True,
        "counting": {"new_project_conditional_physical_dictionary_groups": 1,
                     "new_cognitive_axioms": 0, "global_count_assigned_by_mainline": True,
                     "whole_roadmap_completed": False},
        "scope": {"adopted_3D_finite_band_RWA_family": True,
                  "finite_mass_reference_and_full_output_retained": True,
                  "fixed_motion_preparation_arbitrary_qubit_and_passive_reference": True,
                  "three_alternative_preparations_not_cloning": True,
                  "linear_drive_not_globally_bounded_below": True,
                  "autonomous_controller_or_complete_source_energy_closed": False,
                  "real_SM_material_matching": False,
                  "full_QED_or_counter_rotating_error_certified": False,
                  "actual_endpoint_group_or_spatial_dimension_generated": False,
                  "calibration_not_continuous_motion_simulation": True},
        "historical_sha256": {rel: digest(ROOT/rel) for rel in HISTORY},
        "rational_certificate": rational_certificate(),
        "forced_propagator": force_certificate(),
        "gaussian_and_recoil": gaussian_calibration(),
        "source_reference_parity": source_and_parity_calibration(),
    }


def compare(a, b, path="root"):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys(), path
        for k in a:
            compare(a[k], b[k], path+"."+k)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+f"[{i}]")
    elif isinstance(a, float):
        assert np.isfinite(a) and np.isfinite(b), path
        assert abs(a-b) <= 3e-13*max(1., abs(a), abs(b)), (path, a, b)
    else:
        assert a == b, (path, a, b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with RESULT.open("x", encoding="utf-8") as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write("\n")
    else:
        compare(result, json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"round": 1055, "passed": True,
                      "mode": "exclusive_first_write" if args.write else "read_only",
                      "history": len(HISTORY),
                      "jacobian_error_upper": result["rational_certificate"]["exact"]["total_jacobian_error_upper"],
                      "new_co_lipschitz_lower": result["rational_certificate"]["exact"]["new_co_lipschitz_lower"]}))


if __name__ == "__main__":
    main()
