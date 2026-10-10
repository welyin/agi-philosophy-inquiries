"""1054: same finite-band dipole action, actual probability coordinates.

Default is a read-only replay. --write creates only the first result, exclusively.
Fraction bounds certify the continuous theorem; quadrature/evolution are diagnostics.
No imported author science functions and no installation required.
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
    "archive_370_428/research_note_427.md",
    "archive_467_530/research_note_488.md",
    "archive_467_530/research_note_522.md",
    "archive_467_530/research_note_523.md",
    "archive_1009_1043/research_note_1038.md",
    "archive_1046_/research_note_1046.md",
    "archive_1046_/1046/proof.md",
    "archive_1046_/1046/photon_direction_receiver_results.json",
    "archive_1046_/1046/research_round_1046_checks.json",
    "archive_1046_/research_note_1051.md",
    "archive_1046_/_admission/finite_direction_descent/selection.md",
    "archive_1046_/_admission/optical_coordinates/selection.md",
)
G = .1
OMEGA = 1.5
TIME = .01
RADIUS = .25
I2 = np.eye(2, dtype=complex)
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pauli(v):
    return np.einsum("i,ijk->jk", v, SIGMA)


def rational_certificate():
    t, r = F(1, 100), F(1, 4)
    z = F(30, 7) * t
    leading = F(1, 100) * t*t * F(161, 1728)
    derivative_remainder = F(2, 15) * z**4 / (6*(1-z))
    kappa = leading - derivative_remainder
    probability_remainder = z**4 / (24*(1-z))
    probability_lower = F(1, 100)*t*t*F(49, 96)-probability_remainder
    c_hessian = F(16, 49)*t*t+F(8, 343)*t**3
    local_radius = F(1, 100000)
    jacobian_error = 6*c_hessian*local_radius
    co_lipschitz = F(1, 200000000)
    assert leading == F(161, 1728000000)
    assert derivative_remainder == F(9, 114905000)
    assert kappa == F(589541, 39711168000000)
    assert kappa > F(1, 100000000)
    assert probability_lower > F(1, 3000000)
    assert jacobian_error < F(1, 500000000)
    assert kappa-jacobian_error > co_lipschitz
    values = {
        "r0": r, "t0": t, "z": z, "minus_F_prime_lower": F(161, 1728),
        "leading_slope_lower": leading, "derivative_remainder_upper": derivative_remainder,
        "exact_slope_lower": kappa, "probability_remainder_upper": probability_remainder,
        "exact_probability_lower": probability_lower,
        "hessian_component_bound": c_hessian,
        "local_ball_radius": local_radius,
        "jacobian_perturbation_upper": jacobian_error,
        "co_lipschitz_lower": co_lipschitz,
        "inverse_lipschitz_upper": 1/co_lipschitz,
        "position_instrument_half_diamond_per_length_upper": 2*t/7,
        "direction_t1_gap_lower_at_origin": F(6889, 4000000),
        "direction_t1_gap_lower_in_ball": F(6889, 4000000)-F(8, 7)*local_radius,
    }
    assert values["direction_t1_gap_lower_in_ball"] > 0
    return {"exact": {k: str(v) for k, v in values.items()},
            "decimal": {k: float(v) for k, v in values.items()}}


def radial(n):
    x, w = np.polynomial.legendre.leggauss(n)
    k = 1.5+x/2
    weights = w/2*(8/3)*np.sin(np.pi*(k-1))**4
    raw = float(weights.sum())
    weights /= raw
    return k, weights, raw


def angular_axis(k, r, count=24):
    """Direct mu integral of the physical transverse matrix and its r derivative."""
    mu, weight = np.polynomial.legendre.leggauss(count)
    phase = np.outer(k*r, mu)
    factors = np.array([(1+mu*mu)/4, (1-mu*mu)/2])
    b = np.cos(phase) @ (weight*factors).T
    db = (-k[:, None]*mu[None, :]*np.sin(phase)) @ (weight*factors).T
    return b, db


def angular_closed(k, r):
    """Independent spherical-Bessel identities, used only away from zero."""
    s = k*r
    assert np.min(abs(s)) > .1
    j0 = np.sin(s)/s
    j1_over_s = (np.sin(s)-s*np.cos(s))/s**3
    dj0 = (s*np.cos(s)-np.sin(s))/s**2
    dj1_over_s = (s*s*np.sin(s)-3*np.sin(s)+3*s*np.cos(s))/s**4
    b = np.column_stack([j0-j1_over_s, 2*j1_over_s])
    db = k[:, None]*np.column_stack([dj0-dj1_over_s, 2*dj1_over_s])
    return b, db


def angular_matrix(k, vector, nc=18, nphi=24):
    mu, w = np.polynomial.legendre.leggauss(nc)
    ans = np.zeros((3, 3), complex)
    gram = np.zeros((3, 3))
    aa = np.zeros((2, 2), complex)
    for c, wc in zip(mu, w):
        for phi in 2*np.pi*np.arange(nphi)/nphi:
            n = np.array([np.sqrt(1-c*c)*np.cos(phi), np.sqrt(1-c*c)*np.sin(phi), c])
            p = np.eye(3)-np.outer(n, n)
            weight = wc/(2*nphi)
            ans += weight*np.exp(1j*k*np.dot(n, vector))*p
            gram += weight*p
            for i in range(3):
                for j in range(3):
                    aa += weight*p[i, j]*(SIGMA[i]@SIGMA[j])
    return ans, gram, aa


def star(n):
    k, w, norm = radial(n)
    h = np.diag(np.r_[OMEGA, k]).astype(complex)
    h[0, 1:] = h[1:, 0] = G*np.sqrt(2*w)
    e, v = np.linalg.eigh(h)
    u = (v*np.exp(-1j*TIME*e))@v.conj().T
    # Independent direct Taylor evolution, no spectral diagonalization.
    term = np.eye(len(h), dtype=complex)
    series = term.copy()
    for order in range(1, 19):
        term = (-1j*TIME/order)*(h@term)
        series += term
    assert np.linalg.norm(u-series, 2) < 3e-14
    return k, w, h, u, norm, float(np.linalg.norm(u-series, 2))


def response(r, data):
    k, w, _, u, _, _ = data
    b, db = angular_axis(k, r)
    coeff = np.column_stack([b[:, 0], b[:, 0], b[:, 1]])
    derivative = np.column_stack([db[:, 0], db[:, 0], db[:, 1]])
    init = np.vstack([np.zeros((1, 3)), np.sqrt(w[:, None]/2)*coeff])
    init_prime = np.vstack([np.zeros((1, 3)), np.sqrt(w[:, None]/2)*derivative])
    out = u@init
    out_prime = u@init_prime
    amp = out[0]
    p = float(np.vdot(amp, amp).real/2)
    dp = float(np.vdot(amp, out_prime[0]).real)
    return p, dp, out, coeff


def full_instrument_diagnostic(data):
    """Controlled angular isometry for this fixed mixed source, not all coherent sources.

    Keep the three source-purification labels coherently by stacking the three maps.
    Actual continuous-field Kraus maps and their common domain are proved in proof.md.
    """
    k, w, h, u, _, _ = data
    p, dp, out, b = response(RADIUS, data)
    parts = []
    effect = np.zeros((2, 2), complex)
    for j in range(3):
        dark = np.sqrt(w*np.maximum(2/3-b[:, j]**2/2, 0))*np.exp(-1j*k*TIME)
        # Each j block has excited spin, coupled bright spin, and unused dark spin.
        vj = np.vstack([out[row, j]*SIGMA[j] for row in range(len(k)+1)]
                       + [value*I2 for value in dark])/np.sqrt(2)
        ke = np.zeros_like(vj)
        ke[:2] = vj[:2]
        parts.append((vj, ke, vj-ke))
        effect += ke.conj().T@ke
    v = np.vstack([item[0] for item in parts])
    ke = np.vstack([item[1] for item in parts])
    kg = np.vstack([item[2] for item in parts])
    completeness = float(np.linalg.norm(v.conj().T@v-I2, 2))
    scalar = float(np.linalg.norm(effect-p*I2, 2))
    assert completeness < 3e-13 and scalar < 3e-18
    rng = np.random.default_rng(1054)
    seed = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    rho = seed@seed.conj().T
    rho /= np.trace(rho)
    rho4 = rho.reshape(2, 2, 2, 2)
    reference_before = np.einsum("aras->rs", rho4)
    reference_after = np.zeros((2, 2), complex)
    actual_probability = None
    for index, kb in enumerate([ke, kg]):
        # Trace the actual physical output while preserving the old reference.
        reference_b = np.einsum("oa,arbs,ob->rs", kb, rho4, kb.conj())
        reference_after += reference_b
        if index == 0:
            actual_probability = float(np.trace(reference_b).real)
    ref_error = float(np.linalg.norm(reference_after-reference_before, 2))
    assert ref_error < 3e-13 and abs(actual_probability-p) < 3e-18
    # Source labels are not traced before constructing the coherent purification map.
    # Report a genuine excited-state effect on |+z>, after tracing only that label.
    test = np.diag([1., 0.]).astype(complex)
    exc = sum(item[1][:2]@test@item[1][:2].conj().T for item in parts)
    normalized = exc/p
    luders_error = float(np.sum(abs(np.linalg.eigvalsh(normalized-test)))/2)
    assert luders_error > .6
    comm = np.linalg.norm(h@np.diag(np.r_[1., np.zeros(len(k))])-
                          np.diag(np.r_[1., np.zeros(len(k))])@h, 2)
    assert abs(comm-G*np.sqrt(2)) < 2e-14
    assert np.linalg.eigvalsh(h)[0] > 6/7 and np.linalg.eigvalsh(h)[-1] < 15/7
    return {
        "preparation_gram_diagonal": [1/3]*3,
        "coherent_source_purification_retained": True,
        "isometry_completeness_residual": completeness,
        "scalar_effect_residual": scalar,
        "unknown_reference_probability": actual_probability,
        "reference_marginal_residual": ref_error,
        "conditional_excited_spin_z": float(np.trace(normalized@SIGMA[2]).real),
        "conditional_excited_spin_distance_from_luders": luders_error,
        "terminal_read_commutator_norm": float(comm),
        "finite_quadrature_H_extreme_eigenvalues": [float(np.linalg.eigvalsh(h)[0]), float(np.linalg.eigvalsh(h)[-1])],
        "radial_probability_derivative": dp,
    }


def run():
    cert = rational_certificate()
    data = star(24)
    finer = star(48)
    k, w, _, _, norm, exp_error = data
    b, db = angular_axis(k, RADIUS)
    bc, dbc = angular_closed(k, RADIUS)
    angular_error = float(max(np.max(abs(b-bc)), np.max(abs(db-dbc))))
    assert angular_error < 2e-12
    matrix_errors = []
    gram_error = self_error = 0.
    for vec in [np.array([RADIUS, 0., 0.]), np.array([0., 0., RADIUS]), RADIUS*np.array([1., 2., 3.])/np.sqrt(14)]:
        for kk in [1.1, 1.5, 1.9]:
            mat, gram, aa = angular_matrix(kk, vec)
            bs, _ = angular_closed(np.array([kk]), RADIUS)
            n = vec/RADIUS
            expected = bs[0, 0]*np.eye(3)+(bs[0, 1]-bs[0, 0])*np.outer(n, n)
            matrix_errors.append(float(np.linalg.norm(mat-expected, 2)))
            gram_error = max(gram_error, float(np.linalg.norm(gram-2*np.eye(3)/3, 2)))
            self_error = max(self_error, float(np.linalg.norm(aa-2*I2, 2)))
    assert max(matrix_errors) < 2e-12 and gram_error < 2e-12 and self_error < 2e-12
    p, dp, _, _ = response(RADIUS, data)
    pf, dpf, _, _ = response(RADIUS, finer)
    assert p > cert["decimal"]["exact_probability_lower"]
    assert -dp > cert["decimal"]["exact_slope_lower"]
    f0 = float((2*np.dot(w, b[:, 0])**2+np.dot(w, b[:, 1])**2)/2)
    fp = float(2*np.dot(w, b[:, 0])*np.dot(w, db[:, 0])+np.dot(w, b[:, 1])*np.dot(w, db[:, 1]))
    assert -fp >= float(F(161, 1728))
    assert abs(p-G*G*TIME*TIME*f0) < cert["decimal"]["probability_remainder_upper"]
    assert abs(dp-G*G*TIME*TIME*fp) < cert["decimal"]["derivative_remainder_upper"]
    anchors = RADIUS*np.eye(3)
    origin_prob = np.full(3, p)
    analytic_jacobian = -dp*np.eye(3)
    hstep = 2e-5
    fd_jacobian = np.empty((3, 3))
    def probabilities(x):
        return np.array([response(np.linalg.norm(x-a), data)[0] for a in anchors])
    for j in range(3):
        shift = hstep*np.eye(3)[j]
        fd_jacobian[:, j] = (probabilities(shift)-probabilities(-shift))/(2*hstep)
    jac_error = float(np.linalg.norm(fd_jacobian-analytic_jacobian, 2))
    assert jac_error < 5e-14
    checks = []
    for x in [np.array([1., 2., -1.])*1e-6, np.array([-3., 0., 2.])*1e-6]:
        prob = probabilities(x)
        ratio = np.linalg.norm(prob-origin_prob)/np.linalg.norm(x)
        assert ratio > cert["decimal"]["co_lipschitz_lower"]
        checks.append({"x": x.tolist(), "probabilities": prob.tolist(), "response_per_length": float(ratio)})
    # Countercontrols: the same instrument with one center loses all angular data.
    same_radius = [np.array([RADIUS, 0., 0.]), np.array([0., RADIUS, 0.])]
    same_read = [response(np.linalg.norm(x), data)[0] for x in same_radius]
    assert same_read[0] == same_read[1]
    collinear = np.vstack([-response(abs(a), data)[1]*np.sign(a)*np.array([1., 0., 0.])
                          for a in [-RADIUS, RADIUS, 2*RADIUS]])
    assert np.linalg.matrix_rank(collinear, tol=1e-15) == 1
    instrument = full_instrument_diagnostic(data)
    return {
        "round": 1054, "date": "2026-10-08", "author_scientific_checks_passed": True,
        "new_project_conditional_physical_dictionary_groups": 1,
        "new_cognitive_axioms": 0, "full_roadmap_completed_here": False,
        "parameters": {"g": G, "Omega": OMEGA, "time": TIME, "anchor_radius": RADIUS,
                       "adopted_dimension": 3, "radial_quadratures": [24, 48],
                       "fixed_local_ball_radius": 1e-5},
        "rational_certificate": cert,
        "probability_and_slope": {"p": p, "minus_p_prime": -dp, "refined_p": pf,
                                  "refined_minus_p_prime": -dpf,
                                  "leading_F": f0, "minus_F_prime": -fp,
                                  "full_minus_leading_probability": p-G*G*TIME*TIME*f0,
                                  "full_minus_leading_derivative": dp-G*G*TIME*TIME*fp},
        "independent_diagnostics": {"radial_raw_normalization": norm,
            "same_free_energy_mean": float(np.dot(w, k)),
            "angular_integral_vs_bessel_max_residual": angular_error,
            "full_transverse_matrix_max_residual": max(matrix_errors),
            "source_gram_residual": gram_error, "self_energy_residual": self_error,
            "star_spectral_vs_direct_taylor_residual": exp_error,
            "quadrature_probability_difference": abs(p-pf),
            "quadrature_derivative_difference": abs(dp-dpf),
            "analytic_jacobian": analytic_jacobian.tolist(),
            "finite_difference_jacobian_residual": jac_error,
            "neighborhood_checks": checks},
        "complete_instrument_diagnostic": instrument,
        "deletion_controls": {"one_center_equal_radius_equal_probability": same_read,
                               "collinear_at_origin_jacobian_rank": 1},
        "scope": {
            "three_dimensional_Maxwell_and_fixed_anchor_frame_adopted": True,
            "c_number_receiver_position_parameter_not_quantum_position": True,
            "same_dipole_effective_action_translation_family_not_autonomous_receiver_motion": True,
            "three_alternative_experiments_not_cloned_input_or_free_reset": True,
            "normal_photon_preparation_permissions_adopted": True,
            "short_time_RWA_error_relative_to_full_QED_certified": False,
            "all_unknown_spin_and_passive_reference_exact_in_proof": True,
            "retained_source_purification_uses_coherent_sum": True,
            "terminal_detector_energy_closed": False,
            "actual_displacement_or_half_cost_contract_derived": False,
            "three_dimensions_generated": False,
            "M3A_complete": False,
            "quadrature_differences_are_not_continuous_error_certificates": True,
        },
        "historical_sha256": {rel: digest(ROOT/rel) for rel in HISTORY},
    }


def compare(a, b, path=""):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for key in a:
            compare(a[key], b[key], path+"/"+key)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+"/"+str(i))
    elif isinstance(a, float):
        assert np.isclose(a, b, rtol=5e-10, atol=3e-13), (path, a, b)
    else:
        assert a == b, (path, a, b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    results = run()
    if args.write:
        with RESULT.open("x", encoding="utf-8") as output:
            json.dump(results, output, ensure_ascii=False, indent=2)
            output.write("\n")
    else:
        compare(results, json.loads(RESULT.read_text(encoding="utf-8")))
    print(json.dumps({"round": 1054, "passed": True,
        "mode": "exclusive_first_write" if args.write else "read_only",
        "p": results["probability_and_slope"]["p"],
        "minus_p_prime": results["probability_and_slope"]["minus_p_prime"],
        "exact_slope_lower": results["rational_certificate"]["exact"]["exact_slope_lower"]}))


if __name__ == "__main__":
    main()
