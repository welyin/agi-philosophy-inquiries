"""1030: finite-window two-channel scattering selection, not a QFT completion.

The analytic theorem concerns an orthonormal unit-flux channel compression of a
unitary S.  HEFT supplies a separately declared leading matrix.  No actual SM
channel map or numerical EFT remainder certificate is constructed here.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "finite_window_scattering_selection_results.json"
TOL = 2e-12


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def integral_polynomial(coefficients):
    return sum((F(2)*c/F(i+1) for i, c in enumerate(coefficients)
                if i % 2 == 0), F(0))


def projection_certificate():
    # Divide elastic amplitudes by x*s/v^2 and inelastic by y*s/v^2.
    # t_M/s=-(1-z)/2, u_M/s=-(1+z)/2; t_M is Mandelstam t.
    a_s = [F(1), F(0)]
    a_t = [F(-1, 2), F(1, 2)]
    a_u = [F(-1, 2), F(-1, 2)]
    singlet = [3*a_s[i]+a_t[i]+a_u[i] for i in range(2)]
    elastic_integral = integral_polynomial(singlet)
    # Singlet ww->hh includes sqrt(3), recorded separately exactly.
    transition_integral_over_sqrt3 = integral_polynomial([F(1)])
    # Projection 1/(64*pi); t_energy=s/(16*pi*v^2).
    elastic_factor = elastic_integral/4
    transition_factor_over_sqrt3 = transition_integral_over_sqrt3/4
    assert elastic_factor == 1
    assert transition_factor_over_sqrt3 == F(1, 2)
    return dict(angular_variable="z=cos(theta)",
                energy_variable="t=s/(16*pi*v^2), not Mandelstam t",
                elastic_isosinglet_polynomial=singlet,
                exact_elastic_integral=elastic_integral,
                exact_transition_integral_over_sqrt3=transition_integral_over_sqrt3,
                K_over_t_ww_over_x=elastic_factor,
                K_over_t_wh_over_sqrt3_y=transition_factor_over_sqrt3,
                K_over_t_hh=F(0),
                scope="I=0,J=0 selected coupled block at leading two-derivative, massless/gaugeless order; not every physical SM channel")


def leading_matrix(x, y, t=1):
    x, y, t = map(float, (x, y, t))
    return t*np.array([[x, math.sqrt(3)*y/2],
                       [math.sqrt(3)*y/2, 0.0]])


def analytic_rho(x, y, t=1):
    x, y, t = map(float, (x, y, t))
    assert t >= 0
    return t*(abs(x)+math.sqrt(x*x+3*y*y))/2


def exact_distance(rho):
    # Rationalized expression avoids subtractive loss close to rho=0.
    rho = float(rho)
    return 2*rho*rho/(math.sqrt(1+4*rho*rho)+1)


def polar_completion(k):
    values, vectors = np.linalg.eigh(k)
    phases = (1+2j*values)/np.sqrt(1+4*values*values)
    s = (vectors*phases) @ vectors.conj().T
    return s, (s-np.eye(k.shape[0]))/(2j)


def region_margin_for_R(x, y, R):
    x, y, R = map(F, (x, y, R))
    assert R >= 0
    return R*R-2*R*abs(x)-3*y*y


def region_contains(x, y, R):
    x, y, R = map(F, (x, y, R))
    # The quadratic inequality loses x at R=0.  Keep the linear bound too.
    return R >= 0 and 2*abs(x) <= R and region_margin_for_R(x, y, R) >= 0


def rational_region_certificate():
    delta, t_max = F(1, 15), F(1)
    radius_squared = delta*(1+delta)
    radius = F(4, 15)
    R = 2*radius/t_max
    assert radius*radius == radius_squared
    cases = []
    named = [
        ("flat_origin", F(0), F(0)),
        ("nonzero_curvature_allowed", F(0), F(1, 5)),
        ("elastic_zero_but_joint_rejected", F(0), F(2, 5)),
        ("positive_x_boundary", F(4, 15), F(0)),
        ("negative_x_boundary", F(-4, 15), F(0)),
        ("positive_x_outside", F(1, 3), F(0)),
        ("negative_mixing_allowed", F(0), F(-1, 5)),
    ]
    for name, x, y in named:
        margin = region_margin_for_R(x, y, R)
        rho = analytic_rho(x, y, t_max)
        values = np.linalg.eigvalsh(leading_matrix(x, y, t_max))
        assert abs(max(abs(values))-rho) < TOL
        feasible = region_contains(x, y, R)
        assert feasible == (exact_distance(rho) <= float(delta)+TOL)
        cases.append(dict(name=name, x=x, y=y, a_squared=1-x,
                          b=1-x-y, exact_region_margin=margin,
                          feasible=feasible, rho=rho,
                          minimal_complex_remainder=exact_distance(rho)))
    # A rational grid cross-checks the exact conic against the analytic radical.
    grid_count = 0
    for i in range(-6, 7):
        for j in range(-6, 7):
            x, y = F(i, 15), F(j, 15)
            exact = region_contains(x, y, R)
            radical = analytic_rho(x, y) <= float(radius)+TOL
            assert exact == radical
            grid_count += 1
    assert cases[1]["feasible"] and not cases[2]["feasible"]
    assert region_contains(0, 0, 0)
    assert not region_contains(F(1, 5), 0, 0)
    return dict(delta=delta, t_max=t_max, rho_limit=radius,
                rho_limit_squared=radius_squared, R=R,
                exact_region="2*abs(x)<=R and 2*R*abs(x)+3*y^2<=R^2, intersect x<=1 for real a; R=0 requires x=y=0",
                cases=cases, exact_grid_points=grid_count,
                zero_budget_boundary=dict(origin_allowed=True,
                    nonzero_x=F(1, 5), y=F(0), nonzero_x_allowed=False,
                    quadratic_margin_alone=region_margin_for_R(F(1, 5), 0, 0)),
                interpretation="Necessary for a real certified amplitude; sufficient only for pointwise contraction-matrix existence. The budget is stipulated, not an EFT error estimate.")


def polar_and_spectrum_certificate():
    points = [(F(0), F(1, 5)), (F(0), F(2, 5)),
              (F(1, 5), F(-1, 7)), (F(-2, 7), F(1, 11)),
              (F(0), F(0))]
    rows = []
    max_unitary_error = max_distance_error = max_spectrum_error = 0.0
    for x, y in points:
        k = leading_matrix(x, y)
        rho = analytic_rho(x, y)
        s, amplitude = polar_completion(k)
        actual = float(np.linalg.norm(amplitude-k, 2))
        predicted = exact_distance(rho)
        unitary_error = float(np.linalg.norm(s.conj().T@s-np.eye(2), 2))
        spectral = np.linalg.eigvalsh(k)
        xf, yf = float(x), float(y)
        radical = math.sqrt(xf*xf+3*yf*yf)
        predicted_eigenvalues = np.array([(xf-radical)/2, (xf+radical)/2])
        spectrum_error = float(np.max(abs(spectral-predicted_eigenvalues)))
        assert abs(actual-predicted) < TOL
        assert unitary_error < TOL and spectrum_error < TOL
        max_unitary_error = max(max_unitary_error, unitary_error)
        max_distance_error = max(max_distance_error, abs(actual-predicted))
        max_spectrum_error = max(max_spectrum_error, spectrum_error)
        rows.append(dict(x=x, y=y, exact_characteristic_coefficients=
                         [F(1), -x, -F(3, 4)*y*y],
                         eigenvalues=spectral.tolist(), rho=rho,
                         actual_distance=actual, sharp_distance=predicted,
                         S_transition_probability=float(abs(s[1, 0])**2)))
    # For x=0,y=1/5: rho^2=3/100, polar conversion=4rho^2/(1+4rho^2).
    rho_squared = F(3, 100)
    conversion = 4*rho_squared/(1+4*rho_squared)
    assert conversion == F(3, 28)
    assert abs(rows[0]["S_transition_probability"]-float(conversion)) < TOL
    return dict(points=rows, exact_allowed_example_rho_squared=rho_squared,
                exact_allowed_example_polar_conversion=conversion,
                max_unitarity_residual=max_unitary_error,
                max_sharp_distance_residual=max_distance_error,
                max_spectrum_residual=max_spectrum_error,
                no_QFT_completion_claim=True)


def compressed_unitary_certificate():
    rng = np.random.default_rng(1030)
    contraction_max = -math.inf
    minimum_lower_bound_gap = math.inf
    count = 24
    for i in range(count):
        dimension = 3+i % 4
        raw = rng.normal(size=(dimension, dimension))+1j*rng.normal(size=(dimension, dimension))
        unitary, _ = np.linalg.qr(raw)
        # Orthogonal coordinate embedding realizes an isometric two-channel J.
        s = unitary[:2, :2]
        amplitude = (s-np.eye(2))/(2j)
        k = leading_matrix(F((i % 7)-3, 11), F((i % 9)-4, 13), F(3, 5))
        rho = float(np.linalg.norm(k, 2))
        contraction = float(np.linalg.norm(s, 2))
        gap = float(np.linalg.norm(amplitude-k, 2))-exact_distance(rho)
        assert contraction <= 1+TOL and gap >= -TOL
        contraction_max = max(contraction_max, contraction)
        minimum_lower_bound_gap = min(minimum_lower_bound_gap, gap)
    return dict(seed=1030, samples=count, full_dimensions=[3, 4, 5, 6],
                max_compressed_S_norm=contraction_max,
                minimum_distance_minus_proved_lower_bound=minimum_lower_bound_gap,
                role="Finite calibration of a proved operator-norm bound, not statistical evidence of universality.")


def energy_window_certificate():
    # Define a diagnostic budget delta(t) by delta(1+delta)=t^2*g(t)^2.
    # g=1/5+(t-1/2)^2 has an exact interior minimum on [1/4,1].
    def g(t):
        return F(1, 5)+(t-F(1, 2))**2
    samples = []
    rho_over_t_squared = F(3, 64)  # x=0,y=1/4
    for t in [F(1, 4), F(1, 2), F(1)]:
        ratio = g(t)
        delta = (math.sqrt(1+4*float(t*t*ratio*ratio))-1)/2
        passes = rho_over_t_squared <= ratio*ratio
        samples.append(dict(t=t, g=ratio, g_squared=ratio*ratio,
                            delta=delta, passes=passes))
    assert samples[0]["passes"] and samples[2]["passes"]
    assert not samples[1]["passes"]
    return dict(t_interval=[F(1, 4), F(1)],
                g_polynomial_ascending=[F(9, 20), F(-1), F(1)],
                budget="delta(t)=(sqrt(1+4*t^2*g(t)^2)-1)/2",
                exact_minimizer=F(1, 2), exact_minimum_rho_over_t_limit=F(1, 5),
                test_x=F(0), test_y=F(1, 4),
                test_rho_over_t_squared=rho_over_t_squared,
                samples=samples,
                constant_delta_rule="rho(t) proportional to t, so t_max is strongest only for a shared constant delta (or proved suitable monotonicity).",
                candidate_dependent_budget_rule="If delta=delta(t;a,b), retain d(rho(t;a,b))<=delta(t;a,b); no fixed conic without a common budget.",
                actual_EFT_budget_certified=False)


def curvature_certificate():
    # u=h/v, r/v=sqrt(F).  (r/v)'=a, (r/v)''=b-a^2 at u=0.
    examples = []
    for a, b in [(F(1), F(4, 5)), (F(1), F(3, 5)),
                 (F(3, 5), F(2, 7)), (F(-1), F(1))]:
        root_coefficients = [F(1), a, (b-a*a)/2]
        # Squaring the local sqrt series recovers exactly F through u^2.
        squared = [sum((root_coefficients[j]*root_coefficients[k]
                        for j in range(3) for k in range(3) if j+k == i), F(0))
                   for i in range(3)]
        assert squared == [F(1), 2*a, b]
        examples.append(dict(a=a, b=b, sqrt_F_series=root_coefficients,
                             v_squared_K_pipi=1-a*a,
                             v_squared_K_pih=a*a-b))
    # F=(1+u)^2+c*u^3, c=1: zero curvature at u=0 but not globally flat.
    c, u = F(1), F(1, 5)
    f = (1+u)**2+c*u**3
    fp = 2*(1+u)+3*c*u*u
    fpp = 2+6*c*u
    angular = (4*f-fp*fp)/(4*f*f)
    radial = (fp*fp-2*f*fpp)/(4*f*f)
    assert angular != 0 and radial != 0
    return dict(local_sqrt_checks=examples,
                flat_at_origin_not_globally_flat=dict(cubic_coefficient=c,
                   evaluation_h_over_v=u, F=f, F_prime_u=fp, F_second_u=fpp,
                   v_squared_K_pipi=angular, v_squared_K_pih=radial),
                interpretation="Scalar-target sectional curvature; not spacetime curvature or a derivation of spatial dimension. Reuses old geometric interpretation.")


def quantifier_certificate():
    delta = F(1, 15)
    # K=0, Delta=-i*delta*I lies in the error ball but S expands norms.
    s_norm = 1+2*delta
    assert s_norm > 1
    return dict(K_zero=True, wrong_remainder="-i*delta*I", delta=delta,
                exact_S_norm=s_norm, error_norm=delta,
                disproved_claim="Every matrix in an allowed error ball is unitary or contractive.",
                retained_claim="The ball intersects the contraction-amplitude set iff d(rho)<=delta; any actual certified amplitude must satisfy this.")


def run():
    sources = [
        "archive_1009_/research_note_1029.md",
        "archive_1009_/1029/NEXT.md",
        "archive_1009_/1029/input_dependency_update_v0_18.md",
        "archive_1009_/research_note_1016.md",
        "archive_1009_/research_note_1023.md",
        "archive_554_584/research_note_568.md",
        "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
        "archive_956_989/981/drafts/common_parent_contract_v1.md",
        "archive_1009_/1009/input_dependency_ledger_v0_1.md",
    ]
    result = dict(round=1030, status="scientific_calibration_verified",
        new_calibration_groups=1, cumulative_test_groups=3807,
        new_cognitive_axioms=0, goal_complete=False,
        all_scientific_calibrations_passed=True,
        scope="Finite-window conditional two-channel contraction test of a declared HEFT leading matrix; not an actual SM channel or EFT error certificate.",
        retained_inputs=["unitary full process and orthonormal unit-flux stable asymptotic channel embedding", "declared custodial HEFT two-derivative leading matrix", "an independently justified complex operator-norm remainder for the chosen finite window", "masses, widths, phase-space, equivalence-theorem, higher-order and mapping errors must be included where relevant"],
        not_proved=["physical W/Z/h exact asymptotic channel compression", "Goldstone fields are physical channels", "an actual numerical remainder for P981 or SM", "analytic crossing-symmetric QFT realization of polar matrices", "UV completion or arbitrary-energy extension", "unique a,b or all other HEFT coefficients", "space dimension from scalar-target curvature"],
        projection=projection_certificate(),
        exact_finite_window_region=rational_region_certificate(),
        polar_and_spectrum=polar_and_spectrum_certificate(),
        unitary_compression=compressed_unitary_certificate(),
        nonconstant_budget_window=energy_window_certificate(),
        scalar_target_curvature=curvature_certificate(),
        error_ball_quantifier=quantifier_certificate(),
        primary_sources=[
            dict(url="https://arxiv.org/html/1311.5993", sections="2,3,5; equations 29,34,39,41,43", role="Adopted leading coupled partial-wave normalization, not the new norm-distance proof."),
            dict(url="https://arxiv.org/html/1511.00724", sections="II; warped scalar metric and equation 47", role="Reused scalar geometry; its F is sqrt of our F.")],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        historical_source_sha256={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in sources})
    return encode(result)


def compare(actual, expected, path="root"):
    if isinstance(actual, dict):
        assert isinstance(expected, dict) and actual.keys() == expected.keys(), path
        for key in actual:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(actual, list):
        assert isinstance(expected, list) and len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(actual, float):
        assert isinstance(expected, (float, int)) and math.isclose(actual, expected, rel_tol=1e-10, abs_tol=TOL), (path, actual, expected)
    else:
        assert actual == expected, (path, actual, expected)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with OUT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(result, json.loads(OUT.read_text(encoding="utf-8")))
    print(json.dumps(dict(status="passed", round=1030,
        mode="exclusive_first_write" if args.write else "read_only_recompute_compare",
        exact_grid_points=result["exact_finite_window_region"]["exact_grid_points"],
        random_unitary_compressions=result["unitary_compression"]["samples"],
        polar_conversion=result["polar_and_spectrum"]["exact_allowed_example_polar_conversion"],
        historical_sources=len(result["historical_source_sha256"]), output=str(OUT)), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
