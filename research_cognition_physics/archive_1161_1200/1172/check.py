"""1172: the original complete material's worldline source and weak Einstein law.

This adds no physical operation or point-source detector. Source matching is
an explicit model dictionary. All weak tests are smooth and compactly
supported inside an already certified free leg. Default is read-only
comparison; --write-once refuses to replace results.json.
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
BASE = HERE.parents[1]
ETA = np.diag([1., -1., -1., -1.])
CAPTURE = math.pi/4
SOURCES = {
    'archive_935_955/research_note_941.md': 'a12e2dd67212a77538f850908efc6057e450d4def00dccd59aaa7cfaea346b7f',
    'archive_935_955/research_note_953.md': 'e510568c8812a86a6829be66d6733195858fd7e72c6686e14f060afce3bc2d88',
    'archive_1106_1119/1112/check.py': '6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f',
    'archive_1106_1119/1119/proof.md': 'b4a97154a8302a49c4f4f6788f615909b4849820bfa42e3bf11574c16808008d',
    'archive_1106_1119/1115/proof.md': 'b7b2239946e44cb1aade62d1282525d998297deb6dcadc36b3005dc7c454a8da',
    'archive_1106_1119/_shared/1116_1119/adoption_and_dependencies.md': '1c1395c16cf2392be0c6e7275d660f454113290bd0e23d6892e08b3aa41b4a0e',
    'archive_1120_1158/_shared/1157_1158/adoption_and_dependencies.md': 'f7af93c680768e4e786edf9aad5e1f95c8c6759e69f2ac94137219960bfe95af',
    'archive_1161_/1166/proof.md': '02335c08c46865b438b075a0d8475e143fc113e01eb711f725b96e4290bb0812',
    'archive_1161_/1166/check.py': '932c295d23a7f2031d41f58ffed75d4d0cc799c6fde6e2148b5386816ebf7099',
    'archive_1161_/1166/results.json': 'd1bf0dc11287685e6eb54b0b4c6a9af36c11b7bca82a51af86d8fd57ec3e7cf3',
    'archive_1161_/1166/acceptance.json': '7a660aa0b7d74fee7320addb51f0f124b9979b20dc64825f22279a64b93bc51a',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def norm(a):
    return float(np.linalg.norm(a, 2))


def zero(a, label, tol=4e-10):
    value = norm(a)
    assert value < tol, (label, value)
    return value


def scalar_zero(x, label, tol=4e-10):
    value = float(abs(x))
    assert value < tol, (label, value)
    return value


def source_bindings():
    for name, expected in SOURCES.items():
        assert digest(BASE/name) == expected, name
    receipt = BASE/'archive_1161_/1166/acceptance.json'
    accepted = read_json(receipt)
    assert accepted['status'] == 'ACCEPTED_WITH_SCOPE'
    files = accepted['files_sha256']
    for name, expected in files.items():
        assert digest(receipt.parent/name) == expected, ('1166 accepted core', name)
    return len(files)


def material_group():
    path = BASE/'archive_1106_1119/1112/check.py'
    spec = importlib.util.spec_from_file_location('r1112_for_1172', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    material = module.Material(1)
    mass = material.baseline*np.eye(288)+material.hfk
    eigenvalues, eigenvectors = np.linalg.eigh(mass)
    baseline = 17+2*math.sqrt(19)
    spectrum_error = zero(
        eigenvalues.reshape(-1, 1)-np.repeat(baseline+np.arange(4), 72)[:, None],
        'all 288 eigenvalues, not prepared-code compression')
    assert eigenvalues[0] > 0
    u = np.array([1.25, .75, 0., 0.])
    p = [x*mass for x in u]
    momentum_error = zero(p[0]@p[0]-sum(x@x for x in p[1:])-mass@mass,
                          'same full mass shell', tol=5e-9)
    duration = .5
    displacement = duration*u
    phase_error = zero(sum(ETA[i, i]*displacement[i]*p[i] for i in range(4))
                       -duration*mass, 'same actual free-leg phase')
    joint = np.eye(3, dtype=complex)[:, None, :]
    prepared = material.encode_joint(joint, CAPTURE).reshape(288, 3)
    U = (eigenvectors*np.exp(-1j*duration*eigenvalues)) @ eigenvectors.conj().T
    output = U@prepared
    mask = np.tile(np.repeat(np.arange(8) % 2, 9), 4)
    old_effect = output.conj().T@(mask[:, None]*output)
    probability = math.sin((CAPTURE+duration)/2)**6
    clock_error = zero(old_effect-probability*np.eye(3), 'original complete Mr')
    gram_error = zero(output.conj().T@output-np.eye(3), 'unknown input/R isometry')
    compressed_mass = prepared.conj().T@mass@prepared
    compression_error = zero(compressed_mass-(baseline+1.5)*np.eye(3),
                             'actual prepared source, no new mass-eigenstate permission')
    # A mathematical lapse variation with integrated proper-time change eps:
    # i U(0)^dagger d_eps U(eps) = M. No physical background switch is granted.
    eps = 2e-6
    plus = (eigenvectors*np.exp(-1j*(duration+eps)*eigenvalues)) @ eigenvectors.conj().T
    minus = (eigenvectors*np.exp(-1j*(duration-eps)*eigenvalues)) @ eigenvectors.conj().T
    source_numeric = 1j*U.conj().T@(plus-minus)/(2*eps)
    variation_error = norm(source_numeric-mass)
    truncation_bound = float(eigenvalues[-1]**3*eps**2/6)
    assert variation_error <= truncation_bound+2e-8
    return dict(
        full_material_dimension=288, open_logical_dimension=3,
        baseline=baseline, mass_min=float(eigenvalues[0]), mass_max=float(eigenvalues[-1]),
        spectrum_levels=(baseline+np.arange(4)).tolist(), spectral_multiplicity=72,
        complete_spectrum_error=spectrum_error, complete_mass_shell_error=momentum_error,
        original_phase_error=phase_error, capture_age=CAPTURE, free_duration=duration,
        original_Mr_probability=probability, original_full_effect_error=clock_error,
        all_input_output_gram_error=gram_error,
        prepared_source_mean=baseline+1.5, prepared_source_compression_error=compression_error,
        full_operator_source_variation_error=variation_error,
        central_difference_truncation_bound=truncation_bound,
        source_variation_step=eps, new_physical_operations=0,
    ), mass


def bump(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    inside = abs(y) < 1
    out[inside] = np.exp(1-1/(1-y[inside]**2))
    return out


def bump_derivative(y):
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    inside = abs(y) < 1
    out[inside] = -2*y[inside]*bump(y[inside])/(1-y[inside]**2)**2
    return out


def weak_group(mass):
    x, weights = np.polynomial.legendre.leggauss(128)
    xf, wf = np.polynomial.legendre.leggauss(256)
    j = float(weights@bump(x))
    jf = float(wf@bump(xf))
    quadrature_difference = abs(j-jf)
    assert quadrature_difference < 2e-10
    rho, separation, width, center = .05, .2, .2, CAPTURE+.5
    spatial_left = (rho*j)**3
    # Same physical shape, evaluated in independently translated coordinates.
    shifted_coordinates = separation+rho*x
    shifted_integral = rho*float(weights@bump((shifted_coordinates-separation)/rho))
    spatial_right = shifted_integral*(rho*j)**2
    volume_numeric = spatial_left-spatial_right
    assert abs(volume_numeric) < 1e-16
    line_numeric = jf/j
    assert abs(line_numeric-1) < 2e-10
    assert float(bump(np.array([0.]))[0]) == 1.
    assert float(bump(np.array([-separation/rho]))[0]) == 0.
    # For an asymmetric compact test along the source line, the integral of
    # its total derivative equals its zero endpoint values (analytic identity).
    poly = 1+.3*x+.1*x*x
    derivative_integral = float(weights@(
        bump_derivative(x)*poly+bump(x)*(.3+.2*x)))
    assert abs(derivative_integral) < 2e-9
    boost = np.array([[1.25, .75, 0, 0], [.75, 1.25, 0, 0],
                      [0, 0, 1, 0], [0, 0, 0, 1]])
    lorentz_error = zero(boost.T@ETA@boost-ETA, 'common proper boost')
    u = boost@np.array([1., 0, 0, 0])
    contravariant_source_coefficient = np.outer(u, u)*line_numeric
    lowered_u = ETA@u
    invariant_source = float(lowered_u@contravariant_source_coefficient@lowered_u)
    contraction_error = scalar_zero(invariant_source-line_numeric, 'observer scalar transport')
    # A t'=const cut includes the 1/u'^0 Jacobian from the worldline delta.
    cut_flux = np.outer(u, u)[0, :]/u[0]
    flux_error = scalar_zero(np.linalg.norm(cut_flux-u), 'same P=Mu surface flux')
    coordinates = np.column_stack([center+width*xf, np.zeros((256, 3))])
    transported = coordinates@boost.T
    recovered = transported@np.linalg.inv(boost).T
    event_transport_error = float(np.max(abs(recovered-coordinates)))
    recovered_line = float(wf@bump((recovered[:, 0]-center)/width)/j)
    scalar_zero(recovered_line-line_numeric, 'same compact scalar test after boost', tol=2e-10)
    kappa = 8*math.pi*.001
    eye = np.eye(len(mass))
    lambda_cases = []
    for lam in (-7., 0., 1.e6):
        residual = lam*volume_numeric*eye-kappa*line_numeric*mass
        error = norm(residual+kappa*mass)
        assert error < 2e-9
        lambda_cases.append(dict(Lambda=lam, residual_operator_norm=norm(residual),
                                 error_from_minus_kappa_M=error))
    return dict(
        test=dict(time_center=center, time_half_width=width,
                  spatial_half_width=rho, translated_bump_center=separation,
                  strict_free_open_time_interval=[CAPTURE+.25, CAPTURE+.75],
                  closed_time_support=[center-width, center+width]),
        bump_integral_128=j, bump_integral_256=jf,
        bump_quadrature_difference=quadrature_difference,
        spatial_positive_lobe_integral=spatial_left,
        spatial_negative_lobe_integral=spatial_right,
        four_volume_integral_numeric=volume_numeric,
        worldline_integral_numeric=line_numeric,
        analytic_four_volume_integral=0, analytic_worldline_integral=1,
        analytic_equalities_proved_by_translation_and_normalization=True,
        asymmetric_test_divergence_integral_numeric=derivative_integral,
        common_lorentz_error=lorentz_error,
        transformed_worldline_four_velocity=u.tolist(),
        transformed_source_coefficient=contravariant_source_coefficient.tolist(),
        transported_observer_contraction_error=contraction_error,
        physical_surface_flux_error=flux_error,
        inverse_event_transport_error=event_transport_error,
        lambda_cases=lambda_cases,
        contact_domain_stress_claimed=False,
    )


def finite_scope_group(mass):
    values = np.linalg.eigvalsh(mass)
    minimum, maximum = float(values[0]), float(values[-1])
    # Conditional error responsibility, not a manufactured detector budget.
    G, source_error, geometry_error = .001, .01, .02
    kappa = 8*math.pi*G
    lower = kappa*(minimum-source_error)-geometry_error
    assert lower > 0
    expectation_floor = kappa*minimum
    # One finite collection of normalized weak tests, with fixed coefficients.
    coefficients = [1., .4, 1.6]
    B = max(coefficients)*maximum
    tolerance = 1e-3
    rows = []
    for power in (1, 4, 12):
        candidate_G = tolerance/(8*math.pi*B)*2.**(-power)
        residuals = [8*math.pi*candidate_G*c*maximum for c in coefficients]
        assert 0 < max(residuals) < tolerance
        rows.append(dict(power=power, positive_G=candidate_G,
                         largest_weak_residual=max(residuals),
                         below_fixed_tolerance=True))
    return dict(
        fixed_G=G, fixed_kappa=kappa,
        exact_any_prepared_state_expectation_floor=expectation_floor,
        exact_operator_residual_norm=kappa*maximum,
        hypothetical_source_error=source_error, hypothetical_geometry_error=geometry_error,
        conditional_finite_error_lower_bound=lower,
        actual_source_or_curvature_measurement_precision_certified=False,
        finite_test_coefficients=coefficients, finite_test_source_norm_max=B,
        fixed_residual_tolerance=tolerance,
        arbitrarily_small_G_comparison=rows,
        small_equation_residual_proves_neighboring_Einstein_solution=False,
        G_may_not_be_refitted_per_precision_in_a_fixed_physical_theory=True,
    )


def calculate():
    bound_count = source_bindings()
    material, mass = material_group()
    weak = weak_group(mass)
    finite = finite_scope_group(mass)
    assert source_bindings() == bound_count
    return dict(
        schema='round1172_complete_worldline_weak_einstein_nonentailment_v1',
        status='PASS_WITH_SCOPE', check_groups=3,
        code_sha256=digest(Path(__file__)), source_sha256=SOURCES,
        old1166_accepted_bindings_verified=bound_count,
        group1_original_complete_material=material,
        group2_compact_weak_source_and_covariance=weak,
        group3_fixed_coupling_and_finite_precision=finite,
        scope=dict(
            premise='S_eff plus explicit free-leg worldline monopole source dictionary',
            conclusion='the old common model fails the local weak Einstein necessary equation for every fixed G>0 and constant Lambda',
            complete_contact_stress_not_constructed=True,
            no_new_actual_source_or_metric_instruments=True,
            no_uniform_fixed_precision_exclusion_of_all_positive_G=True,
            new_adopted_axioms=0, new_empirical_calibration=0,
            GR_positive_hard_deliverables_completed=0,
            broader_constructive_goal_completed=False))


def compare(actual, saved):
    if isinstance(saved, dict):
        assert actual.keys() == saved.keys()
        for key in saved:
            compare(actual[key], saved[key])
    elif isinstance(saved, list):
        assert len(actual) == len(saved)
        for a, b in zip(actual, saved):
            compare(a, b)
    elif isinstance(saved, float):
        assert math.isclose(actual, saved, rel_tol=2e-9, abs_tol=3e-10), (actual, saved)
    else:
        assert actual == saved, (actual, saved)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write_once:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        compare(result, read_json(path))
    print(json.dumps(dict(
        status=result['status'], mode='WRITE_ONCE' if args.write_once else 'READ_ONLY_MATCH',
        groups=3, old_accepted_bindings=result['old1166_accepted_bindings_verified'],
        full_dimension=288,
        source_minimum=result['group1_original_complete_material']['mass_min'],
        conditional_bound=result['group3_fixed_coupling_and_finite_precision']['conditional_finite_error_lower_bound'],
        actual_geometry_instruments_certified=False), ensure_ascii=False))


if __name__ == '__main__':
    main()
