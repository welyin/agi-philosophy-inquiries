"""1172 independent review: original material, weak source and quantifiers.

Only the original 1112 material implementation is imported. The author's
1172 check is run as a separate, read-only subprocess. Independent evolution
uses spectral polynomials; compact integration uses y=tanh(r) and Simpson's
rule, rather than the author's diagonalization/Gauss-Legendre calculation.
Default compares review_checks.json. --write-once refuses to overwrite it.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
EXPECTED = {
    'archive_1161_/research_note_1172.md': '2230aa0a93175dc39721ea8575778e09afa9ca9c1b5d299edc1e79a3ee577638',
    'archive_1161_/1172/admission.md': '4d557fc6f657d50331a1cc4f5c0e46b35a839e209a7001c0d3c6c4c11506affe',
    'archive_1161_/1172/proof.md': '36387fcba4e591fad44ccbc7b239ace5a1601d2b96147497d2fc42e98acc4073',
    'archive_1161_/1172/check.py': '51ccd935d6961a4c200efad496fee5f16c638bb2d407e3b51c612e448b3dabe6',
    'archive_1161_/1172/results.json': 'f5c5280aad4ba08b546a8b16a0605eed42c40c6c0dd5749ae6cb64e19d143cdb',
    'archive_1161_/1172/NEXT.md': '5ae1ff4c8f43d34289001d70f7c147d3f0bdbcff75b878632add707e8d9a33f7',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def norm(a):
    return float(np.linalg.norm(a, 2))


def close(a, b, tol=3e-10):
    err = norm(np.atleast_1d(a-b))
    assert err < tol, err
    return err


def bindings(saved):
    for name, expected in EXPECTED.items():
        assert sha(BASE/name) == expected, ('core changed', name)
    sources = saved['source_sha256']
    assert len(sources) == 11
    for name, expected in sources.items():
        assert sha(BASE/name) == expected, ('source changed', name)
    receipt_path = BASE/'archive_1161_/1166/acceptance.json'
    receipt = read_json(receipt_path)
    assert receipt['status'] == 'ACCEPTED_WITH_SCOPE'
    old = receipt['files_sha256']
    assert len(old) == 8
    for name, expected in old.items():
        assert sha(receipt_path.parent/name) == expected, ('accepted old file changed', name)
    return dict(author_core_sha256=EXPECTED, source_sha256=sources,
                old1166_acceptance_sha256=sha(receipt_path),
                old1166_bound_file_sha256=old,
                all_bound_files_unchanged=True)


def original_material(saved):
    path = BASE/'archive_1106_1119/1112/check.py'
    spec = importlib.util.spec_from_file_location('original1112_reviewer1172', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    material = module.Material(1)
    h = material.hfk
    identity = np.eye(288, dtype=complex)
    baseline = 17+2*math.sqrt(19)
    close(material.baseline, baseline)
    mass = baseline*identity+h
    projectors = []
    for k in range(4):
        p = identity.copy()
        for j in range(4):
            if j != k:
                p = p@(h-j*identity)/(k-j)
        projectors.append(p)
    errors = []
    for k, p in enumerate(projectors):
        errors.extend([close(p, p.conj().T), close(p@p, p),
                       close(h@p, k*p), close(np.trace(p), 72.)])
        for j in range(k):
            errors.append(close(p@projectors[j], np.zeros_like(p)))
    errors.append(close(sum(projectors), identity))
    levels = baseline+np.arange(4)

    def evolution(t):
        return sum(np.exp(-1j*t*levels[k])*projectors[k] for k in range(4))

    delta, duration = math.pi/4, .5
    prepared = material.encode_joint(np.eye(3, dtype=complex)[:, None, :], delta).reshape(288, 3)
    out = evolution(duration)@prepared
    # Array axes are original (head=4, memory=8, payload=3, work=3).
    mr = np.indices((4, 8, 3, 3))[1].reshape(-1) % 2
    effect = out.conj().T@(mr[:, None]*out)
    expected_p = math.sin((delta+duration)/2)**6
    effect_error = close(effect, expected_p*np.eye(3))
    gram_error = close(out.conj().T@out, np.eye(3))
    source_error = close(prepared.conj().T@mass@prepared, (baseline+1.5)*np.eye(3))
    saved_error = close(expected_p, saved['group1_original_complete_material']['original_Mr_probability'])
    return dict(method='degree-three spectral projectors on original full H_FK; no 1172 import or eigendecomposition',
                full_dimension=288, actual_unknown_input_dimension=3,
                spectral_multiplicities=[float(np.trace(p).real) for p in projectors],
                maximum_projector_identity_error=max(errors),
                baseline=baseline, mass_levels=levels.tolist(),
                original_capture_age=delta, original_free_duration=duration,
                original_Mr_probability=expected_p,
                complete_Mr_effect_error=effect_error, all_input_isometry_error=gram_error,
                actual_prepared_source_compression_error=source_error,
                saved_clock_probability_difference=saved_error), mass, levels, projectors


def independent_weak_calculation(saved, mass, levels, projectors):
    # beta(tanh r) = exp(-sinh(r)^2), dy=sech(r)^2 dr.
    n, radius = 4096, 3.
    r = np.linspace(-radius, radius, n+1)
    weights = np.ones(n+1)
    weights[1:-1:2] = 4
    weights[2:-1:2] = 2
    weights *= (2*radius/n)/3
    y, jac = np.tanh(r), 1/np.cosh(r)**2
    beta = np.exp(-np.sinh(r)**2)
    J = float(weights@(beta*jac))
    J_tail_bound = float(2*math.exp(-math.sinh(radius)**2)*(1-math.tanh(radius)))
    quadrature_check = close(J, saved['group2_compact_weak_source_and_covariance']['bump_integral_256'], 2e-11)
    # This is a numerical diagnostic; exact normalization/translation are analytic.
    h, rho = .2, .05
    a = beta/(h*J)
    line = float(weights@(a*h*jac))
    close(line, 1.)
    poly = 1+.3*y+.1*y*y
    derivative = -2*y*beta/(1-y*y)**2*poly+beta*(.3+.2*y)
    divergence = float(weights@(derivative*jac))
    close(divergence, 0.)
    one_lobe_volume = (rho*J)**3

    # Independent metric, rather than time-translation, perturbation:
    # g_00(lambda,z(s))=1+lambda*a(s). The proper time is integrated first.
    def proper_time_change(lam):
        z = lam*a
        return float(weights@(h*jac*z/(np.sqrt(1+z)+1)))

    step = 1e-4
    values = [proper_time_change(mult*step) for mult in [-2, -1, 1, 2]]
    phases = [np.exp(-1j*levels*d) for d in values]
    derivative_u = (phases[0]-8*phases[1]+8*phases[2]-phases[3])/(12*step)
    source_levels = 1j*derivative_u
    full_metric_source = sum(source_levels[k]*projectors[k] for k in range(4))
    metric_variation_error = close(full_metric_source, mass/2, 2e-9)

    eta = np.diag([1., -1., -1., -1.])
    u = np.array([1.25, .25, .5, .5])
    L = np.eye(4)
    L[0, 0], L[0, 1:], L[1:, 0] = u[0], u[1:], u[1:]
    L[1:, 1:] += np.outer(u[1:], u[1:])/(u[0]+1)
    lorentz_error = close(L.T@eta@L, eta)
    det_error = close(np.linalg.det(L), 1.)
    ul = eta@u
    tensor = np.outer(u, u)
    scalar_error = close(ul@tensor@ul, 1.)
    # Tilted future spacelike cut, not merely the author's t'=constant cut.
    normal = np.array([5/3, 4/3, 0., 0.])
    close(normal@eta@normal, 1.)
    nl = eta@normal
    cut_denominator = float(nl@u)
    assert cut_denominator > 0
    flux_error = close(nl@tensor/cut_denominator, u)
    kappa = 8*math.pi*.001
    R = -kappa*mass
    residual_norm = norm(R)
    analytic_norm = kappa*levels[-1]
    norm_error = close(residual_norm, analytic_norm)
    norm_saved_error = close(residual_norm, saved['group3_fixed_coupling_and_finite_precision']['exact_operator_residual_norm'])
    return dict(method='tanh-transformed composite Simpson; full spectral-polynomial metric variation; non-collinear boost and tilted cut',
                Simpson_intervals=n, transformed_radius=radius,
                bump_integral=J, transformed_tail_bound=J_tail_bound,
                independent_quadrature_difference=quadrature_check,
                normalized_worldline_integral=line,
                equal_translated_lobe_volume=one_lobe_volume,
                four_volume_zero_basis='exact translation identity, not floating-point tolerance',
                total_derivative_integral=divergence,
                metric_variation_step=step, complete_metric_source_error_from_M_over_2=metric_variation_error,
                proper_time_perturbations=values,
                full_four_velocity=u.tolist(), non_collinear_Lorentz_error=lorentz_error,
                determinant_error=det_error, invariant_source_contraction_error=scalar_error,
                tilted_cut_normal=normal.tolist(), tilted_cut_intersection_factor=cut_denominator,
                complete_flux_coefficient_error=flux_error,
                residual_operator_norm=residual_norm, residual_norm_analytic_error=norm_error,
                residual_norm_saved_difference=norm_saved_error,
                all_prepared_state_expectation_floor=kappa*levels[0]), kappa


def finite_quantifiers(saved, levels, kappa):
    lower = kappa*(levels[0]-.01)-.02
    error = close(lower, saved['group3_fixed_coupling_and_finite_precision']['conditional_finite_error_lower_bound'])
    assert lower > 0
    B = 1.6*levels[-1]
    tol = .001
    witnesses = []
    for factor in [2., 16., 4096.]:
        G = tol/(8*math.pi*B*factor)
        r = 8*math.pi*G*B
        assert 0 < r < tol
        witnesses.append(dict(G=G, maximum_residual=r))
    return dict(conditional_source_error=.01, conditional_full_geometric_error=.02,
                conditional_lower_bound=lower, saved_bound_difference=error,
                actual_instrument_precision_certified=False,
                finite_test_source_bound=B, fixed_tolerance=tol,
                positive_G_below_fixed_tolerance=witnesses,
                zero_source_bound_case='at B=0 and Lambda=0 all selected residuals already vanish',
                no_uniform_exclusion_of_all_positive_G=True,
                small_residual_does_not_establish_Einstein_solution=True)


def text_check():
    count, fences = 0, 0
    for name in EXPECTED:
        if not name.endswith('.md'):
            continue
        path = BASE/name
        text = path.read_text(encoding='utf-8-sig')
        assert not re.search(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', text), name
        num = len(re.findall(r'^\s*\$\$\s*$', text, re.M))
        assert num % 2 == 0, name
        fences += num
        for target in re.findall(r'\[[^\]]*\]\(([^\s)]+)\)', text):
            if re.match(r'^[a-zA-Z]+:', target) or target.startswith('#'):
                continue
            target = target.split('#', 1)[0]
            assert (path.parent/target).resolve().exists(), (name, target)
            count += 1
    return dict(local_links_checked=count, display_math_fences=fences,
                controls_absent=True, display_math_fences_balanced=True,
                image_inspection_performed=False)


def calculate():
    saved = read_json(HERE/'results.json')
    before = bindings(saved)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    run = subprocess.run([sys.executable, '-B', '-X', 'utf8', str(HERE/'check.py')],
                         text=True, encoding='utf-8', capture_output=True, env=env)
    assert run.returncode == 0, (run.returncode, run.stdout, run.stderr)
    stdout = json.loads(run.stdout)
    assert stdout['mode'] == 'READ_ONLY_MATCH' and stdout['groups'] == 3
    material, mass, levels, projectors = original_material(saved)
    weak, kappa = independent_weak_calculation(saved, mass, levels, projectors)
    finite = finite_quantifiers(saved, levels, kappa)
    text = text_check()
    assert bindings(saved) == before
    return dict(schema='round1172_independent_review_v1', date='2026-10-13',
                status='PASS_WITH_SCOPE', blockers=0, reviewer_script_sha256=sha(Path(__file__)),
                provenance=before,
                author_default_run=dict(exit_code=run.returncode, stdout=stdout, stderr=run.stderr),
                independent_material=material, independent_weak_source=weak,
                finite_error_quantifiers=finite, document_checks=text,
                independence=dict(blind_review=False,
                    pre_freeze_input='source identity, full-history scope, expectation interpretation, B=0 and error-quantifier suggestions',
                    author_core_edited_by_reviewer=False, author_wrote_review_assets=False,
                    author_1172_module_imported=False,
                    original1112_module_imported=True),
                scientific_scope=dict(new_adopted_axioms=0, new_actual_operations=0,
                    no_complete_contact_stress=True, no_actual_curvature_instrument=True,
                    conclusion='S_eff plus explicit allowed free-leg WL dictionary does not force the exact weak Einstein necessary equation for any fixed G>0 and constant Lambda',
                    not_all_cognitive_extensions_refuted=True,
                    positive_GR_hard_deliverables_completed=0,
                    constructive_goal_completed=False))


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
    path = HERE/'review_checks.json'
    if args.write_once:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        compare(result, read_json(path))
    print(json.dumps(dict(status=result['status'], blockers=0,
        mode='WRITE_ONCE' if args.write_once else 'READ_ONLY_MATCH',
        author_default_exit_code=result['author_default_run']['exit_code'],
        independent_check_exit_code=0,
        projector_error=result['independent_material']['maximum_projector_identity_error'],
        metric_source_error=result['independent_weak_source']['complete_metric_source_error_from_M_over_2']),
        ensure_ascii=False))


if __name__ == '__main__':
    main()
