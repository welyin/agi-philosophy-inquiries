"""1030: independent finite certificates and limited asset verification.

The verifier checks matrix and algebraic certificates.  It does not certify a
physical stable-channel embedding, EFT remainder, or analytic QFT completion.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import ast
import hashlib
import json
import math
import re
from urllib.parse import unquote

import numpy as np
import finite_window_scattering_selection as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / 'research_note_1030.md'
OUT = HERE / 'research_round_1030_checks.json'
TOL = 3e-12
OWN = ['finite_window_scattering_selection.py', 'finite_window_scattering_selection_results.json',
       'drafts/finite_window_scattering_derivation.md', 'review.md', 'verify_round1030.py',
       'selection_audit.md', 'input_dependency_update_v0_19.md', 'NEXT.md']
HISTORICAL = {
    'archive_1009_/research_note_1029.md', 'archive_1009_/1029/NEXT.md',
    'archive_1009_/1029/input_dependency_update_v0_18.md',
    'archive_1009_/research_note_1016.md', 'archive_1009_/research_note_1023.md',
    'archive_554_584/research_note_568.md',
    'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_1009_/1009/input_dependency_ledger_v0_1.md'}
LINK = re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$', re.M)
EXCLUDED = re.compile(r'\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$', re.S | re.M)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_links(text):
    for match in LINK.finditer(EXCLUDED.sub(lambda m: ' ' * len(m[0]), text)):
        value = (match[1] if match[1] is not None else match[2]).strip()
        target = value[1:value.index('>')] if value.startswith('<') and '>' in value else re.split(r'\s+[\"\']', value, 1)[0]
        local = unquote(target.split('#', 1)[0])
        if local and not re.match(r'^[a-zA-Z]+:', local) and not local.startswith('//'):
            yield target, local


def near(x, y):
    assert math.isclose(float(x), float(y), rel_tol=1e-10, abs_tol=TOL), (x, y)


def k_matrix(x, y, tau=1):
    mixing = math.sqrt(3 * float(y) ** 2 / 4)
    if float(y) < 0:
        mixing = -mixing
    return float(tau) * np.array([[float(x), mixing], [mixing, 0.0]])


def spectral_radius(x, y, tau=1):
    x, y = float(x), float(y)
    return float(tau) * (abs(x) + math.hypot(x, math.sqrt(3) * y)) / 2


def distance(rho):
    return (math.hypot(1, 2 * rho) - 1) / 2


def conic(x, y, radius):
    x, y, radius = F(x), F(y), F(radius)
    assert radius >= 0
    return abs(x) <= radius and 3 * y * y <= 4 * radius * (radius - abs(x))


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads(science.OUT.read_text('utf8')))
    assert fresh['round'] == 1030 and fresh['all_scientific_calibrations_passed']
    assert fresh['new_calibration_groups'] == 1 and fresh['cumulative_test_groups'] == 3807
    assert fresh['new_cognitive_axioms'] == 0 and fresh['goal_complete'] is False
    assert fresh['code_sha256'] == sha(HERE / 'finite_window_scattering_selection.py')
    assert 'an actual numerical remainder for P981 or SM' in fresh['not_proved']
    assert 'physical W/Z/h exact asymptotic channel compression' in fresh['not_proved']
    assert 'Goldstone fields are physical channels' in fresh['not_proved']

    projection = fresh['projection']
    assert list(map(F, projection['elastic_isosinglet_polynomial'])) == [F(2), F(0)]
    assert F(projection['exact_elastic_integral']) == 4
    assert F(projection['exact_transition_integral_over_sqrt3']) == 2
    # Integrating constant 2 and 1 over [-1,1], then comparing 1/(64 pi)
    # with tau=s/(16 pi v^2), fixes the normalization without importing code.
    assert F(projection['K_over_t_ww_over_x']) == F(4, 4) == 1
    assert F(projection['K_over_t_wh_over_sqrt3_y']) == F(2, 4)
    assert F(projection['K_over_t_hh']) == 0

    region = fresh['exact_finite_window_region']
    delta, radius, r = F(region['delta']), F(region['rho_limit']), F(region['R'])
    assert delta == F(1, 15) and radius == F(4, 15) and r == 2 * radius
    assert radius ** 2 == delta * (1 + delta) == F(region['rho_limit_squared'])
    cases = region['cases']; assert len(cases) == 7
    by_name = {}
    for case in cases:
        x, y = F(case['x']), F(case['y'])
        by_name[case['name']] = case
        assert F(case['a_squared']) == 1 - x and F(case['b']) == 1 - x - y
        margin = 4 * radius * (radius - abs(x)) - 3 * y * y
        assert margin == F(case['exact_region_margin'])
        assert bool(case['feasible']) == conic(x, y, radius)
        rho = spectral_radius(x, y)
        near(case['rho'], rho); near(case['minimal_complex_remainder'], distance(rho))
    assert by_name['nonzero_curvature_allowed']['feasible'] is True
    assert by_name['elastic_zero_but_joint_rejected']['feasible'] is False
    assert by_name['positive_x_outside']['feasible'] is False
    assert F(by_name['positive_x_boundary']['exact_region_margin']) == 0
    assert F(by_name['negative_x_boundary']['exact_region_margin']) == 0
    assert F(3, 100) < radius * radius < F(3, 25)
    grid = 0
    for i in range(-6, 7):
        for j in range(-6, 7):
            x, y = F(i, 15), F(j, 15)
            assert conic(x, y, radius) == (spectral_radius(x, y) <= float(radius) + TOL)
            grid += 1
    assert grid == region['exact_grid_points'] == 169
    # The zero-budget degeneration must retain the |x| condition.
    assert conic(0, 0, 0)
    assert not conic(F(1, 5), 0, 0) and not conic(0, F(1, 5), 0)
    zero = region['zero_budget_boundary']
    assert zero['origin_allowed'] is True and zero['nonzero_x_allowed'] is False
    assert F(zero['nonzero_x']) == F(1, 5) and F(zero['y']) == 0
    assert F(zero['quadratic_margin_alone']) == 0

    polar = fresh['polar_and_spectrum']
    assert len(polar['points']) == 5
    for point in polar['points']:
        x, y = F(point['x']), F(point['y'])
        assert list(map(F, point['exact_characteristic_coefficients'])) == [1, -x, -F(3, 4) * y * y]
        k = k_matrix(x, y)
        rho = spectral_radius(x, y)
        # Independent polar algorithm: science diagonalizes K; use SVD of I+2iK.
        raw = np.eye(2) + 2j * k
        left, singular, right = np.linalg.svd(raw)
        s = left @ right
        amplitude = (s - np.eye(2)) / (2j)
        assert np.linalg.norm(s.conj().T @ s - np.eye(2), 2) < TOL
        near(np.linalg.norm(amplitude - k, 2), distance(rho))
        near(max(singular), math.sqrt(1 + 4 * rho * rho))
        near(point['actual_distance'], distance(rho))
        near(point['sharp_distance'], distance(rho))
        near(point['S_transition_probability'], abs(s[1, 0]) ** 2)
        assert np.max(np.abs(np.sort(np.asarray(point['eigenvalues'])) - np.linalg.eigvalsh(k))) < TOL
    rho2 = F(polar['exact_allowed_example_rho_squared'])
    assert rho2 == F(3, 100)
    assert 4 * rho2 / (1 + 4 * rho2) == F(polar['exact_allowed_example_polar_conversion']) == F(3, 28)
    assert polar['no_QFT_completion_claim'] is True
    for name in ('max_unitarity_residual', 'max_sharp_distance_residual', 'max_spectrum_residual'):
        assert 0 <= polar[name] < TOL

    compressed = fresh['unitary_compression']
    assert compressed['samples'] == 24 and compressed['full_dimensions'] == [3, 4, 5, 6]
    rng = np.random.default_rng(1030)
    maximum, minimum_gap = 0.0, math.inf
    for i in range(24):
        n = 3 + i % 4
        raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        full, _ = np.linalg.qr(raw)
        c = full[:2, :2]
        a = (c - np.eye(2)) / (2j)
        k = k_matrix(F(i % 7 - 3, 11), F(i % 9 - 4, 13), F(3, 5))
        gap = np.linalg.norm(a - k, 2) - distance(np.linalg.norm(k, 2))
        assert np.linalg.norm(c, 2) <= 1 + TOL and gap >= -TOL
        # Direct optical inequality; omitted channels contribute positive loss.
        imaginary = (a - a.conj().T) / (2j)
        loss = imaginary - a.conj().T @ a
        assert np.min(np.linalg.eigvalsh(loss)) >= -TOL
        maximum = max(maximum, float(np.linalg.norm(c, 2)))
        minimum_gap = min(minimum_gap, float(gap))
    near(compressed['max_compressed_S_norm'], maximum)
    near(compressed['minimum_distance_minus_proved_lower_bound'], minimum_gap)

    window = fresh['nonconstant_budget_window']
    assert list(map(F, window['g_polynomial_ascending'])) == [F(9, 20), F(-1), F(1)]
    assert F(window['exact_minimizer']) == F(1, 2)
    assert F(window['exact_minimum_rho_over_t_limit']) == F(1, 5)
    test_square = F(window['test_rho_over_t_squared']); assert test_square == F(3, 64)
    statuses = []
    for sample in window['samples']:
        tau = F(sample['t'])
        g = F(1, 5) + (tau - F(1, 2)) ** 2
        assert g == F(sample['g']) and g * g == F(sample['g_squared'])
        near(float(sample['delta']) * (1 + float(sample['delta'])), float(tau * tau * g * g))
        passes = test_square <= g * g
        assert sample['passes'] is passes
        statuses.append(passes)
    assert statuses == [True, False, True]
    assert window['actual_EFT_budget_certified'] is False

    curvature = fresh['scalar_target_curvature']
    assert len(curvature['local_sqrt_checks']) == 4
    for case in curvature['local_sqrt_checks']:
        a, b = F(case['a']), F(case['b'])
        assert list(map(F, case['sqrt_F_series'])) == [1, a, (b - a * a) / 2]
        assert F(case['v_squared_K_pipi']) == 1 - a * a
        assert F(case['v_squared_K_pih']) == a * a - b
    nonlinear = curvature['flat_at_origin_not_globally_flat']
    u, c = F(nonlinear['evaluation_h_over_v']), F(nonlinear['cubic_coefficient'])
    f, fp, fpp = (1 + u) ** 2 + c * u ** 3, 2 * (1 + u) + 3 * c * u ** 2, 2 + 6 * c * u
    assert u == F(1, 5) and c == 1
    assert f == F(nonlinear['F']) and fp == F(nonlinear['F_prime_u']) and fpp == F(nonlinear['F_second_u'])
    assert F(nonlinear['v_squared_K_pipi']) == (4 * f - fp * fp) / (4 * f * f) != 0
    assert F(nonlinear['v_squared_K_pih']) == (fp * fp - 2 * f * fpp) / (4 * f * f) != 0

    quantifier = fresh['error_ball_quantifier']
    d = F(quantifier['delta'])
    assert quantifier['K_zero'] is True and d == F(1, 15)
    assert F(quantifier['error_norm']) == d
    assert F(quantifier['exact_S_norm']) == 1 + 2 * d == F(17, 15) > 1

    assert set(fresh['historical_source_sha256']) == HISTORICAL
    for name, digest in fresh['historical_source_sha256'].items():
        assert sha(BASE / name) == digest, name
    assets = [NOTE] + [HERE / name for name in OWN]
    assert len(assets) == 9 and not any(p.name in {'README.md', 'research_direction.md', 'RESEARCH_STATE.md'} for p in assets)
    links = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == '.py':
            ast.parse(path.read_text('utf-8-sig'))
        if path.suffix == '.md':
            text = path.read_text('utf-8-sig')
            assert text.count('$$') % 2 == 0, path
            for target, local in local_links(text):
                resolved = (path.parent / local.replace('\\', '/')).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (path, target)
                links += 1
    assert all(f'## {n}.' in NOTE.read_text('utf8') for n in range(1, 11))
    assert '独立科学签审通过' in (HERE / 'review.md').read_text('utf8')
    return dict(round=1030, date='2026-10-08', all_delivery_checks_passed=True,
                scientific_result_reproduced=True, new_calibration_groups=1,
                cumulative_research_groups=3807, new_cognitive_axioms=0,
                calibration_kind='finite-window coupled-channel contraction distance',
                named_region_cases=7, exact_region_grid_points=169,
                polar_spectrum_cases=5, full_unitary_compressions=24,
                independent_polar_algorithm='SVD of I+2iK, compared with scientific eigendecomposition',
                zero_budget_requires_both_x_y_zero=True,
                allowed_example_transition_probability='3/28',
                nonconstant_budget_endpoint_status=[True, False, True],
                wrong_error_ball_remainder_S_norm='17/15',
                complete_complex_error_required=True, omitted_channels_allowed=True,
                theorem_is_analytic_not_sample_inference=True,
                pointwise_matrix_completion_only=True,
                actual_SM_stable_channel_bridge_certified=False,
                actual_EFT_remainder_certified=False, analytic_QFT_completion_certified=False,
                scalar_target_curvature_is_not_spacetime_curvature=True,
                goal_completed=False, live_navigation_frozen=False, neighboring_round_frozen=False,
                visual_checks_performed=False, frozen_current_files=9,
                historical_input_files=len(HISTORICAL), local_links_checked=links,
                historical_source_sha256=fresh['historical_source_sha256'],
                source_sha256={str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in assets})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', '--write-receipt', dest='write', action='store_true')
    args = parser.parse_args()
    result = verify(prospective=args.write)
    if args.write:
        with OUT.open('x', encoding='utf8') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    else:
        assert result == json.loads(OUT.read_text('utf8'))
    print(json.dumps({k: v for k, v in result.items() if not k.endswith('sha256')}, ensure_ascii=False, indent=2))
