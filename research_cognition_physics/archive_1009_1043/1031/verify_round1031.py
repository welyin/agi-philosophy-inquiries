"""1031: independent finite certificates and limited delivery verification.

This verifies a leading classical common-parent branch and its finite freedom.
It does not verify an interacting quantum construction or actual EFT errors.
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
import common_parent_coupling_freedom as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / 'research_note_1031.md'
OUT = HERE / 'research_round_1031_checks.json'
TOL = 2e-11
OWN = ['common_parent_coupling_freedom.py', 'common_parent_coupling_freedom_results.json',
       'drafts/common_parent_coupling_freedom_derivation.md', 'review.md', 'verify_round1031.py',
       'selection_audit.md', 'input_dependency_update_v0_20.md', 'NEXT.md']
HISTORICAL = {
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_990_1008/993/common_candidate_v1.md',
    'archive_1009_/1009/input_dependency_ledger_v0_1.md',
    'archive_1009_/research_note_1009.md', 'archive_1009_/research_note_1015.md',
    'archive_1009_/research_note_1024.md', 'archive_1009_/research_note_1027.md',
    'archive_1009_/research_note_1029.md', 'archive_1009_/research_note_1030.md',
    'archive_554_584/research_note_572.md', 'archive_554_584/research_note_573.md',
    'archive_742_763/research_note_753.md'}
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


def maxabs(x):
    return float(np.max(np.abs(x)))


def comm(a, b):
    return a @ b - b @ a


def independent_curvature(g, f, z):
    # Construct Pauli matrices independently rather than use scientific generators.
    pauli = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]],
                      [[1, 0], [0, -1]]], dtype=complex) / 2
    t = np.zeros((3, 3, 3), complex)
    t[:, :2, :2] = pauli
    assert maxabs(comm(t[0], t[1]) - 1j * t[2]) == 0
    a = np.zeros((4, 3, 3), complex)
    a[1:] = f * t
    field = np.zeros((4, 4, 3, 3), complex)
    for i in range(3):
        field[0, i + 1] = z * t[i]
        field[i + 1, 0] = -z * t[i]
    for i in range(1, 4):
        for j in range(1, 4):
            field[i, j] = -1j * g * comm(a[i], a[j])
    # Conformal cancellation reduces the full Euler tensor to this matrix form.
    gauss = sum((comm(a[i], field[i, 0]) for i in range(1, 4)), np.zeros((3, 3), complex))
    euler = []
    acceleration = -2 * g * g * f ** 3
    for i in range(1, 4):
        spatial = sum((comm(a[j], field[j, i]) for j in range(1, 4)), np.zeros((3, 3), complex))
        euler.append(-acceleration * t[i - 1] - 1j * g * spatial)
    assert maxabs(gauss) < TOL and maxabs(euler) < TOL
    assert maxabs(field[:, :, 2, :]) == 0 and maxabs(field[:, :, :, 2]) == 0
    return field


def independent_stress(field, scale):
    metric = np.diag([-1., 1., 1., 1.])
    signs = np.diag(metric)
    contraction = sum(signs[i] * signs[j] * 2 * np.trace(field[i, j] @ field[i, j]).real
                      for i in range(4) for j in range(4))
    tensor = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            tensor[i, j] = sum(signs[k] * 2 * np.trace(field[i, k] @ field[j, k]).real
                               for k in range(4)) - metric[i, j] * contraction / 4
    return tensor / scale ** 2, contraction / scale ** 4


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads(science.OUT.read_text('utf8')))
    assert fresh['round'] == 1031 and fresh['all_scientific_calibrations_passed']
    assert fresh['new_calibration_groups'] == 1 and fresh['cumulative_test_groups'] == 3808
    assert fresh['new_cognitive_axioms'] == 0 and fresh['goal_complete'] is False
    assert fresh['code_sha256'] == sha(HERE / 'common_parent_coupling_freedom.py')
    assert 'Not a full quantum/cognitive countermodel' in fresh['scope']

    su = fresh['su3_structure']
    assert su['embedded_generators'] == [1, 2, 3]
    for key in ('outside_subalgebra_residual', 'matrix_commutator_residual', 'jacobi_residual'):
        assert 0 <= su[key] < TOL
    assert len(su['coupled_representation_cases']) == 10
    for case in su['coupled_representation_cases']:
        near(case['cubic_coefficient'], case['g'])
        near(case['quartic_coefficient'], case['g'] ** 2)
        assert case['covariant_algebra_residual'] < TOL
    assert {case['representation'] for case in su['coupled_representation_cases']} == {'fundamental', 'antifundamental'}

    tensors = fresh['full_field_tensors']
    assert tensors['samples'] == len(tensors['cases']) == 5
    assert tensors['full_matrix_color_equations_per_sample'] == 32
    planck_squares = [2., 2., 3., 2.5, 4.]
    tensor_max = 0.0
    for data, mp2 in zip(tensors['cases'], planck_squares):
        g, f, z, a = (data[k] for k in ('g', 'f', 'f_prime', 'scale_factor'))
        energy = z * z + g * g * f ** 4
        near(data['energy_invariant'], energy)
        field = independent_curvature(g, f, z)
        tensor, field_square = independent_stress(field, a)
        expected = np.diag([1.5, .5, .5, .5]) * energy / (a * a)
        assert maxabs(tensor - expected) < TOL
        rho = tensor[0, 0] / (a * a)
        near(data['source_density'], rho)
        near(data['electric_density'], 1.5 * z * z / a ** 4)
        near(data['magnetic_density'], 1.5 * g * g * f ** 4 / a ** 4)
        near(data['magnetic_fraction'], g * g * f ** 4 / energy)
        near(data['magnetic_fraction'], .5 * (1 + field_square / (4 * rho)))
        # Vacuum energy is kept in the complete Einstein equation with its sign.
        vacuum, bare = data['vacuum_energy'], data['bare_lambda']
        effective = bare + vacuum / mp2
        near(data['lambda_effective'], effective)
        ap2 = energy / (2 * mp2) + effective * a ** 4 / 3
        app = 2 * effective * a ** 3 / 3
        einstein = np.diag([3 * ap2 / a ** 2] + [ap2 / a ** 2 - 2 * app / a] * 3)
        metric = a * a * np.diag([-1., 1., 1., 1.])
        assert maxabs(einstein + bare * metric - (tensor - vacuum * metric) / mp2) < TOL
        for key, value in data.items():
            if key.endswith(('_max', '_residual')):
                assert 0 <= value < TOL, key
                tensor_max = max(tensor_max, value)
    bad = tensors['wrong_acceleration_control']
    assert bad['color_matrix_euler_max'] > .04 and bad['color_euler_max'] > .08
    assert bad['stress_divergence_max'] > .1
    assert bad['gauss_max'] < TOL and bad['bianchi_max'] < TOL

    certificate = fresh['rational_finite_certificate']
    p, t = F(certificate['common_p']), F(certificate['conformal_endpoint'])
    g1, g2 = map(F, certificate['g_interval'])
    assert (p, t, g1, g2) == (1, F(1, 4), 1, 2)
    q = g2 ** 2 * p ** 2 * t ** 4
    assert q == F(certificate['bootstrap_condition_g_squared_p_squared_T_fourth']) == F(1, 64)
    m1 = g1 ** 2 * p ** 2 * t ** 4
    m2 = g2 ** 2 * (p * t * (1 - q / 10)) ** 4 / p ** 2
    assert m1 == F(certificate['endpoint_M1_exact_upper'])
    assert m2 == F(certificate['endpoint_M2_exact_lower'])
    assert m2 - m1 == F(certificate['exact_gap_bound_before_Bernoulli'])
    gap = g2 ** 2 * p ** 2 * t ** 4 * (1 - 4 * q / 10) - m1
    assert gap == F(certificate['simple_strict_gap']) == F(119, 10240) > 0
    assert m2 - m1 >= gap
    lower = 2 * g1 * p ** 2 * t ** 4 * F(637, 640) * F(127, 128)
    upper = 2 * g2 * p ** 2 * t ** 4
    assert F(639, 640) ** 3 >= F(637, 640)
    assert lower == F(certificate['derivative_lower_at_T']) == F(80899, 10485760)
    assert upper == F(certificate['derivative_upper_on_full_interval']) == F(1, 64)
    interval = certificate['finite_constraint_example']
    center, bound = F(interval['center']), F(interval['bound'])
    center_upper = center ** 2 * t ** 4
    margin, radius = bound - center_upper, (bound - center_upper) / (2 * upper)
    assert center_upper == F(interval['center_prediction_upper']) == F(9, 1024)
    assert margin == F(interval['certified_margin']) == F(7, 1024)
    assert upper == F(interval['uniform_Lipschitz_constant'])
    assert radius == F(interval['half_margin_radius']) == F(7, 32)
    assert list(map(F, interval['surviving_interval'])) == [center - radius, center + radius] == [F(41, 32), F(55, 32)]
    assert g1 < center - radius < center + radius < g2
    transport = certificate['conditional_error_transport']
    hypothetical = F(transport['hypothetical_each_prediction_error'])
    assert hypothetical == F(1, 1024)
    assert gap - 2 * hypothetical == F(transport['surviving_prediction_gap']) == F(99, 10240) > 0
    assert transport['actual_EFT_error_certified'] is False

    trajectories = fresh['numerical_trajectories']
    assert trajectories['samples'] == len(trajectories['rows']) == 9
    assert (trajectories['coarse_steps'], trajectories['fine_steps']) == (128, 256)
    assert trajectories['numerical_integration_is_rigorous_error_bound'] is False
    assert F(trajectories['exact_common_proper_time_at_endpoint']) == t + t * t / 4 == F(17, 64)
    scale, adot = F(1) + t / 2, F(1, 2)
    rho = F(3, 2) / scale ** 4
    assert scale == F(9, 8) and F(trajectories['exact_common_density_at_endpoint']) == rho
    numerical_m = []
    max_cross_difference, max_energy_error = 0.0, 0.0
    for index, row in enumerate(trajectories['rows']):
        g, f, z = (row[k] for k in ('g', 'f', 'f_prime'))
        near(g, 1 + index / 8)
        assert float(p * t - F(g) ** 2 * p ** 3 * t ** 5 / 10) - TOL <= f <= float(p * t) + TOL
        assert float(p - F(g) ** 2 * p ** 3 * t ** 4 / 2) - TOL <= z <= float(p) + TOL
        near(z * z + g * g * f ** 4, 1)
        m = g * g * f ** 4
        near(row['magnetic_fraction'], m)
        # Scaling identity eliminates the scaled y to give 2*g*T*f^3*f'.
        derivative = 2 * g * float(t) * f ** 3 * z
        near(row['d_g_M_from_scaling'], derivative)
        assert float(lower) <= derivative <= float(upper)
        near(row['common_a'], scale); near(row['common_a_prime'], adot)
        near(row['common_rho_exact'], rho)
        assert row['coarse_fine_difference'] < 5e-11
        assert row['scaling_difference'] < 5e-12 and row['numerical_energy_residual'] < 5e-12
        max_cross_difference = max(max_cross_difference, row['coarse_fine_difference'], row['scaling_difference'])
        max_energy_error = max(max_energy_error, row['numerical_energy_residual'])
        full = row['full_equations_at_numerical_endpoint']
        for key, value in full.items():
            if key.endswith(('_max', '_residual')):
                assert 0 <= value < TOL
        near(full['source_density'], rho); near(full['magnetic_fraction'], m)
        numerical_m.append(m)
    assert all(a < b for a, b in zip(numerical_m, numerical_m[1:]))
    near(trajectories['endpoint_numerical_gap'], numerical_m[-1] - numerical_m[0])
    assert numerical_m[-1] - numerical_m[0] > float(gap)

    mass = fresh['anomaly_and_fixed_tree_mass']
    charges = list(map(F, mass['hypercharges']))
    multiplicities = mass['multiplicities']
    assert charges == [F(1, 6), F(-2, 3), F(1, 3), F(-1, 2), F(1)]
    assert multiplicities == [6, 3, 3, 2, 1]
    independent_anomalies = [sum(n * y for n, y in zip(multiplicities, charges)),
                            sum(n * y ** 3 for n, y in zip(multiplicities, charges)),
                            2 * charges[0] + charges[1] + charges[2],
                            3 * charges[0] + charges[3], F(2 - 1 - 1)]
    assert independent_anomalies == [0] * 5
    assert len(mass['exact_anomaly_coefficients']) == 5
    assert all(F(value) == 0 for value in mass['exact_anomaly_coefficients'].values())
    assert mass['total_SU2_doublets'] == 3 * (3 + 1) == 12 and mass['SU2_global_parity'] == 0
    flavor = mass['exact_flavor_certificate']
    rotation = [[F(value) for value in row] for row in flavor['mixing_matrix']]
    for i in range(3):
        for j in range(3):
            assert rotation[i][j] == F(i == j) - F((i + 1) * (j + 1), 7)
            assert sum(rotation[i][k] * rotation[j][k] for k in range(3)) == (i == j)
    hu = list(map(F, flavor['Hu_diagonal']))
    assert hu == [F(1, 25), F(4, 25), F(16, 25)] and len(set(hu)) == 3
    hd = [[sum(rotation[i][k] * F((k + 1) ** 2, 49) * rotation[j][k] for k in range(3))
           for j in range(3)] for i in range(3)]
    assert hd == [[F(value) for value in row] for row in flavor['Hd']]
    assert [hd[0][1], hd[0][2], hd[1][2]] == list(map(F, flavor['Hd_upper_off_diagonal'])) == [F(18, 343), F(12, 343), F(6, 343)]
    # Distinct Hu fixes a diagonal commutant, and the complete nonzero Hd graph
    # fixes one scalar.  This exact upper bound supplements numerical rank.
    assert mass['quark_common_Hermitian_commutant_dimension'] == 1
    expected_masses = ([k / (5 * math.sqrt(2)) for k in (1, 2, 4) for _ in range(6)]
                       + [k / (7 * math.sqrt(2)) for k in (1, 2, 3) for _ in range(6)]
                       + [k / (11 * math.sqrt(2)) for k in (1, 2, 3) for _ in range(2)]
                       + [k / 26 for k in (1, 2, 4)])
    assert len(expected_masses) == mass['full_tree_Weyl_mass_rank'] == 45
    assert maxabs(np.sort(expected_masses) - np.sort(mass['mass_singular_values'])) < TOL
    assert maxabs(np.sort(mass['neutral_Takagi_masses']) - np.array([1, 2, 4]) / 26) < TOL
    assert mass['residual_charge_mass_Ward'] < TOL
    assert len(mass['not_identical']) == 4 and 'all-order pole masses' in mass['not_identical']

    ledger = {row['id']: row for row in fresh['common_constraint_ledger']}
    assert len(ledger) == 12
    for key in ('representation_anomaly', 'vector_vertices', 'matter_Noether', 'universal_gravity', 'tree_flavor_mass'):
        assert ledger[key]['status'] == 'exact_preserved'
    for key in ('other_classical_fields', 'finite_same_geometry'):
        assert ledger[key]['status'] == 'exact_solution_branch'
    assert ledger['linear_Higgs_coefficients']['status'] == 'exact_preserved_only'
    assert ledger['high_order_matching']['status'] == 'not_included'
    for key in ('quantum_state_source', 'QCD_running_confinement', 'FUCP_and_instruments'):
        assert ledger[key]['status'] == 'not_proved'
    assert '1030 stable physical channel map' in ledger['linear_Higgs_coefficients']['not_claimed']

    assert set(fresh['historical_source_sha256']) == HISTORICAL
    for name, digest in fresh['historical_source_sha256'].items():
        assert sha(BASE / name) == digest, name
    assets = [NOTE] + [HERE / name for name in OWN]
    assert len(assets) == 9 and not any(path.name in {'README.md', 'research_direction.md', 'RESEARCH_STATE.md'} for path in assets)
    links = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == '.py':
            ast.parse(path.read_text('utf-8-sig'))
        if path.suffix == '.md':
            text = path.read_text('utf-8-sig')
            assert text.count('$$') % 2 == 0, path
            assert not any(ord(c) < 32 and c not in '\n\r\t' for c in text), path
            for target, local in local_links(text):
                resolved = (path.parent / local.replace('\\', '/')).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (path, target)
                links += 1
    assert all(f'## {n}.' in NOTE.read_text('utf8') for n in range(1, 11))
    assert '独立科学签审通过' in (HERE / 'review.md').read_text('utf8')
    return dict(round=1031, date='2026-10-08', all_delivery_checks_passed=True,
                scientific_result_reproduced=True, new_calibration_groups=1,
                cumulative_research_groups=3808, new_cognitive_axioms=0,
                calibration_kind='common leading classical parent coupling freedom',
                full_tensor_samples=5, full_color_components_per_sample=32,
                trajectory_samples=9, coupled_representation_cases=10,
                all_color_matrix_entries_checked=True,
                independent_curvature_and_stress_reconstruction=True,
                complete_Einstein_and_vacuum_sign_check=True,
                maximum_tensor_calibration_residual=tensor_max,
                maximum_RK_and_scaling_cross_difference=max_cross_difference,
                maximum_numerical_energy_residual=max_energy_error,
                numerical_integration_is_rigorous_error_bound=False,
                exact_finite_gap='119/10240', derivative_lower='80899/10485760',
                derivative_upper='1/64', finite_constraint_interval=['41/32', '55/32'],
                conditional_error_gap='99/10240', common_a_endpoint='9/8',
                common_proper_time_endpoint='17/64', fixed_tree_mass_rank=45,
                exact_flavor_commutant_dimension=1, exact_anomaly_coefficients_vanish=True,
                classical_Weyl_zero_is_not_quantum_vacuum=True,
                same_tree_mass_is_not_same_full_quantum_spectrum=True,
                actual_EFT_error_certified=False, full_quantum_parent_certified=False,
                round1030_physical_channel_bridge_certified=False,
                full_cognitive_countermodel_certified=False,
                theorem_is_analytic_not_sample_inference=True,
                goal_completed=False, live_navigation_frozen=False, neighboring_round_frozen=False,
                visual_checks_performed=False, frozen_current_files=9,
                historical_input_files=len(HISTORICAL), local_links_checked=links,
                historical_source_sha256=fresh['historical_source_sha256'],
                source_sha256={str(path.relative_to(ROOT)).replace('\\', '/'): sha(path) for path in assets})


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
    print(json.dumps({key: value for key, value in result.items() if not key.endswith('sha256')}, ensure_ascii=False, indent=2))
