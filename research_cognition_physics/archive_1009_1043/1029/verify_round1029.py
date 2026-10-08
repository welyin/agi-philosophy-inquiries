"""1029: independent finite certificates and limited frozen-asset verification.

Finite coefficient checks do not prove completeness of arbitrary local gauge
generators.  That step is the regular Euler-jet/Noether argument in the note.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import argparse
import ast
import hashlib
import json
import re
from urllib.parse import unquote

import first_order_lift_equivalence as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent / 'research_note_1029.md'
OUT = HERE / 'research_round_1029_checks.json'
OWN = ['first_order_lift_equivalence.py', 'first_order_lift_equivalence_results.json',
       'drafts/first_order_lift_derivation.md', 'review.md', 'verify_round1029.py',
       'selection_audit.md', 'input_dependency_update_v0_18.md', 'NEXT.md']
HISTORICAL = {
    'archive_1009_/research_note_1028.md',
    'archive_1009_/1028/input_dependency_update_v0_17.md',
    'archive_1009_/1028/NEXT.md', 'archive_1009_/research_note_1027.md',
    'archive_1009_/1027/charged_source_selection.py',
    'archive_1009_/research_note_1018.md',
    'archive_342_369/research_note_358.md',
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_990_1008/993/common_candidate_v1.md'}
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


def det(matrix):
    """Leibniz determinant, independent of science's elimination."""
    n = len(matrix)
    total = F(0)
    for perm in permutations(range(n)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        term = F(sign)
        for i, j in enumerate(perm):
            term *= F(matrix[i][j])
        total += term
    return total


def matrix(saved):
    return [[F(x) for x in row] for row in saved]


def scalar_product(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def polynomial(saved):
    result = {}
    for key, value in saved.items():
        powers = ast.literal_eval(key)
        assert isinstance(powers, tuple) and len(powers) == 4
        assert all(isinstance(p, int) and p >= 0 for p in powers)
        result[powers] = F(value)
    return result


def differentiated(poly, *axes):
    result = {}
    for powers, value in poly.items():
        powers = list(powers)
        for axis in axes:
            value *= powers[axis]
            powers[axis] = max(0, powers[axis] - 1)
        if value:
            key = tuple(powers)
            result[key] = result.get(key, F(0)) + value
    return {p: c for p, c in result.items() if c}


def linear_combination(terms):
    result = {}
    for coefficient, poly in terms:
        for powers, value in poly.items():
            result[powers] = result.get(powers, F(0)) + coefficient * value
    return {p: c for p, c in result.items() if c}


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads(science.OUT.read_text('utf8')))
    assert fresh['round'] == 1029 and fresh['all_scientific_calibrations_passed']
    assert fresh['new_calibration_groups'] == 1 and fresh['cumulative_test_groups'] == 3806
    assert fresh['new_cognitive_axioms'] == 0 and fresh['goal_complete'] is False
    assert fresh['code_sha256'] == sha(HERE / 'first_order_lift_equivalence.py')
    assert 'bare off-shell lift uniqueness' in fresh['not_proved']
    assert 'uniqueness or generation of S1' in fresh['not_proved']
    kinetic = fresh['kinetic_principal_blocks']
    assert kinetic['bosonic_velocity_variables'] == 36 + 4 + 1 == 41
    assert kinetic['bosonic_hessian_rank'] == 41 and F(kinetic['bosonic_hessian_determinant']) == 1
    assert kinetic['independent_Weyl_modules_with_internal_and_flavor_multiplicity'] == 45
    weyl = matrix(kinetic['weyl_time_symbol_realified'])
    assert weyl == [[0, 0, -1, 0], [0, 0, 0, -1], [1, 0, 0, 0], [0, 1, 0, 0]]
    assert det(weyl) == F(kinetic['weyl_time_symbol_real_determinant']) == 1
    assert kinetic['weyl_time_symbol_real_rank'] == 4

    pf = fresh['PF_principal_and_constraint_blocks']
    expected = [[F(int(i == j) - 1) if i < 3 and j < 3 else F(2 * int(i == j))
                 for j in range(6)] for i in range(6)]
    assert matrix(pf['spatial_velocity_hessian']) == expected
    assert det(expected) == F(pf['spatial_velocity_determinant']) == -16
    assert pf['spatial_velocity_rank'] == 6
    pivots = pf['constraint_pivots']
    constraints = [{p: F(c) for p, c in row.items()} for row in pf['constraint_polynomials']]
    # Hamiltonian: off-diagonal 2 d_i d_j h_ij, minus all crossed trace terms.
    expected_hamiltonian = {'h01_d01': F(2), 'h02_d02': F(2), 'h12_d12': F(2)}
    for i in range(3):
        for j in range(3):
            if i != j:
                expected_hamiltonian[f'h{i}{i}_d{j}{j}'] = F(-1)
    assert constraints[0] == expected_hamiltonian
    for i in range(3):
        assert constraints[i + 1] == {f'pi{min(i,j)}{max(i,j)}_d{j}': F(1) for j in range(3)}
    pivot_matrix = [[row.get(p, F(0)) for p in pivots] for row in constraints]
    assert pivot_matrix == matrix(pf['constraint_pivot_matrix'])
    assert det(pivot_matrix) == F(pf['constraint_pivot_determinant']) == -1

    symbol_cases = fresh['PF_parameter_symbol']['cases']
    assert len(symbol_cases) == 3
    for case in symbol_cases:
        k = list(map(F, case['k_covector']))
        pairs = [(mu, nu) for mu in range(4) for nu in range(mu, 4)]
        symbol = [[k[mu] * (c == nu) + k[nu] * (c == mu) for c in range(4)] for mu, nu in pairs]
        assert symbol == matrix(case['symbol_matrix'])
        rows = [pairs.index(tuple(p)) for p in case['selected_rows']]
        columns = case['column_order']
        minor = [[symbol[row][c] for c in columns] for row in rows]
        a = columns[0]
        assert det(minor) == F(case['minor_determinant']) == 2 * k[a] ** 4 != 0
        assert case['rank'] == 4
    null_k = list(map(F, symbol_cases[1]['k_covector']))
    assert -null_k[0] ** 2 + sum(x * x for x in null_k[1:]) == 0

    lift = fresh['genuine_higgs_trivial_lift']
    h, hy = list(map(F, lift['h_real'])), list(map(F, lift['partial_y_h']))
    euler = list(map(F, lift['actual_higgs_Euler']))
    force = [F(lift['lambda_h']) * (scalar_product(h, h) - F(lift['v_squared'])) * v for v in h]
    assert force == list(map(F, lift['potential_gradient']))
    assert [-F(acc) - v for acc, v in zip(lift['partial_t_squared_h'], force)] == euler
    ty = matrix(lift['tY'])
    expected_ty = [[0, 0, F(-1, 2), 0], [0, 0, 0, F(-1, 2)],
                   [F(1, 2), 0, 0, 0], [0, F(1, 2), 0, 0]]
    assert ty == expected_ty and all(ty[i][j] == -ty[j][i] for i in range(4) for j in range(4))
    b = scalar_product(h, hy)
    assert b == F(lift['parameter_scalar_b']) == F(5, 14)
    delta = [b * scalar_product(row, euler) for row in ty]
    assert delta == list(map(F, lift['homogeneous_lift_delta_h']))
    assert scalar_product(euler, delta) == F(lift['Euler_contraction']) == 0
    assert scalar_product(h, delta) == F(lift['off_shell_delta_O']) == F(-75, 224)

    transport = fresh['full_E_recoil_transport']
    f, p, q = F(transport['radial_f']), F(transport['radial_p']), F(transport['source_q'])
    potential = F(2, 5) * (f * f - F(9, 4)) ** 2 / 4
    trace = -p * p - 4 * potential
    recoil, e1 = -q * trace / 2, q * trace / 2
    response = f ** 3 * p
    assert potential == F(transport['radial_potential']) == F(121, 2560)
    assert trace == F(transport['full_matter_stress_trace']) == F(-8489, 31360)
    assert recoil == F(transport['free_PF_Euler_linear_recoil_coefficient']) == -e1 != 0
    assert e1 == F(transport['first_interaction_Euler_coefficient'])
    assert response == F(transport['mixed_trivial_operator_O_factor']) == F(125, 224)
    assert F(transport['old_epsilon_M_E0_O']['eps^2']) == response * recoil == F(212225, 1404928)
    inherited = transport['inherited_minus_epsilon_squared_M_E1_and_higher']
    removed = transport['removed_epsilon_M_Efull_O']
    assert F(inherited['eps^2']) == -response * e1 == response * recoil
    assert F(inherited['eps^3*c']) == -response and F(removed['eps^3*c']) == response
    assert set(removed) == {'eps^3*c'}
    contractions = list(map(F, transport['mixed_Euler_contractions']))
    assert contractions == [F(-15, 286), F(15, 286)] and sum(contractions) == 0

    improvement = fresh['identically_conserved_improvement']
    fpoly, opoly = polynomial(improvement['f_polynomial']), polynomial(improvement['O_polynomial'])
    square = {}
    for p1, v1 in fpoly.items():
        for p2, v2 in fpoly.items():
            pp = tuple(x + y for x, y in zip(p1, p2))
            square[pp] = square.get(pp, F(0)) + v1 * v2 / 2
    assert opoly == square
    signs = [-1, 1, 1, 1]
    box = linear_combination([(signs[mu], differentiated(opoly, mu, mu)) for mu in range(4)])
    assert box == polynomial(improvement['box_O_polynomial'])
    current = [[polynomial(poly) for poly in row] for row in improvement['I_upper_polynomials']]
    for mu in range(4):
        for nu in range(4):
            expected_i = linear_combination([(signs[mu] * signs[nu], differentiated(opoly, mu, nu)),
                                            (-signs[mu] if mu == nu else 0, box)])
            assert current[mu][nu] == expected_i == current[nu][mu]
    for nu in range(4):
        assert linear_combination([(1, differentiated(current[mu][nu], mu)) for mu in range(4)]) == {}
    origin = [[poly.get((0, 0, 0, 0), F(0)) for poly in row] for row in current]
    assert origin == matrix(improvement['I_upper_at_origin']) == [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
    assert improvement['exact_divergence_polynomials'] == [{}, {}, {}, {}]
    assert improvement['no_equations_of_motion_used'] is True

    assert set(fresh['historical_source_sha256']) == HISTORICAL
    for name, digest in fresh['historical_source_sha256'].items():
        assert sha(BASE / name) == digest, name
    assets = [NOTE] + [HERE / name for name in OWN]
    assert len(assets) == 9 and not any(p.name in {'README.md', 'research_direction.md', 'RESEARCH_STATE.md'} for p in assets)
    links = 0
    for path in assets:
        assert path.is_file(), path
        if path.suffix == '.md':
            text = path.read_text('utf-8-sig')
            assert text.count('$$') % 2 == 0, path
            for target, local in local_links(text):
                resolved = (path.parent / local.replace('\\', '/')).resolve()
                assert resolved.exists() or (prospective and resolved == OUT), (path, target)
                links += 1
    assert all(f'## {n}.' in NOTE.read_text('utf8') for n in range(1, 11))
    assert '独立科学签审通过' in (HERE / 'review.md').read_text('utf8')
    return dict(round=1029, date='2026-10-08', all_delivery_checks_passed=True,
                scientific_result_reproduced=True, new_calibration_groups=1,
                cumulative_research_groups=3806, new_cognitive_axioms=0,
                calibration_kind='fixed-vertex lift-equivalence principal and transport certificates',
                bosonic_velocity_rank=41, Weyl_time_symbol_real_rank=4,
                PF_spatial_velocity_determinant='-16', PF_constraint_pivot_determinant='-1',
                PF_parameter_symbol_cases=3, null_parameter_symbol_included=True,
                off_shell_nonunique_lift_delta_O='-75/224',
                actual_trace_source='-8489/31360',
                recoil_E0_contribution_retained='212225/1404928',
                identically_conserved_nonzero_improvement=True,
                arbitrary_local_completeness_is_analytic_not_rank_certified=True,
                field_dependent_parameter_transport_is_analytic_not_finitely_classified=True,
                fixed_S0_S1_remain_inputs=True, all_first_order_vertices_classified=False,
                arbitrary_improvement_completion_proved=False,
                quantum_parent_realization_certified=False, goal_completed=False,
                live_navigation_frozen=False, neighboring_round_frozen=False,
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
