"""1286: exact Gaussian-polynomial moments for an autonomous local-field record.
Only standard-library Fraction arithmetic. No truncated propagator or numerical ODE.
Default: read-only saved-result comparison. --write: exclusive first creation.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir() and (p/'scripts').is_dir())
RESEARCH = ROOT/'research_cognition_physics'
ONE = {(0, 0, 0): F(1)}
VARS = [{tuple(int(i == j) for i in range(3)): F(1)} for j in range(3)]
DIMS = (3, 2, 1)


def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, F(0)) + v
    return {k: v for k, v in out.items() if v}


def scale(p, c):
    return {k: v*c for k, v in p.items() if v*c}


def mul(p, q):
    out = {}
    for k, v in p.items():
        for l, w in q.items():
            z = tuple(a+b for a, b in zip(k, l))
            out[z] = out.get(z, F(0)) + v*w
    return {k: v for k, v in out.items() if v}


def derivative(p, axis):
    out = {}
    for k, v in p.items():
        if k[axis]:
            z = list(k)
            z[axis] -= 1
            out[tuple(z)] = k[axis]*v
    return out


def kinetic(p, axis):
    v, d = VARS[axis], DIMS[axis]
    dp = derivative(p, axis)
    return add(scale(mul(v, derivative(dp, axis)), -2),
               mul(add(v, scale(ONE, -d)), dp),
               mul(add(scale(ONE, F(d, 4)), scale(v, F(-1, 8))), p))


def potential_components():
    x, y, z = VARS
    return {'gauge_scalar_interaction': scale(mul(x, y), F(1, 2)),
            'scalar_quartic': scale(mul(y, y), F(119, 1248)),
            'scalar_record_interaction': scale(mul(y, z), F(1, 2))}


def act_h(p):
    return add(*(kinetic(p, i) for i in range(3)),
               *(mul(v, p) for v in potential_components().values()))


@lru_cache(None)
def monomial_moment(k):
    return math.prod(math.prod(d+2*j for j in range(n)) for n, d in zip(k, DIMS))


def expectation(p):
    return sum((v*monomial_moment(k) for k, v in p.items()), F(0))


def inner(p, q):
    return expectation(mul(p, q))


def laguerre(n):
    return {(0, j, 0): F((-1)**j*math.comb(n, j), 2**j*math.factorial(j))
            for j in range(n+1)}


def ceil_root(x):
    assert x >= 0
    n = math.isqrt(x.numerator//x.denominator)
    n += int(n*n != x)
    assert (n == 0 or (n-1)**2 < x) and x <= n*n
    return n


def fs(x):
    x = F(x)
    return f'{x.numerator}/{x.denominator}'


def matrix_strings(m):
    return [[fs(x) for x in row] for row in m]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def calculate():
    p0 = laguerre(3)
    p1 = scale(add(scale(laguerre(0), 4), scale(laguerre(6), 3)), F(1, 5))
    code = [p0, p1]
    gram = [[inner(p, q) for q in code] for p in code]
    hg = [[inner(p, act_h(q)) for q in code] for p in code]
    yg = [[inner(p, mul(VARS[1], q)) for q in code] for p in code]
    energy = F(4561, 78)
    assert gram == [[1, 0], [0, 1]]
    assert hg == [[energy, 0], [0, energy]]
    assert yg == [[14, 0], [0, F(266, 25)]]

    # Check the Laguerre normalization and the stated diagonal energy polynomial
    # from the actual differential Hamiltonian, without importing earlier code.
    basis = [laguerre(n) for n in range(9)]
    for i, p in enumerate(basis):
        en = F(357*i*i + 1683*i + 860, 156)
        assert inner(p, act_h(p)) == en
        assert inner(p, mul(VARS[1], p)) == 4*i+2
        for j, q in enumerate(basis):
            assert inner(p, q) == int(i == j)
            if abs(i-j) > 2:
                assert inner(p, act_h(q)) == 0
            if abs(i-j) > 1:
                assert inner(p, mul(VARS[1], q)) == 0

    ledgers = []
    component_gram = {}
    for name, axis in [('gauge_kinetic', 0), ('scalar_kinetic', 1), ('record_kinetic', 2)]:
        component_gram[name] = [[inner(p, kinetic(q, axis)) for q in code] for p in code]
    for name, v in potential_components().items():
        component_gram[name] = [[inner(p, mul(v, q)) for q in code] for p in code]
    for j in range(2):
        row = {name: m[j][j] for name, m in component_gram.items()}
        row['record_quartic'] = F(0)
        assert sum(row.values(), F(0)) == energy
        assert all(v >= 0 for v in row.values())
        ledgers.append({name: fs(v) for name, v in row.items()})
    assert all(m[0][1] == m[1][0] == 0 for m in component_gram.values())

    moment_results = []
    m4 = []
    expected_ceil = [[1, 79, 12463, 2925311, 968191741],
                     [1, 124, 44613, 20742456, 12107230522]]
    for j, p in enumerate(code):
        powers = [p]
        for _ in range(4):
            powers.append(act_h(powers[-1]))
        norm2 = [inner(v, v) for v in powers]
        ceils = [ceil_root(x) for x in norm2]
        assert ceils == expected_ceil[j]
        bound = sum(math.comb(4, k)*ceils[k]*ceils[4-k] for k in range(5))
        m4.append(bound)
        # Self-adjointness crosscheck at the polynomial/Gaussian level.
        assert inner(powers[0], powers[2]) == norm2[1]
        assert inner(powers[1], powers[3]) == norm2[2]
        moment_results.append({'code': j, 'norm_squared': [fs(x) for x in norm2],
                               'norm_integer_upper': ceils,
                               'polynomial_term_counts': [len(v) for v in powers],
                               'fourth_probability_derivative_upper': bound})
    assert m4 == [4717138248, 56732896010]
    total = sum(m4)
    assert total == 61450034258

    # Characteristic function of N(0,1) gives <r sin r>=exp(-1/2).
    # Prove the conservative exp(-1/2)>1/2 using an exact series upper for e.
    e_upper = sum((F(1, math.factorial(n)) for n in range(9)), F(0)) + F(10, 9*math.factorial(9))
    assert e_upper < 3 < 4
    y_difference = yg[0][0]-yg[1][1]
    second_order_prefactor = y_difference/4
    assert y_difference == F(84, 25) and second_order_prefactor == F(21, 25)
    t = F(1, 10**6)
    point_lower = second_order_prefactor*t*t/2 - F(total, 24)*t**4
    window_coefficient = second_order_prefactor/2 - F(total, 24)*(2*t)**2
    uniform_lower = F(2, 5)*t*t
    assert window_coefficient > F(2, 5)
    assert point_lower > uniform_lower == F(1, 2500000000000)
    assert total*(2*t)**2 < F(12, 25)*1  # Equivalent safe fourth-order margin.

    volume, mp2 = F(1), F(11, 2)
    slope2 = energy/(3*volume*mp2)
    pressure_coefficient = energy/3
    assert slope2 == F(4561, 1287)
    friedmann = 3*mp2*volume*slope2-energy
    acceleration = -slope2+(energy+3*pressure_coefficient)/(6*volume*mp2)
    continuity = -4*energy+3*(energy+pressure_coefficient)
    assert friedmann == acceleration == continuity == 0
    # No implicit adjustment of source in the coupling-deletion comparison.
    deleted_energies = [energy-F(7), energy-F(133, 25)]
    assert deleted_energies[0] != deleted_energies[1]

    sources = [HERE/'proof.md', HERE/'check.py', HERE.parent/'1285/proof.md',
               HERE.parent/'1285/results.json', HERE.parent/'_shared/common_process_boundary_after1285.md',
               RESEARCH/'archive_956_989/research_note_959.md',
               RESEARCH/'archive_1159_1160/1160/proof.md',
               RESEARCH/'archive_585_628/research_note_625.md']
    return {
        'round': 1286, 'status': 'PASS_WITH_SCOPE', 'author_groups': 1, 'new_axioms': 0,
        'method': 'exact Fraction Gaussian radial differential recurrence; no Hamiltonian truncation or propagation quadrature',
        'source_bindings': {p.relative_to(ROOT).as_posix(): sha(p) for p in sources},
        'parameters': {'g': '1/1', 'kappa': '1/1', 'lambda': '119/312', 'lambda_record': '0/1',
                       'coordinate_volume': fs(volume), 'M_P_squared': fs(mp2), 'initial_scale': '1/1',
                       'Gaussian_coordinate_variance': '1/1', 'T_conformal': fs(t)},
        'normal_code': {'polynomials': ['L_3(Y/2)', '(4 L_0(Y/2)+3 L_6(Y/2))/5'],
                        'Gram': matrix_strings(gram), 'compressed_H': matrix_strings(hg),
                        'compressed_Y': matrix_strings(yg), 'Laguerre_basis_exact_checks': 9},
        'energy_ledger_by_basis_code': ledgers,
        'energy_component_grams': {name: matrix_strings(m) for name, m in component_gram.items()},
        'exact_H_norm_certificates': moment_results,
        'finite_record_certificate': {'effect': '(I+cos r)/2', 'difference_order': 'p_0-p_1',
            'second_order_coefficient': '21/25 * exp(-1/2)', 'e_exact_series_upper': fs(e_upper),
            'fourth_derivative_sum_upper': total, 'point_T_strict_lower': fs(point_lower),
            'point_T_lower_decimal_diagnostic': float(point_lower),
            'window': [fs(t), fs(2*t)], 'window_coefficient_lower': fs(window_coefficient),
            'uniform_strict_gap': fs(uniform_lower),
            'deleted_kappa_exact_record_independence': True,
            'deleted_kappa_initial_energies': [fs(x) for x in deleted_energies]},
        'same_source_geometry': {'energy_all_unknown_inputs': fs(energy), 'slope_squared': fs(slope2),
            'scale': 'a(eta)=1+sqrt(4561/1287)*eta',
            'proper_time': 'tau(eta)=eta+sqrt(4561/1287)*eta^2/2',
            'density': '(4561/78)/a^4', 'pressure': '(4561/234)/a^4',
            'Friedmann_residual_coefficient': fs(friedmann),
            'acceleration_residual_coefficient': fs(acceleration),
            'continuity_residual_coefficient': fs(continuity)},
        'scope': {'local_classical_neutral_field_and_interaction_declared': True,
            'normal_orthogonal_equal_mean_energy_code': True,
            'arbitrary_unknown_input_and_reference_exact_joint_isometry': True,
            'common_nonselective_homogeneous_mean_source_for_code': True,
            'quantum_record_marginal_distinguishability': True,
            'complete_record_and_interaction_energy_included': True,
            'autonomous_no_external_switch': True,
            'source_covers_unread_process': True,
            'finite_probability_remainder_rigorous': True,
            'declared_homogeneous_quantum_effective_class': True,
            'full_4D_quantum_invariant_sector_proved': False,
            'actual_terminal_effect_instrument_or_cost_proved': False,
            'strong_copy_or_stable_classical_log': False,
            'conditional_branch_semiclassical_geometry': False,
            'thermal_reference_identity': False,
            'new_axiom_adopted': False,
            'stage_completion_claimed': False}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write:
        with path.open('x', encoding='utf-8', newline='\n') as out:
            json.dump(result, out, ensure_ascii=False, indent=2)
            out.write('\n')
        print('PASS_WITH_SCOPE: EXCLUSIVE_FIRST_WRITE')
    else:
        assert result == json.loads(path.read_text(encoding='utf-8')), 'Saved result mismatch'
        print('PASS_WITH_SCOPE: READ_ONLY_MATCH')
    print('uniform record gap > '+result['finite_record_certificate']['uniform_strict_gap'])
    print('common energy = '+result['same_source_geometry']['energy_all_unknown_inputs'])
    print('results_sha256 = '+sha(path))


if __name__ == '__main__':
    main()
