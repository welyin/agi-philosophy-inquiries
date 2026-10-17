"""1276 author check: exact new counting/transport, no old ODE import or rerun.
Default read-only compare; --write exclusively creates the first results.json.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'research_cognition_physics').is_dir())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def fs(x):
    x = F(x)
    return f'{x.numerator}/{x.denominator}'

def binary(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log(p) - (1 - p) * math.log1p(-p)

def compute():
    baseline = json.loads((HERE / 'source_baseline.json').read_text(encoding='utf-8-sig'))
    bindings = baseline['files_sha256']
    for rel, expected in bindings.items():
        assert sha(ROOT / rel) == expected, ('source changed', rel)
    old_path = ROOT / 'research_cognition_physics/archive_1161_1200/1175/results.json'
    old = json.loads(old_path.read_text(encoding='utf-8-sig'))
    transport = old['group2_autonomous_common_transport']
    encl = old['group3_fixed_menu_and_error_certificate']['numerical_enclosure']
    resources = old['group1_parent_and_full_background_bridge']
    terminal = transport['background']['history'][-1]
    znum = F(str(terminal['z']))
    error = F(str(encl['Bloch_terminal_error_bound']))
    padding = F(1, 10**15)
    zlower, zupper = F(-29, 250), F(-14, 125)
    assert F(str(terminal['t'])) == 2
    assert zlower < znum - error - padding
    assert znum + error + padding < zupper
    assert resources['finite_K_budget_per_branch'] == 0.0001
    n = resources['K_defined_as_power_of_10']
    assert n == 52275242
    assert F(str(resources['resource_constants']['L_full_Fock_bound'])) > 7000000

    # The large K is never allocated. All final budget comparisons are rational.
    eta = F(1, 10000)
    rhalf_upper = eta / (2 * 7000000 * 2)
    assert rhalf_upper < F(1, 2**38)
    h_upper = F(39, 2**38)  # h2(x)<=x(1-log x), log2<1.
    assert h_upper < F(1, 10**8)
    assert n >= 12
    # (3j+1)/10^j decreases: 3(j+1)+1 < 10(3j+1), j>=1.
    assert 10 * (3 * 1 + 1) - (3 * (1 + 1) + 1) > 0
    assert 10 * 3 - 3 > 0
    type_upper = F(37, 10**12)
    assert type_upper < F(1, 10**8)
    b_sharp_rational = eta + h_upper + type_upper
    b_common = eta + F(2, 10**8)
    assert b_sharp_rational < b_common < F(101, 10**6)
    deficit = zupper**2 / 2
    assert deficit == F(98, 15625)
    net_drop = deficit - 2 * b_common
    assert net_drop == F(151799, 25000000)
    assert net_drop > F(3, 500)

    epsilon, nu, ell = F(3, 20), F(7, 20), F(16, 65)
    amin = (1 + epsilon * zupper) / (1 + nu * zupper)
    growth = amin**2 - 1
    assert amin == F(1229, 1201)
    assert growth == F(68040, 1442401) > F(47, 1000)
    # Independent direct enumeration on the original periodic graph.
    region = set(itertools.product(range(32), repeat=3))
    directions = [(a, s) for a in range(3) for s in [-1, 1]]
    boundary = set()
    cut_edges = 0
    for site in region:
        for axis, sign in directions:
            dst = list(site)
            dst[axis] = (dst[axis] + sign) % 65
            if tuple(dst) not in region:
                boundary.add(site)
                cut_edges += 1
    assert len(boundary) == 32**3 - 30**3 == 5768
    assert cut_edges == 6 * 32**2 == 6144
    area0 = cut_edges * ell**2
    assert area0 == F(1572864, 4225)
    assert area0 == 6 * (32 * ell)**2

    # Exact finite type theorem at all m, several moderate K.
    type_cases = 0
    min_scaled_type_probability = None
    for k in [1, 2, 3, 4, 7, 8, 16, 31, 32, 64]:
        for m in range(k + 1):
            prob = F(math.comb(k, m)) * F(m, k)**m * F(k - m, k)**(k - m)
            assert F(1, k + 1) <= prob <= 1
            scaled = (k + 1) * prob
            min_scaled_type_probability = scaled if min_scaled_type_probability is None else min(min_scaled_type_probability, scaled)
            type_cases += 1

    # Diagnostics use actual spectral probabilities, including entangled vectors.
    # They do not infer small variance for an arbitrary state or prove the K-large theorem.
    rng = np.random.default_rng(1276)
    nonproduct = []
    for k in [4, 6, 8]:
        dimension = 2**k
        mvalues = np.array([k - j.bit_count() for j in range(dimension)], dtype=float)
        zvalues = 2 * mvalues / k - 1
        lvalues = np.array([math.log(math.comb(k, int(m))) / k for m in mvalues])
        ghz = np.zeros(dimension, dtype=complex)
        ghz[0] = ghz[-1] = 1 / math.sqrt(2)
        middle = np.array([1.0 if int(m) == k // 2 else 0.0 for m in mvalues], dtype=complex)
        middle /= np.linalg.norm(middle)
        random = rng.normal(size=dimension) + 1j * rng.normal(size=dimension)
        random /= np.linalg.norm(random)
        for name, psi in [('GHZ', ghz), ('Dicke', middle), ('complex_random', random)]:
            probability = np.abs(psi)**2
            ztarget = float(probability @ zvalues)
            r = math.sqrt(max(0.0, float(probability @ ((zvalues - ztarget)**2))))
            lhs = abs(float(probability @ lvalues) - binary((1 + ztarget) / 2))
            rhs = math.log(k + 1) / k + binary(min(1.0, r) / 2)
            assert lhs <= rhs + 2e-14
            mat = psi.reshape(2, -1)
            red = mat @ mat.conj().T
            purity = float(np.trace(red @ red).real)
            assert purity < 0.999
            nonproduct.append({'K': k, 'state': name, 'z': round(ztarget, 15),
                               'r': round(r, 15), 'lhs': round(lhs, 15),
                               'upper': round(rhs, 15), 'one_qubit_purity': round(purity, 15)})
    assert len(nonproduct) == 9
    # A GHZ state has mean z=0 yet no binomial reference entropy; variance matters.
    assert nonproduct[0]['r'] == 1.0
    assert abs(nonproduct[0]['lhs'] - math.log(2)) < 2e-14

    binding_names = ['admission.md', 'source_baseline.json', 'proof.md', 'check.py']
    return {
        'round': 1276,
        'verdict': 'PASS_SCOPED_CONDITIONAL_REFERENCE_COUNT',
        'author_groups': 1,
        'source_bindings': {name: sha(HERE / name) for name in binding_names},
        'inherited_files_verified': len(bindings),
        'inherited_1175': {
            'T': 2, 'log10_K_integer': n,
            'z_saved_decimal': str(terminal['z']),
            'Bloch_terminal_error_bound': str(encl['Bloch_terminal_error_bound']),
            'extra_decimal_padding': fs(padding),
            'certified_z_enclosure': [fs(zlower), fs(zupper)],
            'old_numerical_assumptions_retained': True,
            'old_ODE_rerun': False,
            'exact_BG_preserved_in_comparator': True,
            'unknown_probe_R_uniform_on_original_domain': True,
        },
        'exact_new_bounds': {
            'eta_upper': fs(eta), 'r_half_upper_strict': fs(rhalf_upper),
            'h2_upper_strict': fs(h_upper), 'type_upper_strict': fs(type_upper),
            'per_endpoint_b_strict': fs(b_common),
            'binary_entropy_deficit': fs(deficit),
            'normalized_reference_net_drop_strict': fs(net_drop),
            'normalized_reference_net_drop_target_strict': fs(F(3, 500)),
            'amin_strict': fs(amin), 'relative_area_growth_strict': fs(growth),
            'relative_area_growth_target_strict': fs(F(47, 1000)),
            'N_lower': fs(F(17, 20)), 'N_upper': fs(F(23, 20)),
        },
        'boundary': {'grid': 65, 'region_side': 32, 'distinct_boundary_sites': len(boundary),
                     'cut_edges_counted_once': cut_edges, 'ell': fs(ell), 'A0': fs(area0)},
        'finite_type_exact_checks': {'cases': type_cases, 'min_scaled_probability': fs(min_scaled_type_probability)},
        'nonproduct_diagnostics': nonproduct,
        'scope': {
            'new_conditional_reference_choice': True,
            'actual_state_replaced_or_dephased': False,
            'actual_von_neumann_entropy_claimed': False,
            'full_Gibbs_reference_identity_claimed': False,
            'new_area_instrument_certified': False,
            'heat_temperature_or_Tkk_identity_certified': False,
            'strong_conjecture_refuted': False,
            'entire_time_interval_monotonicity_claimed': False,
            'all_counts_and_resources_remain_finite': True,
            'new_microscopic_ODE_or_large_K_simulation': False,
        },
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true', help='exclusive first result creation')
    args = parser.parse_args()
    actual = compute()
    output = HERE / 'results.json'
    if args.write:
        with output.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(actual, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write('\n')
        mode = 'WRITE_ONCE_PASS'
    else:
        expected = json.loads(output.read_text(encoding='utf-8'))
        assert expected == actual, 'READ_ONLY_MISMATCH'
        mode = 'READ_ONLY_MATCH'
    print(json.dumps({'status': mode, 'verdict': actual['verdict'],
                      'type_cases': actual['finite_type_exact_checks']['cases'],
                      'nonproduct_diagnostics': len(actual['nonproduct_diagnostics']),
                      'net_drop_strict': actual['exact_new_bounds']['normalized_reference_net_drop_strict'],
                      'area_growth_strict': actual['exact_new_bounds']['relative_area_growth_strict']}, ensure_ascii=False))

if __name__ == '__main__':
    main()
