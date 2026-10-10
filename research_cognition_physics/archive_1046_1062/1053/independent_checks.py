"""Independent 1053 review: exact bounds and analytic-angle Slater integrals.

No author module is imported. The integral uses closed angular averages and a
one-dimensional uniform Simpson rule, unlike the author's tensor Gaussian rule.
Integral agreement is a diagnostic, not a rigorous Coulomb-propagation bound.
Default is read-only; --write exclusively creates the review result.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / 'independent_checks_results.json'
RECEIPT = HERE / 'research_round_1053_checks.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integral():
    R, a, limit, intervals = 1000., 500., 32., 8192
    r = np.linspace(0., limit, intervals+1)
    C, D = r*r + R*R+a*a, 2*R*r
    local = 1/np.sqrt(r*r+a*a)
    # (1/2) integral_{-1}^{1} (C+D z)^(-1/2) dz,
    # rationalized to avoid subtracting close square roots.
    remote_mean = 2/(np.sqrt(C+D)+np.sqrt(C-D))
    remote_square = np.empty_like(r)
    remote_square[0] = 1/C[0]
    remote_square[1:] = np.arctanh(D[1:]/C[1:])/D[1:]
    mean_f = local-remote_mean
    mean_f2 = local*local-2*local*remote_mean+remote_square
    weights = np.ones(intervals+1)
    weights[1:-1:2] = 4
    weights[2:-1:2] = 2
    radial_density = 4*r*r*np.exp(-2*r)
    def integrate(value):
        return float((limit/intervals)/3 * np.dot(weights, radial_density*value))
    norm = integrate(np.ones_like(r))
    raw_f, raw_f2 = integrate(mean_f), integrate(mean_f2)
    # The stated positive bound controls the omitted radial probability only.
    # No Simpson truncation estimate is claimed.
    tail_probability = math.exp(-2*limit)*(2*limit**2+2*limit+1)
    assert tail_probability < 4e-25
    return dict(R=R, a=a, method='closed angular averages; 8192-panel Simpson on [0,32]',
                normalization=norm, raw_f=raw_f, raw_f2=raw_f2,
                coefficient_before_Lowdin=a*raw_f,
                omitted_radial_probability=tail_probability,
                certified_quadrature_error=False)


def build():
    receipt_bytes = RECEIPT.read_bytes()
    receipt = json.loads(receipt_bytes)
    assert set(receipt['owned_sha256']) == {
        '../research_note_1053.md', 'proof.md', 'coulomb_newton_source.py',
        'results.json', 'dependency_update.md', 'NEXT.md', 'review.md',
        'sources.md', 'verify_round1053.py'}
    assert len(receipt['historical_sha256']) == 8
    for name, digest in receipt['owned_sha256'].items():
        assert sha(HERE/name) == digest, name
    for name, digest in receipt['historical_sha256'].items():
        assert sha(ROOT/name) == digest, name
    author = json.loads((HERE/'results.json').read_text('utf8'))
    cert = author['rational_certificate']['exact']

    # Independently build a shorter positive Taylor lower bound for exp(20).
    term = F(1)
    lower_e20 = term
    for j in range(1, 20):
        term *= F(20, j)
        lower_e20 += term
    s = F(1, 10**6)
    assert lower_e20*s > 1+20+F(400, 3)
    # All square roots are bounded using rational squares with positive sides.
    assert F(2)/(1-s) < F(100, 49)
    assert (1+s) < 4*(1-s)**3
    assert 84 < F(3025, 36)
    assert F(1, 3) < F(676, 2025)
    assert F(1, 5) < F(179, 400)**2
    a, T, lb = F(500), F(1000, 137), F(22, 7000)
    initial = F(40, 7)/(3*a)+2*s
    natural = initial + T*(F(4, 3)/a*F(55, 6)*F(26, 45)+3/a**2)
    driven = natural + T*F(4, 3)/a*lb*F(26, 45)
    eps = F(22, 7000)*(F(20, 7)+T*F(37, 4))
    bad = eps*eps
    d_lo, delta_loose = F(221, 400), F(27, 250)
    joint = d_lo-delta_loose-bad
    recovery = (d_lo-delta_loose-2*bad)/2
    values = {
        'initial_source_error_strict_upper': initial,
        'natural_source_error_strict_upper': natural,
        'driven_source_error_strict_upper': driven,
        'old_instrument_error_strict_upper': eps,
        'bad_pointer_probability_strict_upper': bad,
        'joint_source_half_test_strict_lower': joint,
        'record_recovery_half_test_strict_lower': recovery,
        'joint_source_half_test_tighter_lower': d_lo-driven-bad,
        'record_recovery_half_test_tighter_lower': (d_lo-driven-2*bad)/2}
    for key, value in values.items():
        assert value == F(cert[key]), key
    assert driven < delta_loose
    assert 2*joint > F(791, 1000)
    assert recovery > F(1733, 10000)
    assert bad < F(49, 1000)
    assert (d_lo-delta_loose-2*F(49, 1000))/2 == F(693, 4000)

    measured = integral()
    expected = author['actual_R1000_orbital_integral']
    differences = {
        'raw_f': abs(measured['raw_f']-expected['raw_right_orbital_f']),
        'raw_f2': abs(measured['raw_f2']-expected['raw_right_orbital_f2']),
        'coefficient': abs(measured['coefficient_before_Lowdin']-
                           expected['coefficient_before_Lowdin'])}
    assert differences['coefficient'] < 1e-9
    assert differences['raw_f2'] < 1e-13
    assert abs(measured['normalization']-1) < 1e-9
    assert RECEIPT.read_bytes() == receipt_bytes
    return {
        'round': 1053, 'date': '2026-10-08', 'all_checks_passed': True,
        'author_functions_imported': False,
        'new_scientific_groups': 0, 'new_cognitive_axioms': 0,
        'author_receipt_sha256': hashlib.sha256(receipt_bytes).hexdigest(),
        'author_assets_checked': 9, 'historical_assets_checked': 8,
        'exact_values': {key: str(value) for key, value in values.items()},
        'independent_positive_Taylor_terms': 20,
        'independent_source_integral': measured,
        'author_integral_differences': differences,
        'task_bounds_depend_on_numerical_quadrature': False,
        'actual_Coulomb_propagation_solved_numerically': False,
        'whole_M4_or_roadmap_complete': False,
        'author_assets_unmodified': True}


def compare(x, y):
    if isinstance(x, dict):
        assert x.keys() == y.keys()
        for key in x:
            compare(x[key], y[key])
    elif isinstance(x, float):
        assert math.isclose(x, y, rel_tol=1e-11, abs_tol=2e-13)
    else:
        assert x == y


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    out = build()
    if args.write:
        with RESULT.open('x', encoding='utf8') as handle:
            json.dump(out, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    else:
        compare(out, json.loads(RESULT.read_text('utf8')))
    print(json.dumps({key: out[key] for key in [
        'round', 'all_checks_passed', 'new_scientific_groups',
        'author_integral_differences', 'author_receipt_sha256']}))


if __name__ == '__main__':
    main()
