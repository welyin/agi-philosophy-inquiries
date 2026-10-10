"""Independent review arithmetic for 1052; no author module is imported.

Decimal input arithmetic/Newton optimization and erfc inversion are separate
from the author's numpy polynomial root and NormalDist implementations.
Default mode only reads; --write exclusively creates this review's result.
This adds no empirical sample, scientific group, or coverage guarantee.
"""
from pathlib import Path
from decimal import Decimal, localcontext
from fractions import Fraction
import argparse
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / 'independent_checks_results.json'
RECEIPT = HERE / 'research_round_1052_checks.json'
D = Decimal


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    receipt_bytes = RECEIPT.read_bytes()
    receipt = json.loads(receipt_bytes)
    assert len(receipt['owned_sha256']) == 9
    assert len(receipt['historical_sha256']) == 5
    for name, digest in receipt['owned_sha256'].items():
        assert sha(HERE / name) == digest, name
    for name, digest in receipt['historical_sha256'].items():
        assert sha(ROOT / name) == digest, name
    data = json.loads((HERE / 'weak_public_inputs.json').read_text('utf8'),
                      parse_float=D)
    transcription = {
        'belle_ii': {
            'R_tau': D('.9675'), 'R_tau_stat': D('.0007'),
            'R_tau_syst': D('.0036'),
            'published_SM_R_tau_rounded': D('.9726'),
            'published_effective_amplitude_ratio': D('.9974'),
            'published_effective_amplitude_ratio_sigma': D('.0019')},
        'atlas': {
            'R_WZ': D('.9990'), 'R_WZ_stat': D('.0022'),
            'R_WZ_syst': D('.0036'), 'external_R_Z': D('1.0009'),
            'external_R_Z_sigma': D('.0028'), 'published_R_W': D('.9995'),
            'published_R_W_stat': D('.0022'),
            'published_R_W_syst': D('.0036'),
            'published_R_W_external': D('.0014'),
            'published_R_W_sigma_rounded': D('.0045')}}
    for block, fields in transcription.items():
        for key, value in fields.items():
            assert data[block][key] == value, (block, key)

    # Solve the two-sided tail probability directly, not with NormalDist.
    lo, hi = 0.0, 5.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > .025:
            lo = mid
        else:
            hi = mid
    q_float = (lo + hi) / 2
    assert abs(math.erfc(q_float / math.sqrt(2)) - .025) < 2e-17
    union_coverage = 1 - 2 * Fraction(1, 40)
    assert union_coverage == Fraction(19, 20)

    with localcontext() as context:
        context.prec = 70
        x, sx, y = D('.9974'), D('.0019'), D('.9995')
        components = [D('.0022'), D('.0036'), D('.0014')]
        sy2 = sum(v * v for v in components)
        sy = sy2.sqrt()
        sx2 = sx * sx
        q = D(str(q_float))
        xi, yi = [x-q*sx, x+q*sx], [y-q*sy, y+q*sy]
        corners = [xx / yy.sqrt() - 1 for xx in xi for yy in yi]
        eta_interval = [min(corners), max(corners)]
        common_interval = [max(xi[0], yi[0].sqrt()),
                           min(xi[1], yi[1].sqrt())]
        assert xi[0] < 1 < xi[1] and yi[0] < 1 < yi[1]
        assert common_interval[0] < 1 < common_interval[1]
        assert eta_interval[0] < 0 < eta_interval[1]

        # Strict convexity on the entire a>=0 domain proves uniqueness here.
        # chi''(a) = 2/sx^2 + (12 a^2 - 4y)/sy^2.
        convexity_floor = 2/sx2 - 4*y/sy2
        assert convexity_floor > 0
        def chi(a):
            return (a-x)**2/sx2 + (a*a-y)**2/sy2
        def grad(a):
            return 2*(a-x)/sx2 + 4*a*(a*a-y)/sy2
        a = D(1)
        for _ in range(16):
            a -= grad(a)/(2/sx2 + (12*a*a-4*y)/sy2)
        assert a > 0 and abs(grad(a)) < D('1e-60')
        fixed_chi, fit_chi = chi(D(1)), chi(a)
        eta = x/y.sqrt()-1
        eta_sigma = (sx2/y + x*x*sy2/(4*y**3)).sqrt()
        sat_a = y.sqrt()
        assert abs(sat_a*(1+eta)-x) < D('1e-65')
        assert abs(sat_a**2-y) < D('1e-65')

        wz, z, sz = D('.9990'), D('1.0009'), D('.0028')
        rw = wz*z.sqrt()
        ext = wz*sz/(2*z.sqrt())
        tau_conversion = (D('.9675')/D('.9726')).sqrt()
        # Allow the last reported decimal in the intermediate quantities.
        half_unit = D('.00005')
        rw_rounding_range = [(wz-half_unit)*(z-half_unit).sqrt(),
                             (wz+half_unit)*(z+half_unit).sqrt()]
        assert rw_rounding_range[0] < y < rw_rounding_range[1]
        sigma_rounding_range = [sum((c+s*half_unit)**2 for c in components).sqrt()
                                for s in [-1, 1]]
        assert max(sigma_rounding_range[0], D('.00445')) < min(
            sigma_rounding_range[1], D('.00455'))
        z_family = [(yy/wz)**2 for yy in [D('.98'), D(1), D('1.02')]]
        assert all(abs(wz*zz.sqrt()-yy) < D('1e-65') for yy, zz in zip(
            [D('.98'), D(1), D('1.02')], z_family))

        independent = {
            'source_conversions': {
                'R_W_from_rounded_WZ_and_Z': float(rw),
                'external_Z_error_propagated': float(ext),
                'W_component_quadrature_sigma': float(sy),
                'tau_amplitude_from_rounded_SM_normalization': float(tau_conversion)},
            'fixed_adopted_SM_branch': {
                'tau_marginal_standardized_residual': float((x-1)/sx),
                'W_marginal_standardized_residual': float((y-1)/sy)},
            'primary_marginal_gaussian_bonferroni': {
                'z_per_two_sided_interval': q_float,
                'tau_amplitude_interval': list(map(float, xi)),
                'W_ratio_interval': list(map(float, yi)),
                'common_amplitude_allowed_interval': list(map(float, common_interval)),
                'relative_transport_residual_point': float(eta),
                'relative_transport_residual_interval': list(map(float, eta_interval))},
            'secondary_independent_gaussian_diagnostics': {
                'fixed_SM_chi_square_two_observables': float(fixed_chi),
                'fixed_SM_nominal_p_two_dof': math.exp(-float(fixed_chi)/2),
                'common_amplitude_best_fit': float(a),
                'common_amplitude_min_chi_square': float(fit_chi),
                'common_amplitude_nominal_p_one_dof': math.erfc(math.sqrt(float(fit_chi)/2)),
                'transport_delta_method_sigma': float(eta_sigma),
                'transport_delta_method_standardized_residual': float(eta/eta_sigma)},
            'freedom_diagnostics': {
                'saturated_amplitude': float(sat_a),
                'saturated_transport_residual': float(eta),
                'saturated_model_remaining_goodness_of_fit_dof': 0}}
        extra = {
            'chi_second_derivative_global_lower_bound': float(convexity_floor),
            'decimal_newton_derivative_abs': str(abs(grad(a))),
            'W_ratio_intermediate_rounding_range': list(map(float, rw_rounding_range)),
            'component_sigma_rounding_range': list(map(float, sigma_rounding_range)),
            'conditional_union_coverage_exact': str(union_coverage),
            'external_Z_family': list(map(float, z_family))}

    author = json.loads((HERE / 'weak_two_scale_results.json').read_text('utf8'))
    differences = []
    def compare(got, reference, label):
        if isinstance(got, dict):
            for key, value in got.items():
                compare(value, reference[key], label+'.'+key)
        elif isinstance(got, list):
            assert len(got) == len(reference), label
            for i, (value, ref) in enumerate(zip(got, reference)):
                compare(value, ref, label+f'[{i}]')
        else:
            diff = abs(got-reference)
            assert diff < 2e-12, (label, got, reference)
            differences.append(diff)
    compare(independent, author, 'result')
    assert RECEIPT.read_bytes() == receipt_bytes
    return {
        'round': 1052, 'date': '2026-10-08', 'all_checks_passed': True,
        'method': 'Independent Decimal(70) Newton and erfc inversion; no author imports',
        'new_scientific_groups': 0, 'new_empirical_groups': 0,
        'source_data_transcription_checked': True,
        'author_assets_checked': 9, 'historical_assets_checked': 5,
        'author_receipt_sha256': hashlib.sha256(receipt_bytes).hexdigest(),
        'numerical_comparisons': len(differences),
        'maximum_author_difference': max(differences),
        'independent_values': independent, 'additional_checks': extra,
        'coverage_assumption_tested_as_fact': False,
        'author_freeze_unmodified': True, 'roadmap_complete': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    out = build()
    if args.write:
        with RESULT.open('x', encoding='utf8') as f:
            json.dump(out, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        assert out == json.loads(RESULT.read_text('utf8'))
    print(json.dumps({k: out[k] for k in ['round', 'all_checks_passed',
        'numerical_comparisons', 'maximum_author_difference',
        'author_receipt_sha256', 'new_scientific_groups']}))


if __name__ == '__main__':
    main()
