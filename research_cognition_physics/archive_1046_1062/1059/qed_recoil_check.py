"""1059: decimal transport of published inputs. Default never writes files."""
from pathlib import Path
from decimal import Decimal as D, localcontext
from statistics import NormalDist
import argparse
import json

HERE = Path(__file__).resolve().parent


def read_inputs():
    return json.loads((HERE / 'public_inputs.json').read_text(encoding='utf-8'))


def polynomial(q, a10, inp):
    x = 1 / (q * D(inp['pi']))
    terms, d_terms, mass_sd_terms = [], [], []
    for row in inp['qed_table_AKN2019']:
        n = row['n']
        if n == 5:
            c = a10 + D(row['mu'])
            u = D(row['u_mu'])
        else:
            c = sum(D(row[key]) for key in ('A1', 'mu', 'tau', 'mixed'))
            # Common mass-ratio errors across orders are never treated as independent.
            u = sum(D(row[key]) for key in ('u_mu', 'u_tau', 'u_mixed'))
        terms.append(c * x**n)
        d_terms.append(-D(n) * c * x**n / q)
        mass_sd_terms.append(u * x**n)
    value = sum(terms) + D(inp['fixed_nonQED']['hadronic']) + D(inp['fixed_nonQED']['weak'])
    return value, sum(d_terms), terms, sum(mass_sd_terms)


def pack(value):
    if isinstance(value, D):
        return format(value, '.25E')
    if isinstance(value, dict):
        return {key: pack(val) for key, val in value.items()}
    if isinstance(value, list):
        return [pack(val) for val in value]
    return value


def calculate(precision=65):
    inp = read_inputs()
    with localcontext() as ctx:
        ctx.prec = precision
        q = {name: D(v['value']) for name, v in inp['alpha_inverse'].items()}
        u = {name: D(v['u']) for name, v in inp['alpha_inverse'].items()}
        z = D(inp['statistics']['normal_quantile_for_two_marginals'])
        assert abs(float(z) - NormalDist().inv_cdf(0.9875)) < 2e-15
        intervals = {name: [q[name]-z*u[name], q[name]+z*u[name]] for name in q}
        gap = intervals['Rb2020'][0] - intervals['Cs2018'][1]
        assert gap > 0
        delta = q['Rb2020'] - q['Cs2018']
        assert delta == D('0.000000160')
        joint_independent_z = delta / (u['Rb2020']**2+u['Cs2018']**2).sqrt()
        joint_gaussian_worst_z = delta / (u['Rb2020']+u['Cs2018'])
        electron = D(inp['electron']['ae'])
        ue = D(inp['electron']['u'])
        versions = {}
        for label, version in inp['theory_versions'].items():
            a10, u10 = D(version['A1_10']), D(version['u_A1_10'])
            if 'setV' in version:
                assert a10 == D(version['setV']) + D(version['fermion_loop'])
                assert u10 == D(version['u_setV']) + D(version['u_fermion_loop'])
            predictions = {}
            for name in q:
                value, derivative, terms, u_mass = polynomial(q[name], a10, inp)
                assert derivative < 0
                x = 1 / (q[name] * D(inp['pi']))
                ua = abs(derivative) * u[name]
                components = [u10*x**5, D(inp['fixed_nonQED']['u_hadronic']), D(inp['fixed_nonQED']['u_weak']), u_mass]
                # Linear upper envelope valid for any covariance of finite-variance errors.
                u_theory_linear = sum(components)
                # Secondary conventional diagnostic; it makes additional independence assumptions.
                u_theory_rss = sum(v*v for v in components).sqrt()
                u_prediction_rss = (ua*ua+u_theory_rss*u_theory_rss).sqrt()
                residual = electron-value
                z_diagnostic = residual / (ue*ue+u_prediction_rss*u_prediction_rss).sqrt()
                # Explicit derivative and a symmetric finite difference agree.
                h = D('0.000001')
                numeric = (polynomial(q[name]+h,a10,inp)[0]-polynomial(q[name]-h,a10,inp)[0])/(2*h)
                assert abs(numeric/derivative-1) < D('1e-15')
                predictions[name] = {
                    'ae': value, 'ae_times_1e12': value*D('1e12'),
                    'qed_terms_times_1e12': [v*D('1e12') for v in terms],
                    'da_d_inverse_alpha': derivative,
                    'alpha_sd_times_1e12': ua*D('1e12'),
                    'known_mass_sd_linear_times_1e12': u_mass*D('1e12'),
                    'five_loop_sd_times_1e12': components[0]*D('1e12'),
                    'theory_sd_linear_upper_times_1e12': u_theory_linear*D('1e12'),
                    'theory_sd_rss_diagnostic_times_1e12': u_theory_rss*D('1e12'),
                    'prediction_sd_rss_diagnostic_times_1e12': u_prediction_rss*D('1e12'),
                    'exp_minus_prediction_times_1e12': residual*D('1e12'),
                    'independent_Gaussian_z_diagnostic_only': z_diagnostic,
                    'alpha_interval_only_prediction_times_1e12': [polynomial(intervals[name][1],a10,inp)[0]*D('1e12'),polynomial(intervals[name][0],a10,inp)[0]*D('1e12')]
                }
            versions[label] = predictions
        main = versions['V2024_main']
        assert main['Cs2018']['ae'] > main['Rb2020']['ae']
        # Exact closed rectangle criterion is not dependent on QED coefficients.
        assert delta > z*(u['Rb2020']+u['Cs2018'])
        assert abs(versions['AKN2019_legacy_check']['Cs2018']['ae_times_1e12']-D('1159652181.606')) < D('.0005')
        result = {
            'all_checks_passed': True, 'scope': 'fixed-version retrospective common-alpha test; mature QED adopted',
            'alpha_closure': {'z_marginal':z, 'intervals':intervals,'center_difference':delta,'positive_interval_gap':gap,
                'independent_Gaussian_z_diagnostic':joint_independent_z,'joint_Gaussian_arbitrary_correlation_min_abs_z':joint_gaussian_worst_z,
                'coverage_needs_marginal_model_not_independence':True,'marginal_normality_does_not_make_difference_Gaussian':True},
            'predictions': versions,
            'AH2025_minus_V2024_times_1e12':{name:(versions['AH2025_sensitivity'][name]['ae']-main[name]['ae'])*D('1e12') for name in q},
            'KNT2020_hadronic_central_shift_times_1e12':(D(inp['hadronic_version_diagnostic']['value'])-D(inp['fixed_nonQED']['hadronic']))*D('1e12'),
            'no_independent_combination_of_two_ae_residuals':True,
            'no_all_order_or_omitted_term_hard_bound':True,
            'mu_tau_coefficient_uncertainty_treated_conservatively':True
        }
        return pack(result)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-first', action='store_true')
    args=parser.parse_args()
    actual=calculate()
    assert actual==calculate(90), 'precision stability failed'
    path=HERE/'results.json'
    if args.write_first:
        with path.open('x',encoding='utf-8') as stream:
            json.dump(actual,stream,ensure_ascii=False,indent=2); stream.write('\n')
    else:
        assert actual==json.loads(path.read_text(encoding='utf-8'))
    print(json.dumps({'all_checks_passed':True,'mode':'exclusive_first_write' if args.write_first else 'read_only',
                     'positive_alpha_gap':actual['alpha_closure']['positive_interval_gap'],
                     'Rba_e':actual['predictions']['V2024_main']['Rb2020']['ae_times_1e12'],
                     'Csa_e':actual['predictions']['V2024_main']['Cs2018']['ae_times_1e12']}))
