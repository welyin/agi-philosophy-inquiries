"""Independent 1059 arithmetic: rational Horner, no author module imports.

Default is read-only. --write-first exclusively creates the review JSON.
The script checks published-input transport, not QED integrals or full errors.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D, localcontext
from statistics import NormalDist
import argparse
import hashlib
import json

BASE = Path(__file__).resolve().parent
FROZEN_RECEIPT_SHA = '81e1963c5bb4a6f9205ef34318d7a26466a9b935d2acd874322701aa24c89c2e'


def dec(x):
    if isinstance(x, F):
        return D(x.numerator) / D(x.denominator)
    return D(x)


def evaluate(q, coefficients, pi):
    x = 1 / (q * pi)
    # Independent exact-rational Horner evaluation of p(x) and p'(x).
    p = F(0)
    dp = F(0)
    for c in reversed([F(0)] + coefficients):
        dp, p = dp*x+p, p*x+c
    return p, -dp*x/q, x


def check():
    receipt_path=BASE/'research_round_1059_checks.json'
    receipt_sha=hashlib.sha256(receipt_path.read_bytes()).hexdigest()
    assert receipt_sha==FROZEN_RECEIPT_SHA
    receipt=json.loads(receipt_path.read_text(encoding='utf-8'))
    own={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()
         for name in receipt['owned_sha256']}
    history={name:hashlib.sha256((BASE.parents[1]/name).read_bytes()).hexdigest()
             for name in receipt['historical_sha256']}
    assert own==receipt['owned_sha256'] and len(own)==9
    assert history==receipt['historical_sha256'] and len(history)==5
    assert receipt['independent_review_complete_at_author_freeze'] is False
    inp = json.loads((BASE/'public_inputs.json').read_text(encoding='utf-8'))
    out = json.loads((BASE/'results.json').read_text(encoding='utf-8'))
    tested = 0
    largest_rounding_ulp = D(0)

    def compare(actual, expected):
        nonlocal tested, largest_rounding_ulp
        target = dec(expected)
        saved = D(actual)
        ulp = D(10) ** saved.adjusted() * D('1e-25') if saved else D('1e-90')
        fraction = abs(saved-target)/ulp
        assert fraction <= D('.501'), (actual, str(target), str(fraction))
        largest_rounding_ulp = max(largest_rounding_ulp, fraction)
        tested += 1

    with localcontext() as ctx:
        ctx.prec = 85
        q = {k:F(v['value']) for k,v in inp['alpha_inverse'].items()}
        u = {k:F(v['u']) for k,v in inp['alpha_inverse'].items()}
        z = F(inp['statistics']['normal_quantile_for_two_marginals'])
        assert abs(float(z)-NormalDist().inv_cdf(.9875)) < 2e-15
        bounds = {k:(q[k]-z*u[k],q[k]+z*u[k]) for k in q}
        gap = bounds['Rb2020'][0]-bounds['Cs2018'][1]
        difference = q['Rb2020']-q['Cs2018']
        assert difference == F(160,10**9) and gap > 0
        closure = out['alpha_closure']
        compare(closure['z_marginal'],z)
        compare(closure['center_difference'],difference)
        compare(closure['positive_interval_gap'],gap)
        for k, ends in bounds.items():
            for actual, expected in zip(closure['intervals'][k],ends):
                compare(actual,expected)
        compare(closure['independent_Gaussian_z_diagnostic'],
                dec(difference)/dec(sum(s*s for s in u.values())).sqrt())
        compare(closure['joint_Gaussian_arbitrary_correlation_min_abs_z'],
                difference/sum(u.values()))
        assert F(1)-2*F(1,40)==F(19,20)

        pi = F(inp['pi'])
        fixed = inp['fixed_nonQED']
        constant = F(fixed['hadronic'])+F(fixed['weak'])
        ae = F(inp['electron']['ae'])
        ue = F(inp['electron']['u'])
        summaries = {}
        central = {}
        for label, version in inp['theory_versions'].items():
            coeff, mass_uncertainty = [], []
            for row in inp['qed_table_AKN2019']:
                if row['n']==5:
                    coeff.append(F(version['A1_10'])+F(row['mu']))
                    mass_uncertainty.append(F(row['u_mu']))
                else:
                    coeff.append(sum(F(row[k]) for k in ('A1','mu','tau','mixed')))
                    mass_uncertainty.append(sum(F(row[k]) for k in ('u_mu','u_tau','u_mixed')))
            if 'setV' in version:
                assert F(version['A1_10'])==F(version['setV'])+F(version['fermion_loop'])
                assert F(version['u_A1_10'])==F(version['u_setV'])+F(version['u_fermion_loop'])
            for name, inverse_alpha in q.items():
                value, derivative, x = evaluate(inverse_alpha,coeff,pi)
                value += constant
                assert derivative < 0
                central[label,name] = value
                ua = abs(derivative)*u[name]
                um = sum(s*x**n for n,s in enumerate(mass_uncertainty,1))
                known = [F(version['u_A1_10'])*x**5,
                         F(fixed['u_hadronic']),F(fixed['u_weak']),um]
                linear = sum(known)
                var_theory = sum(s*s for s in known)
                var_prediction = ua*ua+var_theory
                residual = ae-value
                data = out['predictions'][label][name]
                expected = {
                    'ae': value,
                    'ae_times_1e12':value*10**12,
                    'da_d_inverse_alpha':derivative,
                    'alpha_sd_times_1e12':ua*10**12,
                    'known_mass_sd_linear_times_1e12':um*10**12,
                    'five_loop_sd_times_1e12':known[0]*10**12,
                    'theory_sd_linear_upper_times_1e12':linear*10**12,
                    'theory_sd_rss_diagnostic_times_1e12':dec(var_theory).sqrt()*10**12,
                    'prediction_sd_rss_diagnostic_times_1e12':dec(var_prediction).sqrt()*10**12,
                    'exp_minus_prediction_times_1e12':residual*10**12,
                    'independent_Gaussian_z_diagnostic_only':dec(residual)/dec(ue*ue+var_prediction).sqrt()
                }
                for key,v in expected.items():
                    compare(data[key],v)
                for n,c in enumerate(coeff,1):
                    compare(data['qed_terms_times_1e12'][n-1],c*x**n*10**12)
                for pos,bound in enumerate(reversed(bounds[name])):
                    point = evaluate(bound,coeff,pi)[0]+constant
                    compare(data['alpha_interval_only_prediction_times_1e12'][pos],point*10**12)
                summaries[label+'/'+name]={
                    'ae_times_1e12':format(dec(value)*10**12,'.18f'),
                    'exp_minus_prediction_times_1e12':format(dec(residual)*10**12,'.18f'),
                    'independent_z_diagnostic':format(expected['independent_Gaussian_z_diagnostic_only'],'.15f')}
        for name in q:
            shift = central['AH2025_sensitivity',name]-central['V2024_main',name]
            compare(out['AH2025_minus_V2024_times_1e12'][name],shift*10**12)
        compare(out['KNT2020_hadronic_central_shift_times_1e12'],
                (F(inp['hadronic_version_diagnostic']['value'])-F(fixed['hadronic']))*10**12)
        # Exact cancellation of the common electron, hadronic and weak rows
        # in residual_Rb - residual_Cs; it is not two independent deviations.
        for label in inp['theory_versions']:
            rrb=ae-central[label,'Rb2020']
            rcs=ae-central[label,'Cs2018']
            assert rrb-rcs==central[label,'Cs2018']-central[label,'Rb2020']
        return {
            'all_computational_checks_passed':True,
            'method':'exact rational Horner plus 85-digit square roots; no author code import',
            'numeric_saved_values_checked':tested,
            'maximum_difference_in_last_printed_ulp':format(largest_rounding_ulp,'.12f'),
            'input_asset_hashes':{name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()
                                 for name in ('public_inputs.json','qed_recoil_check.py','results.json')},
            'author_receipt_sha256':receipt_sha,
            'all_author_owned_sha256_verified':own,
            'all_historical_sha256_verified':history,
            'alpha_gap_times_1e9':format(dec(gap)*10**9,'.16f'),
            'bonferroni_nominal_joint_coverage':'at least 0.95 under adopted marginal coverage',
            'unknown_correlation_joint_gaussian_z_lower':str(dec(difference/sum(u.values()))),
            'summaries':summaries,
            'shared_residual_cancellation_checked':True,
            'round_status':'arithmetic review only; scientific scope and frozen hashes signed separately',
            'no_integrals_raw_data_fit_or_all_order_bound':True}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-first',action='store_true')
    args=parser.parse_args()
    result=check()
    target=BASE/'statistical_review_checks.json'
    if args.write_first:
        with target.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
    elif target.exists():
        assert result==json.loads(target.read_text(encoding='utf-8'))
    print(json.dumps({'pass':True,'compared':result['numeric_saved_values_checked'],
                      'max_last_ulp':result['maximum_difference_in_last_printed_ulp'],
                      'mode':'exclusive_first_write' if args.write_first else 'read_only'}))
