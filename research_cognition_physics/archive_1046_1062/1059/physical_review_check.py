"""Independent 1059 physical-numerical review. No author imports or mutations.

Primary numbers are transcribed independently. Exact Fraction arithmetic is used
for finite printed-input transport; square roots/normal diagnostics use Decimal
and math respectively. Default compares the saved review result, if present.
--record is an exclusive first write of this reviewer's result only.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from pathlib import Path
from hashlib import sha256
from math import erfc, sqrt
import argparse
import json

HERE = Path(__file__).resolve().parent
PI = F('3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679')
CAL = {'Cs2018': (F('137.035999046'), F('27e-9')),
       'Rb2020': (F('137.035999206'), F('11e-9'))}
# Each row: universal, muon, tau, mixed; from AKN2019 Table 1.
ROWS = [('0.5','0','0','0'),
        ('-.328478965579193','.519738676e-6','.183790e-8','0'),
        ('1.181241456587','-.737394164e-5','-.658273e-7','.1909e-12'),
        ('-1.912245764','.916197070e-3','.74292e-5','.74687e-6')]
MASS_U = [('0','0','0'), ('.000000024e-6','.000025e-8','0'),
          ('.000000024e-5','.000079e-7','.0001e-12'),
          ('.000000037e-3','.00012e-5','.00028e-6'),
          ('.00039','0','0')]
VERSIONS = {'V2024_main':(F('5.891'),F('.061')),
            'AH2025_sensitivity':(F('6.800')-F('.93042'), F('.128')+F('.00361')),
            'AKN2019_legacy_check':(F('6.737'),F('.159'))}
AE, UE = F('.00115965218059'), F('1.3e-13')
HAD, WH = F('1.693e-12'), F('.03053e-12')
UH, UW = F('.012e-12'), F('.00023e-12')
Z = F('2.241402727604947')


def dec(x):
    return Decimal(x.numerator)/Decimal(x.denominator) if isinstance(x,F) else Decimal(x)


def string(x):
    return format(dec(x), '.19E')


def evaluate(q, a10):
    t = 1/(PI*q)
    coeff = [sum(map(F,row)) for row in ROWS] + [a10-F('.00382')]
    # Horner, separately from author's per-order summation.
    p = F(0)
    for c in reversed(coeff):
        p = (p+c)*t
    value = p+HAD+WH
    derivative = -sum((n+1)*c*t**(n+1) for n,c in enumerate(coeff))/q
    terms = [c*t**(n+1) for n,c in enumerate(coeff)]
    umass = sum(sum(map(F,row))*t**(n+1) for n,row in enumerate(MASS_U))
    return value, derivative, terms, umass


def assert_close(actual, expected, absolute='1e-24'):
    assert abs(dec(actual)-Decimal(expected)) <= Decimal(absolute), (str(actual), expected)


def run():
    with localcontext() as ctx:
        ctx.prec = 70
        author = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
        inputs = json.loads((HERE/'public_inputs.json').read_text(encoding='utf-8'))
        for key,(q,uq) in CAL.items():
            assert F(inputs['alpha_inverse'][key]['value']) == q
            assert F(inputs['alpha_inverse'][key]['u']) == uq
        for row,expected in zip(inputs['qed_table_AKN2019'][:4],ROWS):
            assert [F(row[k]) for k in ('A1','mu','tau','mixed')] == list(map(F,expected))
        assert F(inputs['electron']['ae']) == AE and F(inputs['electron']['u']) == UE
        assert F(inputs['fixed_nonQED']['hadronic']) == HAD
        assert F(inputs['fixed_nonQED']['weak']) == WH
        assert F(inputs['fixed_nonQED']['u_hadronic']) == UH
        assert F(inputs['fixed_nonQED']['u_weak']) == UW
        assert F(inputs['hadronic_version_diagnostic']['value']) == F('1.7030e-12')
        assert F(inputs['hadronic_version_diagnostic']['u']) == F('.0077e-12')
        for name,(a,u) in VERSIONS.items():
            assert F(inputs['theory_versions'][name]['A1_10']) == a
            assert F(inputs['theory_versions'][name]['u_A1_10']) == u
        for row,expected in zip(inputs['qed_table_AKN2019'],MASS_U):
            assert F(row['u_mu']) == F(expected[0])
            if row['n'] < 5:
                assert F(row['u_tau']) == F(expected[1])
                assert F(row['u_mixed']) == F(expected[2])
        assert 'omitted' in inputs['qed_table_AKN2019'][4]['tau']
        assert 'omitted' in inputs['qed_table_AKN2019'][4]['mixed']
        # Bonferroni uses two-sided 0.025 each: marginal normality only.
        normal_tail = erfc(float(Z)/sqrt(2))
        assert abs(normal_tail-.025) < 2e-16
        gap = CAL['Rb2020'][0]-CAL['Cs2018'][0]-Z*(CAL['Rb2020'][1]+CAL['Cs2018'][1])
        assert gap>0
        assert_close(gap,author['alpha_closure']['positive_interval_gap'])
        # L2 triangle bound on standard deviation, not a deterministic error.
        delta = CAL['Rb2020'][0]-CAL['Cs2018'][0]
        z_ind = dec(delta)/dec(CAL['Rb2020'][1]**2+CAL['Cs2018'][1]**2).sqrt()
        z_worst = delta/(CAL['Rb2020'][1]+CAL['Cs2018'][1])
        assert z_worst == F(80,19)
        assert_close(z_ind,author['alpha_closure']['independent_Gaussian_z_diagnostic'])
        predictions = {}
        for version,(a10,u10) in VERSIONS.items():
            predictions[version] = {}
            for name,(q,uq) in CAL.items():
                value,deriv,terms,um = evaluate(q,a10)
                pub=author['predictions'][version][name]
                assert_close(value,pub['ae'],'1e-27')
                assert_close(deriv,pub['da_d_inverse_alpha'],'1e-29')
                for own,old in zip(terms,pub['qed_terms_times_1e12']):
                    assert_close(own*10**12,old,'1e-15')
                t=1/(PI*q)
                u_alpha=abs(deriv)*uq
                u_terms=[u10*t**5,UH,UW,um]
                theory_linear=sum(u_terms)
                theory_rss=dec(sum(x*x for x in u_terms)).sqrt()
                pred_rss=(dec(u_alpha)**2+theory_rss**2).sqrt()
                residual=AE-value
                zz=dec(residual)/(dec(UE)**2+pred_rss**2).sqrt()
                for x,k in [(um,'known_mass_sd_linear_times_1e12'),
                            (u_alpha,'alpha_sd_times_1e12'),
                            (theory_linear,'theory_sd_linear_upper_times_1e12'),
                            (residual,'exp_minus_prediction_times_1e12')]:
                    assert_close(x*10**12,pub[k],'1e-20')
                assert_close(zz,pub['independent_Gaussian_z_diagnostic_only'],'1e-23')
                predictions[version][name]={'ae_times_1e12':string(value*10**12),
                    'residual_times_1e12':string(residual*10**12),
                    'known_mass_sd_times_1e12':string(um*10**12),
                    'theory_sd_linear_times_1e12':string(theory_linear*10**12),
                    'independent_gaussian_z_only':string(zz)}
        # An incomplete Set V substitution is not the declared sensitivity.
        assert VERSIONS['AH2025_sensitivity'][0] == F('5.86958')
        assert VERSIONS['AH2025_sensitivity'][0] != F('6.800')
        assert F('1.7030e-12')-HAD == F('.0100e-12')
        return {'all_checks_passed':True,'review_not_new_science_group':True,
            'method':'independent primary-number transcription, Fraction Horner and derivative; no author import',
            'author_sha256':{p:sha256((HERE/p).read_bytes()).hexdigest() for p in ('public_inputs.json','qed_recoil_check.py','results.json')},
            'normal_two_sided_marginal_tail':normal_tail,
            'alpha_interval_gap_exact':str(gap),'alpha_interval_gap_decimal':string(gap),
            'independent_gaussian_calibration_z':string(z_ind),
            'joint_gaussian_any_correlation_min_z_exact':str(z_worst),
            'predictions':predictions,
            'bounded_scope':{'omitted_order_hard_bound':False,'recoil_fully_qed_independent':False,
                'joint_residual_independence_asserted':False,'higher_dipole_family_excluded':False}}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--record',action='store_true')
    args=parser.parse_args()
    r=run()
    out=HERE/'physical_review_checks.json'
    if args.record:
        with out.open('x',encoding='utf-8') as f:
            json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
    elif out.exists():
        assert r==json.loads(out.read_text(encoding='utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
