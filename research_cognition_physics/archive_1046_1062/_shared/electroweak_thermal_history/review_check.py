"""Independent check of the published central freeze equation, not lattice data.

Uses Newton iteration instead of importing the author's bisection routine.
--freeze exclusively writes the first independent record; default is read-only.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import math

HERE = Path(__file__).resolve().parent
SHARED = HERE.parent
RECEIPT = HERE / 'review_checks.json'
AUTHORS = [
    'electroweak_thermal_history_adoption.md',
    'electroweak_thermal_history/check.py',
    'electroweak_thermal_history/results.json',
    'electroweak_thermal_history/sources.md',
    '../_admission/after1059_asymmetry_history/selection.md',
]


def run(freeze=False):
    saved = json.loads((HERE / 'results.json').read_text(encoding='utf-8'))
    expected = {'fit_slope_per_GeV': '0.83', 'fit_intercept': '-147.7',
                'alpha_H': '0.1015', 'g_star': '106.75',
                'reduced_Planck_mass_GeV_diagnostic': '2.435e18'}
    for key, value in expected.items():
        assert saved['inputs'][key] == value, key
    # Independently rewrite Gamma/T^3 = alpha_H H as exp(aT+b) = beta*T.
    beta = 0.1015 * math.sqrt(math.pi**2 * 106.75 / 90) / 2.435e18
    def residual(t):
        return 0.83*t - 147.7 - math.log(beta*t)
    t = 132.0
    for _ in range(6):
        t -= residual(t) / (0.83 - 1/t)
    delta = t - saved['freeze_central_root_GeV']
    assert residual(130) < 0 < residual(140)
    assert 0.83 - 1/130 > 0
    assert abs(residual(t)) < 1e-12 and abs(delta) < 1e-10
    assert not saved['lattice_or_systematic_error_recomputed']
    assert not saved['baryon_yield_computed']
    record = {
        'method': 'Newton iteration of independently written logarithmic equation',
        'central_root_GeV': t,
        'author_root_difference_GeV': delta,
        'log_residual_at_root': residual(t),
        'log_residual_at_130': residual(130),
        'log_residual_at_140': residual(140),
        'difference_from_published_rounded_estimate_GeV': t-131.7,
        'source_uncertainty_reestimated': False,
        'no_baryon_yield_or_lattice_simulation': True,
        'new_science_groups': 0,
        'author_sha256': {rel: sha256((SHARED/rel).read_bytes()).hexdigest() for rel in AUTHORS},
        'all_checks_passed': True,
    }
    if freeze:
        with RECEIPT.open('x', encoding='utf-8') as stream:
            json.dump(record, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        assert record == json.loads(RECEIPT.read_text(encoding='utf-8'))
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--freeze', action='store_true')
    args = parser.parse_args()
    record = run(args.freeze)
    print(json.dumps({'all_checks_passed': True, 'new_science_groups': 0,
                      'central_root_GeV': record['central_root_GeV'],
                      'author_root_difference_GeV': record['author_root_difference_GeV']}))
