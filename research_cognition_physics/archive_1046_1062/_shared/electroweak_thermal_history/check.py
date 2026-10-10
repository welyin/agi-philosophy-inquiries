"""Read-only central-equation check; not a lattice or baryon-yield calculation.

The optional --write-results exclusively creates the first diagnostic record.
All energies and rates below use hbar=c=kB=1, with GeV as the numerical unit.
No published systematic uncertainty is reinterpreted as a Gaussian fit error.
"""
from pathlib import Path
import argparse
import json
import math

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results.json'
INPUTS = {
    'rate_source': 'arXiv:1404.3565v1, equations (7) and (9)',
    'fit_slope_per_GeV': '0.83',
    'fit_intercept': '-147.7',
    'alpha_H': '0.1015',
    'g_star': '106.75',
    'reduced_Planck_mass_GeV_diagnostic': '2.435e18',
    'primary_fit_window_GeV': ['140', '155'],
    'paper_supported_extension_down_to_GeV': '130',
    'published_freeze_temperature_GeV': '131.7',
    'published_freeze_error_GeV_not_recomputed': '2.3',
    'crossover_maintext_GeV': '159.6',
}


def compute():
    a = float(INPUTS['fit_slope_per_GeV'])
    b = float(INPUTS['fit_intercept'])
    alpha_h = float(INPUTS['alpha_H'])
    gstar = float(INPUTS['g_star'])
    mbar = float(INPUTS['reduced_Planck_mass_GeV_diagnostic'])
    hcoef = math.sqrt(math.pi**2 * gstar / 90) / mbar

    def log_ratio(t):
        # log[(Gamma/T^3)/(alpha_H*H)], where Gamma has energy dimension four.
        return a*t + b - math.log(t) - math.log(alpha_h*hcoef)

    lo, hi = 130.0, 140.0
    assert log_ratio(lo) < 0 < log_ratio(hi)
    derivative_lower = a - 1 / lo
    assert derivative_lower > 0
    for _ in range(64):
        mid = (lo + hi) / 2
        if log_ratio(mid) < 0:
            lo = mid
        else:
            hi = mid
    tstar = (lo + hi) / 2
    hubble = hcoef*tstar**2
    gamma = tstar**4*math.exp(a*tstar+b)
    ratio = gamma/tstar**3/(alpha_h*hubble)
    assert abs(log_ratio(tstar)) < 1e-12
    assert abs(ratio-1) < 1e-12
    assert tstar < 140 < float(INPUTS['crossover_maintext_GeV'])
    return {
        'inputs': INPUTS,
        'natural_unit_energy_dimensions': {'Gamma': 4, 'Gamma_over_T3': 1, 'H': 1},
        'freeze_central_root_GeV': tstar,
        'root_bracket_GeV': [lo, hi],
        'log_ratio_derivative_lower_per_GeV_on_130_140': derivative_lower,
        'log_ratio_at_root': log_ratio(tstar),
        'Gamma_GeV4_at_root': gamma,
        'H_GeV_at_root': hubble,
        'rate_over_alphaH_H_at_root': ratio,
        'rate_over_alphaH_H_at_130': math.exp(log_ratio(130)),
        'rate_over_alphaH_H_at_140': math.exp(log_ratio(140)),
        'root_minus_published_center_GeV_diagnostic': tstar-131.7,
        'crossover_minus_central_freeze_GeV_diagnostic': 159.6-tstar,
        'central_root_uses_paper_extension_below_primary_fit_window': True,
        'lattice_or_systematic_error_recomputed': False,
        'baryon_yield_computed': False,
        'new_science': 0,
        'new_empirical_group': 0,
        'passed': True,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.write_results:
        with RESULTS.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
    else:
        saved = json.loads(RESULTS.read_text(encoding='utf-8'))
        assert saved == result, 'Diagnostic changed; preserve old evidence.'
    print(json.dumps({'passed': True, 'new_science': 0,
                      'freeze_central_root_GeV': result['freeze_central_root_GeV'],
                      'readonly': not args.write_results}))
