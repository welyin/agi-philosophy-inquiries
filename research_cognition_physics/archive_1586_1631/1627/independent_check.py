"""Independent L3 recoil-window check and separate 1998 inclusive comparison."""

import argparse
import json
from itertools import product
from math import atanh, cos, exp, lgamma, log, pi, radians, sqrt
from pathlib import Path

import numpy as np


ENERGY = 188.6
LUMINOSITY = 176.4
ALPHA = 1 / 137.035999
GEV2_TO_PB = 389379338.0
COS_MAX = cos(radians(43.0))
OUTPUT = Path(__file__).with_name('independent_results.json')
# The last three beam energy ranges have published luminosities but no event-by-event weights.
BEAMS = ((188.6, 188.6, 176.0), (191.6, 191.6, 29.5), (195.5, 195.5, 83.9),
         (199.5, 199.5, 81.3), (201.7, 201.7, 34.8), (202.5, 205.5, 74.8),
         (205.5, 207.2, 130.2), (207.2, 209.2, 8.6))


def production_pb(mass, coupling, energy=ENERGY, z_max=COS_MAX):
    s = energy**2
    r = mass**2 / s
    z, weights = np.polynomial.legendre.leggauss(48)
    z *= z_max
    t = -0.5 * (s - mass**2) * (1 - z)
    u = -0.5 * (s - mass**2) * (1 + z)
    # The electron chiral couplings have opposite signs but equal squared modulus.
    g_eff_sq = coupling**2 / 2
    squared_amplitude = 2 * (4 * pi * ALPHA) * g_eff_sq * (
        t / u + u / t + 2 * s * mass**2 / (t * u)
    )
    quadrature = z_max * float(weights @ squared_amplitude) * (1 - r) / (32 * pi * s)
    log_ratio = 2 * atanh(z_max)
    integral = (1 - r) * (2 * log_ratio - 2 * z_max) + 4 * r * log_ratio / (1 - r)
    analytic = ALPHA * g_eff_sq / (2 * s) * integral
    assert abs(quadrature / analytic - 1) < 5e-13
    return quadrature * GEV2_TO_PB


def poisson_cdf_leq(count, mean):
    # Sum from the observed count down. This avoids starting with a tiny P(N=0).
    term = exp(-mean + count * log(mean) - lgamma(count + 1))
    cumulative = term
    for k in range(count, 0, -1):
        term *= k / mean
        cumulative += term
    return cumulative


def upper_mean(count, tail):
    low, high = count, 2 * count + 16.0
    assert poisson_cdf_leq(count, low) >= tail >= poisson_cdf_leq(count, high)
    for _ in range(70):
        mid = (low + high) / 2
        if poisson_cdf_leq(count, mid) > tail:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def inclusive_1998():
    upper = upper_mean(297, 0.05)
    assert abs(poisson_cdf_leq(297, upper) - 0.05) < 5e-13
    corners = []
    for mass, coupling in product((56, 64), (0.62, 0.65)):
        sigma = production_pb(mass, coupling)
        corners.append({'mass_GeV': mass, 'gF': coupling, 'Born_barrel_pb': sigma,
                        'Born_barrel_events_before_decay': LUMINOSITY * sigma})
    born_min = min(x['Born_barrel_events_before_decay'] for x in corners)
    thresholds = []
    for k in (0, 1, 1.10):
        signal_upper = (sqrt(upper) + sqrt(k * 288.9))**2
        thresholds.append({'k': k, 'signal_events_upper': signal_upper,
                           'B_times_transfer_bound_at_weakest_point': signal_upper / born_min,
                           'transfer_bound_for_B_one_third_at_weakest_point':
                           3 * signal_upper / born_min})
    return {'inputs': {'energy_GeV': ENERGY, 'luminosity_pb_inv': LUMINOSITY,
                       'alpha': ALPHA, 'GeV_minus2_in_pb': GEV2_TO_PB,
                       'barrel_abs_cos_max': COS_MAX, 'observed_count': 297,
                       'sm_neutrino_events': 288.9, 'total_expected_events': 292.2},
            'corners': corners, 'poisson_upper_95_mean': upper,
            'thresholds': thresholds}


def calculate():
    rows = []
    for lo, hi, lum in BEAMS:
        # For fixed angle the shape increases with r, hence decreases with s at fixed mass.
        sigma_min = production_pb(56, 0.62, hi, 0.730)
        sigma_max = production_pb(64, 0.65, lo, 0.730)
        rows.append({'energy_low_GeV': lo, 'energy_high_GeV': hi, 'luminosity_pb_inv': lum,
                     'sigma_bounds_pb': [sigma_min, sigma_max],
                     'born_events_bounds': [lum * sigma_min, lum * sigma_max]})
    born_min, born_max = [sum(row['born_events_bounds'][i] for row in rows) for i in range(2)]
    upper = upper_mean(2, 0.05)
    assert abs((1 + upper + upper**2 / 2) * exp(-upper) - 0.05) < 1e-14
    invisible_degrees = 1
    charged_degrees = 2
    branch = invisible_degrees / (invisible_degrees + charged_degrees)
    caps = []
    for sm_cap in (0.0, 2.0, 2.2, 2.42):
        limit = (sqrt(upper) + sqrt(sm_cap))**2
        caps.append({'SM_interfering_events_cap': sm_cap, 'squared_signal_events_cap': limit,
                     'B_times_transfer_at_Born_min': limit / born_min,
                     'transfer_at_Born_min_for_Bnu_one_third': limit / (born_min * branch)})
    return {'combined_2004': {'energy_rows': rows, 'luminosity_sum_pb_inv': sum(i[2] for i in BEAMS),
                              'observed_by_abs_cos_bin': [1, 1, 0, 0],
                              'sm_expectations_by_bin': [0.5, 0.5, 0.4, 0.6],
                              'sm_efficiencies_percent_by_bin': [82, 80, 81, 79],
                              'born_events_bounds': [born_min, born_max], 'U95_observed_2': upper,
                              'Bnu_massless_e_and_electron_neutrino_only': branch,
                              'caps': caps},
            'inclusive_1998_overlapping_not_combined': inclusive_1998()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true', help='compare to the saved independent result')
    args = parser.parse_args()
    result = calculate()
    if args.check:
        assert result == json.loads(OUTPUT.read_text(encoding='utf-8'))
    else:
        OUTPUT.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
