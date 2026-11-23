"""Electromagnetic spectrum tests using published OPAL, L3 and CLEO observables."""
import argparse
import hashlib
import json
from math import atanh, hypot, log, pi, sqrt
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
ALPHA0_INV = 137.03599976
GEV2_TO_PB = 389379338.0
# OPAL table 17, ratios only: energy, fit, upper error, lower error, SM.
OPAL = [
    (188.6, 126.2, 4.2, 3.8, 127.9), (191.6, 138.1, 14.3, 10.6, 127.9),
    (195.5, 129.9, 7.6, 6.4, 127.9), (199.5, 130.3, 8.0, 6.7, 127.8),
    (201.6, 134.0, 12.2, 9.4, 127.8), (204.9, 123.7, 6.8, 5.8, 127.8),
    (206.6, 117.9, 4.7, 4.2, 127.8),
]
# CLEO table VII, endpoints: energy, R, statistical, common, uncorrelated.
CLEO = [(6.964, 3.597, .033, .057, .020), (10.538, 3.591, .003, .067, .049)]
X, WG = np.polynomial.legendre.leggauss(128)


def bisect_zero(fn, lo, hi):
    left = fn(lo)
    assert left * fn(hi) <= 0
    for _ in range(75):
        mid = (lo + hi) / 2
        if fn(mid) * left > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def beta(energy, mass):
    return sqrt(max(0., 1. - (2 * mass / energy) ** 2))


def pair(energy, mass):
    speed = beta(energy, mass)
    return speed * (3 - speed ** 2) / 2


def timelike(energy, mass):
    """Real inverse-alpha decrement per unit electric weight above threshold."""
    assert energy >= 2 * mass > 0
    v = beta(energy, mass)
    return (-8/3 + v*v + v*(3-v*v)*atanh(v)) / (3*pi)


def independent_integral(energy, mass):
    """Resolve both integrable logarithmic poles by a symmetry and squared endpoint maps."""
    pole = (1 - beta(energy, mass)) / 2
    t = (X+1)/2
    total = 0.
    for direction, span in [(-1, pole), (1, .5-pole)]:
        x = pole + direction * span*t*t
        # The factorised form avoids a subtraction near x=pole.
        log_abs = np.log((energy/mass)**2 * span*t*t * np.abs(1-pole-x))
        total += np.dot(WG/2, 2*span*t * x*(1-x)*log_abs)
    return 4*total/pi


def spacelike_interval(low_q2, high_q2, mass):
    x = (X+1)/2
    u = x*(1-x)
    return np.dot(WG, u*np.log1p((high_q2-low_q2)*u/(mass**2+low_q2*u)))/pi


def opal_result():
    energy = 193.2
    single = bisect_zero(lambda m: timelike(energy, m), 10., energy/2)
    minimum = .85*OPAL[0][0]
    maximum = OPAL[-1][0]
    broad = bisect_zero(lambda m: timelike(minimum, m) + timelike(maximum, m), 10., minimum/2)
    at_edges = [timelike(e, broad) for e in [minimum, maximum]]
    assert abs(sum(at_edges)) < 1e-14
    integration_error = max(abs(independent_integral(e, m)-timelike(e, m))
                            for e, m in [(energy, single), (minimum, broad), (maximum, broad)])
    assert integration_error < 2e-8
    return {
        "ratio_fits_energy_center_plus_minus_SM": OPAL,
        "combined_130_207_SM_template_only": [energy, 126.7, 2.4, 2.3, 127.9],
        "single_zero_mass_GeV": single, "pair_kernel_at_zero": pair(energy, single),
        "seven_bin_true_ISR_energy_interval_GeV": [minimum, maximum],
        "balanced_above_threshold_mass_GeV": broad, "real_shift_per_weight_edges": at_edges,
        "real_shifts_per_weight_nominal_energies": [timelike(row[0], broad) for row in OPAL],
        "max_shift_divided_by_smallest_separate_width": at_edges[1]/min(min(row[2:4]) for row in OPAL),
        "minimum_pair_kernel_on_true_ISR_range": pair(minimum, broad),
        "independent_quadrature_max_error": integration_error,
    }


def opal_cut(mass):
    # OPAL table 18, no lepton universality. Extra absolute transfer allowance in each j.
    rows = {}
    for extra in (0., .02):
        fb_low = .763-2*.070-extra
        tot_high = -.080+2*.044+extra
        rows[str(extra)] = fb_low/(fb_low+tot_high)
    energy = 188.6
    born = 4*pi/(3*ALPHA0_INV**2*energy**2)*GEV2_TO_PB
    # Table 5 errors; add table 3 luminosity in quadrature conservatively even if included.
    width = hypot(.4, .1, .0021*22.)
    excess = 22.+2*width-22.2
    factor = rows["0.02"]
    floor = born * pair(energy, mass)*factor/3
    e_vector = sqrt(-.080+2*.044+.02)
    e_axial = sqrt(.763-2*.070-.02)
    electron = np.array([e_vector+e_axial,e_vector-e_axial], dtype=complex)
    optimal = -np.vdot(electron,np.ones(2))/np.vdot(electron,electron)
    amplitude = np.ones(2)+optimal*electron
    assert abs(np.vdot(amplitude,amplitude).real/2-factor)<1e-12
    assert np.vdot(np.ones(2)+(optimal+.17j)*electron,
                   np.ones(2)+(optimal+.17j)*electron).real/2 > factor
    return {
        "electron_S_matrix_jfb_jtot_plus_widths": [.763,.070,-.080,.044],
        "min_photon_fraction_by_extra_absolute_matching_in_each_j": rows,
        "hadronic_energy_GeV_events_obs_stat_sys_lumi_fraction_reference_pb":
            [energy, 4019, 22., .4, .1, .0021, 22.2],
        "working_two_width_excess_over_reference_pb": excess,
        "unit_charge_born_massless_pb": born,
        "mass_GeV": mass, "minimal_color_signal_pb_before_transfer": floor,
        "efficiency_transfer_required_at_reference": excess/floor,
        "tolerable_downward_ordinary_rate_shift_at_unit_transfer_pb": floor-excess,
        "mass_upper_for_minimal_color_at_unit_transfer_and_reference_GeV":
            bisect_zero(lambda m: born*pair(energy,m)*factor/3-excess, 0., energy/2),
    }


def continuum_result():
    low, high = CLEO
    common_sigma = hypot(low[2], high[2], low[4], high[4], high[3]-low[3])
    unaligned_sigma = hypot(low[2],high[2],low[4],high[4],high[3]+low[3])
    diff = high[1]-low[1]
    headroom = diff + 2*common_sigma
    def critical_drift(m):
        assert low[0] <= 2*m <= high[0]
        return pair(high[0],m)/3-headroom
    return {
        "input_ascending_energy_R_stat_common_uncorrelated": CLEO,
        "high_minus_low_R": diff, "one_shared_standardized_systematic_sigma": common_sigma,
        "arbitrary_orientation_common_sigma_upper": unaligned_sigma,
        "two_width_positive_room_before_background_drift": headroom,
        "two_width_room_no_common_orientation": diff+2*unaligned_sigma,
        "max_background_downward_drift_by_mass_unit_transfer_weight_one_third":
            {str(m):critical_drift(m) for m in [3.5,4.,4.5,5.]},
        "masses_tested_above_GeV": low[0]/2,
        "max_mass_by_chosen_background_decline_unit_transfer_GeV":
            {str(d):bisect_zero(lambda m:critical_drift(m)-d,low[0]/2,high[0]/2)
             for d in [0.,.05,.10]},
    }


def calculate():
    opal = opal_result()
    small = spacelike_interval(12.25,3434.,400.)/3
    return {
        "sources": {"OPAL":"https://arxiv.org/pdf/hep-ex/0309053",
                    "L3":"https://arxiv.org/pdf/hep-ex/0002035",
                    "CLEO":"https://arxiv.org/pdf/0706.2813"},
        "OPAL_vacuum_polarization": opal,
        "OPAL_direct_hadronic_rate": opal_cut(opal["balanced_above_threshold_mass_GeV"]),
        "L3_spacelike": {"Q2_GeV2": [12.25,3434.], "measured_difference_width": [3.80,1.29],
                         "all_mass_ceiling_per_unit_weight": log(3434./12.25)/(3*pi),
                         "mass_400_GeV_minimal_color_weight": small,
                         "one_published_width_over_mass_400_signal": 1.29/small},
        "CLEO_direct_continuum": continuum_result(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check = parser.parse_args().check
    result = calculate()
    # Standardize floating-point output for reproducible text across numerical backends.
    def rounded(item):
        if isinstance(item, dict):
            return {k:rounded(v) for k,v in item.items()}
        if isinstance(item, (tuple,list)):
            return [rounded(v) for v in item]
        if isinstance(item, float):
            return round(item,11)
        return item
    document = json.dumps(rounded(result),indent=2,ensure_ascii=False)+"\n"
    path = ROOT/"photon.json"
    if check:
        if path.read_text(encoding="utf-8") != document:
            raise SystemExit("Published photon.json differs")
    else:
        path.write_text(document,encoding="utf-8")
    print("1785 photon",hashlib.sha256(document.encode()).hexdigest()[:12])
