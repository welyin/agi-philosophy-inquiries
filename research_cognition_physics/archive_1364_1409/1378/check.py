"""Finite electroweak scalar polarizations, a threshold inverse and OPAL sensitivity."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def integral(argument):
    """Subtracted scalar shape, positive for positive timelike argument below 4."""
    if argument == 4:
        return 2/9
    if not argument < 4:
        raise ValueError("threshold domain exceeded")
    # Positive half of the even integral after y=2*x-1.
    return float(np.sum(WEIGHTS * NODES**2 *
                         -np.log1p(-argument*(1 - NODES**2)/4)))


NODES, WEIGHTS = np.polynomial.legendre.leggauss(80)
NODES, WEIGHTS = (NODES + 1)/2, WEIGHTS/2


def energy_ratio(z, fraction):
    return 1/fraction if z == 0 else integral(z) / integral(fraction*z)


def inverse_ratio(ratio, fraction):
    lower, upper = 0., 4.
    if not 1/fraction < ratio < energy_ratio(upper, fraction):
        raise ValueError("No unique finite mass in the open interval")
    for _ in range(80):
        mid = (lower + upper)/2
        if energy_ratio(mid, fraction) < ratio:
            lower = mid
        else:
            upper = mid
    return (lower + upper)/2


def real_bounds(ratio, fraction, q_high, eps):
    low = ratio*(1-eps)/(1+eps)
    high = ratio*(1+eps)/(1-eps)
    masses = [q_high/math.sqrt(inverse_ratio(r, fraction)) for r in (high, low)]
    return masses


def trace_powers(twice_t, twice_y, copies):
    m = np.arange(-twice_t, twice_t+1, 2)/2
    y = twice_y/2
    q = m + y
    return (copies*np.asarray([m@m, len(m)*y*y, sum(m*y), q@q, sum(q**4)])).tolist()


def certificate():
    cfg = json.loads(Path(__file__).with_name("inputs.json").read_text(encoding="utf-8"))
    design, opal = cfg["timelike_design"], cfg["opal"]
    low, high = design["energies_GeV"]
    mass = design["reference_mass_GeV"]
    fraction, z = (low/high)**2, (high/mass)**2
    truth = energy_ratio(z, fraction)
    assert math.isclose(high/math.sqrt(inverse_ratio(truth, fraction)), mass)
    assert np.all(np.diff([energy_ratio(s, fraction) for s in np.linspace(0, 3.999, 71)]) > 0)
    series = sum(math.factorial(n)**2 / (n*math.factorial(2*n+1)*(2*n+3))*z**n
                 for n in range(1, 66))
    assert math.isclose(series, integral(z), rel_tol=1e-13)
    moment, mimic = [trace_powers(t, y, copies) for copies, t, y in
                     cfg["group_models"]["multiplets_copies_twiceT_twiceY"]]
    assert moment[:4] == mimic[:4] and moment[4] == 980 and mimic[4] == 56
    # Observable is Delta alpha(-b)-Delta alpha(-a), not an equality
    # between OPAL's log-slope fit and any arbitrary threshold template.
    a, b = opal["positive_minus_t_endpoints_GeV2"]
    factor = moment[3] / (cfg["electromagnetic_normalization"]["inverse_alpha_Thomson"] * 4*math.pi)
    sigma_opal = math.sqrt(sum(x*x for x in opal["uncertainties_stat_syst_theory"]))
    masses = cfg["threshold_comparison"]["masses_GeV"]
    space = []
    for m in masses:
        change = factor*(integral(-a/m**2) - integral(-b/m**2))
        bound = factor*(b-a)/(30*m**2)
        assert 0 < change <= bound*(1+1e-13)
        space.append(dict(M_GeV=m, endpoint_change=change, upper=bound,
                          fraction_of_OPAL_error=change/sigma_opal))
    # Distinct masses can generate the same two-energy ratio when the
    # normalization is free and the first point fixes the weights.
    alternatives = design["mass_mixture_comparison_GeV"]
    ratios = [energy_ratio((high/m)**2, fraction) for m in alternatives]
    weight = (truth-ratios[1])/(ratios[0]-ratios[1])
    assert 0 < weight < 1 and np.isclose(weight*ratios[0] + (1-weight)*ratios[1], truth)
    return dict(
        su2_u1_mixed_Q2_Q4=dict(septet=moment, doublets=mimic),
        timelike_synthetic=dict(energies_GeV=[low, high], mass_GeV=mass,
             F_at_two_energies=[integral(fraction*z),integral(z)], ratio=truth,
             ratio_domain=[1/fraction, energy_ratio(4, fraction)],
             mass_range_at_one_percent_per_extracted_amplitude_GeV=real_bounds(
                 truth, fraction, high, design["relative_error_each_extracted_extra_amplitude"]),
             example_mixture_masses_GeV=alternatives, fraction_of_low_point_from_lighter=weight),
        opal=dict(observed_change=opal["published_delta_alpha"], quoted_SM=opal["paper_SM_prediction"],
                  quoted_SM_error=opal["paper_SM_error"],
                  experimental_and_theory_errors=opal["uncertainties_stat_syst_theory"],
                  combined_size=sigma_opal, examples=space))


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for k in actual:
            compare(actual[k], saved[k])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for a,b in zip(actual, saved):
            compare(a,b)
    else:
        assert math.isclose(actual, saved, rel_tol=2e-11, abs_tol=1e-48), (actual, saved)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = certificate()
    path = Path(__file__).with_name("results.json")
    if args.write:
        with path.open("x", encoding="utf-8") as fp:
            fp.write(json.dumps(report, indent=2)+"\n")
    else:
        compare(report, json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(report))
