"""Bounds for a heavy, minimally coupled two-singlet dipole. No fitting or sampling."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def read(name):
    return json.loads((BASE/name).read_text(encoding="utf-8-sig"))


def g_enclosure(x, n=192):
    """Directed integral bounds. x >= 1, all partitions and summands rational."""
    assert x >= 1

    def primitive(t):
        return t**2/2 - 2*t**3/3 + t**4/4

    lower = upper = Fraction(0)
    for i in range(n):
        left, right = Fraction(i, n), Fraction(i+1, n)
        area = 3*x*(primitive(right)-primitive(left))/2
        lower += area/(1+(x-1)*right)
        upper += area/(1+(x-1)*left)
    return float(lower), float(upper)


def certificate():
    data = read("inputs.json")
    weak = read(data["previous_weak"])
    t_center, t_width = weak["tau"]["g_tau_over_e_leptonic"]
    min_heavy = data["heavy_min_GeV"]
    min_w, max_w = data["W_mass_guard_GeV"]
    x = Fraction(min_heavy, max_w)**2
    gl, gh = g_enclosure(x)
    xn = float(x)
    analytic = -(2*xn**3+5*xn**2-xn)/(4*(1-xn)**3)-3*xn**3*math.log(xn)/(2*(1-xn)**4)
    assert gl < analytic < gh
    light_cap = data["light_mass_guard_eV"]
    mass_data = read(data["previous_masses"])["NO"]
    assert max(mass_data["solar"][2], mass_data["atm_abs"][2]) < light_cap**2
    delta = (light_cap/(min_heavy*1e9))**2
    light_amp = (light_cap/(min_w*1e9))**2/4
    # max joint effect of the two heavy columns departing from saturation,
    # with arbitrary phases. F_max is the infinite-mass value 1/2.
    loss = (gl+.5)/2*math.sqrt(delta) + light_amp
    prefactor = 3/(2*math.pi*data["inverse_alpha_upper"])
    meg = data["MEG_II_observed_upper"]
    stress = data["stress"]

    def bounds(extra_t=0., amp_reduction=0., width_cap=1., dilation=1.):
        t = t_center-t_width-extra_t
        assert 1 < t < math.sqrt(2)
        root_product = math.sqrt(2)*(1-1/t**2)
        amp = (1-amp_reduction)*max(0., gl*root_product-loss)
        lower = prefactor*amp**2/width_cap
        pmax = (math.sqrt(meg*dilation*width_cap/prefactor)/(1-amp_reduction)+loss)/gl
        return dict(t_lower=t, row_product_lower=root_product**2, amplitude_lower=amp,
                    branching_lower=lower, experiment_ratio=lower/(dilation*meg),
                    maximal_sqrt_row_product=pmax,
                    inferred_t_upper=1/math.sqrt(1-pmax/math.sqrt(2)))

    basic = bounds()
    guarded = bounds(stress["tau_e_additional_absolute_error"], stress["destructive_fraction_of_leading_amplitude"],
                     stress["inclusive_width_relative_upper"], stress["limit_dilation"])
    # alpha(0) < 1/137. Even without evaluating a phase, this upper bound
    # suffices for the exact two-width example of 1408.
    prior = read(data["previous_witness"])["full_mass_example"]
    lam = max(prior["X_eigenvalues"])
    witness_br_max = 3/(2*math.pi*137)*(.5*lam+light_amp)**2/(1-lam)**2
    assert basic["experiment_ratio"] > 700 and guarded["experiment_ratio"] > 100
    assert witness_br_max < meg
    return dict(gamma_threshold_enclosure=[gl,gh], gamma_closed_form=analytic,
                conservative_emu_minor_upper=delta, light_amplitude_upper=light_amp,
                central_dictionary=basic, stress_dictionary=guarded,
                two_width_example_branching_upper=witness_br_max,
                minimal_destructive_fraction_baseline=1-math.sqrt(meg/basic["branching_lower"]))


def compare(actual, expected):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys()
        for k in actual:
            compare(actual[k], expected[k])
    else:
        np.testing.assert_allclose(actual, expected, rtol=3e-12, atol=0)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write", action="store_true")
    result = certificate()
    if cli.parse_args().write:
        with (BASE/"results.json").open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        compare(result, read("results.json"))
    print(json.dumps(result, allow_nan=False))
