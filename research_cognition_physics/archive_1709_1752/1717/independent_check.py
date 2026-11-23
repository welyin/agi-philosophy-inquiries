"""Independent weak-charge readout of the five-singlet reciprocal down sector."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


def singlet_fraction(coupling, k, yb, vev, bare):
    """One light family and its two signed heavy paths."""
    mixing = yb / math.sqrt(coupling**2 - yb**2)
    a = coupling * vev / bare
    w = np.sqrt([(1 + k) / 2, (1 - k) / 2])
    matrix = np.array([[0, a*w[0], a*w[1]],
                       [-mixing*w[0], 1, 0],
                       [-mixing*w[1], 0, -1]])
    left, mass, _ = np.linalg.svd(matrix)
    light_mass_squared = mass[-1]**2
    odds = a**2 * ((1-k*k)/(1-light_mass_squared)**2
                   + k*k/(1+mixing**2-light_mass_squared)**2)
    prob = np.sum(left[1:, -1]**2)
    assert light_mass_squared < 1 and abs(odds/(1+odds)-prob) < 2e-14
    assert prob/(1-prob) >= a*a*(1-k*k)*(1-1e-12)
    return float(prob)


def inverse_calibration(proton, electron, cs, z, n):
    h = proton + electron
    rho = 2*h + (z*proton - cs)/n
    x = h/rho
    angle = (1 + electron/rho)/4
    assert rho > 0 and 0 < angle < 1
    vector = 1-4*angle
    np.testing.assert_allclose(
        [rho*(vector+x), -rho*vector, rho*(z*vector-n+(z+2*n)*x)],
        [proton, electron, cs], rtol=2e-14, atol=2e-14)
    return x, rho, angle


def fit_windows(cs, widths, proton, pe, ne, reference, z, n, g_cap, tests):
    """The two Cs determinations are alternative reductions, each used on its own."""
    d_p = proton[0]
    p_err = proton[1]
    c_err = math.hypot(*widths)
    z0 = 1-4*reference
    e0 = -z0+4*(pe-reference)
    means = np.array([z0+d_p, e0, z*z0-n+cs])
    sigmas = np.array([p_err, 4*ne, c_err])
    out = {
        "fixed_reference_pulls_as_marginal_widths": {
            name: [(x-d_p)/p_err, ((z+2*n)*x-cs)/c_err]
            for name, x in tests.items()},
        "rho_one_free_angle_minimum_widths": {
            name: abs((cs-z*d_p)-2*n*x)/(c_err+z*p_err)
            for name, x in tests.items()}}
    for count in (1, 2):
        d = (cs-z*d_p)/(2*n)
        spread = count*(c_err+z*p_err)/(2*n)
        two = [max(0, d-spread), min(1, d+spread)]
        corners = [inverse_calibration(*(means+count*np.array(v)*sigmas), z, n)
                   for v in itertools.product((-1., 1.), repeat=3)]
        low, high = min(a[0] for a in corners), max(a[0] for a in corners)
        assert low < 0 < high < 1
        out[str(count)] = {
            "epsilon_from_p_and_cs_free_angle_fixed_rho": two,
            "g_upper_necessary_from_two_targets": g_cap(two[1]),
            "epsilon_from_all_three_free_rho_and_angle": [0, high],
            "g_upper_necessary_from_three_targets": g_cap(high)}
    return out


def compute():
    weak = read("../1385/inputs.json")
    atom = read("../1462/inputs.json")
    electron = read("../1568/inputs.json")["low_MSbar_targets"]["electron"]
    observed, p_error = weak["proton"]["pv_electron_scattering_only"]
    standard, pbase_error = weak["proton"]["SM"]
    proton = [observed-standard, math.hypot(p_error, pbase_error)]
    z, n = (weak["caesium"][name] for name in ("Z", "N"))
    ref_cs, baseline_error = atom["baseline_Cs"]
    source = read("../1713/results.json")
    yb, _, yd = source["tagged_yukawa_values_desc"]
    kd = -yd/yb
    bare = read("../1713/inputs.json")["mediator_matching_GeV"]
    vev = read("../1492/inputs.json")["v246_GeV"] / math.sqrt(2)
    trials = read("inputs.json")["counterfactual_g_values"]
    tests = {str(g): singlet_fraction(g, kd, yb, vev, bare) for g in trials}
    def g_cap(upper):
        return bare/vev * math.sqrt(upper/((1-upper)*(1-kd*kd)))
    sigma = tests["1"]
    out = {"d_singlet_fractions": tests,
           "fixed_reference_at_g_one": {"proton": standard+sigma, "Cs133": ref_cs+(z+2*n)*sigma,
                                        "neutron_shift": 2*sigma,
                                        "Cs_minus_Z_proton_shift": 2*n*sigma},
           "data": {"Qweak_p_offset_and_width": proton, "electron_MSbar_angle_and_width": electron,
                    "Cs_common_baseline_and_error": [ref_cs, baseline_error]}}
    ref_angle = .23868
    for key, (obs, ex, th) in atom["atomic_branches"].items():
        out[key] = fit_windows(obs-ref_cs, (ex, th, baseline_error), proton, *electron,
                               ref_angle, z, n, g_cap, tests)
    assert out["Dzuba_2012"]["rho_one_free_angle_minimum_widths"]["1"] < 1
    assert 1 < out["Sahoo_2021"]["rho_one_free_angle_minimum_widths"]["1"] < 2
    for key in atom["atomic_branches"]:
        assert all(sigma <= out[key]["1"]["epsilon_from_all_three_free_rho_and_angle"][1]
                   for sigma in tests.values())
    return out


def round_values(x):
    if isinstance(x, dict):
        return {k: round_values(v) for k, v in x.items()}
    if isinstance(x, list):
        return [round_values(v) for v in x]
    return float(f"{x:.11g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = round_values(compute())
    path = HERE/"independent_results.json"
    if args.check:
        assert result == read(path)
        print("1717 independent weak-charge projections verified")
    else:
        path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")

