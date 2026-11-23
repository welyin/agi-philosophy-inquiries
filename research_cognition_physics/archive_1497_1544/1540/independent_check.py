"""Independent mass, Euclidean loop and physical two-body width tests for round 1540."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
J0 = math.pi**2 / 3
SQUARE = math.pi**2 / 6 - math.log(2)**2


def lower_product():
    cfg = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    osc = json.loads((HERE / cfg["oscillation_inputs"]).read_text(encoding="utf-8"))
    widths = cfg["widths"]
    amp, _, amp_err = osc["Daya_Bay_A13"]
    theta = (1 - math.sqrt(1 - amp - widths * amp_err)) / 2
    solar, err = cfg["IO_solar_angle"]
    splitting, lower_err, _ = osc["T2K_mu_disapp_dm13_IO_eV2"]
    block = (1 - theta) * (1 - 2 * (solar + widths * err)) * math.sqrt(splitting - widths * lower_err)
    ratios = [cfg[name] for name in ("leptonic_tau_over_e_HFLAV", "leptonic_mu_over_e_HFLAV")]
    upper_sum = sum((cfg["lower_to_true_ratio_transfer_upper"] / (v - widths * dv))**2
                    for v, dv in ratios)
    # Solve the unexpanded leptonic width inequality by bisection.
    lo, hi = 0.0, 1.0
    for _ in range(80):
        s = (hi + lo) / 2
        if (2 + 2 * s + s*s) / (1 + s*s/4) <= upper_sum:
            lo = s
        else:
            hi = s
    contact = math.sqrt(2) * cfg["G0_over_Gmu_upper"] * cfg["muon_Fermi_constant_GeVm2"] * hi
    prod = (block * 1e-9 * 48 * math.pi**2 /
            (cfg["physical_mass_block_over_two_loop_bound_upper"] * cfg["electron_pole_mass_GeV"]**2 * contact))
    return block, prod


def integral(a, b, c, order):
    # Two ordered half-planes in log(X) and log(Y/X); integrate the symmetric kernel directly.
    u, weight = np.polynomial.legendre.leggauss(order)
    x = np.exp(25 * u)[:, None]
    y = x * np.exp(12.5 * (u - 1))[None, :]
    kernel = 2 / (x + y + c + np.sqrt((x - y)**2 + 2*c*(x+y) + c*c))
    lines = 1/((x+a)*(y+b)) + 1/((x+b)*(y+a))
    return float(312.5 * np.sum(weight[:, None] * weight[None, :] * x*y * kernel * lines))


def phase_space(first, second):
    if first + second >= 1:
        return 0.
    return math.sqrt((1 - (first + second)**2) * (1 - (first - second)**2))


def compute():
    block, product = lower_product()
    author = json.loads((HERE / "results.json").read_text(encoding="utf-8"))["finite_transfers"]
    np.testing.assert_allclose([block, product], [author["B_min_eV"], author["mu_times_g_ee_GeV_lower"]])
    np.testing.assert_allclose(integral(1, 1, 0, 144), J0, atol=2e-7, rtol=0)
    u, w = np.polynomial.legendre.leggauss(64)
    z = (u+1)/2
    np.testing.assert_allclose(np.dot(w, np.log1p(z)/(z*(1+z))), SQUARE, atol=1e-13)
    cases = {}
    for ratio in (2.5, 3., 4.):
        rho = 1 - (1 - ratio**-2) * SQUARE / J0
        actual = integral(1, 1, ratio**2, 144) / J0
        np.testing.assert_allclose(actual, integral(1, 1, ratio**2, 96)/J0, atol=2e-6, rtol=0)
        assert 0 < actual < rho < 1
        beta = math.sqrt(1 - 4/ratio**2)
        geom = product/rho * math.sqrt(beta)/(8*math.pi)
        cases[str(ratio)] = {"rho_upper": rho, "single_mass_integral_ratio": actual,
                             "mu_times_g_ee_GeV_lower": product/rho,
                             "geometric_mean_width_GeV_lower": geom,
                             "total_width_GeV_lower": 2*geom,
                             "parent_mass_GeV_lower_at_one_percent_width": 200*geom}
    # Mixed spectrum: numbers in units of the parent pole mass for phase space.
    fractions = np.array([.2, .35, .4])
    weights = np.array([.2, .3, .5])
    a = (fractions / fractions.max())**2
    c = fractions.max()**-2
    matrix = np.array([[integral(ai, aj, c, 144) for aj in a] for ai in a])
    assert np.linalg.eigvalsh(matrix).min() > 0
    ratio = float((weights @ matrix @ weights)/(J0*np.sum(weights/a)))
    rho = cases["2.5"]["rho_upper"]
    assert ratio < rho
    phases = np.array([[phase_space(fi, fj) for fj in fractions] for fi in fractions])
    beta_weighted = float(weights @ phases @ weights)
    pair_product = product / rho
    parent = 5e6
    g2 = pair_product * math.sqrt(beta_weighted) / parent
    mu2 = pair_product**2/g2
    gamma_ee = parent*g2/(8*math.pi)
    gamma_hh = mu2*beta_weighted/(8*math.pi*parent)
    np.testing.assert_allclose(gamma_ee, gamma_hh)
    assert gamma_ee > cases["2.5"]["geometric_mean_width_GeV_lower"]
    assert phase_space(.55, .55) == 0 and integral(1, 1, 1/.55**2, 144) > 0
    return {"mass_block_eV_lower": block, "baseline_product_GeV_lower": product,
            "square_integral": SQUARE, "open_spectra": cases,
            "mixed_example": {"integral_over_baseline_bound": ratio, "phase_space_weight": beta_weighted,
                              "equal_partial_width_GeV_at_product_floor": gamma_ee},
            "closed_pair_counterexample": True}


def verify(current, saved):
    if isinstance(current, dict):
        assert current.keys() == saved.keys()
        for key, val in current.items():
            verify(val, saved[key])
    else:
        np.testing.assert_allclose(current, saved, rtol=2e-8, atol=1e-11)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    values = compute()
    path = HERE / "independent_results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(values, stream, indent=2)
            stream.write("\n")
    else:
        verify(values, json.loads(path.read_text(encoding="utf-8")))
    print("1540 independent mass, loop and two-channel width tests passed")
