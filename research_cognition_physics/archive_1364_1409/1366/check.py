"""Joint two-channel Higgs comparison in the CMS decay-mode fit."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


DIR = Path(__file__).resolve().parent


def gaussian_contrast(y, covariance, q):
    direction = np.array([q, 1.0])
    inv = np.linalg.inv(covariance)
    normalization = (direction @ inv @ y) / (direction @ inv @ direction)
    residual = y - normalization * direction
    contrast = np.array([1.0, -q])
    z = (contrast @ y) / math.sqrt(contrast @ covariance @ contrast)
    assert normalization > 0 and np.isclose(z * z, residual @ inv @ residual)
    return dict(q=q, best_common_normalization=float(normalization), signed_z=float(z),
                delta_chi2=float(z * z))


def ratio_interval(y, cov, threshold):
    # Fieller region for the linear contrast in the adopted local normal approximation.
    a = y[1] ** 2 - threshold * cov[1, 1]
    b = -2 * y[0] * y[1] + 2 * threshold * cov[0, 1]
    c = y[0] ** 2 - threshold * cov[0, 0]
    assert a > 0 and b * b > 4 * a * c
    bounds = sorted(np.roots([a, b, c]).tolist())
    assert bounds[0] > 0
    return bounds


def run():
    source = json.loads((DIR / "inputs.json").read_text(encoding="utf-8-sig"))
    cms, masses = source["cms"], source["masses"]
    measured_mu, measured_tau = cms["decay_mu"], cms["decay_tau"]
    y = np.array([measured_mu["central"], measured_tau["central"]])
    sigma_mu = (measured_mu["minus"] + measured_mu["plus"]) / 2
    sigma_tau = (measured_tau["minus"] + measured_tau["plus"]) / 2
    rho = cms["printed_correlation"]
    cov = np.array([[sigma_mu ** 2, rho * sigma_mu * sigma_tau],
                    [rho * sigma_mu * sigma_tau, sigma_tau ** 2]])
    mu = masses["muon_MeV"] / 1000
    tau = masses["tau_MeV"] / 1000
    mh = cms["higgs_mass_GeV"]
    ratio_mass = mu / tau
    ratio_tree = ratio_mass ** 2 * ((1 - (2 * mu / mh) ** 2) / (1 - (2 * tau / mh) ** 2)) ** 1.5
    ratio_ref = cms["reference_mumu_percent"] / cms["reference_tautau_percent"]
    q_hat = y[0] / y[1]
    tests = []
    for q in source["relative_width_hypotheses"]:
        test = gaussian_contrast(y, cov, q)
        pulls = []
        for sigma in [measured_mu["minus"], measured_mu["plus"]]:
            for corr in [rho - 0.005, rho + 0.005]:
                variance = sigma ** 2 + q * q * sigma_tau ** 2 - 2 * q * corr * sigma * sigma_tau
                pulls.append(abs(y[0] - q * y[1]) / math.sqrt(variance))
        test["abs_z_sensitivity_at_corners"] = [min(pulls), max(pulls)]
        tests.append(test)
    box_ratio = [(y[0] - measured_mu["minus"]) / (y[1] + measured_tau["plus"]),
                 (y[0] + measured_mu["plus"]) / (y[1] - measured_tau["minus"])]
    return dict(
        conventions="mu, tau order; signed z is a local Gaussian diagnostic",
        low_energy_mass_ratio=ratio_mass, on_shell_tree_width_ratio=ratio_tree,
        cms_reference_partial_width_ratio=ratio_ref, inferred_width_ratio_central=q_hat * ratio_ref,
        ratio_of_signal_strengths=q_hat, relative_abs_coupling_central=math.sqrt(q_hat),
        gaussian_covariance=cov.tolist(), gaussian_tests=tests,
        gaussian_q_delta_chi2_1=ratio_interval(y, cov, 1),
        gaussian_q_delta_chi2_3_841459=ratio_interval(y, cov, 3.841459),
        marginal_one_sigma_box_projected_q=box_ratio,
        union_bound_if_each_nominal_coverage_0_6827=2 * 0.6827 - 1,
    )


def same(actual, expected):
    if isinstance(actual, dict):
        return actual.keys() == expected.keys() and all(same(v, expected[k]) for k, v in actual.items())
    if isinstance(actual, list):
        return len(actual) == len(expected) and all(same(a, b) for a, b in zip(actual, expected))
    if isinstance(actual, float):
        return bool(np.isclose(actual, expected, atol=1e-12, rtol=1e-10))
    return actual == expected


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = run()
    path = DIR / "results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        assert same(result, json.loads(path.read_text(encoding="utf-8")))
        print("PASS: CMS two-lepton rate contrast, covariance profile and conditional bounds")


if __name__ == "__main__":
    main()

