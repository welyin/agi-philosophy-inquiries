"""Check the down-quark phase directly with the full six-state scalar propagator."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
POINTS, WEIGHTS = np.polynomial.legendre.leggauss(96)
POINTS = (POINTS + 1) / 2
WEIGHTS = WEIGHTS / 2


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def named_directions(cfg):
    entry = (HERE / cfg["family"]).resolve()
    spec = importlib.util.spec_from_file_location("actual_family_1521", entry)
    family = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(family)
    frame, mother, right, _, left = family.family(read(cfg["family_inputs"]))
    _, _, kernel_rows = np.linalg.svd(mother, full_matrices=True)
    kernel = kernel_rows[-1].conj()
    np.testing.assert_allclose(mother @ kernel, 0, atol=1e-15)
    ell = left.conj().T @ frame[-1]
    r = right.conj().T @ kernel
    return ell / np.linalg.norm(ell), r / np.linalg.norm(r)


def terms(mu, ell, r):
    phase = ell * r.conj()
    return np.array([[(mu[j] / mu[i] - mu[i] / mu[j]) *
                      (phase[i] * phase[j].conj()).imag for j in range(3)]
                     for i in range(3)]) * np.triu(np.ones((3, 3)), 1)


def uniform_bound(mu, invariant_terms):
    # Applies to any x,y >= 0 for unit a/M, including the diagonal F terms.
    u2 = max(mu) ** 2
    w_floor = 1 / (1 + u2 + u2 ** 2)
    determinant_floor = w_floor ** 2 - u2
    condition = abs(invariant_terms).sum() / abs(invariant_terms.sum())
    return (condition * (1 - w_floor ** 2) + 1 - determinant_floor) / determinant_floor


def woodbury(mu, ell, r, x, y):
    summands = terms(mu, ell, r)
    tree = summands.sum()
    total = leading = 0.0
    for t, weight in zip(POINTS, WEIGHTS):
        ai = 1 + t * mu ** 2
        bi = t * mu
        w = 1 / (ai ** 2 - bi ** 2)
        alpha, beta = ai * w, -bi * w
        aa = (abs(ell) ** 2) @ alpha
        bb = (abs(r) ** 2) @ alpha
        cross = np.vdot(ell, beta * r)
        determinant = (1 + t * x * aa) * (1 + t * y * bb) - t ** 2 * x * y * abs(cross) ** 2
        area = weight * t ** 2 * (1 - t) * x * y
        total += area * (w @ summands @ w) / determinant
        leading += area * tree / ((1 + t * x) * (1 + t * y))
    return float(total), float(leading)


def spectral_kernel(masses):
    # F(s) = s log(s)/(s-1)-1; avoid subtracting two nearly equal numbers.
    value = []
    for s in masses:
        h = s - 1
        if abs(h) < 0.1:
            val = sum((-1) ** (n + 1) * h ** n / (n * (n + 1)) for n in range(1, 32))
        else:
            val = s * math.log1p(h) / h - 1
        value.append(val)
    return np.array(value)


def six_state(mu, ell, r, x, y):
    common = np.diag(1 + mu ** 2)
    mass = np.block([[common + x * np.outer(ell, ell.conj()), np.diag(mu)],
                     [np.diag(mu), common + y * np.outer(r, r.conj())]])
    eigen, modes = np.linalg.eigh(mass)
    assert eigen.min() > 0
    kernel = (modes * spectral_kernel(eigen)) @ modes.conj().T
    return float((kernel[:3, 3:].diagonal() / mu).imag.sum()), float(eigen.min())


def calculate():
    cfg = read("inputs.json")
    assert cfg["gluino_over_M"] == cfg["real_Aeff_over_M"] == 1
    scale = cfg["soft_scale_GeV"]
    # The low-energy SM mass matching length is independently set to the
    # displayed conventional value for this finite six-state completion.
    vev = 246.22
    table = np.array(cfg["down_yukawa_dr_1TeV"])
    mu = table[:, 0] * vev / (math.sqrt(2) * scale)
    ell, r = named_directions(cfg)
    summands = terms(mu, ell, r)
    g3 = read(cfg["top_transport"])["sm_dr_at_one_TeV"]["g3"][0]
    prefactor = g3 ** 2 / (6 * math.pi ** 2)
    author = read("results.json")
    np.testing.assert_allclose(summands.sum(), author["invariant_at_one_TeV"], rtol=1e-12)
    bound = uniform_bound(mu, summands)

    # Enclose the condition number over the whole registered rectangular mass
    # range using independent interval terms; each mass ratio is monotone.
    lower, upper = (table[:, 0] + s * 2 * table[:, 1] for s in (-1, 1))
    factor = ell * r.conj()
    abs_sum, min_sum, max_sum = 0., 0., 0.
    for i, j in ((0, 1), (0, 2), (1, 2)):
        ratios = np.array([lower[j] / upper[i], upper[j] / lower[i]])
        values = (factor[i] * factor[j].conj()).imag * (ratios - 1 / ratios)
        abs_sum += np.max(abs(values))
        min_sum += np.min(values)
        max_sum += np.max(values)
    assert min_sum * max_sum > 0
    box_condition = abs_sum / min(abs(min_sum), abs(max_sum))
    ub2 = (upper[-1] * vev / (math.sqrt(2) * scale)) ** 2
    w2 = (1 + ub2 + ub2 ** 2) ** -2
    floor = w2 - ub2
    box_bound = (box_condition * (1 - w2) + 1 - floor) / floor
    # All three contributions have one sign throughout the ordered mass boxes.
    assert (summands[np.triu_indices(3, 1)] < 0).all() and min(lower[1:] - upper[:-1]) > 0
    sign_specific_bound = max(1 - w2, 1 / floor - 1)

    checks = {}
    for x, y in ((0, .8), (.6, 0), (.001, .001), (.01, .01), (.1, .1),
                 (.27, .41), (1., 1.), (2., 3.)):
        exact, leading = woodbury(mu, ell, r, x, y)
        spectral, eigen_min = six_state(mu, ell, r, x, y)
        np.testing.assert_allclose(spectral, exact, atol=1.8e-10, rtol=2e-8)
        relative = abs(exact / leading - 1) if leading else 0.
        assert relative <= bound
        if x == y:
            numerical = author["theta_at_equal_positive_fractional_mass_squared_splits"]
            entry = numerical[str(int(x) if x == 1 else x)]
            np.testing.assert_allclose(abs(leading) * prefactor, entry, rtol=2e-12)
        checks[f"{x},{y}"] = {
            "full_signed_theta": exact * prefactor,
            "first_LR_signed_theta": leading * prefactor,
            "relative_change": relative,
            "minimum_squark_squared_mass_over_M2": eigen_min,
        }
    return {
        "vev_GeV": vev, "light_masses_GeV": (mu * scale).tolist(),
        "ds_db_sb_contributions": summands[np.triu_indices(3, 1)].tolist(),
        "all_nonnegative_xy_central_relative_bound": bound,
        "two_width_all_nonnegative_xy_relative_bound": box_bound,
        "two_width_condition_upper_bound": box_condition,
        "two_width_one_sign_uniform_relative_bound": sign_specific_bound,
        "examples": checks,
    }


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    else:
        np.testing.assert_allclose(a, b, rtol=2e-8, atol=1e-14)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    modes = cli.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--write-new", action="store_true")
    args = cli.parse_args()
    result = calculate()
    if args.write_new:
        with (HERE / "independent_results.json").open("x", encoding="utf8") as output:
            json.dump(result, output, indent=2)
            output.write("\n")
    else:
        compare(result, read("independent_results.json"))
    print("Checked six squark modes, rank-two resolvent, and uniform phase bound.")
