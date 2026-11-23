"""Three singlets, one Higgs doublet, one-loop mass matching (round 1558)."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOOP = 8 * np.pi**2


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def reconstruct(p, u, ev):
    f = np.zeros((3, 3), complex)
    i, j = np.triu_indices(3)
    f[i, j] = np.array(p["F_real_upper"]) + 1j * np.array(p["F_imag_upper"])
    f[j, i] = f[i, j]
    h = np.diag(p["H_diag"])
    s = complex(p["s_real"], p["s_imag"])
    r = p["r"]
    y = r * (h - 3 * s * f)
    majorana = -u**2 * r**2 * f / (ev * p["t_eV"])
    light = -u**2 * y @ np.linalg.solve(majorana, y.T)
    return y, majorana, light


def takagi(matrix):
    _, masses, vh = np.linalg.svd(matrix)
    axes = vh.conj().T[:, ::-1]
    phases = np.diag(axes.T @ matrix @ axes)
    axes *= np.exp(-0.5j * np.angle(phases))
    np.testing.assert_allclose(axes.T @ matrix @ axes, np.diag(masses[::-1]),
                               rtol=1e-10, atol=2e-15 * masses[0])
    return masses[::-1], axes


def summarize(p, parent, u, ev, mh, mhmax, wmin, boost):
    y, majorana, light = reconstruct(p, u, ev)
    heavy, axes = takagi(majorana)
    masses, light_axes = takagi(light)
    np.testing.assert_allclose(masses[0] / ev * 1000, parent["m0_meV"], rtol=2e-8)
    yy = y @ axes
    norms2 = np.sum(abs(yy)**2, axis=0)
    terms = heavy**2 * norms2 / LOOP
    np.testing.assert_allclose(terms.sum(), np.linalg.norm(y @ majorana.conj().T, "fro")**2 / LOOP)
    exact_min = (heavy**3 @ masses[::-1]) / (LOOP * u**2)
    assert exact_min <= terms.sum() * (1 + 1e-10)
    # Exhibit a physical seesaw with the same full light matrix and all three
    # heavy masses at the minimum, then apply a complex orthogonal 1-2 boost.
    permuted_axes = light_axes.conj()[:, ::-1]
    for angle in (0, boost):
        rotation = np.eye(3, dtype=complex)
        a, b = np.cosh(angle), 1j * np.sinh(angle)
        rotation[:2, :2] = [[a, b], [-b, a]]
        y_free = 1j * (permuted_axes * np.sqrt(masses[::-1])) @ rotation * np.sqrt(heavy) / u
        np.testing.assert_allclose(-u**2 * (y_free / heavy) @ y_free.T, light,
                                   atol=2e-12 * masses[-1], rtol=1e-8)
        analytic = exact_min + sum(heavy[:2]**3) * sum(masses[::-1][:2]) * np.sinh(angle)**2 / (LOOP * u**2)
        np.testing.assert_allclose(np.sum(abs(y_free * heavy)**2) / LOOP, analytic)
    a = (mh**2 / terms.sum())**(1 / 3)
    rescaled = u**2 * (np.sqrt(a) * y) @ np.linalg.solve(a * majorana, np.sqrt(a) * y.T)
    np.testing.assert_allclose(rescaled, -light, rtol=1e-9)
    return dict(heavy_masses_GeV=heavy.tolist(), light_masses_meV=(masses / ev * 1000).tolist(),
                heavy_column_norms_squared=norms2.tolist(), threshold_terms_GeV2=terms.tolist(),
                sum_GeV2=float(terms.sum()), least_positive_threshold_over_Higgs_window=float(wmin * terms.sum() / mhmax**2),
                freely_reoriented_min_GeV2=float(exact_min),
                freely_reoriented_min_over_Higgs_window=float(wmin * exact_min / mhmax**2),
                fixed_point_max_common_response=float(mhmax**2 / terms.sum()),
                common_M_factor_preserving_light_for_central_Higgs=float(a),
                rescaled_heavy_masses_GeV=(a * heavy).tolist(),
                rescaled_Y_spectral_norm=float(np.sqrt(a) * np.linalg.norm(y, 2)))


def compute():
    cfg = read("inputs.json")
    path = HERE / cfg["tensors"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == cfg["tensors_sha256"]
    contract, osc = read(cfg["contract"]), read(cfg["oscillation"])
    mh, sigma = read(cfg["higgs"])["atlas"]["mass_GeV"]
    u = contract["v246_GeV"] / np.sqrt(2)
    ev = cfg["ev_to_GeV"]
    count = cfg["halfwidths"]
    max_mh = mh + count * sigma
    lo_w, hi_w = contract["positive_Higgs_threshold_transport"]
    model = read(cfg["tensors"])["candidates"]
    rows = {name: summarize(row["tensors"], row["predictions"], u, ev, mh, max_mh, lo_w,
                            cfg["boost_test"]) for name, row in model.items()}
    # m_(2) is the middle of the three mass eigenvalues in either ordering;
    # the source atmospheric conventions differ between NO and IO.
    budgets = {}
    for ordering in ("NO", "IO"):
        sg, sgerr = osc[ordering]["solar_gap"]
        ag, agerr = osc[ordering]["atmospheric_gap"]
        solar, atmos = sg - count * sgerr, ag - count * agerr
        if ordering == "NO":
            mmax, middle = np.sqrt(solar + atmos), np.sqrt(solar)
        else:
            mmax, middle = np.sqrt(solar + atmos), np.sqrt(atmos)
        available = LOOP * u**2 * max_mh**2 / (lo_w * ev)
        budgets[ordering] = dict(
            mmax_lower_eV=float(mmax), middle_lower_eV=float(middle),
            M1_upper_GeV=float(np.cbrt(available / (mmax + middle))),
            M2_upper_GeV=float(np.cbrt(available / middle)),
        )
    # Check the squared-mass normalization against the previous independent
    # finite six-state check, which used u=174 GeV instead of v/sqrt(2).
    historical = read("../1556/independent_results.json")["A_charged_first_family_held_out"]
    np.testing.assert_allclose(np.array(rows["release_down_A"]["heavy_masses_GeV"]) *
                               (cfg["original_u_example_GeV"] / u)**2,
                               historical["finite_check"]["right_singular_GeV"], rtol=3e-8)
    return dict(contract=dict(u_GeV=float(u), mH_upper_GeV=max_mh, transfer=[lo_w, hi_w]),
                from_oscillation_only=budgets, candidates=rows)


def check(saved, actual):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            check(saved[key], actual[key])
    else:
        np.testing.assert_allclose(actual, saved, rtol=1e-9, atol=1e-35)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = cli.parse_args()
    out = compute()
    path = HERE / "results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    else:
        check(read("results.json"), out)
    print("1558 verified")
