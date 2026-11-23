"""Recompute the round 1392 anomaly, period and haloscope certificates."""
import argparse
from fractions import Fraction
from functools import reduce
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def anomaly(fields):
    # Each multiplet gives its electric charge, color dimension, 2T and X_L-X_R.
    color = sum(index * axial for charge, dimension, index, axial in fields)
    em = sum(2 * dimension * charge**2 * axial
             for charge, dimension, index, axial in fields)
    return color, em


def wall_number(color, vev_charges):
    stabilizer = reduce(math.gcd, vev_charges)
    assert color % stabilizer == 0
    return abs(color) // stabilizer


def certificate():
    inputs = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    qcd, admx, capp = (inputs[name] for name in ("qcd", "admx", "capp"))
    down = (Fraction(-1, 3), 3, 1, 1)
    lepton = (Fraction(1), 1, 0, 1)
    n, e = anomaly([down])
    n_comp, e_comp = anomaly([down, lepton])
    assert (n, e, n_comp, e_comp) == (1, Fraction(2, 3), 1, Fraction(8, 3))
    # Explicit full G6 bundle coefficients in the integer hypercharge normalization.
    q_integer = int(6 * down[0])
    k1 = down[1] * q_integer**2
    k1_comp = k1 + 36
    assert (24 * n + k1) % 36 == (24 * n_comp + k1_comp) % 36 == 0
    c = qcd["photon_subtraction"]
    dc = qcd["photon_subtraction_width"]
    mf = qcd["m_times_f_GeV2"]
    dmf = sum(qcd["m_times_f_errors_GeV2"])
    alpha = 1 / qcd["inverse_alpha_zero_rounded"]
    photon_ratio = e / n

    def slope(ratio, subtraction, mass_scale):
        return alpha / (2 * math.pi * mass_scale) * (float(ratio) - subtraction)

    center = slope(photon_ratio, c, mf)
    corners = [abs(slope(photon_ratio, sub, scale)) for sub, scale
               in itertools.product((c - dc, c + dc), (mf - dmf, mf + dmf))]
    interval = [min(corners), max(corners)]
    mass = capp["example_mass_microeV"] * 1e-15
    capp_g = [k * mass for k in interval]
    fractions = [(capp["published_approx_limit_per_GeV"] / g)**2 for g in capp_g]
    alternative = slope(e_comp / n_comp, c, mf)
    same_power_fraction = (alternative / center)**2
    assert 0 < same_power_fraction < 1 and center < 0 < alternative
    to_mass = qcd["h_eV_second"] * 1e12
    dfsz = admx["dfs_g_gamma"]
    power_common_f = ((c - float(photon_ratio)) / (2 * dfsz))**2
    frac_common_f_box = [(2 * dfsz / (sub - float(photon_ratio)))**2
                          for sub in (c + dc, c - dc)]
    periods = {
        "one_linear_S": wall_number(n, [1]),
        "same_rescaled_labels": wall_number(7 * n, [7]),
        "two_phase_locked_vevs": wall_number(2, [2, 1])
    }
    assert list(periods.values()) == [1, 1, 2]
    return {
        "N": n, "E_over_N": str(photon_ratio),
        "E_over_N_with_extra_lepton": str(e_comp / n_comp),
        "primitive_period_2pi_units_for_one_S": 1,
        "QCD_vacua_per_primitive_orbit": periods,
        "G6_K1_down_and_with_lepton": [k1, k1_comp],
        "G6_unit_orbit_u_coefficient": [(24 * n + k1) // 36, (24 * n_comp + k1_comp) // 36],
        "g_over_mass_GeV_minus2_signed": center,
        "abs_g_over_mass_theory_rectangle_GeV_minus2": interval,
        "extra_lepton_g_over_mass_signed": alternative,
        "target_fraction_for_full_extra_lepton_equal_halo_power": same_power_fraction,
        "ADMX_nominal_rest_mass_window_microeV": [f * to_mass for f in admx["published_frequency_MHz"]],
        "ADMX_veto_nominal_rest_mass_microeV": [f * to_mass for f in admx["veto_MHz"]],
        "ADMX_common_f_candidate_over_DFSZ_power_center": power_common_f,
        "ADMX_common_f_fraction_matching_DFSZ_theory_rectangle": frac_common_f_box,
        "CAPP_at_example_mass_abs_g_GeV_minus1": capp_g,
        "CAPP_at_example_mass_central_abs_g_GeV_minus1": abs(center) * mass,
        "CAPP_representative_fraction_threshold_working_range": sorted(fractions),
        "CAPP_representative_fraction_threshold_center":
            (capp["published_approx_limit_per_GeV"] / (abs(center) * mass))**2
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = certificate()
    path = HERE / "results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
            f.write("\n")
    else:
        saved = json.loads(path.read_text(encoding="utf-8"))
        assert result.keys() == saved.keys()
        for name, value in result.items():
            if isinstance(value, list):
                assert len(value) == len(saved[name])
                assert all(math.isclose(x, y, rel_tol=2e-14, abs_tol=0)
                           for x, y in zip(value, saved[name]))
            elif isinstance(value, float):
                assert math.isclose(value, saved[name], rel_tol=2e-14, abs_tol=0)
            else:
                assert value == saved[name]
    print(json.dumps(result, ensure_ascii=False))
