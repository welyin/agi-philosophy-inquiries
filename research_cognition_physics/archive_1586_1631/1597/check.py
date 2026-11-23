"""Calibrate a charged-dipole difference with the fine structure of neutral helium-4."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8-sig"))


def envelope(values):
    values = tuple(values)
    return [min(values), max(values)]


def scaled_interval(bounds, slopes):
    return envelope(x * y for x, y in itertools.product(bounds, slopes))


def closure(upper, lower1, lower2):
    residual = upper[0] - lower1[0] - lower2[0]
    radius = 2 * (upper[1] + lower1[1] + lower2[1])
    return {"signed_residual_Hz": residual, "two_width_radius_Hz": radius,
            "minimum_additional_total_Hz": max(0, abs(residual) - radius)}


def remaining_component(norm_squared, fixed, unknown_factor, fixed_factor):
    # Invert x²+y²+(a*x+b*y)² >= norm_squared at given y=fixed.
    c = norm_squared - (1 + fixed_factor ** 2) * fixed ** 2
    if c <= 0:
        return 0.
    linear = unknown_factor * fixed_factor * fixed
    return c / (math.sqrt(linear ** 2 + (1 + unknown_factor ** 2) * c) + linear)


def run():
    inputs = read("inputs.json")
    old = {k: read(v) for k, v in inputs["oldRefs"].items()}
    he = inputs["helium"]
    f01 = he["FOSOF_01_2024_Hz"]
    f02 = he["laser_02_2017_Hz"]
    f12 = he["laser_12_2017_Hz"]
    microwave12 = he["FOSOF_12_2018_Hz"]
    closure_results = {"laser_and_FOSOF01": closure(f02, f01, f12),
                       "FOSOF01_and_FOSOF12": closure(f02, f01, microwave12)}
    assert closure_results["laser_and_FOSOF01"]["minimum_additional_total_Hz"] == 838

    ref, theory_width = he["theory_01_Hz"]
    sensitivity = [ref * a for a in he["sensitivity_exponents"]]
    table_slope = 1.e6 * sum((level0 - level1) * power for power, (level0, level1)
                             in enumerate(he["table_III_intervals_MHz"], start=2))
    assert sensitivity[0] <= table_slope <= sensitivity[1]
    h_ref = he["reference_inverse_alpha"]
    residual = f01[0] - ref
    radius = 2 * (f01[1] + theory_width) + he["theory_print_half_step_Hz"]
    log_bounds = envelope(d / k for d, k in itertools.product(
        [residual - radius, residual + radius], sensitivity))
    assert max(map(abs, log_bounds)) < he["local_abs_log_alpha_max"]
    dip = old["dipole"]
    offset = float(dip["inverse_alpha_offset"])
    unit = dip["offset_units"]
    ref_offset = (h_ref - offset) / unit
    he_offsets = envelope(ref_offset + h_ref * math.expm1(-z) / unit for z in log_bounds)
    trap, trap_error = dip["offsets"]["trap_with_historical_SM"]
    slopes = [unit * s for s in dip["historical_abs_da_per_inverse_alpha"]]
    ae_He = scaled_interval([he_offsets[0] - (trap + 2 * trap_error),
                             he_offsets[1] - (trap - 2 * trap_error)], slopes)
    common_zero_residual = residual - table_slope * math.log(h_ref / (offset + unit * trap))
    assert abs(common_zero_residual) < 2 * theory_width
    assert ae_He[0] < 0 < ae_He[1]

    m_e, m_mu = dip["me_GeV"], dip["mmu_GeV"]
    ratio = (m_e / m_mu) ** 2
    mu_center, mu_width = old["muon"]["WP25_difference_as_published"]
    mu = [1e-11 * (mu_center - 2 * mu_width), 1e-11 * (mu_center + 2 * mu_width)]
    assert mu[0] < 0 < mu[1]
    io = old["osc"]["IO"]
    m2 = math.sqrt(io["atm_abs"][1])
    s13 = io["s13"][2]
    s23 = old["direction"]["theta23_sin2"]["IO"][2]
    a = math.sqrt(s13 / ((1 - s13) * (1 - s23)))
    b = math.sqrt(s23 / (1 - s23))
    bridge = old["bridge_results"]
    row_norm = bridge["normalized_e_row_norm_sq_floor"]
    generic = bridge["generic_inverted_floor_eV"]
    center_ext = inputs["prospective_inverse_alpha_offset"]
    recoil = old["recoil"]["recoil_ae_rectangles_one_two"]
    scenarios = {}
    for case in old["bridge"]["scenarios"][1:]:
        label = case["label"]
        off = bridge["cases"][label]["offdiag_allowance_GeV_minus2"]
        amplitude = off + case["defect_per_m2_GeV_minus2"]
        ae_extra = case["other_ae"]
        mu_near_extra = case["other_amu"] + inputs["muon_alpha_transport_amplitude_near_recoil"]
        mu_wide_extra = case["other_amu"] + inputs["muon_alpha_transport_amplitude_wide_He"]
        windows = {}
        for target_meV in inputs["target_mass_meV"]:
            h = remaining_component(row_norm, target_meV / (1000 * m2), b, a)
            gap_required = amplitude / h
            margin = 2 * m_e ** 2 * gap_required + ae_extra + ratio * mu_near_extra
            blind = [ratio * mu[0] - margin, ratio * mu[1] + margin]
            assert blind[0] < 0 < blind[1]
            # Positive and negative thresholds use the smallest absolute slope.
            inverse_cut = [trap - 2 * trap_error + blind[0] / slopes[0],
                           trap + 2 * trap_error + blind[1] / slopes[0]]
            allowed = center_ext - inverse_cut[1]
            freq_budget = (sensitivity[0] * math.log1p(allowed * unit /
                           (offset + inverse_cut[1] * unit))) if allowed > 0 else 0.
            windows[str(target_meV)] = {"gap_required_GeV_minus2": gap_required,
                "blind_anomaly_interval": blind, "external_inverse_alpha_offset_cuts": inverse_cut,
                "positive_candidate_offset_room_signed": allowed,
                "positive_candidate_total_He_budget_Hz": freq_budget}

        def mass_from_ae(bounds, transport):
            signed = [(bounds[0] - ae_extra - ratio * (mu[1] + transport)) / (2 * m_e ** 2),
                      (bounds[1] + ae_extra - ratio * (mu[0] - transport)) / (2 * m_e ** 2)]
            gap = max(0., signed[0], -signed[1])
            h = min(1., amplitude / gap) if gap > 0 else 1.
            floor = max(generic, m2 * remaining_component(row_norm, h, a, b))
            return {"signed_gap_GeV_minus2": signed, "minimum_absolute_gap_GeV_minus2": gap,
                    "IO_mass_floor_meV": 1000 * floor}

        actual = {"Rb_two_widths": mass_from_ae(recoil["Rb"][1], mu_near_extra),
                  "Cs_two_widths_signed_extension": mass_from_ae(recoil["Cs"][1], mu_near_extra),
                  "He_two_widths": mass_from_ae(ae_He, mu_wide_extra)}
        assert actual["Cs_two_widths_signed_extension"]["signed_gap_GeV_minus2"][1] < 0
        assert actual["Rb_two_widths"]["minimum_absolute_gap_GeV_minus2"] == 0
        scenarios[label] = {"target_windows": windows, "published_inputs": actual}
    return {"frequency_closure": closure_results, "slope_Hz_per_log_alpha_table": table_slope,
            "registered_slope_interval": sensitivity, "helium_two_width_budget_Hz": radius,
            "helium_log_alpha_interval_relative_to_reference": log_bounds,
            "helium_inverse_alpha_offsets": he_offsets, "helium_electron_anomaly_interval": ae_He,
            "zero_both_dipoles_required_atomic_offset_Hz": common_zero_residual,
            "mass_ratio_squared": ratio, "scenarios": scenarios,
            "shift_of_both_inverse_cuts_for_alternative_trap_origin":
                inputs["alternative_trap_offset"] - trap}


def match(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            match(actual[key], saved[key])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for value, reference in zip(actual, saved):
            match(value, reference)
    else:
        assert math.isclose(actual, saved, rel_tol=2.e-10, abs_tol=1.e-23), (actual, saved)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check", action="store_true")
    action.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    out = run()
    if args.write_new:
        with (HERE / "results.json").open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(out, indent=2, allow_nan=False) + "\n")
    elif args.check:
        match(out, read("results.json"))
        print("1597: helium closure, signed charged-dipole intervals and inverse IO targets verified")
    else:
        print(json.dumps(out, indent=2))
