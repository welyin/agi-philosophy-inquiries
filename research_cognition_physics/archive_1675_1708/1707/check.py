"""Tree-calibrated CKM roots, two neutral-B phases and a phase-uniform bound."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


CFG = read("inputs.json")
MODEL = runpy.run_path(str(HERE / CFG["shooting_program"]))
TREE = read(CFG["semileptonic_moduli"])["direct_semileptonic_moduli"]


def error(record):
    return math.hypot(*record[1:])


def band(record, n):
    return [record[0] - n * error(record), record[0] + n * error(record)]


def phase_data(v):
    cd, cs, cb = v[1]
    td, ts, tb = v[2]
    beta = np.angle(-cd * cb.conjugate() / (td * tb.conjugate()))
    beta_s = np.angle(-ts * tb.conjugate() / (cs * cb.conjugate()))
    ratio = (cd * ts) / (td * cs)
    np.testing.assert_allclose(np.exp(1j * (beta + beta_s)), ratio / abs(ratio), atol=2e-15)
    np.testing.assert_allclose(abs(1 - 1 / ratio), abs(v[0, 2]) / abs(cd * ts), atol=2e-15)
    return beta, beta_s


def report_point(gamma, theta):
    v = MODEL["low_ckm"](gamma, theta)
    beta, beta_s = phase_data(v)
    phases = phase_data(np.exp(1j * np.array([1.31, -.07, .41]))[:, None]
                        * v * np.exp(1j * np.array([-.91, -.54, 2.41])))
    np.testing.assert_allclose(phases, [beta, beta_s], atol=2e-15)
    golden, phi = math.sin(2 * beta), -2 * beta_s
    def measured_width(value, record):
        return (value - record[0]) / error(record)
    dist = {"S_d": measured_width(golden, CFG["S_d"])}
    phase_only = {}
    central_decay_phase = math.asin(CFG["S_d"][0] / math.sqrt(1 - CFG["C_d"][0]**2))
    for key, source in CFG["phi_s_rad"].items():
        dist[key] = measured_width(phi, source)
        phase_only[key] = [central_decay_phase - 2 * beta, source[0] - phi]
    return {"gamma_deg": float(gamma), "moduli_us_ub_cb": abs(v[[0, 0, 1], [1, 2, 2]]).tolist(),
            "S_d": golden, "phi_s_rad": phi, "Vtd_over_Vts": abs(v[2, 0] / v[2, 1]),
            "gamma_in_time_integrated": {k: bool(limits[0] <= gamma <= limits[1])
                for k, limits in CFG["time_integrated_gamma_deg"].items() if k != "central"},
            "signed_marginal_widths": dist,
            "separate_offsets_for_central_phase_pair_rad": phase_only}


def bound_on_sine(r_hi, n):
    a_lo, a_hi = band(TREE["Vus"], n)
    _, b_hi = band(TREE["Vcb"], n)
    ub_hi = r_hi * b_hi
    c13_lo = math.sqrt(1 - ub_hi**2)
    s12_hi, s23_hi = a_hi / c13_lo, b_hi / c13_lo
    cd_lo = a_lo * math.sqrt(1 - s23_hi**2) - s23_hi * ub_hi
    ts_over_b_lo = math.sqrt(1 - s12_hi**2) - s12_hi * r_hi
    return r_hi / (cd_lo * ts_over_b_lo)


def phase_uniform_bound(name, ratio_record, n):
    r, down, up = ratio_record
    lo, hi = r - n * math.hypot(*down), r + n * math.hypot(*up)
    upper_sine = bound_on_sine(hi, n)
    assert 0 < upper_sine < 1
    angle_limit = 2 * math.asin(upper_sine)
    theta_min = math.asin(band(CFG["S_d"], n)[0])
    out = {}
    for port, measurement in CFG["phi_s_rad"].items():
        phase_min, phase_max = band(measurement, n)
        # Minimum circular separation from the complete allowed Bd arc to the Bs arc.
        gap = min(theta_min - phase_max, math.pi + theta_min + phase_min)
        assert gap > 0
        r_lower, r_upper = .01, .15
        for _ in range(55):
            mid = (r_lower + r_upper) / 2
            if bound_on_sine(mid, n) < math.sin(gap / 2):
                r_lower = mid
            else:
                r_upper = mid
        out[port] = {"minimum_data_arc_separation_rad": gap,
                     "differential_offset_lower_bound_rad": max(0., gap - angle_limit),
                     "necessary_r_floor_for_common_offset": r_lower}
    return {"name": name, "ratio_interval": [lo, hi], "phase_difference_upper_rad": angle_limit,
            "max_sine_beta_plus_beta_s": upper_sine, "phase_ports": out}


def main():
    frozen = read(CFG["charged_roots"])
    rows = {}
    for kind, branch in frozen["charged_endpoints"].items():
        for k, root in branch.items():
            rows[f"{kind}_k{k}"] = report_point(root["gamma_MZ_deg"], MODEL["angles"](kind))
    for label, root in frozen["down_one_width_lower"].items():
        rows["lower_down_" + label] = report_point(root["gamma_MZ_deg"], MODEL["angles"]("direct"))
    for label, root in frozen["independent_Vub_sensitivity"].items():
        n = float(label.split("_")[-1])
        rows["direct_" + label] = report_point(root["gamma_MZ_deg"], MODEL["angles"]("direct", n))
    new = {}
    for name, (r, _, _) in CFG["ratio_ports"].items():
        dz = (r * TREE["Vcb"][0] - TREE["Vub"][0]) / error(TREE["Vub"])
        theta = MODEL["angles"]("direct", dz)
        for k in (1, 2):
            root = MODEL["solve"](k, "direct", vub_shift=dz)
            row = report_point(root["gamma_MZ_deg"], theta)
            ub, cb = row["moduli_us_ub_cb"][1:]
            assert abs(ub / cb - r) < 1e-15
            row.update({key: root[key] for key in ("mass_invariant_residual",
                       "max_charged_spectrum_relative_error", "large_modulus", "frame_step_norm")})
            new[f"{name}_k{k}"] = row
    outer = {name: {str(n): phase_uniform_bound(name, rec, n) for n in (1, 2)}
             for name, rec in CFG["ratio_ports"].items()}
    assert outer["Bs_low_original"]["2"]["phase_ports"]["CMS_8plus13"]["differential_offset_lower_bound_rad"] > .12
    assert outer["Bs_low_original"]["2"]["phase_ports"]["LHCb"]["differential_offset_lower_bound_rad"] > .09
    return {"R1702": rows, "new_charged_roots": new, "phase_uniform_certificate": outer}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = main()
    out = HERE / "results.json"
    if args.check:
        assert json.loads(out.read_bytes()) == result
        print("1707 verified")
    else:
        out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
