"""Calibrate the archived charged family without direct semileptonic CKM moduli."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


CFG = read("inputs.json")
STRONG = read(CFG["strong_and_inclusive"])
EXCLUSIVE = read(CFG["exclusive"])["direct_semileptonic_moduli"]
CLOCKS = read(CFG["measured_clocks"])["data"]
POINT = runpy.run_path(str(HERE / CFG["charged_model"]))["point"]
AMPLITUDES = runpy.run_path(str(HERE / CFG["incoming_mixing_amplitudes"]))["amplitudes"]


def evaluate(x, step):
    us, ub, cb = np.exp(x)
    row = POINT(HERE, us, ub, cb, step)
    moduli, gamma = row["low_moduli"], row["gamma_deg"]
    _, _, beta, norms = AMPLITUDES(moduli, gamma)
    return dict(ds=row["down_over_s"], moduli=moduli, gamma=gamma,
                tbs_squared=float(norms[1]**2), td_ts=float(norms[0]/norms[1]),
                beta_deg=math.degrees(beta), sin2beta=math.sin(2*beta), bs=row["bottom_over_s"])


def logs(row):
    return np.log([row["ds"], row["gamma"], row["tbs_squared"]])


def jacobian(x, step, base):
    delta = 1.e-4
    return np.column_stack([(logs(evaluate(x + d, step)) - base) / delta
                            for d in np.eye(3) * delta])


def solve(target, x, jac, step):
    row = evaluate(x, step)
    res = logs(row) - target
    for _ in range(12):
        if np.max(np.abs(res)) < 5.e-10:
            break
        change = np.linalg.solve(jac, -res)
        x = x + change
        row = evaluate(x, step)
        new_res = logs(row) - target
        jac = jac + np.outer(new_res - res - jac @ change, change) / (change @ change)
        res = new_res
    assert max(abs(res)) < 5.e-10
    return row, x


def compare(row):
    def pull(value, measured):
        return (value - measured[0]) / math.hypot(*measured[1:])
    v = row["moduli"]
    rho_xi = CLOCKS["HPQCD"]
    md, ms = CLOCKS["delta_md"], CLOCKS["delta_ms"]
    mbd, mbs = CLOCKS["mass_d"], CLOCKS["mass_s"]
    ratio = rho_xi[0] * math.sqrt(md[0]/ms[0] * mbs[0]/mbd[0])
    return {
        "inclusive_B_pulls_ub_cb": [pull(v["Vub"], STRONG["inclusive_ub"]),
                                     pull(v["Vcb"], STRONG["inclusive_cb"])],
        "exclusive_B_pulls_ub_cb": [pull(v["Vub"], EXCLUSIVE["Vub"]),
                                     pull(v["Vcb"], EXCLUSIVE["Vcb"])],
        "Kl3_us_pull": pull(v["Vus"], STRONG["holdout_Kl3"]),
        "golden_S_pull": pull(row["sin2beta"], CLOCKS["golden_S"]),
        "td_ts_using_Bd_and_calibrating_Bs": ratio,
        "td_ts_difference": row["td_ts"] - ratio,
    }


def run():
    ds = 1 / STRONG["QCDSF_s_d"][0]
    gamma = CFG["time_integrated_gamma_deg"][0]
    clock = CLOCKS["delta_ms"][0]
    coeff, error = CFG["mixing_coefficient_ps_inv"]
    def target(g, c, mass_ratio=ds):
        return np.log([mass_ratio, g, clock/c])
    start = np.log([.225, .004, .042])
    step = .25
    base = logs(evaluate(start, step))
    jac = jacobian(start, step, base)
    central, x = solve(target(gamma, coeff), start, jac, step)
    local = jacobian(x, step, logs(central))
    rows = {"central": central}
    second = read(CFG["secondary_tree_angle"])["time_integrated_gamma_deg"]["central"]
    for name, g, c in [
        ("archived_time_integrated_angle", second, coeff),
        ("angle_plus_one", gamma + CFG["time_integrated_gamma_deg"][1], coeff),
        ("coefficient_plus_one", gamma, coeff + error),
        ("coefficient_minus_one", gamma, coeff - error),
    ]:
        rows[name], _ = solve(target(g, c), x, local, step)
    bmw = 2 / (STRONG["BMW_s_mud"][0] * (1 + STRONG["BMW_u_d"][0]))
    rows["BMW_strong_central"], _ = solve(target(gamma, coeff, bmw), x, local, step)
    fine = evaluate(x, step / 2)
    fine_delta = logs(fine) - target(gamma, coeff)
    assert max(abs(fine_delta)) < 2.e-6
    return {
        "rows": {name: {**row, "holdouts": compare(row)} for name, row in rows.items()},
        "fractional_observable_vs_trial_moduli_jacobian": local.tolist(),
        "fractional_inverse": np.linalg.inv(local).tolist(),
        "numerical_singular_values": np.linalg.svd(local, compute_uv=False).tolist(),
        "fine_log_residual_ds_gamma_squared_modulus": fine_delta.tolist(),
    }


def rounded(value):
    if isinstance(value, dict):
        return {k: rounded(v) for k, v in value.items()}
    if isinstance(value, list):
        return [rounded(v) for v in value]
    return float(f"{value:.9g}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = rounded(run())
    path = HERE / "results.json"
    if args.check:
        if result != read("results.json"):
            raise SystemExit("saved result differs")
        print("1765 charged-only solve verified")
    else:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
