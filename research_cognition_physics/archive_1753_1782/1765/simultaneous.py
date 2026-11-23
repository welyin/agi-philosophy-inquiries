"""Complementary charged calibration and a simultaneous marginal-box witness."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
BASE = runpy.run_path(str(HERE / "check.py"))
STRONG = BASE["STRONG"]
GAMMA = BASE["read"](BASE["CFG"]["secondary_tree_angle"])["time_integrated_gamma_deg"]
CLOCK = BASE["CLOCKS"]["delta_ms"]
CS = BASE["CFG"]["mixing_coefficient_ps_inv"]


def pull(value, data):
    return (value - data[0]) / math.hypot(*data[1:])


def inverse_clock(row):
    return CLOCK[0] / row["tbs_squared"]


def project(row, coefficient):
    v = row["moduli"]
    ud = math.sqrt(1-v["Vus"]**2-v["Vub"]**2)
    g = row["gamma"]
    gw = (GAMMA["central"]-GAMMA["68"][0] if g < GAMMA["central"]
          else GAMMA["68"][1]-GAMMA["central"])
    other = BASE["compare"](row)
    readings = {
        "Kl3": other["Kl3_us_pull"], "golden_S": other["golden_S_pull"],
        "Bs_clock": pull(coefficient*row["tbs_squared"], CLOCK),
        "Cs": pull(coefficient, CS),
        "QCDSF": pull(1/row["ds"], STRONG["QCDSF_s_d"]),
        "gamma": (g-GAMMA["central"])/gw,
        "inclusive_ub": other["inclusive_B_pulls_ub_cb"][0],
        "inclusive_cb": other["inclusive_B_pulls_ub_cb"][1],
    }
    return {
        "weak": v, "sd": 1/row["ds"], "gamma": g, "golden_S": row["sin2beta"],
        "Cs": coefficient, "clock": coefficient*row["tbs_squared"],
        "pulls_in_selected_box": readings,
        "largest_selected_marginal_pull": max(map(abs, readings.values())),
        "heldout_exclusive_ub_cb_pulls": other["exclusive_B_pulls_ub_cb"],
        "ud": ud, "ud_pull": pull(ud, STRONG["holdout_nuclear"]),
        "us_over_ud": v["Vus"]/ud, "axis_pull": pull(v["Vus"]/ud, STRONG["holdout_axial"]),
    }


def central_residual(row):
    return np.log([
        row["moduli"]["Vus"]/STRONG["holdout_Kl3"][0],
        row["sin2beta"]/BASE["CLOCKS"]["golden_S"][0],
        row["tbs_squared"]*CS[0]/CLOCK[0],
    ])


def calculate():
    x = np.log([.224, .0038, .0425])
    def solve_at(trial):
        row = BASE["evaluate"](trial, .25)
        return row, central_residual(row)
    row, res = solve_at(x)
    step = 1e-4
    jac = np.column_stack([(solve_at(x + d)[1]-res)/step for d in step*np.eye(3)])
    for _ in range(10):
        if max(abs(res)) < 1e-10:
            break
        dx = np.linalg.solve(jac, -res)
        x += dx
        row, new = solve_at(x)
        jac += np.outer(new-res-jac@dx,dx)/(dx@dx)
        res = new
    assert max(abs(res)) < 1e-10
    central = project(row, CS[0])
    # Chosen after viewing all declared box observables. No statistical combination.
    trial = np.log([.223875, .003835, .04264])
    witness = BASE["evaluate"](trial, .25)
    fine = BASE["evaluate"](trial, .125)
    fixed = inverse_clock(witness)
    joint = project(witness, fixed)
    fine_joint = project(fine, fixed)
    assert central["largest_selected_marginal_pull"] > 1
    assert max(joint["largest_selected_marginal_pull"], fine_joint["largest_selected_marginal_pull"]) < 1
    kl = STRONG["holdout_Kl3"]
    ub = STRONG["inclusive_ub"]
    us_max, ub_max = kl[0]+math.hypot(*kl[1:]), ub[0]+math.hypot(*ub[1:])
    ud_min = math.sqrt(1-us_max**2-ub_max**2)
    ratio_max = us_max/ud_min
    geometry = {
        "us_max": us_max, "ub_max": ub_max, "unitary_ud_min": ud_min,
        "nuclear_ud_max": sum(STRONG["holdout_nuclear"]),
        "unitary_us_over_ud_max": ratio_max,
        "axial_us_over_ud_min": STRONG["holdout_axial"][0]-STRONG["holdout_axial"][1],
    }
    assert ud_min > geometry["nuclear_ud_max"] and ratio_max < geometry["axial_us_over_ud_min"]
    return {"central": central, "simultaneous": joint,
            "fine_witness_max_marginal_pull": fine_joint["largest_selected_marginal_pull"],
            "first_row_incompatibility_at_one_width": geometry}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = BASE["rounded"](calculate())
    path = HERE / "simultaneous_results.json"
    if args.check:
        if result != json.loads(path.read_bytes()):
            raise SystemExit("saved simultaneous comparison differs")
        print("1765 simultaneous calibration and first-row bounds verified")
    else:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
