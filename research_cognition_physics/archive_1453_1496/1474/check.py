"""Replay the septet width, neutral sum and conditional ATLAS rate envelope."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def width_coefficient(mass, mw, vev):
    x = (mw / mass) ** 2
    return mass**2 / (32 * math.pi * vev**2) * math.sqrt(1 - 4*x) * (1 - 4*x + 12*x*x)


def mass_at_coefficient(target, mw, vev):
    low, high = 2 * mw, 10000.
    assert width_coefficient(high, mw, vev) > target > 0
    for _ in range(80):
        middle = (low + high) / 2
        if width_coefficient(middle, mw, vev) < target:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def certificate():
    data = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    study = data["ATLAS_139fb_single_SSWW"]
    limit = 2 * (study["reported_GM_sinH_upper_range"][1] + study["last_digit_outward_allowance"])**2
    width_budget = study["max_template_width_over_mass"]
    y = data["septet_y_witness"]
    x = math.sqrt(1 - y*y)
    mw, vev = data["W_mass_GeV"], data["v_GeV"]
    t = data["spin_T"]
    ladder = np.diag([math.sqrt((t-m)*(t+m+1)) for m in range(-t, t)], k=-1)
    neutral = np.eye(2*t+1)[:, t-data["hypercharge_Y"]]
    squared_ladder = float(np.linalg.norm(ladder @ ladder @ neutral)**2)
    z = squared_ladder * y*y / 8
    h = np.array([x, y/4, math.sqrt(15)*y/4])
    weak_couplings = np.array([x, 4*y, 0.])
    residual = float(np.linalg.norm(weak_couplings - (weak_couplings @ h)*h)**2)
    assert np.allclose([squared_ladder, h @ h, weak_couplings @ h, h[0]/x,
                        h[0]**2 + 16*h[1]**2, residual], [120., 1., 1., 1., 1., z])
    points = []
    for mass in data["M_points_known_tested_GeV"]:
        coeff = width_coefficient(mass, mw, vev)
        beta = math.sqrt(1 - 4*(mw/mass)**2)
        g_squared = (2*mw**2/vev)**2 * z
        explicit = beta*g_squared*(2 + (mass**2/(2*mw**2)-1)**2)/(32*math.pi*mass**2)
        assert math.isclose(explicit, coeff*z, rel_tol=1e-13)
        width_feasible = explicit <= width_budget
        branch_floor = explicit / width_budget
        item = dict(mass_GeV=mass, width_coefficient=coeff, partial_width_over_mass=explicit,
                    narrow_witness_feasible=width_feasible,
                    min_B_in_window=branch_floor if width_feasible else None,
                    min_zB_if_narrow=z * branch_floor if width_feasible else None,
                    max_transfer_for_witness_if_narrow=limit/(z*branch_floor) if width_feasible else None)
        for transfer in data["benchmark_non_signal_transfer"]:
            z_upper = min(width_budget/coeff, math.sqrt(width_budget * limit/(transfer*coeff)))
            item[f"y_upper_at_a_{transfer}"] = math.sqrt(z_upper/15)
        points.append(item)
    c300 = width_coefficient(300., mw, vev)
    example_total = .04
    example_rate = z*z*c300 / example_total
    assert z*c300 < example_total and math.isclose(math.sqrt(example_rate*example_total/c300), z)
    point = data["positive_potential_benchmark"]
    yf = point["y"]
    xf = math.sqrt(1-yf*yf)
    zf = 15*yf*yf
    hv = np.array([xf, yf/4, math.sqrt(15)*yf/4])
    slope = np.array([-1.5*yf/xf, 1., 0.])
    heavy, other = point["other_neutral_masses_GeV"]
    assert heavy == other
    mass, mh = point["charged_mass_GeV"], point["h_mass_GeV"]
    neutral_mass_squared = heavy**2*np.eye(3) + (mh**2-heavy**2)*np.outer(hv,hv)
    curvature_test = float(mass**2 * slope @ np.linalg.solve(neutral_mass_squared, slope))
    residual_min = float(np.linalg.eigvalsh(neutral_mass_squared - mass**2 * np.outer(slope,slope))[0])
    assert curvature_test < 1 and residual_min > 0 and (1.25*y*mass)**2 > mh**2
    cf = width_coefficient(mass, mw, vev)
    positive_benchmark = dict(y=yf, mass_GeV=mass, residual_neutral_strength=zf,
                              rank_one_positivity_test=curvature_test, smallest_remaining_curvature_GeV2=residual_min,
                              min_B=cf*zf/width_budget, min_zB=cf*zf*zf/width_budget,
                              max_transfer_for_witness=limit*width_budget/(cf*zf*zf),
                              needed_total_width_over_mass=cf*zf*zf/limit)
    return dict(conservative_search_upper_kappa_squared=limit, positive_potential_benchmark=positive_benchmark,
                y_if_BR_one_a_one=math.sqrt(limit/15),
                neutral_residual_squared=residual, doubly_charged_kappa_squared=z,
                double_ladder_norm_squared=squared_ladder, points=points,
                witness_mass_crossing_for_rate_envelope_GeV=mass_at_coefficient(limit*width_budget/z**2, mw, vev),
                witness_mass_crossing_for_width_budget_GeV=mass_at_coefficient(width_budget/z, mw, vev),
                total_width_over_mass_needed_at_300_a_one=c300*z*z/limit,
                synthetic_inversion=dict(total_width_over_mass=example_total, production_times_BR=example_rate,
                                         recovered_residual=math.sqrt(example_rate*example_total/c300)))


def compare(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(compare(v,b[k]) for k,v in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(compare(x,y) for x,y in zip(a,b))
    if isinstance(a, float):
        return isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-13)
    return a == b


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    mode = args.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    options = args.parse_args()
    result = certificate()
    path = HERE / "results.json"
    if options.write:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        assert compare(result,json.loads(path.read_text(encoding="utf-8")))
        print("1474: SU(2) ladder, neutral sum, five mass windows and width inversion verified.")
