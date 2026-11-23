"""Calibrated fourth moment, positive inclusive width, and Euclidean bubble curvature."""
import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def read(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))

def bubble_kernel(u, mass_squared, step_squared):
    a = u * step_squared
    return math.log1p(a*a / ((mass_squared + a)*(mass_squared + 3*a)))

def bubble_interval(mass, scale, v, segments):
    # The integrand increases on [0, 1/2] and is symmetric.
    values = []
    for i in range(segments + 1):
        x = i / (2 * segments)
        values.append(bubble_kernel(x*(1-x), mass*mass, scale*scale))
    factor = mass**4 / (8 * math.pi**2 * v*v * segments)
    return [factor * math.fsum(values[:-1]), factor * math.fsum(values[1:])]

def compute():
    cfg = read("inputs.json")
    archive = read(cfg["radial_result"])["bosonic_fourth_root_with_f_equal_v_no_neutral"]
    base = read(cfg["radial_inputs"])
    mass = [archive[k]["Higgs_matching_boundary_GeV"] for k in ("box_low","center","box_high")]
    v, dv = base["matching"]["v_GeV"]
    h, dh = base["atlas"]["mass_GeV"]
    br = cfg["invisible_upper"] / cfg["production_times_acceptance_floor"]
    width_limit = cfg["visible_width_ceiling_GeV"] * br / (1 - br)
    weighted_beta = (8 * math.pi * (v+dv)**2 * (h+dh) * width_limit
                     / (mass[0]**4 * cfg["relative_invisible_width_floor"]))
    light_limits = {str(x): weighted_beta / math.sqrt(1-x*x)
                    for x in cfg["visible_mass_fraction_of_threshold"]}
    rows = {}
    for scale in cfg["euclidean_step_GeV"]:
        minimum = bubble_interval(mass[0],scale,v+dv,cfg["segments"])
        one = bubble_interval(mass[1],scale,v,cfg["segments"])
        four = [4*y for y in bubble_interval(mass[1]/math.sqrt(2),scale,v,cfg["segments"])]
        assert four[0] > one[1] > minimum[0]
        rows[str(scale)] = {"floor_GeV2":minimum[0],
                            "one_central_GeV2":one,"four_central_GeV2":four,
                            "floor_over_step_squared":minimum[0]/scale**2}
    return {"fourth_root_GeV":mass,"width_upper_GeV":width_limit,
            "fourth_moment_weighted_open_velocity_upper":weighted_beta,
            "fraction_upper_for_masses_below_x_times_threshold":light_limits,
            "euclidean_second_differences":rows}

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    options = p.parse_args()
    data = compute()
    if options.write_new:
        with (ROOT/"results.json").open("x",encoding="utf-8") as dest:
            dest.write(json.dumps(data,indent=2,allow_nan=False)+"\n")
    else:
        assert data == read("results.json")
        assert 0 < data["fourth_moment_weighted_open_velocity_upper"] < 2e-6
        print("1570 inclusive invisible limit and sharp one-singlet Euclidean envelope checked")
