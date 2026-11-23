"""Bound how much of the GW signed mass moment actual weakly charged neutrinos supply."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def invisible_loss(x):
    if x >= 0.25:
        return 1.0
    return 1.0 - (1.0 - x) * math.sqrt(1.0 - 4.0 * x)


def inverse_loss(loss):
    lo, hi = 0.0, 0.25
    for _ in range(70):
        mid = (lo + hi) / 2.0
        if invisible_loss(mid) < loss:
            lo = mid
        else:
            hi = mid
    return hi


def calculate(settings):
    base = read(settings["baseline_data"])
    old = read(settings["baseline_results"])
    v = sum(base["matching"]["v_GeV"])
    h = sum(base["atlas"]["mass_GeV"])
    mass_term = 8 * math.pi**2 * v**2 * h**2
    trace_term = -min(old["SM_trace_over_v4"]["marginal_box"]) * v**4
    assert math.isclose((mass_term + trace_term)**0.25,
                        old["bosonic_fourth_root_with_f_equal_v_no_neutral"]["box_high"]["Higgs_matching_boundary_GeV"])
    invisible = settings["lep_invisible"]
    loss = 3 - invisible["retained_lower"] + invisible["extra_yield_and_transfer_upper"]
    assert 0 < loss <= 0.1000000001
    x_max = inverse_loss(loss)
    z = settings["z_mass_upper_GeV"]
    transfer = settings["neutrino_running_over_pole_upper"]
    neutral = 4 * (transfer * z)**4 * x_max**2
    flavor_floor = settings["charged_flavor_floor_GeV"]
    named = mass_term + trace_term + neutral
    charged_cap = (named / 2)**0.25
    min_negative = 2 * flavor_floor**4 - named

    frontier = []
    for q in settings["charged_pole_over_running_tests"]:
        needed = 2 * (flavor_floor / q)**4
        for kappa in settings["radial_kappa_tests"]:
            available = mass_term/kappa**2 + trace_term + neutral
            frontier.append({"pole_over_running": q, "kappa_V": kappa,
                             "extra_negative_or_remainder_min_GeV4": max(0, needed - available)})
    return {
        "loss_budget": loss,
        "maximum_single_pole_GeV": z * math.sqrt(x_max),
        "active_neutral_signed_moment_cap_GeV4": neutral,
        "loose_neutral_cap_GeV4": 4 * (transfer*z)**4 * (loss/3)**2,
        "mass_term_GeV4": mass_term, "charged_SM_trace_deficit_cap_GeV4": trace_term,
        "radial_charged_mass_cap_GeV": charged_cap,
        "pole_over_running_min_if_only_registered": flavor_floor / charged_cap,
        "extra_negative_min_pure_doublet_GeV4": min_negative,
        "one_colorless_Dirac_mass_equivalent_GeV": (min_negative / 4)**0.25,
        "pure_vev_dilution_kappa_V_max": math.sqrt(mass_term / (2*flavor_floor**4 - trace_term - neutral)),
        "frontier": frontier
    }


def compare(saved, computed):
    if isinstance(computed, dict):
        assert saved.keys() == computed.keys()
        for key, value in computed.items():
            compare(saved[key], value)
    elif isinstance(computed, list):
        assert len(saved) == len(computed)
        for left, right in zip(saved, computed):
            compare(left, right)
    else:
        assert math.isclose(saved, computed, rel_tol=2e-12, abs_tol=1e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    options = parser.parse_args()
    results = calculate(read("inputs.json"))
    if options.write_new:
        with (HERE / "results.json").open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(results, indent=2, allow_nan=False) + "\n")
    elif options.check:
        compare(read("results.json"), results)
        print("1589: Higgs mass, Z width and charged flavor frontiers passed")
    else:
        print(json.dumps(results, indent=2))
