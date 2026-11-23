"""Compare axial kaon and pion rates with the three-pass flavour trajectory."""
import argparse
import json
import math
from pathlib import Path
import runpy


HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


# Purely leptonic radiation-inclusive K/pi axial ratios. Both analyses use
# overlapping decay rates. The vector alternative includes pion beta.
RATIOS = {
    "RM123_axial": (0.23135, 0.00046),  # arXiv:1904.08731, eq. 108
    "RBC_axial": (0.23154, math.hypot(.00028,.00015,.00045,.00065)),  # 2211.12865 eq. 6.11
    "Seng_vector": (0.22908, 0.00087),  # arXiv:2107.14708
}
MODULE = runpy.run_path(str(HERE / "../1732/check.py"))
MODULI = read("../1702/inputs.json")["direct_semileptonic_moduli"]
POINT = read("../1732/inputs.json")["exploratory_weak_offsets"]["registered_with_lower_us"]
ETA = read("../1730/results.json")["max_all_light_singlet_probability_at_Higgs_budget"]


def measured(key, shift=0.0):
    center, *errors = MODULI[key]
    return center + shift * math.hypot(*errors)


def required_loss(strange, beauty, ratio):
    return 1.0 - beauty**2 - strange**2 * (1.0 + 1.0/ratio**2)


def box(ratio, error, halfwidth):
    return required_loss(measured("Vus",halfwidth), measured("Vub",halfwidth),
                         ratio - error*halfwidth)


def meeting_width(center, error):
    lo, hi = 0.0, 6.0
    if box(center,error,0) <= ETA:
        return 0.0
    assert box(center,error,hi) < ETA
    for _ in range(55):
        mid = (lo+hi)/2
        if box(center,error,mid) > ETA:
            lo = mid
        else:
            hi = mid
    return hi


def compute():
    beauty = measured("Vub",POINT["ub"])
    strange = measured("Vus",POINT["us"])
    down = math.sqrt(1 - strange**2 - beauty**2)
    original_ratio = strange/down
    outputs = {}
    for name,(ratio,error) in RATIOS.items():
        # With three unitary light flavours, the rate ratio and independent
        # semileptonic beauty modulus set strange; no nuclear beta input.
        new_strange = ratio * math.sqrt((1-beauty**2)/(1+ratio**2))
        weak = dict(POINT,us=(new_strange-measured("Vus"))/math.hypot(*MODULI["Vus"][1:]))
        MODULE["WEAK"][name] = weak
        row = MODULE["solve"](3,name,MODULE["CONFIG"]["log_step"])
        neutral = MODULE["neutral_tests"]({name:row})[name]
        loss = required_loss(strange,beauty,ratio)
        # Equality is attainable in the uniform-projection envelope. The
        # specific mass model is allowed far less attenuation.
        if 0 < loss < 1:
            undiluted = [strange/ratio, strange, beauty]
            assert abs(sum((v/math.sqrt(1-loss))**2 for v in undiluted)-1) < 2e-14
        outputs[name] = {
            "central_ratio_and_marginal_error":[ratio,error],
            "ratio_residual_original_in_reported_widths":(original_ratio-ratio)/error,
            "required_rate_rescaling_original":(ratio/original_ratio)**2-1,
            "lower_bound_on_uniform_loss_original": loss,
            "lower_bound_on_uniform_loss_over_one_width_box":box(ratio,error,1),
            "minimum_common_marginal_halfwidth":meeting_width(ratio,error),
            "unitary_recalibration_from_ratio":{"Vus":new_strange,"Vus_in_original_widths":weak["us"],
                "mass_widths_d_s_b": row["down_strange_bottom_marginal_widths"],
                "gamma_MZ_deg":row["gamma_MZ_deg"],"neutral_fixed_gain":neutral}
        }
    assert all(outputs[k]["lower_bound_on_uniform_loss_over_one_width_box"]>ETA for k in ("RM123_axial","RBC_axial"))
    assert abs(outputs["Seng_vector"]["ratio_residual_original_in_reported_widths"]) < 1
    return MODULE["PARENT"]["tidy"]({"weak_Vus":strange,"weak_Vub":beauty,"weak_Vud":down,
                                    "weak_Vus_over_Vud":original_ratio,"curvature_loss_budget":ETA,
                                    "experimental_comparisons":outputs})


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    args = cli.parse_args()
    result = compute()
    path = HERE / "independent_results.json"
    if args.check:
        assert result == json.loads(path.read_bytes())
        print("1734 verified axial-rate boxes, attainable unitary bound and three full recalibrated trajectories")
    else:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n",encoding="utf8")
