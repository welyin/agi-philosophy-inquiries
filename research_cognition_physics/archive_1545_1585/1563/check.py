"""Replay the kaon shape and rate envelopes. No experimental likelihood is fitted."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def endpoints(fn, *ranges):
    vals = [fn(*point) for point in itertools.product(*ranges)]
    return [min(vals), max(vals)]


def compute(data):
    def box(key):
        value, error = data[key]
        w = data["width_multiplier"] * error
        return [value - w, value + w]

    def val(key):
        return data[key][0]

    divide = lambda a, b: a / b
    B = endpoints(divide, box("FK_over_Fpi_isosymmetric"), box("fplus_neutral"))
    Bmid = val("FK_over_Fpi_isosymmetric") / val("fplus_neutral")
    Z = endpoints(lambda a, u, p: a * u / p, box("inclusive_axial_product"),
                  box("nuclear_vector_product"), box("neutral_electron_vector_product"))
    right = endpoints(lambda z, b: 0.5 * math.log(z / b), Z, B)
    iso = box("charged_to_isosymmetric_F_squared_minus_one")
    c0 = endpoints(lambda b, i, d: math.log(b * math.sqrt(1 + i) + d),
                   B, iso, data["neutral_CT_additive_work_range"])
    c0mid = math.log(Bmid * math.sqrt(1 + val("charged_to_isosymmetric_F_squared_minus_one")))
    cplus = endpoints(lambda x, y: x + y, c0, data["charged_minus_neutral_QCD_logC_work_range"])
    kin = data["kinematic_GeV"]
    m_mu = kin["m_mu"]
    dk0 = kin["M_K0"] ** 2 - kin["M_piplus"] ** 2
    dkp = kin["M_Kplus"] ** 2 - kin["M_pi0"] ** 2
    m_diff = endpoints(lambda s, u: s - u, box("working_ms_msbar_2GeV_GeV"),
                       box("working_mu_msbar_2GeV_GeV"))
    common = m_mu * (val("working_ms_msbar_2GeV_GeV") - val("working_mu_msbar_2GeV_GeV"))
    eps0 = data["benchmark_real_epsilon"]
    fit = data["work_abs_fitted_logC_residual"]
    neutral = box("shape_KTeV_neutral")
    charged = box("shape_NA48_2_charged")
    eps = endpoints(lambda obs, pred, rem, m: (obs - pred - rem) * m_mu * m / dk0,
                    neutral, c0, [-fit, fit], m_diff)
    axial = 0.2754
    electron = axial * val("nuclear_vector_product") / Bmid
    charge_shift = -0.010
    assert box("inclusive_axial_product")[0] <= axial <= box("inclusive_axial_product")[1]
    assert box("neutral_electron_vector_product")[0] <= electron <= box("neutral_electron_vector_product")[1]
    assert neutral[0] <= c0mid <= neutral[1]
    assert charged[0] <= c0mid + charge_shift <= charged[1]
    tmax = max((kin["M_K0"] - kin["M_piplus"]) ** 2, (kin["M_Kplus"] - kin["M_pi0"]) ** 2)
    zmax = eps0 * tmax / (m_mu * min(m_diff))
    log_diff = endpoints(lambda a, b: a - b, neutral, charged)
    diff_slope = (dk0 - dkp) / common
    return {
        "isospin_symmetric_B": B, "rate_Z": Z, "B_center": Bmid,
        "half_log_Z_over_B": right,
        "neutral_QCD_logC_center_at_zero_remainder": c0mid,
        "neutral_QCD_logC": c0, "charged_QCD_logC": cplus,
        "neutral_shape": neutral, "charged_shape": charged,
        "measured_neutral_minus_charged_logC": log_diff,
        "epsilon_neutral_work_interval": eps,
        "delta_K0pi_and_Kpluspi_GeV2": [dk0, dkp],
        "scalar_slopes_neutral_and_charged": [dk0 / common, dkp / common],
        "scalar_slope_difference": diff_slope,
        "logC_difference_lever_at_0p001": eps0 * diff_slope,
        "worst_abs_physical_log_shape_remainder_at_0p001": zmax * zmax / (2 * (1 - zmax)),
        "joint_zero_scalar_standard_weak_witness": {"axial": axial, "KLe3": electron,
            "neutral_logC": c0mid, "charged_logC": c0mid + charge_shift},
        "same_rate_double_ratio_right_current_witness": {"xi_d": 0.1, "xi_s": 0.1,
            "common_vector": 1, "common_axial": 0.9 / 1.1}
    }


def rounded(obj):
    if isinstance(obj, float):
        return float(format(obj, ".12g"))
    if isinstance(obj, list):
        return [rounded(x) for x in obj]
    if isinstance(obj, dict):
        return {k: rounded(v) for k, v in obj.items()}
    return obj


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write-new", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    out = rounded(compute(json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))))
    target = HERE / "results.json"
    if args.write_new:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(out, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
    elif json.loads(target.read_text(encoding="utf-8")) != out:
        raise SystemExit("Saved certificate differs")
    print(json.dumps(out, ensure_ascii=True))
