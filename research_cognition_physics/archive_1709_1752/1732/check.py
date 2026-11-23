"""Solve all three down masses using genuine one-loop shooting for each port length."""
import argparse
from functools import lru_cache
import json
import math
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
CONFIG = json.loads((HERE / "inputs.json").read_bytes())
PARENT = runpy.run_path(str(HERE / CONFIG["parent"]))
DOWN = PARENT["DOWN"]
REGISTERED = {name: PARENT["PRIOR"][name] for name in ("us", "ub", "cb")}
WEAK = {"central": {}, "registered": REGISTERED, **CONFIG["exploratory_weak_offsets"]}


@lru_cache(maxsize=64)
def trial(strange, weak_set, step):
    # R1727 accepts standard-width offsets. Start and solve in absolute Yukawa units.
    offset = (strange - DOWN[1, 0]) / DOWN[1, 1]
    weak = WEAK[weak_set]
    return PARENT["exchange_tau_for_bottom"](
        {**weak, "strange": offset}, True, CONFIG["seed_gamma_deg"], step
    )


def solve(k, weak_set, step):
    def evaluate(strange):
        row = trial(strange, weak_set, step)
        return math.log(row["x"]) - k * math.log(row["abs_z"]), row

    lo, hi = CONFIG["strange_absolute_bracket"]
    flo, _ = evaluate(lo)
    fhi, _ = evaluate(hi)
    assert flo > 0 > fhi
    for _ in range(32):
        strange = (lo * fhi - hi * flo) / (fhi - flo)
        residual, row = evaluate(strange)
        if abs(residual) < 2e-10:
            break
        if residual > 0:
            lo, flo = strange, residual
        else:
            hi, fhi = strange, residual
    else:
        raise ValueError("port length did not close")
    c = row["Vud_10PeV"]
    predicted = row["MZ_down_yukawas"]
    ratio = math.prod(row["high_e_mu"]) * row["tau_10PeV"] / math.prod(row["high_down_yukawas"])
    assert abs(ratio / c**(2*k + 1) - 1) < 2e-8
    assert abs(row["return_relative"]) < 2e-9
    assert abs(predicted[1] / strange - 1) < 1e-12
    assert abs(row["tau_MZ"] / PARENT["LEPTONS"][2, 0] - 1) < 2e-10
    widths = (predicted - DOWN[:,0]) / DOWN[:,1]
    gamma_status = PARENT["observed"](row)
    return {
        "k": k, "positive_x": row["x"], "Vud_10PeV": c, "high_mass_product_ratio": ratio,
        "down_strange_bottom_MZ": predicted, "down_strange_bottom_marginal_widths": widths.tolist(),
        "gamma_MZ_deg": row["gamma_MZ_deg"], "gamma_in_time_integrated_68": gamma_status["in_time_integrated_gamma_68"],
        "gamma_in_time_integrated_95": gamma_status["in_time_integrated_gamma_95"], "equation_residual": residual
    }


def neutral_tests(solutions):
    neutral = runpy.run_path(str(HERE / CONFIG["neutral_transport"]))
    combined = {}
    for name, charged in solutions.items():
        wd, ws, wb = charged["down_strange_bottom_marginal_widths"]
        weak = {key: WEAK[name].get(key, 0.0) for key in ("us", "ub", "cb")}
        pt = neutral["endpoint"]("explicit_open", CONFIG["log_step"],
                                 explicit_offsets={**weak, "strange":ws, "bottom":wb})
        assert abs(pt["row"]["gamma_MZ_deg"] - charged["gamma_MZ_deg"]) < 1e-7
        assert abs(pt["row"]["x"] - charged["positive_x"]) < 1e-8
        obs = neutral["observables"](pt["C"], pt["b"], neutral["FROZEN"][0])
        combined[name] = {
            "up_monomials": neutral["CFG"]["up_monomials"], "frozen_raw_gain": obs["rho"],
            "charged_mass_widths": [wd, ws, wb], "charged_gamma_in_68": charged["gamma_in_time_integrated_68"],
            "neutral_angles": obs["angle_readings"].tolist(), "angles_in_box": obs["angles_in_box"],
            "dm31_and_32_eV2": [obs["dm31_eV2"],obs["dm32_eV2"]],
            "dm31_and_32_in_box": obs["dm31_and_32_in_box_at_central_solar"],
            "delta_deg": obs["delta_deg"], "SK_T2K_FC68_90": obs["SK_T2K_FC68_90"],
            "mbb_meV": obs["mbb_meV"], "mbeta_meV": obs["mbeta_meV"]
        }
    return combined


def compute():
    step = CONFIG["log_step"]
    results = {name: [solve(k, name, step) for k in CONFIG["path_lengths"]]
               for name in ("central", "registered")}
    finer = {}
    for name, rows in results.items():
        row = next(r for r in rows if r["k"] == 3)
        half = solve(3, name, step/2)
        width_change = max(abs(a-b) for a,b in zip(row["down_strange_bottom_marginal_widths"],
                                                  half["down_strange_bottom_marginal_widths"]))
        gamma_change = abs(half["gamma_MZ_deg"] - row["gamma_MZ_deg"])
        assert width_change < 1e-5 and gamma_change < 1e-5
        finer[name] = {"max_width_change":width_change, "gamma_deg_change":gamma_change}
    k3 = {name:next(r for r in rows if r["k"]==3) for name,rows in results.items()}
    exploratory = {name:solve(3, name, step) for name in CONFIG["exploratory_weak_offsets"]}
    neutral = neutral_tests({**k3, **exploratory})
    witness = neutral["registered_with_lower_us"]
    assert max(abs(v) for v in witness["charged_mass_widths"]) < 1 and witness["charged_gamma_in_68"]
    assert all(witness["angles_in_box"]) and all(witness["dm31_and_32_in_box"])
    assert list(witness["SK_T2K_FC68_90"]) == [False, True]
    assert not all(neutral["charged_interior"]["dm31_and_32_in_box"])
    return PARENT["tidy"]({"solutions": results, "half_step_k3": finer, "exploratory_k3":exploratory,
                           "same_point_neutral_fixed_gain_solar_only":neutral})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    computed = compute()
    path = HERE / "results.json"
    if args.check:
        assert computed == json.loads(path.read_bytes())
        print("1732 checked three down masses, CKM phase, port-length competitors and determinant identity")
    else:
        with path.open("w", encoding="utf8") as stream:
            json.dump(computed, stream, indent=2, allow_nan=False)
            stream.write("\n")
