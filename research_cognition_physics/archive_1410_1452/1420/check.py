"""Solar production weights under a signed adiabatic matter potential, round 1420."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def compute():
    inputs = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
    reactor = json.loads((HERE / inputs["reactor_inputs"]).read_text(encoding="utf-8"))
    sno = inputs["sno"]
    juno = reactor["juno_primary"]
    daya = reactor["daya_bay"]
    bands = []
    for widths in (0, *inputs["input_envelopes"]):
        p_max = sno["P"] + widths * (sno["stat"] + sno["syst"])
        solar_sin2_min = juno["sin_sq_theta12"] - widths * juno["one_sigma"]
        c_max = 1 - 2 * solar_sin2_min
        d_max = daya["sin_sq_2theta13"] + widths * daya["one_sigma"]
        z_max = (1 - math.sqrt(1 - d_max)) / 2
        flat_floor = z_max**2 + (1 - z_max)**2 / 2
        k_min = 2 * (flat_floor - p_max) / (1 - z_max)**2
        assert 0 < k_min and k_min * (1 + 2 * c_max) > c_max**2
        q_res = k_min / c_max
        q_align = (c_max + q_res) / (1 + c_max)
        c_min = 1 - 2 * (juno["sin_sq_theta12"] + widths * juno["one_sigma"])
        opposite_floor = flat_floor + (1 - z_max)**2 * c_min**2 / 2
        bands.append(dict(widths=widths, sno_upper=p_max, threshold_floor=flat_floor,
                          mismatch_margin_without_theta12=flat_floor - p_max,
                          opposite_floor=opposite_floor, k_min=k_min, aligned_floor=q_align,
                          beyond_resonance_floor=q_res, c_max=c_max))
    central = bands[0]
    c = central["c_max"]
    w = -central["beyond_resonance_floor"]
    effective_a = c - math.sqrt(1 - c*c) * w / math.sqrt(1-w*w)
    assert abs((c-effective_a)/math.hypot(c-effective_a, math.sqrt(1-c*c)) - w) < 1e-12

    # Same recoil distributions for two opposite signed total electron amplitudes.
    def kernel(left, right, y, mass_over_energy):
        return left**2 + right**2 * (1-y)**2 - left*right*mass_over_energy*y
    max_diff = 0.
    for left, y, ratio, anti in itertools.product((.73, -.27), (.1, .3, .7), (.02, .05), (False, True)):
        right = .23
        if anti:
            left, right = right, left
        max_diff = max(max_diff, abs(kernel(left, right, y, ratio) - kernel(-left, -right, y, ratio)))
    assert max_diff < 1e-14
    return {"round":1420, "bands":bands, "single_density_example_signed_A":effective_a,
            "electron_recoil_mirror_max_error":max_diff}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = HERE / "results.json"
    result = compute()
    if args.write:
        with output.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        expected = json.loads(output.read_text(encoding="utf-8"))
        def rounded(obj):
            if isinstance(obj, float):
                return round(obj, 11)
            if isinstance(obj, list):
                return [rounded(v) for v in obj]
            if isinstance(obj, dict):
                return {k: rounded(v) for k, v in obj.items()}
            return obj
        assert rounded(expected) == rounded(result)
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
