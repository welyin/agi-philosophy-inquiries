"""One-loop elimination without a weak-angle input; marginal boxes, not likelihoods."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def interval(center, halfwidth):
    return center - halfwidth, center + halfwidth


def hull(values):
    values = list(values)
    return [min(values), max(values)]


def common_boxes(inp, widths):
    a, aerr = inp["A_em_MSbar"]
    alpha_s, serr = inp["alpha_s_lattice"]
    strong = interval(alpha_s, widths * serr)
    ref, ref_error = inp["alpha_s_reference_in_dispersion"]
    assert strong[0] >= ref - ref_error and strong[1] <= ref + ref_error
    available_A = interval(a, widths * aerr + inp["extra_A_for_alpha_s"] + inp["extra_A_for_scheme"])
    a3 = [1 / strong[1], 1 / strong[0]]
    bands = {name: interval(value, widths * error + inp["pole_transfer_allowance"])
             for name, (value, error) in inp["targets"].items()}
    joint = [max(x[0] for x in bands.values()), min(x[1] for x in bands.values())]
    assert joint[1] >= joint[0]
    return available_A, a3, bands, joint


def thresholds(inp, susy_mass):
    lower = inp["SM_beta_Y_2_3"]
    upper = inp["MSSM_beta_Y_2_3"] if susy_mass else lower
    low_time = math.log(susy_mass / inp["MZ_GeV"]) / (2 * math.pi) if susy_mass else 0.
    return lower, upper, low_time


def evaluate(inp, A, a3, k, susy_mass, shifts=(0., 0., 0.)):
    low, high, tlo = thresholds(inp, susy_mass)
    c, d = high[1] - high[2], high[0] - k * high[1]
    c0, d0 = low[1] - low[2], low[0] - k * low[1]
    ey, ew, ec = shifts  # additive at the low scale, physical hypercharge normalization
    xi = ey - k * ew - (d / c) * (ew - ec)
    denom = 1 + k + d / c
    s = (1 + (d / c) * a3 / A + (d * c0 / c - d0) * tlo / A - xi / A) / denom
    thi = (A * s - a3 - c0 * tlo - ew + ec) / c
    unified = a3 - low[2] * tlo - high[2] * thi - ec
    assert thi >= 0 and unified > 0
    actual = [k * unified, unified, unified]
    back = [x + b0 * tlo + b1 * thi + correction for x, b0, b1, correction in
            zip(actual, low, high, shifts)]
    assert max(abs(a - b) for a, b in zip(back, [A * (1 - s), A * s, a3])) < 3e-13
    return s, inp["MZ_GeV"] * math.exp(2 * math.pi * (tlo + thi)), unified


def k_from_target(inp, A, a3, s, susy_mass):
    low, high, tlo = thresholds(inp, susy_mass)
    thi = (A * s - a3 - (low[1] - low[2]) * tlo) / (high[1] - high[2])
    unified = a3 - low[2] * tlo - high[2] * thi
    assert thi > 0 and unified > 0
    return (A * (1 - s) - low[0] * tlo - high[0] * thi) / unified


def analyze_case(inp, k, mass, a, a3, targets, joint):
    predicted = [evaluate(inp, A, a_s, k, mass) for A, a_s in itertools.product(a, a3)]
    low, high, tlo = thresholds(inp, mass)
    slope = (high[0] - k * high[1]) / (high[1] - high[2])
    mismatch = []
    closest_witness = None
    for A, strong, target in itertools.product(a, a3, joint):
        s0 = evaluate(inp, A, strong, k, mass)[0]
        # A common offset is invisible; a shift on the weak factor alone saturates the span bound.
        delta_w = A * (1 + k + slope) * (target - s0) / (k + slope)
        reconstructed = evaluate(inp, A, strong, k, mass, (0., delta_w, 0.))[0]
        assert abs(reconstructed - target) < 2e-15
        mismatch.append(delta_w)
        if closest_witness is None or abs(delta_w) < abs(closest_witness[-1]):
            closest_witness = [A, strong, target, delta_w]
    missing = max(min(mismatch), -max(mismatch), 0.)
    if not missing:
        closest_witness = None
    ranges = list(map(hull, zip(*predicted)))
    def meets(band):
        return max(ranges[0][0], band[0]) <= min(ranges[0][1], band[1])
    return {"angle": ranges[0], "unification_GeV": ranges[1], "inverse_alpha_U": ranges[2],
            "individual_intersections": {name: meets(band) for name, band in targets.items()},
            "joint_intersection": meets(joint), "min_nonuniversal_inverse_spread": missing,
            "if_missing_one_weak_shift_A_a3_target_delta2": closest_witness}


def compute():
    inp = read("inputs.json")
    cases = inp["cases"]
    result = {"center": {}, "boxes": {}}
    for name, (n, d, mass) in cases.items():
        result["center"][name] = evaluate(inp, inp["A_em_MSbar"][0], 1 / inp["alpha_s_lattice"][0], n/d, mass)
    for widths in (1, 2):
        A, a3, bands, joint = common_boxes(inp, widths)
        result["boxes"][str(widths)] = {
            "A_em": A, "a3": a3, "separate_target_pole_angles": bands, "common_pole_angle": joint,
            "allowed_k": {}, "cases": {}}
        out = result["boxes"][str(widths)]
        ec, *e_errors = inp["original_E158_pole"]
        original_ee = interval(ec, widths * sum(e_errors) + inp["pole_transfer_allowance"])
        original_joint = [max(original_ee[0], bands["proton"][0]), min(original_ee[1], bands["proton"][1])]
        out["original_E158_with_proton_joint"] = original_joint
        out["original_E158_joint_compatible"] = {
            name: analyze_case(inp, n/d, mass, A, a3, bands, original_joint)["joint_intersection"]
            for name, (n, d, mass) in cases.items()}
        for spectrum, mass in (("SM", 0), ("MSSM_1TeV", 1000)):
            out["allowed_k"][spectrum] = hull(
                k_from_target(inp, ae, strong, angle, mass)
                for ae, strong, angle in itertools.product(A, a3, joint))
        for name, (n, d, mass) in cases.items():
            out["cases"][name] = analyze_case(inp, n/d, mass, A, a3, bands, joint)
    # With a free independent Y threshold, a target value fixes it for any positive k and scale.
    A, a3 = inp["A_em_MSbar"][0], 1 / inp["alpha_s_lattice"][0]
    s = .2315
    for k in (1.25, 5/3):
        eq_k = k_from_target(inp, A, a3, s, 0)
        unified = evaluate(inp, A, a3, eq_k, 0)[2]
        result.setdefault("same_two_targets_by_Y_threshold", []).append({
            "k": k, "delta_Y": (eq_k - k) * unified, "reconstructed_s":
            evaluate(inp, A, a3, k, 0, ((eq_k - k)*unified, 0., 0.))[0]})
    assert result["boxes"]["2"]["cases"]["orbifold_SM_32_25"]["joint_intersection"]
    assert not result["boxes"]["2"]["cases"]["canonical_SM"]["joint_intersection"]
    return result


def rounded(x):
    if isinstance(x, dict):
        return {k: rounded(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return list(map(rounded, x))
    return float(format(x, ".12g")) if isinstance(x, float) else x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--check", action="store_true")
    flags.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = rounded(compute())
    if args.write_new:
        with (HERE/"results.json").open("x", encoding="utf-8") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    elif read("results.json") != result:
        raise SystemExit("1568 gauge/target certificate changed")
    print("1568 reconstructed the three gauge couplings and both separate weak-angle boxes")
