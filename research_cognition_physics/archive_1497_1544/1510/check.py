"""Exact interval elimination for ground-state neutrinoless double beta decay."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def common_range(center, radius):
    return max(c - d for c, d in zip(center, radius)), min(c + d for c, d in zip(center, radius))


def solve_box(labels, ratio, width, rho):
    """Light mass x (eV), nuisance y = -m_e eta_H (eV), centers r_i x."""
    lower = [(1 - rho) * r for r in ratio]
    upper = [(1 + rho) * r for r in ratio]
    common_low, common_high = common_range(ratio, [rho * r for r in ratio])
    if common_low <= common_high:
        return {"bound_eV": None, "common_ratio_witness": (common_low + common_high) / 2}
    candidates = []
    for i, j in itertools.permutations(range(len(labels)), 2):
        if lower[i] > upper[j]:
            bound = (width[i] + width[j]) / (lower[i] - upper[j])
            candidates.append((bound, i, j))
    x, i, j = min(candidates)
    midpoint = [x * r for r in ratio]
    radii = [d + rho * mid for d, mid in zip(width, midpoint)]
    lo, hi = common_range(midpoint, radii)
    y = (lo + hi) / 2
    witness = [max(l, min(y / x, u)) for l, u in zip(lower, upper)]
    ratios = [((r * x - y) / d) ** 2 for r, d in zip(witness, width)]
    assert math.isclose(max(ratios), 1, rel_tol=2e-11)
    assert math.isclose(ratios[i], 1, rel_tol=2e-11) and math.isclose(ratios[j], 1, rel_tol=2e-11)
    # A direct feasibility check at either side, including all remaining isotopes.
    for factor, expected in [(0.999999, True), (1.000001, False)]:
        centers = [factor * c for c in midpoint]
        lo, hi = common_range(centers, [d + rho * c for c, d in zip(centers, width)])
        assert (lo <= hi) == expected
    return {"bound_eV": x, "binding_pair": [labels[i], labels[j]], "minus_m_e_etaH_eV": y,
            "witness_L_over_H": dict(zip(labels, witness)), "rate_over_limit": dict(zip(labels, ratios))}


def evaluate(data, model, labels, target_io):
    ratio = [model["light"][i] / model["heavy"][i] for i in labels]
    width = []
    for iso in labels:
        row = data["isotopes"][iso]
        amplitude = 1 / (model["gA"]**2 * math.sqrt(row["G_per_1e15_year"] * 1e-15 * row["T_year"]))
        width.append(data["m_e_eV"] * amplitude / model["heavy"][iso])
    limits = [(i, d / r) for i, d, r in zip(labels, width, ratio)]
    fixed = solve_box(labels, ratio, width, 0.0)
    x = fixed["bound_eV"]
    crit = (max(ratio) - min(ratio)) / (max(ratio) + min(ratio))
    blind = (1 + crit) * min(ratio)
    assert all(abs(blind / r - 1) <= crit + 1e-14 for r in ratio)
    return {"single_mechanism_eV": dict(limits), "fixed": fixed, "rho_blind": crit,
            "blind_ratio_at_threshold": blind, "worst_case_over_ratio_box": {
                str(rho): solve_box(labels, ratio, width, rho) for rho in data["differential_ratio_half_widths"]},
            "common_half_life_factor_for_IO_floor": (x / target_io)**2,
            "at_katrin_minimax_max_rate_over_limit": (data["katrin_upper_eV"] / x)**2}


def isolated_improvement(data, model, mass_targets):
    labels = ["Ge76", "Te130", "Xe136", "Mo100"]
    widths, ratios = {}, {}
    for i in labels:
        row = data["isotopes"][i]
        scale = math.sqrt(row["T_year"] * row["G_per_1e15_year"] * 1e-15)
        widths[i] = data["m_e_eV"] / (model["gA"]**2 * model["heavy"][i] * scale)
        ratios[i] = model["light"][i] / model["heavy"][i]
    base = labels[:-1]
    floor = min(widths[i] / abs(ratios[i] - ratios["Mo100"]) for i in base)
    answer = {}
    for name, x in mass_targets.items():
        lo, hi = common_range([x * ratios[i] for i in base], [widths[i] for i in base])
        assert lo <= hi
        center = x * ratios["Mo100"]
        gap = max(lo - center, center - hi, 0.0)
        answer[name] = data["isotopes"]["Mo100"]["T_year"] * (widths["Mo100"] / gap)**2 if gap else None
    return {"bound_eV_at_perfect_Mo": floor, "Mo_half_life_needed_year": answer}


def calculate():
    data = read(HERE / "inputs.json")
    target_io = read(HERE / data["rank2_reference"])["IO"]["m_bb_rank2_window"][0]
    labels = ["Ge76", "Te130", "Xe136"]
    result = {"io_rank2_central_floor_eV": target_io}
    for name, model in data["matrix_models"].items():
        result[name] = evaluate(data, model, labels, target_io)
    model = data["matrix_models"]["IBM2"]
    result["IBM2_plus_Mo"] = evaluate(data, model, labels + ["Mo100"], target_io)
    # Rectangular illustration only: neither independence nor a confidence distribution is known.
    lo, hi = [], []
    for i in labels + ["Mo100"]:
        l, h = model["light"][i], model["heavy"][i]
        dl, dh = model["light_error_TableVI"][i], model["heavy_error_TableVI"][i]
        lo.append((l - dl) / (h + dh))
        hi.append((l + dl) / (h - dh))
    assert max(lo) < min(hi)
    result["independent_matrix_error_box_common_ratio"] = [max(lo), min(hi)]
    result["only_Mo_improved"] = isolated_improvement(data, model, {"KATRIN": data["katrin_upper_eV"], "rank2_IO_floor": target_io})
    return result


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-16), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = calculate()
    path = HERE / "results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as handle:
            json.dump(out, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        compare(out, read(path))
        print("1510 four isotope and nuclear-input checks passed.")
