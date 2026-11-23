"""Independent residual-line, spectrum-normalization and pole tests. Uses only local Python/NumPy."""
import argparse
import itertools
import json
from math import acos, cos, degrees, pi, sqrt
from pathlib import Path

import numpy as np

OUT = Path(__file__).with_name("independent_results.json")
GRID = [cos(k * pi / 14) ** 2 for k in range(7)] + [0.0]


def row(fast, slow):
    p3 = (1 - sqrt(1 - fast)) / 2
    return np.array([(1 - slow) * (1 - p3), slow * (1 - p3), p3])


def ray_angle(p):
    return degrees(acos(sqrt(p)))


def pole_response(asym):
    ratio = asym / (1 + sqrt(1 - asym**2))
    qed = 1 + 3 / (4 * pi * 128)
    return 2 / ((1 + ratio**2) * qed)


def window(solar, error):
    rows = [row(a, x) for a, x in itertools.product(
        (0.0851 - 0.0048, 0.0851 + 0.0048), (solar - 2 * error, solar + 2 * error))]
    p_low, p_high = np.min(rows, axis=0), np.max(rows, axis=0)
    result = []
    for low, high in zip(p_low, p_high):
        gap = min(abs(ray_angle(g) - ray_angle(min(max(g, low), high))) for g in GRID)
        candidates = [(1 - g / max(g, low), g) for g in GRID if 0 < g <= high]
        if candidates:
            loss, g = min(candidates)
            result.append(dict(fs_gap_deg=gap, grid=g, missing=loss,
                               active_spectral_loss=loss / (1 - g)))
        else:
            result.append(dict(fs_gap_deg=gap, grid=None, missing=None))
    best = min((item for item in result if item["missing"] is not None),
               key=lambda item: item["missing"])
    ell = best["missing"]
    strongest = min(item["active_spectral_loss"] for item in result
                    if item["missing"] is not None)
    upper = 2 + (1 - ell)**2
    physical_upper = 2 + (1 - strongest)**2
    u_low = 2.9963 - 2 * 0.0074
    r_low = 5.943 - 0.032
    response = pole_response(0.1513 - 0.0042)
    return {
        "columns": result,
        "normalized_p_at_n1_extremum": min(rows, key=lambda p: p[0]).tolist(),
        "min_missing_e": ell, "U_upper": upper, "U_upper_active_line": physical_upper,
        "species_gap": u_low - upper,
        "pole_remainder_required_no_extra_invisible": r_low / (response * upper),
        "additional_width_over_charged_lepton_at_1p02": r_low - 1.02 * response * upper,
        "Rinv": r_low, "r1_max": response,
    }


def common_representation():
    p_low = np.min([row(fast, solar) for fast, solar in itertools.product(
        (0.0803, 0.0899), (0.2908, 0.3164))], axis=0)
    p_high = np.max([row(fast, solar) for fast, solar in itertools.product(
        (0.0803, 0.0899), (0.2908, 0.3164))], axis=0)
    upper = max(max(p_high), 1 - p_low[0], 1 - p_low[1],
                (1 - p_low[2]) * max(0.624, 1 - 0.402))
    # Opposite representations with identical group generators and an isolated muon.
    s12sq, s23sq = 0.3036, 0.514
    s13sq = row(0.0851, s12sq)[2]
    incoherent = s12sq * (1 - s23sq) + (1 - s12sq) * s23sq * s13sq
    interference = 2 * sqrt(s12sq * (1 - s12sq) * s23sq * (1 - s23sq) * s13sq)
    assert abs(GRID[5] - incoherent) < interference
    return {"all_nine_upper": upper, "identical_actions_min_rotation_deg":
            ray_angle(upper) - ray_angle(GRID[1]),
            "other_irrep_mu1": GRID[5], "other_irrep_abs_delta_deg":
            degrees(acos((GRID[5] - incoherent) / interference))}


def witness(stats):
    p = np.array(stats["normalized_p_at_n1_extremum"])
    g = stats["columns"][0]["grid"]
    target = stats["columns"][0]["active_spectral_loss"]
    e, other, third, sterile = np.eye(4)
    n = sqrt(g) * e + sqrt(1 - g) * other
    t = (e - sqrt(g) * n) / sqrt(1 - g)
    heavy = sqrt(target) * t + sqrt(1 - target) * sterile
    light = sqrt(1 - target) * t - sqrt(target) * sterile
    c, s = np.sqrt(p[1:] / sum(p[1:]))
    mixing = np.column_stack((n, c * light + s * third, -s * light + c * third, heavy))
    assert np.max(abs(mixing.T @ mixing - np.eye(4))) < 1e-12
    N = mixing[:3, :3]
    gram = N @ N.T
    assert np.max(abs(N[0]**2 / gram[0, 0] - p)) < 1e-12
    U = np.trace(gram @ gram)
    assert abs(stats["U_upper_active_line"] - U) < 1e-12
    charged = 2 * np.outer(e, e) - np.eye(4)
    neutral = 2 * np.outer(n, n) - np.eye(4)
    generator = neutral @ charged
    assert np.linalg.norm(np.linalg.matrix_power(generator, 14) - np.eye(4)) < 1e-12
    for power in (1, 2, 7):
        assert np.linalg.norm(np.linalg.matrix_power(generator, power) - np.eye(4)) > 0.1
    assert np.max(abs(charged @ generator @ charged - generator.T)) < 1e-12
    return dict(mixing=mixing.tolist(), gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
                U=U, missing_e=1 - gram[0, 0],
                gmu_over_g0=sqrt(gram[0, 0] * gram[1, 1]))


def formatted(data):
    if isinstance(data, dict):
        return {k: formatted(v) for k, v in data.items()}
    if isinstance(data, list):
        return [formatted(v) for v in data]
    if isinstance(data, (float, np.floating)):
        return round(float(data), 10)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_mutually_exclusive_group()
    command.add_argument("--write", action="store_true")
    command.add_argument("--check", action="store_true")
    options = parser.parse_args()
    no = window(0.3036, 0.0064)
    io = window(0.307, 0.013)
    items = dict(normal=no, inverted_historical_solar=io,
                 faithful_comparison=common_representation(),
                 conditioned_electron_loss_upper_LEP=1 - sqrt(2.9963 - 0.0148 - 2),
                 reactor_and_active_line_witness=witness(no))
    assert 0.0847 < no["min_missing_e"] < 0.085
    assert io["pole_remainder_required_no_extra_invisible"] > 1.02
    result = formatted(items)
    if options.write:
        with OUT.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    else:
        if result != json.loads(OUT.read_text(encoding="utf-8")):
            raise SystemExit("Independent archive mismatch")
        print("1616 independent reproduction passed")


if __name__ == "__main__":
    main()
