"""Independent heavy neutral pole and neutral-current orientation tests."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def projection_checks(z):
    # Independent saturation via the two orthogonal rows of the top current.
    tilts = {}
    for top_cap, ratio_cap, label in [
        (1.03, 1.01 / .94, "nominal"),
        (1.05, 1.01 / .94 + .02, "registered")
    ]:
        slope = z * (.98 if label == "registered" else 1)
        length = math.hypot(1, slope)
        atlas = math.sqrt(length ** 2 - top_cap ** 2)
        # For fixed |t/a|, positive t/a uses more of the light top projection.
        light_row = np.array([slope, ratio_cap - 1])
        light_row /= np.linalg.norm(light_row)
        orthogonal_row = np.array([-light_row[1], light_row[0]])
        cms = float(orthogonal_row @ [1, slope])
        angle = math.atan2(slope, 1) - math.acos(top_cap / length)
        assert abs(math.cos(angle) + slope*math.sin(angle) - top_cap) < 1e-12
        assert math.cos(angle) >= (.88 if label == "registered" else .90)
        tilts[label] = {"atlas_even_top_min": atlas, "cms_even_top_min": cms}
    return tilts


def heavy_rates(z, max_mt, min_m, low_y, transfer):
    beta = math.sqrt(1 - 4 * (max_mt/min_m)**2)
    scalar = np.array([1., z, 0.])
    pseudo = np.array([0., 0., z])
    width_matrix = beta**3*np.outer(scalar, scalar) + beta*np.outer(pseudo, pseudo)
    eig, directions = np.linalg.eigh(width_matrix)
    residual = float(eig[0] + eig[1])
    assert np.allclose(residual, min(beta**3*(1+z*z), beta*z*z))
    # Put the most strongly coupled direction into the unregulated light slot.
    heavy_rows = directions[:, :2].T
    assert abs(np.trace(heavy_rows @ width_matrix @ heavy_rows.T)-residual) < 1e-12
    # Exact alignment, with an arbitrary physical heavy-plane rotation.
    rotations = np.array([[0., 3/5, 4/5], [0., -4/5, 3/5]])
    aligned = float(np.trace(rotations @ width_matrix @ rotations.T))
    assert abs(aligned - z*z*(beta+beta**3)) < 1e-12
    coeff = transfer*3*low_y**2/(16*math.pi)
    return {"beta_min": beta, "summed_partial_width_over_mass_min": coeff*residual,
            "aligned_sum_min": coeff*aligned,
            "cp_odd_pole_min_in_conserving_model": coeff*beta*z*z}


def compute():
    previous = read(HERE.parent / "1586" / "results.json")
    p = previous["second_doublet"]["vev_share_upper"]
    z = math.sqrt(1/(p*p)-1)
    mt_cap = max(read(HERE.parent / "1535" / "inputs.json")["kinetic_top_mass_GeV"])
    # Antusch et al. MZ inputs, PDG2024 row.
    y_mid, y_err = read(HERE.parent / "1386" / "inputs.json")["yukawa"]["up"][-1]
    y_floor = y_mid - 2*y_err
    scan_low, atlas_high = .4, 3.49
    predictions = {}
    for label, slope in [("same_scale", z), ("two_percent_slope_floor", .98*z)]:
        tanmax = 1/slope
        assert scan_low < tanmax < atlas_high
        predictions[label] = {
            "z_min": slope, "tan_beta_max": tanmax,
            "covered_excluded_intersection": [scan_low, tanmax],
            "ratio_limit_to_projected_tan_beta": atlas_high / tanmax,
            "tt_in_400_to_1000_GeV_window": heavy_rates(slope, mt_cap, 400, y_floor, .8)
        }
    return {"source_top_doublet_share_max": p,
            "actual_heavy_search_tan_beta_lower_limit_at_400_GeV": atlas_high,
            "corresponding_minimum_p_within_scanned_benchmark":
            atlas_high/math.hypot(1,atlas_high),
            "test_two_width_marginal_witnesses": projection_checks(z),
            "conditional_spectra": predictions}


def near(left, right):
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(near(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(near(a,b) for a,b in zip(left,right))
    return math.isclose(left, right, rel_tol=2e-12, abs_tol=2e-14)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    output = HERE / "independent_results.json"
    result = compute()
    if args.write_new:
        with output.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
    else:
        assert near(result, read(output)), "independent archived result changed"
    print(json.dumps(result, ensure_ascii=True))
