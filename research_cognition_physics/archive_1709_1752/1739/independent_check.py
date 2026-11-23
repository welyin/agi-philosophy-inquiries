"""Reuse the five registered color poles. Run with --check to verify the saved result."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def previous(number):
    return json.loads((HERE.parent / str(number) / "results.json").read_text(encoding="utf-8"))


def test_rotation(index, theta, phi, scalar, full_y, poles):
    """The external Higgs source stays fixed while the heavy mass rows rotate."""
    c, s = math.cos(theta), math.sin(theta)
    z = np.exp(1j * phi)
    rot = np.eye(3, dtype=complex)
    rot[np.ix_([0, index], [0, index])] = [[c, -s * z.conjugate()], [s * z, c]]
    hv = 246.22 / math.sqrt(2)
    old_heavy = np.diag([0, *poles])[1:]
    new_heavy = old_heavy @ rot.conj().T
    old = np.vstack([full_y * hv, old_heavy])
    new = np.vstack([full_y * hv, new_heavy])
    ratio = np.linalg.det(new) / np.linalg.det(old)
    ratio_from_light = (full_y @ rot[:, 0]) / full_y[0]
    assert abs(ratio - ratio_from_light) < 2e-12
    assert np.allclose(np.linalg.svd(new_heavy, compute_uv=False), poles[::-1])
    assert abs(abs(ratio) - scalar) < 2e-12
    phase = float(np.angle(ratio))
    repaired = new.copy()
    repaired[index, :] *= np.exp(-1j * phase)
    repaired_ratio = np.linalg.det(repaired) / np.linalg.det(old)
    assert abs(repaired_ratio - scalar) < 2e-12
    transferred = abs((full_y @ rot)[index]) ** 2 - full_y[index] ** 2
    budget = (1 - scalar**2) * full_y[0] ** 2
    assert math.isclose(transferred, budget, abs_tol=3e-20)
    kappa = math.sqrt((poles[1] / poles[0]) ** 2 - 1)
    # Ordinary tagged d^c overlaps the old (light, low, high) lines as (1,0,-k)/R.
    ordinary_rate_ratio = abs(rot[0, 0] - kappa * rot[2, 0]) ** 2
    m = poles[index - 1]
    q0 = poles[0] * math.exp(-0.75)
    delta_higgs = -3 * m**2 / (8 * math.pi**2) * (math.log((m/q0)**2) - 1) * budget
    return dict(angle_rad=theta, interference_phase_rad=phi, old_new_light_overlap=c,
                new_mass_over_old=float(abs(ratio)), delta_arg_full_color_det_rad=phase,
                heavy_row_repair_phase_rad=-phase,
                repaired_full_det_over_old_real=float(repaired_ratio.real),
                norm_of_new_row_to_old_light_GeV=m*abs(s),
                higgs_weight_gained_by_changed_pole=budget,
                delta_pure_fermion_higgs_c_GeV2=delta_higgs,
                tagged_ordinary_d_to_lepton_rate_ratio=float(ordinary_rate_ratio))


def calculate():
    old = previous(1720)["central_conditional_combined_point"]
    tag = previous(1712)["cases"]["interior_all_one_width"]
    yd, ys, yb = tag["upper_yukawas"]
    f = math.sqrt(tag["upper_probabilities"][0])
    s = math.sqrt(1 - f*f)
    r = old["R"]
    kappa = math.sqrt(r*r - 1)
    g = yb * r / kappa
    assert math.isclose(g, old["g"], rel_tol=5e-11)
    poles = old["five_color_poles_GeV"][::4]
    assert poles[1] == r*poles[0]
    low_yd = g * math.sqrt(1 - (yd/yb)**2)
    low_ys = g * math.sqrt(1 - (ys/yb)**2)
    higgs = np.array([yd, low_yd, yd/kappa])
    min_changed_col = s*low_ys - f*yd
    assert ys < min_changed_col < yb
    allowed_y = 2*f*yd/s
    ratio_low = low_yd/yd
    small_angle = 2*math.atan((f - 1)/(ratio_low + math.sqrt(ratio_low**2 + s*s)))
    low = test_rotation(1, small_angle, 0, f, higgs, poles)
    high = test_rotation(2, math.acos(f), math.acos(-s/(2*kappa*f)), f, higgs, poles)
    edm_limit = previous(1389)["conditional_edm"]["theta_limit_rad"]
    assert math.isclose(abs(high["delta_arg_full_color_det_rad"]),
                        2*math.asin(s/(2*kappa*f)), abs_tol=2e-14)
    assert math.isclose(high["tagged_ordinary_d_to_lepton_rate_ratio"],
                        1+kappa*kappa*s*s, rel_tol=1e-12)
    high["conditional_strong_angle_over_old_EDM_budget"] = (
        abs(high["delta_arg_full_color_det_rad"])/edm_limit)
    high["tagged_ordinary_d_probability"] = high["tagged_ordinary_d_to_lepton_rate_ratio"]/r**2
    return dict(reference="1712 interior_all_one_width before 1721 threshold retune; 1720 R=800",
                y_light_dsb=[yd,ys,yb], f= f, R=r, pole_GeV=poles,
                old_strong_angle_budget_rad=edm_limit,
                any_low_pole_at_required_geometric_overlap=dict(
                    min_old_low_Y=low_ys, required_old_low_Y_upper=allowed_y,
                    lower_bound_second_light_yukawa=min_changed_col,
                    bound_over_old_strange=min_changed_col/ys,
                    necessary_upper_mixed_probability_if_strange_fixed=((ys+yd)/low_ys)**2,
                    required_mixed_probability=s*s),
                low_real_interference=low, high_same_overlap_and_modulus=high)


def rounded(obj):
    if isinstance(obj, dict):
        return {key: rounded(value) for key, value in obj.items()}
    if isinstance(obj, list):
        return [rounded(value) for value in obj]
    if isinstance(obj, float):
        return float(f"{obj:.10g}")
    return obj


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    out = rounded(calculate())
    saved = HERE / "independent_results.json"
    if args.check:
        assert out == json.loads(saved.read_text(encoding="utf-8"))
        print("1739 independent check OK: two old low poles, high-pole phase, current and curvature.")
    else:
        saved.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

