"""Two-source plane geometry and retained reactor-data checks for round 1371."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def small_reactor_root(amplitude):
    return amplitude / (2 * (1 + math.sqrt(1 - amplitude)))


def solar(x):
    return 1 - 2 / (3 * (1 - x))


def cp_coefficient(x):
    return (1 - 5 * x) / (2 * math.sqrt(2 * x * (1 - 3 * x)))


def orthogonal_angle(electron_probability):
    return math.acos(math.sqrt(electron_probability))


def test_sources():
    null = np.array([2.0, -1.0, 1.0]) / math.sqrt(6)
    test_kernel = np.array([[1.4 + 0.3j, -0.9 + 0.2j], [-0.9 + 0.2j, 0.1 - 1.3j]])
    assert abs(np.linalg.det(test_kernel)) > 0.1
    rephase = np.diag(np.exp(1j * np.array([0.32, -0.17, 0.51])))
    change = np.array([[1 + 1j, 0.8], [0.2, 1.4 - 0.2j]])
    probes = []
    for n in [0, 1, 3, 2 + 0.7j]:
        source = np.array([[0, 1], [1, n], [1, n - 2]], dtype=complex)
        mass = source @ test_kernel @ source.T
        singular = np.linalg.svd(mass, compute_uv=False)
        assert singular[1] > 1e-2 and np.linalg.norm(mass @ null) < 1e-12
        moved = source @ change
        kernel_moved = np.linalg.solve(change, test_kernel) @ np.linalg.inv(change).T
        assert np.allclose(mass, moved @ kernel_moved @ moved.T)
        phased = rephase @ mass @ rephase.T
        assert np.linalg.norm(phased @ (rephase.conj() @ null)) < 1e-12
        probes.append(float(singular[1]))
    # Rotate the full source relative to the fixed, physically named charged ports.
    angle = 0.3
    rotate = np.array([[math.cos(angle), -math.sin(angle), 0],
                       [math.sin(angle), math.cos(angle), 0], [0, 0, 1]])
    rotated = rotate @ mass @ rotate.T
    assert np.allclose(np.linalg.svd(rotated, compute_uv=False), singular)
    other_e_fraction = float((rotate @ null)[0] ** 2)
    return dict(second_singular_value_checks=probes, rotated_zero_e_fraction=other_e_fraction)


def run():
    data = json.loads((HERE / "inputs.json").read_text(encoding="utf-8-sig"))
    a = data["daya_bay"]["sin_sq_2theta13"]
    sa = data["daya_bay"]["one_sigma"]
    j = data["juno_primary"]["sin_sq_theta12"]
    sj = data["juno_primary"]["one_sigma"]
    x = small_reactor_root(a)
    xlo, xhi = small_reactor_root(a - sa), small_reactor_root(a + sa)
    y = solar(x)
    ylo, yhi = solar(xhi), solar(xlo)
    closest = max(0, ylo - (j + sj), (j - sj) - yhi)
    sigma_y = sa / (6 * math.sqrt(1 - a) * (1 - x) ** 2)
    # Second moment only: covariance of estimators need not vanish.
    error_lo, error_hi = abs(sj - sigma_y), sj + sigma_y
    coeff = cp_coefficient(x)
    lower23 = 0.5 * math.atan(coeff)
    min_tilt = orthogonal_angle(xhi) - orthogonal_angle(2 / 3)
    assert closest > 0 and min_tilt > 0
    assert math.isclose(math.sin(2 * math.asin(math.sqrt(2 / 3))) ** 2, 8 / 9)
    # Only the normal-order mass splittings are inherited; no fitted angles enter.
    prior = json.loads((HERE.parent / "1370" / "inputs.json").read_text(encoding="utf-8-sig"))
    delta, atm = prior["NO"]["solar"][0], prior["NO"]["atm_abs"][0]
    p2 = 1 / 3 - x
    beta = math.sqrt(p2 * delta + x * atm)
    term2, term3 = p2 * math.sqrt(delta), x * math.sqrt(atm)
    return dict(
        beta_common_normal_eV=beta, majorana_ee_normal_eV=[abs(term2-term3), term2+term3],
        reactor_s13_squared=x, reactor_one_sigma=[xlo, xhi],
        no_solar_central=y, no_solar_from_reactor_one_sigma=[ylo, yhi],
        juno_one_sigma=[j - sj, j + sj], distance_of_one_sigma_intervals=closest,
        primary_center_difference=y - j, primary_difference_in_juno_quoted_widths=(y - j) / sj,
        solar_linearized_reactor_error=sigma_y,
        unknown_correlation_linearized_rms_bounds=[error_lo, error_hi],
        center_distance_over_possible_linearized_rms=[(y - j) / error_hi, (y - j) / error_lo],
        io_predicted_reactor_sin_sq_2theta13=8 / 9,
        min_io_tilt_degrees_small_branch_over_reactor_box=math.degrees(min_tilt),
        no_cp_cos_delta_over_minus_cot_2theta23=coeff,
        no_allowed_theta23_degrees=[math.degrees(lower23), 90 - math.degrees(lower23)],
        source_kernel_checks=test_sources()
    )


def same(data, stored):
    if isinstance(data, dict):
        return data.keys() == stored.keys() and all(same(v, stored[k]) for k, v in data.items())
    if isinstance(data, list):
        return len(data) == len(stored) and all(same(a, b) for a, b in zip(data, stored))
    return math.isclose(data, stored, rel_tol=1e-10, abs_tol=1e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    values = run()
    path = HERE / "results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as out:
            json.dump(values, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        assert same(values, json.loads(path.read_text(encoding="utf-8")))
        print("PASS: common plane, ordering, reactor sum rule, tilted-plane and CP checks.")

