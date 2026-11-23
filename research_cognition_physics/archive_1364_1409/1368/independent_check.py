"""Invert the left Gram matrix and certify extrema over both running-mass boxes."""
import argparse
import json
from math import asin, atan, atan2, cos, isclose, sin, sqrt
from pathlib import Path

import numpy as np


def interval_ratio(a, b):
    return (a[0]-a[1])/(b[0]+b[1]), (a[0]+a[1])/(b[0]-b[1])


def branch_factor(r, e):
    return (e + 1/e - r - 1/r)/(r - 1/r)**2


def weight(r, e, angle):
    k = branch_factor(r, e)
    return sqrt(1+k*(r/cos(angle)**2+1/(r*sin(angle)**2)))


def roots_at_four(r, e):
    # Quadratic in tan(theta)^2; use the root product near the small solution.
    linear = 15/branch_factor(r, e) - r - 1/r
    large = (linear + sqrt(linear**2-4))/(2*r)
    return [atan(1/(r*sqrt(large))), atan(sqrt(large))]


def weight_window(rbox, ebox, angular_window):
    lower_angle = min(max(atan(1/sqrt(rbox[1])), angular_window[0]), angular_window[1])
    lower = weight(rbox[1], ebox[0], lower_angle)
    upper = max(weight(rbox[0], ebox[1], t) for t in angular_window)
    return [lower, upper]


def left_matrix(r, theta):
    # Prescribe the physical light and heavy eigenvectors of Y_d Y_d^dagger.
    left = np.array([[cos(theta), sin(theta)], [-sin(theta), cos(theta)]])
    gram = (left * [1/r, r]) @ left.T
    aa = sqrt(gram[0, 0])
    return np.array([[0., aa], [1/aa, gram[1, 0]/aa]]), left


def gauge_freedom(r, e, theta_d, cabibbo):
    mat, down = left_matrix(r, theta_d)
    wk = weight(r, e, theta_d)
    other = mat.copy()
    other[1, 1] *= -wk
    sv = np.linalg.svd(other, compute_uv=False)
    assert isclose(sv[0]/sv[1], e, rel_tol=4e-14)
    # A completely specified 3x3 CKM, including leakage, is allowed by arbitrary up basis.
    def rotation(i, j, angle, phase=0):
        u = np.eye(3, dtype=complex)
        a = np.exp(1j*phase)
        u[i,i] = u[j,j] = cos(angle)
        u[i,j], u[j,i] = sin(angle)*a, -sin(angle)*a.conjugate()
        return u
    target = (rotation(1,2,0.18) @ rotation(0,2,0.12,0.71)
              @ rotation(0,1,cabibbo))
    ud = np.eye(3, dtype=complex)
    ud[:2,:2] = down
    up = ud @ target.conj().T
    up_yukawa = up * [0.022, 0.8, 9.]
    vals, measured_up = np.linalg.eigh(up_yukawa @ up_yukawa.conj().T)
    mixed = measured_up.conj().T @ ud
    assert np.allclose(vals, np.array([0.022,0.8,9.])**2)
    assert np.allclose(abs(mixed), abs(target))
    projected = measured_up[:2,0]
    ray_distance = atan2(abs(projected[1]), abs(projected[0]))
    assert abs(ray_distance-abs(cabibbo-theta_d)) < 1e-7
    return {"weight": wk, "physical_first_row_moduli": abs(mixed[0]).tolist(),
            "projected_angle": ray_distance}


def examine(down, strange, electron, muon, up, charm, cab, is_sine):
    r, e = strange[0]/down[0], muon[0]/electron[0]
    rbox, ebox = interval_ratio(strange, down), interval_ratio(muon, electron)
    center_c = asin(cab[0]) if is_sine else cab[0]
    cbox = [cab[0]-cab[1],cab[0]+cab[1]]
    if is_sine:
        cbox = [asin(t) for t in cbox]
    delta = atan(sqrt(up[0]/charm[0]))
    delta_max = atan(sqrt(interval_ratio(up, charm)[1]))
    sym_angles = [cbox[0]-delta_max, cbox[1]+delta_max]
    four_angles = roots_at_four(r, e)
    distances = [abs(t-center_c) for t in four_angles]
    # Over the boxes, the favorable mass pair gives the widest superlevel region.
    closest_roots = roots_at_four(rbox[0], ebox[1])
    min_required = min(cbox[0]-closest_roots[0],closest_roots[1]-cbox[1])
    assert min_required > delta_max > 0
    windows = [weight_window(rbox, ebox, angles) for angles in (cbox,sym_angles)]
    assert all(3 < low < high < 4 for low, high in windows)
    for theta in four_angles + [center_c]:
        mat, left = left_matrix(r, theta)
        assert np.allclose(left.T @ mat @ mat.T @ left, np.diag([1/r,r]))
        lepton = mat.copy()
        lepton[1,1] *= weight(r,e,theta)
        s = np.linalg.svd(lepton, compute_uv=False)
        assert isclose(s[0]/s[1], e)
    return {
        "central_cabibbo_radians": center_c, "central_angle_from_up_zero_texture": delta,
        "aligned_central_weight": weight(r,e,center_c),
        "symmetric_central_weight_window": weight_window((r,r),(e,e),
                                                         [center_c-delta,center_c+delta]),
        "aligned_full_box_weight": windows[0], "symmetric_full_box_weight": windows[1],
        "weight_4_two_down_angles": four_angles, "central_required_up_misalignments": distances,
        "full_box_minimum_up_misalignment_for_weight_4": min_required,
        "full_box_symmetric_up_maximum": delta_max,
        "arbitrary_up_w4_demo": gauge_freedom(r,e,four_angles[0],center_c),
        "projected_alignment_with_nonzero_Vub_demo": gauge_freedom(r,e,center_c,center_c)
    }


def run():
    # Dimensionless masses, limited to the source's quoted precision.
    return {
        "source": "https://arxiv.org/html/2510.01312v2", "source_sections": "Eqs. (4),(5), Table 2",
        "MZ": examine((1.54e-5,.02e-5),(3.06e-4,.04e-4),(2.77713e-6,.00036e-6),
                      (5.85042e-4,.00075e-4),(7.04e-6,.15e-6),(3.56e-3,.06e-3),
                      (.2251,.0008),True),
        "10^16_GeV": examine((.65e-5,.01e-5),(1.29e-4,.02e-4),(2.6935e-6,.0091e-6),
                             (5.6745e-4,.0192e-4),(2.87e-6,.07e-6),(1.45e-3,.03e-3),
                             (.22708,.00082),False)
    }


def portable(value):
    if isinstance(value, dict):
        return {k: portable(v) for k,v in value.items()}
    if isinstance(value, list):
        return [portable(v) for v in value]
    return round(value, 8) if isinstance(value, float) else value


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true", help="verify frozen results (default)")
    options = parser.parse_args()
    path = Path(__file__).with_name("independent_results.json")
    result = portable(run())
    if options.write:
        with path.open("x", encoding="utf-8") as output:
            json.dump(result,output,ensure_ascii=False,indent=2)
            output.write("\n")
    else:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    print("1368 independent mass-box, SVD and three-generation checks passed.")
