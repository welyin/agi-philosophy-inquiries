"""Independent two-family fit using determinant-one matrices and Hermitian spectra."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


# Order: d, s, e, mu. Central values and marginal HPD half-widths,
# in their actual Yukawa units. Source: Antusch–Hinze–Saad v2, Eq. 4 and Table 2.
DATA = {
    "M_Z": [(1.54e-5, .02e-5), (3.06e-4, .04e-4),
            (2.77713e-6, .00036e-6), (5.85042e-4, .00075e-4)],
    "1e16_GeV": [(.65e-5, .01e-5), (1.29e-4, .02e-4),
                 (2.6935e-6, .0091e-6), (5.6745e-4, .0192e-4)],
}


def shape(ratio):
    return ratio + 1 / ratio


def lower_weight(values):
    down, strange, electron, muon = values
    r, s = strange / down, muon / electron
    assert s >= r > 1
    # For r>1, sqrt(r+1/r-2)=(r-1)/sqrt(r).
    return (s - 1) / (r - 1) * math.sqrt(r / s)


def checked_ratio(matrix):
    h = matrix @ matrix.conj().T
    small, large = np.linalg.eigvalsh(h)
    return math.sqrt(large / small)


def witnesses(r, s, weight):
    # Normalize |det D|=1 and solve for the two possible first row norms.
    c2 = (shape(s) - shape(r)) / (weight**2 - 1)
    total = shape(r) - c2
    assert c2 >= 0 and total >= 2
    larger = (total + math.sqrt((total - 2) * (total + 2))) / 2
    angles = []
    for a2 in (1 / larger, larger):
        a, b = math.sqrt(a2) * np.exp(.41j), np.exp(-.22j) / math.sqrt(a2)
        c = math.sqrt(c2) * np.exp(.58j)
        down = np.array([[0, a], [b, c]])
        lepton = down.copy()
        lepton[1, 1] *= -weight
        assert abs(checked_ratio(down) / r - 1) < 1e-11
        for test in (lepton, lepton.T):
            assert abs(checked_ratio(test) / s - 1) < 1e-11
        _, left_vectors = np.linalg.eigh(down @ down.conj().T)
        angle = math.degrees(math.asin(abs(left_vectors[0, 1])))
        assert abs(abs(left_vectors[0, 1])**2 - (a2 - 1 / r) / (r - 1 / r)) < 1e-12
        # A free up-sector left basis realizes either trial CKM angle.
        for theta in (0., .225):
            v = np.array([[math.cos(theta), math.sin(theta)],
                          [-math.sin(theta), math.cos(theta)]])
            u_up = left_vectors @ v.conj().T
            assert np.allclose(u_up.conj().T @ left_vectors, v)
        angles.append(round(angle, 8))
    return angles


def get_summary():
    summary = {}
    for scale, rows in DATA.items():
        center = [x for x, _ in rows]
        corners = list(itertools.product(*[(x - e, x + e) for x, e in rows]))
        thresholds = [lower_weight(c) for c in corners]
        r, s = center[1] / center[0], center[3] / center[2]
        max_r = max(c[1] / c[0] for c in corners)
        min_s = min(c[3] / c[2] for c in corners)
        highest_shape = 2 + 9 * (shape(max_r) - 2)
        highest_ratio_three = (highest_shape + math.sqrt(highest_shape**2 - 4)) / 2
        assert highest_ratio_three < min_s and 3 < min(thresholds) <= max(thresholds) < 4
        summary[scale] = {
            "center_down_lepton_ratios": [round(r, 9), round(s, 9)],
            "minimum_weight_at_center": round(lower_weight(center), 9),
            "weight_threshold_extremes_marginal_box": [round(min(thresholds), 9), round(max(thresholds), 9)],
            "C3_max_allowed_ratio_across_box": round(highest_ratio_three, 9),
            "min_lepton_ratio_across_box": round(min_s, 9),
            "C4_down_left_angles_degrees_at_center": witnesses(r, s, 4),
        }
    return summary


def compare_saved(actual, saved):
    assert actual.keys() == saved.keys()
    for scale, findings in actual.items():
        assert findings.keys() == saved[scale].keys()
        for name, value in findings.items():
            # Witnesses are rounded to eight decimal places or finer in the file.
            observed = np.asarray(value)
            reference = np.asarray(saved[scale][name])
            assert observed.shape == reference.shape, (scale, name)
            np.testing.assert_allclose(observed, reference, rtol=0, atol=1.1e-8,
                                       err_msg=f"{scale}: {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = get_summary()
    path = Path(__file__).with_name("independent_results.json")
    if args.write:
        with path.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2) + "\n")
    else:
        compare_saved(result, json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(result))
