"""Round 1368. Finite two-family texture with a calibrated weak left port."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


DIR = Path(__file__).resolve().parent


def p_of(r):
    return r + 1 / r


def hierarchy(p):
    return (p + math.sqrt(p * p - 4)) / 2


def t_of(r, theta):
    a_over_b = math.cos(theta) ** 2 / r + r * math.sin(theta) ** 2
    return p_of(r) - p_of(a_over_b)


def weight(r, lepton, theta):
    return math.sqrt(1 + (p_of(lepton) - p_of(r)) / t_of(r, theta))


def lepton_for_four(r, theta):
    return hierarchy(p_of(r) + 15 * t_of(r, theta))


def angles_at_four(r, lepton):
    s = p_of(r) - (p_of(lepton) - p_of(r)) / 15
    upper = hierarchy(s)
    return [math.asin(math.sqrt((a - 1 / r) / (r - 1 / r)))
            for a in [1 / upper, upper]]


def extremum_angles(r, angles):
    theta_opt = math.atan(1 / math.sqrt(r))
    return [*angles, float(np.clip(theta_opt, *angles))]


def box_response(r, lepton, angles):
    # Inverse weight decreases with r; predicted charged-lepton hierarchy increases with r.
    min_w = min(weight(r[1], lepton[0], a) for a in extremum_angles(r[1], angles))
    max_w = max(weight(r[0], lepton[1], a) for a in angles)
    min_lepton = min(lepton_for_four(r[0], a) for a in angles)
    max_lepton = max(lepton_for_four(r[1], a) for a in extremum_angles(r[1], angles))
    return dict(weight=[min_w, max_w], lepton_for_four=[min_lepton, max_lepton])


def conservative_four_angles(r, lepton):
    # Monotone enclosures of the two angle branches before excluding a weak port.
    s_lo = (16 * p_of(r[0]) - p_of(lepton[1])) / 15
    s_hi = (16 * p_of(r[1]) - p_of(lepton[0])) / 15
    z_lo, z_hi = hierarchy(s_lo), hierarchy(s_hi)
    out = []
    for z0, z1 in [(1 / z_hi, 1 / z_lo), (z_lo, z_hi)]:
        sin_sq_lo = (z0 - 1 / r[0]) / (r[1] - 1 / r[1])
        sin_sq_hi = (z1 - 1 / r[1]) / (r[0] - 1 / r[0])
        out.append([math.asin(math.sqrt(max(0, sin_sq_lo))),
                    math.asin(math.sqrt(min(1, sin_sq_hi)))])
    return out


def ratio_with_box(central, error, numerator, denominator):
    return (central[numerator] / central[denominator],
            [(central[numerator] - error[numerator]) / (central[denominator] + error[denominator]),
             (central[numerator] + error[numerator]) / (central[denominator] - error[denominator])])


def angle_distance(angles, cb):
    return max(0.0, cb[0] - angles[1], angles[0] - cb[1])


def run():
    params = json.loads((DIR / "inputs.json").read_text(encoding="utf-8-sig"))
    inherited = json.loads((DIR / params["mass_data"]).read_text(encoding="utf-8-sig"))
    results = {}
    for name, data in inherited["points"].items():
        central = np.array(data["central"]) * inherited["units"]
        error = np.array(data["half_width"]) * inherited["units"]
        r, rb = ratio_with_box(central, error, 1, 0)
        lepton, lb = ratio_with_box(central, error, 3, 2)
        item = params[name]
        u, du = item["up"]
        c, dc = item["charm"]
        cab, delta_cab = item["cabibbo"]
        transform = math.asin if item["angle_type"] == "sine" else float
        cb = [transform(cab - delta_cab), transform(cab + delta_cab)]
        cab = transform(cab)
        up_angle = math.atan(math.sqrt(u / c))
        up_box = [math.atan(math.sqrt((u - du) / (c + dc))),
                  math.atan(math.sqrt((u + du) / (c - dc)))]
        allowed_angles = [cb[0] - up_box[1], cb[1] + up_box[1]]
        assert allowed_angles[0] > 0 and allowed_angles[1] < math.pi / 2
        theta_four = angles_at_four(r, lepton)
        angle_boxes = conservative_four_angles(rb, lb)
        min_distance = min(angle_distance(pair, cb) for pair in angle_boxes)
        assert min_distance > up_box[1]
        # Recover both weak angles from left eigenvectors of actual lower-mass states.
        for theta in theta_four:
            ratio_a_b = math.cos(theta) ** 2 / r + r * math.sin(theta) ** 2
            a = math.sqrt(ratio_a_b)
            mat = np.array([[0, a], [1 / a, math.sqrt(t_of(r, theta))]])
            eig, vec = np.linalg.eigh(mat @ mat.T)
            assert np.isclose(eig[1] / eig[0], r ** 2)
            assert np.isclose(math.acos(abs(vec[0, 0])), theta)
            assert np.isclose(lepton_for_four(r, theta), lepton)
        result = dict(
            cabibbo_radians=cab, symmetric_up_angle_radians=up_angle, up_angle_box_radians=up_box,
            input_down_ratio_box=rb, input_lepton_ratio_box=lb,
            aligned=dict(inferred_w=weight(r, lepton, cab),
                         lepton_predicted_w4=lepton_for_four(r, cab), box=box_response(rb, lb, cb)),
            symmetric_up=dict(permitted_angle_box_radians=allowed_angles,
                              box=box_response(rb, lb, allowed_angles)),
            w4_angles_deg=[math.degrees(a) for a in theta_four],
            w4_min_delta_deg_central=min(math.degrees(abs(a - cab)) for a in theta_four),
            w4_angle_enclosures_deg=[[math.degrees(x) for x in b] for b in angle_boxes],
            w4_min_delta_deg_necessary_over_box=math.degrees(min_distance),
        )
        assert result["symmetric_up"]["box"]["weight"][1] < 4
        results[name] = result
    return results


def same(actual, expected):
    if isinstance(actual, dict):
        return actual.keys() == expected.keys() and all(same(v, expected[k]) for k, v in actual.items())
    if isinstance(actual, list):
        return len(actual) == len(expected) and all(same(a, b) for a, b in zip(actual, expected))
    return bool(np.isclose(actual, expected, rtol=1e-10, atol=1e-12))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    output = run()
    path = DIR / "results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as target:
            json.dump(output, target, indent=2, allow_nan=False)
            target.write("\n")
    else:
        assert same(output, json.loads(path.read_text(encoding="utf-8")))
        print("PASS: Cabibbo projection, analytic mass-ratio bounds and both four-weight branches.")

