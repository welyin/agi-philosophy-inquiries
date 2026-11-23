"""Independent two-source bound in the up-quark basis, from nonnegative Gram terms."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE.parent / relative).read_text(encoding="utf-8-sig"))


def box(row, widths):
    return row[0] - widths * sum(row[1:]), row[0] + widths * sum(row[1:])


def triangle_box(tree, belle, widths):
    k, f = (box(tree[key], widths) for key in ("K_product", "K_form"))
    us = k[0] / f[1], k[1] / f[0]
    cb = box(belle, widths)
    ratio = box(tree["ub_over_cb"], widths)
    ub = ratio[0] * cb[0], ratio[1] * cb[1]

    def cosine(s):
        return math.sqrt(1 - s[1] ** 2), math.sqrt(1 - s[0] ** 2)

    c13 = cosine(ub)
    s12 = us[0] / c13[1], us[1] / c13[0]
    s23 = cb[0] / c13[1], cb[1] / c13[0]
    c12, c23 = cosine(s12), cosine(s23)
    # Bounds on the moduli of each of the three columns, valid for EVERY Dirac phase.
    d_lower = np.array([c12[0] * c13[0],
                        s12[0] * c23[0] - c12[1] * s23[1] * ub[1],
                        s12[0] * s23[0] - c12[1] * c23[1] * ub[1]])
    s_upper = np.array([us[1], c12[1] * c23[1] + s12[1] * s23[1] * ub[1],
                        c12[1] * s23[1] + s12[1] * c23[1] * ub[1]])
    b_upper = np.array([ub[1], cb[1], c23[1] * c13[1]])
    assert np.all(d_lower > 0)
    s_over_d, b_over_d = s_upper / d_lower, b_upper / d_lower
    # In the last row factor out the common s23 before bounding the ratio.
    t = c23[1] * c13[1] * ratio[1]
    denom = s12[0] - c12[1] * t
    assert denom > 0
    s_over_d[2] = min(s_over_d[2], (c12[1] + s12[1] * t) / denom)
    b_over_d[2] = min(b_over_d[2], b_upper[2] / (s23[0] * denom))
    return np.ones(3), s_over_d, b_over_d


def certify(d_lower, s_upper, b_upper, squared_yukawas):
    d_max, s_min, b_min = squared_yukawas
    # The entries of D times |V_id|² share ONE budget, rather than three.
    def projected_cap(left, right):
        return d_max * np.max(left * right / d_lower**2)

    a = s_min - projected_cap(s_upper, s_upper)
    b = b_min - projected_cap(b_upper, b_upper)
    c = projected_cap(s_upper, b_upper)
    assert b > a and a > 0
    lower = a - 2 * c**2 / (b - a + math.hypot(b - a, 2 * c))
    return {"second_positive_source_norm_floor": float(lower),
            "fraction_of_y_s_squared_floor": float(lower / s_min),
            "s_diagonal_cap": float(s_min - a),
            "b_diagonal_cap": float(b_min - b),
            "sb_cap": float(c)}


def calculate():
    mass = read("1386/inputs.json")["yukawa"]["down"]
    tree = read("1381/inputs.json")
    belle = read("1412/inputs.json")["BelleII_Vcb_BGL"]
    answer = {}
    for widths in (0, 1, 2):
        directions = triangle_box(tree, belle, widths)
        squared = [box(row, widths)[index] ** 2 for row, index in zip(mass, (1, 0, 0))]
        answer[f"widths_{widths}_all_delta"] = certify(*directions, squared)
        ratios = [squared[0] / squared[1], 1., squared[2] / box(mass[1], widths)[1] ** 2]
        answer[f"widths_{widths}_all_delta"]["uniform_fraction_of_actual_y_s_squared"] = (
            certify(*directions, ratios)["second_positive_source_norm_floor"])

    # Reconstruct two simultaneous positive spectra even in the CP-even branch.
    us = tree["K_product"][0] / tree["K_form"][0]
    cb = belle[0]
    s13 = cb * tree["ub_over_cb"][0]
    c13 = math.sqrt(1 - s13**2)
    s12, s23 = us/c13, cb/c13
    c12, c23 = math.sqrt(1 - s12**2), math.sqrt(1 - s23**2)
    spectrum = np.array([row[0] ** 2 for row in mass])
    phase_checks = []
    for delta in (0, 1.2, 3.1):
        phase = np.exp(1j * delta)
        V = np.array([[c12*c13, s12*c13, s13/phase],
                      [-s12*c23-c12*s23*s13*phase, c12*c23-s12*s23*s13*phase, s23*c13],
                      [s12*s23-c12*c23*s13*phase, -c12*s23-s12*c23*s13*phase, c23*c13]])
        H = (V * spectrum) @ V.conj().T
        assert np.linalg.norm(V @ V.conj().T - np.eye(3)) < 1e-14
        D = spectrum[0] * np.eye(3)
        parts = [np.outer(V[:, i], V[:, i].conj()) * (spectrum[i]-spectrum[0]) for i in (1, 2)]
        assert np.linalg.norm(H-D-sum(parts), 2) < 1e-18
        fixed = certify(abs(V[:, 0]), abs(V[:, 1]), abs(V[:, 2]), spectrum)
        assert fixed["second_positive_source_norm_floor"] <= np.linalg.norm(parts[0], 2)
        phase_checks.append({"delta_radians_example": delta, **fixed})
    assert answer["widths_2_all_delta"]["fraction_of_y_s_squared_floor"] > 0.7
    answer["central_two_source_construction_smaller_norm"] = float(spectrum[1]-spectrum[0])
    answer["fixed_central_profiles"] = phase_checks
    return answer


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    dest = HERE / "independent_results.json"
    if args.check:
        saved = json.loads(dest.read_text(encoding="utf-8"))
        def close(found, expected):
            if isinstance(found, dict):
                return found.keys() == expected.keys() and all(close(v, expected[k]) for k, v in found.items())
            if isinstance(found, list):
                return len(found) == len(expected) and all(close(a, b) for a, b in zip(found, expected))
            return math.isclose(found, expected, rel_tol=1e-11, abs_tol=1e-18)
        assert close(result, saved)
    else:
        dest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(result["widths_2_all_delta"])
