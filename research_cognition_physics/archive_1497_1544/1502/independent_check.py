"""Independent FX audit. Inputs are the registered running masses and decay data."""
import argparse
import cmath
import json
from math import asin, atan, degrees, pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load_json(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8-sig"))


def endpoint(row, width, sign):
    return row[0] + sign * width * sum(row[1:])


def square_cabibbo(p, q, v):
    return ((q + p) / (1 + p) - v * v * p) / (1 + q)


def constrained_beta(masses, width):
    data = load_json("../1381/inputs.json")
    belle = load_json("../1412/inputs.json")["BelleII_Vcb_BGL"]
    q = endpoint(masses[0], width, -1) / endpoint(masses[1], width, 1)
    q_high = endpoint(masses[0], width, 1) / endpoint(masses[1], width, -1)
    k = endpoint(data["K_product"], width, 1) / endpoint(data["K_form"], width, -1)
    v = endpoint(belle, width, 1)
    d = k * k * (1 + q)
    b = 1 - v * v - d
    # Stable small root of v^2 p^2 - b p + (d-q) = 0.
    if d < q:
        raise ValueError("No low-p solution at the registered endpoints")
    p = 2 * (d - q) / (b + sqrt(b * b - 4 * v * v * (d - q)))
    assert 0 < p < 1 and abs(square_cabibbo(p, q, v) - k * k) < 1e-14
    c = sqrt(1 - v * v * (1 + p))
    t = sqrt(p / q) / c
    bounded = 2 * t / (1 + t * t)
    high_branch_bound = 2 * sqrt(q_high) / (1 + q_high)
    assert t < 1 and high_branch_bound < bounded
    measured_lower = 0.717 - width * (0.013 + 0.008)
    return {"q_min": q, "k_max": k, "v_max": v, "low_p_max": p,
            "all_p_abs_sin2beta_upper": bounded, "p_ge_1_upper": high_branch_bound,
            "LHCb_Run2_lower": measured_lower,
            "required_additive_budget": max(0., measured_lower - bounded),
            "required_phase_offset_degrees_if_only_phase": degrees(asin(measured_lower) - asin(bounded))}


def direct(q, u, v, sign):
    d = sqrt(q)
    cu = 1 / sqrt(1 + u*u)
    cd = 1 / sqrt(1 + d*d)
    sf = v / cu
    c = sqrt(1 - sf*sf)
    ru = np.array([[cu, -cu*u, 0.], [cu*u, cu, 0.], [0., 0., 1.]])
    rd = np.array([[cd, -cd*d, 0.], [cd*d, cd, 0.], [0., 0., 1.]])
    middle = np.array([[-1j*sign, 0., 0.], [0., c, sf], [0., -sf, c]])
    mat = ru.T @ middle @ rd
    def invariant(i, j, k, l):
        return mat[i,j]*mat[k,l]*mat[i,l].conjugate()*mat[k,j].conjugate()
    j = invariant(0,0,1,1).imag
    alpha = cmath.phase(-mat[2,0]*mat[2,2].conjugate()/(mat[0,0]*mat[0,2].conjugate()))
    beta = cmath.phase(-mat[1,0]*mat[1,2].conjugate()/(mat[2,0]*mat[2,2].conjugate()))
    gamma = cmath.phase(-mat[0,0]*mat[0,2].conjugate()/(mat[1,0]*mat[1,2].conjugate()))
    assert abs(alpha - sign*(pi/2-atan(u*d*c))) < 1e-14
    assert abs(beta - sign*atan(u/(d*c))) < 1e-14
    assert abs(gamma + alpha + beta - sign*pi) < 1e-14
    assert abs(abs(mat[0,1])**2 - square_cabibbo(u*u, q, v)) < 1e-14
    return {"J": j, "alpha_mod_180_deg": degrees(alpha) % 180,
            "beta_signed_deg": degrees(beta), "gamma_mod_180_deg": degrees(gamma) % 180,
            "cabibbo": abs(mat[0,1])}


def compute():
    masses = load_json("../1386/inputs.json")["yukawa"]["down"]
    high = load_json("../1397/inputs.json")["high_scale"]["down"]
    u = load_json("../1381/inputs.json")["ub_over_cb"][0]
    v = load_json("../1412/inputs.json")["BelleII_Vcb_BGL"][0]
    return {"central_MZ": [direct(masses[0][0]/masses[1][0], u, v, sign) for sign in (1, -1)],
            "MZ": {str(n): constrained_beta(masses, n) for n in (1,2)},
            "high_conditional_transfer": {str(n): constrained_beta(high, n) for n in (1,2)},
            "LHCb_S_Run2": [0.717, 0.013, 0.008],
            "BelleII_alpha_published_branch_deg": [92.6, 4.7, 4.5],
            "mixed_Z3_color_anomaly_twoT_mod3_minimal_quarks": 2}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    results = compute()
    output = HERE / "independent_results.json"
    if args.check:
        assert json.loads(output.read_text(encoding="utf-8-sig")) == results
        print("1502: independent CKM invariants and joint K / B-phase bounds agree.")
    else:
        output.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
