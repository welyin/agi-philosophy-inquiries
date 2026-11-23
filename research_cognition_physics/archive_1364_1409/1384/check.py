"""Permutation-invariant CKM/PMNS compatibility with an isolated reactor frequency."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def bounds(entry, multiplier):
    delta = multiplier * sum(entry[1:])
    return entry[0] - delta, entry[0] + delta


def amplitude(p):
    return 4 * p * (1 - p)


def small_root(f):
    return f / (2 * (1 + math.sqrt(1-f)))


def folded_angle(f):
    return 0.5 * math.asin(math.sqrt(f))


def rotation_witness(a, b):
    """SO(3) with diagonal (a,b,a+b-1), saturating the two-tag bound."""
    c = a+b-1
    x, y = math.sqrt((1-b)/(1-c)), math.sqrt((1-a)/(1-c))
    axis = np.array([x, y, 0.0])
    cross = np.array([[0, 0, y], [0, 0, -x], [-y, x, 0]])
    return c*np.eye(3) + (1-c)*np.outer(axis, axis) + math.sqrt(1-c*c)*cross


def mixing(s12, s13, s23, phase):
    c12, c13, c23 = [math.sqrt(1-s*s) for s in (s12, s13, s23)]
    eta = np.exp(1j*phase)
    return np.array([[c12*c13, s12*c13, s13/eta],
                     [-s12*c23-c12*s23*s13*eta, c12*c23-s12*s23*s13*eta, s23*c13],
                     [s12*s23-c12*c23*s13*eta, -c12*s23-s12*c23*s13*eta, c23*c13]])


def frame_witness(qdata, db, juno):
    a = qdata["ud_beta"][0]
    c = qdata["F_times_cb"][0]/qdata["F_example"]
    b = c*qdata["ub_over_cb"][0]
    s12q = math.sqrt((1-a*a-b*b)/(1-b*b))
    v = mixing(s12q, b, c/math.sqrt(1-b*b), 1.14)
    p = small_root(db["sin_sq_2theta13"])
    u = mixing(math.sqrt(juno["sin_sq_theta12"]), math.sqrt(p), math.sqrt(0.5), 0.5)
    eye = np.eye(3)
    common_weak_family_map = eye
    for (x, y, expected) in ((eye, v, v), (eye, u, u)):
        assert np.linalg.norm(x.conj().T @ common_weak_family_map @ y - expected) < 1e-13
        assert np.linalg.norm(expected.conj().T @ expected - eye) < 1e-13
    # Any nondegenerate positive squared masses fix the corresponding mass eigenbases.
    return {"weak_family_map": "I_3", "moduli_quarks": abs(v).tolist(),
            "moduli_leptons": abs(u).tolist(), "reactor_amplitude": amplitude(abs(u[0,2])**2)}


def evaluate():
    cfg = read(ROOT / "inputs.json")
    q = read(ROOT / cfg["quark_inputs"])
    l = read(ROOT / cfg["reactor_inputs"])
    db, ju = l["daya_bay"], l["juno_primary"]
    out = {}
    for k in cfg["error_multipliers"]:
        amin, amax = bounds(q["ud_beta"], k)
        cmax = bounds(q["F_times_cb"], k)[1]/q["F_lower"]
        bmax = cmax*bounds(q["ub_over_cb"], k)[1]
        t = cmax**2+bmax**2
        d = 1-amax**2-t
        assert 0 < t < d < .5
        assert (1-amin**2+cmax**2) < .5
        small_ceiling, other_floor = amplitude(t), amplitude(d)
        flo, fhi = bounds([db["sin_sq_2theta13"], db["one_sigma"]], k)
        sep = min(flo-small_ceiling, other_floor-fhi)
        angle = min(folded_angle(flo)-folded_angle(small_ceiling),
                    folded_angle(other_floor)-folded_angle(fhi))
        tlo, thi = small_root(flo), small_root(fhi)
        smin, smax = bounds([ju["sin_sq_theta12"], ju["one_sigma"]], k)
        max_e = math.sqrt(max((1-tlo)*(1-smin), (1-tlo)*smax, thi))
        csmin = bounds(q["charm_s"], k)[0]
        floor_diag = amin+csmin-1
        rot = rotation_witness(amin, csmin)
        assert np.allclose(rot.T@rot, np.eye(3)) and np.allclose(np.diag(rot), [amin,csmin,floor_diag])
        assert sep > 0 and floor_diag > max_e and angle > 0
        out[str(k)] = {
            "Vud_interval":[amin,amax], "Vub_Vcb_caps":[bmax,cmax],
            "folded_probability_five_max":t, "folded_probability_four_min":d,
            "quark_amplitude_five_ceiling_four_floor":[small_ceiling,other_floor],
            "reactor_amplitude":[flo,fhi],
            "small_and_large_modulus_roots":[[math.sqrt(tlo),math.sqrt(thi)],
                                                  [math.sqrt(1-thi),math.sqrt(1-tlo)]],
            "distance_to_quark_amplitude_bands_at_least":sep,
            "projective_frame_turn_sum_rad_at_least":angle,
            "alternative_fit_also_separate":small_ceiling < cfg["extra"]["effective_frequency_daya_bay_amplitude"][0]-k*cfg["extra"]["effective_frequency_daya_bay_amplitude"][1]
            and cfg["extra"]["effective_frequency_daya_bay_amplitude"][0]+k*cfg["extra"]["effective_frequency_daya_bay_amplitude"][1] < other_floor,
            "independent_two_quark_tags_plus_solar":{"quark_row_or_column_max_at_least":floor_diag,
                                                      "electron_row_max_at_most":max_e,
                                                      "diagonal_inequality_saturation":np.diag(rot).tolist()}}
    out["counterexample_to_rejecting_a_common_gauge_current"] = frame_witness(q, db, ju)
    return out


def same(x, y):
    if isinstance(x, dict):
        return x.keys() == y.keys() and all(same(x[k], y[k]) for k in x)
    if isinstance(x, list):
        return len(x) == len(y) and all(same(i,j) for i,j in zip(x,y))
    if isinstance(x, float):
        return math.isclose(x, y, rel_tol=1e-11, abs_tol=1e-14)
    return x == y


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = evaluate()
    path = ROOT / "results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(data, stream, indent=2)
            stream.write("\n")
    else:
        assert same(data, read(path))
    print("1384 OK")
