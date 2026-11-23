"""Direct three-family rotations, unitarity-triangle reconstruction and a CP holdout."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rotate(r, t, vcb, short, phase_sign):
    c0, s0 = math.cos(math.pi/14), math.sin(math.pi/14)
    p = r*r/(1+r*r)
    q = t*t/(1+t*t)
    cp = (q-p*c0*c0-(1-p)*s0*s0)/(2*c0*s0*math.sqrt(p*(1-p)))
    assert abs(cp) <= 1
    ray = np.array([math.sqrt(p), math.sqrt(1-p)*complex(cp, phase_sign*math.sqrt(1-cp*cp))])
    sine = vcb * math.sqrt(1+r*r)
    cosine = math.sqrt(1-sine*sine) * (1 if short else -1)
    rotor = np.eye(3, dtype=complex)
    rotor[:2, :2] += (cosine-1)*np.outer(ray, ray.conj())
    rotor[:2, 2] = sine*ray
    rotor[2, :2] = -sine*ray.conj()
    rotor[2, 2] = cosine
    seed = np.array([[c0, -s0, 0], [s0, c0, 0], [0, 0, 1]])
    v = rotor @ seed
    np.testing.assert_allclose(v @ v.T.conj(), np.eye(3), atol=2e-15)
    np.testing.assert_allclose(abs(v[0,2]/v[1,2]), r, atol=2e-15)
    np.testing.assert_allclose(abs(v[2,0]/v[2,1]), t, atol=2e-15)
    return v


def readout(v):
    u, c, t = v[:, 0]*v[:, 2].conj()
    gamma = np.angle(-u/c)
    beta = np.angle(-c/t)
    weak_s = math.sin(2*beta)
    cos_beta = (abs(c)**2+abs(t)**2-abs(u)**2)/(2*abs(c*t))
    assert abs(math.cos(beta)-cos_beta) < 2e-12
    assert abs(u+c+t) < 1e-15
    return {"ud": abs(v[0,0]), "us": abs(v[0,1]), "ub": abs(v[0,2]),
            "us_over_ud": abs(v[0,1]/v[0,0]), "gamma_degrees": math.degrees(gamma),
            "beta_degrees": math.degrees(beta), "sin_2beta": weak_s,
            "J_magnitude": abs((v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj()).imag)}


def calc():
    old = read(HERE.parent / "1381" / "inputs.json")
    belle = read(HERE.parent / "1412" / "inputs.json")["BelleII_Vcb_BGL"]
    old_beta_high = old["ud_beta"][0]+2*old["ud_beta"][1]
    old_k_high = (old["K_product"][0]+2*old["K_product"][1])/(old["K_form"][0]-2*old["K_form"][1])
    vmax = belle[0]+2*sum(belle[1:])
    ratio_max = old["ub_over_cb"][0]+2*sum(old["ub_over_cb"][1:])
    sine_max = vmax*math.sqrt(1+ratio_max**2)
    dmax = ratio_max*vmax*sine_max/(1+math.sqrt(1-sine_max*sine_max))
    assert math.cos(math.pi/14)-dmax > old_beta_high
    s_obs, s_error = .717, math.hypot(.013,.008)
    cases = {
        "central": (.083, .2052, .04057),
        "old_nuclear_kaon_boundary": (.094, .2095312, .0442),
        "updated_nuclear_and_KS_example": (.09, .2118, .04057)}
    result = {}
    for name, (r,t,b) in cases.items():
        both = []
        for short in (True, False):
            a = readout(rotate(r, t, b, short, 1))
            rev = readout(rotate(r, t, b, short, -1))
            assert abs(a["sin_2beta"]+rev["sin_2beta"]) < 2e-13
            both.append(a if a["gamma_degrees"] > 0 else rev)
        result[name] = {"r_t_b": [r,t,b], "short_arc": both[0], "long_arc": both[1]}
    edge = result["old_nuclear_kaon_boundary"]["long_arc"]
    assert edge["ud"] < old_beta_high and edge["us"] < old_k_high
    assert abs(edge["sin_2beta"]-s_obs) < 2*s_error
    modern = result["updated_nuclear_and_KS_example"]["long_arc"]
    assert abs(modern["sin_2beta"]-s_obs) < s_error
    assert abs(modern["ud"]-.97359) < 2*.00033
    assert 62.6 < modern["gamma_degrees"] < 74.9
    # Illustrative matching inputs within the cited separate lattice intervals.
    assert abs(.967*modern["us"]-.2170) < 2*.0009
    assert abs(1.192*modern["us_over_ud"]-.27599) < 2*(.00029+.00024)
    new_nuclear_high = .97359+2*.00033
    floors = [math.sqrt(1-u*u-old_k_high**2) for u in (old_beta_high, new_nuclear_high)]
    assert floors[1] > ratio_max*vmax > floors[0]
    return {"historical_K_upper": old_k_high, "old_and_revised_nuclear_with_old_K_min_ub": floors,
            "ub_max_in_inputs": ratio_max*vmax, "short_arc_tree_only_max_ud_displacement": dmax,
            "S_psi_KS_direct": [s_obs, .013, .008], "two_component_quadrature_error": s_error,
            "cases": result}


def canonical(obj):
    if isinstance(obj, dict):
        return {key: canonical(v) for key, v in obj.items()}
    if isinstance(obj, list):
        return [canonical(v) for v in obj]
    return round(float(obj), 11)


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    value = canonical(calc())
    target = HERE / "independent_results.json"
    if args.write:
        target.write_text(json.dumps(value, indent=2)+"\n", encoding="utf-8")
    else:
        assert read(target) == value
    print("1614 independent: exact CKM, delta-B holdout and row-norm budget verified")


if __name__ == "__main__":
    main()
