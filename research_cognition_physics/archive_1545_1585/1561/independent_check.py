"""Signed cycle invariant, compressed CKM data and a spectral competitor."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def ckm(theta12, theta13, theta23, phase):
    su, sb, sc = np.sin([theta12, theta13, theta23])
    cu, cb, cc = np.cos([theta12, theta13, theta23])
    z = np.exp(1j*phase)
    return np.array([[cu*cb, su*cb, sb/z],
                     [-su*cc-cu*sc*sb*z, cu*cc-su*sc*sb*z, sc*cb],
                     [su*sc-cu*cc*sb*z, -cu*sc-su*cc*sb*z, cc*cb]])


def spectrum_factor(masses):
    d, s, b = masses
    return (1-(d/s)**2)*(1-(s/b)**2)*((b/d)**2-1)


def cycle(gram):
    triple = gram[0,1]*gram[1,2]*gram[2,0]
    return float(triple.imag/np.linalg.det(gram).real)


def compute():
    cfg = read(HERE/"inputs.json")
    high = read(HERE/cfg["charged_target"])["high_scale"]
    original = read(HERE/cfg["shared_tensors"])
    two = cfg["widths"]
    down_bounds = np.array(high["down"]) @ [[1,1],[-two,two]]
    d, b = down_bounds[0,0], down_bounds[2,1]
    r = b/d
    best_s = math.sqrt(b*d)
    phi_max = (1-1/r)**2*(r**2-1)
    assert math.isclose(phi_max, spectrum_factor([d, best_s, b]))
    rows = high["angles_theta12_theta13_theta23_delta"]
    th13, th23 = [x[0]+two*x[1] for x in rows[1:3]]
    max_j = math.sin(th13)*math.cos(th13)**2*math.sin(th23)*math.cos(th23)/2
    x, y = math.sin(th13), math.sin(th23)
    max_jac = phi_max*max_j
    edge = ckm(math.pi/4, th13, th23, math.pi/2)
    gram_edge = (edge*np.square([d,best_s,b]))@edge.conj().T
    assert math.isclose(abs(cycle(gram_edge)), max_jac, rel_tol=2e-11)
    cfg_tree = read(HERE/"../1381/inputs.json")
    q0, *q_err = cfg_tree["ub_over_cb"]
    q_hi = q0 + two*sum(q_err)
    output = dict(target_two_widths=dict(lightest_low=d, heaviest_high=b,
        unconstrained_optimal_middle=best_s, cp_invariant_ceiling=max_jac,
        Vub_upper=x, Vcb_upper_bound=y, tree_Vub_over_Vcb_upper=q_hi), branches={})
    for key in cfg["keys"]:
        vals = dict(zip(original["labels"], original["candidates"][key]["values"]))
        masses = np.array([vals[k] for k in ("yd","ys","yb")])
        angles = [vals[k] for k in ("theta12_q","theta13_q","theta23_q","delta_q")]
        v = ckm(*angles)
        j = (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
        gram = (v*masses**2)@v.conj().T
        direct, spectral = cycle(gram), j*spectrum_factor(masses)
        assert math.isclose(abs(direct),abs(spectral),rel_tol=1e-10)
        # A finite change of every complex up-aligned endpoint preserves the cycle.
        gain = np.diag(np.array([1.8,.4,2.7])*np.exp(1j*np.array([.5,-.75,1.3])))
        assert math.isclose(direct,cycle(gain@gram@gain.conj().T),rel_tol=1e-10)
        needed_product = 2*abs(direct)/phi_max
        need_ratio = needed_product/y**2
        # One down-aligned gain repairs only the lightest mass and preserves all CKM entries.
        repaired = np.array([high["down"][0][0], *masses[1:]])
        down_gain = (v*(repaired/masses))@v.conj().T
        target = (v*repaired**2)@v.conj().T
        assert np.allclose(down_gain@gram@down_gain.conj().T, target,rtol=1e-10,atol=1e-20)
        output["branches"][key] = dict(signed_invariant=direct, J=float(j),
            smallest_required_Vub_times_Vcb=needed_product,
            Vub_required_given_Vcb_upper=needed_product/y,
            Vub_over_Vcb_required_given_Vcb_upper=need_ratio,
            ratio_transport_required=need_ratio/q_hi,
            residual_against_loose_target=abs(direct)-max_jac,
            down_aligned_gain=repaired[0]/masses[0],repaired_signed_invariant=cycle(target))
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = compute()
    target = HERE/"independent_results.json"
    if args.write_new:
        with target.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
    else:
        prev = read(target)
        assert prev.keys() == result.keys()
        for key in result["target_two_widths"]:
            assert np.isclose(prev["target_two_widths"][key], result["target_two_widths"][key])
        for branch, row in result["branches"].items():
            for key, value in row.items():
                assert math.isclose(value, prev["branches"][branch][key], rel_tol=1e-10, abs_tol=1e-13)
        print("1561 independent signed weak cycle, uncalibrated angles and spectrum verified")
