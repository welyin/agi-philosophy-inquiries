"""Compare a tagged right-handed down line with any rank-two neutrino null line."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def interval(row, widths):
    step = widths * sum(row[1:])
    return row[0] - step, row[0] + step


def ckm(s12, s13, s23, delta):
    c12, c13, c23 = (math.sqrt(1 - s*s) for s in (s12, s13, s23))
    z = np.exp(1j * delta)
    return np.array([[c12*c13, s12*c13, s13/z],
                     [-s12*c23-c12*s23*s13*z, c12*c23-s12*s23*s13*z, s23*c13],
                     [s12*s23-c12*c23*s13*z, -c12*s23-s12*c23*s13*z, c23*c13]])


def right_tags(tree, belle, mass, gamma, widths):
    k, f = [interval(tree[key], widths) for key in ("K_product", "K_form")]
    cb, r = interval(belle, widths), interval(tree["ub_over_cb"], widths)
    us, ub = (k[0]/f[1], k[1]/f[0]), (cb[0]*r[0], cb[1]*r[1])
    c13lo = math.sqrt(1-ub[1]**2)
    s12hi, s23hi = us[1]/c13lo, cb[1]/c13lo
    c12lo = math.sqrt(1-s12hi**2)
    c23lo = math.sqrt(1-s23hi**2)
    angles = (gamma["central_degrees"]-widths*gamma["lower_error"]-gamma["additional_degrees"],
              gamma["central_degrees"]+widths*gamma["upper_error"]+gamma["additional_degrees"])
    delta_pad = math.asin(s23hi*ub[1]/(us[0]*c23lo))
    min_sine = min(abs(math.sin(math.radians(angles[0])-delta_pad)),
                   abs(math.sin(math.radians(angles[1])+delta_pad)))
    assert delta_pad < math.radians(angles[0]) < math.radians(angles[1]) < math.pi/2-delta_pad
    bound = {}
    for key, min_den in (("all_quark_phases", us[0]-r[1]),
                         ("also_charged_B_gamma", max(us[0]-r[1], us[0]*min_sine))):
        assert min_den > 0
        # Last row: Vtd/s23 = s12-c12*t*exp(i delta), t=c23*c13*r.
        d_lower = (c12lo*c13lo, us[0]*c23lo-s23hi*ub[1], cb[0]*min_den)
        s_over_d = (us[1]/d_lower[0], 1/d_lower[1],
                     (1+s12hi*r[1])/min_den)
        b_over_d = (ub[1]/d_lower[0], cb[1]/d_lower[1], 1/d_lower[2])
        yd, ys, yb = (interval(row, widths) for row in mass)
        bound[key] = {row: 1/math.sqrt(1+(sd*yd[1]/ys[0])**2+(bd*yd[1]/yb[0])**2)
                      for row, sd, bd in zip(("u","c","t"), s_over_d, b_over_d)}
    return bound


def neutrino_bounds(order, nu, atm):
    _, x0, x1 = nu["s12"]
    _, z0, z1 = nu["s13"]
    _, v0, v1 = atm
    if order == "IO":
        moduli = np.sqrt([z1, v1*(1-z0), (1-v0)*(1-z0)])
        first_nonzero = math.sqrt(nu["atm_abs"][1]-nu["solar"][2])
    else:
        moduli = [math.sqrt((1-x0)*(1-z0)),
                  math.sqrt(x1*(1-v0))+math.sqrt((1-x0)*v1*z1),
                  math.sqrt(x1*v1)+math.sqrt((1-x0)*(1-v0)*z1)]
        first_nonzero = math.sqrt(nu["solar"][1])
    return float(max(moduli)), first_nonzero


def finite_witness(V, y, nu, atm, order):
    # One relative family phase convention is realized in the charged mass tensors.
    lepton_delta = math.radians(-75 if order == "NO" else -90)
    unitary = ckm(math.sqrt(nu["s12"][0]), math.sqrt(nu["s13"][0]), math.sqrt(atm[0]), lepton_delta)
    i = 0 if order == "NO" else 2
    n = unitary[:, i]
    phases = np.exp(-1j*np.angle(V[2]*n))
    U = phases[:, None]*unitary
    n = U[:, i]
    masses = (np.sqrt([0, nu["solar"][0], nu["atm_abs"][0]]) if order == "NO" else
              np.sqrt([nu["atm_abs"][0]-nu["solar"][0], nu["atm_abs"][0], 0]))
    N = (U*masses) @ U.T
    Y = V*y
    left_column = Y@n
    left_column /= np.linalg.norm(left_column)
    w = V[2].conj()/y
    w /= np.linalg.norm(w)
    cosine = float(abs(np.vdot(w, n)))
    # Rotate w onto n in their span. This matches the two bases if the relative
    # charged mass rotation is freed; it preserves all singular values.
    n1 = n*np.exp(-1j*np.angle(np.vdot(w,n)))
    sine = math.sqrt(1-cosine*cosine)
    v = (n1-cosine*w)/sine
    rot = (np.eye(3)+(cosine-1)*(np.outer(w,w.conj())+np.outer(v,v.conj()))
           +sine*(np.outer(v,w.conj())-np.outer(w,v.conj())))
    assert np.linalg.norm(rot.conj().T@rot-np.eye(3)) < 1e-12
    assert np.linalg.norm(rot@w-n1) < 1e-12 and np.linalg.norm(N@n.conj()) < 1e-14
    Ue = rot.conj().T
    common_neutral = Ue @ N @ Ue.T
    assert np.linalg.norm(common_neutral @ w.conj()) < 1e-14
    assert np.linalg.norm(Ue.conj().T@(Ue@U)-U) < 1e-12
    assert np.max(abs(np.linalg.svd(N, compute_uv=False)-np.sort(masses)[::-1])) < 1e-13
    return {"quark_to_neutral_angle_deg_with_pure_left_t": math.degrees(math.acos(cosine)),
            "left_t_impurity_deg_if_bases_align_and_null_is_exact": math.degrees(math.acos(abs(left_column[2]))),
            "left_t_leakage_probability_if_bases_align": float(1-abs(left_column[2])**2)}


def calculate():
    cfg = read(HERE/"inputs.json")
    def old(key):
        return read(HERE.parent/cfg[key])
    tree, belle = old("tree"), old("belle")["BelleII_Vcb_BGL"]
    gamma, mass = old("gamma")["chargedB"], old("down_and_scheme")["yukawa"]["down"]
    atm, osc = old("atmospheric_angles")["theta23_sin2"], old("oscillations")
    tags = right_tags(tree, belle, mass, gamma, cfg["charged_widths"])
    y = np.array([m[0] for m in mass])
    cb, r = belle[0], tree["ub_over_cb"][0]
    us = tree["K_product"][0]/tree["K_form"][0]
    s13 = cb*r
    c13 = math.sqrt(1-s13**2)
    s12, s23 = us/c13, cb/c13
    gam = math.radians(gamma["central_degrees"])
    delta = gam+math.asin(math.sqrt(1-s12*s12)*s23*s13*math.sin(gam)/(s12*math.sqrt(1-s23*s23)))
    V = ckm(s12, s13, s23, delta)
    outcome = {"lower_abs_d_amplitude_by_up_port":tags, "orders":{}}
    for order in ("NO", "IO"):
        largest, mass_gap = neutrino_bounds(order, osc[order], atm[order])
        results = {}
        for key, rows in tags.items():
            theta = max(0, math.acos(largest)-math.acos(min(rows.values())))
            angles = cfg["matching_angle_examples_deg"]
            perturbations = {str(a): 1000*mass_gap*math.sin(max(0, theta-math.radians(a))) for a in angles}
            condition_inverse = interval(mass[0], 2)[0]/interval(mass[2], 2)[1]
            results[key] = {"projective_gap_deg": math.degrees(theta),
                            "source_error_lower_meV_vs_matching_degrees":perturbations,
                            "left_leak_angle_deg_necessary_if_only_purity_is_relaxed":
                            math.degrees(math.asin(condition_inverse*math.sin(theta)))}
        assert results["all_quark_phases"]["projective_gap_deg"] > 1
        outcome["orders"][order] = {"largest_null_component_upper":largest,
                                     "smallest_nonzero_mass_lower_meV":1000*mass_gap,
                                     "necessary":results,
                                     "central_counterexamples":finite_witness(V,y,osc[order],atm[order],order)}
    return outcome


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    dest = HERE/"results.json"
    if args.check:
        def same(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(same(v,b[k]) for k,v in a.items())
            return math.isclose(a,b,rel_tol=2e-10,abs_tol=1e-12)
        assert same(result,read(dest))
    else:
        with dest.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    print({order:result["orders"][order]["necessary"]["all_quark_phases"]["projective_gap_deg"] for order in ("NO","IO")})
