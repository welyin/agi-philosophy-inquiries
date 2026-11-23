"""Reproduce the tagged-plane completions and the independent charged-lepton holdout."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
IN = json.loads((HERE / "inputs.json").read_text(encoding="utf-8-sig"))
NAMES = ("u", "c", "t")


def load(name):
    return json.loads((HERE / IN[name]).read_text(encoding="utf-8-sig"))


def mixing(angles):
    t12, t13, t23, phase = angles
    s, r, t = np.sin([t12, t13, t23])
    c, v, w = np.cos([t12, t13, t23])
    z = np.exp(1j * phase)
    return np.array([[c*v, s*v, r*z.conjugate()],
                     [-s*w-c*t*r*z, c*w-s*t*r*z, t*v],
                     [s*t-c*w*r*z, -c*t-s*w*r*z, w*v]])


def fixed_plane(d, leptons, row):
    D = np.diag(d)
    b = row.conj() / d
    b = b / np.linalg.norm(b)
    P = np.outer(b, b.conj())
    base = D @ (np.eye(3) - P)
    roots, U = np.linalg.eigh(base @ base.conj().T)
    roots[0] = 0.0
    s1, s2 = np.sqrt(roots[1:])
    lo = max(leptons[0]/s1, leptons[1]/s2)
    hi = min(leptons[1]/s1, leptons[2]/s2)
    assert lo < 1 < hi
    weights = [-np.prod(a - leptons**2) / np.prod(a - np.delete(roots, i))
               for i, a in enumerate(roots)]
    assert min(weights) > 0
    col = U @ np.sqrt(weights)
    E = base + np.outer(col, b.conj())
    assert np.allclose(np.linalg.svd(E, compute_uv=False)[::-1], leptons,
                       rtol=1e-8, atol=1e-15)
    assert np.linalg.norm((E-D)@(np.eye(3)-P)) < 1e-15
    return {"compressed_singular": [s1, s2], "gain_interval": [lo, hi],
            "k1_leptons": np.linalg.svd(E, compute_uv=False)[::-1].tolist(),
            "squared_mass_weights_of_added_column": weights}


def bracket(fn, left, right):
    lv = fn(left)
    if lv * fn(right) >= 0:
        raise ValueError("Root not bracketed")
    for _ in range(70):
        mid = (left + right)/2
        value = fn(mid)
        if value == 0:
            return mid
        if lv * value > 0:
            left, lv = mid, value
        else:
            right = mid
    return (left + right)/2


def endpoint_preserving(d, mu, tau, row):
    # Calibrate only the muon and tau. The electron is an output.
    a, s, _ = (d / d[2])**2
    q = (mu/tau)**2
    p = np.abs(row)**2

    def secular(x):
        light = p[0]/((1-x*a)*(q-x*a)) + p[1]/((1-x*s)*(q-x*s))
        return (x-q)*(1-x)*light - p[2]

    if secular(0.5) <= 0:
        return []  # A full no-solution proof, when used, is in uniform_bounds.
    out = []
    for ends in [(q, 0.5), (0.5, 1.)]:
        x = bracket(secular, *ends)
        h = 1 / np.sum(p/(1-x*np.array([a,s,1.])))
        k = math.sqrt(x)*tau/d[2]
        N = 1 / np.linalg.norm(row.conj()/d)
        gain = math.sqrt(1 + h*(tau/(k*N))**2)
        b = N*row.conj()/d
        P = np.outer(b, b.conj())
        E = k*np.diag(d)@(np.eye(3)+(gain-1)*P)
        _, singular, right = np.linalg.svd(E)
        assert np.allclose(singular[:2], [tau,mu], rtol=1e-7)
        q_e = (x**3*a*s + h*x*x*(p[0]*s+p[1]*a+p[2]*a*s))/q
        assert abs(math.sqrt(q_e)*tau - singular[2]) < 1e-12
        out.append({"x": x, "h": h, "k": k, "column_amplitude_gain": gain,
                    "electron_prediction": singular[2],
                    "tag_tau_probability": abs(right[0] @ b)**2})
    return out


def uniform_bounds(m, ang):
    width = IN["margin_widths"]
    lo, hi = m[:,0]-width*m[:,1], m[:,0]+width*m[:,1]
    left, right = ang[:,0]-width*ang[:,1], ang[:,0]+width*ang[:,1]
    s12, s13, s23 = np.sin(left[:3])
    c12, c13, c23 = np.cos(right[:3])
    S12, S13, S23 = np.sin(right[:3])
    pd_min = np.array([(c12*c13)**2, (s12*c23-S23*S13)**2,
                       (s12*s23-S13)**2])
    assert min(s12*c23-S23*S13, s12*s23-S13) > 0
    bottom_tag_cap = (hi[0]/lo[2])**2 / min(pd_min)
    x_cap = 1/(1-bottom_tag_cap)
    qlo, qhi = (lo[4]/hi[5])**2, (hi[4]/lo[5])**2
    slo, shi = (lo[1]/hi[2])**2, (hi[1]/lo[2])**2
    electron_upper = hi[3]/lo[5]
    assert x_cap*shi < qlo  # negative update cannot give the muon.
    assert lo[0] > hi[3] and hi[1] < lo[4] and hi[2] < lo[5]
    assert lo[2]*math.sqrt(1-bottom_tag_cap) > hi[4]
    ptb_lower = (c13*c23)**2
    top_secular_upper = (1-ptb_lower)/(4*(1-shi)*(qlo-shi))
    assert top_secular_upper < ptb_lower
    holdout = {}
    for j, name in enumerate(NAMES[:2]):
        l2 = pd_min[j]*slo
        # sqrt(q(1-q)l2) on [q,1/2], and h >= q-x*s on [1/2,1].
        lower_low = math.sqrt(qlo*(1-qhi)*l2)
        lower_high = 0.5*math.sqrt((1-shi/qlo)*l2)
        holdout[name] = {"min_e_over_tau_for_x_below_half": lower_low,
                         "min_e_over_tau_for_x_above_half": lower_high}
    assert holdout["u"]["min_e_over_tau_for_x_below_half"] > electron_upper
    assert min(h["min_e_over_tau_for_x_above_half"] for h in holdout.values()) > electron_upper
    xmax_c = electron_upper * math.sqrt(2*qhi/(pd_min[1]*slo))
    io_max = IN["neutral_zero_column_maximum"]["IO"]
    k1_tau_probability = 1-(hi[2]/lo[5])**2
    k_io_min = lo[5]/hi[2]*math.sqrt(1-io_max**2)
    k1_io_angle = math.degrees(math.acos(io_max)-math.acos(math.sqrt(k1_tau_probability)))
    assert k_io_min > 1 and k1_io_angle > 0
    pmin_tau = 1-xmax_c
    angles = {order: math.degrees(math.acos(maximum)
                                  - math.acos(math.sqrt(pmin_tau)))
              for order, maximum in IN["neutral_zero_column_maximum"].items()}
    return {"bottom_probability_max_on_tag_all": bottom_tag_cap,
            "absolute_x_max_same_plane": x_cap, "q_range":[qlo,qhi],
            "down_strange_squared_range":[slo,shi],
            "top_left_bottom_probability_lower":ptb_lower,
            "top_secular_rhs_upper":top_secular_upper,
            "electron_over_tau_data_upper":electron_upper,"electron_bounds":holdout,
            "charm_x_upper_after_e":xmax_c,"charm_k_upper_after_e":math.sqrt(xmax_c)*hi[5]/lo[2],
            "protected_tau_probability_lower":pmin_tau,
            "minimum_ray_transfer_degrees":angles,
            "general_k1_tau_probability_min":k1_tau_probability,
            "general_IO_k_min_zero_transfer":k_io_min,
            "general_k1_IO_minimum_transfer_degrees":k1_io_angle}


def compute():
    m = np.array(load("down_and_leptons_file")["scales"]["1e16"])
    ang = np.array(load("high_ckm_file")["high_scale"]["angles_theta12_theta13_theta23_delta"])
    d, light = m[:3,0], m[3:,0]
    V = mixing(ang[:,0])
    plane = {j: fixed_plane(d,light,row) for j,row in zip(NAMES,V)}
    strong = {j: endpoint_preserving(d,*light[1:],row) for j,row in zip(NAMES,V)}
    # Existence only: choose b Yukawa at its upper edge. Calibrate strange mass to e.
    witness_d = d.copy()
    witness_d[2] += IN["margin_widths"]*m[2,1]

    def difference(strange):
        test = witness_d.copy()
        test[1] = strange
        return endpoint_preserving(test,*light[1:],V[1])[0]["electron_prediction"]-light[0]

    width = IN["margin_widths"]*m[1,1]
    witness_d[1] = bracket(difference,d[1]-width,d[1]+width)
    witness = endpoint_preserving(witness_d,*light[1:],V[1])[0]
    return {"comparison_scale_GeV":IN["scale_GeV"],"central_fixed_input":plane,
            "central_two_mass_calibration":strong,
            "charm_existence_witness":{"down":witness_d.tolist(),**witness},
            "all_phases_and_closed_two_width_box":uniform_bounds(m,ang)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check = parser.parse_args().check
    rendered = json.dumps(compute(), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    path = HERE/"results.json"
    if check:
        assert path.read_text(encoding="utf-8-sig") == rendered
        print("1601: fixed-plane completions, electron holdout and two-ordering transfer checked")
    else:
        path.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
