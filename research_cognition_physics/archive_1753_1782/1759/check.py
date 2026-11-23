"""Exact two-circle reconstruction at frozen left currents; finite VL doublet check."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def triangle(left, vector, axial):
    """R in a convention with L positive real; either orientation is allowed."""
    real = (vector**2 - axial**2) / (4 * left)
    square = ((vector + axial)**2 - 4 * left**2) * (4 * left**2 - (vector - axial)**2)
    assert square > -3e-14
    return complex(real, math.sqrt(max(0., square)) / (4 * left))


def arg_ratio(left, right):
    return float(np.angle((left - right) / (left + right), deg=True))


def charged_ckm(row):
    s13 = row["ub"]
    c13 = math.sqrt(1 - s13**2)
    s12, s23 = row["us"] / c13, row["cb"] / c13
    c12, c23 = math.sqrt(1 - s12**2), math.sqrt(1 - s23**2)
    gamma = math.radians(row["gamma_deg"])
    k = c12 * s23 * s13 / (s12 * c23)
    phase = np.exp(1j * (gamma + math.asin(k * math.sin(gamma))))
    return np.array([[c12*c13, s12*c13, s13*phase.conjugate()],
                     [-s12*c23-c12*s23*s13*phase, c12*c23-s12*s23*s13*phase, s23*c13],
                     [s12*s23-c12*c23*s13*phase, -c12*s23-s12*c23*s13*phase, c23*c13]])


def rotation(column):
    light, last = column[:3], column[3].real
    unitary = np.eye(4, dtype=complex)
    unitary[:3, :3] -= np.outer(light, light.conjugate()) / (1 + last)
    unitary[:3, 3], unitary[3, :3], unitary[3, 3] = light, -light.conjugate(), last
    return unitary


def embed(ckm, rd, rs, p):
    fu = p["right_up_doublet_weight"]
    a = np.array([math.sqrt(fu), 0, 0, math.sqrt(1-fu)])
    light = np.array([rd.conjugate(), rs.conjugate(), 0]) / a[0]
    b = np.append(light, math.sqrt(1-np.vdot(light, light).real))
    mu, md = p["illustrative_light_quark_masses_GeV"]
    md = np.append(md, p["upper_mass_B_GeV"])
    common = np.linalg.norm(md*b)
    mt = math.sqrt(common**2 - (mu[0]*a[0])**2) / a[3]
    mu = np.append(mu, mt)
    base = np.eye(4, dtype=complex)
    base[:3, :3] = ckm
    vl = rotation(mu*a/common) @ base @ rotation(md*b/common).conj().T
    vr = np.outer(a, b.conjugate())
    mass_up = np.diag(mu) @ rotation(a)
    mass_down = vl @ np.diag(md) @ rotation(b)
    res = np.max(abs(mass_up[:, 3]-mass_down[:, 3]))
    assert res < 1e-9 and np.max(abs(vl @ vl.conj().T-np.eye(4))) < 4e-15
    assert np.allclose(vr[0, :2], [rd, rs], atol=1e-15)
    left = vl[0, :2]
    lt = vl[2, 1].conjugate()*vl[2, 0]
    lc = vl[1, 1].conjugate()*vl[1, 0]
    return {"masses_T_B_GeV": [mt, md[3]], "common_bare_column_residual_GeV": res,
            "max_change_light_CKM": np.max(abs(vl[:3, :3]-ckm)), "abs_left_uB": abs(vl[0, 3]),
            "abs_Zsd_twice_isospin": abs(b[0]*b[1]), "delta_Qp": 2*fu-abs(b[0])**2,
            "delta_QCs": 188*fu-211*abs(b[0])**2,
            "top_lambda_re_im_and_charm_lambda_real": [lt.real, lt.imag, lc.real],
            "physical_nuclear_vector_K_vector_axial_ratio": [abs(left[0]+rd), abs(left[1]+rs),
                                                           abs((left[1]-rs)/(left[0]-rd))]}


def lower_endpoint(c, s, meas, width):
    # Minimize the two current moduli simultaneously under the three rate bands.
    (u, eu), (v, ev), (ratio, er) = meas
    ud = min(c, u+width*eu)
    ad = 2*c-ud
    amin = max(s, (ratio-width*er)*ad, 2*s-(v+width*ev))
    vec = max(v-width*ev, 2*s-amin)
    assert vec <= v+width*ev+3e-15 and amin <= (ratio+width*er)*ad+3e-15
    rs = triangle(s, vec, amin)
    rd = ud-c
    return rd, rs, [ud, vec, amin/ad]


def threshold(test, left, right):
    assert not test(left) and test(right)
    for _ in range(56):
        mid = (left+right)/2
        if test(mid):
            right = mid
        else:
            left = mid
    return right


def one_row(row, p, meas, caps, fu_cap, rare):
    V = charged_ckm(row)
    c, s = V[0, :2].real
    (u, du), (v, dv), (ratio, dr) = meas
    def gap(w):
        return v-w*dv+(ratio-w*dr)*(2*c-u-w*du)-2*s
    wreal = threshold(lambda w: gap(w) <= 0, 0, 1)
    # Neutron phase may have either sign in its two-marginal-error interval.
    d0, dwidth = p["neutron_phase_minus_180_deg"]
    phi = math.radians(abs(d0)+p["neutron_phase_marginal_width"]*dwidth)
    amin = 2*c-u
    amax = -u*math.cos(phi)+math.sqrt((2*c)**2-(u*math.sin(phi))**2)
    low_s, high_s = [triangle(s, v, ratio*a) for a in (amin, amax)]
    im_d_cap = u*amax*math.sin(phi)/(2*c)
    assert (ratio*amax)**2 < v**2+4*s**2
    dcoef, scoef = p["illustrative_PQ_matched_neutron_coefficients_e_cm"]
    edm = scoef*s*low_s.imag-dcoef*c*im_d_cap
    base_rd, base_rs, _ = lower_endpoint(c, s, meas, 0)
    def minimum(w):
        rd, rs, _ = lower_endpoint(c, s, meas, w)
        return abs(rd*rs)
    widths = [0, .25, .5, 1, 2, 2.5]
    floors = {}
    for w in widths:
        rd, rs, vals = lower_endpoint(c, s, meas, w)
        floors[str(w)] = {"Rd": rd, "Rs_re_im": [rs.real, rs.imag], "attaining_rates": vals,
                          "min_abs_RdRs": abs(rd*rs), "min_abs_Fsd_using_PV": abs(rd*rs)/fu_cap}
    # Standard left top/charm central at the *same* original CKM. For real negative
    # currents the new contribution aligns with the negative real SM amplitude.
    lam = rare["lambda_reference"]
    lt, lc = V[2, 1].conjugate()*V[2, 0], V[1, 1].conjugate()*V[1, 0]
    z0 = lt*rare["Xtop"][0]/lam**5 + lc.real/lam * rare["Pcharm"][0]
    assert z0.real < 0
    shell = caps["two_width_measured_annulus_total_amplitude"][1]
    h = caps["normalizations"]["h"]
    signed_Fcap = (math.sqrt(shell**2-z0.imag**2)+z0.real)/h
    assert signed_Fcap > 0
    cap_widths = {k: threshold(lambda w: minimum(w) <= b["right_W_product_cap"], 1, 2.8)
                  for k, b in caps["branches"].items()}
    signed_width = threshold(lambda w: minimum(w) <= signed_Fcap*fu_cap, 1, 2.8)
    # Concrete two-width charged fit and almost-decoupled joint case.
    def real_witness(w):
        rd, rs, vals = lower_endpoint(c, s, meas, w)
        assert abs(rs.imag) < 1e-7
        rs = rs.real
        bound = embed(V, rd, rs, p)
        fu = p["right_up_doublet_weight"]
        re, im, lcre = bound["top_lambda_re_im_and_charm_lambda_real"]
        z = complex(re, im)*rare["Xtop"][0]/lam**5 + lcre/lam*rare["Pcharm"][0]-h*rd*rs/fu
        br = rare["kplus_e11"][0]*(1+rare["delta_EM"])*abs(z)**2
        kl = rare["kL_e11"][0]*z.imag**2
        return {"Rd": rd, "Rs": rs, "rates_nominal": vals, "finite_doublet": bound,
                "central_tree_top_charm_plus_Z_BR_times_1e11": [br, kl]}
    pion, epi = p["pion_beta_experimental_normalization"]
    return {"left_ud_us": [c, s], "critical_real_charged_width": wreal,
            "neutron_two_width_exact_rate_centers": {
                "axial_d_interval": [amin, amax], "abs_ImRd_upper": im_d_cap,
                "abs_ImRs_interval": [low_s.imag, high_s.imag],
                "abs_strange_axial_vector_phase_deg": [abs(arg_ratio(s, low_s)), abs(arg_ratio(s, high_s))],
                "edm_benchmark_lower_e_cm": edm, "benchmark_over_EDM_limit": edm/p["edm_90_upper_e_cm"]},
            "sharp_charged_product": floors, "necessary_width_for_old_product_caps": cap_widths,
            "same_point_left_rare_baseline": [z0.real, z0.imag], "real_negative_Fsd_cap_at_central_matching": signed_Fcap,
            "necessary_width_for_positive_Zsd_at_central_matching": signed_width,
            "pion_beta_pull_fit_and_left_only": [(u-pion)/epi, (c-pion)/epi],
            "complex_center_doublet_witness": embed(V, base_rd, base_rs, p),
            "real_two_width_witness": real_witness(2), "nonzero_2p70_width": real_witness(2.70),
            "near_decoupled_2p8_width": real_witness(2.8),
            "pure_left_required_width": max(abs(c-u)/du, abs(s-v)/dv, abs(s/c-ratio)/dr)}


def main():
    p = read("inputs.json")
    raw = read(p["charged_measurements"])
    meas = [(raw[k][0], math.hypot(*raw[k][1:])) for k in ("holdout_nuclear", "holdout_Kl3", "holdout_axial")]
    parents = read(p["parent_model"])
    caps = read(p["neutral_caps"])
    fu_cap = read("../1460/results.json")["weak_charges_Fu_upper"]
    rare = read(p["neutral_normalization"])
    return {key: one_row(parents[key], p, meas, caps, fu_cap, rare) for key in p["parent_rows"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    def clean(obj):
        if isinstance(obj, dict):
            return {k: clean(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [clean(v) for v in obj]
        return float(f"{obj:.11g}") if isinstance(obj, (float, np.floating)) else obj
    result = json.dumps(clean(main()), ensure_ascii=False, indent=2, allow_nan=False)+"\n"
    destination = HERE/"results.json"
    if args.write:
        destination.write_text(result, encoding="utf-8")
    elif destination.read_text(encoding="utf-8") != result:
        raise SystemExit("1759 result mismatch")
    else:
        print("1759 OK: charged phase, sharp product, neutron/pion, full 4x4 and rare K")
