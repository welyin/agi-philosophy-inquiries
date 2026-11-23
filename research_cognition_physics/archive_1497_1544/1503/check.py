"""Replay the two-phase texture, the physical PQ norm, and strong/flavour/photon holds."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def read(p):
    return json.loads((BASE/p).read_text(encoding="utf-8-sig"))


def row_rotation(t):
    cosine = 1/math.sqrt(1+t*t)
    return np.array([[cosine, t*cosine, 0], [-t*cosine, cosine, 0], [0, 0, 1.]])


def pq_strength(vev, scale):
    wu, wd = vev["wu_GeV"], vev["wd_GeV"]
    sum_w2 = vev["v_GeV"]**2/2
    vcharges = np.array([0., 1., -1.])
    w2 = np.array([sum_w2-wu**2-wd**2, wu**2, wd**2])
    mean = np.dot(w2, vcharges)/sum_w2
    return math.sqrt(2*(np.dot(w2, (vcharges-mean)**2)
                     +scale**2*(vev["xu_GeV"]**2+vev["xd_GeV"]**2)))


def kaon_width(mass, mk, mpi, form, fv):
    pion, ax = (mpi/mk)**2, (mass/mk)**2
    assert mass+mpi < mk
    phase_space = math.sqrt((1-pion-ax)**2-4*pion*ax)
    return mk**3/(16*math.pi*fv**2)*(1-pion)**2*phase_space*form**2


def solve():
    cfg = read("inputs.json")
    texture = read(cfg["texture"])
    y = read("../1386/inputs.json")["yukawa"]
    old_photon = read(cfg["old_photon"])["qcd"]
    recast = read(cfg["kaon_search"])["massless_reinterpretation"]
    kaon, edm = read(cfg["kaon_constants"]), read(cfg["edm"])["edm"]
    z = texture["z3"]
    charge = {"Q1": -1, "Q23": 0, "uR": 0, "dR": 0, "H0": 0,
              "Hu": 1, "Hd": -1, "Xu": 1, "Xd": -1, "A": 0}
    assert all((z[key]+q)%3 == 0 for key,q in charge.items())
    assert -charge["Q1"]-charge["Hu"]+charge["uR"] == 0
    assert -charge["Q1"]+charge["Hd"]+charge["dR"] == 0
    assert charge["Xu"]+charge["Xd"] == 0
    anomaly = 2*charge["Q1"]
    electric_anomaly = 2*3*charge["Q1"]*((2/3)**2+(-1/3)**2)
    w = texture["witness"]
    td = math.sqrt(w["down_yukawa"]/w["strange_yukawa"])
    tu = w["up_tangent"]
    ou = row_rotation(tu)
    s23 = w["Vcb"]*math.sqrt(1+tu**2)
    c23 = math.sqrt(1-s23**2)
    od = np.array([[1., 0, 0], [0, c23, s23], [0, -s23, c23]])@row_rotation(td)
    qs = np.diag([-1., 0, 0])
    kd = (od.T@qs@od)/abs(anomaly)
    ku = (ou.T@qs@ou)/abs(anomaly)
    yu = ou@np.diag([-y["up"][0][0], y["up"][1][0], y["up"][2][0]])@ou.T
    yd = od@np.diag([-w["down_yukawa"], w["strange_yukawa"], y["down"][2][0]])@od.T
    # Trace the two CP vacua and multiple bare angles; the same light masses and CKM persist.
    max_error = 0.
    j_values = []
    for sign in (1, -1):
        current_j = []
        for bare in (0., .83, math.pi):
            for t in (.17, .58, .5*(bare+sign*math.pi/2)):
                up_phase = np.diag([np.exp(-1j*t), 1, 1])
                dn_phase = np.diag([np.exp(1j*(sign*math.pi/2-t)), 1, 1])
                ck = (up_phase@ou).conj().T@(dn_phase@od)
                m_arg = np.angle(np.linalg.det(up_phase@yu)*np.linalg.det(dn_phase@yd))
                residual = abs(np.exp(1j*m_arg)-np.exp(1j*(sign*math.pi/2-2*t)))
                max_error = max(max_error, residual)
                if math.isclose(2*t,bare+sign*math.pi/2):
                    assert abs(np.exp(1j*(bare+m_arg))-1) < 1e-12
                current_j.append((ck[0,0]*ck[1,1]*ck[0,1].conjugate()*ck[1,0].conjugate()).imag)
        assert max(current_j)-min(current_j) < 1e-15
        j_values.append(current_j[0])
    assert max_error < 1e-12 and j_values[0]*j_values[1] < 0
    root, err1, err2 = cfg["chisus_4throot_MeV"]
    n = cfg["window_widths"]
    chiq = np.array([root-n*(err1+err2),root,root+n*(err1+err2)])*1e-3
    anomaly_ratio = electric_anomaly/anomaly
    alpha = 1/cfg["alpha_inverse"]
    slopes = {}
    for label, (c, err) in {"lattice": cfg["photon_lattice"],
                    "historical": [old_photon["photon_subtraction"],old_photon["photon_subtraction_width"]]}.items():
        assert c-n*err > anomaly_ratio
        slopes[label] = [alpha/(2*np.pi)*(c-n*err-anomaly_ratio)/chiq[2]**2,
                         alpha/(2*np.pi)*(c-anomaly_ratio)/chiq[1]**2,
                         alpha/(2*np.pi)*(c+n*err-anomaly_ratio)/chiq[0]**2]
    ratio_lo = math.sqrt((y["down"][0][0]-n*y["down"][0][1])/(y["down"][1][0]+n*y["down"][1][1]))
    ratio_hi = math.sqrt((y["down"][0][0]+n*y["down"][0][1])/(y["down"][1][0]-n*y["down"][1][1]))
    fv_limit = recast["FV_GeV_lower"]
    kbox = [t/(2*(1+t*t)) for t in [ratio_lo,ratio_hi]]
    fa_min = fv_limit*kbox[0]/2
    fit_fa = pq_strength(texture["vev_example"], 1)/2
    light_mass = chiq[1]**2/fit_fa
    fv = 2*fit_fa/abs(kd[0,1])
    width = kaon_width(light_mass,kaon["MK_GeV"],kaon["Mpi_GeV"],cfg["kaon_fplus0"][0],fv)
    gamma_k = kaon["hbar_GeV_s"]/kaon["tauK_s"]
    big = pq_strength(texture["vev_example"], 5e9)/2
    assert big > fv_limit*kbox[1]/2 and width/gamma_k > 1e5
    edm_slope = edm["abs_slope_e_fm"]
    theta_cap = edm["limit_90_ecm"]/(edm["fm_in_cm"]*(edm_slope[0]-sum(edm_slope[1:])))
    pin_bounds = {}
    for deg in cfg["weak_phase_tolerance_degrees"]:
        phase_margin = math.pi/4-1.5*(math.radians(deg)+theta_cap)
        bound = 2*(1+cfg["relative_torque_tolerance"])*theta_cap/(3*math.sin(phase_margin))
        pin_bounds[str(deg)] = {"amplitude_over_chi_at_most": bound,
                                "amplitude_GeV4_using_chi_upper": bound*chiq[2]**4,
                                "kappa_at_30_GeV_at_most": bound*chiq[2]**4/(2*30.**4)}
    return {"PQ_color_N": anomaly, "PQ_EM_E": electric_anomaly, "E_over_N": anomaly_ratio,
            "Z3_color_instanton_exponent_mod3": (2*z["Q1"])%3, "walls_per_fixed_weak_branch": 2,
            "weak_J_on_two_orbits": j_values, "largest_determinant_phase_error": max_error,
            "low_energy_absolute_photon_slope_GeV_minus2_lo_center_hi": slopes,
            "kV_abs_sd_two_widths": kbox, "witness_abs_kV_uc_sd": [abs(ku[0,1]), abs(kd[0,1])],
            "kaon_conditioned_fa_lower_GeV": {str(r):r*fa_min for r in cfg["kaon_source_transfer_examples"]},
            "unit_transfer_m_upper_microeV": 1e15*chiq[2]**2/fa_min,
            "unit_transfer_lattice_g_upper_GeV_neg1": alpha/(2*np.pi)*
                (cfg["photon_lattice"][0]+n*cfg["photon_lattice"][1]-anomaly_ratio)/fa_min,
            "1502_small_vev": {"physical_fa_GeV": fit_fa, "mass_keV_central":light_mass*1e6,
                "Fsd_GeV": fv, "K_to_pi_a_partial_width_GeV": width,
                "observed_total_width_GeV":gamma_k, "partial_over_observed_total":width/gamma_k,
                "max_amplitude_retention_even_using_whole_K_width": math.sqrt(gamma_k/width)},
            "enlarged_singlets_factor_5e9_fa_GeV": big,
            "small_angle_theta_edm_cap_at_one_slope_width":theta_cap,
            "single_real_cubic_pin_sharp_torque_test":pin_bounds}


def close(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(close(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-24)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    flags=p.parse_args()
    out=json.loads(json.dumps(solve(), default=float))
    file=BASE/"results.json"
    if flags.write_new:
        with file.open("x",encoding="utf-8") as f:
            f.write(json.dumps(out,indent=2)+"\n")
    elif not close(out,json.loads(file.read_text(encoding="utf-8"))):
        raise ValueError("Round 1503 results changed")
    print("1503 PQ, weak CKM, kaon width, photon slope and cubic torque checks passed.")

