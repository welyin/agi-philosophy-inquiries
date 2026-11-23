"""Reproduce the three-flavor bound without importing the author certificate."""
import argparse
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
M_E, M_MU, M_TAU_LO = 0.00051099895, 0.1056583755, 1.775
ALPHA, HBAR = 0.0072973525205055605, 6.582119569e-25


def edge_bound(m_light, m_heavy, life_s, br, efficiency):
    phase = (1 - (m_light / m_heavy) ** 2) ** 3
    return sqrt(4 * HBAR * br / (life_s * ALPHA * m_heavy**3 * phase
                                 * efficiency * m_light * m_heavy))


def branching_bound(efficiency):
    e_mu = edge_bound(M_E, M_MU, 2.1969803e-6, 1.5e-13, efficiency)
    e_tau = edge_bound(M_E, M_TAU_LO, 288e-15, 3.3e-8, efficiency)
    mu_tau = edge_bound(M_MU, M_TAU_LO, 288e-15, 4.4e-8, efficiency)
    path = sqrt(5) * (e_mu + e_tau) + sqrt(8) * mu_tau
    return {"b_e_mu_GeV-2": e_mu, "b_e_tau_GeV-2": e_tau,
            "b_mu_tau_GeV-2": mu_tau, "path_GeV-2": path,
            "max_negative_ae": 2 * M_E**2 * path}


def matrix_examples():
    m = np.diag([M_E, M_MU])
    scale = 1.e-8
    z_e, z_mu = scale * (.3 + .05j), scale * (1.7 - .05j)
    target = z_e * z_mu
    p = target.real
    eig = np.sort(np.roots([1., -2., p / scale**2])) * scale
    product = M_E * M_MU * (target - p)
    c12 = sqrt(abs(product)) * np.exp(.73j)
    c21 = product / c12
    c = np.array([[M_E*z_e, c12], [c21, M_MU*z_mu]])
    x = c @ np.diag([1/M_E, 1/M_MU])
    u1 = (x - eig[1]*np.eye(2)) @ m / (eig[0] - eig[1])
    u2 = m - u1
    assert all(np.linalg.matrix_rank(u, tol=1.e-12) == 1 for u in (u1, u2))
    assert np.linalg.norm(eig[0]*u1 + eig[1]*u2 - c) < 1.e-23
    assert abs(np.linalg.norm([c12, c21])**2 - 2*abs(product)) < 1.e-31
    assert np.linalg.norm(x.conj().T @ x - x @ x.conj().T) / scale**2 > .01
    # Four independent physical channels for three flavors. Opposite mass
    # contributions need no negative spectral gain.
    masses = np.array([M_E, M_MU, M_TAU_LO])
    us = np.array([[-2*M_E,0,0],[3*M_E,0,0],[0,M_MU,0],[0,0,M_TAU_LO]])
    gains = np.array([3e-6, 1e-6, 5e-8, 1e-8])
    assert np.allclose(us.sum(axis=0), masses, rtol=1.e-15, atol=1.e-18)
    delta = 2 * masses * (gains @ us)
    return {"complex_example_gains_GeV-2": eig.tolist(),
            "minimum_offdiag_norm_sq_GeV-2": 2*abs(product),
            "four_channel_anomalies": delta.tolist()}


def report():
    src = json.loads((ROOT.parent/"1532"/"independent_results.json").read_text(encoding="utf-8"))
    edm = src["ten_percent_alpha_budget"]["d_in_1e_minus30"]
    branch = branching_bound(.8)
    one = branching_bound(1.)
    cs_width = (166 - 2*15) - (46 + 2*27)
    cs = cs_width * 1.e-9 * 8.3e-6
    stress = branch["max_negative_ae"] / .8 + 1e-13
    # Use the full range of the QED coefficients and a separate slope budget.
    q0, q1 = 137.035998, 137.036
    drift = sum(n*c / (pi**n*q0**(n+1))
                for n,c in [(2,.35),(3,2),(4,3),(5,12)]) + 1e-8
    slope = .5 / (pi*q1**2) - drift
    assert slope > 8.3e-6 and stress < cs
    return {"round":1591, "three_flavor_eta_1":one, "three_flavor_eta_08":branch,
            "conservative_Cs_two_width_min_abs":cs, "historic_derivative_lower":slope,
            "with_20pct_amplitude_cancellation_and_1e-13_diagonal":stress,
            "Cs_min_separation_with_stress":cs-stress,
            "critical_common_efficiency_without_extras":(one["max_negative_ae"]/cs)**2,
            "critical_common_efficiency_with_extras":
                (one["max_negative_ae"]/(.8*(cs-1e-13)))**2,
            "two_flavor_trace_muon_WP25_e_lower":-1.64e-9*(M_E/M_MU)**2,
            "two_flavor_trace_muon_edm_e_cm":[-edm[1]*1e-30*M_MU/M_E,
                                             -edm[0]*1e-30*M_MU/M_E],
            "examples":matrix_examples()}


def check_tree(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for k in actual:
            check_tree(actual[k], saved[k])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for a,b in zip(actual,saved):
            check_tree(a,b)
    elif isinstance(actual,(int,float)):
        assert np.isclose(actual,saved,rtol=2.e-11,atol=1.e-45)
    else:
        assert actual == saved


if __name__ == "__main__":
    cmd = argparse.ArgumentParser()
    cmd.add_argument("--write",action="store_true")
    cmd.add_argument("--check",action="store_true")
    args=cmd.parse_args()
    result=report()
    path=ROOT/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8",newline="\n") as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write("\n")
    else:
        check_tree(result,json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(result,ensure_ascii=False))
