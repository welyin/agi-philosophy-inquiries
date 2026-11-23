"""Solve the weighted trace directly, and transport all nine Yukawa singular values."""
import argparse
import json
from itertools import product
from math import exp, log, pi, sqrt
from pathlib import Path

import numpy as np


LOOP = 16*pi**2
GAUGE = np.array([.461228, .65096, 1.2123])
GAUGE_ERR = np.array([.000026, .00004, .0046])
B = np.array([-41/10, 19/6, 7]) * 2/LOOP   # inverse squares
TRACES = np.array([[11/10, 9/10], [3/2, 1/2], [2, 0]])
MZ = 91.188
Y = np.array([[7.04e-6, .00356, .967], [1.54e-5, .000306, .01630],
              [2.77713e-6, .000585042, .0099378]])
Y_ERR = np.array([[.15e-6, .00006, .004], [.02e-5, .000004, .00009],
                  [.00036e-6, .000000075, .0000014]])


def solve_trace(g):
    # Three unknowns: q weight*normalization, lepton weight*normalization and log scale.
    q, ell, t = np.linalg.solve(np.column_stack((TRACES, -B)), 1/g**2)
    xs = 1/g**2 + t*B
    assert q>0 and ell>0 and np.allclose(TRACES@[q,ell], xs)
    return t, 3*q/ell, sqrt(xs[1]/xs[2])


def minimum_flavour_weight(d,e):
    return (sqrt(e)-1/sqrt(e))/(sqrt(d)-1/sqrt(d))


def envelope(tmax):
    upper = Y + Y_ERR
    lower = Y - Y_ERR
    u0 = np.sum(upper[:2]**2)
    l0 = np.sum(upper[2]**2)
    g3_min_at_mz = GAUGE[2]-GAUGE_ERR[2]
    beta = B[2]*g3_min_at_mz**2
    e_integral = 7/beta * (1-(1+beta*tmax)**(-1/7))
    denominator = 1-9*u0*e_integral/LOOP
    lowest_g2_squared = 1/(1/(GAUGE[1]-GAUGE_ERR[1])**2 + B[1]*tmax)
    lmax = l0*exp((5*.001+12)*tmax/LOOP)
    assert 2*.001 < 4.5*lowest_g2_squared
    assert denominator>0 and u0/denominator < 2 and lmax < .001
    delta_log_r = -log(denominator)/6
    r_limit = upper[1,1]/lower[1,0]*exp(delta_log_r)
    e_limit = lower[2,1]/upper[2,0]
    return {"comparison_denominator": denominator, "bootstrap_U_ceiling": u0/denominator,
            "bootstrap_L_ceiling": lmax, "integrated_log_down_error": delta_log_r,
            "down_ratio_upper":r_limit, "lepton_ratio_lower":e_limit,
            "required_weight_lower":minimum_flavour_weight(r_limit,e_limit)}


def rotation(i,j,sine,phase=0):
    matrix = np.eye(3,dtype=complex)
    matrix[i,i] = matrix[j,j] = sqrt(1-sine**2)
    matrix[i,j] = sine*np.exp(-1j*phase)
    matrix[j,i] = -matrix[i,j].conjugate()
    return matrix


def running_rates(t, yukawas):
    g1,g2,g3 = 1/(1/GAUGE**2+t*B)
    up, down, leptons = yukawas @ yukawas.conj().transpose(0,2,1)
    common_trace = np.trace(3*(up+down)+leptons).real
    costs = np.array([17*g1/20+9*g2/4+8*g3, g1/4+9*g2/4+8*g3, 9*(g1+g2)/4])
    directions = 1.5*np.array([up-down,down-up,leptons])
    scalar = (common_trace-costs)[:,None,None]*np.eye(3)
    return ((directions+scalar) @ yukawas)/LOOP


def transport(tstar, nstep):
    ckm = rotation(1,2,.04193) @ rotation(0,2,.00370,1.139) @ rotation(0,1,.2251)
    matrices = np.array([np.diag(Y[0]), ckm @ np.diag(Y[1]), np.diag(Y[2])])
    h = tstar/nstep
    for i in range(nstep):
        t = i*h
        k1 = running_rates(t,matrices)
        k2 = running_rates(t+h/2,matrices+h*k1/2)
        k3 = running_rates(t+h/2,matrices+h*k2/2)
        k4 = running_rates(t+h,matrices+h*k3)
        matrices += h*(k1+2*k2+2*k3+k4)/6
    u, d, e = np.sort(np.linalg.svd(matrices, compute_uv=False),axis=1)
    return np.array([d[1]/d[0],e[1]/e[0],u[1]/u[0]])


def two_tabulated_points():
    # These are tabulated RGE outputs, not measurements at two new accelerators.
    cases = [(1e9,[.509051,.57639,.6872],[.000035,.00002,.0008]),
             (1e12,[.534498,.55157,.6052],[.000040,.00002,.0006])]
    checks = []
    for scale,c,err in cases:
        raw = np.array(c) + np.array(list(product([-1,1],repeat=3)))*err
        xx = 1/raw**2
        residual = xx @ np.array([5/3,-3,4/3])
        checks.append({"scale_GeV":scale, "trace_residual_box":[min(residual),max(residual)],
                       "min_colour_to_weak": min(raw[:,2]/raw[:,1])})
    assert checks[0]["trace_residual_box"][0]>0>checks[1]["trace_residual_box"][1]
    return checks


def compute():
    t0,w0,ratio = solve_trace(GAUGE)
    corners = np.array([solve_trace(GAUGE + np.array(s)*GAUGE_ERR)
                        for s in product([-1,1],repeat=3)])
    lo, hi = corners.min(axis=0), corners.max(axis=0)
    env = envelope(hi[0])
    fine = transport(t0,256)
    coarse = transport(t0,128)
    assert max(abs(1-fine/coarse)) < 2e-10
    assert fine[0] <= env["down_ratio_upper"] and fine[1] >= env["lepton_ratio_lower"]
    assert lo[2] > 1 and hi[1] < 3 < env["required_weight_lower"]
    return {"inputs":"https://arxiv.org/html/2510.01312v2",
            "one_loop":"https://arxiv.org/pdf/hep-ph/0501272#page=40",
            "central_log_scale": t0, "central_GeV":MZ*exp(t0),
            "central_gauge_weight":w0,"central_colour_to_weak":ratio,
            "log_scale_box":[lo[0],hi[0]],"gauge_weight_box":[lo[1],hi[1]],
            "colour_to_weak_box":[lo[2],hi[2]], "flavour_envelope":env,
            "central_same_scale_down_lepton_up_ratios":fine.tolist(),
            "central_required_flavour_weight": minimum_flavour_weight(fine[0],fine[1]),
            "split_factor_lower":env["required_weight_lower"]/hi[1],
            "tabulated_RGE_crosscheck":two_tabulated_points()}


def rounded(value):
    if isinstance(value, dict):
        return {k: rounded(v) for k,v in value.items()}
    if isinstance(value, list):
        return [rounded(v) for v in value]
    return float(format(value,".9g")) if isinstance(value,float) else value


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true",help="verify saved results (default)")
    parser.add_argument("--write",action="store_true")
    options=parser.parse_args()
    file=Path(__file__).with_name("independent_results.json")
    result=rounded(compute())
    if options.write:
        with file.open("x",encoding="utf-8") as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write("\n")
    else:
        assert result==json.loads(file.read_text(encoding="utf-8"))
    print("1369 independent trace, nine Yukawas and full-box exclusion verified.")
