"""One-loop Higgs double-zero test with the exact two-singlet control spectrum of 1607.

Low renormalized SM inputs fixed at 200 GeV; step matching, no Yukawas
for heavy mass eigenstates. This is not a pole-mass prediction for the extension.
"""
import argparse
import json
import math
from pathlib import Path

from check import BETA, LOOP, critical_y_squared

HERE = Path(__file__).resolve().parent
DIRAC_DOWN = (4/9, 0, 2/3)


def squares(t, couplings, gates):
    extra_log = sum(max(0, t-threshold) for threshold in gates)
    return tuple(1/(1/g**2-2*(b*t+db*extra_log)/LOOP)
                 for g,b,db in zip(couplings, BETA, DIRAC_DOWN))


def slopes(t, q, lam, couplings, gates):
    gy, gw, gs = squares(t, couplings, gates)
    beta_q = q*(9*q-17*gy/6-9*gw/2-16*gs)
    beta_lam = (24*lam**2+(12*q-9*gw-3*gy)*lam+
                6*(critical_y_squared(gy,gw)**2-q**2))
    return beta_q/LOOP, beta_lam/LOOP


def integrate(t, couplings, gates, total_steps):
    x,y,_ = squares(t, couplings, gates)
    q, lam = critical_y_squared(x,y), 0.
    high = t
    # Split exactly at each pole. No integrator step crosses a threshold.
    for low in sorted([a for a in gates if 0<a<t]+[0.], reverse=True):
        count = max(8, math.ceil(total_steps*(high-low)/t))
        dt = (low-high)/count
        for i in range(count):
            now = high+i*dt
            a,b = slopes(now,q,lam,couplings,gates)
            c,d = slopes(now+dt/2,q+dt*a/2,lam+dt*b/2,couplings,gates)
            e,f = slopes(now+dt/2,q+dt*c/2,lam+dt*d/2,couplings,gates)
            g,h = slopes(now+dt,q+dt*e,lam+dt*f,couplings,gates)
            q += dt*(a+2*c+2*e+g)/6
            lam += dt*(b+2*d+2*f+h)/6
        high = low
    return math.sqrt(q),lam


def root(cfg, masses, steps):
    ref = cfg["reference"]
    unit = cfg["matching_scale"]
    g = tuple(ref[key] for key in ("gY","g2","g3"))
    gates = [math.log(m/unit) for m in masses]
    lo, hi = [math.log(m/unit) for m in cfg["window"]]
    assert integrate(lo,g,gates,steps)[1] < ref["lambda"] < integrate(hi,g,gates,steps)[1]
    for _ in range(40):
        t = (lo+hi)/2
        if integrate(t,g,gates,steps)[1] > ref["lambda"]:
            hi = t
        else:
            lo = t
    t = (lo+hi)/2
    yt, lam = integrate(t,g,gates,steps)
    x, y, _ = squares(t,g,gates)
    q = critical_y_squared(x,y)
    dq, beta = slopes(t,q,0.,g,gates)
    assert abs(beta)<1e-15 and abs(lam-ref["lambda"])<1e-10
    # Along the upward flow, local positivity is a necessary minimum check.
    delta = 1e-5
    beta_above = slopes(t+delta,q+delta*dq,0.,g,gates)[1]
    beta_below = slopes(t-delta,q-delta*dq,0.,g,gates)[1]
    curvature_proxy = (beta_above-beta_below)/(2*delta)
    assert curvature_proxy>0
    return {"yt_at_200":yt,"double_zero_field_GeV":unit*math.exp(t),
            "d_beta_lambda_dt_at_zero":curvature_proxy}


def evaluate():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf-8-sig"))
    old = json.loads((HERE.parent/"1607"/"independent_results.json").read_text(encoding="utf-8"))
    masses = old["two_weak_singlet_control"]["gauge_basis_five_masses_GeV"][-2:]
    pure = root(cfg,[],512)
    heavy = root(cfg,masses,512)
    coarse = root(cfg,masses,256)
    assert abs(heavy["yt_at_200"]-coarse["yt_at_200"])<2e-9
    assert abs(heavy["double_zero_field_GeV"]/coarse["double_zero_field_GeV"]-1)<3e-8
    assert .8993<pure["yt_at_200"]<.8995 and .9104<heavy["yt_at_200"]<.9107
    logsum = sum(math.log(3e18/m) for m in masses)
    return {"heavy_masses_GeV":masses, "pure_SM_proxy":pure, "two_singlet_proxy":heavy,
            "added_inverse_gauge_squares_at_3e18":[-2*c*logsum/LOOP for c in DIRAC_DOWN],
            "weak_yukawa_shift":heavy["yt_at_200"]-pure["yt_at_200"],
            "formal_field_ratio":heavy["double_zero_field_GeV"]/pure["double_zero_field_GeV"]}


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args = cli.parse_args()
    data = evaluate()
    output = HERE/"joint_threshold_results.json"
    if args.check:
        assert data == json.loads(output.read_text(encoding="utf-8"))
    else:
        with output.open("x",encoding="utf-8") as file:
            file.write(json.dumps(data,indent=2,allow_nan=False)+"\n")
    print("1608 cross-sector check: electroweak-singlet heavy flavor thresholds included")
