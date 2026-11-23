"""Test the frozen messenger orientation with correlated, inclusive weak currents."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE/name).read_bytes())


def upper_crossing(func, target, end):
    lo, hi = 0., end
    assert func(lo) <= target < func(hi)
    for _ in range(64):
        mid = (lo+hi)/2
        if func(mid) <= target:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def stable(value):
    if isinstance(value, dict):
        return {k: stable(v) for k, v in value.items()}
    if isinstance(value, (list, np.ndarray)):
        return [stable(v) for v in value]
    return float(format(value, ".12g"))


def compute():
    cfg = load("inputs.json")
    old = runpy.run_path(str(HERE/cfg["charged_model"]))
    e, _ = old["OLD"]["charged_basis"](old["read"](old["INPUT"]["charged"]))
    q = math.sqrt(7)/5
    np.testing.assert_allclose(math.sqrt(q), load(cfg["static_bridge"])["Hq_over_H1"], atol=1.e-11)
    w = abs(e[2, :])**2
    assert 0 < w[0] < w[1] < w[2] < 1
    lim = 1-q

    def coords(f):
        u = q*f/(1-q-f)
        return {"contrast":f, "tL_over_H1": math.sqrt(u),
                "tL_over_g":math.sqrt((1+q)*u), "S":(1-f*w)/(1+u)}

    def ratios(f):
        s = 1-f*w
        return np.sqrt([s[2]/s[1], s[1]/s[0]])

    hf = cfg["hflav"]
    target = np.array(hf["ratios"])
    sig = np.array(hf["uncertainties"])
    rho = hf["correlation"]
    cov = np.outer(sig,sig) * np.array([[1.,rho],[rho,1.]])

    def chi(f):
        diff = ratios(f)-target
        return float(diff @ np.linalg.solve(cov, diff))

    baseline = chi(0.)
    cut = cfg["diagnostic"]["hflav_delta_chi2"]
    joint = upper_crossing(chi, baseline+cut, lim)
    pi = cfg["pienu"]
    pion_max = pi["ratio"]+cfg["diagnostic"]["pienu_marginal_widths"]*pi["uncertainty"]
    pion = upper_crossing(lambda f: 1/ratios(f)[1], pion_max, lim)
    assert joint < pion
    def report(f):
        out = coords(f)
        ratios_at_f = ratios(f)
        out.update(HFLAV_tau_mu_and_mu_e=ratios_at_f, PIENU_e_mu=1/ratios_at_f[1],
                   HFLAV_chi2=chi(f))
        # Exact row-sum rule, also tested away from zero.
        st, sm, se = out["S"][::-1]
        np.testing.assert_allclose((st-se)*(w[1]-w[0]), (sm-se)*(w[2]-w[0]), atol=3.e-15)
        return out

    points = []
    for k in cfg["diagnostic"]["illustrative_tL_over_g"]:
        u = k*k/(1+q)
        f = (1-q)*u/(q+u)
        row = report(f)
        # Direct null projectors for three massive spectator doublets in this frame.
        D = np.array([1+q,1+q,q*(1+q)])
        direct = e.conj().T @ np.diag(D/(D+k*k)) @ e
        np.testing.assert_allclose(direct.diagonal(), row["S"], atol=1.e-14)
        points.append(row)
    # Preserve the tree-level source/sink products on all legs but align left
    # couplings with the square roots of the unequal heavy denominators.
    D = np.array([1+q,1+q,q*(1+q)])
    left = .05*np.sqrt(D/D[0])
    right = .02*np.sqrt(D[0]/D)
    np.testing.assert_allclose(left*right, np.full(3,.001))
    balanced = D/(D+left**2)
    np.testing.assert_allclose(balanced, [balanced[0]]*3)
    return stable(dict(charged_weights=w, universal_difference_ratio=(w[2]-w[0])/(w[1]-w[0]),
                       HFLAV_covariance=cov, HFLAV_min_chi2=baseline, HFLAV_upper=report(joint),
                       PIENU_upper=report(pion), joint_rectangle_and_correlated_upper=report(min(pion,joint)),
                       sample=points, competitor_left=left, competitor_right=right,
                       competitor_equal_residues=balanced))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args=parser.parse_args()
    result = compute()
    if args.check:
        assert result == load("results.json")
        print("PASS: frozen charged basis, exact row identity, correlated HFLAV, PIENU and competitor")
    else:
        with (HERE/"results.json").open("x", encoding="utf8") as file:
            file.write(json.dumps(result, indent=2)+"\n")
