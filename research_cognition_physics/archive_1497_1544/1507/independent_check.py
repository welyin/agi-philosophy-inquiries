"""Independent photon continuum and trace checks for round 1507."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def continuum(energy, counts, bg, lumi, mu_nb, eff, isr, r, stat, syst):
    raw = (counts-bg)/(1000*lumi*mu_nb*eff*isr)
    assert abs(r-raw) < .002
    low, high = r-2*(stat+syst), r+2*(stat+syst)
    k_min, k_max = .95, 1.25
    def charges(d, anomaly):
        if anomaly:
            return (1+1/d)/2, (-1+1/d)/2
        return 2/3, -1/3
    def norm(d, anomaly):
        u, q = charges(d, anomaly)
        return d*(u*u+2*q*q)
    fits = {}
    for anomaly in (False, True):
        label = "anomaly_linked" if anomaly else "fixed"
        compatible = [d for d in range(1, 31)
                      if norm(d, anomaly)*k_min <= high and norm(d, anomaly)*k_max >= low]
        assert compatible == [3]
        assert norm(4, anomaly)*k_min > high and norm(2, anomaly)*k_max < low
        fits[label] = {"integers_checked": compatible, "required_K": {
            str(d): [low/norm(d,anomaly), high/norm(d,anomaly)] for d in (2,3,4)}}
        for d in (1, 2, 3, 4, 7, 11):
            u,q = charges(d, True)
            assert abs(d*(u*u-q*q)-1) < 2e-15
    return {"GeV":energy, "published":r, "from_tabulated_counts":raw,
            "twice_linear_width": [low,high], "fits":fits}


def trace_ratio(matrices):
    c = sum(a @ a for a in matrices)
    traces = np.array([[np.trace(a @ b).real for a in matrices] for b in matrices])
    return float(np.sum(traces*traces)/np.trace(c @ c).real)


def basis():
    es = np.eye(3, dtype=complex)
    hs = []
    for i,j in ((0,1),(0,2),(1,2)):
        m = np.outer(es[i], es[j])
        hs.extend([(m+m.T)/math.sqrt(2), 1j*(m-m.T)/math.sqrt(2)])
    return hs + [np.diag([1,-1,0])/math.sqrt(2),np.diag([1,1,-2])/math.sqrt(6)]


def audit():
    photons = [continuum(2.4,96627,2331,3.415,15.079,.6729,1.204,2.260,.008,.042),
               continuum(2.8,83802,2075,3.753,11.078,.7225,1.219,2.233,.008,.055)]
    hs = basis()
    threshold = (3-math.sqrt(3))/2
    optimum = 2*math.sqrt(3)-3
    assert abs(trace_ratio(hs)-3/8) < 5e-15
    values = {}
    for t in (0, optimum, 1, 8):
        ratio = trace_ratio([hs[0],hs[1],hs[6],math.sqrt(t)*hs[7]])
        assert abs(ratio-(6+2*t*t)/(9+2*t+t*t)) < 2e-14
        assert ratio >= threshold-2e-15
        values[str(t)] = ratio
    diagonal_bound = 1.5
    for angle in (0,.08,.49,1.3):
        a = math.cos(angle)*hs[6]+math.sin(angle)*hs[7]
        b = -math.sin(angle)*hs[6]+math.cos(angle)*hs[7]
        for t in (0,.13,1,8):
            ratio = trace_ratio([a,math.sqrt(t)*b])
            formula = (1+t*t)/((1+t*t)/2+t/3)
            assert abs(ratio-formula) < 3e-14 and ratio >= diagonal_bound-2e-15
    # Add a ninth (traceful) generator. Invariant positive kinetic terms
    # for the simple and Abelian factors allow independent couplings.
    y = .37
    roots = np.roots([3-y,-16*y,24-64*y])
    t_u3 = min(roots.real)
    assert t_u3 > 0 and not np.any(roots.imag)
    r_u3 = trace_ratio(hs + [math.sqrt(t_u3/3)*np.eye(3)])
    assert abs(r_u3-y)<2e-14
    enlarged = hs + [math.sqrt(t_u3/3)*np.eye(3)]
    c = sum(a @ a for a in enlarged)
    x_u3 = sum(np.linalg.norm(a@b-b@a)**2 for a in enlarged for b in enlarged)/np.trace(c @ c).real
    assert abs(x_u3-2.25*(8/(8+t_u3))**2)<3e-14 and 2.02 < x_u3 < 2.32
    return {"BESIII_table_II":"https://doi.org/10.1103/PhysRevLett.128.062004",
            "conditional_inclusive_window":[.95,1.25], "photon_witnesses":photons,
            "traceless_proper_subgroup_min":threshold,
            "max_extra_jet_transfer_at_two_widths":threshold-(.37+2*(.04+.06)),
            "reducible_weak2_plus_abelian_checks":values, "diagonal_min":diagonal_bound,
            "U3_counterexample":{"relative_abelian_squared_coupling":t_u3, "y":r_u3,
                                  "ALEPH_self_coupling_analogue":float(x_u3)}}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
