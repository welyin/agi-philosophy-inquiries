"""Independent replay of the two-layer test, using LEP fits without flavour universality."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
CENTRE = np.array([91.1876, 41.4807, 20.8038, 20.7842, 20.7644])
ERROR = np.array([.0021, .0325, .0497, .0335, .0448])
CORRELATION = np.array([
    [1, -.0507, .0783, 0, .0019],
    [-.0507, 1, .1138, .1391, .0987],
    [.0783, .1138, 1, .0694, .0464],
    [0, .1391, .0694, 1, .0696],
    [.0019, .0987, .0464, .0696, 1],
])
UNIT_NB = 389379.37
REFERENCE = 1.99060


def line_shape(row):
    mz, had_nb, re, rmu, rtau = row
    inv_e = math.sqrt(12 * math.pi * re * UNIT_NB / (mz**2 * had_nb))
    missed = inv_e - re - 1 - re / rmu - re / rtau
    return np.array([missed / REFERENCE, re / rmu, re / rtau, inv_e])


def charged_bound(measured):
    # a=2 sin²(theta_eff) in [0.44,0.48]; rate=eta*tree with unknown eta.
    dlo, dhi = measured[0] / 1.005, measured[1] / .995
    qmax, distance = .44**2 + .56**2, .52
    def shift(d):
        return (d - 1) * qmax / (math.sqrt(distance**2 + (d - 1)*qmax) + distance)
    assert dlo < 1 < dhi
    return [1 + shift(dlo), 1 + shift(dhi)]


def algebra_rank(s):
    # Singular powers rank from 1,t,t² in the three blocks.
    t = np.square(s)
    return int(np.linalg.matrix_rank(np.column_stack([t**k for k in range(3)]),
                                     tol=1e-12)) * 3


def run():
    observed = line_shape(CENTRE)
    jac = np.column_stack([(line_shape(CENTRE + 1e-4 * h) -
                             line_shape(CENTRE - 1e-4 * h)) / 2e-4
                            for h in np.diag(ERROR)])
    sigma = np.sqrt(np.diag(jac @ CORRELATION @ jac.T))
    corners = np.array([line_shape(CENTRE + 2 * ERROR * sg)
                        for sg in np.array(list(itertools.product([-1, 1], repeat=5)))])
    lower, upper = corners.min(axis=0), corners.max(axis=0)
    lep_mu, lep_tau = [charged_bound([lower[i], upper[i]]) for i in (1, 2)]
    z_coefficient = .98 * 1.1663787e-5 * 91.18**3 / (6 * math.sqrt(2) * math.pi * 2.51)
    limits = np.array([2.62e-7, 4.5e-6, 5.6e-6])
    ab, ac, bc = np.sqrt(limits / z_coefficient)
    off = np.array([[0, ab, ac], [ab, 0, bc], [ac, bc, 0]])
    eigen_box = [float(np.linalg.eigvalsh(np.diag([1, lep_mu[0], lep_tau[0]])-off)[0]),
                 float(np.linalg.eigvalsh(np.diag([1, lep_mu[1], lep_tau[1]])+off)[-1])]
    # W-only replay, using different matrix bound calculation from algebraic moment test.
    weak = [(np.array([1, (.9995-.009)/1.002, (.975-.064)/1.002]), -1),
            (np.array([1, (.9995+.009)/.998, (.975+.064)/.998]), 1)]
    w_box = [float(np.linalg.eigvalsh(np.diag(d)+sg*off)[0 if sg < 0 else -1])
             for d, sg in weak]
    # LFV upper cap belongs to otherwise unclassified width; fixed standard electron response.
    tr_min = (lower[0] * REFERENCE - upper[3] * sum(limits)) / (REFERENCE + .00042)
    tr_max = upper[0] * REFERENCE / (REFERENCE - .00042)
    # Exact neutrino-only contraction implies every s^4 >= tr G² - 2.
    pure_lower = max(tr_min - 2, 0)**.25
    rows = {}
    for label, c, eps in [("general", 1., .001),
                          ("contraction", ((2.9963 - 1)/2)**.25, .0002)]:
        base = np.array([1., c, c])
        split = np.array([1., c+eps, (2*c**4-(c+eps)**4)**.25])
        assert abs(sum(base**4)-sum(split**4)) < 2e-14
        if label == "contraction":
            assert max(split) <= 1
        assert algebra_rank(split) == 9
        assert algebra_rank(base) == (3 if c == 1 else 6)
        # Initial charged e/mu superposition, final sum of antisymmetric e/mu POVM
        # over the two charge sectors. For c=1 the baseline has zero dark rate.
        t = math.pi / eps
        def probability(s):
            charged = (np.cos(s[0]*t) - np.cos(s[1]*t))/2
            neutral = (np.sin(s[0]*t) - np.sin(s[1]*t))/2
            return float(charged**2 + neutral**2)
        if c == 1:
            assert probability(split) > 1-1e-12 and probability(base) < 1e-12
        rows[label] = dict(equal_base=base.tolist(), split=split.tolist(),
                           equal_trace=float(sum(split**4)), ranks=[algebra_rank(base),9],
                           rate_gram_diag=(split**2).tolist(),
                           norm_difference=float(max(abs(split-base))),
                           integrated_area=t, dark_prob=[probability(base),probability(split)])
    return dict(lep_B13_central=observed.tolist(), correlated_linearized_sigma=sigma.tolist(),
                two_width_coordinate_box=[lower.tolist(),upper.tolist()],
                trace_G_squared_box=[tr_min,tr_max],
                contraction_only_lower_singular=pure_lower,
                charged_ratio_to_diagonal={"mu":lep_mu,"tau":lep_tau},
                LFV_plus_charged_G_spectrum=eigen_box,
                singular_spectrum=np.sqrt(eigen_box).tolist(),
                LFV_plus_W_G_spectrum=w_box,
                witness=rows)


def normalize(value):
    if isinstance(value, dict):
        return {k:normalize(v) for k,v in value.items()}
    if isinstance(value, list):
        return [normalize(v) for v in value]
    return round(float(value), 9)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--write-new", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = cli.parse_args()
    calculated = normalize(run())
    dest = ROOT / "independent_results.json"
    if args.write_new:
        with dest.open("x", encoding="utf-8") as f:
            f.write(json.dumps(calculated, indent=2, ensure_ascii=False) + "\n")
    elif args.check:
        assert json.loads(dest.read_text(encoding="utf-8")) == calculated
        print("1473 independent check passed")
    else:
        print(json.dumps(calculated, indent=2))
