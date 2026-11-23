"""Neutral-current counterpart of the protected charged-current window."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


MZ = 91.1876
LEP_MEAN, LEP_ERROR = 2.9963, 0.0074
ROOT = Path(__file__).resolve().parent
nodes, weights = np.polynomial.legendre.leggauss(96)
nodes, weights = (nodes + 1) / 2, weights / 2


def integrated(q):
    if q <= 0:
        return 1.0
    if q >= 1:
        return 0.0
    u = math.sqrt(1 - q)
    if q > 0.8:
        # This positive form stays accurate at threshold.
        return float(3 * u**7 * np.dot(weights, nodes**2 *
                     ((1 - nodes**2) / (1 - u*u*nodes**2))**2))
    a = math.log1p(u) - 0.5 * math.log(q)
    return u**3 - 6*q*(a - u) + 1.5*q*(u - q*a)


def single_window(q):
    if q <= 0:
        return 1.0
    if q >= 1:
        return 0.0
    x = (math.sqrt(q*q + 48*q) - q) / 6
    return (1 - x)**1.5 * (1 - q/x)**2


def min_budget(target, charge_factor):
    lo, hi = 0.0, MZ / (2 * math.sqrt(charge_factor))
    for _ in range(65):
        mid = (lo + hi) / 2
        if integrated(4 * charge_factor * (mid / MZ)**2) > target:
            lo = mid
        else:
            hi = mid
    return hi


def takagi(matrix):
    if np.isrealobj(matrix):
        ev, u = np.linalg.eigh(matrix)
        idx = np.argsort(abs(ev))
        return abs(ev[idx]), u[:, idx] * np.where(ev[idx] < 0, 1j, 1)
    _, _, vh = np.linalg.svd(matrix)
    u = vh.conj().T[:, ::-1]
    u *= np.exp(-0.5j * np.angle(np.diag(u.T @ matrix @ u)))
    diagonal = u.T @ matrix @ u
    masses = np.diag(diagonal).real
    assert np.max(abs(diagonal - np.diag(masses))) < 1e-7
    return np.maximum(masses, 0), u


def width(masses, current):
    x = (masses / MZ)**2
    a, b = x[:, None], x[None, :]
    open_pair = masses[:, None] + masses[None, :] < MZ
    root = np.sqrt(np.maximum((1-a-b)**2 - 4*a*b, 0.0))
    pieces = root * ((1 - (a+b)/2 - (a-b)**2/2) * abs(current)**2 -
                     3 * np.sqrt(a*b) * (current**2).real)
    assert np.min(pieces[open_pair], initial=0) > -1e-12
    return float(np.sum(pieces[open_pair]))


def example(seed, p, mirrors, singlets):
    rng = np.random.default_rng(seed)
    n = p + 2*mirrors + singlets
    plus, minus = p+mirrors, p+2*mirrors
    def random_matrix(rows, cols):
        return rng.normal(size=(rows, cols)) + 1j*rng.normal(size=(rows, cols))
    # A weak basis aligned to ker M0 on the positive-charge doublets.
    m = np.zeros((n, n), complex)
    m[p:plus, plus:minus] = 90 * random_matrix(mirrors, mirrors)
    m[:minus, minus:] = 3 * random_matrix(minus, singlets)
    m[minus:, minus:] = 60 * random_matrix(singlets, singlets)
    m += m.T.copy()
    budgets = np.sort(np.linalg.svd(m[:, :p], compute_uv=False))
    charges = np.repeat([1.0, -1.0, 0.0], [plus, mirrors, singlets])
    masses, u = takagi(m)
    current = u.conj().T @ (charges[:, None] * u)
    c = 2 if mirrors else 1
    for mu in np.linspace(2, MZ/2, 24):
        keep = masses <= mu
        protected = np.maximum(1 - c*(budgets/mu)**2, 0)
        assert np.sum(abs(current[np.ix_(keep, keep)])**2) + 1e-9 >= protected @ protected
    exact = width(masses, current)
    bound = sum(integrated(4*c*(b/MZ)**2) for b in budgets)
    assert exact + 1e-9 >= bound > 0
    return exact, bound


def run():
    main = json.loads((ROOT / "results.json").read_bytes())
    cfg = json.loads((ROOT / "inputs.json").read_bytes())
    v = json.loads((ROOT / cfg["old_width_and_vev"]).read_bytes())["v_GeV"]
    cols = json.loads((ROOT / cfg["old_loop_norm"]).read_bytes())["Nc"]
    for y in cfg["weak_budget_examples"] + [4*math.pi/math.sqrt(cols)]:
        mass_scale = v*y/math.sqrt(2)
        charged = main["optimized_extra_pair_rate_over_one_lepton"][str(round(y,8))]
        for energy, outcome in charged.items():
            moment = 8*(mass_scale/float(energy))**2
            predicted = cols*(cfg["family_index"]-cfg["known_families"])*single_window(moment)
            assert abs(predicted - outcome["relative_rate"]) < 1e-8
    for q in np.linspace(0.01, 0.99, 43):
        u = math.sqrt(1 - q)
        quad = 3*u**7 * np.dot(weights, nodes**2 *
                               ((1-nodes**2)/(1-u*u*nodes**2))**2)
        assert abs(integrated(q) - quad) < 1e-11
        assert integrated(q) >= single_window(q) - 1e-14
    for mass in [0.0, 1.0, 22.0, 44.0]:
        x = (mass/MZ)**2
        m, u = takagi(np.array([[0.0, mass], [mass, 0.0]]))
        assert abs(width(m, u.conj().T @ np.diag([1, 0]) @ u) -
                   math.sqrt(1-4*x)*(1-x)) < 1e-12
        assert abs(width(np.array([mass]), np.ones((1,1))) - (1-4*x)**1.5) < 1e-12
    upper = LEP_MEAN + 2*LEP_ERROR
    table = []
    for factor in [1.0, 0.98]:
        for c in [1, 2]:
            table.append(dict(width_matching_lower_factor=factor, charge_factor=c,
                              common_budget_min_GeV=min_budget(upper/(4*factor), c),
                              fourth_min_if_three_null_GeV=min_budget(upper/factor-3, c)))
    branches = []
    for bare in [0.0, 6000.0]:
        m, u = takagi(np.array([[0., 60.], [60., bare]]))
        branches.append(dict(singlet_bare_GeV=bare, extra_masses_GeV=m.tolist(),
                             light_total_with_three_massless=3+width(m, u.conj().T @ np.diag([1,0]) @ u)))
    sample = [example(100+k, 4, k%3, 5) for k in range(75)]
    return dict(tree_level=True, pole_mass_GeV=MZ, LEP_reanalysis=[LEP_MEAN, LEP_ERROR],
                illustrative_upper=upper, bounds=table, competitors=branches,
                random_complex_models=len(sample), minimum_width_margin=min(x-y for x,y in sample))


def clean(x):
    if isinstance(x, dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x, list):
        return [clean(v) for v in x]
    return float(f"{x:.11g}") if isinstance(x, float) else x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = clean(run())
    payload = json.dumps(result, indent=2) + "\n"
    path = ROOT / "independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes()) == result
    else:
        path.write_text(payload, encoding="utf-8")
    print("1776 neutral spectral bound", hashlib.sha256(payload.encode()).hexdigest()[:12])
