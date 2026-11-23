"""Exact Weyl-field budgets and a common-column quark counterexample."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import json
import math
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
Q = (3, 2, F(1, 6))
L = (1, 2, F(-1, 2))
U = (3, 1, F(2, 3))
D = (3, 1, F(-1, 3))
E = (1, 1, F(-1))
SM = Counter({(Q, 1): 3, (L, 1): 3, (U, -1): 3, (D, -1): 3, (E, -1): 3})
BASE = tuple(map(F, CFG["SM"]))


def field_table(fields):
    beta = list(map(F, CFG["bosons_plus_one_Higgs"]))
    anomalies = [F(0) for _ in range(5)]
    balance = Counter()
    charge_weight = F(0)
    doublets = 0
    for ((d, w, y), hand), n in fields.items():
        assert n >= 0 and (6*y - (1 if w == 2 else 0)*3 + 2*(d == 3)) % 6 == 0
        beta[0] += F(n*d*(w == 2), 3)
        beta[1] += F(n*w*(d == 3), 3)
        beta[2] += F(2*n*d*w, 3)*y**2
        anomalies[0] += hand*n*w*(d == 3)
        anomalies[1] += hand*n*w*(d == 3)*y/2
        anomalies[2] += hand*n*d*(w == 2)*y/2
        anomalies[3] += hand*n*w*d*y
        anomalies[4] += hand*n*w*d*y**3
        doublets += n*d*(w == 2)
        for t in ((F(-1, 2), F(1, 2)) if w == 2 else (0,)):
            charge = y+t
            if charge:
                balance[d, charge] += hand*n
            charge_weight += F(2*n*d, 3)*charge**2
    assert not any(balance.values()) and not any(anomalies) and doublets % 2 == 0
    assert beta[0]+beta[2] == charge_weight - 7
    return tuple(value-base for value, base in zip(beta, BASE))


def tables():
    competitors = {}
    for name, field in (("vector_lepton_doublet", L), ("charged_lepton_singlet", E),
                        ("colour_down_singlet", D)):
        table = SM.copy()
        table[field, 1] += 1
        table[field, -1] += 1
        competitors[name] = table
    mirror = SM.copy()
    for field in (U, D, E):
        mirror[field, -1] -= 1
    for field in (Q, L, E):
        mirror[field, -1] += 1
    mirror[E, 1] += 1
    competitors["replace_right_quark_and_lepton_with_doublets_add_charged_singlet"] = mirror
    return competitors


def unique(w2, w3, wY):
    return (w2 < F(4, 3) and (w2 < F(2, 3) or wY < F(2, 3))
            and wY < F(4, 3) and (w3 < F(2, 3) or wY < F(4, 9)))


def mirror_mass():
    model = runpy.run_path(str(HERE / CFG["colour_threshold_model"]))
    yu, yd, V, _, _ = model["threshold"]()
    v = model["VEV"]
    up = yu*v/math.sqrt(2)
    c = abs(V[0, 0])
    down = yd*v/math.sqrt(2)*np.array([c, 1, math.sqrt((1+c)/2)])
    # Retain the finite physical down poles from that example for d and b.
    for k, f in ((0, c), (2, math.sqrt((1+c)/2))):
        down[k] = model["pole"](yd[k]*v/math.sqrt(2), f)["light_GeV"]
    assert up[0] < down[0] and up[-1] > down[-1]
    diff = np.diag(1/up**2) - (V/down**2) @ V.conj().T
    left, right = 0., math.pi/2
    def direction(t):
        return np.array([math.cos(t), 0., math.sin(t)])
    assert diff[0, 0].real > 0 and diff[2, 2].real < 0
    for _ in range(90):
        middle = (left+right)/2
        line = direction(middle)
        if (line @ diff @ line).real > 0:
            left = middle
        else:
            right = middle
    line = direction((left+right)/2)
    column = line/np.linalg.norm(line/up)
    # In each mass frame the image of a unit right vector is precisely column.
    # Complete that vector into a unitary, retaining the last physical column.
    moments = []
    roots = []
    for values, frame in ((up, np.eye(3)), (down, V)):
        state = (frame.conj().T @ column)/values
        assert abs(np.linalg.norm(state)-1) < 2e-12
        _, _, vh = np.linalg.svd(state.conj()[None, :], full_matrices=True)
        basis = np.column_stack((vh[1:].conj().T, state))
        mass = (frame*values) @ basis
        np.testing.assert_allclose(mass[:, -1], column, atol=2e-13)
        np.testing.assert_allclose(mass@mass.conj().T, (frame*values**2)@frame.conj().T,
                                   rtol=2e-12, atol=2e-12)
        spectrum = np.linalg.svd(mass, compute_uv=False)[::-1]
        np.testing.assert_allclose(spectrum, values, rtol=3e-11)
        roots.append(spectrum.tolist())
        moments.append(state)
    wr = np.outer(moments[0], moments[1].conj())
    assert abs(np.linalg.norm(wr)**2 - 1) < 1e-11
    return dict(masses_GeV=roots, common_bare_column_GeV=column.tolist(),
                CKM_abs=np.abs(V).tolist(), right_W_strength_by_flavour=np.abs(wr**2).tolist(),
                total_W_left=float(np.linalg.norm(V)**2), total_W_right=float(np.linalg.norm(wr)**2))


def neutral_sign():
    old = json.loads((HERE / CFG["lepton_Z"]).read_bytes())
    conv = old["conventions"]
    values = old["lep_sld"]["absolute_gL"].values()
    smin = conv["weak_mixing_sin_squared_bounds"][0]
    nmax = conv["common_Z_normalization_bounds"][1]
    residual = conv["lepton_specific_absolute_budget"]
    lower = [2*(smin+(x-CFG["Z_marginal_widths"]*e)/nmax-residual)
             for x, e in values]
    positive = [2*(x-.2312) for x, _ in values]
    assert sum(lower) > 2 and all(0 < x < 1 for x in positive)
    real, sigma = CFG["LEP_SLD_table_7_7_electron_right_absolute"]
    right_upper = real+CFG["Z_marginal_widths"]*sigma
    mirror_lower = conv["common_Z_normalization_bounds"][0]*(
        .5-conv["weak_mixing_sin_squared_bounds"][1])-residual
    branches = [(qL, qR, 3, lR) for qL, qR, lR in product(range(4), repeat=3)
                if 3*(qL+qR)+3+lR == 12 and qL-qR == 3-lR]
    assert branches == [(1, 1, 3, 3), (3, 0, 3, 0)] and mirror_lower > right_upper
    return dict(negative_branch_ordinary_lower=lower, negative_branch_sum=sum(lower),
                positive_branch_exotic_upper_probabilities_without_ordinary=positive,
                branches_after_lepton_left_and_tight_slopes=branches,
                right_electron_all_doublet_lower=mirror_lower, right_electron_data_upper=right_upper)


def run():
    assert field_table(SM) == (0, 0, 0)
    corners = {name: field_table(spectrum) for name, spectrum in tables().items()}
    width_grid = [F(0), F(1, 3), F(4, 9), F(2, 3), F(1), F(4, 3), F(2)]
    count = 0
    for widths in product(width_grid, repeat=3):
        witness = any(all(abs(v) <= w for v, w in zip(row, widths)) for row in corners.values())
        assert witness != unique(*widths)
        count += 1
    endpoint = json.loads((HERE / CFG["ordinary_left_window"]).read_bytes())
    q_strength = endpoint["windows_L"]["ordinary_visible_window"]["visible_visible"]
    assert q_strength > 2
    return dict(SM_slopes=CFG["SM"], corners={k: list(map(str, v)) for k, v in corners.items()},
                closed_rectangles_checked=count, projected_quark_left_strength_1775=q_strength,
                common_column_witness=mirror_mass(), neutral_current=neutral_sign())


def tidy(x):
    if isinstance(x, dict):
        return {k: tidy(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [tidy(v) for v in x]
    if isinstance(x, float):
        return float(f"{x:.10g}")
    return x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    result = tidy(run())
    output = HERE / "results.json"
    if parser.parse_args().check:
        assert json.loads(output.read_bytes()) == result
        print("1781 verified: 343 closed windows, four anomaly-free boundaries, six masses, CKM and signed Z.")
    else:
        output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf8")
