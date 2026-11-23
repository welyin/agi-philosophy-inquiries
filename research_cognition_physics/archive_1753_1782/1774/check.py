"""Coherent spectator-pair preparation, finite mass filters, and a hidden line."""
import argparse
from itertools import combinations
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())


def pairs(n):
    return list(combinations(range(n), 2))


def exterior(matrix):
    return np.array([[np.linalg.det(matrix[np.ix_(r, c)])
                      for c in pairs(matrix.shape[1])] for r in pairs(matrix.shape[0])])


def attempts(prob):
    return math.ceil(math.log(CFG["tail_miss_probability"]) / math.log1p(-prob))


def pair_filter(h, target):
    """Interferometric success factors on the full two-particle spectrum."""
    spec = np.array([h[i] + h[j] for i, j in pairs(len(h))])
    other = [k for k in range(len(spec)) if k != target]
    t = math.pi / np.abs(spec[target] - spec[other])

    # Evaluate phases as differences before multiplication to avoid cancellation.
    def stable_poly(values):
        return np.prod((1-np.exp(-1j*(values[:, None]-spec[other])*t))/2, axis=1)

    exact = stable_poly(spec)
    assert np.max(np.abs(exact - np.eye(len(spec))[target])) < 1e-12
    gap = min(np.diff(np.sort(spec)))
    err = gap * CFG["pair_spectrum_absolute_error_over_minimum_gap"]
    bound = err * t.sum()/2
    perturbed = max(np.max(np.abs(stable_poly(spec+sign*err)-exact)) for sign in (-1, 1))
    assert perturbed <= bound*(1+1e-7)
    return dict(pair_pulse_area=t.sum(), single_slot_summed_area=2*t.sum(),
                minimum_pair_gap=gap, absolute_calibration_error=err,
                operator_error_bound=bound, perturbed_test_error=perturbed)


def spectrum_frame(y):
    u, s, _ = np.linalg.svd(y)
    return s[::-1]**2, u[:, ::-1]


def dark(v, theta):
    s, c = math.sin(theta), math.cos(theta)
    v0 = np.eye(4, dtype=complex)
    v0[:3, :3] = v
    rotation = np.eye(4)
    rotation[2:, 2:] = [[c, s], [-s, c]]
    ve = v0 @ rotation
    out = ve.conj().T
    idx = pairs(4).index((2, 3))
    amplitudes = exterior(out)[:, 0]
    hidden = sum(abs(amplitudes[j])**2 for j, p in enumerate(pairs(4)) if 3 in p)
    assert abs(amplitudes[idx]) < 1e-14
    np.testing.assert_allclose(hidden, s*s*sum(abs(v[:2,2])**2), atol=1e-16)
    flip_u = np.diag([-1, 1, 1, 1])
    flip_b = np.diag([1, 1, -1, 1])
    combined = out @ flip_u @ ve @ flip_b @ out
    detected = abs(exterior(combined)[idx, 0])**2
    expected = 4*math.sin(2*theta)**2*abs(v[0,2]*v[1,2])**2
    np.testing.assert_allclose(detected, expected, atol=1e-19, rtol=1e-10)
    return dict(theta_rad=theta, single_pass_any_fourth_probability=hidden,
                single_pass_b4_probability=0., known_phase_three_pass_b4_probability=expected,
                three_pass_independent_trials_95=attempts(expected) if theta else None)


def run():
    model = runpy.run_path(str(HERE / CFG["charged_snapshot"]))
    matrices, ref = model[CFG["snapshot_function"]]()
    hu, left_u = spectrum_frame(matrices[0])
    hd, left_d = spectrum_frame(matrices[1])
    v = left_u.conj().T @ left_d
    np.testing.assert_allclose(abs(v[0,1]),ref["Vus_at_upper_neutral_cluster"],rtol=1e-9)
    ev = exterior(v.conj().T)
    np.testing.assert_allclose(ev.conj().T @ ev,np.eye(3),atol=1e-13)
    prob = abs(ev[:, 0])**2
    np.testing.assert_allclose(prob,abs(v[2, ::-1])**2,atol=1e-14)
    classical = []
    for i, j in pairs(3):
        classical.append(abs(v[0,i]*v[1,j])**2 + abs(v[0,j]*v[1,i])**2)
    rows = []
    output_filters = {}
    for j, key in enumerate(CFG["labels"]):
        filt = pair_filter(hd,j)
        delta = filt["operator_error_bound"]
        good = prob[j]*(1-delta)**2
        bad = (1-prob[j])*delta**2
        filt["targeted_preparation_wrong_pair_fraction_bound"] = bad/(good+bad)
        output_filters[key] = filt
        rows.append(dict(pair=key, coherent_probability=prob[j],
                         distinguishable_unordered_different_probability=classical[j],
                         trials_for_95_percent_at_least_once=attempts(prob[j])))
    # Column and row order deliberately put the vanishing DFT minor at the b4 readout.
    idx = np.arange(4)
    dft = (1j**np.outer(idx,idx)/2)[np.ix_([0,2,1,3],[1,3,0,2])]
    assert np.min(abs(dft)) == .5 and abs(exterior(dft.conj().T)[-1,0]) == 0
    uniform = 3*min(prob)
    np.testing.assert_allclose(prob*(min(prob)/prob)/uniform, np.ones(3)/3)
    assert 0 < min(prob) < .001
    return dict(snapshot="1711 same 7171 GeV upper cluster, model Yukawa running",
                CKM_abs=np.abs(v).tolist(), source="u wedge c", conditional_ports=rows,
                uniform_three_pair_one_pass_acceptance=uniform,
                classical_same_flavour_probability=1-sum(classical),
                eigenvalues_up=hu.tolist(),eigenvalues_down=hd.tolist(),
                source_mass_filter=pair_filter(hu,0),
                output_mass_filters=output_filters,
                dark_rotations=[dark(v,t) for t in CFG["fourth_rotations_rad"]],
                DFT4=dict(all_single_probabilities=.25, coherent_b4=0,
                          distinguishable_b4=.125))


def rounded(x):
    if isinstance(x, dict):
        return {key:rounded(v) for key,v in x.items()}
    if isinstance(x, list):
        return [rounded(v) for v in x]
    if isinstance(x, (float,np.floating)):
        return float(f"{x:.9g}")
    return x


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    check = ap.parse_args().check
    output = rounded(run())
    path = HERE/"results.json"
    if check:
        assert json.loads(path.read_bytes()) == output
        print("1774 verified: 3 mass pairs, 4D dark port, known-label coherent rescue")
    else:
        path.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n", encoding="utf8")
