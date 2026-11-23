"""Electroweak breaking in the R1713 five-singlet mass matrix and the kaon axial/vector test."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())


def read(path):
    return json.loads((HERE / path).read_bytes())


def kaon_test(data, pion):
    def interval(k):
        d, s, a, _ = data["rectangle"](k, pion)
        return max(s[0]/d[1], a[0]), min(s[1]/d[0], a[1])

    lo, hi = 0., 4.
    for _ in range(64):
        t = (lo + hi) / 2
        lower, upper = interval(t)
        if lower > upper:
            lo = t
        else:
            hi = t
    d, s, a, _ = data["rectangle"](1., pion)
    ratio = [s[0]/d[1], s[1]/d[0]]
    if pion:
        assert max(ratio[0], a[0]) < min(ratio[1], a[1])
    else:
        assert ratio[1] < a[0]
    return {"vector_ratio_one_width": ratio, "axial_ratio_one_width": a,
            "first_overlap_source_width_multiplier": hi}


def diagonalize(coupling, a, w, link, singular, framed, bare, vev):
    gain = singular[-1]
    kappa = gain / math.sqrt(coupling**2 - gain**2)
    x = vev * coupling
    full = np.block([[np.zeros((3, 3)), x*w.conj().T],
                     [-bare*kappa*w@link.conj().T, bare*a]])
    left, mass, right_h = np.linalg.svd(full)
    light = left[:, ::-1][:, :3]
    masses = mass[::-1][:3]
    # framed consists of light-down eigenvectors in the common ordinary left space.
    relative = abs(framed.conj().T @ light[:3])
    ell = relative.diagonal()
    neutral = light[:3].conj().T @ light[:3]
    np.testing.assert_allclose(relative, np.diag(ell), atol=2e-9)
    np.testing.assert_allclose(neutral, np.diag(ell**2), atol=6e-14)
    # With all v-dependence in the ordinary left rows, the Higgs vertex
    # in the exact light mass basis is the neutral-current projector times m.
    vr = right_h.conj().T[:, ::-1][:, :3]
    top = np.zeros_like(full)
    top[:3] = full[:3]
    higgs = light.conj().T @ top @ vr
    np.testing.assert_allclose(higgs, neutral * masses, atol=3e-9)
    # Independent symmetric logarithmic derivative of *physical* light masses.
    delta = 3e-4
    varied = [np.linalg.svd(full + np.expm1(t)*top, compute_uv=False)[::-1][:3]
              for t in (-delta, delta)]
    derivative = np.log(varied[1] / varied[0]) / (2*delta)
    np.testing.assert_allclose(derivative, ell**2, rtol=2e-6)
    heavy = left[:, :len(a)]
    weights = np.sum(abs(heavy[:3])**2, axis=0)
    np.testing.assert_allclose(weights.sum(), np.sum(1 - ell**2), atol=4e-14)
    tau_heavy = (CFG["Higgs_mass_GeV"] / (2*mass[:len(a)]))**2
    assert np.max(tau_heavy) < 1e-4
    shape = (singular/gain)**2
    q = (masses/bare)**2
    formula = 1/np.sqrt(1 + (x/bare)**2 *
                         ((1-shape)/(1-q)**2 + shape/(1+kappa**2-q)**2))
    np.testing.assert_allclose(ell, formula, atol=5e-15, rtol=5e-15)
    return {"g": coupling, "kappa": kappa, "doublet_projections_d_s_b": ell.tolist(),
            "finite_vev_masses_over_leading": (masses/(vev*singular)).tolist(),
            "neutral_current_left_shifts_d_s_b": ((1-ell**2)/2).tolist(),
            "physical_higgs_diagonal_over_m_over_v": (ell**2).tolist(),
            "new_heavy_ggH_amplitude_over_one_infinite_mass_SM_fundamental": [
                float(weights.sum()), float(np.sum(weights/(1-tau_heavy)))],
            "largest_offdiagonal_Z_projector": float(abs(neutral - np.diag(neutral.diagonal())).max())}


def compute():
    src = runpy.run_path(str(HERE / CFG["mass_model"]))
    data = runpy.run_path(str(HERE / CFG["data"]))
    y, _, _ = src["actual_target"]()
    link, _ = src["generated_alignment"](y)
    eigenframe, values, _ = np.linalg.svd(y)
    singular, framed = values[::-1], eigenframe[:, ::-1]
    a, w = src["dilation"](y @ link, singular[-1])
    bare = read(CFG["matching_scale_and_inputs"])["mediator_matching_GeV"]
    vev = read(CFG["vev"])["v246_GeV"]/math.sqrt(2)
    trials = [diagonalize(g, a, w, link, singular, framed, bare, vev)
              for g in CFG["counterfactual_g_values"]]
    baseline = trials[0]
    point = read(CFG["charged_point"])["cases"]["interior_all_one_width"]
    us = point["quark_neutrino"]["Vus_MZ"]
    ub = data["limits"](data["CKM"]["Vub"], point["input_widths"]["Vub"])[1]
    core = np.array([math.sqrt(1-us**2-ub**2), us, ub])
    moduli = core * baseline["doublet_projections_d_s_b"]
    # For a pure left weak current its absolute normalization drops out of the two ratios.
    right = read("../1713/independent_results.json")["branches"]["nuclear_beta_b_minus"]
    L = right["one_width_real_competing_current_witness"]["left_normalization"]
    candidate_g = bare/vev * math.sqrt(1/L**2 - 1)
    alternative = diagonalize(candidate_g, a, w, link, singular, framed, bare, vev)
    return {"bare_mass_GeV": bare, "example_vev_over_sqrt_two_GeV": vev,
            "cases": trials, "frozen_row_without_running_new_operators": moduli.tolist(),
            "frozen_row_defect": 1 - float(moduli @ moduli),
            "auxiliary_g_for_common_normalization_limit_at_kappa_zero": candidate_g,
            "alternative_tree_projection": alternative,
            "nuclear_light_rates": kaon_test(data, False), "pion_light_rates": kaon_test(data, True)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = compute()
    result = json.loads(json.dumps(data), parse_float=lambda n: float(f"{float(n):.10g}"))
    path = HERE / "results.json"
    if args.check:
        assert json.loads(path.read_bytes()) == result
        print("1717 verified exact light W/Z projections and the independent axial/vector width test")
    else:
        with path.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
