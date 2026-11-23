"""Test two visible charged-lepton lines against d and s proton-decay ports."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
# Table VIII of arXiv:2111.01608. Mixed-chirality W0,W1 magnitudes and errors.
# The pion isospin conversion is applied after loading the charged-pion results.
FORMS = {
    "pi_e": ([.159, .015, .020, .025], [.169, .014, .018, .029]),
    "pi_mu": ([.161, .015, .020, .026], [.170, .014, .018, .028]),
    "K_e": ([.0854, .0057, .0055, .0090], [.0258, .0038, .0003, .0004]),
    "K_mu": ([.0860, .0056, .0055, .0091], [.0259, .0038, .0004, .0002]),
}
YEARS = {"pi_e": 2.4e34, "pi_mu": 1.6e34, "K_e": 1e33, "K_mu": 3.6e33}
MP, MPI, MK, ME, MMU = .938272, .134977, .497611, .000511, .105658


def flavour():
    old = runpy.run_path(str(HERE / "../1715/check.py"))
    _, point = old["charged_background"](.1)
    D = old["unpack"](point["down_texture_10PeV"])
    E = old["unpack"](point["charged_texture_10PeV"]).T
    # Spectral projectors independent of both singular-vector phase conventions.
    _, down_basis = np.linalg.eigh(D @ D.conj().T)
    _, lepton_basis = np.linalg.eigh(E @ E.conj().T)
    weighted = np.diag([2, 1, 1]) @ down_basis[:, :2]
    gram = weighted.conj().T @ weighted
    np.testing.assert_allclose(gram, np.eye(2) + 3*np.outer(
        down_basis[0, :2].conj(), down_basis[0, :2]), atol=1e-14)
    eig = np.linalg.eigvalsh(gram)
    np.testing.assert_allclose(eig[0], 1, atol=2e-14)
    # When tau is chosen along the largest left singular vector, the lower
    # bound on the sum of the four light couplings is attained exactly.
    left = np.linalg.svd(weighted, full_matrices=True)[0]
    np.testing.assert_allclose(np.linalg.norm(left[:, 1:].conj().T @ weighted)**2, 1)
    rows = {}
    for j, meson in enumerate(("pi", "K")):
        for k, lepton in enumerate(("e", "mu")):
            terms = lepton_basis[:, k].conj() * weighted[:, j]
            lengths = abs(terms)
            rows[meson + "_" + lepton] = {
                "tagged": abs(sum(terms)),
                "min_over_common_family_phases": max(0, 2*max(lengths)-sum(lengths)),
                "max_over_common_family_phases": sum(lengths)}
    return gram, eig, rows


def coefficients(key, conservative):
    meson, lepton = key.split("_")
    m, ml = (MPI if meson == "pi" else MK), (ME if lepton == "e" else MMU)
    e = (MP**2 + ml**2 - m**2)/(2*MP)
    p = math.sqrt(e**2-ml**2)
    w0, w1 = FORMS[key]
    a0 = w0[0] - conservative*sum(w0[1:])
    a1 = w1[0] + conservative*sum(w1[1:])
    form = (a0 - ml/MP*a1)/(math.sqrt(2) if meson == "pi" else 1)
    assert form > 0
    k = p*(e-ml)/(8*math.pi*MP) * form**2
    upper = 6.582119569e-25/(31557600*YEARS[key])
    return math.sqrt(upper/k)


def scenario(gram, tagged, conservative, g2, crossing):
    caps = {key: coefficients(key, conservative) for key in YEARS}
    # Sum positive constraints with weights (1,w) on the two columns. Every
    # tested w separately certifies its numeric bound. Optimality is not needed.
    a = caps["pi_e"]**2 + caps["pi_mu"]**2
    b = caps["K_e"]**2 + caps["K_mu"]**2
    def common_bound(w):
        s = np.diag([1., math.sqrt(w)])
        floor = np.linalg.eigvalsh(s @ gram @ s)[0]
        return math.sqrt((a + w*b)/floor)
    candidates = np.geomspace(.025, 40., 2049)
    best_weight = min(candidates, key=common_bound)
    upper = common_bound(best_weight)
    fixed = {key: math.sqrt(g2*row["min_over_common_family_phases"]/caps[key])
             for key, row in tagged.items()}
    return {"each_low_energy_coupling_upper_GeVm2": caps,
            "minimal_MX_GeV_common_phases_unit_transfer_by_mode": fixed,
            "common_coefficient_upper_any_flavour_rotation_unweighted_GeVm2": common_bound(1.),
            "positive_K_over_pi_weight_used": best_weight,
            "common_coefficient_upper_any_flavour_rotation_GeVm2": upper,
            "arbitrary_lepton_frame_MX_floor_at_unit_transfer_GeV": math.sqrt(g2/upper),
            "max_transfer_at_gauge_crossing_any_lepton_frame": upper*crossing**2/g2}


def results():
    gram, eig, tagged = flavour()
    gauge = json.loads((HERE / "../1720/results.json").read_bytes())
    g2 = 4*math.pi/gauge["central_conditional_combined_point"]["common_inverse_coupling"]
    crossing = gauge["baseline_unification_GeV"]
    return {"two_column_gram_eigenvalues": eig.tolist(),
            "min_visible_e_mu_d_s_squared_sum_arbitrary_frame": 1,
            "fixed_tags_amplitudes": tagged, "separate_partial_lifetime_limits_year": YEARS,
            "lattice_central": scenario(gram, tagged, 0, g2, crossing),
            "each_lattice_error_one_width": scenario(gram, tagged, 1, g2, crossing),
            "g_squared": g2, "crossing_GeV": crossing,
            "sources": ["https://arxiv.org/pdf/hep-ph/0403286",
                        "https://arxiv.org/pdf/2111.01608",
                        "https://arxiv.org/abs/2010.16098",
                        "https://arxiv.org/abs/hep-ex/0502026",
                        "https://arxiv.org/abs/2208.13188"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    path = HERE / "independent_results.json"
    data = json.loads(json.dumps(results()), parse_float=lambda z: float(f"{float(z):.11g}"))
    if args.check:
        assert data == json.loads(path.read_bytes())
        print("1725 independent check passed")
    else:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
