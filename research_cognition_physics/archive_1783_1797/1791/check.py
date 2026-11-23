"""Leading inert-scalar loop, common scalar pole spectrum; no fit packages."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent


def kernels(masses, invariant, order):
    """dB0/ds, B0(s)-B0(0), equal-on-shell C0 (positive Feynman integral)."""
    nodes, w = leggauss(order)
    u, w = (nodes + 1) / 2, w / 2
    xy = u[:, None]
    yz = (1 - xy) * u[None, :]
    zx = 1 - xy - yz
    measure = (w * (1 - u))[:, None] * w[None, :]
    shape = xy * yz + yz * zx + zx * xy
    out = []
    for m in masses:
        z = u * (1 - u)
        out.append([
            np.dot(w, z / (m**2 - invariant * z)),
            -np.dot(w, np.log1p(-invariant * z / m**2)),
            np.sum(measure / (m**2 - invariant * shape)),
        ])
    return np.array(out).T


def cut(x, data, order):
    # Reduced helicity partial waves; the common (32 pi)^-2 factor is restored.
    z, w = leggauss(order)
    mass, h = data["charged_mass_GeV"], data["higgs_mass_GeV"]
    s, v = data["photon_Ecm_GeV"] ** 2, data["vev_GeV"]
    beta, bh = np.sqrt(1 - 4 * mass**2 / s), np.sqrt(1 - 4 * h**2 / s)
    q = beta * bh * s / 2
    k = s / 2 - h**2
    amp_photon = 8 * np.pi * data["alpha_em"] * beta**2 * (1 - z**2) / (1 - beta**2 * z**2)
    amp_hh = (v * x) ** 2 * (1 / (k - q * z) + 1 / (k + q * z))
    d20, p2 = np.sqrt(3 / 8) * (1 - z**2), (3 * z**2 - 1) / 2
    partial_photon = np.dot(w, d20 * amp_photon) / (32 * np.pi)
    partial_higgs = np.dot(w, p2 * amp_hh) / (32 * np.pi)
    return float(2 * beta * partial_photon * partial_higgs)


def evaluate(data, order):
    v, mh = data["vev_GeV"], data["higgs_mass_GeV"]
    charged, odd = data["charged_mass_GeV"], data["odd_neutral_mass_GeV"]
    even = data["even_neutral_mass_GeV"]
    masses = np.array([charged, odd, even[0], charged, odd, even[1]])
    n = np.array([2, 1, 1, 2, 1, 1])
    sign = np.array([1, 1, 1, -1, -1, -1])
    offsets = 2 * (masses**2 - charged**2) / v**2
    bp, db, tri = kernels(masses, mh**2, order)
    weight = n * bp
    a = np.dot(weight, sign**2)
    b = np.dot(weight, offsets * sign)
    center = -b / a
    coeff = 32 * np.pi**2
    tree = 3 * mh**2 / v
    cases = []
    for x in [center - data["portal_half_separation"], center + data["portal_half_separation"]]:
        eta = offsets + x * sign
        g = eta * v
        residue_defect = np.dot(weight, g**2) / coeff
        effective = np.sum(n * g**3 / masses**2) / coeff
        triangle = np.sum(n * g**3 * tri) / (16 * np.pi**2)
        pole_mass_term = 3 * np.dot(n * g**2, db) / (v * coeff)
        vertex_bubble = -3 * np.dot(n * eta * g, db) / coeff
        min_portals = np.array([min(eta[:3]), min(eta[3:])])
        margin = data["inert_quartic"] - np.dot(min_portals, min_portals) / (4 * (mh**2 / (2*v**2)))
        cases.append(dict(eta_charged=[float(x), float(-x)], bare_mu2_GeV2=(charged**2 - v**2/2*np.array([x, -x])).tolist(),
                          all_radial_portals=eta.tolist(), stability_margin=float(margin),
                          pole_residue_defect=float(residue_defect), charged_one_H_portal_sum=float(g[0] + g[3]),
                          curvature_inert_shift_over_tree=float(effective/tree),
                          equal_pole_triangle_over_tree=float(triangle/tree),
                          equal_pole_inert_vertex_over_tree_before_common_v_match=float(triangle/tree - 1.5*residue_defect),
                          pole_mass_term_GeV=float(pole_mass_term), vertex_bubble_GeV=float(vertex_bubble),
                          charged_cut_Im_a2=float(cut(x, data, order))))
    p, r = cases
    assert abs(p["pole_residue_defect"] - r["pole_residue_defect"]) < 1e-14
    for c in cases:
        assert c["stability_margin"] > 0 and c["charged_one_H_portal_sum"] == 0
        assert abs(c["pole_mass_term_GeV"] + c["vertex_bubble_GeV"]) < 1e-11
        assert c["charged_cut_Im_a2"] > 0
    assert np.isclose(r["charged_cut_Im_a2"]/p["charged_cut_Im_a2"], (r["eta_charged"][0]/p["eta_charged"][0])**2)
    weight_limit = p["pole_residue_defect"] ** 1.5 * np.sqrt(coeff) / tree
    cubic_weight = 1 / (np.sqrt(n) * bp**1.5)
    uniform_bounds = dict(curvature=float(weight_limit * max(cubic_weight / masses**2)),
                          pole_triangle=float(weight_limit * max(2 * cubic_weight * tri)))
    assert all(abs(c["equal_pole_triangle_over_tree"]) <= uniform_bounds["pole_triangle"] for c in cases)
    return dict(portal_symmetry_center=float(center), masses_C_A_H_GeV=masses.tolist(),
                all_portal_absolute_bounds_over_SM_tree_at_same_residue=uniform_bounds,
                cases=cases, curvature_difference_over_SM_tree=r["curvature_inert_shift_over_tree"]-p["curvature_inert_shift_over_tree"],
                three_pole_residue_difference_over_SM_tree=r["equal_pole_triangle_over_tree"]-p["equal_pole_triangle_over_tree"],
                cut_ratio_second_over_first=r["charged_cut_Im_a2"]/p["charged_cut_Im_a2"])


def rounded(x):
    if isinstance(x, dict):
        return {k: rounded(v) for k, v in x.items()}
    if isinstance(x, list):
        return [rounded(v) for v in x]
    return float(f"{x:.11g}") if isinstance(x, float) else x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = HERE / "inputs.json"
    data = json.loads(source.read_text(encoding="utf-8"))
    answer = rounded(evaluate(data, data["quadrature_order"]))
    assert answer == rounded(evaluate(data, data["quadrature_order"] * 2))
    answer["inputs_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    dest = HERE / "results.json"
    result = json.dumps(answer, indent=2, ensure_ascii=False) + "\n"
    if args.check:
        assert json.loads(dest.read_text(encoding="utf-8")) == answer
        print("1791 OK", hashlib.sha256(result.encode()).hexdigest()[:12])
    else:
        dest.write_text(result, encoding="utf-8")


if __name__ == "__main__":
    main()
