"""Check the two-doublet spectrum and independent T -> th,tZ amplitudes."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOP, HIGGS, Z, VEV = 172.5, 125.09, 91.1876, 246.22
LEFT_WEAK_LOSS = .001


def gammas():
    eye, zero = np.eye(2), np.zeros((2, 2))
    matrices = [np.block([[eye, zero], [zero, -eye]])]
    pauli = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    matrices.extend(np.block([[zero, x], [-x, zero]]) for x in pauli)
    return matrices


GAMMA = gammas()
IDENTITY = np.eye(4)
G5 = 1j * GAMMA[0] @ GAMMA[1] @ GAMMA[2] @ GAMMA[3]
LEFT, RIGHT = (IDENTITY - G5) / 2, (IDENTITY + G5) / 2


def phase(m, daughter):
    return math.sqrt((m*m - (TOP + daughter)**2) * (m*m - (TOP - daughter)**2))


def ratio_formula(m):
    top2, mass2, z2, h2 = TOP**2, m**2, Z**2, HIGGS**2
    hs = (mass2 + top2) * (mass2 + top2 - h2) + 4 * mass2 * top2
    zs = (mass2 - top2)**2 + z2 * (mass2 + top2 - 2*z2)
    return phase(m, HIGGS) / phase(m, Z) * hs / zs


def spin_trace_width(m, daughter, vertices):
    momentum = phase(m, daughter) / (2*m)
    e_top, e_boson = math.hypot(momentum, TOP), math.hypot(momentum, daughter)
    initial = m * (GAMMA[0] + IDENTITY)
    final = e_top * GAMMA[0] - momentum * GAMMA[3] + TOP * IDENTITY
    ops = vertices
    if len(vertices) == 4:
        # Boson momentum is along -z. Contract covariant polarization vectors.
        ops = [vertices[1], vertices[2], (momentum * vertices[0] + e_boson * vertices[3]) / daughter]
    summed = 0.
    for op in ops:
        adjoint = GAMMA[0] @ op.conj().T @ GAMMA[0]
        total = np.trace(final @ op @ initial @ adjoint) / 2
        assert abs(total.imag) < 1e-9
        summed += total.real
    return momentum * summed / (8 * math.pi * m*m)


def svd_spectrum(bottom_mass):
    kappa = 1 - LEFT_WEAK_LOSS * bottom_mass**2 / TOP**2
    partner = math.sqrt(TOP**2 + (bottom_mass**2 - TOP**2) / kappa)
    naked_y = math.sqrt(TOP**2 + partner**2 - bottom_mass**2)
    s = TOP * partner / (bottom_mass * naked_y)
    c = math.sqrt(1 - s*s)
    unit_radial = np.array([[naked_y*s, 0], [naked_y*c, 0]])
    matrix = unit_radial + np.diag([0, bottom_mass])
    l, masses, rt = np.linalg.svd(matrix)
    l, masses, r = l[:, ::-1], masses[::-1], rt.T[:, ::-1]
    assert np.allclose(masses, [TOP, partner])
    coupling = l.T @ unit_radial @ r / VEV
    projection = r.T @ np.diag([1, 0]) @ r
    assert np.allclose(coupling * VEV, np.diag(masses) @ projection)
    assert math.isclose(kappa, projection[0, 0], rel_tol=2e-13)
    assert math.isclose(1 - l[0, 0]**2, LEFT_WEAK_LOSS, rel_tol=1e-11)
    assert np.allclose(l.T @ np.eye(2) @ l, np.eye(2))
    right_isospin = r.T @ np.diag([0, .5]) @ r
    scalar = coupling[0, 1] * RIGHT + coupling[1, 0] * LEFT
    vector = [(2*Z/VEV) * right_isospin[0, 1] * g @ RIGHT for g in GAMMA]
    width_h = spin_trace_width(partner, HIGGS, [scalar])
    width_z = spin_trace_width(partner, Z, vector)
    ratio = ratio_formula(partner)
    assert math.isclose(width_h / width_z, ratio, rel_tol=5e-13)
    # Orthogonal change of the two unmarked doublet labels leaves right-handed Gram fixed.
    turn = np.array([[.6, -.8], [.8, .6]])
    assert np.allclose((turn @ matrix).T @ (turn @ matrix), matrix.T @ matrix)
    assert np.linalg.norm(turn @ l[:, 0] - l[:, 0]) > .5
    return {"B_GeV": bottom_mass, "T_GeV": partner, "site_weight_s_squared": s*s,
            "tree_top_kappa": kappa, "conditional_tth_rate": kappa*kappa,
            "right_isospin_top": right_isospin[0, 0],
            "Gamma_th_tree_GeV": width_h, "Gamma_tZ_tree_GeV": width_z,
            "partial_width_ratio_th_over_tZ": ratio}


def results():
    assert abs(ratio_formula(1e6) - 1) < 4e-7
    return {"case_1200": svd_spectrum(1200.), "case_1600": svd_spectrum(1600.),
            "ATLAS_working_positive_kappa_low": math.sqrt(.81-.19),
            "large_mass_ratio_check_1e6_GeV": ratio_formula(1e6)}


def compare(a, b):
    if isinstance(a, dict):
        assert set(a) == set(b)
        for name in a:
            compare(a[name], b[name])
    else:
        assert math.isclose(a, b, rel_tol=3e-10, abs_tol=1e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = results()
    path = HERE / "independent_results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(result))
