"""Custodial self-energy from the physical five-singlet spectrum, with sharp CKM envelopes."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE / relative).read_bytes())


def physical_endpoint():
    structure = runpy.run_path(str(HERE / "../1713/check.py"))
    down, _, _ = structure["actual_target"]()
    link, _ = structure["generated_alignment"](down)
    ordinary, values, _ = np.linalg.svd(down)
    scale = values[0]
    signs, embed = structure["dilation"](down @ link, scale)
    return ordinary[:, ::-1], link, signs, embed, scale


def charged_frame():
    ck = read("../1702/inputs.json")["direct_semileptonic_moduli"]
    point = read("../1712/results.json")["cases"]["interior_all_one_width"]
    pulls = point["input_widths"]
    us, ub, cb = [ck[name][0] + pulls[name] * math.hypot(*ck[name][1:]) for name in ("Vus", "Vub", "Vcb")]
    c13 = math.sqrt(1 - ub**2)
    builder = runpy.run_path(str(HERE / "../1702/check.py"))["low_ckm"]
    return builder(point["quark_neutrino"]["quarks"]["gamma_deg"], np.arcsin([us/c13, ub, cb/c13]))


def spectrum(g, bare, vev, endpoint, chiral=False):
    basis, link, a, w, t = endpoint
    kappa = 0. if chiral else t / math.sqrt(g**2 - t**2)
    mass = np.block([[np.zeros((3, 3)), g*vev/math.sqrt(2)*w.conj().T],
                     [-bare*kappa*w@link.conj().T, bare*a]])
    left, masses, _ = np.linalg.svd(mass)
    return masses[::-1], basis.conj().T @ left[:3, ::-1]


def calculate():
    cfg = read("../1642/inputs.json")
    base = runpy.run_path(str(HERE / "../1642/check.py"))
    ker = base["kernels"]
    norm = 3/(16*math.pi**2*cfg["vev_GeV"]**2)
    up = cfg["up_GeV"]
    mass0 = read("../1713/inputs.json")["mediator_matching_GeV"]
    endpoint = physical_endpoint()
    ck = charged_frame()
    vev = cfg["vev_GeV"]
    permutations = list(itertools.permutations(range(3)))

    def delta(g, bare, chiral=False):
        masses, weak = spectrum(g, bare, vev, endpoint, chiral)
        neutral = weak.conj().T @ weak
        np.testing.assert_allclose(neutral @ neutral, neutral, atol=4e-15)
        block = np.array([[ker(u, d)[0] for d in masses] for u in up])
        z = sum(abs(neutral[i, j])**2 * ker(masses[i], masses[j])[0]
                for j in range(8) for i in range(j))
        full = (np.sum(abs(ck @ weak)**2 * block) - z - np.sum(abs(ck)**2 * block[:, :3])) * norm
        if chiral:
            # Degenerate light states can rotate freely. The subtraction still uses
            # only their common zero mass and is independent of that rotation.
            return full
        operators = [(weak * row) @ weak.conj().T for row in block]
        # H_i includes a sum over every heavy mass, even within a degenerate cluster.
        diagonal = np.array([h.diagonal().real for h in operators]) - block[:, :3]
        assert max(np.linalg.norm(h-np.diag(h.diagonal()), 2) for h in operators) < 3e-5
        resolved = (np.sum(abs(ck)**2*diagonal) - z)*norm
        np.testing.assert_allclose(full, resolved, atol=3e-12)
        extremes = np.array([(sum(diagonal[i,p[i]] for i in range(3)) - z)*norm for p in permutations])
        assert extremes.min()-1e-13 <= full <= extremes.max()+1e-13
        return full, extremes

    def closed(g, bare):
        x = g * vev / math.sqrt(2)
        h2 = bare**2 + x*x
        eps = x*x / h2
        masses2 = np.square(up)
        log_sum = np.sum(masses2 * h2 / (h2 - masses2) * np.log(h2 / masses2))
        return norm * (3*eps**2*h2 - 2*eps*log_sum)

    def point(g):
        finite, ext = delta(g, mass0)
        simple = closed(g, mass0)
        np.testing.assert_allclose(delta(g, mass0, chiral=True), simple, atol=8e-12)
        return {"coupling":g, "full_delta_rho": finite, "full_T": finite*cfg["alpha_inverse"],
                "all_orientations_delta_rho_range": [float(ext.min()), float(ext.max())],
                "zero_down_mass_analytic_delta_rho": simple}

    lo, hi = 1., 3.
    for _ in range(45):
        mid = (lo+hi)/2
        if delta(mid, mass0)[0] > 0:
            hi = mid
        else:
            lo = mid
    root = (hi+lo)/2
    examples = [point(g) for g in (.05, 1., 1.88297623, 2., 3.)]
    assert examples[1]["all_orientations_delta_rho_range"][1] < 0
    assert examples[-1]["all_orientations_delta_rho_range"][0] > 0
    np.testing.assert_allclose(delta(root, mass0)[0], 0, atol=3e-12)
    return {"bare_M_GeV": mass0, "reference_up_GeV": up, "vev_GeV": vev,
            "examples": examples, "zero_of_full_rho_in_registered_CKM_frame_g": root,
            "fixed_g1_double_and_quadruple_M_chiral_delta_rho": [closed(1., mass0*i) for i in (1,2,4)],
            "rho_with_same_chiral_weak_projection_and_scaled_M_and_g":
                [closed(i, mass0*i) for i in (1,2,3)]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    value = json.loads(json.dumps(calculate()), parse_float=lambda x:float(f"{float(x):.10g}"))
    out = HERE / "independent_results.json"
    if args.check:
        assert json.loads(out.read_bytes()) == value
        print("1719 verified electroweak vacuum polarization, full mass spectrum, sharp unitary envelopes and chiral limits")
    else:
        with out.open("x", encoding="utf8") as stream:
            stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
