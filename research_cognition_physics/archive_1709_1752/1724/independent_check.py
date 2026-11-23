"""An independent Bs mixing readout of the tagged X radial scalar current."""
import argparse
import cmath
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
HAD = {125.25: (163., [-.30, .40, -.14, -.29]),
       1000.: (1000., [-.37, .51, -.17, -.33])}


def vertices():
    old = runpy.run_path(str(HERE / "../1715/check.py"))
    flow, rep = old["charged_background"](.1)
    down = old["unpack"](rep["down_texture_10PeV"])
    ul, y, vh = np.linalg.svd(down)
    ul, y, vr = ul[:, ::-1], y[::-1], vh.conj().T[:, ::-1]
    omega = cmath.exp(-2j * math.pi / 3)
    long = omega * down[1, 2].imag / omega.imag
    assert abs((down[1, 2] - long).imag) < 1e-14
    x = np.zeros_like(down)
    x[1, 2] = long
    g = ul.conj().T @ x @ vr / y[2]
    np.testing.assert_allclose(g, np.outer(ul[1].conj(), vr[2]) * long / y[2])
    top = ul[2, 1].conjugate() * ul[2, 2]
    remove = top.conjugate() / abs(top)
    v = flow["low_ckm"](rep["gamma_MZ_deg"], flow["angles"]("direct"))
    beta = cmath.phase(-v[2, 1] * v[2, 2].conjugate() / (v[1, 1] * v[1, 2].conjugate()))
    return g[1, 2] * remove, g[2, 1].conjugate() * remove, beta


def contractions(mass, alpha_mz, factors):
    mu, had = HAD[mass]
    inverse = 1 / alpha_mz + (23 / 3 * math.log(163. / 91.1876)
                              + 7 * math.log(mu / 163.)) / (2 * math.pi)
    a = 1 / (4 * math.pi * inverse)
    log = 2 * math.log(mass / mu)
    coeff = np.array([-1.5 * a, 1 - a, 1 + (4.5 - 3 * log) * a,
                      (.125 - log / 12) * a])
    weighted = coeff * had * factors
    return sum(weighted[:2]), sum(weighted[2:])


def amplitude(mass, r, l):
    lr, sll = contractions(mass, .118, np.ones(4))
    return -(sll * (r**2 + l**2) / 2 + lr * r * l) / mass**2


def robust_lower(mass, r, l):
    # Use the signed dominant RR term. The sum of the other two moduli
    # bounds every cancellation even if their phases are completely free.
    had_box = itertools.product((.8, 1.2), repeat=4)
    contractions_box = [contractions(mass, alpha, fac)
                        for fac in had_box for alpha in (.116, .120)]
    # Each expression is multi-affine in alpha(mu) and the four hadronic factors.
    lo_s = min(-s for _, s in contractions_box)
    hi_s = max(abs(s) for _, s in contractions_box)
    hi_lr = max(abs(lr) for lr, _ in contractions_box)
    center = cmath.phase(r**2)
    width = math.radians(10.)
    assert 0 < center - width < center + width < math.pi
    sine = min(math.sin(center - width), math.sin(center + width))
    lower = (lo_s / 2 * (.9 * abs(r))**2 * sine
             - hi_s / 2 * (1.1 * abs(l))**2
             - hi_lr * 1.1**2 * abs(r * l)) / mass**2
    assert lower > 0
    return lower


def pair(z):
    return [z.real, z.imag]


def results():
    r, l, beta = vertices()
    freq, hbar = 17.7683, 6.582119569e-13
    freq_high = freq + 2 * math.hypot(.0051, .0032)
    mix_high = freq_high * hbar / math.sqrt(4 - .02**2)
    # A budget on decay/absorptive phase and independent non-top mixing
    # perpendicular to the top current. No SM box magnitude is necessary.
    phase = abs(-.039 + 2 * beta) + 2 * math.hypot(.022, .006) + .01
    bound = mix_high * (math.sin(phase) + .01)
    tables = {}
    for mass in HAD:
        a = amplitude(mass, r, l)
        lower = robust_lower(mass, r, l)
        tables[str(mass)] = {"central_M12star_GeV_per_gain_squared": pair(a),
                             "robust_imag_lower_GeV_per_gain_squared": lower,
                             "central_coeff_gain_upper": math.sqrt(bound / a.imag),
                             "robust_gain_upper": math.sqrt(bound / lower)}
    return {"top_rephased_vertices_R_L_per_gain": [pair(r), pair(l)],
            "absolute_chirality_ratio": abs(r/l), "SM_phi_s_from_low_quark_point": -2 * beta,
            "upper_total_M12_GeV": mix_high, "phase_offset_upper_rad": phase,
            "perpendicular_budget_GeV": bound, "poles_GeV": tables,
            "assumptions": {"absolute_Gamma12_over_M12": .02,
                            "other_decay_and_absorptive_phase_rad": .01,
                            "other_perpendicular_mixing_over_total_upper": .01,
                            "vertex_modulus_multipliers": [.9, 1.1],
                            "right_vertex_phase_halfwidth_degrees": 5,
                            "hadronic_element_modulus_multipliers": [.8, 1.2],
                            "alpha_s_MZ_interval": [.116, .120]},
            "source_hadronic_NDR_tables_1_2_and_Wilson": "https://arxiv.org/pdf/1303.3723",
            "source_oscillation_frequency": "https://arxiv.org/abs/2104.04421",
            "source_decay_phase": "https://arxiv.org/abs/2308.01468"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    out = results()
    path = HERE / "independent_results.json"
    if args.check:
        assert out == json.loads(path.read_bytes())
        print("1724 independent check passed")
    else:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n", encoding="utf-8")
