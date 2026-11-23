"""Rebuild the conditional two-branch Bd/Bs frequency-ratio certificate."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np


def window(center, error):
    return (center - error, center + error)


def quotient(x, y):
    assert min(x) >= 0 and min(y) > 0
    return (x[0] / y[1], x[1] / y[0])


def product(x, y):
    values = [i * j for i in x for j in y]
    return (min(values), max(values))


def summed(x, y):
    return (x[0] + y[0], x[1] + y[1])


def complement(s):
    return math.sqrt(1 - s * s)


def components(k, w, rho):
    s13 = rho * w
    c13 = complement(s13)
    s23, s12 = w / c13, k / c13
    c23, c12 = complement(s23), complement(s12)
    h = c12 * s23 * s13 / (s12 * c23)
    return s12, c12, s23, c23, s13, c13, h


def ratio(k, w, rho, gamma):
    s, c, sb, cb, st, ct, h = components(k, w, rho)
    delta = gamma + math.asin(h * math.sin(gamma))
    x = rho * cb * ct
    cross = 2 * s * c * x * math.cos(delta)
    square_ratio = (s*s + (c*x)**2 - cross) / (c*c + (s*x)**2 + cross)
    phase = cmath.exp(1j * delta)
    matrix = np.array([
        [c*ct, s*ct, st/phase],
        [-s*cb-c*sb*st*phase, c*cb-s*sb*st*phase, sb*ct],
        [s*sb-c*cb*st*phase, -c*sb-s*cb*st*phase, cb*ct],
    ])
    return square_ratio, matrix


def certificate():
    k = quotient(window(.2157, .0006), window(.9704, .0032))
    w_max = (.0351 + .0019 + .0018 + .0008) / .80
    rho = window(.083, .008)
    xi = window(1.206, .024)
    mass = quotient(window(5279.58, .43), window(5366.90, .51))
    short_range = quotient(window(1., .01), window(1., .01))
    scale = product(quotient(mass, product(xi, xi)), short_range)
    obs = quotient(window(.5050, .0031), window(17.7683, .0083))
    gamma_degrees = (62.6, 74.9)
    st_max = rho[1] * w_max
    ct_min = complement(st_max)
    sb_max = w_max / ct_min
    cb_min = complement(sb_max)
    s = (k[0], k[1] / ct_min)
    c = (complement(s[1]), complement(s[0]))
    x = (rho[0] * cb_min * ct_min, rho[1])
    h_max = c[1] * sb_max * st_max / (s[0] * cb_min)
    pad = math.degrees(math.asin(h_max))

    def branch(offset):
        delta = (gamma_degrees[0] + offset - pad, gamma_degrees[1] + offset + pad)
        assert ((0 <= delta[0] < delta[1] <= 90) or
                (180 <= delta[0] < delta[1] <= 270))
        cosines = tuple(sorted(math.cos(math.radians(v)) for v in delta))
        cross = product(product(product((2, 2), s), product(c, x)), cosines)
        cx, sx = product(c, x), product(s, x)
        numerator = summed(summed(product(s, s), product(cx, cx)), (-cross[1], -cross[0]))
        denominator = summed(summed(product(c, c), product(sx, sx)), cross)
        return product(quotient(numerator, denominator), scale)

    allowed, mirror = branch(0), branch(180)
    central_k, central_w = .2157 / .9704, .0351 / .93
    central_mass = 5279.58 / 5366.90
    central_scale = central_mass / 1.206**2
    observed_center = .5050 / 17.7683
    # One simultaneous witness; also supply a fixed round-number choice prior to inversion.
    fixed = {}
    for gamma in (65., 68.7, 245., 248.7):
        value, _ = ratio(central_k, central_w, .083, math.radians(gamma))
        fixed[str(gamma)] = value * central_scale
    s0, c0, _, cb, _, ct, h = components(central_k, central_w, .083)
    x0 = .083 * cb * ct
    z = observed_center / central_scale
    cosine = (s0**2 + (c0*x0)**2 - z*(c0**2 + (s0*x0)**2)) / (2*s0*c0*x0*(1+z))
    delta = math.acos(cosine)
    gamma = delta - cmath.phase(1 + h * cmath.exp(1j * delta))
    value, _ = ratio(central_k, central_w, .083, gamma)
    assert gamma_degrees[0] < math.degrees(gamma) < gamma_degrees[1]
    assert abs(value * central_scale - observed_center) < 1.e-14
    assert obs[0] < fixed["65.0"] < obs[1]
    assert mirror[0] > obs[1]
    signed_lower = k[0]**2 / (1-k[0]**2) * scale[0]
    assert signed_lower > obs[1]
    rng = np.random.default_rng(1351)
    err_ratio, err_phase, err_unitary = 0., 0., 0.
    for n in range(256):
        offset = 180 * (n % 2)
        k_i, w_i, rho_i = rng.uniform(*k), rng.uniform(.00001, w_max), rng.uniform(*rho)
        g = math.radians(rng.uniform(*gamma_degrees) + offset)
        v, m = ratio(k_i, w_i, rho_i, g)
        ck = -m[0,0]*m[0,2].conj() / (m[1,0]*m[1,2].conj())
        err_phase = max(err_phase, abs(cmath.phase(ck * cmath.exp(-1j*g))))
        err_ratio = max(err_ratio, abs(v - abs(m[2,0]/m[2,1])**2))
        err_unitary = max(err_unitary, float(np.linalg.norm(m@m.conj().T - np.eye(3))))
        band = mirror if offset else allowed
        assert band[0] <= v*scale[0] <= v*scale[1] <= band[1]
    assert max(err_ratio, err_phase, err_unitary) < 1.e-12
    return {
        "round": 1351,
        "conditional_deterministic_boxes": {
            "Vus": k, "Vcb_upper": w_max, "Vub_over_Vcb": rho, "xi": xi,
            "gamma_degrees_first_branch": gamma_degrees, "second_offset_degrees": 180,
            "per_flavour_residual_max": .01, "mass_ratio": mass, "h_upper": h_max,
            "abs_delta_minus_gamma_degrees_upper": pad, "scale_with_residuals": scale,
        },
        "Bd_over_Bs": {
            "tree_branch_enclosure": allowed, "mirror_enclosure": mirror,
            "cos_delta_nonpositive_simple_lower": signed_lower,
            "measured_ratio_enclosure": obs, "measured_central": observed_center,
            "central_inputs_by_gamma_degrees": fixed,
            "simultaneous_witness_all_other_central_gamma_deg": math.degrees(gamma),
            "witness_ratio": value * central_scale,
            "mirror_extra_abs_Cd_over_Cs_necessary_enclosure": quotient(obs, mirror),
            "mirror_central_248p7deg_exact_recovery_factor": observed_center / fixed["248.7"],
        },
        "seeded_256_checks": {
            "max_CKM_unitarity_norm": err_unitary, "max_gamma_wrapped_error": err_phase,
            "max_ratio_error": err_ratio,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = certificate()
    if args.write:
        with Path(__file__).with_name("results.json").open("x", encoding="utf-8") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)
            f.write("\n")
    print(json.dumps(result, indent=2, ensure_ascii=False))
