"""Reproduce the photon/ZZ diagnostic via a scalar triangle integral and a 2D fit."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np


def scalar_triangle(tau):
    # Feynman triangular domain, including the scalar contact graph normalization.
    nodes, weights = np.polynomial.legendre.leggauss(24)
    x = ((nodes + 1) / 2)[:, None]
    y = (1 - x) * ((nodes + 1) / 2)[None, :]
    area = (1 - x) * weights[:, None] * weights[None, :] / 4
    assert 0 < tau < 1
    val = float(np.sum(area * 8 * x*y / (1 - 4*tau*x*y)))
    closed = (math.asin(math.sqrt(tau))**2 - tau) / tau**2
    assert math.isclose(val, closed, rel_tol=3e-12)
    return val


def threshold(tau):
    # acosh branch for masses below half of the Higgs mass.
    if tau <= 1:
        return math.asin(math.sqrt(tau))**2
    return -(math.acosh(math.sqrt(tau)) - .5j*math.pi)**2


def charged_sector():
    # Independent transcription of Antusch, Hinze & Saad, v2, equation (4).
    # Weak scale running Yukawa eigenvalues, GUT-normalized g1.
    v, g1, g2, mh = 248.401, .461228, .65096, 125.38
    mw, mz = v*g2/2, v*math.sqrt(g2*g2 + .6*g1*g1)/2
    yukawas = [
        (.00000704, 3, 2/3), (.00356, 3, 2/3), (.967, 3, 2/3),
        (.0000154, 3, -1/3), (.000306, 3, -1/3), (.0163, 3, -1/3),
        (.00000277713, 1, -1), (.000585042, 1, -1), (.0099378, 1, -1),
    ]
    tw = (mh / (2*mw))**2
    aw = -(2 + 3/tw + 3*(2/tw - 1/tw**2)*threshold(tw))
    af = 0j
    fermion_mass_moment = 0.
    for yukawa, multiplicity, charge in yukawas:
        mass = v*yukawa/math.sqrt(2)
        fermion_mass_moment += 4*multiplicity*mass**4
        t = (mh / (2*mass))**2
        af += multiplicity * charge**2 * 2/t * (1 + (1 - 1/t)*threshold(t))
    scalar_mass = ((8*math.pi**2 * v*v * mh*mh - 6*mw**4 - 3*mz**4
                    + fermion_mass_moment) / 2)**.25
    a0 = scalar_triangle((mh / (2*scalar_mass))**2)
    return scalar_mass, aw, af, a0


def gaussian(q, zz_error):
    # CMS Figure 2b and B.5b. Positive common normalization is free.
    obs = np.array([1.13, .97])
    errors = np.array([.09, zz_error])
    cov = np.outer(errors, errors) * np.array([[1., .20], [.20, 1.]])
    direction = np.array([q, 1.])
    inv = np.linalg.inv(cov)
    norm = max(0., (direction @ inv @ obs) / (direction @ inv @ direction))
    residual = obs - norm*direction
    z = math.sqrt(residual @ inv @ residual)
    difference = obs[0] - q*obs[1]
    assert norm > 0 and np.isclose(z*z, difference**2 / (cov[0, 0] + q*q*cov[1, 1] - 2*q*cov[0, 1]))
    return z


def certificate():
    mass, w, f, a0 = charged_sector()
    base = w + f
    result = {
        "M_complex_scalar_GeV": mass, "A_spin1": w, "A_all_spin_half_real": f.real,
        "A_all_spin_half_imag": f.imag, "A_one_complex_scalar_unit_charge": a0,
        "ratio_CMS_center": 1.13/.97, "scalar_to_heavy_limit": 3*a0,
    }
    for q in (0, 1, 2):
        amp_shift = q*q*a0
        rate = abs(base + amp_shift)**2 / abs(base)**2
        reflected = -2*base.real - amp_shift
        assert np.isclose(abs(base+reflected), abs(base+amp_shift))
        result[f"charge_{q}_rate"] = rate
        result[f"charge_{q}_reflected_shift"] = reflected
        result[f"charge_{q}_positive_common_norm_z"] = gaussian(rate, .115)
        z_end = [gaussian(rate, e) for e in (.11, .12)]
        result[f"charge_{q}_smallest_z_endpoint"] = min(z_end)
        result[f"charge_{q}_largest_z_endpoint"] = max(z_end)
    # A real-loop extraction from a rate admits two branches.
    rho = result["ratio_CMS_center"]
    delta = cmath.sqrt(rho*abs(base)**2 - base.imag**2).real
    result["central_rate_real_shift_minus"] = -base.real - delta
    result["central_rate_real_shift_plus"] = -base.real + delta
    return {key: round(float(val), 10) for key, val in result.items()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    actual = certificate()
    path = Path(__file__).with_name("independent_results.json")
    if args.write:
        with path.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(actual, indent=2) + "\n")
    else:
        stored = json.loads(path.read_text(encoding="utf-8"))
        assert stored.keys() == actual.keys()
        for key, value in stored.items():
            assert math.isclose(value, actual[key], rel_tol=1e-10, abs_tol=2e-9), key
    print(json.dumps(actual))
