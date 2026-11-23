"""Test tree-level mass, scalar, gauge, and Yukawa tensor relations in the same basis."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent


def tensor_test(rng, v):
    # Singlets use real canonical fields. Each orthogonal neutral doublet direction
    # has a gauge-invariant completion by Re or Im of H_1^dagger H_i.
    nd, ns = 3, 2
    dimension = nd + ns
    mass = np.sort(rng.uniform(15, 400, dimension))
    rotation = np.linalg.qr(rng.normal(size=(dimension, dimension)))[0]
    hessian = (rotation * mass**2) @ rotation.T
    jac = np.array([v] * nd + [1.] * ns)
    potential = hessian / np.outer(jac, jac)
    r = Polynomial([0, 1])
    invariant = (r*r-v*v)/2
    # Differentiate the stable polynomial (q^T C q)/2 along the vacuum ray.
    # Here q=(X,r D_2,...,S_1,...); q on the ray is (X,0,...).
    gradients = [potential[j, 0] * (r if j < nd else 1.) * invariant
                 for j in range(dimension)]
    cubic = rotation.T @ np.array([p.deriv(2)(v) for p in gradients])
    quartic = rotation.T @ np.array([p.deriv(3)(v) for p in gradients])
    kappa = rotation[0]
    contact = rotation[:nd].T @ rotation[:nd]
    b = mass**2 * kappa
    dimensionless_error = max(abs(v*cubic - b - 2*contact@b).max(),
                              abs(v*v*quartic - 3*contact@b).max()) / max(mass)**2
    # Two arbitrary complex Yukawa matrix entries, including an off-diagonal one.
    direct = rng.normal(size=(2, nd)) + 1j*rng.normal(size=(2, nd))
    y = direct @ rotation[:nd]
    np.testing.assert_allclose(y @ (v*cubic), 3*y@b, rtol=1e-10)
    np.testing.assert_allclose(y @ kappa, direct[:, 0], atol=2e-14)
    assert dimensionless_error < 1e-11
    return dimensionless_error


def interval(mass2, light_max2, p):
    # Extremize a projector bilinear over the spectrum and singlet orientation.
    c = math.sqrt(1-p)
    assert 0 < p < 1 and light_max2 < min(mass2, 4*c*c*mass2)
    b = c*(2*mass2-light_max2)
    y = light_max2**2 / (2*(math.sqrt(b*b + p*light_max2**2) + b))
    mu = light_max2/2 + c*y
    assert 0 <= mu <= light_max2
    return [1-2*p/3, 1+p*y/(3*c*mass2)], mu


def check_spectral_extrema(mass2, light_max2, p):
    bounds, mu = interval(mass2, light_max2, p)
    recorded = []
    for mean in (0., mu):
        # Light weight at each endpoint has that mean and maximal variance.
        weights = np.array([1-p, p*(1-mean/light_max2), p*mean/light_max2])
        n = np.sqrt(weights)
        e = np.array([1., 0., 0.])
        b = n*np.array([mass2, 0., light_max2])
        hp, bp = e-n*(e@n), b-n*(b@n)
        spectrum = np.linalg.eigvalsh((np.outer(hp, bp)+np.outer(bp, hp))/2)
        recorded.append([1-2*spectrum[-1]/(3*n[0]*mass2),
                         1-2*spectrum[0]/(3*n[0]*mass2)])
    np.testing.assert_allclose([recorded[0][0], recorded[1][1]], bounds, atol=2e-14)
    return {"strength_cap": p, "projected_cubic_interval": bounds,
            "upper_extremizing_light_mean_mass_squared_GeV2": mu}


def run():
    data = json.loads((HERE/"inputs.json").read_bytes())
    recoil = json.loads((HERE/data["OPAL_table_projection"]).read_bytes())["recoil"]
    window = next(row for row in recoil if row["upper_edge_GeV"] == data["light_upper_GeV"])
    m = data["h_mass_GeV"]
    rng = np.random.default_rng(data["seed"])
    errors = [tensor_test(rng, data["v_GeV"]) for _ in range(128)]
    assert max(errors) < 1e-11
    cases = {}
    for key in ["two_width_upper", "all_stat_corr_two_width_upper"]:
        cases[key] = check_spectral_extrema(m*m, data["light_upper_GeV"]**2, window[key])
    ex = data["examples"]
    def two_pole_ratio(weight, other_mass):
        return 1 + 2*weight*((other_mass/m)**2-1)/3
    return {"round": 1790, "rotated_potentials_and_complex_yukawa_tests": len(errors),
            "OPAL_response_conditioned": cases, "pure_doublet_projected_ratio": 1.,
            "light_singlet_mixing_example": two_pole_ratio(ex["second_weight"], ex["second_mass_GeV"]),
            "unseen_heavy_singlet_mixing_example": two_pole_ratio(ex["heavy_weight"], ex["heavy_mass_GeV"])}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    body = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    target = HERE/"results.json"
    if args.check:
        assert json.loads(body) == json.loads(target.read_bytes())
    else:
        target.write_text(body, encoding="utf8")
    print("1790", hashlib.sha256(body.encode()).hexdigest()[:12])
