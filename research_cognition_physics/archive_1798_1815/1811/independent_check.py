"""Heavy-light couplings of the radial SU(5) adjoint for two vectorlike fives."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE / relative).read_bytes())


def unitary(n, rng):
    q, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n)))
    return q


def check_geometry(cosines, rng):
    full = unitary(5, rng)
    basis_t, alternate = full[:, :2], full[:, 2:4]
    basis_d = basis_t*np.array(cosines) + alternate*np.sqrt(1-np.square(cosines))
    subspace_t = basis_t @ basis_t.conj().T
    subspace_d = basis_d @ basis_d.conj().T
    mass_t = unitary(2, rng) @ np.diag([9., 3.]) @ unitary(2, rng) @ basis_t.conj().T
    mass_d = unitary(2, rng) @ np.diag([2., .7]) @ unitary(2, rng) @ basis_d.conj().T
    adjoint = (mass_t-mass_d)/5
    bare = (3*mass_t+2*mass_d)/5
    for mass, weight in ((mass_t, 2), (mass_d, -3)):
        np.testing.assert_allclose(bare+weight*adjoint, mass, atol=1e-14)
    y_d = -3*adjoint @ (np.eye(5)-subspace_d)
    y_t = 2*adjoint @ (np.eye(5)-subspace_t)
    def spectrum(matrix):
        return np.linalg.svd(matrix, compute_uv=False)
    # For an exactly overlapping heavy direction use the cross-multiplied equality.
    lhs = np.prod(spectrum(y_d))*np.prod(spectrum(mass_d))
    rhs = 9/4*np.prod(spectrum(y_t))*np.prod(spectrum(mass_t))
    np.testing.assert_allclose(lhs, rhs, atol=8e-13)
    expected = (1-np.square(cosines)).prod()
    np.testing.assert_allclose((np.prod(spectrum(y_d))/
        (.6**2*np.prod(spectrum(mass_t))))**2, expected, atol=3e-13)
    if max(cosines) == 0:
        np.testing.assert_allclose(spectrum(y_d), .6*spectrum(mass_t))
        np.testing.assert_allclose(spectrum(y_t), .4*spectrum(mass_d))
    # A gauge-singlet mass allowed to scale along the same fluctuation kills both vertices.
    for mass, projection in ((mass_t, subspace_t), (mass_d, subspace_d)):
        np.testing.assert_allclose(mass @ (np.eye(5)-projection), 0, atol=1e-13)


def gauge_mass_gains(sample):
    # Central witnesses are pairwise degenerate within each species. Integrate
    # d log m / d log mu = -6 sum C_i g_i^2/(16 pi^2), starting at its threshold.
    q = sample["matching_scale_GeV"]
    doublet = np.array([.6, 1., 0.])
    triplet = 1-doublet
    b = np.array([4.1, -19/6, -7.])+4*(doublet+triplet)/3
    inv = np.array(sample["inverse_alphas_at_matching"])
    casimirs = {"D": np.array([3/20, .75, 0.]), "T": np.array([1/15, 0., 4/3])}
    gains = {}
    area = np.zeros(3)
    for mass, step, label in [(sample["M_T_GeV"][0], triplet, "T"),
                              (sample["M_D_GeV"][0], doublet, "D")]:
        assert mass < q
        dt = math.log(q/mass)/(2*math.pi)
        area += np.log1p(b*dt/inv)/b
        gains[label] = math.exp(-3*(casimirs[label] @ area))
        inv = inv + b*dt
        b -= 4*step/3
        q = mass
    return gains


def run():
    rng = np.random.default_rng(1811)
    cases = [[0., 0.], [.24, .83], [1., .31], [1., 1.]]
    for angles in cases:
        for _ in range(7):
            check_geometry(angles, rng)
    saved = read("../1807/independent_results.json")
    out = {"complex_matrix_orientations_checked": 7*len(cases)}
    for key in ("common_MZ", "independent_lattice"):
        box = saved[key]["two_marginal_box"]
        tmax = box["log_scale"][1]
        kmin = 3*math.pi*box["K"][0]
        kmax = 3*math.pi*box["K"][1]
        # Outward enclosure, extrema in log scale and K need not occur at the same input.
        bounds = [math.exp(kmin)/math.cosh(2*tmax-kmin),
                  math.exp(kmax)*math.cosh(2*tmax-kmin)]
        assert kmin > tmax and kmax < 2*box["log_scale"][0]
        t, k = saved[key]["central"]["log_scale"], 3*math.pi*saved[key]["central"]["K"]
        def ratio(ld, lt):
            assert abs(sum(ld)-sum(lt)-k) < 1e-13
            return sum(math.exp(-2*x) for x in lt)/sum(math.exp(-2*x) for x in ld)
        least = ratio([t, k-t], [0., 0.])
        largest = ratio([t, t], [2*t-k, 0.])
        np.testing.assert_allclose([least, largest],
            [math.exp(k)/math.cosh(2*t-k), math.exp(k)*math.cosh(2*t-k)])
        sample = saved[key]["two_split_complete_pairs"]
        gains = gauge_mass_gains(sample)
        d = np.array(sample["M_D_GeV"])
        trip = np.array(sample["M_T_GeV"])
        inv_mass = sample["matching_scale_GeV"]
        assert min(trip) < inv_mass
        rd = (.6*trip/inv_mass)**2
        rt = (.4*d/inv_mass)**2
        # Set the radial vev and scalar pole to Q for illustration; these are extra inputs.
        # Both species are internally degenerate in these examples.
        factor_d = (1-(d[0]/inv_mass)**2)**2
        factor_t = (1-(trip[0]/inv_mass)**2)**2
        out[key] = {
            "nonzero_portal_singular_product_D_over_T_range":
                [2.25*math.exp(kmin), 2.25*math.exp(kmax)],
            "exact_screening_frobenius_squared_D_over_T_enclosure": [2.25*b for b in bounds],
            "mass_square_ratio_extrema_fixed_central_gauge": [least, largest],
            "gauge_only_MS_mass_gains_from_each_pole_to_matching": gains,
            "singular_product_ratio_multiplier_if_gauge_only_mass_running":
                (gains["T"]/gains["D"])**2,
            "central_example_sigma_vev_and_mass_equal_to_matching_scale": {
                "mass_GeV": inv_mass, "sigma_D_light_portal_squared_sum": float(sum(rd)),
                "sigma_T_light_portal_squared_sum": float(sum(rt)),
                "summed_partial_width_over_mass_D": float(2*sum(rd)*factor_d/(8*math.pi)),
                "summed_partial_width_over_mass_T": float(3*sum(rt)*factor_t/(8*math.pi))}}
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    answer = run()
    # Proven ratios span 25 orders of magnitude. Retain ten significant figures.
    raw = json.dumps(json.loads(json.dumps(answer), parse_float=lambda x: float(f"{float(x):.10g}")),
                     indent=2, allow_nan=False) + "\n"
    target = HERE/"independent_results.json"
    if parser.parse_args().check:
        assert json.loads(raw) == json.loads(target.read_bytes())
    else:
        target.write_text(raw, encoding="utf-8")
    print("1811 independent", hashlib.sha256(raw.encode()).hexdigest()[:12])
