"""Finite strong-angle response to the gauge-allowed reverse doublet column."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def load(path):
    return json.loads((ROOT / path).read_bytes())


def decode(data):
    array = np.array(data)
    return array[..., 0] + 1j * array[..., 1]


def saved_matrix_case(witness, budget, vev):
    m = np.array(witness["physical_masses_GeV"])
    baseline = decode(witness["weak_basis_mass_GeV"])
    # Enforce the analytic zeros in the rounded published matrix before taking ratios.
    assert max(np.max(abs(baseline[:4, 4])), abs(baseline[4, 3]),
               np.max(abs(baseline[3, :3]))) < 1e-8
    mass = baseline.copy()
    mass[:4, 4] = 0
    mass[4, 3] = 0
    mass[3, :3] = 0
    r = decode(witness["right_D"])
    s = decode(witness["left_S"])
    ell = decode(witness["left_D"])
    assert max(abs(s[:2])) < 1e-12
    kernel = np.linalg.inv(mass)[4, 3]
    spectral = np.sum(r.conj() * s / m)
    heavy_factor = np.prod(1 - (m[2] / m[3:])**2)
    via_pole = r[2].conjugate() * s[2] / m[2] * heavy_factor
    np.testing.assert_allclose([kernel, spectral], via_pole, rtol=2e-9, atol=2e-13)
    assert abs(np.sum(ell.conj() * s)) < 1e-10
    assert abs(np.sum((m/m[-1])**2 * ell.conj() * s)) < 1e-10
    disk = math.sin(budget) / abs(kernel)
    angle_errors = []
    for target in (0.1j, -1.3 + 0.6j, -0.96, 1j*budget):
        varied = mass.copy()
        varied[3, 4] = target / kernel
        ratio = np.linalg.det(varied) / np.linalg.det(mass)
        np.testing.assert_allclose(ratio, 1+target, rtol=2e-12, atol=3e-13)
        angle_errors.append(abs(np.angle(ratio / (1+target))))
    assert max(angle_errors) < 3e-12
    return {"physical_masses_GeV": m.tolist(), "abs_r_b": abs(r[2]), "abs_s_b": abs(s[2]),
            "heavy_factor": heavy_factor, "abs_inverse_entry_per_GeV": abs(kernel),
            "phase_uniform_radius_GeV": disk, "phase_uniform_radius_eV": disk*1e9,
            "phase_uniform_yukawa_radius": math.sqrt(2)*disk/vev,
            "response_to_imaginary_w_for_abs_E_1eV_rad": math.atan(abs(kernel)*1e-9)}


def random_checks(cfg):
    rng = np.random.default_rng(cfg["random_seed"])
    def rand(shape):
        return rng.normal(size=shape) + 1j*rng.normal(size=shape)
    case_count = 0
    for n in (1, 2, 3):
        for _ in range(cfg["random_trials_per_singlet_count"]):
            a = 4*np.eye(3) + rand((3, 3))
            b = rand((3, n))
            c = rand((1, 3))
            d = rand((1, n))
            s = 4*np.eye(n) + rand((n, n))
            f = 2.7 + 0.6j
            e = rand((n, 1))
            def matrix(t, reverse):
                return np.block([[t*a, t*b, np.zeros((3, 1))],
                                 [np.zeros((n, 3)), s, t*reverse],
                                 [t*c, t*d, np.array([[f]])]])
            zero = np.zeros((n, 1))
            m0 = matrix(1, zero)
            k = -(d - c @ np.linalg.solve(a, b)) @ np.linalg.inv(s) / f
            np.testing.assert_allclose(k[0], np.linalg.inv(m0)[-1, 3:3+n], rtol=1e-11, atol=1e-12)
            w = (k@e).item()
            for t in (0.5, 0.9, 1.1):
                ratio = np.linalg.det(matrix(t,e)) / np.linalg.det(matrix(t,zero))
                np.testing.assert_allclose(ratio, 1+t*t*w, rtol=1e-11, atol=1e-12)
            step = 1e-5
            derivative = np.angle(np.linalg.det(matrix(np.exp(step),e)) /
                                  np.linalg.det(matrix(np.exp(-step),e))) / (2*step)
            np.testing.assert_allclose(derivative, 2*w.imag / abs(1+w)**2, atol=2e-8, rtol=5e-8)
            # Exhibit the tangent giving the maximum absolute angle for a disk radius <1.
            for rho in (r for r in cfg["radius_fractions"] if r < 1):
                target = -rho**2 + 1j*rho*math.sqrt(1-rho**2)
                reverse = k.conj().T * target / np.vdot(k,k).real
                np.testing.assert_allclose(np.linalg.norm(reverse)*np.linalg.norm(k), rho)
                ratio = np.linalg.det(matrix(1,reverse))/np.linalg.det(m0)
                np.testing.assert_allclose(np.angle(ratio), math.asin(rho), atol=1e-12)
            # A disk extending beyond the origin contains a negative real determinant ratio.
            rho = max(cfg["radius_fractions"])
            reverse = -rho*k.conj().T / np.vdot(k,k).real
            ratio = np.linalg.det(matrix(1,reverse))/np.linalg.det(m0)
            np.testing.assert_allclose(ratio, 1-rho, atol=1e-12)
            for phase in cfg["finite_probe_complex_angles"]:
                rotated = m0.copy()
                rotated[-1,-1] *= np.exp(1j*phase)
                complete = np.exp(1j*phase)*np.linalg.det(rotated)/np.linalg.det(m0)
                np.testing.assert_allclose(complete, np.exp(2j*phase), atol=1e-12)
            case_count += 1
    return case_count


def normalized(value):
    if isinstance(value, dict):
        return {key: normalized(item) for key, item in value.items()}
    if isinstance(value, list):
        return [normalized(item) for item in value]
    if isinstance(value, float):
        return float(f"{value:.11g}")
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    cfg = load("inputs.json")
    witnesses = load(cfg["pole_matrices"])
    budget = load(cfg["neutron_dictionary"])["theta_max_under_neutron_dictionary"]
    vev = load(cfg["mass_inputs"])["witness"]["Higgs_vev_GeV"]
    out = {"round": 1637, "conditional_angle_budget_rad": budget,
           "random_complex_block_cases": random_checks(cfg)}
    for name in ("common_witness","high_mass_witness"):
        out[name] = saved_matrix_case(witnesses[name], budget, vev)
    result = normalized(out)
    if args.check:
        assert result == load("results.json")
    else:
        (ROOT/"results.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    for name in ("common_witness", "high_mass_witness"):
        print(name, result[name]["abs_inverse_entry_per_GeV"], result[name]["phase_uniform_radius_eV"])


if __name__ == "__main__":
    main()
