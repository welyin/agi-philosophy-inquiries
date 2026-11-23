"""Physical sensitivity of a shared lepton-leg matching and a sparse two-pole orbit."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def angle_bound(k, theta, charged_ratio_squared):
    mass_part = charged_ratio_squared * (k*k-1) / (2*(1-charged_ratio_squared*k*k))
    return min(math.pi/2, 2*math.atan(k*math.tan(theta/2)) + math.asin(min(1, mass_part)))


def threshold(theta, required_angle, ratio):
    lo, hi = 1.0, 5.0
    assert angle_bound(lo, theta, ratio) < required_angle < angle_bound(hi, theta, ratio)
    for _ in range(60):
        mid = (lo + hi) / 2
        if angle_bound(mid, theta, ratio) < required_angle:
            lo = mid
        else:
            hi = mid
    return lo


def shared_leg_checks(theta, ratio):
    rng = np.random.default_rng(1787)
    e = np.eye(3)[:, 0]
    b = np.array([math.cos(theta), math.sin(theta), 0])
    h = np.diag([ratio, 1, 280.0])
    for _ in range(128):
        matrix = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
        unitary = np.linalg.qr(matrix)[0]
        k = rng.uniform(1.01, 4.9)
        transform = (unitary * [1.0, math.sqrt(k), k]) @ unitary.conj().T
        zero = np.linalg.solve(transform, b)
        zero /= np.linalg.norm(zero)
        _, vectors = np.linalg.eigh(transform.conj().T @ h @ transform)
        angle = math.acos(min(1, abs(np.vdot(zero, vectors[:, 0]))))
        assert angle < angle_bound(k, theta, ratio) + 1e-11
        massless_e = np.linalg.solve(transform, e)
        assert np.linalg.norm((transform.conj().T @ np.diag([0, 1, 280.0]) @ transform) @ massless_e) < 1e-10
    return 128


def run():
    stage = HERE.parent
    prior = json.loads((stage/"1784/results.json").read_bytes())
    oscillation = json.loads((stage/"1782/results.json").read_bytes())["experimental_rectangles"]
    bounds = prior["class_bounds"]
    theta = math.asin(bounds["first_axis_tilt_sine"])
    # Conditional electron/muon mass cap rounded upward for the common endpoint.
    ratio = 0.005**2
    d21_max = oscillation["solar_frequency_eV2"][1]
    juno = json.loads((stage/"1782/inputs.json").read_bytes())["juno"]
    central, *errors = juno["NO_positive_d31_eV2"]
    m3max = math.sqrt(central + 2*max(errors))
    m2max = math.sqrt(d21_max)
    x3max = float(oscillation["isolated_weight_roots"][0][1])
    amplitude = math.sqrt(m3max*x3max)
    sine_needed = (math.sqrt(amplitude**2 + .003) - amplitude) / math.sqrt(m2max+m3max)
    angles = {"inverted_ordering": math.asin(bounds["IO_required_tilt_sine_min"]),
              "three_meV_at_normal_ordering": math.asin(sine_needed)}
    witness = prior["two_source_witness_02"]
    coeffs = np.array(witness["raw_real_D_coefficients_before_common_scale"])
    charges = np.array([2, 0, 1])
    u, v = coeffs[0]*.5**charges, coeffs[1]*.5**(charges+2)
    null_partner = u-v/2
    heavy = np.array([[1., 1.], [1., 0.]])
    m0 = 2e7
    # An illustrative vacuum convention; only Yukawa rows depend on this choice.
    vu_over_sqrt2 = 170.0
    mass_scale = witness["common_mass_scale_eV"] * 1e-9
    physical = math.sqrt(m0*mass_scale)
    lights = np.outer(null_partner, v)+np.outer(v, null_partner)
    eig, directions = np.linalg.eigh(heavy)
    masses = m0*np.abs(eig)
    logs = np.log(masses.max()/masses)
    yukawa_cap = 4*math.pi
    finite_leg_size = yukawa_cap / (32*math.pi**2)
    one_loop_kappa = math.exp(yukawa_cap*max(logs)/(16*math.pi**2)) / (1-finite_leg_size)**2
    examples = {}
    for boost in (.001, 1, 1000):
        dirac = np.column_stack((boost*null_partner+v/(2*boost), v/boost))
        assert dirac[0, 1] == 0.
        recovered = dirac @ np.linalg.solve(heavy, dirac.T)
        assert np.linalg.norm(recovered-lights, 2) < 3e-9
        y = physical/vu_over_sqrt2 * dirac @ directions
        gram = y @ y.T
        log_gram = (y*logs) @ y.T / (16*math.pi**2)
        examples[str(boost)] = {"frobenius_yukawa_squared": float(np.trace(gram)),
                                "largest_eigenvalue_split_log_gram": float(np.linalg.eigvalsh(log_gram)[-1]),
                                "electron_row_yukawa_squared": float(gram[0, 0])}
    return {"conditional_me_over_mmu_cap": math.sqrt(ratio), "baseline_tilt_sine": math.sin(theta),
            "required_tilt_sines": {k: math.sin(v) for k, v in angles.items()},
            "necessary_effective_common_leg_condition_numbers":
                {k: threshold(theta, v, ratio) for k, v in angles.items()},
            "shared_kinetic_matrix_trials": shared_leg_checks(theta, ratio),
            "split_pole_example_truncated_common_leg_kappa_upper": one_loop_kappa,
            "tilt_sine_if_effective_common_leg_has_this_kappa": math.sin(angle_bound(one_loop_kappa, theta, ratio)),
            "illustrative_heavy_masses_GeV": masses.tolist(), "vu_over_sqrt2_GeV": vu_over_sqrt2,
            "fixed_light_matrix_fixed_heavy_pole_orbit": examples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    out = HERE/"independent_results.json"
    if args.check:
        assert json.loads(out.read_bytes()) == result
    else:
        out.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print("1787 common lepton metric", result["necessary_effective_common_leg_condition_numbers"])
