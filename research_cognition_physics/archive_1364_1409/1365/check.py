"""Round 1365. Run with no arguments to verify the saved numerical results."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
HIGGS_MASS = 8.0
VEV = 20.0
UP_MASSES = np.array([0.5, 2.0, 9.0])


def scalar_loop(mass):
    tau = (HIGGS_MASS / (2 * mass)) ** 2
    if tau <= 1:
        f = math.asin(math.sqrt(tau)) ** 2
    else:
        rapidity = 2 * math.acosh(math.sqrt(tau))
        f = -0.25 * (rapidity - 1j * math.pi) ** 2
    return 1.5 * (tau + (tau - 1) * f) / tau**2


def cycle(matrix):
    return (matrix[0, 1] * matrix[1, 2] * matrix[2, 0]).imag


def polar_data(matrix):
    left, masses, right_dagger = np.linalg.svd(matrix)
    order = np.argsort(masses)
    return left[:, order], masses[order], right_dagger.conj().T[:, order]


def widths(masses, kappas, reference):
    assert np.allclose(kappas.imag, 0, atol=1e-10)
    up = sum(scalar_loop(mass) for mass in UP_MASSES)
    down = sum(k.real * scalar_loop(m) for m, k in zip(masses, kappas))
    amplitude = (up + down) / reference
    return {
        "bb_width_over_SM": float(kappas[2].real ** 2),
        "gg_width_over_SM": float(abs(amplitude) ** 2),
        "gg_amplitude_over_SM": [float(amplitude.real), float(amplitude.imag)],
    }


def evaluate():
    matrix = np.array([[1, 0, 0, 0], [-1, 2, 0, 0], [1, -1, 3, 0],
                       [-1+2j, -2+1j, -2+2j, 4]], dtype=complex)
    left, masses, right = polar_data(matrix)
    light_v = left[:3, :3]
    projector = np.diag([1, 1, 1, 0])
    # Light charged currents fix all diagonal scalar responses in this model.
    deficits = 1 - np.sum(abs(light_v)**2, axis=0)
    kappas = np.r_[1 - deficits, sum(deficits)]
    assert np.allclose(kappas, np.diag(left.conj().T @ projector @ left))
    assert np.allclose(left.conj().T @ projector @ matrix @ right,
                       (left.conj().T @ projector @ left) @ np.diag(masses))
    # The real upper block of M M^dagger supplies one simple root.
    h0 = (matrix @ matrix.conj().T)[:3, :3]
    heavy_projector = np.outer(left[:3, 3], left[:3, 3].conj())
    plus = cycle(h0 + heavy_projector)
    minus = cycle(h0 - heavy_projector)
    linear, quadratic = (plus - minus) / 2, (plus + minus) / 2
    heavy_squares = [masses[3]**2, masses[3]**2 - linear / quadratic]
    assert abs(cycle(h0)) < 1e-12 and abs(quadratic) > 1e-7
    assert max(masses[:3]) < HIGGS_MASS / 2 < min(np.sqrt(heavy_squares))
    reference = sum(scalar_loop(m) for m in [*UP_MASSES, *masses[:3]])
    predictions = []
    for square in heavy_squares:
        branch_masses = np.r_[masses[:3], math.sqrt(square)]
        assert abs(cycle(h0 + (square - masses[3]**2) * heavy_projector)) < 1e-11
        predictions.append(widths(branch_masses, kappas, reference))
    assert predictions[0]["bb_width_over_SM"] == predictions[1]["bb_width_over_SM"]
    separation = abs(predictions[0]["gg_width_over_SM"] - predictions[1]["gg_width_over_SM"])
    assert separation > 1e-5
    # Local EFT mass coefficient (H^dagger H - v^2/2) Dbar_L D_R / Lambda.
    # Choose v^2 times its real coefficient to be one mass unit.
    def eft_mass(background):
        value = matrix.copy()
        value[:3] *= background / VEV
        value[3, 3] += (background**2 - VEV**2) / (2 * VEV**2)
        return value

    assert np.array_equal(eft_mass(VEV), matrix)
    step = VEV * 1e-5
    derivative = (eft_mass(VEV+step) - eft_mass(VEV-step)) * VEV / (2*step)
    eft_kappas = np.diag(left.conj().T @ derivative @ right) / masses
    assert np.allclose(eft_kappas, kappas + 4 * abs(left[3])**2 / masses**2)
    eft_rates = widths(masses, eft_kappas, reference)
    assert abs(eft_rates["bb_width_over_SM"] - predictions[0]["bb_width_over_SM"]) > 1e-3
    return {"scope": "synthetic common mass units, tree-level Yukawa and one-loop hgg",
            "scalar_mass": HIGGS_MASS, "vev": VEV, "up_masses": UP_MASSES.tolist(),
            "light_down_masses": masses[:3].tolist(),
            "heavy_mass_squared_branches": [float(x) for x in heavy_squares],
            "light_column_deficits": deficits.tolist(), "diagonal_kappas": kappas.tolist(),
            "predictions": predictions, "absolute_gg_width_ratio_separation": separation,
            "EFT_same_vacuum_first_branch": eft_rates}


def check_saved(actual, saved):
    if isinstance(actual, list) and actual and isinstance(actual[0], dict):
        assert len(actual) == len(saved)
        for a, s in zip(actual, saved):
            check_saved(a, s)
    elif isinstance(actual, dict):
        assert set(actual) == set(saved)
        for key, value in actual.items():
            check_saved(value, saved[key])
    elif isinstance(actual, str):
        assert actual == saved
    else:
        np.testing.assert_allclose(actual, saved, rtol=2e-8, atol=1e-10)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = evaluate()
    path = HERE / "results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as output:
            json.dump(data, output, ensure_ascii=False, indent=2, allow_nan=False)
            output.write("\n")
    else:
        check_saved(data, json.loads(path.read_text(encoding="utf-8")))
    print("Round 1365 verified; normalized gg width separation =", data["absolute_gg_width_ratio_separation"])
