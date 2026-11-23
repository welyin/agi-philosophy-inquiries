"""Round 1375. Canonical kinetic terms and charge-preserving scalar vacua."""
import argparse
import json
from pathlib import Path

import numpy as np

DIRECTORY = Path(__file__).resolve().parent


def rho_and_weights(isospin, hypercharge):
    charged = 2 * (isospin * (isospin + 1) - hypercharge**2)
    neutral = 4 * hypercharge**2
    return charged, neutral


def calculate(data):
    y = data["synthetic"]["y"]
    v = data["synthetic"]["v_in_chosen_units"]
    x = np.sqrt(1 - y**2)
    a = np.array([x, y / 4, np.sqrt(15) * y / 4])
    direction = np.array([x, 4 * y, 0])
    # Complete the light pole to a real orthogonal mixing basis.
    second = np.array([a[1], -a[0], 0])
    second /= np.linalg.norm(second)
    poles = np.stack([a, second, np.cross(a, second)])
    eigs = data["synthetic"]["positive_hessian_eigenvalues"]
    hessian = poles.T @ np.diag(eigs) @ poles
    v_doublet, v_septet = x * v, y * v / 4
    jac = np.diag([v_doublet, v_septet, 1])
    # On real neutral displacements, the sextic polynomial follows v_D**6.
    lock_grad = np.array([-6 * v_septet / v_doublet, 1, 0])
    epsilon = data["synthetic"]["alignment_penalty"]
    remainder = hessian - epsilon * np.outer(lock_grad, lock_grad)
    coefficients = np.linalg.solve(jac, remainder) @ np.linalg.inv(jac)
    assert np.linalg.eigvalsh(coefficients)[0] > 0
    rebuilt = jac.T @ coefficients @ jac + epsilon * np.outer(lock_grad, lock_grad)
    assert np.allclose(rebuilt, hessian) and np.allclose(hessian @ a, eigs[0] * a)

    rep_weights = [rho_and_weights(t, h) for _, t, h in data["representations"]]
    assert rep_weights == [(1, 1), (16, 16), (0, 0)]
    linear = poles @ direction
    kf, kv, k2v = a[0] / x, linear[0], a[0]**2 + 16 * a[1]**2
    assert np.allclose([kf, kv, k2v], 1)
    assert np.isclose(linear[1:] @ linear[1:], 15 * y**2)
    two_field_vals = [x**2 + sign * 4 * y**2 for sign in (1, -1)]
    assert all(not np.isclose(z, 1) for z in two_field_vals)
    n, m = 1, 0
    pell = []
    for _ in range(data["pell_witness_count"]):
        n, m = 2 * n + 3 * m, n + 2 * m
        t, h = (n - 1) / 2, m / 2
        w, z = rho_and_weights(t, h)
        assert n*n - 3*m*m == 1 and w == z > 0
        assert h <= t and (t - h).is_integer()
        pell.append([t, h])
    return {
        "rho_one_complex_representations_T_absY": pell,
        "W_Z_mass_weights": rep_weights,
        "doublet_plus_septet_kV_at_kf_one": two_field_vals,
        "h_directions_doublet_septet_singlet": a.tolist(),
        "kappa_f_V_2V_contact": [float(kf), float(kv), float(k2v)],
        "other_neutral_poles_sum_kV_squared": float(linear[1:] @ linear[1:]),
        "lower_bound_on_larger_other_kV_squared": 15 * y**2 / 2,
        "constructed_hessian_eigenvalues": np.linalg.eigvalsh(rebuilt).tolist(),
        "positive_radial_coefficients": bool(np.all(np.linalg.eigvalsh(coefficients) > 0)),
        "kinetic_counterexample_scope": "synthetic tree-level couplings; extra pole masses and visibility not inferred"
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate(json.loads((DIRECTORY / "inputs.json").read_text(encoding="utf-8")))
    path = DIRECTORY / "results.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    else:
        saved = json.loads(path.read_text(encoding="utf-8"))
        for k, value in result.items():
            if isinstance(value, (list, float, int)):
                np.testing.assert_allclose(value, saved[k], rtol=1e-10, atol=1e-11)
            else:
                assert value == saved[k]
        print("1375: ordinary Yukawa, kinetic weights, equal light couplings and extra pole sum verified")


if __name__ == "__main__":
    main()
