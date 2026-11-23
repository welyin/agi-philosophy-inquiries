"""Recompute the joint projection and the heavy-scale decoupling example."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(relative):
    return json.loads((HERE / relative).read_text(encoding="utf-8"))


def quartet(left):
    return float((left[0, 0] * left[1, 1] * left[0, 1].conj()
                  * left[1, 0].conj()).imag)


def calculate():
    inputs = read("inputs.json")
    base = read(inputs["texture"])
    rare = read(inputs["neutral_rate"])
    old = read(inputs["old_witness"])
    m = np.array(base["ordinary_mass"])
    row = np.array(base["mixing_row_real"]) + 1j * np.array(base["mixing_row_imag"])
    mass = base["heavy_mass"]
    norm2 = float(np.vdot(row, row).real + mass**2)
    z = m @ row.conj()
    scaled = norm2 * (m @ m.T) - np.outer(z, z.conj())
    odd_numerator = (scaled[0, 1] * scaled[1, 2] * scaled[2, 0]).imag
    assert norm2.is_integer() and odd_numerator == round(odd_numerator) != 0
    eff = scaled / norm2
    light_m2, w = np.linalg.eigh(eff)
    slopes = row @ m.T @ w / norm2
    p_slope = float(abs(slopes[0] * slopes[1]))
    j_limit = quartet(w)
    assert min(np.diff(light_m2)) > 0 and p_slope > 0 and abs(j_limit) > 1e-6
    output = []
    for t in inputs["heavy_scales"]:
        matrix = np.block([[m, np.zeros((3, 1))],
                           [t * row[None, :], np.array([[t * mass]])]])
        left, singular, _ = np.linalg.svd(matrix)
        left = left[:, np.argsort(singular)]
        p = np.outer(left[-1].conj(), left[-1])
        assert np.allclose(p, np.eye(4) - left[:3].conj().T @ left[:3])
        assert np.allclose(p @ p, p)
        p_abs = float(abs(p[0, 1]))
        output.append({"t": t, "abs_P_ds": p_abs, "t_squared_abs_P_ds": t*t*p_abs,
                       "weak_quartet": quartet(left), "tree_Z_shape": p_abs**2,
                       "tree_axion_shape": p_abs**2 / t**2})
    assert abs(output[-1]["t_squared_abs_P_ds"] / p_slope - 1) < 1e-4
    assert abs(output[-1]["weak_quartet"] / j_limit - 1) < 1e-4
    bound = inputs["massless_reinterpretation"]["FV_GeV_lower"]
    p_nu = rare["projection_abs_upper"]
    p_old = old["branches"]["proton_null"]["ds_projector"]
    return {"rare_projection_cap": p_nu, "crossing_fa_GeV": bound * p_nu / 2,
            "two_body_cap_at_example": 2 * inputs["example_fa_GeV"] / bound,
            "conditional_two_body_fa_floor_for_old_projection_GeV": bound * p_old / 2,
            "norm_squared": norm2, "odd_numerator": int(odd_numerator),
            "odd_denominator": int(norm2)**3, "light_m2_limit": light_m2.tolist(),
            "limit_t_squared_abs_P": p_slope, "limit_weak_quartet": j_limit,
            "samples": output}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    answer = calculate()
    if args.check:
        def compare(a, b):
            if isinstance(a, dict):
                assert a.keys() == b.keys()
                for key in a:
                    compare(a[key], b[key])
            elif isinstance(a, list):
                assert len(a) == len(b)
                for left, right in zip(a, b):
                    compare(left, right)
            else:
                assert np.isclose(a, b, rtol=2e-9, atol=1e-22), (a, b)
        compare(answer, read("results.json"))
        print("1393 verified")
    else:
        print(json.dumps(answer, indent=2, allow_nan=False))
