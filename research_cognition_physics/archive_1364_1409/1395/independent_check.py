"""Independent upwards one-loop integration and boundary-sensitivity test."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REFERENCE = HERE.parent / "1369" / "inputs.json"
GAUGE_BETA = np.array([4.1, -19 / 6, -7.0])
GAUGE_CHARGE = np.array([[17 / 20, 9 / 4, 8], [1 / 4, 9 / 4, 8], [9 / 4, 9 / 4, 0]])
YUKAWA_WEIGHT = np.array([[4.5, 1.5, 1], [1.5, 4.5, 1], [3.0, 3.0, 2.5]])
LOOP = 16 * math.pi**2


def evolve(log_span, couplings, initial, steps, top_only=False):
    """Integrate log Yukawas and their forward log-to-log Jacobian."""
    n = len(initial)
    weight = YUKAWA_WEIGHT[:n, :n]
    charge = GAUGE_CHARGE[:n, :]
    assert (n == 1) == top_only
    start = np.r_[np.log(initial), np.eye(n).ravel()]
    step = log_span / steps

    def rhs(t, values):
        y2 = np.exp(2 * values[:n])
        g2 = couplings**2 / (1 - 2 * GAUGE_BETA * couplings**2 * t / LOOP)
        dlog = (weight @ y2 - charge @ g2) / LOOP
        jac = values[n:].reshape((n, n))
        jac_dot = ((2 * weight * y2) @ jac) / LOOP
        return np.r_[dlog, jac_dot.ravel()]

    state = start
    for j in range(steps):
        t = j * step
        k1 = rhs(t, state)
        k2 = rhs(t + step / 2, state + step * k1 / 2)
        k3 = rhs(t + step / 2, state + step * k2 / 2)
        k4 = rhs(t + step, state + step * k3)
        state = state + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return np.exp(state[:n]), state[n:].reshape((n, n))


def results():
    data = json.loads(REFERENCE.read_text(encoding="utf-8"))
    g, dg = np.array([data[name] for name in ("g1_GUT", "g2", "g3")]).T
    y, dy = np.array([data[name] for name in ("yt", "yb", "ytau")]).T
    records = {}
    for high in (10**4, 10**10, 10**16):
        span = math.log(high / data["mu_ref_GeV"])
        uv, jac = evolve(span, g, y[:1], 512, top_only=True)
        slope = 1 / jac[0, 0]
        # For the exact affine inverse-square flow S=1-y_IR^2/y_star^2.
        ceiling = y[0] / math.sqrt(1 - slope)
        full_uv, full_jac = evolve(span, g, y, 512)
        uv_half, jac_half = evolve(span, g, y, 256)
        assert max(np.max(abs(uv_half - full_uv)), np.max(abs(jac_half - full_jac))) < 3e-9
        low_g = g - 2 * dg
        high_g = g + 2 * dg
        _, lower_jac = evolve(span, low_g, [y[0] + 2 * dy[0]], 512, top_only=True)
        _, upper_jac = evolve(span, high_g, [y[0] - 2 * dy[0]], 512, top_only=True)
        lower, upper = 1 / lower_jac[0, 0], 1 / upper_jac[0, 0]
        assert .01 < lower <= slope <= upper < 1
        if high == 10**16:
            # Perturb each infrared coordinate in turn and check the UV Jacobian.
            bump = 1e-4
            for k in range(3):
                direction = np.zeros(3)
                direction[k] = bump
                right, _ = evolve(span, g, y * np.exp(direction), 256)
                left, _ = evolve(span, g, y * np.exp(-direction), 256)
                numerical = (np.log(right) - np.log(left)) / (2 * bump)
                assert np.allclose(numerical, full_jac[:, k], atol=4e-8, rtol=3e-7)
        records[str(high)] = {
            "top_only_y_UV": float(uv[0]),
            "central_log_sensitivity_ir_to_uv": slope,
            "top_only_sensitivity_exact_extrema_in_two_width_input_box": [lower, upper],
            "central_top_only_large_boundary_ceiling": ceiling,
            "central_y_IR_needed_for_at_most_one_percent_sensitivity": math.sqrt(.99) * ceiling,
            "aligned_three_Yukawa_uv": full_uv.tolist(),
            "central_log_sensitivity_fixed_other_uv_yukawas": np.linalg.inv(full_jac)[0, 0],
            "central_log_sensitivity_fixed_other_ir_yukawas": 1 / full_jac[0, 0]
        }
    return records


def check(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            check(actual[key], saved[key])
    else:
        assert np.allclose(actual, saved, atol=5e-11, rtol=2e-9)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    record = results()
    path = HERE / "independent_results.json"
    if args.write:
        with path.open("x", encoding="utf-8") as f:
            json.dump(record, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        check(record, json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(record, allow_nan=False))
