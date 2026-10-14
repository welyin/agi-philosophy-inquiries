"""After-1194 candidate screen; diagnostic only, not round 1195.
Run: python -B -X utf8 check.py [--output path]
No finite-K quantum, instrument, interval or GR certificate is claimed.
"""
from pathlib import Path
import argparse
import json
import math
import numpy as np


def rhs(z, delta):
    sx, sy, sz, tx, ty, tz, tau = z
    b = .15 * (2 - .5 * tx) + delta * tz
    return np.array([
        -2 * b * sy, 2 * b * sx + sz, -sy,
        -2 * delta * sz * ty,
        2 * delta * sz * tx + (1 + .15 * sz) * tz,
        -(1 + .15 * sz) * ty, 1 + .15 * sz,
    ])


def integrate(T, delta, steps):
    z = np.array([0., 0., 1., 1., 0., 0., 0.])
    dt = T / steps
    for _ in range(steps):
        k1 = rhs(z, delta)
        k2 = rhs(z + dt * k1 / 2, delta)
        k3 = rhs(z + dt * k2 / 2, delta)
        k4 = rhs(z + dt * k3, delta)
        z += dt * (k1 + 2*k2 + 2*k3 + k4) / 6
    return {
        "state": z.tolist(),
        "probe_probability": (1 + math.sin(.2 * z[6])) / 2,
        "reaction_Gt_probability": (1 + z[5]) / 2,
        "bloch_norm_errors": [
            abs(float(np.dot(z[:3], z[:3])) - 1),
            abs(float(np.dot(z[3:6], z[3:6])) - 1),
        ],
    }


def main():
    rows = []
    for T in [.5, math.pi/2]:
        branches = {}
        for d in [0., .2]:
            coarse = integrate(T, d, 4000)
            fine = integrate(T, d, 8000)
            fine["step_doubling_state_difference"] = float(
                np.max(np.abs(np.array(coarse["state"]) - fine["state"])))
            assert fine["step_doubling_state_difference"] < 1e-10
            assert max(fine["bloch_norm_errors"]) < 1e-10
            branches[str(d)] = fine
        gap = abs(branches["0.0"]["probe_probability"]
                  - branches["0.2"]["probe_probability"])
        # A screening decision, not a rigorous bound on the physical device.
        rows.append({
            "T": T, "branches": branches,
            "probe_probability_gap_diagnostic": gap,
            "old_two_branch_instrument_budget": 4e-5,
            "diagnostic_gap_below_old_instrument_budget": bool(gap < 4e-5),
            "reaction_Gt_probability_difference":
                branches["0.2"]["reaction_Gt_probability"]
                - branches["0.0"]["reaction_Gt_probability"],
        })
    assert all(r["diagnostic_gap_below_old_instrument_budget"] for r in rows)
    return {
        "schema": "after1194_unaccepted_candidate_screen_v1",
        "status": "DIAGNOSTIC_REPRODUCED_NOT_PHYSICAL_CERTIFICATE",
        "method": "RK4 at 4000 and 8000 steps, NumPy float64",
        "certified_new_research_rounds": 0,
        "new_scientific_check_groups": 0,
        "assumptions": {
            "model": "Hartree/Bloch comparator for proposed reaction code",
            "epsilon": .15, "Omega": .5,
            "d_central": 0, "delta_compared": [0., .2],
            "initial_state": [0, 0, 1, 1, 0, 0, 0],
            "old_probe_readout": "(1+sin(0.2*tau))/2",
        },
        "rows": rows,
        "limits": [
            "Step doubling is a numerical consistency diagnostic, not interval certification.",
            "No full finite-K quantum evolution or actual readout is simulated.",
            "The larger reaction-readout gap does not replace the original probe readout.",
            "No tuning of time, precision or resources to force a passing witness.",
            "Deduplication against 1040 and 1178 is documented in README and scope_review.",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = main()
    text = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(text)
    print(text)
