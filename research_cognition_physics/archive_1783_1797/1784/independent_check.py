"""Uniform bounds for the rank-changing SM dimension-five running below matching."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
TWO_LOOP = (16 * math.pi**2)**2


def prior_inputs(round_number):
    return json.loads((HERE.parent / str(round_number) / "inputs.json").read_bytes())


def oscillation_floors():
    juno = prior_inputs(1782)["juno"]
    widths = juno["marginal_widths"]
    solar, err = juno["d21_eV2"]
    d_solar = (widths + juno["extra_solar_width_for_ordering_and_external_constraint"]) * err
    normal, *normal_errors = juno["NO_positive_d31_eV2"]
    inverted, *inverted_errors = juno["IO_absolute_d31_eV2"]
    # Use the larger asymmetric error in either direction to give conservative windows.
    max_atm = normal + widths * max(normal_errors)
    min_atm = inverted - widths * max(inverted_errors)
    return {"NO": math.sqrt((solar - d_solar) / max_atm),
            "IO": math.sqrt(min_atm / (min_atm + solar + d_solar))}


def random_checks():
    # Deliberately large exposures make the limiting mass ratios resolvable in float64.
    rng = np.random.default_rng(1784)
    for _ in range(256):
        p = rng.uniform(0, 1, (5, 3))
        radius = np.ptp(p, axis=1) / 2
        centered = p - (p.max(axis=1) - radius)[:, None]
        gram = centered.T @ centered / 5
        epsilon = np.mean(radius**2)
        kernel = np.exp(-gram)
        vectors = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
        for rank in (1, 2):
            columns = vectors[:, :rank]
            before = columns @ columns.T
            after = before * kernel
            _, masses, right = np.linalg.svd(after)
            eta = np.expm1(epsilon)
            assert np.linalg.norm(after - before, 2) <= eta * masses[0] + 1e-12
            if rank == 1:
                assert masses[1] <= np.tanh(epsilon) * masses[0] + 1e-12
            else:
                assert masses[2] <= eta * masses[0] + 1e-12
                ratio = masses[1] / masses[0]
                if ratio > eta:
                    null = np.linalg.svd(before)[2][-1].conj()
                    observed = right[-1].conj()
                    sin_angle = np.linalg.norm(observed - null * np.vdot(null, observed))
                    assert sin_angle <= 2 * eta / (ratio - eta) + 1e-12
    for epsilon in (0.001, 0.3, 1.4):
        block = np.exp(-epsilon * np.array([[1, -1], [-1, 1]])) / 2
        masses = np.linalg.svd(block, compute_uv=False)
        assert abs(masses[1] / masses[0] - math.tanh(epsilon)) < 1e-14
    return {"complex_trials_per_rank": 256, "sharp_extremizers": 3}


def run():
    weak_input = prior_inputs(1369)
    weak = weak_input["mu_ref_GeV"]
    floors = oscillation_floors()
    windows = []
    for scale in (1e14, 1e16, 1e19):
        length = math.log(scale / weak)
        span_cap = 0.02
        epsilon = span_cap**4 * length / (2 * TWO_LOOP)
        eta = math.expm1(epsilon)
        row = {"scale_GeV": scale, "sqrt_yukawa_squared_span_cap": span_cap,
               "epsilon_cap": epsilon, "rank_one_middle_ratio_cap": math.tanh(epsilon),
               "rank_two_light_ratio_cap": eta}
        for ordering, floor in floors.items():
            row[ordering] = {"observed_to_rank_one_bound": floor / math.tanh(epsilon),
                             "rank_two_null_rotation_sine_cap": 2 * eta / (floor - eta),
                             "required_effective_sqrt_span":
                                 (2 * TWO_LOOP * math.atanh(floor) / length)**0.25}
        windows.append(row)
    return {"weak_reference_GeV": weak, "weak_tau_msbar": weak_input["ytau"],
            "middle_over_heaviest_floors": floors, "windows": windows,
            "verification": random_checks()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    dest = HERE / "independent_results.json"
    if args.check:
        assert json.loads(dest.read_bytes()) == result
    else:
        dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("1784 two-loop rank bounds verified", result["middle_over_heaviest_floors"])
