"""Group-independent three-clock readout of the live two-leg electric energy."""
from pathlib import Path
from itertools import product
import argparse
import json
import math

import numpy as np


def generate():
    durations = [math.pi, 2 * math.pi / 5, math.pi / 3]
    branches = {"SU3": ([1, 2.5], [3, 6]),
                "SO3": ([0, 1, 3], [1, 3, 5]),
                "mu3_SO3_local_clock": ([0, 1, 3], [1, 3, 5])}
    out = {}
    for label, (energies, counts) in branches.items():
        numeric = np.array(energies)
        frequencies = np.array(counts) / 9
        histories = {}
        for word in product((1, -1), repeat=3):
            amplitudes = np.ones(len(energies), complex)
            for sign, tau in zip(word, durations):
                amplitudes *= (1 + sign * np.exp(-1j * numeric * tau)) / 2
            histories["".join("+" if s == 1 else "-" for s in word)] = float(
                np.dot(frequencies, abs(amplitudes) ** 2))
        assert np.isclose(sum(histories.values()), 1)
        target = 0 if label == "SU3" else 1/9
        assert abs(histories["+++"] - target) < 1e-14
        out[label] = {"all_eight_records": histories, "joint_with_path_minus": histories["+++"] * .75}
    error = .02
    return {"clock_times_in_units_hbar_per_kappa": durations, "branches": out,
            "max_phase_error_per_mode_per_clock": error,
            "positive_mode_acceptance_max": math.sin(error / 2)**2,
            "zero_mode_acceptance_min": math.cos(error / 2)**6}


if __name__ == "__main__":
    opts = argparse.ArgumentParser()
    opts.add_argument("--write", action="store_true")
    args = opts.parse_args()
    record = generate()
    path = Path(__file__).with_name("clock_delivery_results.json")
    if args.write:
        with path.open("x", encoding="utf8") as out:
            out.write(json.dumps(record, indent=2) + "\n")
    else:
        saved = json.loads(path.read_text(encoding="utf8"))
        for group, row in record["branches"].items():
            expected = saved["branches"][group]
            for key, val in row["all_eight_records"].items():
                assert abs(val-expected["all_eight_records"][key]) < 1e-12
            assert abs(row["joint_with_path_minus"]-expected["joint_with_path_minus"]) < 1e-12
        for key in record.keys() - {"branches"}:
            np.testing.assert_allclose(saved[key], record[key], atol=1e-13)
    print({name: v["joint_with_path_minus"] for name, v in record["branches"].items()})
