"""Two-site calibration of a stable common qubit drift using independent probes."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


DIRECTORY = Path(__file__).resolve().parent
VISIBILITY = math.sqrt(3) / 2
ALPHA = math.atan2(4, 3)


def angle_distance(angle):
    return abs(math.remainder(angle, math.tau))


def nearest_correction(phase, budget):
    return min(range(-budget, budget + 1), key=lambda k: (angle_distance(phase + k * ALPHA), abs(k), k))


def drift(vector, angle):
    # Last axis is the transported B qubit, the first axis is any partner.
    return vector * np.array([1, np.exp(1j * angle)])


def density(vector):
    flat = vector.ravel()
    return np.outer(flat, flat.conj())


def trace_distance(first, second):
    return float(np.abs(np.linalg.eigvalsh(first - second)).sum() / 2)


def simulate_local_readouts(rng, n, phase):
    # Born law: P(a,b) = [1 + a*b*s*cos(phase)]/4 for a,b in {-1,1}.
    alice = 2 * rng.integers(0, 2, size=n) - 1
    mean_product = VISIBILITY * math.cos(phase)
    parity = np.where(rng.random(n) < (1 + mean_product) / 2, 1, -1)
    bob = alice * parity
    return alice, bob


def choose_word(local, received, budget):
    means = [float(np.mean(a * b)) for a, b in zip(local, received, strict=True)]
    x = means[0] / VISIBILITY
    y = (3 * means[0] - 5 * means[1]) / (4 * VISIBILITY)
    return nearest_correction(math.atan2(y, x), budget)


def radius_certificate(budget):
    orbit = json.loads((DIRECTORY / "relational_gate_orbit_results.json").read_text(encoding="utf-8"))
    upper = orbit["results"]["rational_phase"][str(budget)]["covering_radius_bounds_times_1e12"][1]
    radius = 2 * math.asin(upper / (10**12 * VISIBILITY))
    words = sorted((k * ALPHA) % math.tau for k in range(-budget, budget + 1))
    gaps = np.diff(words + [words[0] + math.tau])
    assert float(max(gaps)) / 2 <= radius
    return radius


def minimum_probe_count(separation, failure):
    fidelity = 1 - (VISIBILITY * math.sin(separation / 2))**2
    count = math.ceil(math.log(4 * failure * (1 - failure)) / math.log(fidelity))
    def error(n):
        return (1 - math.sqrt(1 - fidelity**n)) / 2

    assert error(count) <= failure < error(count - 1)
    return count


def run():
    n, budget, failure = 100000, 256, 0.01
    t = math.sqrt(2 * math.log(4 / failure) / n)
    angular_stat_error = math.asin(math.sqrt(5) * t / VISIBILITY)
    angular_grid_error = radius_certificate(budget)
    angular_bound = angular_stat_error + angular_grid_error
    data_bound = math.sin(angular_bound / 2)
    rng = np.random.default_rng(71923)
    reference = np.diag([math.sqrt(3) / 2, 0.5])
    private_data = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
    private_data /= np.linalg.norm(private_data)
    x_basis = np.array([[1, 1], [1, -1]]) / math.sqrt(2)
    examples = []
    for theta in (0., 0.025, 0.4, 1.75, -2.9, math.pi):
        records = []
        for offset in (0., ALPHA):
            phase = theta + offset
            readout_prob = np.abs(x_basis @ drift(reference, phase) @ x_basis.T)**2
            signs = np.array([1, -1])
            expected = (1 + np.outer(signs, signs) * VISIBILITY * math.cos(phase)) / 4
            assert np.allclose(readout_prob, expected)
            records.append(simulate_local_readouts(rng, n, phase))
        alice, bob = zip(*records)
        word = choose_word(alice, bob, budget)
        assert word == choose_word(bob, alice, budget)
        drifted = drift(private_data, theta)
        corrected = drift(drifted, word * ALPHA)
        restored = drift(corrected, -word * ALPHA)
        residual = angle_distance(theta + word * ALPHA)
        measured = trace_distance(density(corrected), density(private_data))
        assert np.allclose(restored, drifted)
        assert measured <= math.sin(residual / 2) + 1e-12
        assert residual <= angular_bound  # deterministic regression on these seeded trials
        examples.append({"drift": round(theta, 8), "agreed_word": word,
                         "residual_radians": round(residual, 8),
                         "trace_distance_with_partner": round(measured, 8)})

    opposite = drift(reference, math.pi)
    assert np.allclose(drift(opposite, math.pi), reference)
    mismatch_floor = trace_distance(density(reference), density(opposite)) / 2
    first_pulse = trace_distance(density(reference), density(drift(reference, ALPHA)))
    assert math.isclose(mismatch_floor, VISIBILITY / 2)
    assert math.isclose(first_pulse, math.sqrt(3 / 20))
    return {"per_setting_probes": n, "total_probes": 2*n, "local_bits_exchanged": 4*n,
            "budget_K": budget, "failure_bound": failure,
            "angle_statistical_bound": round(angular_stat_error, 10),
            "angle_grid_bound": round(angular_grid_error, 10),
            "all_partner_trace_bound": round(data_bound, 10),
            "minimum_probes_for_two_angles_0_and_0.06_at_99pct": minimum_probe_count(0.06, failure),
            "unshared_drift_error_floor": round(mismatch_floor, 10),
            "first_pulse_distance": round(first_pulse, 10), "trials": examples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    result = run()
    output = Path(__file__).with_name("shared_drift_calibration_results.json")
    if parser.parse_args().check:
        assert json.loads(output.read_text(encoding="utf-8")) == result
    else:
        with output.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
