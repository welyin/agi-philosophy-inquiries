"""Checked predictions for repeated weak-source and clock records (round 1350)."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np


def power(matrix, exponent):
    return np.linalg.matrix_power(matrix, exponent)


def root(matrix):
    vals, vecs = np.linalg.eigh(matrix)
    assert vals.min() >= -1e-12
    return (vecs * np.sqrt(np.maximum(0, vals))) @ vecs.conj().T


def clock_ports(effect):
    values, vectors = np.linalg.eigh(effect)
    assert values.min() > -1e-12 and values.max() < 1 + 1e-12
    energies = 2 * np.arccos(np.sqrt(np.clip(values, 0, 1)))
    timed = (vectors * np.exp(-1j * energies)) @ vectors.conj().T
    eye = np.eye(len(effect))
    return energies, ((eye + timed) / 2, (eye - timed) / 2)


def direct_record(source, ports, word):
    joint = source.copy()
    for port in word:
        joint = ports[port] @ joint
    return joint.conj().T @ joint


def readings(moment, preparations):
    return {label: float(np.vdot(a, moment @ a).real) for label, a in preparations.items()}


def model(effect, source, preparations):
    energies, ports = clock_ports(effect)
    records = []
    for depth in range(5):
        measured = direct_record(source, ports, (0,) * depth)
        np.testing.assert_allclose(measured, source.conj().T @ power(effect, depth) @ source, atol=1e-12)
        total = sum((direct_record(source, ports, word) for word in product((0, 1), repeat=depth)))
        np.testing.assert_allclose(total, source.conj().T @ source, atol=1e-12)
        records.append(readings(measured, preparations))
    return {"clock_energy_phases": energies.tolist(), "unconditional_records": records}


def calculate():
    unit = np.eye(2)
    e1, e2 = unit
    states = {"first": e1, "second": e2, "plus": (e1 + e2) / np.sqrt(2),
              "phase": (e1 + 1j * e2) / np.sqrt(2)}
    first = np.array([[.65, .10 + .15j], [.10 - .15j, .35]])
    leakage = np.diag([.33, .27])
    second = first @ first + leakage @ leakage
    lower = second @ np.linalg.solve(first, second)
    upper = second - (first - second) @ np.linalg.solve(unit - first, first - second)
    lower_block = leakage @ np.linalg.solve(first, leakage)
    upper_block = unit - leakage @ np.linalg.solve(unit - first, leakage)
    assert np.linalg.eigvalsh(upper_block - lower_block).min() > 0
    source = np.vstack((unit, np.zeros_like(unit)))
    contracts = {}
    for label, block in (("low", lower_block), ("high", upper_block)):
        effect = np.block([[first, leakage], [leakage, block]])
        contracts[label] = model(effect, source, states)
        for depth, moment in enumerate((unit, first, second, {"low": lower, "high": upper}[label])):
            np.testing.assert_allclose(source.T @ power(effect, depth) @ source, moment, atol=1e-12)
    # A prescribed noncentral third moment, realized by an interior physical filter.
    gap = root(upper - lower)
    fraction = np.array([[.6, .1j], [-.1j, .4]])
    requested = lower + gap @ fraction @ gap
    block = lower_block + np.linalg.solve(leakage, requested - lower) @ np.linalg.inv(leakage)
    effect = np.block([[first, leakage], [leakage, block]])
    assert np.linalg.eigvalsh(effect).min() > 0 and np.linalg.eigvalsh(effect).max() < 1
    np.testing.assert_allclose(source.T @ power(effect, 3) @ source, requested, atol=1e-12)

    flat = model(first, unit, states)
    gram = np.block([[unit, first], [first, second]])
    meter_error = .001
    operator_error = (1 + 2 * np.sqrt(2)) * meter_error
    gram_spectrum = np.linalg.eigvalsh(gram)
    assert gram_spectrum.min() > 2 * operator_error
    separation = {k: contracts["high"]["unconditional_records"][3][k]
                  - contracts["low"]["unconditional_records"][3][k] for k in states}
    assert min(separation.values()) > 2 * meter_error
    # Finite-error illustration with independently warranted complete source C0=I.
    beta = operator_error
    leakage_upper = beta + 2 * beta
    holdout_radius = {k: (k * (k - 1) / 2) * leakage_upper + k * beta for k in (3, 4)}
    return {"round": 1350, "time_units": "tau/hbar = 1", "source_identity": True,
            "readout_abs_error": meter_error, "tomography_operator_error": operator_error,
            "flat": flat, "nonflat": contracts, "gram_eigenvalues": gram_spectrum.tolist(),
            "rank_certified_at_stated_error": int(sum(gram_spectrum > 2 * operator_error)),
            "third_record_separation": separation, "flat_prediction_radius_exact_full_source": holdout_radius}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = Path(__file__).with_name("results.json")
    raw = json.dumps(calculate(), indent=2, ensure_ascii=False, allow_nan=False)
    if args.write:
        with path.open("x", encoding="utf8") as out:
            out.write(raw + "\n")
    else:
        def compare(actual, expected):
            if isinstance(actual, dict):
                assert actual.keys() == expected.keys()
                for name in actual:
                    compare(actual[name], expected[name])
            elif isinstance(actual, list):
                assert len(actual) == len(expected)
                for a, b in zip(actual, expected):
                    compare(a, b)
            elif isinstance(actual, float):
                np.testing.assert_allclose(actual, expected, atol=5e-8, rtol=1e-10)
            else:
                assert actual == expected
        compare(json.loads(raw), json.loads(path.read_text(encoding="utf8")))
    values = json.loads(raw)
    print(json.dumps({k: values[k] for k in ("gram_eigenvalues", "third_record_separation",
                                            "flat_prediction_radius_exact_full_source")}))


if __name__ == "__main__":
    main()
