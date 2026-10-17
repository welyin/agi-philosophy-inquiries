"""Unnumbered source-gate admission diagnostic after round 1265.

Default: recompute, compare results.json, and do not write.
--write / --write-once: create results.json exclusively on first archival run.
No Fock propagation or autonomous geometry is certified here.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def repository_root() -> Path:
    for candidate in HERE.parents:
        if (candidate / "README.md").is_file() and (
            candidate / "research_cognition_physics" / "archive_1234_" / "1253" / "proof.md"
        ).is_file():
            expected = candidate / "research_cognition_physics" / "archive_1234_" / "_admission" / "source_gate_after1265"
            if expected.resolve() != HERE:
                raise RuntimeError("Script is outside its declared admission directory.")
            return candidate
    raise RuntimeError("Repository markers not found from the script location.")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calculate() -> dict:
    root = repository_root()
    dim = 6
    a = np.zeros((dim, dim), dtype=complex)
    for n in range(1, dim):
        a[n-1, n] = np.sqrt(n)
    D = 0.5j * (a.conj().T @ a.conj().T - a @ a)
    e0 = np.eye(dim, dtype=complex)[:, 0]
    e1 = np.eye(dim, dtype=complex)[:, 1]
    d20 = float(np.vdot(D @ e0, D @ e0).real)
    d21 = float(np.vdot(D @ e1, D @ e1).real)
    assert abs(d20 - 0.5) < 1e-14
    assert abs(d21 - 1.5) < 1e-14
    assert abs(np.vdot(e1, D @ e1)) < 1e-14
    assert abs(np.vdot(e0, D @ e0)) < 1e-14
    # Only |0> -> |2> and |1> -> |3> are used: exact listed moments,
    # not truncated propagation. D changes total N by +/-2, so P1 D P1=0.
    # On N=1, the pair ladder commutator is N + n_real/2.
    real_modes = 1728
    D2_one = F(real_modes, 2) + 1
    assert D2_one == 865
    gamma, T = F(1, 1000), F(1, 64)
    assert gamma < 26
    assert 2 * 32 * T <= 1
    pair_failure = D2_one * gamma**2 * T**2
    full_state_upper = 30 * gamma * T
    assert D2_one < 30**2
    assert pair_failure == F(173, 819200000)
    assert full_state_upper == F(3, 6400)
    # Full zero-point energy is included; no normal-order subtraction.
    E_total_lower = F(26 * 865)
    integrated_A_lower = F(5, 6) * gamma * E_total_lower * T**2
    assert integrated_A_lower == F(2249, 491520)
    assert integrated_A_lower > F(1, 250)
    # Commuting initial constraint identity, not sharp a and Pa quantum data.
    Hg, Ec = F(-100), F(1)
    E_A, E_B = F(20), F(25)
    Pth_A, Pth_B = -Hg - Ec - E_A, -Hg - Ec - E_B
    assert Ec + Pth_A + E_A == Ec + Pth_B + E_B == -Hg
    dependencies = [
        Path(__file__).resolve(),
        HERE / "proof.md",
        root / "research_cognition_physics/archive_1234_/1253/proof.md",
        root / "research_cognition_physics/archive_1234_/1249/canonical.md",
        root / "research_cognition_physics/archive_1234_/1252/independent_review.md",
    ]
    hashes = {p.relative_to(root).as_posix(): sha256(p) for p in dependencies}
    return {
        "schema_version": 1,
        "status": "PASS_UNNUMBERED_COMPARISON_ONLY",
        "numbered_round": None,
        "new_scientific_groups": 0,
        "bound_files_sha256": hashes,
        "original_temporary_provenance_sha256": {
            "clock_source_gate_after1265.md": "b4eb32a37c21af2293466208acef4f44ded2a36c6e5482237cff5ef00f68e863",
            "clock_source_gate_after1265.py": "8867e98a86ca0423aab05ca31ee93e7db10790950993992b9864a6dbb8b5a601",
        },
        "real_mode_single_moment": {
            "vacuum_D2": d20,
            "one_particle_D2": d21,
            "authoritative_exact_values": ["1/2", "3/2"],
            "float_comparison_absolute_tolerance": 1e-13,
        },
        "full_real_modes": real_modes,
        "P1_D_P1_zero": True,
        "P1_D2_P1_coefficient": str(D2_one),
        "constant_dilation_comparison": {
            "gamma": str(gamma),
            "T": str(T),
            "pair_failure_upper": str(pair_failure),
            "full_state_error_strict_upper": str(full_state_upper),
            "minus_integral_of_volume_times_pressure_strict_lower": str(integrated_A_lower),
        },
        "initial_energy_compensation": {
            "delta_KG": str(E_B - E_A),
            "delta_clock": str(Pth_B - Pth_A),
            "delta_full_material_energy": "0",
        },
        "scope": [
            "prescribed constant gamma comparison, not autonomous feedback",
            "not a finite self-consistent geometry certificate",
            "not a propagation calculation with a Fock cutoff",
            "one-particle projection is not invariant under the full Hamiltonian",
            "fixed original mass weights alone do not supply a new geometry-source contrast",
            "stop pair engineering; no numbered scientific addition",
        ],
    }


def compare(expected, actual, label="results"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"{label}: keys differ")
        for key in expected:
            compare(expected[key], actual[key], f"{label}.{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"{label}: sequence differs")
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{label}[{i}]")
    elif isinstance(expected, float):
        if not isinstance(actual, (int, float)) or not math.isclose(
            expected, actual, rel_tol=0.0, abs_tol=1e-13
        ):
            raise AssertionError(f"{label}: float differs")
    elif expected != actual:
        raise AssertionError(f"{label}: value differs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", "--write-once", dest="write_once", action="store_true")
    args = parser.parse_args()
    result = calculate()
    result_path = HERE / "results.json"
    if args.write_once:
        with result_path.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        receipt = "FIRST_WRITE_EXCLUSIVE"
    else:
        with result_path.open("r", encoding="utf-8") as handle:
            saved = json.load(handle)
        compare(saved, result)
        receipt = "READ_ONLY_MATCH"
    print(json.dumps({
        "status": result["status"],
        "receipt": receipt,
        "new_scientific_groups": 0,
        "results": str(result_path),
        "results_sha256": sha256(result_path),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
