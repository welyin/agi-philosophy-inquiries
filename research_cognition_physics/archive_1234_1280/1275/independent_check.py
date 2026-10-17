"""1275 independent finite-reference two-time audit.

Uses only the inherited S/T swap relations and the three-site spectral
projectors, not author code or results. The 42-dimensional dressed block is
exact for these effects: all 720 flag blocks give the same trace.
Lueders operations and the Gibbs input are mathematical test contracts.
No physical preparation, intermediate readout, or free measurement is granted.
Default: recompute and compare read-only. --write: exclusive first creation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULT = HERE / "independent_results.json"
TOL = 3e-12


def need(condition, message):
    if not bool(condition):
        raise AssertionError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fraction(value):
    return {"numerator": value.numerator, "denominator": value.denominator}


def calculation():
    beta = np.log(2.5) / 2.0
    chain = np.array([[0., 1., 0.], [1., 0., 1.], [0., 1., 0.]])
    b = chain / np.sqrt(2.)
    plus = (b @ b + b) / 2
    minus = (b @ b - b) / 2
    zero = np.eye(3) - b @ b
    need(np.linalg.norm(plus + minus + zero - np.eye(3)) < TOL,
         "three-site spectral completeness")
    e0 = np.diag([1., 0., 0.])

    # B (2+5) tensor W (0,2e), with e=1.  S=P4 tensor |1><1|.
    # T is the two-dimensional singlet-embedding image in P6 tensor |0><0|.
    h0 = np.repeat([4., 4., 6., 6., 6., 6., 6.], 2)
    h0 += np.tile([0., 2.], 7)
    p4 = np.diag([1.] * 4 + [0.] * 10)
    s = np.zeros((14, 14))
    t = np.zeros((14, 14))
    s[1, 1] = s[3, 3] = 1.
    t[4, 4] = t[6, 6] = 1.
    q = p4 - s + t
    delta = s - t
    p = np.kron(np.eye(3), q) + np.kron(e0, delta)
    identity = np.eye(42)
    effects = [p, identity - p]
    need(np.linalg.norm(p @ p - p) < TOL, "dressed area projector")
    need(np.linalg.norm(q @ delta + t) < TOL, "Q(S-T)=-T")
    need(np.linalg.norm(delta @ delta - s - t) < TOL, "(S-T)^2=S+T")
    rd_values = np.exp(-beta * h0)
    rd_values /= rd_values.sum()
    rd = np.diag(rd_values)
    w_s, w_t = np.trace(rd @ s), np.trace(rd @ t)
    need(abs(w_s - 1 / 7) < TOL and abs(w_t - 1 / 7) < TOL,
         "equal resonant reference weights")

    target = np.array([[3 / 7, 1 / 14], [1 / 14, 3 / 7]])
    isolated = np.eye(2) / 2
    cases = []
    for g in [0.25, 1.0, 2.0]:
        need(4.0 > np.sqrt(2.) * g, "inherited positive-spectrum range")
        z = beta * np.sqrt(2.) * g
        zc = 1 + 2 * np.cosh(z)
        rc = (np.exp(-z) * plus + np.exp(z) * minus + zero) / zc
        rho = np.kron(rc, rd)
        h = np.kron(g * chain, np.eye(14)) + np.kron(np.eye(3), np.diag(h0))
        delay = np.pi / (2 * np.sqrt(2.) * g)
        theta = np.sqrt(2.) * g * delay
        uc = np.exp(-1j * theta) * plus + np.exp(1j * theta) * minus + zero
        u = np.kron(uc, np.diag(np.exp(-1j * h0 * delay)))
        evolved = u.conj().T @ p @ u
        correlation = np.trace(rho @ evolved @ p)
        joint = np.array([
            [np.trace(second @ u @ first @ rho @ first @ u.conj().T @ second).real
             for second in effects] for first in effects
        ])
        tv = float(np.abs(joint - isolated).sum() / 2)
        post = sum(effect @ rho @ effect for effect in effects)
        # This positive energy change is why mathematical Lueders insertion
        # must not be silently identified with a free original operation.
        measurement_energy = float(np.trace(h @ (post - rho)).real)
        expected_energy = 2 / 7 * np.sqrt(2.) * g * np.sinh(z) / zc
        expected_imag = -np.sinh(z) / (14 * zc)
        residuals = {
            "gibbs_stationarity": float(np.linalg.norm(h @ rho - rho @ h)),
            "unitarity": float(np.linalg.norm(u.conj().T @ u - identity)),
            "joint_table": float(np.max(np.abs(joint - target))),
            "correlation_real": float(abs(correlation.real - 3 / 7)),
            "correlation_imag": float(abs(correlation.imag - expected_imag)),
            "joint_first_marginal": float(np.max(np.abs(joint.sum(axis=1) - 0.5))),
            "joint_second_marginal": float(np.max(np.abs(joint.sum(axis=0) - 0.5))),
            "tv": abs(tv - 1 / 7),
            "measurement_energy": abs(measurement_energy - expected_energy),
        }
        need(max(residuals.values()) < TOL, "independent 42-dimensional check")
        need(measurement_energy > 0, "first ideal read is not free")
        cases.append({
            "e": 1.0, "g": g, "beta_star": float(beta),
            "delay": float(delay), "joint": joint.tolist(),
            "correlation_real": float(correlation.real),
            "correlation_imag": float(correlation.imag), "tv": tv,
            "ideal_first_read_total_energy_change": measurement_energy,
            "residuals": residuals,
        })

    # The exact margin is a statement about joint distributions, not an
    # instrument construction. Any physical error budget must be supplied.
    epsilon = Fraction(1, 100)
    remaining = Fraction(1, 7) - 2 * epsilon
    need(remaining == Fraction(43, 350), "triangle-inequality margin")
    sources = {
        "../1237/proof.md": digest(HERE.parent / "1237" / "proof.md"),
        "../1274/proof.md": digest(HERE.parent / "1274" / "proof.md"),
        "independent_check.py": digest(Path(__file__)),
    }
    return {
        "status": "PASS_WITH_SCOPE", "source_hashes": sources,
        "method": "independent S/T blocks and exact three-site spectral projectors",
        "author_code_imported": False,
        "exact_joint_table": [[fraction(Fraction(3, 7)), fraction(Fraction(1, 14))],
                              [fraction(Fraction(1, 14)), fraction(Fraction(3, 7))]],
        "exact_tv": fraction(Fraction(1, 7)),
        "each_joint_distribution_tv_error_if_supplied": fraction(epsilon),
        "remaining_tv_margin": fraction(remaining),
        "cases": cases,
        "scope": {
            "mathematical_reference_and_readout_contract_only": True,
            "actual_reference_preparation_certified": False,
            "actual_two_time_instrument_certified": False,
            "full_gibbs_kms_refuted": False,
            "physical_clock_or_horizon_identity_certified": False,
            "geometry_dynamics_selected": False,
        },
    }


def compare(old, new, location="root"):
    if isinstance(new, dict):
        need(isinstance(old, dict) and old.keys() == new.keys(), location + " keys")
        for key in new:
            compare(old[key], new[key], location + "." + key)
    elif isinstance(new, list):
        need(isinstance(old, list) and len(old) == len(new), location + " length")
        for i, value in enumerate(new):
            compare(old[i], value, location + "." + str(i))
    elif isinstance(new, float):
        need(abs(float(old) - new) <= TOL, location + " numeric")
    else:
        need(old == new, location + " value")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = calculation()
    if args.write:
        with RESULT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        compare(json.loads(RESULT.read_text(encoding="utf-8")), result)
    print(json.dumps({"status": result["status"], "cases": len(result["cases"]),
                      "exact_tv": result["exact_tv"],
                      "actual_instrument_certified": False}, ensure_ascii=False))
