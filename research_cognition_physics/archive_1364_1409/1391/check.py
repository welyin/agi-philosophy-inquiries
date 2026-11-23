"""Recompute the single-carrier two-phase witness. Run with --check or --write-new."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def quartet(left):
    return float((left[0, 0] * left[1, 1] * left[0, 1].conjugate()
                  * left[1, 0].conjugate()).imag)


def observe(upper, row, scale, phi):
    mass = np.vstack((upper, scale * np.exp(1j * phi) * row))
    left, singular, right_h = np.linalg.svd(mass)
    order = np.argsort(singular)
    left = left[:, order]
    singular = singular[order]
    right = right_h.conj().T[:, order]
    weak = left[:3, :]
    neutral = weak.conj().T @ weak
    singlet = left[3:, :].conj().T @ left[3:, :]
    projected = np.diag([1, 1, 1, 0])
    higgs = left.conj().T @ projected @ mass @ right
    pq = left.conj().T @ (np.eye(4) - projected) @ mass @ right
    diagonal = np.diag(singular)
    assert np.allclose(higgs, neutral @ diagonal, atol=2e-10)
    assert np.allclose(pq, singlet @ diagonal, atol=2e-10)
    assert np.allclose(higgs + pq, diagonal, atol=2e-10)
    assert np.linalg.matrix_rank(singlet, tol=1e-8) == 1
    ds2 = abs(singlet[0, 1])**2
    assert np.isclose(ds2, singlet[0, 0].real * singlet[1, 1].real, atol=1e-20, rtol=1e-10)
    gram = mass @ mass.conj().T
    cycle = (gram[0, 1] * gram[1, 3] * gram[3, 0]).imag
    return {"masses": singular.tolist(), "J_weak": quartet(weak),
            "cycle_over_scale_sq": float(cycle / scale**2),
            "light_cabibbo_current_singlet_trace": float(np.trace(singlet[:3, :3]).real),
            "abs_ds": float(abs(singlet[0, 1])), "arg_det": float(np.angle(np.linalg.det(mass)))}


def calculate(data):
    model = data["model"]
    upper = np.column_stack((model["m"], model["c"]))
    delta = math.acos(model["cos_delta"])
    phase = np.exp(1j * delta)
    row = np.r_[np.array(model["bS"]) + phase * np.array(model["bT"]), model["heavy_yukawa"]]
    det_factor = 4 - 0.5 * phase
    phi_star = -float(np.angle(det_factor))
    records = {}
    for scale in model["scalings"]:
        record = observe(upper, row, scale, phi_star)
        assert abs(record["arg_det"]) < 1e-12
        assert abs(record["cycle_over_scale_sq"] + 2 * math.sin(delta)) < 1e-12
        assert abs(record["J_weak"]) > 1.e-5
        for phi in model["phi_checks"]:
            alternate = observe(upper, row, scale, phi)
            assert np.allclose(record["masses"], alternate["masses"])
            assert abs(alternate["J_weak"] - record["J_weak"]) < 1e-10
            angle_error = np.angle(np.exp(1j * (alternate["arg_det"] - phi - np.angle(det_factor))))
            assert abs(angle_error) < 1.e-12
        records[str(scale)] = record
    cp_partner = observe(upper, row.conj(), 1, -phi_star)
    assert np.isclose(cp_partner["J_weak"], -records["1"]["J_weak"])
    # Isolate the light eigensystem as the complete heavy row is scaled to infinity.
    coupling = upper @ row.conj()
    light_gram = upper @ upper.T - np.outer(coupling, coupling.conj()) / np.vdot(row, row).real
    _, light_left = np.linalg.eigh(light_gram)
    limiting_J = quartet(light_left)
    assert abs(limiting_J - records["1000"]["J_weak"]) < 1e-7

    v_s, v_t = model["radial_vev_S"], model["radial_vev_T"]
    decay_sq = 2 * (v_s**2 + v_t**2)
    relative_sq = 2 * v_s**2 * v_t**2 / (v_s**2 + v_t**2)
    weight = v_t**2 / (v_s**2 + v_t**2)
    slope = float((-0.5j * phase / det_factor).imag - weight)
    k_delta = 2 * model["positive_lock"] * (v_s * v_t * math.sin(delta))**2
    chi = model["chi_example"]
    angle_hessian = chi * np.outer([1, slope], [1, slope]) + np.diag([0, k_delta])
    kinetic_inv_root = np.diag([1/math.sqrt(decay_sq), 1/math.sqrt(relative_sq)])
    mass_sq = np.linalg.eigvalsh(kinetic_inv_root @ angle_hessian @ kinetic_inv_root)
    assert np.isclose(np.linalg.det(angle_hessian), chi * k_delta) and min(mass_sq) > 0
    epsilon = model["energy_spread_example"]
    radius = math.sqrt(epsilon / (model["positive_lock"] * (v_s * v_t)**2))
    cos_interval = [model["cos_delta"] - radius, model["cos_delta"] + radius]
    kaon = data["kaon"]
    m_k, m_pi = kaon["mass_GeV"], kaon["pion_mass_GeV"]
    coefficient = (m_k**2 - m_pi**2)**3 / (64 * math.pi * m_k**3)
    scale_limit = math.sqrt(coefficient * kaon["lifetime_s"] /
                            (kaon["hbar_GeV_s"] * kaon["massless_line_B90"]))
    return {"delta":delta, "phi_relaxed":phi_star, "pq_color_anomaly":1, "em_over_color":2/3,
            "samples":records, "conjugate_J":cp_partner["J_weak"], "decoupling_J":limiting_J,
            "example_mass_sq":mass_sq.tolist(),"stiff_estimate_chi_over_F2":chi/decay_sq,
            "cos_interval_given_spread":cos_interval,
            "min_abs_sin_given_spread":math.sqrt(1-max(abs(x) for x in cos_interval)**2),
            "massless_K_F_over_ds_f0_lower_GeV":scale_limit}


def equivalent(observed, stored):
    if isinstance(observed, dict):
        return observed.keys() == stored.keys() and all(equivalent(v, stored[k]) for k, v in observed.items())
    return bool(np.allclose(observed, stored, rtol=1e-8, atol=1e-12))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    values = calculate(json.loads((ROOT / "inputs.json").read_text(encoding="utf-8-sig")))
    result_path = ROOT / "results.json"
    if args.write_new:
        with result_path.open("x", encoding="utf-8") as target:
            json.dump(values, target, ensure_ascii=False, indent=2, allow_nan=False)
            target.write("\n")
    else:
        assert equivalent(values, json.loads(result_path.read_text(encoding="utf-8-sig")))
        print("1391 reproducibility check passed")
