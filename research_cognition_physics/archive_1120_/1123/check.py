"""1123: physical-generator defect closure at fixed local CAR support.

Default execution is read-only.  --write-results creates the initial artifact.
Finite generator matrices check algebra, not a discretized Dirac spectrum.
The physical unbounded-domain and spatial-support claims are analytic.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOL = 3e-11
THETA = math.pi / 3
ETA = 1 / 5
RADIUS = 2.0
MASS = 0.7
SIGMA_SPACE = math.sqrt(11 / 3) / RADIUS
SIGMA_TIME = math.sqrt(11) / RADIUS


def norm(a):
    return float(np.linalg.norm(a, 2))


def check_close(a, b, label):
    error = norm(np.atleast_2d(np.asarray(a) - np.asarray(b)))
    assert error <= TOL, (label, error)
    return error


def car(n):
    out = []
    for j in range(n):
        a = np.zeros((1 << n, 1 << n), dtype=complex)
        for state in range(1 << n):
            if state & (1 << j):
                a[state ^ (1 << j), state] = (-1) ** (
                    state & ((1 << j) - 1)
                ).bit_count()
        out.append(a)
    return out


def q(a, cs):
    result = np.zeros_like(cs[0])
    for j in range(len(cs)):
        for k in range(len(cs)):
            if a[j, k] != 0:
                result += a[j, k] * cs[j].conj().T @ cs[k]
    return result


def phase(h, theta):
    values, vectors = np.linalg.eigh(h)
    return (vectors * np.exp(-1j * theta * values)) @ vectors.conj().T


def mode_mix(cs, i=0, j=1):
    return phase(cs[i].conj().T @ cs[j] + cs[j].conj().T @ cs[i], THETA)


def f_parity(n):
    return np.diag([(-1) ** (s & 3).bit_count() for s in range(1 << n)])


def kraus(u, z):
    return [math.sqrt(1 - ETA) * u, math.sqrt(ETA) * z @ u]


def heisenberg(ks, a):
    return sum((k.conj().T @ a @ k for k in ks), np.zeros_like(a))


def schrodinger(ks, rho):
    return sum((k @ rho @ k.conj().T for k in ks), np.zeros_like(rho))


def projector(v):
    return np.outer(v, v.conj())


def half_trace(a):
    return float(np.abs(np.linalg.eigvalsh((a + a.conj().T) / 2)).sum() / 2)


def cross(n, sigma, first=(0, 1), second=(2, 3)):
    a = np.zeros((n, n), dtype=complex)
    for left, right in zip(first, second):
        a[left, right] = a[right, left] = sigma
    return a


def closure():
    cs = car(6)
    identity = np.eye(64)
    n_f = cs[0].conj().T @ cs[0] + cs[1].conj().T @ cs[1]
    u = mode_mix(cs)
    z = f_parity(6)
    ks = kraus(u, z)
    tp = check_close(sum(k.conj().T @ k for k in ks), identity, "TP")
    gauge = max(norm(k @ n_f - n_f @ k) for k in ks)
    assert gauge < TOL
    cases = []
    for label, mu, sigma in [
        ("time_Dirac_jet", MASS, SIGMA_TIME),
        ("space_Dirac_jet", 0.0, SIGMA_SPACE),
    ]:
        p = np.diag([mu, mu, 0.7, -0.2, 1.3, -0.9]).astype(complex)
        c = cross(6, sigma)
        p += c
        # These real exterior couplings ensure the four-mode block is NOT
        # invariant under this algebra-test generator.
        p[2, 4] = p[4, 2] = 0.37
        p[3, 5] = p[5, 3] = -0.29
        p[2, 3] = 0.11j
        p[3, 2] = -0.11j
        full_p, source = q(p, cs), q(c, cs)
        exterior = full_p - mu * n_f - source
        exterior_commutator = max(norm(exterior @ k - k @ exterior) for k in ks)
        assert exterior_commutator < TOL
        actual = heisenberg(ks, full_p) - full_p
        bounded = heisenberg(ks, source) - source
        ideal = u.conj().T @ full_p @ u - full_p
        residual = check_close(actual, bounded, "unbounded-exterior cancellation")
        source_norm = check_close(norm(source), 2 * sigma, "full CAR source norm")
        target_norm = check_close(
            norm(ideal), 4 * sigma * abs(math.sin(THETA / 2)), "target defect"
        )
        sharp = check_close(norm(actual - ideal), 4 * sigma * ETA, "sharp bound")
        cases.append({
            "label": label,
            "mu": mu,
            "sigma": sigma,
            "source_norm": norm(source),
            "actual_error": norm(actual - ideal),
            "bound": 4 * sigma * ETA,
            "exterior_leakage_norm": norm(p[4:6, 2:4]),
            "exterior_commutator_residual": exterior_commutator,
            "closure_residual": residual,
            "norm_residuals": [source_norm, target_norm, sharp],
        })
    return {"hilbert_dimension": 64, "TP_residual": tp,
            "Kraus_number_commutator": gauge, "cases": cases}


def sharp_and_code():
    cs = car(4)
    u = mode_mix(cs)
    z = f_parity(4)
    ks = kraus(u, z)
    source = q(cross(4, SIGMA_SPACE), cs)
    actual_error = heisenberg(ks, source) - u.conj().T @ source @ u
    source_formula = check_close(
        actual_error, -2 * ETA * u.conj().T @ source @ u, "exact source error"
    )
    code = np.eye(16)[:, [1, 2]]
    code_errors = []
    for i in range(2):
        for j in range(2):
            matrix_unit = np.outer(code[:, i], code[:, j].conj())
            code_errors.append(norm(
                schrodinger(ks, matrix_unit) - u @ matrix_unit @ u.conj().T
            ))
    assert max(code_errors) < TOL

    # One fermionic reference; both summands have even TOTAL parity.
    cs_ref = car(3)
    u_ref = mode_mix(cs_ref)
    z_ref = f_parity(3)
    psi = (np.eye(8)[:, 0] + np.eye(8)[:, 5]) / math.sqrt(2)
    total_parity = np.diag([(-1) ** s.bit_count() for s in range(8)])
    parity_residual = check_close(total_parity @ psi, psi, "even reference input")
    desired = u_ref @ psi
    wrong = z_ref @ desired
    overlap = abs(np.vdot(desired, wrong))
    assert overlap < TOL
    actual = schrodinger(kraus(u_ref, z_ref), projector(psi))
    trace_error = half_trace(actual - projector(desired))
    check_close(trace_error, ETA, "diamond lower witness")
    return {
        "fixed_theta": THETA,
        "mixture_probability": ETA,
        "half_diamond_exact": ETA,
        "half_diamond_justification": "Convex upper bound eta plus explicit orthogonal-output witness.",
        "reference_half_trace_witness": trace_error,
        "reference_total_parity_residual": parity_residual,
        "reference_output_overlap": float(overlap),
        "N_F_one_code_matrix_unit_errors": code_errors,
        "full_source_error": norm(actual_error),
        "linear_bound": 4 * SIGMA_SPACE * ETA,
        "source_formula_residual": source_formula,
    }


def expanded_code():
    cs = car(4)
    u = mode_mix(cs)
    z = f_parity(4)
    ks = kraus(u, z)
    source = q(cross(4, SIGMA_SPACE), cs)
    eigenvalues, vectors = np.linalg.eigh(source)
    outputs = vectors[:, [-1, 0]]
    code = u.conj().T @ outputs
    check_close(code.conj().T @ code, np.eye(2), "expanded code isometry")
    defect_error = heisenberg(ks, source) - u.conj().T @ source @ u
    compressed = code.conj().T @ defect_error @ code
    check_close(
        compressed, np.diag([-4 * SIGMA_SPACE * ETA, 4 * SIGMA_SPACE * ETA]),
        "expanded code source error",
    )
    traces, expectations = [], []
    for j in range(2):
        rho = projector(code[:, j])
        output_difference = schrodinger(ks, rho) - u @ rho @ u.conj().T
        traces.append(half_trace(output_difference))
        expectations.append(float(np.trace(source @ output_difference).real))
        check_close(traces[-1], ETA, "expanded code diamond witness")
        check_close(
            expectations[-1], np.vdot(code[:, j], defect_error @ code[:, j]).real,
            "bounded source duality",
        )
    total_parity = np.diag([(-1) ** s.bit_count() for s in range(16)])
    parity_residual = check_close(total_parity @ code, code, "even expanded code")
    return {
        "code_dimension": 2,
        "retained_output_modes": 4,
        "total_even_parity_residual": parity_residual,
        "half_diamond_exact": ETA,
        "endpoint_half_trace_errors": traces,
        "endpoint_source_errors": expectations,
        "compressed_source_error_norm": norm(compressed),
        "bound": 4 * SIGMA_SPACE * ETA,
        "interpretation": "Finite CAR contract demonstration, not a certified old preparation or readout.",
    }


def hidden_mode():
    # First two modes are visible F.  Modes 2,3 are hidden g, and 4,5
    # are its actual-generator leakage directions.  Orthogonal supports
    # make the separate continuous field algebras have precisely this form.
    cs = car(6)
    u_f = mode_mix(cs, 0, 1)
    u_g = mode_mix(cs, 2, 3)
    combined = u_g @ u_f
    visible_n = cs[0].conj().T @ cs[0] + cs[1].conj().T @ cs[1]
    assert norm(visible_n @ u_g - u_g @ visible_n) < TOL
    vacuum_f = (np.eye(64) - cs[0].conj().T @ cs[0]) @ (
        np.eye(64) - cs[1].conj().T @ cs[1]
    )
    # Explicit complete matrix units of CAR(F), including odd ones.
    creation_words = [
        np.eye(64), cs[0].conj().T, cs[1].conj().T,
        cs[0].conj().T @ cs[1].conj().T,
    ]
    endpoint_errors = []
    for wi in creation_words:
        for wj in creation_words:
            a = wi @ vacuum_f @ wj.conj().T
            endpoint_errors.append(norm(
                combined.conj().T @ a @ combined - u_f.conj().T @ a @ u_f
            ))
    assert max(endpoint_errors) < TOL
    cases = []
    for scale in [1.0, 3.0]:
        source = q(cross(6, scale * SIGMA_SPACE, (2, 3), (4, 5)), cs)
        actual = combined.conj().T @ source @ combined - source
        ideal = u_f.conj().T @ source @ u_f - source
        difference = actual - ideal
        expected = 4 * scale * SIGMA_SPACE * abs(math.sin(THETA / 2))
        check_close(norm(difference), expected, "hidden defect scaling")
        cases.append({"relative_scale": scale, "sigma": scale * SIGMA_SPACE,
                      "hidden_defect_norm": norm(difference), "formula": expected})
    ratio = cases[1]["hidden_defect_norm"] / cases[0]["hidden_defect_norm"]
    check_close(ratio, 3.0, "hidden mode ratio")
    return {
        "visible_CAR_matrix_units_checked": len(endpoint_errors),
        "visible_endpoint_max_residual": max(endpoint_errors),
        "fixed_angle": THETA,
        "cases": cases,
        "defect_ratio": ratio,
        "interpretation": "Actual same-ball smooth scaling is proved analytically; these matrices test its CAR identity.",
    }


def run():
    return {
        "round": 1123,
        "new_axioms": 0,
        "scientific_calibration": 0,
        "tolerance": TOL,
        "physical_mode_diagnostic": {
            "radius": RADIUS, "mass": MASS,
            "sigma_space": SIGMA_SPACE, "sigma_time": SIGMA_TIME,
            "profile": "Reused 1121 H0^1 intersect H^2 radial diagnostic; smooth H1 approximation required.",
        },
        "groups": {
            "full_CAR_closure": closure(),
            "sharp_constant_and_invisible_code": sharp_and_code(),
            "expanded_output_code_certificate": expanded_code(),
            "same_field_hidden_mode": hidden_mode(),
        },
        "analytic_only": [
            "Actual positive-energy Dirac representation and unbounded physical generator domains.",
            "Strong commutation of the residual with the entire local CAR algebra.",
            "Arbitrary passive references and representation multiplicity.",
            "All smooth same-ball hidden modes and unbounded spatial scaling family.",
        ],
        "not_certified": [
            "Actual finite-duration controller, recoil resource, or point Hamiltonian implementation.",
            "Old preparation/readout permissions for all expanded modes.",
            "Local stress, source Hessian, higher energy moments, or complete four-momentum distributions.",
        ],
    }


def compare(actual, saved, path="root"):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys(), path
        for key in actual:
            compare(actual[key], saved[key], path + "." + key)
    elif isinstance(actual, list):
        assert len(actual) == len(saved), path
        for i, (left, right) in enumerate(zip(actual, saved)):
            compare(left, right, path + f"[{i}]")
    elif isinstance(actual, float):
        assert math.isclose(actual, saved, rel_tol=1e-9, abs_tol=1e-11), path
    else:
        assert actual == saved, path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write_results:
        if path.exists():
            raise SystemExit("Refusing to overwrite existing results.json")
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        action = "initial results written"
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
        action = "SAVED_RESULTS_MATCH; read-only"
    print(f"1123 PASS: 4 groups; {action}; new axioms 0; scientific calibration 0")


if __name__ == "__main__":
    main()
