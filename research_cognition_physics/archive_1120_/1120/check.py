"""1120: necessary complete interaction action for the original 1111 window.

Default: read-only calculation and saved-result comparison. --write-results
writes only this round's results.json. Numeric checks do not certify physical
controls, local stress, a lapse operation, or a model of all inherited axioms.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / "archive_1106_1119" / "1111" / "check.py"
RESULTS = HERE / "results.json"
TOL = 2e-10


def residual(a, label):
    r = float(np.linalg.norm(a))
    if r > TOL:
        raise AssertionError(f"{label}: {r}")
    return r


def evolve(h, t):
    vals, vecs = np.linalg.eigh(h)
    return (vecs * np.exp(-1j * t * vals)) @ vecs.conj().T


def halfwidth(h):
    vals = np.linalg.eigvalsh(h)
    return float((vals[-1] - vals[0]) / 2)


def assert_close(x, y, label):
    return residual(np.array(x) - np.array(y), label)


def pure_distance(a, b):
    return math.sqrt(max(0.0, 1.0 - abs(np.vdot(a, b)) ** 2))


def calculate():
    spec = importlib.util.spec_from_file_location("round1111_reused", OLD)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    m = old.model()
    gamma, hd, hc = m["gamma"], m["hout"], m["hc"]
    pc = np.kron(np.eye(8), np.kron(np.diag([0, 0, 1]), np.eye(3)))
    pd = gamma @ np.kron(np.eye(4), pc) @ gamma.conj().T
    delta = math.pi / 8
    # Reuse the actual old ready columns and the proved Gamma intertwiner.
    tau = math.pi - delta
    ws = gamma @ (np.kron(evolve(hc + m["baseline"] * np.eye(4), tau),
                         np.eye(72)) @ (m["phase"] * m["ready"]))
    r_intertwine = residual(hd @ gamma - gamma @ m["hin"], "old Gamma intertwiner")
    r_unitary = residual(gamma.conj().T @ gamma - np.eye(288), "Gamma unitary")
    r_source = residual(pd @ ws - ws @ np.diag([0, 1]), "actual reachable source")
    r_iso = residual(ws.conj().T @ ws - np.eye(2), "original window input isometry")
    x, y = np.eye(2, dtype=complex).T
    v1 = np.kron(ws[:, 1], x)
    s1 = np.kron(ws[:, 1] - pd @ ws[:, 1], x) + np.kron(pd @ ws[:, 1], y)
    r_orthogonal = residual(np.vdot(v1, s1), "same-task orthogonal material output")
    r_y = residual(s1 - np.kron(ws[:, 1], y), "old controlled material y")

    # Exact logical source x material reduction after retaining the old clock
    # and dressing. This is a reduction of this action, not a new apparatus.
    identity = np.eye(4, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    S = np.zeros((4, 4), dtype=complex)
    S[:2, :2] = np.eye(2)
    S[2:, 2:] = X
    D = 4 * (identity - S)
    J = identity[:, [0, 2]]
    vq0, vq1 = J[:, 0], J[:, 1]
    r_boundary = residual(evolve(D, delta) - S, "full logical boundary identity")
    r_involution = residual(S @ S - identity, "involution")
    q = halfwidth(D)
    assert_close(q, 4.0, "sharp average halfwidth")
    py4 = np.diag([0, 1, 0, 1]).astype(complex)
    cross_star = float(np.linalg.norm(py4 @ D @ (identity - py4), 2))
    assert_close(cross_star, 4.0, "sharp material-crossing budget")
    expected = [float(np.vdot(v, D @ v).real) for v in (vq0, vq1)]
    assert_close(expected, [0.0, 4.0], "source-branch energy contrast")
    shift_checks = []
    for c in (-9.0, -4.0, 2.5):
        centered = D + c * identity
        contrast = float(np.vdot(vq1, centered @ vq1).real - np.vdot(vq0, centered @ vq0).real)
        shift_checks.append({"common_scalar": c, "halfwidth": halfwidth(centered),
                             "branch_contrast": contrast})
        assert_close(halfwidth(centered), 4.0, "scalar invariant action")
        assert_close(contrast, 4.0, "scalar invariant contrast")
    hmin = float(np.linalg.eigvalsh(hc)[0] + m["baseline"] + 10)
    if hmin <= 0:
        raise AssertionError("old positive energy")

    epsilon_rows = []
    for eps in (0.0, 0.001, 0.01, 0.1, 0.5, 1.0):
        bound = math.acos(math.sqrt(eps))
        de = (bound / delta) * (identity - S)
        state = evolve(de, delta) @ vq1
        p_y = float(abs(state[3]) ** 2)
        assert_close(p_y, 1 - eps, "probability bound saturation")
        assert_close(delta * halfwidth(de), bound, "complete action for saturation")
        epsilon_rows.append({"epsilon_port_probability": eps, "minimum_action": bound,
                             "minimum_average_halfwidth": bound / delta,
                             "mathematical_witness_y_probability": p_y})

    # Auxiliary witness deliberately changes a0 to a1. It does not claim a
    # complete original channel realization after discarding the auxiliary.
    a = np.eye(8, dtype=complex)
    src, dst = a[:, 4], a[:, 7]  # q1,x,a0 -> q1,y,a1
    active = np.outer(src, src.conj()) + np.outer(dst, dst.conj())
    flip = np.outer(src, dst.conj()) + np.outer(dst, src.conj())
    daux = 4 * (active - flip)
    auxout = evolve(daux, delta) @ src
    aux_res = residual(auxout - dst, "auxiliary may retain changed record")
    assert_close(delta * halfwidth(daux), math.pi / 2, "full auxiliary action")

    # Deterministic, finite noncommuting pulse examples exercise the complete
    # inequality. The analytic FS proof, not this grid, covers arbitrary K(t).
    rng = np.random.default_rng(1120)
    href = np.kron(np.eye(4), np.diag([0.0, 0.3]))
    py = np.diag([0, 0, 1, 1, 0, 0, 1, 1]).astype(complex)
    residual(href @ py - py @ href, "reference preserves original material projection")
    identity8 = np.eye(8, dtype=complex)
    aux_cross = float(np.linalg.norm(py @ daux @ (identity8 - py), 2))
    assert_close(aux_cross, 4.0, "auxiliary-assisted crossing budget")
    diagonal_rows = []
    for amount in (4.0, 1000.0):
        diagonal = amount * py
        cross = float(np.linalg.norm(py @ diagonal @ (identity8 - py), 2))
        state = evolve(href + diagonal, delta) @ src
        probability = float(np.vdot(state, py @ state).real)
        assert_close(cross, 0.0, "diagonal correction no crossing")
        assert_close(probability, 0.0, "diagonal energy cannot complete material task")
        diagonal_rows.append({"correction_halfwidth": halfwidth(diagonal),
                              "crossing_norm": cross, "y_probability": probability})
    pulse_rows = []
    for budget in (0.2, 0.7, 1.2, 1.5):
        psi = src.copy()
        dt = delta / 3
        accounted = 0.0
        cross_accounted = 0.0
        for _ in range(3):
            raw = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
            hermitian = (raw + raw.conj().T) / 2
            correction = hermitian * (budget / delta) / halfwidth(hermitian)
            accounted += dt * halfwidth(correction)
            cross_accounted += dt * float(np.linalg.norm(py @ correction @ (identity8 - py), 2))
            psi = evolve(href + correction, dt) @ psi
        p_y = float(np.vdot(psi, py @ psi).real)
        upper = math.sin(accounted) ** 2
        cross_upper = math.sin(cross_accounted) ** 2
        if p_y > upper + TOL:
            raise AssertionError("complete noncommuting-pulse action bound")
        if p_y > cross_upper + TOL or cross_accounted > accounted + TOL:
            raise AssertionError("complete material-crossing bound")
        pulse_rows.append({"action": accounted, "y_probability": p_y,
                           "material_crossing_action": cross_accounted,
                           "success_upper_bound": upper, "crossing_success_upper_bound": cross_upper,
                           "bound_slack": upper - p_y, "crossing_bound_slack": cross_upper - p_y})

    response_rows = []
    max_gram_residual = 0.0
    for branch in (0, 1):
        multiple = 2 * branch + 1
        db = multiple * D
        residual(evolve(db, delta) - S, "same original boundary for each branch")
        for eta in (-0.2, -0.1, 0.01, 0.1, 0.2):
            n = 1 + eta
            ua, ub = evolve(db, n * delta) @ J, S @ J
            angle = multiple * math.pi * eta / 2
            z = np.exp(-1j * angle) * math.cos(angle)
            gram = ub.conj().T @ ua
            rg = residual(gram - np.diag([1, z]), "full source Gram identity")
            max_gram_residual = max(max_gram_residual, rg)
            rdist = pure_distance(ua[:, 1], ub[:, 1])
            exact = abs(math.sin(angle))
            assert_close(rdist, exact, "arbitrary-reference task distance witness")
            wrong = float(abs(ua[2, 1]) ** 2)
            assert_close(wrong, exact ** 2, "same material wrong port")
            # For all p, |p+(1-p)z|^2 = |z|^2+p^2(1-|z|^2).
            # This algebraic coefficient identity proves the minimum at p=0.
            residual(z.real - abs(z) ** 2, "Gram convex-hull coefficient")
            response_rows.append({"branch": branch, "N": n, "action_at_nominal_window":
                                  delta * halfwidth(db), "source_branch_energy_contrast":
                                  float(np.vdot(vq1, db @ vq1).real),
                                  "wrong_port_probability": wrong,
                                  "exact_complete_task_half_diamond": exact})

    files = {"source": Path(__file__), "old_source": OLD,
             "proof": HERE / "proof.md", "selection_audit": HERE / "selection_audit.md",
             "next": HERE / "NEXT.md"}
    return {
        "schema": "round1120_window_source_budget_v1",
        "round": 1120,
        "scientific_calibration_increment": 0,
        "delta": delta,
        "sha256": {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in files.items()},
        "group1_original_1111_material_and_reachable_input": {
            "old_apparatus_dimension": 288, "original_product_mother_dimension": 1152,
            "gamma_intertwiner_residual": r_intertwine,
            "gamma_unitary_residual": r_unitary, "same_source_code_residual": r_source,
            "same_source_isometry_residual": r_iso, "orthogonality_residual": r_orthogonal,
            "material_y_identity_residual": r_y},
        "group2_sharp_same_window_mathematical_refinement": {
            "logical_action_dimension": 4, "involution_residual": r_involution,
            "boundary_residual": r_boundary, "extra_halfwidth": q,
            "material_crossing_norm": cross_star,
            "non_scalar_action": delta * q, "reference_energy_minimum": hmin,
            "correction_energy_expectations_q0_q1": expected,
            "scalar_invariance_checks": shift_checks},
        "group3_finite_error_necessary_bounds": epsilon_rows,
        "group4_complete_auxiliary_and_noncommuting_candidates": {
            "auxiliary_changes_and_remains_retained_residual": aux_res,
            "auxiliary_complete_action": delta * halfwidth(daux),
            "auxiliary_crossing_action": delta * aux_cross,
            "diagonal_energy_cannot_replace_material_crossing": diagonal_rows,
            "noncommuting_pulse_rows": pulse_rows},
        "group5_time_response_not_fixed_by_one_boundary": {
            "maximum_source_gram_residual": max_gram_residual, "rows": response_rows},
        "scope": {
            "universal_quantifiers_supported_by_analytic_proof": True,
            "new_actual_time_variation_permission_proved": False,
            "actual_internal_controller_or_local_source_constructed": False,
            "K_star_identified_with_original_H_energy": False,
            "whole_axiom_countermodel_claimed": False,
            "stage_completion_claimed": False}}


def compare(saved, actual, path="root"):
    if isinstance(actual, dict):
        if not isinstance(saved, dict) or saved.keys() != actual.keys():
            raise AssertionError(path + ": keys")
        for k in actual:
            compare(saved[k], actual[k], path + "." + k)
    elif isinstance(actual, list):
        if not isinstance(saved, list) or len(saved) != len(actual):
            raise AssertionError(path + ": length")
        for i, (s, a) in enumerate(zip(saved, actual)):
            compare(s, a, f"{path}[{i}]")
    elif isinstance(actual, float):
        if not math.isclose(saved, actual, rel_tol=3e-8, abs_tol=3e-9):
            raise AssertionError(f"{path}: {saved} != {actual}")
    elif saved != actual:
        raise AssertionError(f"{path}: {saved} != {actual}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = calculate()
    if args.write_results:
        RESULTS.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print("WROTE results.json")
    else:
        compare(json.loads(RESULTS.read_text(encoding="utf-8")), result)
        print("SAVED_RESULTS_MATCH")
    print("PASS: 5 groups; same 1111 task, finite-error action and complete source response")
    print("Scientific calibration: 0; physical internal implementation and stage completion not claimed")


if __name__ == "__main__":
    main()
