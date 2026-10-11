"""1110 finite-box four-momentum bridge. Default mode is read-only.

Only --write writes results.json in this new round directory.
No plotting, external dependencies, runtime installation, or dense large-bank W.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results.json"
ATOL = 1e-9
RTOL = 1e-8


def maxabs(a):
    return float(np.max(np.abs(a))) if np.size(a) else 0.0


def opnorm(a):
    gram = a.conj().T @ a
    return float(np.sqrt(max(0.0, np.linalg.eigvalsh((gram + gram.conj().T) / 2)[-1])))


def must_small(value, label, tol=1e-10):
    if not np.isfinite(value) or value > tol:
        raise AssertionError(f"{label}: residual {value} > {tol}")
    return float(value)


def must_le(value, bound, label, tol=1e-10):
    if not np.isfinite(value) or value > bound + tol:
        raise AssertionError(f"{label}: {value} > {bound}")
    return float(value)


def rotation(d, i, j, angle):
    a = np.eye(d, dtype=complex)
    c, s = np.cos(angle), np.sin(angle)
    a[i, i] = a[j, j] = c
    a[i, j], a[j, i] = -s, s
    return a


def dense_unitaries():
    d = 3
    k = np.arange(d)
    f = np.exp(2j * np.pi * np.outer(k, k) / d) / np.sqrt(d)
    u = f @ np.diag(np.exp(1j * np.array([0.17, -0.31, 0.23])))
    u = u @ rotation(d, 0, 1, 0.37) @ rotation(d, 1, 2, -0.22)
    v = rotation(d, 0, 2, 0.29)
    v = v @ np.diag(np.exp(1j * np.array([0.19, 0.43, -0.47]))) @ f.conj().T
    return u, v


def labels(d, m):
    return np.array(list(itertools.product(range(m + 2), repeat=d)), dtype=int)


def indices(ns, d, m):
    ns = np.asarray(ns, dtype=int)
    if ns.ndim == 1:
        ns = ns[None, :]
    return np.ravel_multi_index(ns.T, (m + 2,) * d)


def bank(d, m, p, q0):
    ns = labels(d, m)
    q = q0[None, :] + ns @ p
    eta = np.zeros(len(ns), dtype=complex)
    inside = np.all((ns >= 1) & (ns <= m), axis=1)
    eta[inside] = m ** (-d / 2)
    return ns, q, eta


def full_w(u, m):
    d = len(u)
    nr = (m + 2) ** d
    w = np.eye(d * nr, dtype=complex)
    occupied = np.zeros(d * nr, dtype=bool)
    complete_count = 0
    for k in itertools.product(range(1, m + 2), repeat=d):
        k = np.array(k, dtype=int)
        ids = []
        for i in range(d):
            n = k.copy()
            n[i] -= 1
            ids.append(i * nr + int(indices(n, d, m)[0]))
        w[np.ix_(ids, ids)] = u
        occupied[ids] = True
        complete_count += 1
    return w, occupied, complete_count


def rectangular_v(u, m):
    d = len(u)
    nr = (m + 2) ** d
    inside = np.array(list(itertools.product(range(1, m + 1), repeat=d)), dtype=int)
    eta = np.zeros(nr, dtype=complex)
    eta[indices(inside, d, m)] = m ** (-d / 2)
    v = np.zeros((d, nr, d), dtype=complex)
    for i in range(d):
        for j in range(d):
            shifted = inside.copy()
            shifted[:, i] += 1
            shifted[:, j] -= 1
            if np.min(shifted) < 0 or np.max(shifted) > m + 1:
                raise AssertionError("shift escaped the finite bank")
            v[j, indices(shifted, d, m), i] = u[j, i] * m ** (-d / 2)
    v0 = np.einsum("ji,r->jri", u, eta)
    return v.reshape(d * nr, d), v0.reshape(d * nr, d), eta


def momentum_diagonals(p, q):
    return (p[:, None, :] + q[None, :, :]).reshape(-1, 4)


def moment_compressions(v, p, q, eta):
    d, nr = len(p), len(q)
    mean = np.abs(eta) ** 2 @ q
    sys_diags = np.broadcast_to(p[:, None, :], (d, nr, 4)).reshape(-1, 4)
    ref_diags = np.broadcast_to(q[None, :, :], (d, nr, 4)).reshape(-1, 4)
    da, dr, balance = [], [], []
    for mu in range(4):
        a = v.conj().T @ (sys_diags[:, mu, None] * v) - np.diag(p[:, mu])
        r = v.conj().T @ (ref_diags[:, mu, None] * v) - mean[mu] * np.eye(d)
        da.append(a)
        dr.append(r)
        balance.append(maxabs(a + r))
    return np.array(da), np.array(dr), mean, balance


def run():
    u, v = dense_unitaries()
    p = np.array([[3.0, 0.0, 0.0, 0.0],
                  [4.0, 1.0, 0.0, 1.0],
                  [5.0, 0.0, 2.0, -1.0]])
    q0 = np.array([7.0, 0.5, -0.5, 0.25])
    groups = {}

    # 1. Complete finite mother-space construction, including every residual block.
    m = 2
    w_u, occupied, count = full_w(u, m)
    w_v, _, _ = full_w(v, m)
    w_uv, _, _ = full_w(u @ v, m)
    w_ud, _, _ = full_w(u.conj().T, m)
    ident = np.eye(len(w_u), dtype=complex)
    nr = (m + 2) ** 3
    ve, v0, eta = rectangular_v(u, m)
    input_embedding = np.einsum("ji,r->jri", np.eye(3), eta).reshape(3 * nr, 3)
    residual_columns = np.flatnonzero(~occupied)
    residual_identity = maxabs(w_u[:, residual_columns] - ident[:, residual_columns])
    checks = {
        "unitarity": maxabs(w_u.conj().T @ w_u - ident),
        "representation": maxabs(w_uv - w_u @ w_v),
        "inverse": maxabs(w_ud - w_u.conj().T),
        "residual_identity": residual_identity,
        "finite_vs_shift_isometry": maxabs(w_u @ input_embedding - ve),
    }
    for key, val in checks.items():
        must_small(val, "finite block " + key)
    groups["01_finite_blocks_and_representation"] = {
        "status": "PASS", "d": 3, "M": m, "mother_dimension": len(w_u),
        "complete_blocks": count, "residual_basis_vectors": len(residual_columns),
        "residuals": checks,
    }

    # 2. Non-collinear material momenta; exact algebra is checked on the full mother.
    ns, q, eta = bank(3, m, p, q0)
    mass2_p = p[:, 0] ** 2 - np.sum(p[:, 1:] ** 2, axis=1)
    mass2_q = q[:, 0] ** 2 - np.sum(q[:, 1:] ** 2, axis=1)
    if np.min(mass2_p) <= 0 or np.min(mass2_q) <= 0 or np.min(q[:, 0]) <= 0:
        raise AssertionError("non-positive material mass or future direction")
    total = momentum_diagonals(p, q)
    comm = [maxabs((total[:, mu, None] - total[None, :, mu]) * w_u) for mu in range(4)]
    for mu, val in enumerate(comm):
        must_small(val, f"full four-momentum commutator {mu}")
    bare = [maxabs(np.diag(p[:, mu]) @ u - u @ np.diag(p[:, mu])) for mu in range(4)]
    if min(bare) <= 1e-3:
        raise AssertionError("test U should be non-conserving in all four bare components")
    groups["02_positive_noncollinear_four_momentum"] = {
        "status": "PASS", "p_i": p.tolist(), "q0": q0.tolist(),
        "min_system_mass_squared": float(np.min(mass2_p)),
        "min_bank_mass_squared": float(np.min(mass2_q)),
        "spatial_catalog_rank": int(np.linalg.matrix_rank(p[1:, 1:] - p[0, 1:])),
        "full_commutator_residuals": comm, "bare_U_commutator_norm_entries": bare,
    }

    # 3. Rectangular isometries and an independent exact cross-Gram identity.
    rows = []
    for large_m in [2, 4, 8, 16]:
        vm, ideal, et = rectangular_v(u, large_m)
        alpha = (1 - 1 / large_m) ** 2
        diagonal_u = np.diag(np.diag(u))
        gram_formula = alpha * np.eye(3) + (1 - alpha) * u.conj().T @ diagonal_u
        cross = ideal.conj().T @ vm
        diff_sq = (1 - alpha) * (2 * np.eye(3) - u.conj().T @ diagonal_u - diagonal_u.conj().T @ u)
        exact_norm = float(np.sqrt(max(0.0, np.linalg.eigvalsh((diff_sq + diff_sq.conj().T) / 2)[-1])))
        err = opnorm(vm - ideal)
        isores = maxabs(vm.conj().T @ vm - np.eye(3))
        gramres = maxabs(cross - gram_formula)
        must_small(isores, "rectangular isometry")
        must_small(gramres, "independent cross Gram")
        must_small(abs(err - exact_norm), "independent exact norm")
        bound = 2 * np.sqrt(2 / large_m)
        must_le(err, bound, "full-output uniform bound")
        rows.append({
            "M": large_m, "reference_dimension": (large_m + 2) ** 3,
            "isometry_rows": len(vm), "full_output_operator_error": err,
            "independent_gram_error": exact_norm, "proved_operator_bound": float(bound),
            "proved_half_diamond_bound": float(min(1.0, bound)),
            "isometry_residual": isores, "cross_gram_residual": gramres,
        })
    groups["03_rectangular_full_output_error"] = {
        "status": "PASS", "dense_U": "DFT3 times phases and two real plane rotations",
        "rows": rows, "no_numerical_diamond_optimization_claimed": True,
    }

    # 4. Whole-input operator recoil balance; state closeness does not erase exchange.
    recoil_rows = []
    expected_probe = p[0] - np.abs(u[:, 0]) ** 2 @ p
    for recoil_m in [4, 16]:
        _, qr, et = bank(3, recoil_m, p, q0)
        vr, ideal, _ = rectangular_v(u, recoil_m)
        da, dr, mean, bal = moment_compressions(vr, p, qr, et)
        for val in bal:
            must_small(val, "all-input momentum exchange", tol=2e-10)
        actual_probe = np.real(dr[:, 0, 0])
        must_small(maxabs(actual_probe - expected_probe), "probe exact recoil", tol=2e-10)
        tensor = vr.reshape(3, len(qr), 3)
        fidelity = float(sum(abs(np.vdot(et, tensor[j, :, 0])) ** 2 for j in range(3)))
        if not 0 <= fidelity <= 1 + 1e-10:
            raise AssertionError("reference fidelity outside range")
        recoil_rows.append({
            "M": recoil_m, "all_input_operator_balance_residuals": bal,
            "probe_reference_delta_four_momentum": actual_probe.tolist(),
            "probe_reference_fidelity_with_initial": fidelity,
            "initial_reference_mean": mean.tolist(),
            "full_output_operator_error": opnorm(vr - ideal),
        })
    if np.linalg.norm(expected_probe) < 0.1:
        raise AssertionError("recoil example accidentally trivial")
    groups["04_reference_recoil_and_all_input_balance"] = {
        "status": "PASS", "rows": recoil_rows,
        "same_nonzero_probe_recoil_at_all_M": expected_probe.tolist(),
    }

    # 5. Same material family under a fixed Lorentz transport, not a finite carrier boost.
    beta = 0.6
    gamma = 1 / np.sqrt(1 - beta ** 2)
    lam = np.eye(4)
    lam[0, 0] = lam[1, 1] = gamma
    lam[0, 1] = lam[1, 0] = gamma * beta
    metric = np.diag([1.0, -1.0, -1.0, -1.0])
    pp, qq0 = p @ lam.T, q0 @ lam.T
    _, qp, etap = bank(3, m, pp, qq0)
    qp_expected = q @ lam.T
    totalp = momentum_diagonals(pp, qp)
    commp = [maxabs((totalp[:, mu, None] - totalp[None, :, mu]) * w_u) for mu in range(4)]
    weight = np.abs(eta) ** 2
    mean = weight @ q
    centered = q - mean
    cov = (centered.T * weight) @ centered
    meanp = np.abs(etap) ** 2 @ qp
    cp = qp - meanp
    covp = (cp.T * weight) @ cp
    mean_formula = q0 + (m + 1) / 2 * np.sum(p, axis=0)
    cov_formula = (m * m - 1) / 12 * p.T @ p
    lres = {
        "metric": maxabs(lam.T @ metric @ lam - metric),
        "catalog": maxabs(qp - qp_expected),
        "mass_squared": maxabs(qp[:, 0] ** 2 - np.sum(qp[:, 1:] ** 2, axis=1) - mass2_q),
        "mean_formula": maxabs(mean - mean_formula),
        "covariance_formula": maxabs(cov - cov_formula),
        "mean_transport": maxabs(meanp - lam @ mean),
        "covariance_transport": maxabs(covp - lam @ cov @ lam.T),
    }
    for key, val in lres.items():
        must_small(val, "Lorentz " + key, tol=1e-9)
    for val in commp:
        must_small(val, "boosted conservation", tol=1e-9)
    groups["05_lorentz_transport_and_resource_moments"] = {
        "status": "PASS", "beta": beta, "residuals": lres,
        "boosted_commutator_residuals": commp,
        "initial_reference_mean": mean.tolist(), "transported_mean": meanp.tolist(),
        "initial_energy_variance": float(cov[0, 0]),
        "transported_energy_variance": float(covp[0, 0]),
    }

    # 6. A complete Stinespring instrument, with a positive-mass degenerate report.
    damping = 0.35
    ud = np.eye(4, dtype=complex)
    c, s = np.sqrt(1 - damping), np.sqrt(damping)
    ud[np.ix_([1, 2], [1, 2])] = np.array([[c, s], [-s, c]])
    ps = np.array([[3.0, 0.0, 0.0, 0.0], [5.0, 1.0, 0.0, 1.0]])
    report_p = np.array([2.0, 0.0, 0.0, 0.0])
    pb = np.repeat(ps, 2, axis=0) + report_p
    inst_m = 4
    _, qb, etb = bank(4, inst_m, pb, q0)
    vb, vb0, _ = rectangular_v(ud, inst_m)
    blank = np.zeros((4, 2), dtype=complex)
    blank[0, 0] = blank[2, 1] = 1
    inst = vb @ blank
    ideal_inst = vb0 @ blank
    it = inst.reshape(4, len(qb), 2)
    target = ud @ blank
    choi = np.einsum("ari,brj->aibj", it, it.conj()).reshape(8, 8) / 2
    tc = np.outer(target.reshape(-1), target.reshape(-1).conj()) / 2
    # Explicitly classicalize the report only for this instrument diagnostic.
    for a in range(4):
        for b in range(4):
            if a % 2 != b % 2:
                choi[2*a:2*a+2, 2*b:2*b+2] = 0
                tc[2*a:2*a+2, 2*b:2*b+2] = 0
    tp = np.trace(choi.reshape(4, 2, 4, 2), axis1=0, axis2=2)
    norm_res = maxabs(tp - np.eye(2) / 2)
    must_small(norm_res, "instrument normalization")
    err_inst = opnorm(inst - ideal_inst)
    must_le(err_inst, 2 * np.sqrt(2 / inst_m), "instrument full-output bound")
    report_one = float(np.sum(np.abs(it[[1, 3], :, 1]) ** 2))
    must_small(abs(report_one - damping), "instrument report probability")
    delta_choi = (choi - tc + (choi - tc).conj().T) / 2
    choi_half_trace = float(np.sum(np.abs(np.linalg.eigvalsh(delta_choi))) / 2)
    _, _, _, balances = moment_compressions(vb, pb, qb, etb)
    for val in balances:
        must_small(val, "instrument whole-material balance", tol=1e-9)
    groups["06_stinespring_instrument_all_outputs"] = {
        "status": "PASS", "expanded_system_dimension": 4, "M": inst_m,
        "reference_dimension": len(qb), "fresh_input_dimension": 2,
        "damping_parameter": damping, "report_one_on_excited_input": report_one,
        "normalization_residual": norm_res,
        "full_isometry_residual": maxabs(inst.conj().T @ inst - np.eye(2)),
        "full_output_operator_error": err_inst,
        "normalized_choi_half_trace_diagnostic": choi_half_trace,
        "choi_diagnostic_is_not_a_computed_diamond_norm": True,
        "all_input_momentum_balance_residuals": balances,
    }

    # 7. Two nodes with fresh banks; the first bank is retained and can be correlated.
    uh = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=complex) / np.sqrt(2)
    u2 = rotation(2, 0, 1, 0.41) @ np.diag(np.exp(1j * np.array([0.2, -0.3])))
    p2 = np.array([[3.0, 0.0, 0.0, 0.0], [5.0, 1.0, 0.5, -0.5]])
    q02 = np.array([9.0, 0.25, -0.125, 0.375])
    m1, m2 = 4, 5
    _, q1, e1 = bank(2, m1, p2, q0)
    _, q2, e2 = bank(2, m2, p2, q02)
    a1, i1, _ = rectangular_v(uh, m1)
    a2, i2, _ = rectangular_v(u2, m2)
    at1, at2 = a1.reshape(2, len(q1), 2), a2.reshape(2, len(q2), 2)
    combined = np.einsum("bqa,ari->brqi", at2, at1).reshape(2 * len(q1) * len(q2), 2)
    combined_ideal = np.einsum("ji,r,q->jrqi", u2 @ uh, e1, e2).reshape(combined.shape)
    err1, err2 = opnorm(a1 - i1), opnorm(a2 - i2)
    err_combined = opnorm(combined - combined_ideal)
    must_le(err_combined, err1 + err2, "complete two-step hybrid bound")
    iso2 = maxabs(combined.conj().T @ combined - np.eye(2))
    must_small(iso2, "two-step isometry")
    mean1, mean2 = np.abs(e1) ** 2 @ q1, np.abs(e2) ** 2 @ q2
    total2 = (p2[:, None, None, :] + q1[None, :, None, :] + q2[None, None, :, :]).reshape(-1, 4)
    bal2 = []
    for mu in range(4):
        actual = combined.conj().T @ (total2[:, mu, None] * combined)
        expected = np.diag(p2[:, mu]) + (mean1[mu] + mean2[mu]) * np.eye(2)
        bal2.append(maxabs(actual - expected))
        must_small(bal2[-1], "two-step full four-momentum", tol=1e-9)
    groups["07_two_step_fresh_reference_composition"] = {
        "status": "PASS", "M_values": [m1, m2],
        "reference_dimensions": [len(q1), len(q2)],
        "output_dimension_including_both_banks": len(combined),
        "first_bank_retained": True,
        "node_operator_errors": [err1, err2],
        "composed_operator_error": err_combined,
        "hybrid_bound_from_actual_node_operator_errors": err1 + err2,
        "proved_uniform_hybrid_bound": float(2*np.sqrt(2/m1) + 2*np.sqrt(2/m2)),
        "isometry_residual": iso2, "full_four_momentum_balance_residuals": bal2,
    }

    return {
        "schema": "research_1110_finite_box_four_momentum_v1",
        "status": "PASS", "groups_count": len(groups), "scientific_increment": 0,
        "scope": {
            "ideal_fucp_axioms_remain_exact": True,
            "finite_physical_implementations_are_approximate_on_original_system": True,
            "full_resources_and_arbitrary_passive_reference_in_proof": True,
            "four_momentum_conserved_exactly_by_finite_block_definition": True,
            "local_port_preparation_and_control_permissions_are_conditional": True,
            "does_not_close_round_1106_or_full_joint_model": True,
            "does_not_claim_optimal_resource_cost": True,
        },
        "groups": groups,
    }


def compare(actual, saved, path="root"):
    if isinstance(actual, dict):
        if not isinstance(saved, dict) or set(actual) != set(saved):
            raise AssertionError(f"{path}: dictionary keys differ")
        for key in actual:
            compare(actual[key], saved[key], path + "." + key)
    elif isinstance(actual, (list, tuple)):
        if not isinstance(saved, (list, tuple)) or len(actual) != len(saved):
            raise AssertionError(f"{path}: sequence differs")
        for i, (a, b) in enumerate(zip(actual, saved)):
            compare(a, b, f"{path}[{i}]")
    elif isinstance(actual, bool) or actual is None or isinstance(actual, str):
        if type(actual) is not type(saved) or actual != saved:
            raise AssertionError(f"{path}: exact value differs")
    elif isinstance(actual, (float, int)):
        if isinstance(saved, bool) or not isinstance(saved, (float, int)):
            raise AssertionError(f"{path}: numeric type differs")
        if not np.isfinite(actual) or not np.isfinite(saved):
            raise AssertionError(f"{path}: non-finite numeric value")
        if abs(actual - saved) > ATOL + RTOL * max(abs(actual), abs(saved)):
            raise AssertionError(f"{path}: {actual} != {saved}")
    else:
        raise TypeError(f"{path}: unsupported type {type(actual)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write this round's results.json")
    args = parser.parse_args()
    result = run()
    if args.write:
        RESULTS.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        action = "results.json written"
    else:
        if not RESULTS.exists():
            raise FileNotFoundError("results.json missing; only the explicit --write mode creates it")
        saved = json.loads(RESULTS.read_text(encoding="utf-8"))
        compare(result, saved)
        action = "saved results match (read-only)"
    for key, group in result["groups"].items():
        print(f"PASS {key}")
    print(f"PASS: {result['groups_count']} groups; {action}; scientific increment 0")


if __name__ == "__main__":
    main()