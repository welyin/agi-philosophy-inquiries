"""1213 author checks. Default replays read-only; --write creates absent results.

Continuous-spectrum and domain statements are proved in proof.md. The matrix
check is restricted to explicitly supported input/paired packet sectors and
does not identify a finite truncated clock with the full half-line clock.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def frac(x):
    return str(x)


def entropy(rho):
    vals = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    require(float(vals.min()) > -2e-12, "non-positive density matrix")
    vals = vals[vals > 1e-14]
    return float(-np.dot(vals, np.log(vals)))


def marginal(rho, dims, keep):
    a = rho.reshape(tuple(dims) * 2)
    remaining = list(range(len(dims)))
    for original in sorted(set(remaining) - set(keep), reverse=True):
        ix = remaining.index(original)
        a = np.trace(a, axis1=ix, axis2=ix + len(remaining))
        remaining.pop(ix)
    d = math.prod(dims[i] for i in keep)
    return a.reshape(d, d)


def compute():
    source_paths = [
        "archive_1201_/1202/results.json",
        "archive_1201_/1205/proof.md",
        "archive_1201_/1209/proof.md",
        "archive_1201_/1211/proof.md",
        "archive_1201_/1211/results.json",
    ]
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
              for p in source_paths}
    hashes["archive_1201_/1213/check.py"] = hashlib.sha256(
        Path(__file__).read_bytes()).hexdigest()
    old = json.loads((ROOT / source_paths[0]).read_text(encoding="utf-8"))
    prep = old["preparation"]
    den = int(prep["root_denominator"])
    boxes = prep["root_boxes"]
    # Accept the exact inherited schema without recomputing the spectrum.
    if isinstance(boxes, dict):
        boxes = list(boxes.values())
    box = [(F(int(x[0]), den), F(int(x[1]), den)) for x in boxes]
    p_exact = [F(x) for x in old["one_marginal_rational"]]
    require(sum(p_exact) == 1 and all(x > 0 for x in p_exact), "native r")
    gap_lo = [F(0)] + [box[i][0] - box[0][1] for i in range(1, 4)]
    gap_hi = [F(0)] + [box[i][1] - box[0][0] for i in range(1, 4)]
    gaps = np.array([float((a + b) / 2) for a, b in zip(gap_lo, gap_hi)])
    require(gap_hi[3] < F(547, 500), "width bound")
    require(gap_lo[1] > 0, "positive first gap")
    p = np.array([float(x) for x in p_exact])
    groups = {}
    groups["native_material_and_resource_scope"] = {
        "pass": True,
        "material": "unchanged h, energy basis; no same-h ground inventory",
        "root_boxes": [[frac(a), frac(b)] for a, b in box],
        "gap_bounds": [[frac(a), frac(b)] for a, b in zip(gap_lo, gap_hi)],
        "clock": "L2(R,dp), P^2/(2M); infinite dimension; pure normal Gaussian",
        "record": "four-dimensional, zero Hamiltonian, independent pure blank",
        "operation": "new global W/G, not native QQ or established local-clock menu",
    }

    # Exact integration for normalized constant energy packets; no quadrature.
    a, b, shift = F(1, 4), F(3, 4), F(1, 2)
    mean = (a + b) / 2
    second = (a*a + a*b + b*b) / 3
    shifted_mean = (a + b + 2*shift) / 2
    shifted_second = ((a+shift)**2 + (a+shift)*(b+shift) + (b+shift)**2) / 3
    tt_dagger_prob = max(F(0), b - max(a, shift)) / (b-a)
    require(shifted_mean == mean + shift, "packet first moment")
    require(shifted_second == second + 2*shift*mean + shift**2, "packet second moment")
    require(tt_dagger_prob == F(1, 2), "threshold projection")
    groups["half_line_isometry_packet_identities"] = {
        "pass": True, "norm_Tf_squared": "1", "norm_Tdagger_Tf_squared": "1",
        "norm_T_Tdagger_f_squared": frac(tt_dagger_prob),
        "mean_before_after": [frac(mean), frac(shifted_mean)],
        "second_before_after": [frac(second), frac(shifted_second)],
        "scope": "exact packet-integral diagnostics; general isometry/domain proof is analytic",
    }

    # Clock packets: normalized intervals [base+gap_i,base+gap_i+ell].
    # All eight are disjoint. On the 14 included fixed/paired sectors,
    # Q = diag(base+gap_system+gap_packet) + u*I, 0<=u<=ell.
    # W preserves the same internal u exactly. Its full 128-dimensional
    # identity completion here is only a diagnostic surrogate off these sectors.
    bases, ell = [F(2), F(4)], F(1, 1000)
    for i in range(3):
        require(gap_hi[i] + ell < gap_lo[i+1], "disjoint shifted packets")
    dc, dk, ds = 8, 4, 4
    dim = ds*dc*dk

    def index(s, c, k):
        return (s*dc + c)*dk + k

    def formal_energy(s, c):
        packet_i = c % 4
        z = [0, 0, 0]
        if s:
            z[s-1] += 1
        if packet_i:
            z[packet_i-1] += 1
        return (bases[c//4], *z)

    perm = np.arange(dim)
    pairs = []
    supported = []
    for base_ix in range(2):
        supported.append(index(0, base_ix*4, 0))
        for i in range(1, 4):
            x, y = index(i, base_ix*4, 0), index(0, base_ix*4+i, i)
            require(formal_energy(i, base_ix*4) == formal_energy(0, base_ix*4+i),
                    "formal total energy mismatch")
            perm[x], perm[y] = y, x
            pairs.append([x, y])
            supported.extend([x, y])
    require(len(set(supported)) == 14, "supported-sector dimension")
    require(np.array_equal(perm[perm], np.arange(dim)), "involution")
    W = np.eye(dim)[perm]
    ec = np.array([float(b0 + ell/2) + gaps[i] for b0 in bases for i in range(4)])
    hs = np.repeat(gaps, dc*dk)
    hc = np.tile(np.repeat(ec, dk), ds)
    q = hs + hc
    q_error = float(np.max(np.abs(q[perm] - q)))
    require(q_error < 3e-15, "numerical compressed energy mismatch")
    groups["supported_swap_and_exact_energy_labels"] = {
        "pass": True, "diagnostic_matrix_dimension": dim,
        "actual_included_fixed_or_paired_sectors": 14, "swap_pairs": pairs,
        "formal_energy_label_equalities": 6,
        "packet_internal_energy_u_is_unchanged": True,
        "compressed_energy_residual": q_error,
        "scope": "coincides with continuous W on the 14 supported sectors only; no cutoff theorem",
    }

    # All-input code isometry, then one reproducible unknown-SR diagnostic.
    ccoeff = np.array([math.sqrt(1/3), 1j*math.sqrt(2/3)])
    A = np.zeros((dim, ds), dtype=complex)
    for i in range(4):
        for b0 in range(2):
            A[index(i, 4*b0, 0), i] = ccoeff[b0]
    J = W @ A
    target_excited = np.array([ix for ix in range(dim) if ix // (dc*dk) > 0])
    isometry_error = float(np.linalg.norm(J.conj().T @ J - np.eye(ds), 2))
    reset_error = float(np.linalg.norm(J[target_excited], 2))
    require(isometry_error < 1e-14 and reset_error == 0, "code reset/isometry")
    rng = np.random.default_rng(1213)
    phi = rng.normal(size=(4, 3)) + 1j*rng.normal(size=(4, 3))
    phi /= np.linalg.norm(phi)
    output = J @ phi
    recovered = J.conj().T @ output
    r_before = phi.T @ phi.conj()
    r_after = output.T @ output.conj()
    recovery_error = float(np.linalg.norm(recovered - phi))
    r_error = float(np.linalg.norm(r_before - r_after))
    require(max(recovery_error, r_error) < 1e-14, "unknown reference recovery")
    groups["unknown_SR_and_complete_record"] = {
        "pass": True, "reference_dimension": 3, "rng_seed": 1213,
        "code_isometry_residual": isometry_error, "excited_output_residual": reset_error,
        "joint_SR_recovery_residual": recovery_error, "R_marginal_residual": r_error,
        "scope": "analytic formula covers every SR; numerical sample and code matrix cross-check it",
    }

    tau = 1.0
    G = math.pi / (2*tau) * (np.eye(dim) - W)
    H = np.diag(q) + G
    evals, evecs = np.linalg.eigh(H)
    U = (evecs * np.exp(-1j*tau*evals)) @ evecs.conj().T
    ideal_U = np.exp(-1j*tau*q)[:, None] * W
    endpoint_error = float(np.linalg.norm(U @ A - ideal_U @ A, 2))
    min_G = float(np.linalg.eigvalsh(G).min())
    require(min_G > -1e-12 and float(evals.min()) >= 0, "positive parent")
    require(endpoint_error < 2e-13, "actual endpoint phase")
    raw_W_error = float(np.linalg.norm(U @ A - W @ A, 2))
    require(raw_W_error > .1, "free phases cannot be erased")
    # Complete revival of populations at 2*tau, hence no released/stable product.
    U2 = (evecs * np.exp(-2j*tau*evals)) @ evecs.conj().T
    revival_error = float(np.linalg.norm(U2 @ A - np.exp(-2j*tau*q)[:, None]*A, 2))
    require(revival_error < 2e-13, "revival")
    groups["positive_parent_endpoint_and_nonrelease"] = {
        "pass": True, "tau": tau, "G_norm": math.pi/tau,
        "min_parent_eigenvalue_diagnostic": float(evals.min()),
        "endpoint_error_with_free_phase": endpoint_error,
        "endpoint_error_if_free_phase_omitted": raw_W_error,
        "two_tau_free_revival_error": revival_error,
        "interaction_switched_off": False, "released_or_persistent_cold_output_certified": False,
    }

    M, pc, s = F(1), F(2), F(1, 2)
    mu = (pc**2 + s**2)/(2*M)
    nu = (pc**4 + 6*pc**2*s**2 + 3*s**4)/(4*M**2)
    lo = sum(pi*gi for pi, gi in zip(p_exact, gap_lo))
    hi = sum(pi*gi for pi, gi in zip(p_exact, gap_hi))
    lo2 = sum(pi*gi**2 for pi, gi in zip(p_exact, gap_lo))
    hi2 = sum(pi*gi**2 for pi, gi in zip(p_exact, gap_hi))
    clock2lo, clock2hi = nu + 2*mu*lo + lo2, nu + 2*mu*hi + hi2
    require(0 < lo <= hi < gap_hi[3], "native energy interval")
    require(clock2hi <= nu+2*mu*gap_hi[3]+gap_hi[3]**2, "second moment")
    rho_in = A @ np.diag(p) @ A.conj().T
    rho_out = ideal_U @ rho_in @ ideal_U.conj().T
    e_change = float(np.trace(np.diag(hc) @ (rho_out-rho_in)).real)
    g_in = float(np.trace(G @ rho_in).real)
    g_out = float(np.trace(G @ rho_out).real)
    energy_error = abs(e_change - float(np.dot(p, gaps)))
    record_out = marginal(rho_out, [4, 8, 4], [2])
    record_error = float(np.linalg.norm(record_out - np.diag(p)))
    joint_ck_entropy = entropy(marginal(rho_out, [4, 8, 4], [1, 2]))
    native_entropy = float(-np.dot(p, np.log(p)))
    require(max(energy_error, abs(g_in-g_out), abs(g_in-math.pi*(1-p[0])/2),
                record_error, abs(joint_ck_entropy-native_entropy)) < 3e-13, "full ledger")
    groups["finite_moments_native_energy_and_entropy_ledger"] = {
        "pass": True, "Gaussian_parameters_M_pc_s": [frac(M), frac(pc), frac(s)],
        "Gaussian_HC_mean": frac(mu), "Gaussian_HC_second": frac(nu),
        "native_energy_transferred_to_clock_interval": [frac(lo), frac(hi)],
        "native_energy_transferred_diagnostic": float(np.dot(p, gaps)),
        "native_Gaussian_clock_second_out_interval": [frac(clock2lo), frac(clock2hi)],
        "native_Gaussian_clock_second_out_diagnostic": float((clock2lo+clock2hi)/2),
        "full_G_mean_in_out": [g_in, g_out], "record_entropy_native": native_entropy,
        "CK_entropy_native": joint_ck_entropy, "clock_energy_residual": energy_error,
        "pure_K_returned": False, "clock_returned": False,
    }

    cdf = []
    for gap in gaps[1:]:
        bound = math.sqrt(2*float(M)*gap)
        prob = .5*(math.erf((bound-float(pc))/(math.sqrt(2)*float(s)))
                   - math.erf((-bound-float(pc))/(math.sqrt(2)*float(s))))
        cdf.append(prob)
    mixed_F = .25 + p[0]/4*sum(cdf)
    upper_F = F(1, 4) + F(3, 4)*p_exact[0]
    require(mixed_F <= float(upper_F) < 1, "mixed blank deletion")
    # Scalar upward shift on a finite ladder cannot be an isometry on all inputs.
    truncated_T = np.zeros((7, 7))
    for k in range(6):
        truncated_T[k+1, k] = 1
    defect = np.eye(7) - truncated_T.T @ truncated_T
    require(np.array_equal(np.diag(defect), [0, 0, 0, 0, 0, 0, 1]), "ceiling defect")
    groups["resource_deletion_and_finite_ceiling"] = {
        "pass": True, "same_W_mixed_K_ground_probability": mixed_F,
        "Gaussian_energy_below_native_gaps": cdf,
        "mixed_K_ground_probability_exact_upper": frac(upper_F),
        "no_energy_shift_loss_for_excited_i": "Delta_i > 0 (formal exact)",
        "finite_ladder_TdaggerT_top_defect": 1,
        "Gaussian_shift_without_K_isometry": False,
        "scope": "same-W deletion, not a universal ban on other pure-clock protocols",
    }

    # The derivative squared behaves as x^(-5/2) near each shifted threshold.
    # This exact exponent is an analytic non-domain witness, not a grid test.
    threshold_power = F(-1, 4)
    derivative_squared_power = 2*(threshold_power-1)
    require(2*threshold_power > -1 and derivative_squared_power < -1,
            "normal state but no finite X2")
    groups["threshold_domain_and_claim_boundaries"] = {
        "pass": True, "threshold_wavefunction_power": frac(threshold_power),
        "threshold_probability_power": frac(2*threshold_power),
        "derivative_squared_power": frac(derivative_squared_power),
        "shifted_Gaussian_normal_L2": True,
        "shifted_Gaussian_in_DX_for_positive_gap": False,
        "continuous_spectrum_proof_replaced_by_matrix": False,
        "old_QQ_or_local_SR_permission_proved": False,
        "all_resources_finite_claim": False,
        "third_law_or_Masanes_Oppenheim_refuted": False,
        "weak_energy_only_endpoint_unattainability_refuted": True,
    }
    require(len(groups) == 8 and all(g["pass"] for g in groups.values()), "groups")
    return {
        "round": 1213, "status": "AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE",
        "all_checks_passed": True, "check_groups": 8,
        "source_sha256": hashes, "groups": groups,
        "new_adopted_axioms": 0,
        "scientific_proof": "proof.md; continuous spectrum/strong commutation are analytic",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    target = HERE / "results.json"
    if args.write and target.exists():
        raise SystemExit("Refusing --write: results.json already exists; use read-only replay.")
    computed = compute()
    if args.write:
        with target.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(computed, indent=2, ensure_ascii=False, sort_keys=True)+"\n")
    stored = json.loads(target.read_text(encoding="utf-8"))
    match = stored == computed
    require(match, "results.json differs; default never overwrites results")
    print(json.dumps({"status": "PASS", "all_checks_passed": True,
                      "results_match": match, "check_groups": 8,
                      "mode": "write_and_replay" if args.write else "read_only_replay"},
                     ensure_ascii=False))


if __name__ == "__main__":
    main()
