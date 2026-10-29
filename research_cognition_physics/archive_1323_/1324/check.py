"""1324: finite coherent reference, full joint evolution and source accounting.

Default / --check is strictly read-only. --write-new creates results exclusively.
No previous scientific implementation is imported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = STAGE.parent.parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_result(actual, saved, path="result"):
    """Verify every saved field recursively, including secondary diagnostics."""
    if isinstance(actual, dict):
        assert isinstance(saved, dict) and actual.keys() == saved.keys(), path
        for k in actual:
            same_result(actual[k], saved[k], path + "." + k)
    elif isinstance(actual, list):
        assert isinstance(saved, list) and len(actual) == len(saved), path
        for k, (a, b) in enumerate(zip(actual, saved)):
            same_result(a, b, f"{path}[{k}]")
    elif isinstance(actual, float):
        assert isinstance(saved, (int, float)) and not isinstance(saved, bool), path
        assert np.isfinite(saved) and np.isclose(
            actual, saved, rtol=2e-10, atol=2e-11
        ), path
    else:
        assert type(actual) is type(saved) and actual == saved, path


def opnorm(a):
    """Tall-column norm without forming a large output-space matrix."""
    a = np.asarray(a, complex)
    return float(np.sqrt(max(0.0, np.linalg.eigvalsh(a.conj().T @ a)[-1])))


def evolve(h, t):
    w, q = np.linalg.eigh(h)
    return (q * np.exp(-1j * t * w)) @ q.conj().T


def moments(L, nu, phase=0.0):
    k = np.arange(L + 1, dtype=float)
    logs = k * math.log(nu) - np.array([math.lgamma(j + 1) for j in range(L + 1)])
    logs -= logs.max()
    p = np.exp(logs)
    p /= p.sum()
    overlaps = np.sqrt(p[:-1] * p[1:])
    s = float(overlaps.sum())
    a = float((np.sqrt(k[1:]) * overlaps).sum())
    mu = float((k * p).sum())
    c = float((p / (k + 1)).sum())
    ds = float((((2 * k[1:] - 1) * overlaps).sum() - 2 * mu * s) / np.sqrt(nu))
    # The finite top occupation controls the ordinary annihilation residual.
    bare_eta2 = float(1 + nu * p[-1] * (1 - p[-1]))
    normalized_eta2 = float(max(1 - s * s, c))
    coeff = np.sqrt(p) * np.exp(1j * phase * k)
    return dict(p=p, coeff=coeff, s=s, a=a, mu=mu, c=c, ds=ds,
                eta2=normalized_eta2, bare_eta2=bare_eta2,
                top=float(p[-1]), vacuum=float(p[0]))


def full_parent(L, M, normalized=True, er=1.0, el=2.0, omega=1.0):
    """All one-fermion states in K=N_b+n_L <= L+1; no dark states removed."""
    n = len(M)
    basis = []
    for K in range(L + 2):
        for m in range(3):
            N = K - int(m != 0)
            if N >= 0:
                for n0 in range(N + 1):
                    for f in range(n):
                        basis.append((m, n0, N - n0, f))
    index = {x: j for j, x in enumerate(basis)}
    d = len(basis)
    V = np.zeros((d, d), complex)
    for (m, n0, n1, f), j in index.items():
        if m:
            continue
        N = n0 + n1
        if not N:
            continue
        for alpha, occ in enumerate((n0, n1)):
            if not occ:
                continue
            scale = np.sqrt(occ) / (np.sqrt(N) if normalized else 1)
            after = (n0 - int(alpha == 0), n1 - int(alpha == 1))
            for h in range(n):
                i = index[(alpha + 1, *after, h)]
                V[i, j] += scale * M[h, f]
                V[j, i] += scale * M[h, f].conjugate()
    nb = np.array([x[1] + x[2] for x in basis], float)
    nl = np.array([float(x[0] != 0) for x in basis])
    h0 = er * (1 - nl) + el * nl + omega * nb
    return basis, V, h0, nb, nl


def full_embedding(basis, L, nu, phase, u, n):
    """Embed every matter/flavor input with the same actual finite reference."""
    data = moments(L, nu, phase)
    J = np.zeros((len(basis), 3 * n), complex)
    for row, (m, n0, n1, f) in enumerate(basis):
        N = n0 + n1
        if N <= L:
            direction = np.sqrt(math.comb(N, n0)) * u[0] ** n0 * u[1] ** n1
            J[row, m * n + f] = data["coeff"][N] * direction
    return J, data


def effective(M, b):
    n = len(M)
    h = np.zeros((3 * n, 3 * n), complex)
    for alpha in range(2):
        lo = slice((alpha + 1) * n, (alpha + 2) * n)
        h[lo, :n] = b[alpha] * M
        h[:n, lo] = b[alpha].conjugate() * M.conj().T
    return h


class ReachableColumns:
    """Exact u=e0 invariant hull of all R,L0,L1 inputs, with n flavor suppressed.

    n=3,4,5 and M=g I_n are exact tensor products of these columns with I_n.
    Includes the highest right occupation L+1 and the transverse bright/dark block.
    """

    def __init__(self, L, g, normalized):
        self.L, self.g, self.normalized = L, g, normalized
        self.basis = [(0, 0, 0)]
        raw_blocks = []
        for K in range(1, L + 2):
            r = (0, K, 0)
            l = (1, K - 1, 0)
            self.basis.extend([r, l])
            raw_blocks.append((K, r, [(l, 1.0)]))
            r = (0, K - 1, 1)
            ls = [((2, K - 1, 0), 1 / np.sqrt(K))]
            if K > 1:
                ls.append(((1, K - 2, 1), np.sqrt((K - 1) / K)))
            self.basis.append(r)
            self.basis.extend(x[0] for x in ls)
            raw_blocks.append((K, r, ls))
        self.index = {x: j for j, x in enumerate(self.basis)}
        assert len(self.index) == len(self.basis)
        self.blocks = [
            (self.index[r], np.array([self.index[l] for l, w in ls]),
             np.array([w for l, w in ls]), g * (1 if normalized else np.sqrt(K)))
            for K, r, ls in raw_blocks
        ]
        self.nb = np.array([b[1] + b[2] for b in self.basis], float)
        self.nl = np.array([float(b[0] != 0) for b in self.basis])
        self.h0 = 1 + self.nb + self.nl

    def embed(self, coeff):
        J = np.zeros((len(self.basis), 3), complex)
        for N, z in enumerate(coeff):
            for m in range(3):
                J[self.index[(m, N, 0)], m] = z
        return J

    def apply(self, X):
        Y = np.zeros_like(X, dtype=complex)
        for r, ls, w, rate in self.blocks:
            Y[r] = rate * (w @ X[ls])
            Y[ls] += rate * w[:, None] * X[r]
        return Y

    def evolve(self, X, t):
        Y = X.copy()
        for r, ls, w, rate in self.blocks:
            c, s = np.cos(t * rate), np.sin(t * rate)
            bright = w @ X[ls]
            Y[r] = c * X[r] - 1j * s * bright
            Y[ls] = X[ls] + (c - 1) * w[:, None] * bright - 1j * s * w[:, None] * X[r]
        return Y

    def total_apply(self, X):
        return self.h0[:, None] * X + self.apply(X)


def mean_matrix(X, diagonal):
    return X.conj().T @ (diagonal[:, None] * X)


def reference_trace_distance(basis, state, coeff):
    """Trace distance of the actual reference from its initial pure state.

    At most three reference vectors occur; compute exact nonzero eigenvalues
    in their four-vector Gram span, rather than a large density matrix.
    """
    refs = sorted({(b[1], b[2]) for b in basis})
    pos = {x: i for i, x in enumerate(refs)}
    vectors = np.zeros((len(refs), 4), complex)
    for (m, n0, n1), z in zip(basis, state):
        vectors[pos[(n0, n1)], m] += z
    for N, z in enumerate(coeff):
        vectors[pos[(N, 0)], 3] = z
    G = vectors.conj().T @ vectors
    w, q = np.linalg.eigh(G)
    root = (q * np.sqrt(np.maximum(w, 0))) @ q.conj().T
    delta = root @ np.diag([1., 1., 1., -1.]) @ root
    return float(np.abs(np.linalg.eigvalsh(delta)).sum() / 2)


def compute():
    checks = []

    def close(name, actual, expected=0.0, tol=3e-9):
        residual = float(np.linalg.norm(np.asarray(actual) - np.asarray(expected)))
        assert residual <= tol, (name, residual, tol)
        checks.append(dict(name=name, kind="equality", residual=residual,
                           tolerance=tol, passed=True))

    def bound(name, value, upper, tol=3e-10):
        assert value <= upper + tol, (name, value, upper)
        checks.append(dict(name=name, kind="upper_bound", value=float(value),
                           upper=float(upper), tolerance=tol, passed=True))

    rng = np.random.default_rng(1324)
    # Non-scalar complex flavor matrix and a genuinely complex SU(2) direction.
    u = np.array([np.sqrt(.37), np.sqrt(.63) * np.exp(.61j)])
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
    r, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
    M = q @ np.diag([.04, .09, .13]) @ r.conj().T
    small_results = []
    for L, nu in ((2, .8), (4, 2.7)):
        for normed in (True, False):
            tag = f"full.L{L}.normalized{normed}"
            basis, V, h0, nb, nl = full_parent(L, M, normed)
            J, data = full_embedding(basis, L, nu, .37, u, 3)
            b = np.exp(.37j) * (data["s"] if normed else data["a"]) * u
            h = effective(M, b)
            eta = opnorm(M) * np.sqrt(data["eta2"] if normed else data["bare_eta2"])
            close(tag + ".input_isometry", J.conj().T @ J, np.eye(9))
            close(tag + ".hermitian", V, V.conj().T)
            close(tag + ".all_K_commute", (nb + nl)[:, None] * V, V * (nb + nl)[None, :])
            close(tag + ".resonant_H0", h0[:, None] * V, V * h0[None, :])
            close(tag + ".short_time_generator", J.conj().T @ V @ J, h)
            close(tag + ".full_input_residual", opnorm(V @ J - J @ h), eta)
            if normed:
                close(tag + ".parent_norm", opnorm(V), opnorm(M))
            t = 1.41
            U = evolve(V, t)
            X, Y = U @ J, J @ evolve(h, t)
            close(tag + ".actual_isometry", X.conj().T @ X, np.eye(9))
            bound(tag + ".joint_Duhamel", opnorm(X - Y), t * eta)
            H = np.diag(h0) + V
            close(tag + ".full_energy", X.conj().T @ H @ X, J.conj().T @ H @ J)
            close(tag + ".number_ledger", mean_matrix(X, nb + nl), mean_matrix(J, nb + nl))
            close(tag + ".free_energy", mean_matrix(X, h0), mean_matrix(J, h0))
            physical = evolve(H, t) @ J
            close(tag + ".physical_interaction_picture",
                  physical, np.exp(-1j * t * h0)[:, None] * X, 1e-8)
            partner = rng.normal(size=(9, 2)) + 1j * rng.normal(size=(9, 2))
            partner /= np.linalg.norm(partner)
            actual, ideal = X @ partner, Y @ partner
            bound(tag + ".unknown_partner_joint_vector",
                  np.linalg.norm(actual - ideal), t * eta)
            overlap = np.vdot(actual, ideal)
            distance = np.sqrt(max(0., 1 - abs(overlap) ** 2))
            bound(tag + ".unknown_partner_full_trace_distance", distance, min(1., t * eta))
            # Independent physical sources; no retuning along resonance.
            for key, D in (("er", 1 - nl), ("el", nl), ("omega", nb)):
                plus, minus = dict(er=1., el=2., omega=1.), dict(er=1., el=2., omega=1.)
                plus[key] += 1e-5
                minus[key] -= 1e-5
                hp = full_parent(L, M, normed, **plus)[2]
                hm = full_parent(L, M, normed, **minus)[2]
                close(tag + ".fixed_other_sources_" + key, (hp - hm) / 2e-5, D, 1e-8)
            small_results.append(dict(L=L, nu=nu, normalized=normed, dimension=len(basis),
                                      matter_flavor_dimension=9, eta=eta,
                                      joint_norm=opnorm(X - Y), time=t))

    # Cross-check the O(L) reachable construction against the entire occupation
    # Hamiltonian, including transverse L1 and its dark partner.
    for normed in (True, False):
        L, nu, g, t = 4, 2.1, .17, 2.3
        model = ReachableColumns(L, g, normed)
        data = moments(L, nu, .43)
        J = model.embed(data["coeff"])
        full_basis, V, _, _, _ = full_parent(L, np.array([[g]]), normed)
        embedding = np.zeros((len(full_basis), len(model.basis)), complex)
        pos = {x: i for i, x in enumerate(full_basis)}
        for j, x in enumerate(model.basis):
            embedding[pos[(*x, 0)], j] = 1
        Jf, _ = full_embedding(full_basis, L, nu, .43, np.array([1., 0.]), 1)
        tag = f"reachable.normalized{normed}"
        close(tag + ".initial_embedding", embedding @ J, Jf)
        close(tag + ".invariant_hull",
              V @ embedding, embedding @ model.apply(np.eye(len(model.basis))))
        close(tag + ".actual_columns", embedding @ model.evolve(J, t), evolve(V, t) @ Jf)

    L, nu0, nu1 = 400, 64., 256.
    phase = .29
    d0, d1 = moments(L, nu0, phase), moments(L, nu1, phase)
    calibrated_amplitude = .1
    gnorm = calibrated_amplitude / d0["s"]
    gbare = calibrated_amplitude / d0["a"]
    time = np.pi / (2 * calibrated_amplitude)
    numerical_allowance = 1e-9
    witness = []
    response = []
    for nu, data in ((nu0, d0), (nu1, d1)):
        k = np.arange(L + 1, dtype=float)
        close(f"moments.nu{nu}.a_top", data["a"], np.sqrt(nu) * (1 - data["top"]))
        close(f"moments.nu{nu}.mu_top", data["mu"], nu * (1 - data["top"]))
        close(f"moments.nu{nu}.bare_residual_identity",
              1 + data["mu"] - data["a"] ** 2, data["bare_eta2"])
        bound(f"moments.nu{nu}.downward_variance_nonnegative",
              -(1 - data["vacuum"] - data["s"] ** 2), 0)
        eps_deriv = 1e-3
        rr = np.sqrt(nu)
        finite_deriv = (moments(L, (rr + eps_deriv) ** 2)["s"]
                       - moments(L, (rr - eps_deriv) ** 2)["s"]) / (2 * eps_deriv)
        close(f"moments.nu{nu}.preparation_derivative", finite_deriv, data["ds"], 2e-9)
        tail_log = -nu + (L + 1) * (1 + np.log(nu / (L + 1)))
        response.append(dict(nu=nu, s=data["s"], a=data["a"], c=data["c"],
                             mean_number=data["mu"], top_probability=data["top"],
                             vacuum_probability=data["vacuum"],
                             normalized_preparation_derivative=data["ds"],
                             poisson_chernoff_tail_upper=float(np.exp(tail_log))))
        for normed, g in ((True, gnorm), (False, gbare)):
            tag = f"large.nu{nu}.normalized{normed}"
            model = ReachableColumns(L, g, normed)
            J = model.embed(data["coeff"])
            amp = data["s"] if normed else data["a"]
            h = effective(np.array([[g]]), np.exp(1j * phase) * amp * np.array([1., 0.]))
            eta = g * np.sqrt(data["eta2"] if normed else data["bare_eta2"])
            residual = model.apply(J) - J @ h
            close(tag + ".all_input_isometry", J.conj().T @ J, np.eye(3))
            close(tag + ".calibration_generator", J.conj().T @ model.apply(J), h)
            close(tag + ".full_residual_eta", opnorm(residual), eta)
            X, Y = model.evolve(J, time), J @ evolve(h, time)
            close(tag + ".actual_all_input_isometry", X.conj().T @ X, np.eye(3))
            joint_norm = opnorm(X - Y)
            bound(tag + ".full_joint_Duhamel", joint_norm, time * eta)
            close(tag + ".exact_number_ledger",
                  mean_matrix(X, model.nb + model.nl),
                  mean_matrix(J, model.nb + model.nl), 1e-8)
            close(tag + ".free_energy_ledger",
                  mean_matrix(X, model.h0), mean_matrix(J, model.h0), 1e-8)
            close(tag + ".interaction_energy_ledger",
                  X.conj().T @ model.apply(X), J.conj().T @ model.apply(J))
            close(tag + ".total_energy_ledger",
                  X.conj().T @ model.total_apply(X), J.conj().T @ model.total_apply(J), 1e-8)
            prob = float(mean_matrix(X, model.nl)[0, 0].real)
            prob_eff = float(mean_matrix(Y, model.nl)[0, 0].real)
            direct = ((1 - data["vacuum"]) * np.sin(g * time) ** 2 if normed
                      else float(np.dot(data["p"], np.sin(g * time * np.sqrt(k)) ** 2)))
            close(tag + ".actual_probability_from_full_columns", prob, direct)
            close(tag + ".effective_probability", prob_eff, np.sin(g * amp * time) ** 2)
            bound(tag + ".probability_difference", abs(prob - prob_eff), time * eta)
            true_change = float(mean_matrix(X, model.nb)[0, 0].real - data["mu"])
            close(tag + ".actual_reference_consumption", true_change, -prob, 1e-8)
            bound(tag + ".bounded_ledger_reconstruction", abs(true_change + prob_eff), time * eta)
            # The product reference has unchanged number while the actual one
            # loses a quantum. This is recorded, not hidden by trace-distance.
            close(tag + ".comparison_reference_mean",
                  mean_matrix(Y, model.nb), data["mu"] * np.eye(3), 1e-8)
            actual_ref_distance = reference_trace_distance(model.basis, X[:, 0], data["coeff"])
            bound(tag + ".reference_poststate_distance", actual_ref_distance, time * eta)
            full_trace = float(np.sqrt(max(0., 1 - abs(np.vdot(X[:, 0], Y[:, 0])) ** 2)))
            bound(tag + ".actual_complete_poststate_distance", full_trace, min(1., time * eta))
            witness.append(dict(nu=nu, normalized=normed, coupling=g, time=float(time),
                                reachable_dimension=len(model.basis), eta=float(eta),
                                joint_norm=joint_norm, joint_error_bound=float(min(1., time * eta)),
                                actual_probability=prob, effective_probability=prob_eff,
                                actual_reference_mean_change=true_change,
                                actual_reference_trace_distance=actual_ref_distance,
                                full_joint_pure_trace_distance=full_trace,
                                effective_interval=[max(0., prob_eff - time * eta - numerical_allowance),
                                                    min(1., prob_eff + time * eta + numerical_allowance)]))

    held_norm = next(x for x in witness if x["nu"] == nu1 and x["normalized"])
    held_bare = next(x for x in witness if x["nu"] == nu1 and not x["normalized"])
    gap = held_norm["effective_interval"][0] - held_bare["effective_interval"][1]
    bound("heldout.strict_interval_separation", .70, gap)
    close("calibration.normalized", gnorm * d0["s"], calibrated_amplitude)
    close("calibration.ordinary", gbare * d0["a"], calibrated_amplitude)

    # Phase averaging removes the short-time reduced generator, but does not
    # remove the parent conversion. Calculate the mixed reference explicitly.
    Lp, nup, gp, tp = 3, 1.7, .23, 2.1
    basis, V, _, nb, nl = full_parent(Lp, np.array([[gp]]), True)
    U = evolve(V, tp)
    pos = {x: i for i, x in enumerate(basis)}
    data = moments(Lp, nup)
    hmix = np.zeros((3, 3), complex)
    rate_mix = 0.
    rho_matter = np.zeros((3, 3), complex)
    refs = sorted({(x[1], x[2]) for x in basis})
    for N, weight in enumerate(data["p"]):
        JN = np.zeros((len(basis), 3), complex)
        for m in range(3):
            JN[pos[(m, N, 0, 0)], m] = 1
        hmix += weight * JN.conj().T @ V @ JN
        state = U @ JN[:, 0]
        rate_mix += weight * float(np.vdot(state, nl * state).real)
        for n0, n1 in refs:
            vec = np.array([state[pos[(m, n0, n1, 0)]] if (m, n0, n1, 0) in pos else 0
                            for m in range(3)])
            rho_matter += weight * np.outer(vec, vec.conj())
    close("phase_averaged.short_time_generator_zero", hmix)
    close("phase_averaged.nonzero_actual_transfer",
          rate_mix, (1 - data["vacuum"]) * np.sin(gp * tp) ** 2)
    bound("phase_averaged.transfer_positive", .1, rate_mix)
    close("phase_averaged.no_reduced_offdiagonal", rho_matter - np.diag(np.diag(rho_matter)))

    # The conversion/mixing preparation response is one common scalar multiplier,
    # not a new choice of M. Check a non-scalar complex matrix.
    ratio = d1["s"] / d0["s"]
    close("matrix_response.same_left_gram",
          (d1["s"] * M) @ (d1["s"] * M).conj().T,
          ratio ** 2 * (d0["s"] * M) @ (d0["s"] * M).conj().T)
    close("matrix_response.same_right_gram",
          (d1["s"] * M).conj().T @ (d1["s"] * M),
          ratio ** 2 * (d0["s"] * M).conj().T @ (d0["s"] * M))
    close("matrix_response.every_singular_value",
          np.linalg.svd(d1["s"] * M, compute_uv=False),
          ratio * np.linalg.svd(d0["s"] * M, compute_uv=False))

    return dict(
        round=1324, status="PASS",
        scope="finite coherent preparation; full joint matter-reference-partner error and exact source ledger",
        new_axioms=0, checks=checks, count=len(checks),
        max_equality_residual=max(x["residual"] for x in checks if x["kind"] == "equality"),
        small_full_parent=small_results, finite_moments=response,
        calibration=dict(kind="exact short-time reduced generator h=J_dagger V J",
                         nu=nu0, normalized_coupling=gnorm, ordinary_coupling=gbare,
                         amplitude=calibrated_amplitude, L=L, time=float(time),
                         no_heldout_refitting=True),
        large_joint_witness=witness,
        heldout_probability_interval_gap=gap,
        numerical_interval_allowance=numerical_allowance,
        phase_averaged_reference=dict(L=Lp, nu=nup, coupling=gp, time=tp,
                                      actual_probability=rate_mix,
                                      short_time_generator_norm=opnorm(hmix)),
        common_preparation_multiplier=ratio,
        flavor_extension=dict(families=[3, 4, 5],
                              method="exact tensor identity for M=g I_n; no additional numerical family tests",
                              general_complex_flavor_test_n=3),
        limitations=[
            "coherent preparation and internal phase resources are declared inputs",
            "general-N theorem reused; implementation contract extended only to fixed K<=401",
            "normalized and ordinary candidates do not both satisfy the 1323 candidate balance contract",
            "reference number/energy uses the exact ledger, never an unbounded-observable trace bound",
            "radial derivative is preparation response, not a Higgs functional source or particle pole",
            "no complete dynamic gauge/Higgs device or UV completion certified",
            "finite probability predictions do not give free preparation, clock or readout permissions",
        ],
    )


def run():
    result = compute()
    paths = [HERE / "proof.md", HERE / "check.py",
             STAGE / "_admission/round1324_admission.md",
             STAGE / "1323/proof.md", STAGE / "1323/check.py", STAGE / "1323/results.json"]
    result["source_sha256"] = {
        p.relative_to(ROOT).as_posix(): digest(p) for p in paths
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write_new:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        same_result(result, json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(dict(round=1324, mode="write-new" if args.write_new else "read-only",
                          passed=True, checks=result["count"],
                          max_equality_residual=result["max_equality_residual"],
                          heldout_gap=result["heldout_probability_interval_gap"])))


if __name__ == "__main__":
    main()
