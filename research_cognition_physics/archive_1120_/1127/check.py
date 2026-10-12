"""1127: finite algebra diagnostics, not a continuum Dirac-vacuum cutoff.

Default: recompute, assert discriminating checks, read/compare results.json.
Only --write-results registers results; no other paths are written.
Requires Python 3 and NumPy. All matrices retain every CAR mode.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOL = 2e-10


def norm(a):
    return float(np.linalg.norm(a))


def expi(a):
    """exp(-i a), for a Hermitian finite matrix."""
    values, vectors = np.linalg.eigh(a)
    return (vectors * np.exp(-1j * values)) @ vectors.conj().T


def unitary(rng, n):
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    return q @ np.diag(np.diag(r) / np.abs(np.diag(r)))


def annihilators(n):
    operators = []
    for j in range(n):
        a = np.zeros((2**n, 2**n), complex)
        for bits in range(2**n):
            if (bits >> j) & 1:
                sign = (-1) ** ((bits & ((1 << j) - 1)).bit_count())
                a[bits ^ (1 << j), bits] = sign
        operators.append(a)
    return operators


def creation(v, ann):
    return sum(v[j] * ann[j].conj().T for j in range(len(ann)))


def exterior_lift(u):
    """Full fermionic exterior algebra, all 2**n states."""
    n = len(u)
    g = np.zeros((2**n, 2**n), complex)
    g[0, 0] = 1
    for k in range(1, n + 1):
        subsets = list(itertools.combinations(range(n), k))
        for left in subsets:
            i = sum(1 << j for j in left)
            for right in subsets:
                j = sum(1 << q for q in right)
                g[i, j] = np.linalg.det(u[np.ix_(left, right)])
    return g


def covariance(state, ann):
    images = [a @ state for a in ann]
    # C_ij = <c_j^dagger c_i>, consistent with C -> U C U^dagger.
    return np.array([[np.vdot(images[j], images[i])
                      for j in range(len(ann))] for i in range(len(ann))])


def factor_target(state, n=6, target=(2, 3)):
    """CAR permutation including its parity, then target x remainder."""
    rest = tuple(j for j in range(n) if j not in target)
    order = target + rest
    result = np.zeros((2**len(target), 2**len(rest)), complex)
    for bits in range(2**n):
        occupied = [j for j in order if (bits >> j) & 1]
        inversions = sum(a > b for i, a in enumerate(occupied)
                         for b in occupied[i + 1:])
        t = sum(((bits >> j) & 1) << k for k, j in enumerate(target))
        r = sum(((bits >> j) & 1) << k for k, j in enumerate(rest))
        result[t, r] = (-1) ** inversions * state[bits]
    return result


def fixture():
    rng = np.random.default_rng(1127)
    n = 6
    ann = annihilators(n)
    q = unitary(rng, n)
    spectrum = np.array([-2.7, -1.9, -1.1, 0.8, 1.6, 3.2])
    h = (q * spectrum) @ q.conj().T
    s = q[:, :3] @ q[:, :3].conj().T
    vacuum = np.zeros(2**n, complex)
    vacuum[0] = 1
    for v in q[:, :3].T[::-1]:
        vacuum = creation(v, ann) @ vacuum
    f = np.eye(n, dtype=complex)[:, :2]
    empty_mask = np.array([float((i & 3) == 0) for i in range(2**n)])
    base_raw = empty_mask * vacuum
    probability = float(np.vdot(base_raw, base_raw).real)
    base = base_raw / np.sqrt(probability)
    j = np.column_stack([creation(f[:, i], ann) @ base for i in range(2)])
    swap = np.eye(n, dtype=complex)
    for a, b in ((0, 2), (1, 3)):
        swap[a, a] = swap[b, b] = 0
        swap[a, b] = swap[b, a] = -1j
    extra = np.eye(n, dtype=complex)
    extra[2:, 2:] = unitary(rng, 4)
    u = swap @ extra
    lift = exterior_lift(u)
    return dict(rng=rng, ann=ann, q=q, h=h, spectrum=spectrum, s=s,
                vacuum=vacuum, f=f, base=base, probability=probability,
                j=j, swap=swap, u=u, lift=lift)


def group_car(d):
    ann, u, lift, j, base = (d[k] for k in ("ann", "u", "lift", "j", "base"))
    identity = np.eye(64)
    car_error = max(norm(a @ b.conj().T + b.conj().T @ a -
                         (identity if i == k else 0))
                    for i, a in enumerate(ann) for k, b in enumerate(ann))
    lift_error = max(norm(lift @ a.conj().T @ lift.conj().T -
                          creation(u[:, i], ann)) for i, a in enumerate(ann))
    out = lift @ j
    xi = lift @ base
    target_creators = [ann[k].conj().T for k in (2, 3)]
    columns_error = max(norm(out[:, i] + 1j * target_creators[i] @ xi)
                        for i in range(2))
    units_error = 0.0
    for i in range(2):
        for k in range(2):
            expected = np.zeros((2, 2), complex)
            expected[i, k] = 1
            observed = out.conj().T @ target_creators[i] @ ann[k + 2] @ out
            units_error = max(units_error, norm(observed - expected))
    tensors = np.stack([factor_target(out[:, i]) for i in range(2)])
    coeff = np.array([[1, 0.2j, -0.3], [0.1 - 0.4j, 0.7, 0.6j]], complex)
    coeff /= norm(coeff)
    psi = np.einsum("qte,qr->ter", tensors, coeff)
    rho_target_r = np.einsum("ter,ues->trus", psi, psi.conj()).reshape(12, 12)
    ideal = np.zeros((4, 3), complex)
    ideal[1, :] = coeff[0]
    ideal[2, :] = coeff[1]
    qr_error = norm(rho_target_r - np.outer(ideal.ravel(), ideal.ravel().conj()))
    environment_error = norm(tensors[0].conj().T @ tensors[0] -
                             tensors[1].conj().T @ tensors[1])
    environment_cross = norm(tensors[0].T @ tensors[1].conj())
    wrong_out = exterior_lift(d["swap"]) @ j
    # Both gates hand off the code, but they differ on the retained environment.
    full_trace_distance = np.sqrt(max(0.0, 1 - abs(np.vdot(out[:, 0], wrong_out[:, 0]))**2))
    correlation = norm(d["s"][:2, 2:])
    errors = [car_error, lift_error, norm(j.conj().T @ j - np.eye(2)),
              columns_error, units_error, qr_error, environment_error, environment_cross]
    assert max(errors) < TOL
    assert correlation > 0.1 and full_trace_distance > 0.1
    assert norm(d["s"] @ u - u @ d["s"]) > 0.1
    return dict(modes=6, fock_dimension=64, reference_dimension=3,
                car_error=car_error, lift_error=lift_error,
                source_orthogonality_error=errors[2], column_error=columns_error,
                all_matrix_units_error=units_error, qr_error=qr_error,
                environment_independence_error=environment_error,
                environment_cross_error=environment_cross,
                vacuum_source_remainder_correlation=correlation,
                deleting_remainder_full_trace_distance=float(full_trace_distance),
                source_empty_probability=d["probability"])


def group_dirac():
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1, -1]).astype(complex)
    alpha = [np.kron(sx, z) for z in (sx, sy, sz)]
    beta = np.kron(sz, np.eye(2))
    momenta = [(0.7, -0.2, 1.1), (-0.4, 0.9, 0.3), (1.5, 0.1, -0.8)]
    durations = [0.13, 0.29, 0.37, 0.21]
    areas = (np.pi / 2) * np.array([0.11, 0.34, 0.39, 0.16])
    m = 0.8
    seed = np.array([1, 0.2j, -0.4, 0.3 + 0.1j], complex)
    seed /= norm(seed)
    exact_error = 0.0
    beta_difference = []
    negative_weights = []
    for p in momenta:
        h = sum(x * a for x, a in zip(p, alpha)) + m * beta
        hfull = np.kron(h, np.eye(2))
        v = np.kron(np.eye(4), sx)
        wrong_v = np.kron(beta, sx)
        u = np.eye(8, dtype=complex)
        wrong = u.copy()
        for dt, angle in zip(durations, areas):
            u = expi(dt * hfull + angle * v) @ u
            wrong = expi(dt * hfull + angle * wrong_v) @ wrong
        source = np.kron(seed, np.array([1, 0]))
        target = -1j * np.kron(expi(sum(durations) * h) @ seed, np.array([0, 1]))
        exact_error = max(exact_error, norm(u @ source - target))
        beta_difference.append(norm(wrong @ source - target))
        epsilon = h / np.sqrt(sum(x*x for x in p) + m*m)
        negative_weights.append(float(np.vdot(seed, (np.eye(4) - epsilon) @ seed).real / 2))
    assert exact_error < TOL
    assert min(beta_difference) > 0.1
    assert all(0 < s < 1 for s in negative_weights)
    return dict(momentum_blocks=3, time_steps=4,
                pulse_area=float(sum(areas)), a0_handoff_error=exact_error,
                beta_mass_mixing_difference_min=min(beta_difference),
                beta_mass_mixing_difference_max=max(beta_difference),
                retained_negative_spectrum_weights=negative_weights,
                scope="Finite spectral algebra only; no spatial-grid or continuum-vacuum claim.")


def group_energy(d):
    s, f, u, h, q, ann = (d[k] for k in ("s", "f", "u", "h", "q", "ann"))
    e = np.eye(6) - f @ f.conj().T
    r = np.eye(2) - f.conj().T @ s @ f
    c0 = e @ (s + s @ f @ np.linalg.solve(r, f.conj().T @ s)) @ e
    conditional_error = norm(c0 - covariance(d["base"], ann))
    probability_error = abs(d["probability"] - np.linalg.det(r).real)
    amplitudes = np.array([np.sqrt(0.37), np.sqrt(0.63) * np.exp(0.71j)])
    rho = np.outer(amplitudes, amplitudes.conj())
    source = d["j"] @ amplitudes
    csource = c0 + f @ rho @ f.conj().T
    code_covariance_error = norm(csource - covariance(source, ann))
    delta = csource - s
    legs = np.column_stack((f, s @ f))
    leg_q, _ = np.linalg.qr(legs)
    rank_residual = norm((np.eye(6) - leg_q @ leg_q.conj().T) @ delta)
    omega = (q * np.abs(d["spectrum"])) @ q.conj().T
    root_omega = (q * np.sqrt(np.abs(d["spectrum"]))) @ q.conj().T
    epsilon = np.eye(6) - 2 * s
    comm = epsilon @ u - u @ epsilon
    weighted_hs_energy = norm(root_omega @ comm)**2 / 4
    svac_out = u @ s @ u.conj().T
    trace_energy = float(np.trace(h @ (svac_out - s)).real)
    hf = sum(h[i, j] * ann[i].conj().T @ ann[j]
             for i in range(6) for j in range(6)) - np.trace(h @ s) * np.eye(64)
    vacuum_out = d["lift"] @ d["vacuum"]
    fock_energy = float(np.vdot(vacuum_out, hf @ vacuum_out).real)
    code_out = d["lift"] @ source
    code_energy = float(np.vdot(code_out, hf @ code_out).real)
    correction = float(np.trace(h @ u @ delta @ u.conj().T).real)
    prepared_energy = float(np.vdot(source, hf @ source).real)
    errors = [conditional_error, probability_error, code_covariance_error, rank_residual,
              abs(weighted_hs_energy - trace_energy), abs(trace_energy - fock_energy),
              abs(code_energy - fock_energy - correction),
              norm(covariance(code_out, ann) - u @ csource @ u.conj().T)]
    unweighted_wrong_energy = norm(comm)**2 / 4
    assert max(errors) < TOL
    assert min(weighted_hs_energy, prepared_energy, code_energy) >= -TOL
    assert abs(unweighted_wrong_energy - weighted_hs_energy) > 0.1
    assert np.linalg.matrix_rank(delta, tol=1e-9) <= 4
    return dict(conditional_covariance_error=conditional_error,
                empty_probability_determinant_error=probability_error,
                coherent_code_covariance_error=code_covariance_error,
                finite_rank=int(np.linalg.matrix_rank(delta, tol=1e-9)),
                finite_rank_legs_error=rank_residual,
                weighted_hs_energy=weighted_hs_energy, trace_energy=trace_energy,
                full_fock_vacuum_energy=fock_energy,
                deleting_energy_weight_error=abs(unweighted_wrong_energy - weighted_hs_energy),
                prepared_code_energy=prepared_energy, output_code_energy=code_energy,
                finite_rank_energy_correction=correction,
                output_energy_formula_error=errors[6], output_covariance_error=errors[7])


def compare(actual, expected, path="results"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for k in actual:
            compare(actual[k], expected[k], path + "." + k)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, path + f"[{i}]")
    elif isinstance(actual, float):
        assert np.isfinite(actual) and np.isclose(actual, expected, atol=1e-10, rtol=1e-10), path
    else:
        assert actual == expected, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    d = fixture()
    result = dict(round=1127, schema=1, status="PASS",
                  scope="Three finite algebra diagnostics; continuum claims are proved in proof.md.",
                  car=group_car(d), dirac=group_dirac(), energy=group_energy(d))
    path = HERE / "results.json"
    if args.write_results:
        path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print("PASS: three diagnostics; results.json registered explicitly.")
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
        print("PASS: three diagnostics; recorded results match; no files written.")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
