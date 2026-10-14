"""Read-only small-K cross-check of mean_field_bound.md; not a research round.

NumPy only. --record creates the result exactly once; the default independently
recomputes and compares. RK4/eigendecomposition checks are numerical diagnostics,
not interval certificates and not substitutes for the analytical inequalities.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RESULT = HERE / "mean_field_bound_checks.json"
PROOF = HERE / "mean_field_bound.md"
PROOF_SHA = "11ac64ce824d3448dec97548efe3dda768bcda30f2e48e1cd88bfe511fcd6b3f"
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                 [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)
TIME = 0.65
FINE_STEPS = 1024
COARSE_STEPS = 512


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(a: np.ndarray) -> float:
    return float(np.linalg.norm(a, 2))


def cpair(a: np.ndarray) -> list:
    return np.stack((a.real, a.imag), axis=-1).tolist()


def rand_hermitian(rng: np.random.Generator, d: int, size: float) -> np.ndarray:
    z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    z = (z + z.conj().T) / 2
    return size * z / norm(z)


def embedded_one(a: np.ndarray, site: int, count: int) -> np.ndarray:
    size = 2 ** count
    columns = np.arange(size)
    shift = count - site - 1
    bit = (columns >> shift) & 1
    cleared = columns & ~(1 << shift)
    out = np.zeros((size, size), dtype=complex)
    for v in range(2):
        out[cleared | (v << shift), columns] = a[v, bit]
    return out


def embedded_two(a: np.ndarray, i: int, j: int, count: int) -> np.ndarray:
    size = 2 ** count
    columns = np.arange(size)
    si, sj = count - i - 1, count - j - 1
    pair = 2 * ((columns >> si) & 1) + ((columns >> sj) & 1)
    cleared = columns & ~((1 << si) | (1 << sj))
    out = np.zeros((size, size), dtype=complex)
    for v in range(4):
        rows = cleared | ((v // 2) << si) | ((v % 2) << sj)
        out[rows, columns] = a[v, pair]
    return out


def tensor_power(vector: np.ndarray, count: int) -> np.ndarray:
    value = np.ones(1, dtype=complex)
    for _ in range(count):
        value = np.kron(value, vector)
    return value


def bg_hamiltonian(h: np.ndarray, wmat: np.ndarray, count: int) -> np.ndarray:
    out = sum(embedded_one(h, i, count) for i in range(count))
    for i in range(count):
        for j in range(i + 1, count):
            out += embedded_two(wmat, i, j, count) / count
    return out


def evolve(h: np.ndarray, state: np.ndarray, t: float) -> np.ndarray:
    values, basis = np.linalg.eigh(h)
    phases = np.exp(-1j * t * values)
    if state.ndim == 1:
        return basis @ (phases * (basis.conj().T @ state))
    return basis @ (phases[:, None] * (basis.conj().T @ state))


def hartree_operator(wmat: np.ndarray, phi: np.ndarray) -> np.ndarray:
    return np.einsum("a,iajb,b->ij", phi.conj(),
                     wmat.reshape(2, 2, 2, 2), phi)


def integrate(h, wmat, phi0, hp, aa, bb, steps):
    """Evolve Hartree phi and its common probe comparator simultaneously."""
    dt = TIME / steps
    state = np.concatenate((phi0, I2.ravel()))

    def rhs(y):
        phi, up = y[:2], y[2:].reshape(2, 2)
        wh = hartree_operator(wmat, phi)
        heff = hp.copy()
        for a, b in zip(aa, bb):
            heff += a * float(np.vdot(phi, b @ phi).real)
        return np.concatenate((-1j * (h + wh) @ phi,
                               (-1j * heff @ up).ravel()))

    samples = [(0.0, state[:2].copy(), state[2:].reshape(2, 2).copy())]
    for n in range(1, steps + 1):
        k1 = rhs(state)
        k2 = rhs(state + dt * k1 / 2)
        k3 = rhs(state + dt * k2 / 2)
        k4 = rhs(state + dt * k3)
        state += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        if n % (steps // 4) == 0:
            samples.append((n * dt, state[:2].copy(),
                            state[2:].reshape(2, 2).copy()))
    return samples


def prepared_check(hk, hp, aa, bbar, initial, up, final_bg, t):
    dim = len(initial)
    full = np.kron(I2, hk) + np.kron(hp, np.eye(dim))
    for a, b in zip(aa, bbar):
        full += np.kron(a, b)
    preparation = np.kron(I2, initial[:, None])
    actual = evolve(full, preparation, t)
    comparison = np.kron(up, final_bg[:, None])
    delta = actual - comparison
    error = norm(delta)
    # A Bell input explicitly retains a reference; the operator norm is
    # stronger, uniform over every input and every ancillary reference.
    zero, one = I2[:, 0], I2[:, 1]
    bell_delta = (np.kron(delta[:, 0], zero)
                  + np.kron(delta[:, 1], one)) / math.sqrt(2)
    probe = np.array([0.8, 0.6j])
    joint = (actual @ probe).reshape(2, dim)
    b_actual = sum(np.vdot(row, bbar[0] @ row) for row in joint).real
    b_bg = np.vdot(final_bg, bbar[0] @ final_bg).real
    return {
        "prepared_isometry_error": error,
        "bell_reference_vector_error": float(np.linalg.norm(bell_delta)),
        "all_input_gram_defect": norm(actual.conj().T @ actual - I2),
        "total_energy_compression_drift":
            norm(actual.conj().T @ full @ actual
                 - preparation.conj().T @ full @ preparation),
        "actual_backreaction_B0_difference": float(b_actual - b_bg),
    }, actual


def random_model(seed):
    rng = np.random.default_rng(seed)
    h = rand_hermitian(rng, 2, 0.4)
    raw = rand_hermitian(rng, 4, 1.0)
    wmat = (raw + SWAP @ raw @ SWAP) / 2
    wmat *= 0.18 / norm(wmat)
    bb = [rand_hermitian(rng, 2, 1.0), rand_hermitian(rng, 2, 0.8)]
    aa = [0.06 * Z, 0.04 * X]
    hp = 0.17 * X + 0.11 * Y
    phi0 = rng.normal(size=2) + 1j * rng.normal(size=2)
    phi0 /= np.linalg.norm(phi0)
    return h, wmat, aa, bb, hp, phi0


def main_checks():
    assert sha(PROOF) == PROOF_SHA, "The completed proof changed."
    random_rows = []
    correction_rows = []
    integrity = []
    max_rk_difference = 0.0
    for seed in (731, 1947):
        h, wmat, aa, bb, hp, phi0 = random_model(seed)
        w = norm(wmat)
        strength = sum(norm(a) * norm(b) for a, b in zip(aa, bb))
        fine = integrate(h, wmat, phi0, hp, aa, bb, FINE_STEPS)
        coarse = integrate(h, wmat, phi0, hp, aa, bb, COARSE_STEPS)
        convergence = max(float(np.linalg.norm(f[1] - c[1]))
                          + norm(f[2] - c[2]) for f, c in zip(fine, coarse))
        max_rk_difference = max(max_rk_difference, convergence)
        integrity.append({
            "seed": seed, "w": w, "L": strength,
            "h": cpair(h), "W": cpair(wmat), "phi0": cpair(phi0),
            "B": [cpair(b) for b in bb],
            "swap_error": norm(SWAP @ wmat @ SWAP - wmat),
            "B_noncommutator": norm(bb[0] @ bb[1] - bb[1] @ bb[0]),
            "RK4_coarse_fine_difference": convergence,
        })
        for count in range(2, 8):
            initial = tensor_power(phi0, count)
            hk = bg_hamiltonian(h, wmat, count)
            bars = [sum(embedded_one(b, i, count) for i in range(count))
                    / count for b in bb]
            values, basis = np.linalg.eigh(hk)
            coefficients = basis.conj().T @ initial
            audits = []
            for t, phi, up in fine:
                bg = basis @ (np.exp(-1j * t * values) * coefficients)
                p = np.outer(phi, phi.conj())
                q = I2 - p
                q1 = embedded_one(q, 0, count)
                alpha = float(np.vdot(bg, q1 @ bg).real)
                mf1 = embedded_one(h + hartree_operator(wmat, phi), 0, count)
                lmat = hk - mf1
                dot_alpha = float((1j * np.vdot(
                    bg, (lmat @ q1 - q1 @ lmat) @ bg)).real)
                frows = []
                for b, bar in zip(bb, bars):
                    beta = norm(b)
                    mean = float(np.vdot(phi, b @ phi).real)
                    ss = sum(embedded_one(p @ b @ q, i, count)
                             for i in range(count)) / count
                    dd = sum(embedded_one(q @ (b - mean * I2) @ q, i, count)
                             for i in range(count)) / count
                    fluctuation = float(np.linalg.norm((bar - mean * np.eye(len(bg))) @ bg))
                    bound = 4 * beta * math.exp(6 * w * t) / math.sqrt(count)
                    comm = ss @ ss.conj().T - ss.conj().T @ ss
                    comm_margin = float(np.linalg.eigvalsh(
                        beta ** 2 / count * np.eye(len(bg)) - comm).min())
                    identity_error = norm(bar - mean * np.eye(len(bg))
                                          - ss - ss.conj().T - dd)
                    assert fluctuation <= bound + 1e-10
                    assert comm_margin >= -1e-10
                    assert identity_error < 1e-9
                    frows.append({
                        "fluctuation": fluctuation, "bound": bound,
                        "raising_commutator_margin": comm_margin,
                        "decomposition_error": identity_error,
                    })
                alpha_bound = math.expm1(12 * w * t) / count
                dot_bound = 12 * w * (max(alpha, 0.0) + 1 / count)
                assert alpha <= alpha_bound + 1e-10 and alpha >= -1e-10
                assert abs(dot_alpha) <= dot_bound + 1e-10
                audits.append({
                    "time": t, "alpha": alpha, "alpha_bound": alpha_bound,
                    "abs_dot_alpha": abs(dot_alpha), "dot_bound": dot_bound,
                    "fluctuations": frows,
                })
            t, phi, up = fine[-1]
            bg = basis @ (np.exp(-1j * t * values) * coefficients)
            whole, _ = prepared_check(hk, hp, aa, bars, initial, up, bg, t)
            bound = 4 * strength * t * math.exp(6 * w * t) / math.sqrt(count)
            whole["analytic_bound"] = bound
            overlap = abs(np.vdot(tensor_power(phi, count), bg))
            whole["phase_optimized_BG_product_distance"] = math.sqrt(max(0, 2 - 2 * overlap))
            assert whole["prepared_isometry_error"] <= bound + 1e-9
            assert whole["bell_reference_vector_error"] <= whole["prepared_isometry_error"] + 1e-10
            assert whole["all_input_gram_defect"] < 1e-9
            assert whole["total_energy_compression_drift"] < 1e-9
            random_rows.append({"seed": seed, "K": count, "background_samples": audits,
                                "complete_output": whole})

            # The exact corrected background is retained in BOTH comparators.
            dmat = 0.09 * I2 + 0.05 * X + 0.03 * Z
            d = norm(dmat)
            hk_corr = bg_hamiltonian(h + dmat / count, wmat, count)
            bg_corr = evolve(hk_corr, initial, t)
            fine_corr = integrate(h + dmat / count, wmat, phi0, hp, aa, bb, FINE_STEPS)
            _, phi_corr, up_corr = fine_corr[-1]
            phi_error = float(np.linalg.norm(phi_corr - phi))
            phi_bound = d * t * math.exp(2 * w * t) / count
            b_errors = [abs(float(np.vdot(phi_corr, b @ phi_corr).real
                                   - np.vdot(phi, b @ phi).real)) for b in bb]
            b_bounds = [2 * norm(b) * phi_bound for b in bb]
            out_corr, actual_corr = prepared_check(
                hk_corr, hp, aa, bars, initial, up_corr, bg_corr, t)
            out_uncorrected = norm(actual_corr - np.kron(up, bg_corr[:, None]))
            correction_budget = strength * d * t * t * math.exp(2 * w * t) / count
            exact_correction = hk_corr - hk
            assert phi_error <= phi_bound + 1e-10
            assert all(a <= b + 1e-10 for a, b in zip(b_errors, b_bounds))
            assert out_corr["prepared_isometry_error"] <= bound + 1e-9
            assert out_uncorrected <= bound + correction_budget + 1e-9
            correction_rows.append({
                "seed": seed, "K": count, "D_norm": d,
                "full_background_correction_norm": norm(exact_correction),
                "Hartree_vector_difference": phi_error, "Hartree_bound": phi_bound,
                "B_differences": b_errors, "B_bounds": b_bounds,
                "corrected_comparator_complete_error": out_corr["prepared_isometry_error"],
                "unmodified_probe_with_exact_corrected_BG_error": out_uncorrected,
                "main_bound": bound, "additional_probe_budget": correction_budget,
                "actual_corrected_energy_drift": out_corr["total_energy_compression_drift"],
            })
    assert max_rk_difference < 1e-9

    # Solvable one-axis twisting: mean one-body excitation fraction goes to
    # zero, but the phase-optimized WHOLE state does not approach the product.
    product_control = []
    for count in (2, 3, 4, 5, 6, 7, 32, 128, 512):
        amplitude = 0j
        for j in range(count + 1):
            weight = math.exp(math.lgamma(count + 1) - math.lgamma(j + 1)
                              - math.lgamma(count - j + 1) - count * math.log(2))
            energy = ((count - 2 * j) ** 2 - count) / (2 * count)
            amplitude += weight * np.exp(-1j * energy)
        alpha = (1 - math.cos(2 / count) ** (count - 1)) / 2
        product_control.append({
            "K": count, "time": 1.0, "W_norm": 1.0, "alpha_exact": alpha,
            "whole_product_overlap_absolute": abs(amplitude),
            "phase_optimized_whole_product_distance":
                math.sqrt(max(0, 2 - 2 * abs(amplitude))),
        })
    assert product_control[-1]["phase_optimized_whole_product_distance"] > 0.5
    assert product_control[-1]["alpha_exact"] < 0.003

    return {
        "schema": "mean_field_bound_small_K_numerical_crosscheck_v1",
        "status": "PASS_NUMERICAL_CROSSCHECK_ONLY",
        "proof_sha256": sha(PROOF), "checker_sha256": sha(Path(__file__)),
        "numpy_version": np.__version__,
        "settings": {
            "cluster_dimension": 2, "probe_dimension": 2,
            "explicit_reference_dimension": 2, "K_values": list(range(2, 8)),
            "time": TIME, "RK4_fine_steps": FINE_STEPS,
            "RK4_coarse_steps": COARSE_STEPS,
            "unknown_probe_reference_coverage":
                "full prepared-isometry operator norm, plus explicit Bell input",
            "numerical_status":
                "floating eigendecomposition and RK4; no interval certification",
        },
        "random_models": integrity,
        "background_and_complete_output": random_rows,
        "finite_K_self_term_correction": correction_rows,
        "whole_product_nonconvergence_control": {
            "H": "K^-1 sum_(r<s) Z_r Z_s; initial |+>^K; h=0",
            "Hartree": "phi(t)=|+>; b_Z=0",
            "rows": product_control,
            "limiting_overlap_absolute": 2 ** (-0.25),
            "limiting_phase_optimized_distance": math.sqrt(2 - 2 * 2 ** (-0.25)),
        },
        "model_mapping_audit": {
            "status": "generic finite tensor-copy hypotheses checked; concrete spatial candidate not certified here",
            "necessary_author_data": [
                "explicit K-independent cluster h and swap-symmetric W",
                "separate same-copy diagonal term D/K from r<s interaction",
                "even local CAR factors and fixed spatial support",
                "strict bounds for w=||W|| and L=sum ||A_a||||B_a||",
                "fixed physical finite domain and product preparation resource",
                "actual clock/endpoint/signal menu and common geometric comparison",
            ],
            "not_claimed": [
                "all-background product convergence",
                "arbitrary initial background/reference classicalization",
                "actual CAR or spatial transport construction from qubit diagnostic",
                "numerical proof of the analytical mean-field theorem",
                "new research round or accepted physical result",
            ],
        },
    }


def compare(actual, saved, path="root"):
    if isinstance(actual, dict):
        assert set(actual) == set(saved), path
        for key in actual:
            compare(actual[key], saved[key], path + "." + key)
    elif isinstance(actual, list):
        assert len(actual) == len(saved), path
        for n, (a, b) in enumerate(zip(actual, saved)):
            compare(a, b, path + f"[{n}]")
    elif isinstance(actual, float):
        assert math.isclose(actual, saved, rel_tol=2e-8, abs_tol=2e-10), (path, actual, saved)
    else:
        assert actual == saved, (path, actual, saved)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", action="store_true",
                        help="Create the result once; never overwrite.")
    args = parser.parse_args()
    if args.record and RESULT.exists():
        raise SystemExit("Refusing to overwrite existing result.")
    data = main_checks()
    if args.record:
        with RESULT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
        mode = "RECORDED_ONCE"
    else:
        compare(data, json.loads(RESULT.read_text(encoding="utf-8")))
        mode = "READ_ONLY_MATCH"
    print(json.dumps({
        "status": mode, "random_model_K_cases": len(data["background_and_complete_output"]),
        "finite_K_correction_cases": len(data["finite_K_self_term_correction"]),
        "max_prepared_isometry_error": max(
            x["complete_output"]["prepared_isometry_error"]
            for x in data["background_and_complete_output"]),
        "whole_product_control_K512_distance":
            data["whole_product_nonconvergence_control"]["rows"][-1]["phase_optimized_whole_product_distance"],
        "proof_sha256": data["proof_sha256"],
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
