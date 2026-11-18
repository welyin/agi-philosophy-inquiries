"""Construct a common autonomous generator for three prescribed finite-dimensional rays."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def adj(matrix):
    return matrix.conj().T


def common_populations(rays):
    size = rays.shape[0]
    projectors = [np.outer(v, v.conj()) for v in rays.T]
    first, second = [p - projectors[0] for p in projectors[1:]]
    residual = np.eye(size, dtype=complex)
    basis = []
    for dim in range(size, 1, -1):
        _, vectors = np.linalg.eigh(adj(residual) @ first @ residual)
        frame = residual @ vectors
        indices = np.arange(dim)
        phases = 2*np.pi*np.outer(indices, indices)/dim
        def direction(k):
            return frame @ np.exp(1j * k) / np.sqrt(dim)
        def sign(k):
            v = direction(k)
            return float(np.vdot(v, second @ v).real)
        values = [sign(k) for k in phases]
        lo, hi = phases[np.argmin(values)], phases[np.argmax(values)]
        if min(abs(s) for s in values) <= 1e-13:
            v = direction(phases[np.argmin(np.abs(values))])
        else:
            assert min(values) < 0 < max(values)
            for _ in range(60):
                mid = (lo + hi) / 2
                if sign(mid) < 0:
                    lo = mid
                else:
                    hi = mid
            v = direction((lo+hi)/2)
        basis.append(v)
        # Complete the chosen vector in the current residual coordinates.
        q, _ = np.linalg.qr((adj(residual) @ v)[:, None], mode="complete")
        residual = residual @ q[:, 1:]
    frame = np.column_stack(basis+[residual[:, 0]])
    coordinates = adj(frame) @ rays
    assert np.allclose(np.abs(coordinates), np.abs(coordinates[:, :1]), atol=2e-7)
    return frame, coordinates


def construct(rays, count, ratio):
    frame, v = common_populations(rays)
    phases = np.angle(v[:, 1] * v[:, 0].conj())
    branches = np.arange(count)
    angles = phases[:, None] + 2*np.pi*branches
    # Minimize Euclidean component errors, including modes of zero weight.
    errors = abs(np.exp(1j*ratio*angles)*v[:, 0, None] - v[:, 2, None])
    selected = angles[np.arange(v.shape[0]), np.argmin(errors, axis=1)]
    advance = lambda time: frame @ (np.exp(1j*time*selected) * v[:, 0])
    first_error = np.linalg.norm(advance(1)-rays[:, 1])
    second_error = np.linalg.norm(advance(ratio)-rays[:, 2])
    assert first_error < 2e-7
    spread = np.ptp(selected)
    angles_fs = [np.arccos(np.clip(abs(np.vdot(rays[:, j], rays[:, j+1])), 0, 1))
                 for j in range(2)]
    return dict(first_error=float(first_error), last_error=float(second_error),
                bandwidth_tau_over_hbar=float(spread), target_angle_gap=float(abs(np.diff(angles_fs)[0])))


def run():
    rng = np.random.default_rng(1342)
    ratio = 1 + np.sqrt(2)
    samples = []
    for n in (2, 3, 4, 7):
        for _ in range(5):
            rays = rng.normal(size=(n, 3))+1j*rng.normal(size=(n, 3))
            rays /= np.linalg.norm(rays, axis=0)
            thin = construct(rays, 32, ratio)
            fine = construct(rays, 2048, ratio)
            assert fine["last_error"] <= thin["last_error"]+1e-10
            assert fine["last_error"] < .003
            samples.append(fine)
    # Deterministic non-orbit at equal steps: orthogonal at first, identical at last.
    extreme = np.eye(3, dtype=complex)[:, [0, 1, 1]]
    extreme_fit = construct(extreme, 2048, ratio)
    assert extreme_fit["last_error"] < .003
    rng2 = np.random.default_rng(4132)
    max_violation = 0.
    for _ in range(150):
        e = rng2.uniform(-4., 4., 5)
        p = rng2.dirichlet(np.ones(5))
        t, s = rng2.uniform(0, 2, size=2)
        fs = lambda time: math.acos(np.clip(abs(p @ np.exp(-1j*e*time)), 0, 1))
        excess = abs(fs(t)-fs(s)) - np.ptp(e)*abs(t-s)/2
        max_violation = max(max_violation, excess)
    assert max_violation < 1e-12
    cut = 2*(math.acos(.87)-math.acos(.96))
    return dict(round=1342, ratio=float(ratio), target_count=len(samples),
                max_second_error=max(z["last_error"] for z in samples),
                max_first_error=max(z["first_error"] for z in samples),
                extreme=extreme_fit, energy_time_required_from_1341=cut,
                finite_band_tests=150,max_sampled_band_violation=max_violation)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = Path(__file__).with_name("results.json")
    result = run()
    if args.write:
        with out.open("x", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
            f.write("\n")
    else:
        saved = json.loads(out.read_text(encoding="utf-8"))
        assert result.keys() == saved.keys()
        for name in result:
            if name == "extreme":
                # The completing basis can move within a degenerate eigenspace.
                assert result[name]["last_error"] < .003 and saved[name]["last_error"] < .003
            elif name.startswith("max_"):
                assert abs(result[name]-saved[name]) < .003
            else:
                assert np.isclose(result[name], saved[name])
    print(json.dumps(result))
