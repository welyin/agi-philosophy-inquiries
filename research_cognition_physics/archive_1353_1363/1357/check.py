"""Round 1357: single-loop matter echoes and a pure-field magnetic pulse."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def unitary(hamiltonian, time):
    values, modes = np.linalg.eigh(hamiltonian)
    return (modes * np.exp(-1j * time * values)) @ modes.conj().T


def mean(op, state):
    return float(np.vdot(state, op @ state).real)


class FluxLoop:
    def __init__(self, cutoff, alternate=False):
        self.labels = [(n, m) for n in range(-cutoff, cutoff + 1) for m in range(3)]
        self.indices = {label: index for index, label in enumerate(self.labels)}
        self.size = len(self.labels)
        self.eye = np.eye(self.size, dtype=complex)
        self.vacuum = self.eye[:, self.indices[0, 0]]
        self.energy = np.diag([n*n + ((n % 3 != 0) + ((n + m) % 3 != 0)
                                    if alternate else int(m != 0))
                               for n, m in self.labels]).astype(complex)
        self.locked = np.zeros_like(self.eye)
        self.discrete = np.zeros_like(self.eye)
        for j, (n, m) in enumerate(self.labels):
            if n < cutoff:
                self.locked[self.indices[n+1, (m-1) % 3], j] = 1
            self.discrete[self.indices[n, (m+1) % 3], j] = 1
        self.magnetic = self.discrete + self.discrete.conj().T
        # Controlled finite-time reversible exchange on the only two electric
        # levels used for matter interrogation; use identity off these levels.
        self.exchange = self.eye.copy()
        for m in range(3):
            i, j = self.indices[0, m], self.indices[1, (m-1) % 3]
            self.exchange[i, i] = self.exchange[j, j] = 0
            self.exchange[i, j] = self.exchange[j, i] = 1

    def compress(self, d):
        return np.column_stack([self.eye[:,i] for i,(n,m) in enumerate(self.labels)
                                if (n+m) % 3 == d])


def echo(loop, pulse_on, time):
    state = unitary(loop.magnetic, math.pi/3) @ loop.vacuum if pulse_on else loop.vacuum
    energy0 = mean(loop.energy, state)
    h = np.array([[1, 1], [1, -1]])/math.sqrt(2)
    split = np.kron(h, loop.eye)
    ctrl = np.kron(np.diag([1, 0]), loop.eye) + np.kron(np.diag([0, 1]), loop.exchange)
    total_e = np.kron(np.eye(2), loop.energy)
    path_start = split @ np.kron([1., 0.], state)
    after_shift = ctrl @ path_start
    after_wait = np.kron(np.eye(2), unitary(loop.energy, time)) @ after_shift
    after_unshift = ctrl.conj().T @ after_wait
    bins = (split @ after_unshift).reshape(2, -1)
    probabilities = np.sum(np.abs(bins)**2, axis=1)
    assert abs(probabilities.sum() - 1) < 1.e-12
    for i,(n, _) in enumerate(loop.labels):
        if n:
            assert np.max(np.abs(bins[:,i])) < 1.e-12
    work = [energy0, mean(total_e, after_shift)-energy0,
            mean(total_e, after_unshift)-mean(total_e, after_wait)]
    assert abs(sum(work)-energy0) < 1.e-12
    final_e = [mean(loop.energy, v) for v in bins]
    assert abs(sum(final_e)-energy0) < 1.e-12
    # Record the full subnormalized flux states on the three occupied labels.
    flux = [loop.indices[0,m] for m in range(3)]
    densities = [np.outer(v[flux], v[flux].conj()) for v in bins]
    split_complex = lambda x: {"real": x.real.tolist(), "imag": x.imag.tolist()}
    return {"plus_minus": probabilities.tolist(), "work_in_units_kappa": work,
            "final_branch_weighted_energy": final_e,
            "flux_branch_densities": [split_complex(rho) for rho in densities]}


def test_compression():
    loop = FluxLoop(3)
    emb = loop.compress(0)
    projector = emb @ emb.conj().T
    record = np.array([[0, 1j], [-1j, .4]])
    h = (np.kron(loop.energy, np.eye(2)) +
         .31*np.kron(loop.locked + loop.locked.conj().T, record))
    joint_emb = np.kron(emb, np.eye(2))
    h_eff = joint_emb.conj().T @ h @ joint_emb
    u, v = unitary(h, .73), unitary(h_eff, .73)
    residual = np.linalg.norm(u @ joint_emb - joint_emb @ v)
    assert residual < 1.e-12
    # The 2-level record can be entangled with a passive system.
    rng = np.random.default_rng(1357)
    psi = rng.normal(size=(h_eff.shape[0], 4))+1j*rng.normal(size=(h_eff.shape[0],4))
    psi /= np.linalg.norm(psi)
    assert np.linalg.norm(u @ joint_emb @ psi - joint_emb @ v @ psi) < 1.e-12
    assert np.linalg.norm(projector @ loop.discrete @ loop.vacuum) < 1.e-12
    # Retain all d labels: this is exactly one rotor with a three-state memory.
    perm = np.column_stack([loop.eye[:,loop.indices[n,(d-n)%3]]
                            for n in range(-3,4) for d in range(3)])
    shift = np.diag(np.ones(6), -1)
    d_shift = np.roll(np.eye(3), 1, axis=0)
    assert np.allclose(perm.T @ loop.locked @ perm, np.kron(shift, np.eye(3)))
    assert np.allclose(perm.T @ loop.discrete @ perm, np.kron(np.eye(7), d_shift))
    return float(residual)


def calculate():
    loop = FluxLoop(3)
    twin = FluxLoop(3, alternate=True)
    times = [math.pi/2, math.pi]
    scenarios = {}
    for name, model, pulse, expected in [
            ("no_pulse", loop, False, [1., 0.]),
            ("pure_gauge_pulse", loop, True, [1/3, 4/9]),
            ("alternative_off_sector_energy", twin, True, [1., 0.])]:
        rows = [echo(model, pulse, tau) for tau in times]
        assert np.allclose([r["plus_minus"][1] for r in rows], expected)
        scenarios[name] = rows
    for model in [loop, twin]:
        assert np.allclose([model.energy[model.indices[0,m], model.indices[0,m]]
                            for m in range(3)], [0., 1., 1.])
    e = loop.compress(0)
    assert np.allclose(e.T @ loop.energy @ e, e.T @ twin.energy @ e)
    assert np.allclose(scenarios["pure_gauge_pulse"][0]["work_in_units_kappa"],
                       [8/9, 1/3, -1/3])
    assert np.allclose(scenarios["alternative_off_sector_energy"][0]["work_in_units_kappa"],
                       [8/9, 1., -1.])
    radius = .04
    lower = 1/3 - radius
    single_gap_floor = 4*lower*(1-lower) - radius
    assert single_gap_floor > 4/9 + radius
    return {"dimensionless_waits": times, "scenarios": scenarios,
            "certified_absolute_probability_error": radius,
            "single_gap_second_lower_including_error": single_gap_floor,
            "compression_residual": test_compression()}


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            compare(actual[key], saved[key])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for a, b in zip(actual, saved):
            compare(a, b)
    else:
        assert np.isclose(actual, saved, rtol=1.e-10, atol=1.e-12)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = calculate()
    dest = Path(__file__).with_name("results.json")
    if args.write:
        with dest.open("x", encoding="utf-8") as out:
            json.dump(result, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(result, json.loads(dest.read_text(encoding="utf-8")))
    print({key: [row["plus_minus"][1] for row in values]
           for key, values in result["scenarios"].items()})
    print("unitary compression residual", result["compression_residual"])
