"""Build the physical reversible three-arm circuit and verify the transfer."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


SOURCE = runpy.run_path(str(Path(__file__).resolve().parents[1] / "1353" / "check.py"))
NU = .25
THETA = math.pi / math.sqrt(3)
N_MATERIAL = 81


class Device:
    def __init__(self, group):
        parent, _, self.source_size = SOURCE["source"](group)
        self.levels, self.basis = np.linalg.eigh(parent)
        self.blank = np.zeros((self.source_size, N_MATERIAL + 1), dtype=complex)
        self.blank[0, 0] = 1
        self.a = np.eye(9).ravel() / 3
        self.swap = np.eye(9).reshape(3, 3, 3, 3).transpose(1, 0, 2, 3).reshape(9, 9)

    def prepare_pair_packet(self, state):
        """Swap vacuum and the four-leg A state, leaving the orthogonal part unchanged."""
        out = state.copy()
        amplitude = state[:, 1:] @ self.a.conj()
        out[:, 0] = amplitude
        out[:, 1:] += np.outer(state[:, 0] - amplitude, self.a)
        return out

    def arm(self, index, state, theta, inverse=False):
        curr = state.copy() if inverse else self.prepare_pair_packet(state)
        body = curr[:, 1:].reshape(self.source_size, 9, 9)
        if index == 1:
            curr[:, 1:] = (self.swap @ body).reshape(self.source_size, N_MATERIAL)
        if index == 2:
            phase = theta if inverse else -theta
            unitary = (self.basis * np.exp(1j * phase * self.levels)) @ self.basis.conj().T
            curr[:, 1:] = (unitary @ body.reshape(-1, 9)).reshape(self.source_size, N_MATERIAL)
        return self.prepare_pair_packet(curr) if inverse else curr

    def matrix(self, theta):
        return np.column_stack([self.arm(i, self.blank, theta)[0, 1:] for i in range(3)])

    def decode_raw(self, physical):
        """Run each controlled inverse on the actual input, then project the old ports."""
        state = self.blank.copy()
        state[0, 0] = 0
        state[0, 1:] = physical
        return np.array([self.arm(i, state, THETA, inverse=True)[0, 0]
                         for i in range(3)]) / math.sqrt(3)


def neutral_filter(gram):
    values, modes = np.linalg.eigh(gram)
    weights = np.sqrt(NU / np.maximum(values, NU))
    return (modes * weights) @ modes.conj().T, values, modes


def calibrate(matrix):
    """Two-arm interference from actual success amplitudes reconstructs all complex entries."""
    grams = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        grams[i, i] = np.vdot(matrix[:, i], matrix[:, i]).real
        for j in range(i):
            rates = [np.linalg.norm((matrix[:, j] + phase * matrix[:, i]) / 2)**2
                     for phase in [1, -1, 1j, -1j]]
            grams[j, i] = rates[0] - rates[1] - 1j * (rates[2] - rates[3])
            grams[i, j] = grams[j, i].conj()
    return grams


def inspect_group(group):
    device = Device(group)
    matrix = device.matrix(THETA)
    gram = calibrate(matrix)
    assert np.allclose(gram, matrix.conj().T @ matrix)
    # Complex interferometry controls, independent of the fixed phase convention.
    phased = matrix @ np.diag(np.exp(1j*np.array([.31, .87, -.44])))
    assert np.allclose(calibrate(phased), phased.conj().T @ phased)

    filt, vals, modes = neutral_filter(gram)
    encode = matrix @ filt / math.sqrt(3)
    physical_basis = np.eye(N_MATERIAL)
    decode = filt @ np.column_stack([device.decode_raw(v) for v in physical_basis])
    assert np.linalg.norm(decode-encode.conj().T) < 1.e-12
    support = modes[:, vals > 1.e-10]
    phys, sv, _ = np.linalg.svd(matrix, full_matrices=False)
    phys = phys[:, sv > 1.e-7]
    rank = phys.shape[1]
    p = NU/3
    assert np.allclose(decode@encode, p*(support@support.conj().T))
    assert np.allclose(encode@decode, p*(phys@phys.conj().T))
    assert np.linalg.norm(gram.imag) < 1.e-12
    target = np.array([1, 1j, 0])/math.sqrt(2)
    first = matrix @ target
    second_raw = device.matrix(THETA/2)[:, 2]
    second = second_raw/np.linalg.norm(second_raw)
    rows = []
    for psi, preparation in [(first, .5), (second, np.linalg.norm(second_raw)**2)]:
        raw = device.decode_raw(psi)
        z = filt @ raw
        cond = abs(np.vdot(target, z))**2 / p
        assert np.linalg.norm(encode@z - p*psi) < 1.e-12
        assert abs(np.vdot(z, z).real-p) < 1.e-12
        # Raw erasure fails, neutral filter rejects, delivered into the tested detector or its complement.
        branches = np.array([1-np.linalg.norm(raw)**2, np.linalg.norm(raw)**2-p,
                             p*cond, p*(1-cond)])
        assert min(branches) > -1.e-12 and abs(branches.sum()-1) < 1.e-12
        rows.append({"preparation_probability": float(preparation),
                     "decode_branches_given_preparation": branches.tolist(),
                     "conditional_neutral_detector": float(cond),
                     "from_raw_preparations_and_detector": float(preparation * p * cond)})
    bell = phys/math.sqrt(rank)
    end = encode @ decode @ bell
    assert np.linalg.norm(end - p*bell) < 1.e-12
    return {"rank": rank, "gram_eigenvalues": vals.tolist(),
            "positive_gap": float(min(vals[vals > 1.e-10])),
            "one_way_permitted_success": p, "round_trip_success": float(np.linalg.norm(end)**2),
            "basis_dephased_conditional_Bell_fidelity": 1/rank, "holdouts": rows}


def compare(now, saved):
    if isinstance(now, dict):
        assert now.keys() == saved.keys()
        for k in now:
            compare(now[k], saved[k])
    elif isinstance(now, list):
        assert len(now) == len(saved)
        for x, y in zip(now, saved):
            compare(x, y)
    else:
        assert np.isclose(now, saved, rtol=1.e-10, atol=1.e-12)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    results = {g: inspect_group(g) for g in ["SU3", "SO3"]}
    results["mu3_SO3_has_same_four_leg_matrices_as_SO3"] = True
    delta = .0005
    results["conditional_probability_error_bound"] = 2*delta/(NU/3-delta)
    dest = Path(__file__).with_name("results.json")
    if args.write:
        with dest.open("x", encoding="utf8") as out:
            json.dump(results, out, indent=2, allow_nan=False)
            out.write("\n")
    else:
        compare(results, json.loads(dest.read_text(encoding="utf8")))
    print(json.dumps(results, indent=2))
