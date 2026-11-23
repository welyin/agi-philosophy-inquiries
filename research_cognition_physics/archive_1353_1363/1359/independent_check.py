"""Independent audit of the two number-changing vertices of round 1359."""
import argparse
from itertools import permutations, product
import json
from pathlib import Path

import numpy as np


def signed_orientations():
    for triple in permutations(range(3)):
        inversions = sum(i > j for i, j in [(triple[0], triple[1]),
                                            (triple[0], triple[2]),
                                            (triple[1], triple[2])])
        yield triple, (-1) ** inversions


def build_vertices():
    # A basis is (colour a,b,f ; spin a,b,f), C is (dual pair, f ; spin f).
    a_index = list(product(range(3), range(3), range(3), range(2), range(2), range(2)))
    c_index = list(product(range(3), range(3), range(2)))
    v = np.zeros((18, 216))
    w = np.zeros((2, 18))
    # Independently construct initial colour-volume/singlet and volume/triplet codes.
    singlet_code = np.zeros((216, 2))
    triplet_code = np.zeros_like(singlet_code)
    for (i, j, k), sign in signed_orientations():
        for spin in range(2):
            for x, y, spin_sign in [(0, 1, 1), (1, 0, -1)]:
                for ell in range(3):
                    v[c_index.index((k, ell, spin)), a_index.index((i, j, ell, x, y, spin))] = sign * spin_sign / 2
                singlet_code[a_index.index((i, j, k, x, y, spin)), spin] = sign * spin_sign / np.sqrt(12)
            triplet_code[a_index.index((i, j, k, 0, 0, spin)), spin] = sign / np.sqrt(6)
    for k, ell, spin in c_index:
        w[spin, c_index.index((k, ell, spin))] = float(k == ell) / np.sqrt(3)
    assert np.allclose(v @ v.T, np.eye(18))
    assert np.allclose(w @ w.T, np.eye(2))
    assert np.allclose(w @ v, singlet_code.T)
    assert np.allclose(v @ triplet_code, 0)
    return v, w, singlet_code, triplet_code


def spectrum(v, w, g1=1., g2=1.):
    zero = np.zeros
    hamiltonian = np.block([
        [zero((216, 216)), g1 * v.T, zero((216, 2))],
        [g1 * v, zero((18, 18)), g2 * w.T],
        [zero((2, 216)), g2 * w, zero((2, 2))],
    ])
    return np.linalg.eigh(hamiltonian)


def symmetry_residuals(v, w):
    rng = np.random.default_rng(1359)
    generators = {}
    candidate = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    unitary = np.linalg.qr(candidate)[0]
    unitary[:, 0] *= np.linalg.det(unitary).conjugate()
    generators["SU3"] = unitary
    rot = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    rot[:, 0] *= np.linalg.det(rot)
    generators["SO3"] = rot
    generators["Z3_times_SO3"] = np.exp(2j * np.pi / 3) * rot
    residuals = {}
    for name, matrix in generators.items():
        a = np.kron(np.kron(np.kron(matrix, matrix), matrix), np.eye(8))
        c = np.kron(np.kron(matrix.conj(), matrix), np.eye(2))
        residuals[name] = float(max(np.linalg.norm(v @ a - c @ v), np.linalg.norm(w @ c - w)))
        assert residuals[name] < 1.e-12
    return residuals


def run():
    v, w, code, triplet = build_vertices()
    values, vectors = spectrum(v, w)
    levels = [-np.sqrt(2), -1, 0, 1, np.sqrt(2)]
    multiplicities = [int(np.isclose(values, level).sum()) for level in levels]
    assert multiplicities == [2, 16, 200, 16, 2]
    # A Bell state of the input f spin and a separately identified two-state partner.
    bell = np.eye(2) / np.sqrt(2)
    input_state = code @ bell
    records = {}
    for label, omega_t, phi in [("complete", np.pi, np.pi / 3),
                                ("half", np.pi / 2, 2 * np.pi / 3)]:
        evolved_from_a = (vectors * np.exp(-1j * values * omega_t / np.sqrt(2))) @ vectors[:216, :].T
        transfer = evolved_from_a[-2:, :]
        h_effect = transfer.conj().T @ transfer
        expected_p = np.sin(omega_t / 2) ** 4
        assert np.allclose(h_effect, expected_p * (code @ code.T))
        assert np.allclose(transfer @ triplet, 0)
        evolved = evolved_from_a @ input_state
        a, c, received = np.split(evolved, [216, 234], axis=0)
        occupation = [float(np.linalg.norm(block) ** 2) for block in [a, c, received]]
        dark_joint_amplitude = (received - np.exp(1j * phi) * received) / 2
        dark = float(np.linalg.norm(dark_joint_amplitude) ** 2)
        path_logged = float(np.linalg.norm(received) ** 2 / 2)
        # Completely dephase the spin without tagging the path. The path probabilities survive.
        spin_logged_dark = sum(np.linalg.norm(dark_joint_amplitude[j, :]) ** 2 for j in range(2))
        bell_probability = float(abs(np.vdot(bell, received)) ** 2 / expected_p)
        logged_bell_probability = float(sum(abs(np.vdot(bell[j], received[j])) ** 2
                                            for j in range(2)) / expected_p)
        assert np.isclose(dark, 1/4 if label == "complete" else 3/16)
        assert np.isclose(dark, spin_logged_dark)
        assert np.isclose(bell_probability, 1) and np.isclose(logged_bell_probability, 1/2)
        records[label] = {"occupations_a_c_h": occupation, "h_dark": dark,
                          "h_dark_if_path_orthogonally_logged": path_logged,
                          "h_dark_if_spin_logged": float(spin_logged_dark),
                          "bell_given_h": bell_probability,
                          "bell_given_h_if_spin_logged": logged_bell_probability}
    different, basis = spectrum(v, w, 2, 1)
    t_star = np.pi / np.sqrt(5)
    received = (basis[-2:, :] * np.exp(-1j * different * t_star)) @ basis[:216, :].T @ input_state
    reduced_peak = float(np.linalg.norm(received) ** 2)
    assert np.isclose(reduced_peak, 16/25)
    return {"energies_in_coupling_units": levels, "multiplicities": multiplicities,
            "symmetry_intertwiner_frobenius_residual": symmetry_residuals(v, w),
            "path_spin_records": records, "unequal_2_to_1_peak": reduced_peak}


def verify(current, saved):
    if isinstance(current, dict):
        assert current.keys() == saved.keys()
        for key, item in current.items():
            verify(item, saved[key])
    else:
        assert np.allclose(current, saved, atol=4.e-12, rtol=4.e-12), (current, saved)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write", action="store_true")
    options = cli.parse_args()
    record = run()
    destination = Path(__file__).with_name("independent_results.json")
    if options.write:
        with destination.open("x", encoding="utf-8") as stream:
            json.dump(record, stream, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        verify(record, json.loads(destination.read_text(encoding="utf-8")))
    print("PASS", record["multiplicities"], record["path_spin_records"])
