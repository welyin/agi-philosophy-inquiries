"""Independent neutral-pointer verification using explicit matter generators and pointer H."""
import argparse
import json
from pathlib import Path

import numpy as np


def su3():
    matrices = []
    for j, k in [(0, 1), (0, 2), (1, 2)]:
        sym = np.zeros((3, 3), complex)
        skew = sym.copy()
        sym[j, k] = sym[k, j] = 1
        skew[j, k], skew[k, j] = -1j, 1j
        matrices.extend([sym, skew])
    matrices.extend([np.diag([1., -1., 0.]), np.diag([1., 1., -2.]) / np.sqrt(3)])
    return np.array(matrices) / 2


def so3():
    rotations = np.zeros((3, 3, 3), complex)
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        rotations[a, b, c], rotations[a, c, b] = -1j, 1j
    return rotations


def pointer_hamiltonian(generators, first, second, gated_phase=False):
    sources = []
    for generator in generators:
        left = np.kron(generator, np.eye(3))
        right = np.kron(np.eye(3), generator)
        second_weight = second
        if gated_phase:
            second_weight += (1j - 1) * first * second
        sources.append(first * left + second_weight * right)
    offdiag = np.concatenate(sources, axis=0)
    size = 9 + offdiag.shape[0]
    hamiltonian = np.zeros((size, size), complex)
    hamiltonian[9:, :9] = offdiag
    hamiltonian[:9, 9:] = offdiag.conj().T
    return hamiltonian, offdiag


def response(H, eigen_casimir):
    energies, directions = np.linalg.eigh(H / np.sqrt(eigen_casimir))
    start = directions.conj().T[:, :9]
    # Each column labels one orthogonal untouched conjugate-partner word in vec(I_9)/3.
    def evaluate(tau):
        purified_output = (directions * np.exp(-1j * tau * energies)) @ start / 3
        return float(np.sum(abs(purified_output[9:, :])**2))
    return evaluate


def intervals(p, k, eta, extra_time, read_error):
    delta = k * np.arcsin(2 * eta) / 2 + extra_time
    radius = np.sqrt(2) * delta + read_error
    return [float(max(0, p - radius)), float(min(1, p + radius))]


def compute():
    expected_casimir = {"SU3": 4 / 3, "SO3": 2.}
    scale_error, extra_time, read_error = 0.0025, 0.001, 0.005
    checks = {}
    for name, generators in [("SU3", su3()), ("SO3", so3())]:
        d = len(generators)
        C = expected_casimir[name]
        assert np.allclose(sum(T @ T for T in generators), C * np.eye(3))
        assert np.allclose(np.einsum("aij,bji->ab", generators, generators), 3 * C / d * np.eye(d))
        h_one, _ = pointer_hamiltonian(generators, 1, 0)
        h_other, _ = pointer_hamiltonian(generators, 0, 1)
        both, L = pointer_hamiltonian(generators, 1, 1)
        counter, L_counter = pointer_hamiltonian(generators, 1, 1, gated_phase=True)
        assert np.allclose(L_counter.conj().T @ L_counter, 2 * C * np.eye(9))
        for i,j,h in [(1,0,h_one),(0,1,h_other)]:
            h_cal, _ = pointer_hamiltonian(generators, i, j, gated_phase=True)
            assert np.array_equal(h_cal, h)
        calc_one, calc_other, calc_both, calc_counter = [response(h, C) for h in [h_one,h_other,both,counter]]
        eig = np.linalg.eigvalsh(L.conj().T @ L / C)
        assert abs(np.mean(eig) - 2) < 1e-12
        fourth = np.mean(eig**2) / 3
        assert abs(fourth - (4 + 4 / d) / 3) < 1e-12
        for tau in [0.137, np.pi/4, np.pi/2, np.pi, 4.23]:
            assert abs(calc_one(tau) - np.sin(tau)**2) < 1e-13
            assert abs(calc_other(tau) - calc_one(tau)) < 1e-13
            ratios = np.maximum(0., eig)
            assert abs(calc_both(tau) - np.mean(np.sin(tau * np.sqrt(ratios))**2)) < 1e-13
            assert abs(calc_counter(tau) - np.sin(np.sqrt(2) * tau)**2) < 1e-13
        windows = {}
        for k in [2, 4]:
            central = calc_both(k * np.pi/4)
            windows[str(k)] = dict(tau_over_pi=k/4, probability=central,
                                   bounded_readout=intervals(central, k, scale_error, extra_time, read_error),
                                   nonadditive_probability=calc_counter(k*np.pi/4))
        # An altered joint source keeps both one-site marginals I/3 but changes the diagonal Casimir law.
        singlet = np.eye(3, dtype=complex).reshape(-1) / np.sqrt(3)
        swap = np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9)
        rho_corr = (np.eye(9) + swap) / 12 if name == "SU3" else np.outer(singlet, singlet.conj())
        assert np.allclose(rho_corr.reshape(3,3,3,3).trace(axis1=1, axis2=3), np.eye(3)/3)
        assert np.allclose(rho_corr.reshape(3,3,3,3).trace(axis1=0, axis2=2), np.eye(3)/3)
        for T in generators:
            total = np.kron(T, np.eye(3)) + np.kron(np.eye(3), T)
            assert np.allclose(total @ rho_corr, rho_corr @ total)
        normalized_gram = L.conj().T @ L / C
        corr_leading = float(np.trace(rho_corr @ normalized_gram).real)
        checks[name] = dict(pointer_matter_dimension=both.shape[0], conjugate_partner_dimension=9,
                            single_probability_pi_over_4=calc_one(np.pi/4),
                            normalized_total_casimir_spectrum=eig.tolist(),
                            fourth_coefficient=float(fourth), correlated_source_quadratic=corr_leading,
                            holdouts=windows)
    assert checks["SO3"]["correlated_source_quadratic"] < 1e-12
    for k in ["2", "4"]:
        assert checks["SU3"]["holdouts"][k]["bounded_readout"][0] > checks["SO3"]["holdouts"][k]["bounded_readout"][1]
    return dict(calibration_probability_error=scale_error, extra_pulse_tau_error=extra_time,
                final_probability_error=read_error, generators=checks)


def compare(actual, stored):
    if isinstance(actual, dict):
        assert actual.keys() == stored.keys()
        for k in actual:
            compare(actual[k], stored[k])
    elif isinstance(actual, list):
        assert len(actual) == len(stored)
        for left, right in zip(actual, stored):
            compare(left, right)
    else:
        assert np.isclose(actual, stored, atol=2e-12, rtol=2e-12), (actual,stored)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    result = compute()
    path = Path(__file__).with_name("independent_results.json")
    if parser.parse_args().check:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
    else:
        with path.open("x", encoding="utf-8") as output:
            json.dump(result, output, indent=2, allow_nan=False)
            output.write("\n")
    print("PASS", result["generators"]["SU3"]["holdouts"]["4"]["probability"],
          result["generators"]["SO3"]["holdouts"]["4"]["probability"])
