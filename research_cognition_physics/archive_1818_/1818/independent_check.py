"""Check the public two-round instrument and an energy-conserving feedback reset."""

import argparse
import hashlib
import itertools
import json
import math
import runpy
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE.parents[2] / "猜想" / "知行闭环的分化与物质生成猜想" / "闭环分化验证"
RESULT = HERE / "independent_results.json"
I = np.eye(2)
E = np.diag([0.0, 1.0])
GROUND = I - E
SWAP = np.eye(4)[[0, 2, 1, 3]]


def close(a, b):
    np.testing.assert_allclose(a, b, atol=2e-12, rtol=0)


def measured_kraus(old, basis, phase):
    """Obtain the two conditional W maps by running the original two-reader gate."""
    output = np.zeros((2, 2, 2), dtype=complex)
    for column in range(2):
        source = np.zeros((2, 2, 2, 2), dtype=complex)
        source[basis, column, 0, 0] = 1
        full = old["issue_reading"](source, phase, 0, (2, 3))
        for r in range(2):
            output[r, :, column] = full[basis, :, r, r]
            close(full[basis, :, r, 1-r], 0)
    close(sum(k.conj().T @ k for k in output), I)
    return output


def reduced(rho, keep, dims):
    arr = rho.reshape(tuple(dims) * 2)
    size = len(dims)
    for i in reversed(range(len(dims))):
        if i not in keep:
            arr = np.trace(arr, axis1=i, axis2=i + size)
            size -= 1
    product = math.prod(dims[i] for i in keep)
    return arr.reshape(product, product)


def entropy(rho):
    eig = np.linalg.eigvalsh(rho)
    eig = eig[eig > 1e-14]
    return float(-np.dot(eig, np.log(eig)))


def reset(angle, r, excited):
    # Port order W, B (finite thermal bath), A (ready buffer).
    bath = np.diag([1-excited, excited])
    partial_swap = math.cos(angle) * np.eye(4) - 1j * math.sin(angle) * SWAP
    swap_wa = np.eye(8).reshape((2,)*6).transpose(2, 1, 0, 3, 4, 5).reshape(8, 8)
    u = swap_wa @ np.kron(partial_swap, I)
    energy = np.kron(E, np.eye(4)) + np.kron(I, np.kron(E, I) + np.kron(I, E))
    close(energy @ u, u @ energy)
    initial = np.kron(np.diag(I[r]), np.kron(bath, GROUND))
    out = u @ initial @ u.conj().T
    w, b, a = [reduced(out, [j], (2, 2, 2)) for j in range(3)]
    close(w, GROUND)
    qb = float(b[1, 1].real - excited)
    ea = float(a[1, 1].real)
    close(qb, math.sin(angle)**2 * (r-excited))
    close(qb + ea, r)
    beta = math.log((1-excited)/excited)
    delta_entropy = entropy(b) - entropy(bath)
    relative_entropy = -entropy(b) - np.trace(b @ np.diag(np.log(np.diag(bath)))).real
    close(beta * qb, delta_entropy + relative_entropy)
    mutual_information = entropy(a) + entropy(b) - entropy(out)
    close(mutual_information, entropy(a) + delta_entropy)
    return out, dict(heat_into_bath_gap_units=qb, energy_into_buffer_gap_units=ea,
                     bath_entropy_change_nats=delta_entropy,
                     bath_relative_entropy_nats=relative_entropy,
                     bath_buffer_mutual_information_nats=mutual_information)


def record_test(old):
    old["run"]()
    histories = 0
    # Four input matrix units span all W operators, hence test the complete instrument
    # including an arbitrary untouched reference tensored with this map.
    phase1, phase2 = 0.317, -0.481
    angles = (0.0, math.pi/4, math.pi/2)
    for basis in range(2):
        first = measured_kraus(old, basis, phase1)
        observable = old["X"] if basis == 0 else old["Y"]
        first_zero = measured_kraus(old, basis, 0.0)
        close(first_zero[0].conj().T @ first_zero[0], (I + observable)/2)
        for r, t in itertools.product(range(2), repeat=2):
            next_read = measured_kraus(old, r, phase2)[t]
            original = next_read @ old["H"] @ (old["X"] if r else I) @ first[r]
            row = first[r][r, :]
            for angle in angles:
                complete, _ = reset(angle, r, 0.25)
                output = reduced(complete, [0], (2, 2, 2))
                later = next_read @ old["H"]
                for i, j in itertools.product(range(2), repeat=2):
                    unit = np.outer(I[i], I[j])
                    expected = original @ unit @ original.conj().T
                    actual = (row @ unit @ row.conj()) * (later @ output @ later.conj().T)
                    close(actual, expected)
                histories += 1
    return histories


def ladder(dim, phase):
    """Invariant signed effect on probe and finite equally spaced reference."""
    mode = np.exp(-1j * phase)
    reference = np.sqrt(2/(dim+1)) * np.sin(np.pi*np.arange(1, dim+1)/(dim+1))
    obs = np.zeros((2*dim, 2*dim), dtype=complex)
    for n in range(dim-1):
        obs[n+1, dim+n] = mode
        obs[dim+n, n+1] = mode.conjugate()
    h = np.diag(np.concatenate([np.arange(dim), 1 + np.arange(dim)]))
    close(obs @ h, h @ obs)
    assert np.linalg.norm(obs, 2) <= 1 + 2e-12
    lift = np.kron(I, reference[:, None])
    target = np.array([[0, mode], [mode.conjugate(), 0]])
    contrast = math.cos(math.pi/(dim+1))
    close(lift.conj().T @ obs @ lift, contrast * target)
    path = np.diag(np.ones(dim-1), 1) + np.diag(np.ones(dim-1), -1)
    close(np.linalg.eigvalsh(path/2)[-1], contrast)
    return (1-contrast)/2


def rounded(value):
    if isinstance(value, dict):
        return {key: rounded(item) for key, item in value.items()}
    if isinstance(value, list):
        return [rounded(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return round(float(value), 12) + 0.0
    return value


def run():
    sys.path.insert(0, str(EXPERIMENT))
    try:
        old = runpy.run_path(str(EXPERIMENT / "public_choice_quantum_return.py"))
    finally:
        sys.path.pop(0)
    checked = record_test(old)
    ledger = {}
    for name, angle in (("no_bath_swap", 0), ("half_probability", math.pi/4),
                        ("full_bath_swap", math.pi/2)):
        ledger[name] = {str(r): reset(angle, r, 0.25)[1] for r in range(2)}
    efforts = {}
    for dim in (1, 2, 4, 8, 16):
        x, y = [ladder(dim, phi) for phi in (0.0, math.pi/2)]
        close(x, y)
        efforts[str(dim)] = x
    return rounded(dict(source_sha256=hashlib.sha256(
        (EXPERIMENT / "public_choice_quantum_return.py").read_bytes()).hexdigest(),
        beta_times_gap=math.log(3), initial_thermal_population=0.25,
        conditional_maps_checked_on_four_matrix_units=checked,
        reservoir_and_buffer=ledger,
        optimal_single_binary_readout_error_by_ladder_size=efforts))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    payload = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.check:
        if result != json.loads(RESULT.read_text(encoding="utf-8")):
            raise SystemExit("Saved results differ from recomputation.")
    else:
        RESULT.write_text(payload, encoding="utf-8")
    print("1818 independent", hashlib.sha256(payload.encode()).hexdigest()[:12])
