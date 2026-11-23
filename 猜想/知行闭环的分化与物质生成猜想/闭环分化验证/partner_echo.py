"""Active joint-correlation test with one fixed local return channel."""

import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])
Y = np.array([[0., -1j], [1j, 0.]])
S = np.diag([1., 1j])


def pair(phase=1., p=.75):
    ket = np.array([sqrt(p), 0., 0., phase * sqrt(1 - p)])
    return np.outer(ket, ket.conj())


def local(state, kraus):
    operators = [np.kron(I, k) for k in kraus]
    return sum(k @ state @ k.conj().T for k in operators)


def mean(state, first, second):
    value = np.trace(state @ np.kron(first, second))
    assert abs(value.imag) < 1e-12
    return value.real.item()


def marginals(state):
    indices = state.reshape(2, 2, 2, 2)
    return indices.trace(axis1=1, axis2=3), indices.trace(axis1=0, axis2=2)


def midpoint_limit(state):
    m, x, z = (mean(state, X, observable) for observable in (I, X, Z))
    return max(abs(m), np.hypot(x, z))


def real_return_checks():
    rng = np.random.default_rng(16017)
    draws = 256
    for _ in range(draws):
        raw = rng.normal(size=(4, 4))
        state = raw @ raw.T
        state /= np.trace(state)
        isometry, _ = np.linalg.qr(rng.normal(size=(6, 2)))
        maps = np.split(isometry, 3)
        assert np.allclose(sum(k.T @ k for k in maps), I)
        measured = mean(local(state, maps), X, X)
        assert abs(measured) <= midpoint_limit(state) + 1e-12
        # Directly test the best admissible binary read on this state.
        m, x, z = (mean(state, X, op) for op in (I, X, Z))
        radius = np.hypot(x, z)
        if abs(m) >= radius:
            best = np.sign(m) * I
        else:
            best = (x * X + z * Z) / radius
        assert np.max(np.abs(np.linalg.eigvalsh(best))) <= 1 + 1e-12
        assert np.isclose(mean(state, X, best), midpoint_limit(state))
    return draws


def check_controls(initial, s):
    bypass = local(initial, [I])
    assert np.isclose(mean(bypass, X, X), s)
    projectors = [np.diag([1., 0.]), np.diag([0., 1.])]
    reprepare_z = local(initial, projectors)
    assert np.isclose(mean(reprepare_z, Z, Z), 1)
    assert np.isclose(mean(reprepare_z, X, X), 0)
    assert np.isclose(np.trace(reprepare_z @ initial), .625)

    # A hidden shared switch chooses I or Z on both legs. No input label is
    # sent, yet dropping the switch from the resource ledger fakes the echo.
    middles = [local(initial, [gate]) for gate in (I, Z)]
    averaged = (middles[0] + middles[1]) / 2
    returned = sum(local(mid, [gate]) for mid, gate in zip(middles, (I, Z))) / 2
    assert np.isclose(midpoint_limit(averaged), 0)
    assert np.isclose(sum(midpoint_limit(mid) for mid in middles) / 2, s)
    assert np.allclose(returned, initial)
    flagged = sum(np.kron(mid, flag) for mid, flag in zip(middles, projectors)) / 2
    difference = flagged - np.kron(averaged, I / 2)
    dependence = np.sum(abs(np.linalg.eigvalsh(difference))) / 2
    assert np.isclose(dependence, s / 2)
    return {"dephasing_final_fidelity": .625, "one_hidden_bit_fakes_echo": True,
            "hidden_register_trace_distance_from_independence": round(float(dependence), 12)}


def check_unknown_extension():
    # Matrix units span all inputs, including the coefficient blocks of
    # arbitrary correlated partners. No assumption about their state is used.
    joint = np.kron(I, S)
    checks = 0
    for row in range(4):
        for column in range(4):
            unit = np.zeros((4, 4), dtype=complex)
            unit[row, column] = 1
            mid = joint @ unit @ joint.conj().T
            assert np.array_equal(joint.conj().T @ mid @ joint, unit)
            checks += 1
    return checks


def run():
    initial = pair()
    mid = local(initial, [S])
    final = local(mid, [S.conj().T])
    expected = np.diag([.75, .25])
    for state in (initial, mid, final):
        for side in marginals(state):
            assert np.allclose(side, expected)
    s = sqrt(3) / 2
    assert np.isclose(midpoint_limit(mid), 0)
    assert np.isclose(mean(mid, X, Y), s)  # Diagnostic outside the X/Z-only test.
    assert np.isclose(mean(final, Z, Z), 1)
    assert np.isclose(mean(final, X, X), s)
    assert np.allclose(final, initial)

    # A real orthogonal rotation genuinely moves the joint pure state while
    # fixing maximally mixed marginals. It moves the nondegenerate ones.
    rot = np.array([[.6, -.8], [.8, .6]])
    rotated_bell = local(pair(p=.5), [rot])
    assert not np.allclose(rotated_bell, pair(p=.5))
    assert all(np.allclose(marginal, I / 2) for marginal in marginals(rotated_bell))
    assert not np.allclose(marginals(local(initial, [rot]))[1], expected)
    return {"real_channel_and_sharp_effect_checks": real_return_checks(),
            "unknown_partner_matrix_units": check_unknown_extension(),
            "midpoint_XI_XX_XZ": [round(mean(mid, X, op), 12) for op in (I, X, Z)],
            "return_ZZ": round(mean(final, Z, Z), 12), "return_XX": round(s, 12),
            "degenerate_real_rotation": True, "controls": check_controls(initial, s)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    options = parser.parse_args()
    results = run()
    target = Path(__file__).with_name("partner_echo_results.json")
    if options.check:
        assert json.loads(target.read_text(encoding="utf-8")) == results
    else:
        with target.open("x", encoding="utf-8") as out:
            out.write(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results))
