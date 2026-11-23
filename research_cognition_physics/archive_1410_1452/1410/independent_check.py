"""Test cycle spectral projectors, invariant metric, missing ports and multiplicity."""
import argparse
import json
import math
from collections import deque
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
FILE = HERE / "independent_results.json"


def action(m):
    identity = np.eye(m + 1)
    cycle = identity[:, [*range(1, m), 0, m]]
    swap = identity[:, [m, *range(1, m), 0]]
    common = np.ones_like(identity) / (m + 1)
    return cycle, swap, identity - common


def closure(m):
    start = tuple(range(m + 1))
    generators = (tuple([*range(1, m), 0, m]), tuple([m, *range(1, m), 0]))
    todo, elements = deque([start]), {start}
    while todo:
        current = todo.popleft()
        for perm in generators:
            candidate = tuple(perm[i] for i in current)
            if candidate not in elements:
                elements.add(candidate)
                todo.append(candidate)
    return elements


def cycle_read(m):
    cycle, swap, contrast = action(m)
    eye = np.eye(m + 1)
    zero = (eye - swap) @ eye[:, -1] / math.sqrt(2)
    powers = [np.linalg.matrix_power(cycle, j) for j in range(m)]
    moments = [np.vdot(zero, gate @ zero) for gate in powers]
    assert np.allclose(moments, [1] + [.5]*(m-1))
    projectors = []
    for k in range(m):
        coefficients = np.exp(-2j*np.pi*k*np.arange(m)/m) / m
        proj = sum(c * element for c, element in zip(coefficients, powers)) @ contrast
        assert np.allclose(proj @ proj, proj) and np.allclose(proj, proj.conj().T)
        assert np.isclose(np.trace(proj), 1)
        projectors.append(proj)
    observed = np.array([np.vdot(zero, p @ zero).real for p in projectors])
    expected = [((m+1)/(2*m))] + [1/(2*m)]*(m-1)
    assert np.allclose(observed, expected)
    return np.array(projectors), zero, observed


def arbitrary_complex_fringe():
    # An actual qutrit with the same cyclic spectrum and a different pure prepared ray.
    p = np.array([.5, .38, .12])
    phases = np.exp(2j*np.pi*np.arange(3)/3)
    vector = np.sqrt(p)*np.exp(1j*np.array([.2, .1, -.4]))
    gate = np.diag(phases)
    reflection = np.eye(3) - 2*np.outer(vector, vector.conj())
    fringe = np.vdot(vector, gate @ vector)
    assert abs(np.trace(gate @ reflection) + 2*fringe) < 1e-12
    read = [(1+2*fringe.real)/3,
            (1-fringe.real)/3+fringe.imag/math.sqrt(3),
            (1-fringe.real)/3-fringe.imag/math.sqrt(3)]
    assert np.allclose(read, p)
    return {"complex_fringe": [fringe.real, fringe.imag], "three_probabilities": read}


def invariant_metric():
    m = 3
    group = closure(m)
    trial = np.arange(1, 17).reshape(4, 4) + 1j*np.diag([1., 3., 4., 9.])
    positive = trial.conj().T @ trial + np.eye(4)
    averaged = sum(positive[np.ix_(g, g)] for g in group) / len(group)
    a = averaged[0, 0].real - averaged[0, 1].real
    b = averaged[0, 1].real
    assert a > 0 and a + 4*b > 0
    assert np.allclose(averaged, a*np.eye(4) + b*np.ones((4, 4)))
    cycle, swap, contrast = action(m)
    unequal = np.diag([1., 1., 1., 4.])
    assert np.allclose(cycle.T @ unequal @ cycle, unequal)
    assert np.allclose(contrast @ averaged @ contrast, a*contrast)
    return {"full_group_example_metric_a": a, "full_group_example_metric_b": b,
            "unequal_weight_swap_metric_defect": np.linalg.norm(swap.T @ unequal @ swap - unequal, 2)}


def two_copy_read(m, efficiencies):
    projectors, zero, ideal = cycle_read(m)
    state = np.kron(zero, [1., 0.])
    observed = []
    for proj, eff in zip(projectors, efficiencies):
        detector_copy = np.sqrt([eff, 1-eff])
        detector = np.kron(proj, np.outer(detector_copy, detector_copy))
        observed.append(np.vdot(state, detector @ state).real)
    assert np.allclose(observed, ideal * efficiencies)
    return np.array(observed)


def disguised(m):
    # One registered readout per distinct phase; missing phases give no light signal.
    efficiencies = np.zeros(m)
    efficiencies[:3] = [4/(m+1), 1., 1.]
    read = two_copy_read(m, efficiencies)
    total = sum(read)
    assert np.isclose(total, 3/m)
    conditional = read[:3] / total
    assert np.allclose(conditional, [2/3, 1/6, 1/6])
    return {"accepted_rate": total, "missing_rate": 1-total,
            "three_conditionals": list(conditional)}


def boundary(loss, channel):
    _, _, ideal = cycle_read(3)
    p = ideal[channel]
    low = max(0., (p-loss)/(1-loss))
    high = min(1., p/(1-loss))
    for target in (low, high):
        accepted = ideal.copy()
        accepted[channel] = (1-loss)*target
        others = [j for j in range(3) if j != channel]
        factor = (1-loss-accepted[channel]) / sum(ideal[others])
        accepted[others] *= factor
        checked = two_copy_read(3, np.clip(accepted/ideal, 0., 1.))
        assert np.isclose(sum(checked), 1-loss)
        assert np.isclose(checked[channel]/sum(checked), target)
    return [low, high]


def gate_errors():
    m = 3
    projectors, odd, _ = cycle_read(m)
    cycle, swap, _ = action(m)
    spectator = np.ones(m+1)/math.sqrt(m+1)
    # Both test states stay inside contrast; finite eigenvector residuals are measurable.
    fixed = projectors[0][:, -1].real
    fixed = fixed / np.linalg.norm(fixed)
    even = fixed - np.dot(fixed, odd)*odd
    even /= np.linalg.norm(even)
    second = projectors[1][:, 0]
    second = second / np.linalg.norm(second)
    phi, theta = .11, .07
    prepared = math.cos(phi)*odd + math.sin(phi)*even
    readout = math.cos(theta)*fixed + math.sin(theta)*second
    assert abs(np.vdot(spectator, prepared)) + abs(np.vdot(spectator, readout)) < 1e-12
    odd_res = np.linalg.norm(swap @ prepared + prepared)/2
    cycle_res = np.linalg.norm(cycle @ readout-readout)
    gap = 2*math.sin(math.pi/m)
    measured = abs(np.vdot(readout, prepared))**2
    angle_shift = abs(math.acos(math.sqrt(measured))-math.acos(math.sqrt((m+1)/(2*m))))
    bound = math.asin(odd_res) + math.asin(cycle_res/gap)
    assert math.isclose(bound, phi+theta, abs_tol=1e-12)
    assert angle_shift <= bound
    return {"prepared_even_norm": float(odd_res), "cycle_eigen_residual": float(cycle_res),
            "observed_probability": float(measured), "angle_shift_rad": angle_shift,
            "angle_bound_rad": bound}


def run():
    dimensions = (2, 3, 4, 5, 9)
    output = {str(m): list(cycle_read(m)[2]) for m in dimensions}
    groups = {str(m): len(closure(m)) for m in (2, 3, 4)}
    assert all(count == math.factorial(int(m)+1) for m, count in groups.items())
    # An allowed bare odd Majorana block under the transposition has nonzero mass.
    _, z, _ = cycle_read(3)
    _, swap, contrast = action(3)
    odd = np.outer(z, z)
    even = contrast - odd
    massive = 2*odd + 3*even
    assert np.allclose(swap.T @ massive @ swap, massive)
    assert np.allclose(massive @ z, 2*z)
    return {"cycle_projectors": output, "permutation_orders": groups, "other_reflection": arbitrary_complex_fringe(),
            "metric": invariant_metric(), "finite_gate_residuals": gate_errors(),
            "conditional_e_at_1pct_loss": boundary(.01, 0),
            "conditional_mu_at_25pct_loss": boundary(.25, 1),
            "complete_fake_with_exact_target_column": {str(m): disguised(m) for m in (4, 5, 9)},
            "odd_mass_allowed_by_swap_alone": float(np.vdot(z, massive @ z).real)}


def matches(a, b):
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(matches(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(matches(x, y) for x, y in zip(a, b))
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=2e-12)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flags = parser.add_mutually_exclusive_group()
    flags.add_argument("--write", action="store_true")
    flags.add_argument("--check", action="store_true")
    data = run()
    if parser.parse_args().write:
        FILE.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    else:
        assert matches(data, json.loads(FILE.read_text(encoding="utf-8")))
        print("1410 cycle orders, invariant metric, sharp lost-strength bounds and countermodels verified")
