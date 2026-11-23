"""Four-leg preparations and a six-permutation readout, round 1353."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def generators(group):
    eye = np.eye(3, dtype=complex)
    pairs = [(0, 1), (0, 2), (1, 2)]
    imaginary = [-1j*(np.outer(eye[i], eye[j]) - np.outer(eye[j], eye[i]))
                 for i, j in pairs]
    if group == "SO3":
        return imaginary
    real = [np.outer(eye[i], eye[j]) + np.outer(eye[j], eye[i]) for i, j in pairs]
    diagonal = [np.diag([1, -1, 0]), np.diag([1, 1, -2])/math.sqrt(3)]
    return [x/2 for x in real + imaginary + diagonal]


def readout(vectors):
    """Expand the actual six permutation matrix elements without an 81^3 array."""
    gram = np.array([[np.vdot(u, v) for v in vectors] for u in vectors])
    terms = []
    for p in itertools.permutations(range(3)):
        inversions = sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3))
        terms.append((-1)**inversions * np.prod([gram[i, p[i]] for i in range(3)]))
    result = sum(terms)/6
    assert abs(result.imag) < 1.e-12
    assert abs(result - np.linalg.det(gram)/6) < 1.e-12
    return max(0., float(result.real))


def packet(operator):
    # Basis ordered as the nine active pairs followed by their nine partners.
    return operator.reshape(-1)/3


def source(group):
    local = generators(group)
    cr = np.trace(sum(t@t for t in local)).real/3
    total = [(np.kron(t, np.eye(3)) + np.kron(np.eye(3), t))/np.sqrt(cr) for t in local]
    casimir = sum(t @ t for t in total)
    d = 1 + len(total)
    parent = np.zeros((9*d, 9*d), dtype=complex)
    for index, charge in enumerate(total, 1):
        parent[:9, 9*index:9*(index+1)] = charge
        parent[9*index:9*(index+1), :9] = charge
    return parent, casimir, d


def setting(group, theta, extra):
    parent, casimir, d = source(group)
    drift = np.kron(np.eye(d), casimir)
    assert np.linalg.norm(drift@parent - parent@drift) < 1.e-12
    eigenvalues, basis = np.linalg.eigh(parent + extra*drift)
    # All source outcomes, with two partners retained via vectorization.
    columns = (basis * np.exp(-1j*theta*eigenvalues)) @ basis[:9, :].conj().T
    outputs = [packet(columns[9*i:9*(i+1)]) for i in range(d)]
    assert abs(sum(np.vdot(x, x).real for x in outputs)-1) < 1.e-12
    swap = np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9)
    p_return = float(np.vdot(outputs[0], outputs[0]).real)
    p_click = readout([packet(np.eye(9)), packet(swap), outputs[0]])
    ideal = 0.
    if group == "SO3":
        f2 = math.cos(math.sqrt(3)*theta)*np.exp(-3j*theta*extra)
        ideal = 10*abs(1-f2)**2/729
    assert abs(p_click-ideal) < 1.e-12
    branches = {"source_returns_antisymmetric": p_click,
                "source_returns_complement": max(0., p_return-p_click),
                "source_excited": max(0., 1-p_return)}
    assert abs(sum(branches.values())-1) < 1.e-12
    return branches


def calculate():
    times = [math.pi/math.sqrt(3), math.pi/(2*math.sqrt(3))]
    output = {"dimensionless_times": times, "rates": {}, "alternative_rates": {}}
    for group in ("SU3", "SO3"):
        a = [setting(group, t, 0.) for t in times]
        b = [setting(group, t, 1/math.sqrt(3)) for t in times]
        for ideal, shifted in zip(a, b):
            assert abs(ideal["source_excited"]-shifted["source_excited"]) < 1.e-12
        output["rates"][group] = a
        output["alternative_rates"][group] = b
    # On two active legs the scalar centre commutes with every generator;
    # on four legs it cancels. The same physical source matrices then apply.
    zeta = np.exp(2j*math.pi/3)
    central_two = (zeta**2)*np.eye(9)
    for t in generators("SO3"):
        charge = np.kron(t, np.eye(3)) + np.kron(np.eye(3), t)
        assert np.linalg.norm(central_two@charge - charge@central_two) < 1.e-12
    assert abs(zeta**2*zeta.conjugate()**2 - 1) < 1.e-12

    identity = packet(np.eye(9))
    swap = packet(np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9))
    zero = np.array([1., 0.])
    tag = np.sqrt([17/27, 10/27])
    false_positive = readout([np.kron(identity, zero), np.kron(swap, zero),
                              np.kron(identity, tag)])
    assert abs(false_positive - 40/729) < 1.e-12
    p_min = 10/729
    systematic = .002
    margin = p_min/2-systematic
    repetitions = math.ceil(math.log(400)/(2*margin**2))
    output["controls"] = {"same_colour_extra_tag_false_positive": false_positive,
                          "mu3_SO3_same_as_SO3_for_both_readouts": True,
                          "absolute_systematic_budget": systematic,
                          "repetitions_per_time_four_tail_bound_0p01": repetitions,
                          "four_tail_bound": 4*math.exp(-2*repetitions*margin**2)}
    return output


def compare(current, saved):
    if isinstance(current, dict):
        assert current.keys() == saved.keys()
        for name in current:
            compare(current[name], saved[name])
    elif isinstance(current, list):
        assert len(current) == len(saved)
        for x, y in zip(current, saved):
            compare(x, y)
    else:
        assert np.isclose(current, saved, rtol=1.e-11, atol=1.e-12)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = calculate()
    path = Path(__file__).with_name("results.json")
    if args.write:
        with path.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        compare(result, json.loads(path.read_text(encoding="utf8")))
    print(json.dumps(result, indent=2))
