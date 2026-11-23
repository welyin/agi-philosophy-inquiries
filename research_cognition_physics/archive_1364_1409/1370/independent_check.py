"""Check absolute-mass observables from triangle area and Takagi witnesses."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np


# NuFIT 6.0 Table 1, lower panel (IC24 with SK); mass-squared units are eV^2.
# Each tuple gives the marginal low endpoint, best fit, high endpoint.
SOLAR = (.275, .308, .345)
SPLITTING = (6.92e-5, 7.49e-5, 8.05e-5)
REACTOR = {"NO": (.02030, .02215, .02388), "IO": (.02060, .02231, .02409)}
ATMOSPHERIC = {"NO": (.002451, .002513, .002578), "IO": (.002421, .002484, .002547)}


def parameters(order, s, t, d, atmospheric):
    offsets = [0, d, atmospheric] if order == "NO" else [atmospheric - d, atmospheric, 0]
    weights = [(1 - s) * (1 - t), s * (1 - t), t]
    return np.array(weights), np.array(offsets)


def zero_interval(weights, offsets):
    # Heron's discriminant = (sum sides^2)^2 - 2 sum sides^4.
    # sides_i^2 = line_i * m0^2 + intercept_i.
    line = weights**2
    intercept = line * offsets
    gram = np.ones((3, 3)) - 2 * np.eye(3)
    poly = [line @ gram @ line, 2 * line @ gram @ intercept, intercept @ gram @ intercept]
    roots = np.roots(poly)
    assert poly[0] < 0
    positive = [np.sqrt(x.real) for x in roots if abs(x.imag) < 1e-13 and x.real > 0]
    return sorted(positive)


def two_phases(sides):
    # Compose sides 1 and 2 to length 3, then reverse their resultant.
    a, b, c = sides
    cos_alpha = (c*c - a*a - b*b) / (2*a*b)
    assert abs(cos_alpha) <= 1
    alpha = np.arccos(cos_alpha)
    beta = np.angle(-(a + b * np.exp(1j * alpha)))
    return np.array([0., alpha, beta])


def witness(weights, offsets, x, phases):
    # An orthonormal electron row is completed using a complex conjugated
    # Householder basis; only this row enters the two measured observables.
    row = np.sqrt(weights)
    vec = np.asarray([1., 0., 0.]) - row
    orthogonal = np.eye(3) - 2 * np.outer(vec, vec) / np.dot(vec, vec)
    u = orthogonal @ np.diag(np.exp(-.5j * phases))
    assert np.allclose(u.conj().T @ u, np.eye(3))
    m = np.sqrt(offsets + x*x)
    target = u.conj() @ np.diag(m) @ u.conj().T
    base = u.conj() @ np.diag(np.sqrt(offsets)) @ u.conj().T
    assert np.linalg.matrix_rank(base, tol=1e-13) == 2
    assert np.allclose(m*m - np.min(m*m), offsets)
    assert abs(np.linalg.norm(target - base, ord=2) - x) < 2e-13
    assert abs((target @ target.conj().T)[0, 0].real - np.dot(weights, m*m)) < 1e-14
    mee = abs(target[0, 0])
    assert abs(mee - abs(np.dot(weights * m, np.exp(1j*phases)))) < 1e-13
    return mee


def run():
    result = {}
    for order in ("NO", "IO"):
        s, t, d, a = SOLAR[1], REACTOR[order][1], SPLITTING[1], ATMOSPHERIC[order][1]
        weights, offsets = parameters(order, s, t, d, a)
        sides = weights * np.sqrt(offsets)
        r = {
            "rank2_mbeta_meV": 1000 * np.sqrt(weights @ offsets),
            "rank2_mee_interval_meV": 1000 * np.array([2 * max(sides) - sum(sides), sum(sides)]),
        }
        for x in (0., .013, .1):
            witness(weights, offsets, x, np.array([0., 1., 2.]))
        ranges = zero_interval(weights, offsets)
        if order == "NO":
            corners = itertools.product(*(pair[::2] for pair in
                                        (SOLAR, REACTOR[order], SPLITTING, ATMOSPHERIC[order])))
            end_pairs = np.array([zero_interval(*parameters(order, *row)) for row in corners])
            x = sum(ranges) / 2
            phase = two_phases(weights * np.sqrt(offsets + x*x))
            assert witness(weights, offsets, x, phase) < 1e-13
            r.update({
                "zero_thresholds_meV": 1000 * np.array(ranges),
                "low_threshold_box_meV": 1000 * np.array([min(end_pairs[:, 0]), max(end_pairs[:, 0])]),
                "high_threshold_box_meV": 1000 * np.array([min(end_pairs[:, 1]), max(end_pairs[:, 1])]),
                "interior_zero_witness_m0_meV": 1000*x,
            })
        else:
            assert not ranges
            # Uniform non-cancellation bound for any m0 >= 0 and *all*
            # independent input choices in the reported marginal box.
            smax, tmax, dmax, amin = SOLAR[-1], REACTOR[order][-1], SPLITTING[-1], ATMOSPHERIC[order][0]
            margin = (1 - tmax) * ((1 - smax) * np.sqrt(1 - dmax/amin) - smax) - tmax
            assert margin > 0
            r["uniform_mee_floor_meV"] = 1000 * margin * np.sqrt(amin)
        result[order] = {name: np.around(value, 8).tolist() for name, value in r.items()}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    actual = run()
    target = Path(__file__).with_name("independent_results.json")
    if args.write:
        with target.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(actual, indent=2) + "\n")
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert actual.keys() == saved.keys()
        for order in actual:
            assert actual[order].keys() == saved[order].keys()
            for key, val in actual[order].items():
                assert np.shape(val) == np.shape(saved[order][key])
                np.testing.assert_allclose(val, saved[order][key], atol=1.1e-8, rtol=0)
    print(json.dumps(actual))
