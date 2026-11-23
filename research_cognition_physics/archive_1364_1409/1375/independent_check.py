"""Independent generators, symmetric-six-spin source, and scalar kinetic readouts."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def gauge_mass_coefficients(two_t, two_y):
    # Remove the gauge couplings. Complex neutral VEV = v / sqrt(2).
    t, y = two_t/2, two_y/2
    m = np.arange(-t, t+1)
    ladder = np.diag(np.sqrt((t-m[:-1]) * (t+m[:-1]+1)), -1)
    t1 = (ladder + ladder.T)/2
    idx, = np.flatnonzero(m == -y)
    return 4*(t1 @ t1)[idx, idx], 4*m[idx]**2


def locking_tensor():
    # Six spins 1/2. J=3 is the completely symmetric representation.
    words = list(itertools.product((0, 1), repeat=6))
    sym = np.zeros((64, 7))
    for i, word in enumerate(words):
        k = sum(word)
        sym[i, k] = 1 / math.sqrt(math.comb(6, k))
    assert np.allclose(sym.T @ sym, np.eye(7))
    def product(v):
        conjugate = np.array([v[1].conjugate(), -v[0].conjugate()])
        state = conjugate
        for _ in range(5):
            state = np.kron(v, state)
        return state
    h = np.array([0., 1.])
    projected = sym.T @ product(h)
    assert np.flatnonzero(projected).tolist() == [5]  # m=3-5=-2
    c, s = np.cos(.31), np.sin(.31)
    g = np.array([[c*np.exp(.2j), s*np.exp(.7j)], [-s*np.exp(-.7j), c*np.exp(-.2j)]])
    assert np.allclose(g.conj().T @ g, np.eye(2)) and np.isclose(np.linalg.det(g), 1)
    six = np.ones((1, 1), dtype=complex)
    for _ in range(6):
        six = np.kron(six, g)
    assert np.allclose(sym.T @ product(g @ h), sym.T @ six @ sym @ projected)
    return float(projected[5])


def scalar_model(y, sign):
    x = math.sqrt(1 - y*y)
    values = np.array([x, y/4, 0.])
    h = np.array([x, y/4, sign*math.sqrt(15)*y/4])
    assert np.isclose(np.linalg.norm(h), 1)
    coeff = np.array([gauge_mass_coefficients(1, 1), gauge_mass_coefficients(6, 4), [0., 0.]])
    np.testing.assert_allclose(coeff, [[1, 1], [16, 16], [0, 0]])
    gauge_mass = lambda fields: fields**2 @ coeff
    # A unit finite step exactly extracts the derivatives of the quadratic mass functions.
    zero, plus, minus = [gauge_mass(values + step*h) for step in (0, 1, -1)]
    np.testing.assert_allclose(zero, [1, 1])
    np.testing.assert_allclose((plus-minus)/4, [1, 1])
    np.testing.assert_allclose((plus+minus-2*zero)/2, [1, 1])
    assert np.isclose(h[0]/x, 1)
    v_current = values * coeff[:, 0]
    remainder = v_current - h
    norm = np.linalg.norm(remainder)
    axis = remainder/norm
    other = np.cross(h, axis)
    hidden1 = (axis + other)/math.sqrt(2)
    hidden2 = (axis - other)/math.sqrt(2)
    frame = np.column_stack((h, hidden1, hidden2))
    np.testing.assert_allclose(frame.T @ frame, np.eye(3), atol=1e-14)
    powers = (frame.T @ v_current)**2
    np.testing.assert_allclose(powers, [1, 7.5*y*y, 7.5*y*y])
    # Small positive locking term; the neutral CP-even invariant Jacobian
    # can supply the rest of a chosen strictly positive mass Hessian.
    target = frame @ np.diag([1., 2., 4.]) @ frame.T
    lock_direction = np.array([-6*values[1]/x, 1., 0.])
    lock = .03*np.outer(lock_direction, lock_direction)
    jacobian = np.diag([x, values[1], 1.])   # H†H, X†X, singlet
    inv = np.linalg.inv(jacobian)
    k = inv.T @ (target-lock) @ inv
    assert min(np.linalg.eigvalsh(k)) > 0
    np.testing.assert_allclose(lock + jacobian.T @ k @ jacobian, target)
    return [float(q) for q in powers]


def certificate():
    # Brute force independent check of the irreducible-complex-multiplet Pell sequence.
    hits = []
    for twice_t in range(1, 100):
        for twice_y in range(1, twice_t + 1):
            if (twice_t - twice_y) % 2 == 0 and twice_t*(twice_t+2) == 3*twice_y**2:
                hits.append([twice_t, twice_y])
    assert hits == [[1, 1], [6, 4], [25, 15], [96, 56]]
    for sign in (-1, 1):
        for y in (.1, .35, .6):
            powers = scalar_model(y, sign)
            assert abs(sum(powers[1:]) - 15*y*y) < 1e-12
    return {
        "twice_T_abs_twice_Y_under_100": hits,
        "J3_m_minus2_source_amplitude": round(locking_tensor(), 10),
        "y_035_split_neutral_kappa_V_squared": np.round(scalar_model(.35, 1), 10).tolist(),
        "y_035_if_concentrated_max_extra_abs_kappa_V": round(math.sqrt(15)*.35, 10),
        "verified_vev_fractions": [.1, .35, .6],
        "checked_both_singlet_signs": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = certificate()
    path = Path(__file__).with_name("independent_results.json")
    if args.write:
        with path.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(result, indent=2) + "\n")
    else:
        saved = json.loads(path.read_text(encoding="utf-8"))
        assert saved.keys() == result.keys()
        for key, value in result.items():
            assert np.shape(value) == np.shape(saved[key])
            np.testing.assert_allclose(value, saved[key], rtol=0, atol=2e-10)
    print(json.dumps(result))
