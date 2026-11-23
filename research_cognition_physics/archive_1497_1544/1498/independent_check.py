"""Single flavour-current checks with tree angles and an archival Wilson-coefficient box."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
PAIRS = [(0, 1), (0, 2), (1, 2)]


def read_prior(number):
    return json.loads((BASE.parent / str(number) / "inputs.json").read_text(encoding="utf-8-sig"))


def two_widths(data):
    return data[0] - 2*sum(data[1:]), data[0] + 2*sum(data[1:])


def tree_data():
    charged = read_prior(1381)
    rate = charged["K_product"]
    form = charged["K_form"]
    cb = read_prior(1412)["BelleII_Vcb_BGL"]
    ratio = charged["ub_over_cb"]
    phase = read_prior(1490)["chargedB"]
    gamma = phase["central_degrees"]
    us = rate[0]/form[0]
    x, y = two_widths(rate)
    a, b = two_widths(form)
    bounds = ((x/b, y/a), two_widths(cb), two_widths(ratio),
              (gamma-2*phase["lower_error"]-phase["additional_degrees"],
               gamma+2*phase["upper_error"]+phase["additional_degrees"]))
    return [us, cb[0], ratio[0], gamma], bounds


def ckm(us, cb, ratio, gamma):
    s13 = cb*ratio
    c13 = math.sqrt(1-s13*s13)
    s12, s23 = us/c13, cb/c13
    c12, c23 = math.sqrt(1-s12*s12), math.sqrt(1-s23*s23)
    phase = math.radians(gamma)
    delta = phase + math.asin(c12*s23*s13*math.sin(phase)/(s12*c23))
    e = np.exp(1j*delta)
    return np.array([[c12*c13, us, s13*e.conjugate()],
                     [-s12*c23-c12*s23*s13*e, c12*c23-s12*s23*s13*e, cb],
                     [s12*s23-c12*c23*s13*e, -c12*s23-s12*c23*s13*e, c23*c13]])


def envelope(bounds):
    (ul, uh), (cl, ch), (rl, rh), (gl, gh) = bounds
    blo, bhi = cl*rl, ch*rh
    c13 = math.sqrt(1-bhi*bhi)
    s12, s23 = uh/c13, ch/c13
    c12, c23 = math.sqrt(1-s12*s12), math.sqrt(1-s23*s23)
    low = np.array([[c12*c13, ul, blo],
                    [ul*c23-s23*bhi, c12*c23-s12*s23*bhi, cl]])
    high = np.array([[1., uh, bhi],
                     [s12+s23*bhi, 1., ch]])
    gamma_sin = min(math.sin(math.radians(t)) for t in (gl, gh))
    j = np.prod(low[:, [0, 2]])*gamma_sin
    h = j/max(np.prod(high, axis=0))
    return low, high, h


def projected_bound(weights, upper, height):
    matrix = np.zeros((3, 3))
    for (i, j), weight in zip(PAIRS, weights):
        matrix[i, j] = matrix[j, i] = weight
    spectral = max(np.sum(matrix, axis=0))
    entry = sum(matrix[i, j]*(upper[0, i]*upper[1, j]+upper[0, j]*upper[1, i])
                for i, j in PAIRS)
    # The operator norm also bounds the same matrix element.
    return height/(weights[3]+min(spectral, entry)+2*spectral*height)


def rank_one_certificate(weights, low, high):
    """Invert necessary inequalities, one possible largest down component at a time."""
    matrix = np.zeros((3, 3))
    for (i, j), weight in zip(PAIRS, weights):
        matrix[i, j] = matrix[j, i] = weight

    def possible(scale):
        for k in range(3):
            ck = math.sqrt(max(1/3, 1-3*scale**2*sum(matrix[k]**2)))
            orth = (high @ matrix[k])*scale/ck
            if np.prod(np.maximum(0, low[:, k]*ck-orth)) <= weights[3]*scale:
                return True
        return False

    a, b = 0., 1.
    while not possible(b):
        b *= 2
    for _ in range(65):
        mid = (a+b)/2
        if possible(mid):
            b = mid
        else:
            a = mid
    return a


def diagonal_minimum(z):
    trial = []
    for i, j, k in itertools.permutations(range(3)):
        x = min(1., max(0., -(z[i]*z[k].conjugate()).real/abs(z[k])**2))
        q = np.zeros(3)
        q[i], q[k] = 1, x
        trial.append((abs(z@q), q.tolist()))
    return min(trial)


def audit():
    central, bounds = tree_data()
    v = ckm(*central)
    assert np.linalg.norm(v@v.conj().T - np.eye(3)) < 1e-13
    low, high, height_lo = envelope(bounds)
    z = v[0]*v[1].conjugate()
    j = abs((z[0]*z[1].conjugate()).imag)
    height = j/max(abs(z))
    assert math.isclose(j, abs(z[0]*z[2])*math.sin(math.radians(central[3])), rel_tol=1e-12)
    min_value, minimizing_q = diagonal_minimum(z)
    assert math.isclose(height, min_value, rel_tol=1e-12)
    caps = [math.hypot(9.6e-13, 4.4e-15), 2.3e-11, 1.1e-9, 7.2e-13]
    weights = np.sqrt(2*np.array(caps))
    # The two flavour-diagonal commuting generators give phases differing by pi/2 in uc.
    lhs = np.array([z[:2].real, z[:2].imag])
    p, r = np.linalg.solve(lhs, [(-z[2].imag), z[2].real])
    q1 = np.array([0., 0., 1.])
    q2 = np.array([p, r, 0.])
    q2 -= min(q2)
    m1, m2 = z@q1, z@q2
    assert abs(m2-1j*m1) < 5e-15 and abs(m1*m1 + m2*m2) < 4e-17
    assert abs(m1)**2+abs(m2)**2 > 0
    # Test the rank-one and general estimates against actual matrices at box endpoints.
    rng = np.random.default_rng(1498)
    for point in itertools.product(*bounds):
        for gamma_shift in (0, 180):
            matrix = ckm(*point[:3], point[3]+gamma_shift)
            true = abs(matrix[:2])
            assert np.all(true >= low-2e-15) and np.all(true <= high+2e-15)
            sides = matrix[0]*matrix[1].conj()
            assert abs((sides[0]*sides[1].conj()).imag)/max(abs(sides)) >= height_lo-1e-15
            x = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
            q = x+x.conj().T
            q /= np.ptp(np.linalg.eigvalsh(q))
            currents = [abs(q[i, j0]) for i, j0 in PAIRS]+[abs((matrix@q@matrix.conj().T)[0, 1])]
            assert max(currents/weights) >= projected_bound(weights, high, height_lo)
    return {
        "channels": ["K0", "Bd", "Bs", "D0"],
        "tree_central_us_cb_ratio_gamma_deg": central, "tree_boxes": bounds,
        "tree_j_abs_central": j, "triangle_height_central": height,
        "full_tree_box_height_floor": height_lo, "central_saturating_diagonal": minimizing_q,
        "unit_spectral_span_unweighted_floor_full_box": projected_bound(np.ones(4), high, height_lo),
        "phase_free_rank_one_unweighted_floor": rank_one_certificate(np.ones(4), low, high),
        "archival_UTfit_table_IV_LL_caps_GeV_neg2": caps,
        "conditional_scale_F0_neg_half_GeV": {
            "rank_one_full_box_no_phase": rank_one_certificate(weights, low, high),
            "rank_one_central": rank_one_certificate(weights, abs(v[:2]), abs(v[:2])),
            "general_hermitian_unit_span_full_box_with_gamma": projected_bound(weights, high, height_lo),
            "general_hermitian_unit_span_central": projected_bound(weights, abs(v[:2]), height)},
        "two_distinct_positive_residue_currents": {
            "charge_1": q1.tolist(), "charge_2": q2.tolist(),
            "single_squared_coefficient_magnitude": abs(m1)**2,
            "combined_squared_coefficient_magnitude": abs(m1*m1+m2*m2),
            "untwisted_flavour_preserving_exchange": abs(m1)**2+abs(m2)**2}}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
