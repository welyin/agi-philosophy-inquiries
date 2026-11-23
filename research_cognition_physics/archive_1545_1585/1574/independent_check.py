"""Independent exact charge count, radial coupling witnesses and OPAL spacelike test."""
import argparse
import json
import math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
Q2_LO, Q2_HI = 1.81, 6.07
ALPHA = 1 / 137.035999
NODES, WEIGHTS = np.polynomial.legendre.leggauss(96)
X = (1 + NODES) / 2
P = X * (1 - X)


def charge_trace(weak_dim, y):
    """Trace Q² over left Weyl components, with 6Y=y."""
    return sum((F(y, 6) + F(2*i - weak_dim + 1, 2))**2 for i in range(weak_dim))


def sum_rule():
    multiplets = []
    for y in [-1]*6 + [2]*3:
        multiplets.extend([(2, y), (1, 3-y), (1, -3-y)])
    multiplets.append((2, 0))
    electric = Counter()
    y_classes = [F(0) for _ in range(3)]
    spectrum = Counter()
    for d, y in multiplets:
        y_classes[y % 3] += charge_trace(d, y)
        for i in range(d):
            electric[F(y, 6) + F(2*i-d+1, 2)] += 1
            spectrum[(y % 3, 2*i-d+1)] += 1
    # A whole exchanged torus: Laurent coefficients give nine deficit doublets.
    pos_minus_neg = {j: spectrum[1, j] - spectrum[2, j] for j in (-1, 0, 1)}
    assert pos_minus_neg == {-1: -9, 0: 18, 1: -9}
    assert electric == Counter({-k: v for k, v in electric.items()})
    assert sum(d*y**3 for d, y in multiplets) == 0
    assert sum(d*y for d, y in multiplets) == 0
    assert sum(y for d, y in multiplets if d == 2) == 0
    assert sum(d == 2 for d, y in multiplets) == 10
    assert y_classes == [F(1,2), F(11,2), F(11,2)]
    assert min(charge_trace(2,y) for y in range(-20,21) if y % 3) == F(5,9)
    return {
        "signed_weak_weights": pos_minus_neg,
        "weyl_components": sum(d for d, y in multiplets),
        "sum_Q_squared_by_y_mod_3": [float(s) for s in y_classes],
        "sum_Q_squared_sharp": str(sum(y_classes)),
        "extra_screening_bY_plus_b2_sharp": str(F(2,3)*sum(y_classes)),
        "chiral_spectrum_as_d2_y_multiplicity": [[d, y, n] for (d,y), n in sorted(Counter(multiplets).items())],
    }


def delta_for_dirac(mass, lower_response=False):
    """One-loop on-shell Delta alpha change, unit charge and no hidden interactions."""
    if lower_response:
        value = math.log(Q2_HI/Q2_LO)*Q2_LO*P / (mass**2 + Q2_LO*P)
    else:
        value = np.log1p((Q2_HI-Q2_LO)*P / (mass**2 + Q2_LO*P))
    return float(ALPHA / math.pi * np.dot(WEIGHTS, P*value))


def threshold(cap, lower_response):
    lo, hi = 0.0, 3.0
    for _ in range(65):
        mid = (lo+hi)/2
        if (23/4)*delta_for_dirac(mid,lower_response) > cap:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


def spacelike():
    widths = {}
    for n in (1,2):
        cap = (440+n*(58+43+30)-202)*1e-5
        widths[n] = {
            "extra_end_to_end_budget": cap,
            "common_upper_mass_boundary_GeV_endpoint": threshold(cap, False),
            "common_upper_mass_boundary_GeV_lowest_log_slope": threshold(cap, True),
        }
    light = 23/4*delta_for_dirac(0)
    assert math.isclose(light, (23/4)*ALPHA/(3*math.pi)*math.log(Q2_HI/Q2_LO), rel_tol=1e-12)
    assert widths[2]["common_upper_mass_boundary_GeV_lowest_log_slope"] > 0
    return {
        "OPAL_abs_gevsquared": [Q2_LO, Q2_HI],
        "OPAL_end_to_end": .00440,
        "OPAL_component_errors": [.00058, .00043, .00030],
        "OPAL_subtracted_leptonic_delta": .00202,
        "extra_massless_lower_delta": light,
        "extra_100_GeV_lower_delta_if_all_below_that_mass": (23/4)*delta_for_dirac(100),
        "working_rectangles": widths,
    }


def radial():
    # A chiral doublet component and singlet, plus vectorlike copies of each.
    # Ordered rows D0L,DvL,SvL and columns S0R,SvR,DvR.
    mixing = math.sqrt(3)
    m0 = np.array([[0.,0.,0.],[0.,0.,mixing],[0.,mixing,0.]])
    y = np.eye(3)
    m = m0+y
    svals = np.linalg.svd(m, compute_uv=False)
    assert np.linalg.matrix_rank(m0) == 2
    assert abs(np.trace(np.linalg.solve(m,y))) < 1e-14
    assert svals[-1]**2 <= np.linalg.norm(y, "fro")**2
    # Phi_K(x) = x f_K(x²), finite Taylor jet at x=1: f_K(t) approximates t^-1/2.
    from numpy.polynomial import Polynomial
    x = Polynomial([0.,1.])
    tests = {}
    for order in range(1,6):
        f, c = Polynomial([1.]), F(1)
        for j in range(1, order+1):
            c *= F(1-2*j, 2*j)
            f += float(c)*(x*x-1)**j
        phi = x*f
        jet = [float(phi.deriv(j)(1.)) for j in range(order+1)]
        assert phi(0.) == 0 and jet == [1.]+[0.]*order
        tests[order] = jet
    # Exact trace bound saturated by a protected mass below two invariant masses.
    m0_eq = np.diag([0.,3.,4.])
    y_eq = np.diag([2.,0.,0.])
    assert min(np.linalg.svd(m0_eq+y_eq, compute_uv=False))**2 == np.linalg.norm(y_eq, "fro")**2
    return {
        "invariant_mix_over_radial_entry": mixing,
        "normalized_singular_masses": svals.tolist(),
        "soft_log_det_response": round(float(np.trace(np.linalg.solve(m,y))), 12),
        "all_flavour_coupling_norm_squared": float(np.linalg.norm(y, "fro")**2),
        "covariant_mass_replacement_derivatives": tests,
    }


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write", action="store_true")
    args = parser.parse_args()
    records = {"round": 1574, "charge": sum_rule(), "spacelike": spacelike(), "higgs": radial()}
    output = HERE / "independent_results.json"
    payload = json.dumps(records, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.write:
        with output.open("x", encoding="utf-8") as target:
            target.write(payload)
    elif json.loads(output.read_text(encoding="utf-8")) != json.loads(payload):
        raise SystemExit("independent results mismatch")
    print("1574 independent check passed")


if __name__ == "__main__":
    main()
