"""Round 1452: tree minors, threshold sums, and the two Abelian currents."""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
PI = math.pi


def chain(epsilon):
    sites = list(range(-8, 9))
    edge = np.zeros((16, 17))
    for row, child in enumerate(site for site in sites if site):
        parent = child - (1 if child > 0 else -1)
        edge[row, sites.index(parent)] = -epsilon
        edge[row, sites.index(child)] = 1
    ideal = np.array([epsilon**abs(site) for site in sites])
    signed = [(-1)**i * np.linalg.det(np.delete(edge, i, axis=1))
              for i in range(17)]
    assert np.allclose(np.array(signed) / signed[8], ideal, atol=1e-14)
    masses = np.sqrt(np.linalg.eigvalsh(edge @ edge.T))
    assert math.isclose(math.prod(masses), np.linalg.norm(ideal), rel_tol=1e-12)
    return masses, np.linalg.norm(ideal)


def abelian_coefficients():
    # Original light Weyl field plus all the telescoping heavy fields have
    # the charge trace of a neutral-F SM family plus full-vector pairs at k.
    ff = Fraction(0)
    yf = Fraction(0)
    extra_yy = Fraction(0)
    species = [(3, Fraction(-2, 3), (8, 4, 0)),
               (3, Fraction(1, 3), (7, 5, 3)),
               (1, Fraction(1), (8, 5, 3))]
    for size, hyper, powers in species:
        for power in powers:
            for k in range(1, power + 1):
                weight = Fraction(4, 3) * size
                yf += weight * hyper * k
                ff += weight * k*k
                extra_yy += weight * hyper**2
    for charge in range(-8, 9):
        scalar = Fraction(2, 3)  # weak doublet, single complex field
        ff += scalar * charge**2
        yf += scalar * Fraction(1, 2) * charge
        if charge:
            extra_yy += scalar / 4
    ff += Fraction(1, 3)  # charge -1 complex flavon
    return [Fraction(41, 6) + extra_yy, yf, ff]


def inverse_couplings(f, end, epsilon, source):
    # Order is physical Y (not GUT normalized), SU(2), SU(3).
    start, one, two, three = source
    inv = 4 * PI / np.array([0.6 * one**2, two**2, three**2])
    inv -= np.array([41 / 6, -19 / 6, -7.0]) * math.log(end/start)/(2*PI)
    fermion = max(0.0, math.log(end/max(f, start)))
    inv -= np.array([148 / 3, 0, 18]) * fermion/(2*PI)
    masses, norm = chain(epsilon)
    tower = sum(max(0.0, math.log(end/max(f/epsilon*m, start))) for m in masses)
    if start <= min(masses)*f/epsilon and end >= max(masses)*f/epsilon:
        assert math.isclose(tower, 16*math.log(end*epsilon/f) - math.log(norm))
    return inv - np.array([1/6, 1/6, 0.0]) * tower/(2*PI)


def result():
    b_yy, b_yf, b_ff = abelian_coefficients()
    assert (b_yy, b_yf, b_ff) == (Fraction(353, 6), Fraction(56, 3), Fraction(7225, 3))
    # Antusch et al., v2 Table 2, 2024-PDG input; Eq. 4 for the weak boundary.
    table = {1e12: (1e12, .534498, .55157, .6052),
             1e16: (1e16, .575240, .52296, .5310)}
    weak = (91.1876, .461228, .65096, 1.2123)
    families = {}
    for epsilon in (.20, .22):
        masses, norm = chain(epsilon)
        lo, hi = math.log(1e8), math.log(1e15)
        for _ in range(95):
            mid = (lo+hi)/2
            # Enforce alpha <= .30 for each of the three SM group factors.
            worst = min(inverse_couplings(math.exp(mid), 1e16, epsilon, weak))
            if worst >= 1/.30:
                hi = mid
            else:
                lo = mid
        output = {"root_norm": norm, "heavy_singular_extrema": [min(masses), max(masses)],
                  "MZ_to_1e16_all_alpha_le_point30_floor_f_GeV": math.exp(hi)}
        for f, end in [(1e12, 1e16), (1e16, 1e18)]:
            inv = inverse_couplings(f, end, epsilon, table[f])
            output[f"formal_inverse_alpha_Y_2_3_from_{f:g}_to_{end:g}"] = inv.tolist()
        f, end = 1e16, 1e18
        match = 2*f/epsilon
        assert match > max(masses)*f/epsilon
        ell = math.log(end/match)/(8*PI**2)
        yy = inverse_couplings(f, end, epsilon, table[f])[0]/(4*PI)
        cross = -float(b_yf)*ell
        critical = (float(b_ff)*ell + cross**2/yy)**(-.5)
        def kinetic(gauge):
            matrix = [[yy, cross], [cross, gauge**-2-float(b_ff)*ell]]
            return np.linalg.eigvalsh(matrix).tolist()
        assert kinetic(.05)[0] > 0 and kinetic(.15)[0] < 0
        output["upper_segment_zero_cross_matching"] = {
            "Q0_GeV": match, "gF_positive_kinetic_upper_strict": critical,
            "eigenvalues_for_gF_point05": kinetic(.05),
            "formal_eigenvalues_for_gF_point15": kinetic(.15)}
        families[str(epsilon)] = output
    return {"basis_physical_Y_F_beta_matrix_entries": [str(b_yy), str(b_yf), str(b_ff)],
            "two_epsilon_checks": families}


def compatible(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(compatible(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(compatible(x, y) for x, y in zip(a, b))
    if isinstance(a, (float, int)) and isinstance(b, (float, int)):
        return math.isclose(a, b, rel_tol=2e-11, abs_tol=2e-12)
    return a == b


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    file = ROOT / "independent_results.json"
    calculation = result()
    if args.write:
        with file.open("x", encoding="utf-8") as output:
            output.write(json.dumps(calculation, indent=2, allow_nan=False) + "\n")
    else:
        assert compatible(calculation, json.loads(file.read_text(encoding="utf-8")))
    print("PASS independent cofactor, threshold and mixed-current checks")
