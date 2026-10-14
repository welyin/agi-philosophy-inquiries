"""Independent exact checks of an unnumbered conserved-quadrupole screen.

Uses differential-polynomial identities and Cartesian Taylor jets, not the
author's implementation. Default mode compares a saved review record and never
writes. --write creates that record once, after all three author files exist.
No physical switching, source preparation, detector or GR claim is certified.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
AUTHOR_NAMES = ("proof.md", "check.py", "results.json")
ZERO = (0, 0, 0)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plus(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, F(0)) + v
        if not out[k]:
            del out[k]
    return out


def scale(a, c):
    return {k: c*v for k, v in a.items() if c*v}


def times(a, b, degree=4):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x+y for x, y in zip(ka, kb))
            if sum(k) <= degree:
                out[k] = out.get(k, F(0)) + va*vb
    return {k: v for k, v in out.items() if v}


def radius_power_jet(R, exponent):
    """Taylor polynomial of (R^2+2Rz+x^2+y^2+z^2)^exponent."""
    u = {(0, 0, 1): 2/R, (2, 0, 0): 1/R**2,
         (0, 2, 0): 1/R**2, (0, 0, 2): 1/R**2}
    out, power, coefficient = {ZERO: F(1)}, {ZERO: F(1)}, F(1)
    for n in range(1, 5):
        power = times(power, u)
        coefficient *= (exponent-n+1)/n
        out = plus(out, scale(power, coefficient))
    assert exponent in (F(1, 2), F(-1, 2))
    return scale(out, R if exponent > 0 else 1/R)


def D_squared(jet):
    return (24*jet.get((4, 0, 0), 0)
            + 24*jet.get((0, 4, 0), 0)
            - 8*jet.get((2, 2, 0), 0))


def xy_hessian(jet):
    return 2*(jet.get((2, 0, 0), 0)+jet.get((0, 2, 0), 0))


def integral_product(a, b):
    return sum((x*y/F(i+j+1) for i, x in enumerate(a)
                for j, y in enumerate(b)), F(0))


def evaluate():
    # Differential symbols are (dt, dx, dy, dz). T is symmetric; all entries
    # operate on the same scalar test distribution, so derivatives commute.
    dt, dx, dy, dz = [{tuple(int(j == i) for j in range(4)): F(1)}
                      for i in range(4)]
    deriv = (dt, dx, dy, dz)
    T = [[{} for _ in range(4)] for _ in range(4)]
    T[0][0] = plus(times(dx, dx), scale(times(dy, dy), -1))
    T[0][1] = T[1][0] = scale(times(dt, dx), -1)
    T[0][2] = T[2][0] = times(dt, dy)
    T[1][1] = times(dt, dt)
    T[2][2] = scale(times(dt, dt), -1)
    div = []
    for nu in range(4):
        value = {}
        for mu in range(4):
            value = plus(value, times(deriv[mu], T[mu][nu]))
        div.append(value)
    assert all(not v for v in div)
    assert not plus(T[1][1], plus(T[2][2], T[3][3]))

    # Author-independent Cartesian derivative reconstruction. All kernels in
    # this calculation have their common 1/pi removed, exactly.
    R, radius, duration = F(3), F(1, 4), F(1)
    radial_integral = lambda power: sum(
        (F((-1)**j*comb(4, j), power+2*j+1) for j in range(5)), F(0))
    moment_one = radius**2*radial_integral(4)/radial_integral(2)
    assert moment_one == 3*radius**2/13
    m2 = 2*moment_one
    inverse = radius_power_jet(R, F(-1, 2))
    length = radius_power_jet(R, F(1, 2))
    G1 = scale(inverse, F(1, 4))
    G2 = plus(scale(length, F(-1, 8)), scale(inverse, -m2/24))
    coeff_ab = D_squared(G1)/2
    coeff_deriv = 2*xy_hessian(G1)+D_squared(G2)/2
    assert coeff_ab == F(3, 2)/R**5
    assert coeff_deriv == -F(3, 4)/R**3-m2/(4*R**5)
    assert R > 2*radius and 0 < duration < R-2*radius
    # In general C_ab/(-C_a'b')=6/(3R^2+m2)<=2/R^2.
    ratio = coeff_ab/(-coeff_deriv)
    assert ratio == 6/(3*R**2+m2) <= 2/R**2
    assert 9/duration**2 > ratio  # pi>3 gives strict Poincare separation.

    # Concrete compact C^2 profile p=u^3(1-u)^3. Not mislabeled C-infinity.
    p = [F(0), F(0), F(0), F(1), F(-3), F(3), F(-1)]
    dp = [i*p[i] for i in range(1, len(p))]
    ip2, idp2 = integral_product(p, p), integral_product(dp, dp)
    assert ip2 == F(1, 12012) and idp2 == F(1, 770)
    pairing_pi = coeff_ab*duration*ip2 + coeff_deriv*idp2/duration
    assert pairing_pi < 0
    # Since 3<pi<22/7 and pairing_pi<0, these are correctly ordered bounds.
    pairing_lo, pairing_hi = pairing_pi/3, pairing_pi/F(22, 7)
    assert pairing_lo < pairing_hi < 0
    return {
        "status": "UNNUMBERED_MATHEMATICAL_SCREEN_REVIEW",
        "groups": 3,
        "blockers": 0,
        "group1_conservation": {
            "method": "commuting four-variable differential polynomial",
            "divergence_polynomials_identically_zero": len(div),
            "spatial_trace_identically_zero": True,
            "does_not_establish_positive_complete_matter_source": True,
        },
        "group2_full_constraint_kernel": {
            "method": "Cartesian Taylor jets through degree four, exact fractions",
            "R": str(R), "radius_each": str(radius), "duration": str(duration),
            "second_moment_each": str(moment_one), "m2": str(m2),
            "pi_times_ab_coefficient": str(coeff_ab),
            "pi_times_derivative_coefficient": str(coeff_deriv),
            "coefficient_ratio": str(ratio),
            "full_momentum_terms_included": True,
            "poincare_sign_is_strict": True,
        },
        "group3_finite_profile": {
            "regularity": "compact C2 profile; distributional conservation exact",
            "integral_p_squared": str(ip2), "integral_pprime_squared": str(idp2),
            "pi_times_pairing": str(pairing_pi),
            "pairing_lower_exact": str(pairing_lo),
            "pairing_upper_exact": str(pairing_hi),
            "pairing_interval_float": [float(pairing_lo), float(pairing_hi)],
        },
        "scope": [
            "No actual Q(t), controller stress or recoil realization is proved.",
            "No signal probability, quantum apparatus or physical xi selection is proved.",
            "This is not a numbered research round or an Einstein derivation.",
            "Whole conservative source differences need not separately have positive energy.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = evaluate()
    out["author_sha256"] = {name: sha(HERE/name) for name in AUTHOR_NAMES}
    out["reviewer_code_sha256"] = sha(Path(__file__))
    target = HERE/"review_checks.json"
    if args.write:
        with target.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(out, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print("WRITE_ONCE_PASS; unnumbered; groups=3; blockers=0")
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        if saved != out:
            raise AssertionError("Review result or bound source has changed")
        print("READ_ONLY_MATCH; unnumbered; groups=3; blockers=0")


if __name__ == "__main__":
    main()
