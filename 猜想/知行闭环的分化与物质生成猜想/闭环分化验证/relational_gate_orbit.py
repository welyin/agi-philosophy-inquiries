"""Certified phase coverage and exact reversal at fixed Schmidt weights 3/4, 1/4."""

import argparse
import json
from fractions import Fraction
from functools import cmp_to_key
from math import isqrt
from pathlib import Path


ONE = (Fraction(1), Fraction(0))
GATE = (Fraction(3, 5), Fraction(4, 5))
QUARTER_TURN = (Fraction(0), Fraction(1))
PRECISION = 10 ** 12


def multiply(left, right):
    x, y = left
    u, v = right
    return (x * u - y * v, x * v + y * u)


def inverse(point):
    return point[0], -point[1]


def finite_words(generator, budget):
    points = {ONE}
    forward = ONE
    for _ in range(budget):
        forward = multiply(forward, generator)
        assert multiply(forward, inverse(forward)) == ONE
        points.update((forward, inverse(forward)))
    return points


def angular_order(p, q):
    # First sort by semicircle, then by the oriented area within it.
    def lower(point):
        x, y = point
        return y < 0 or (y == 0 and x < 0)

    if lower(p) != lower(q):
        return 1 if lower(p) else -1
    cross = p[0] * q[1] - p[1] * q[0]
    return int(cross < 0) - int(cross > 0)


def sqrt_bounds(value):
    """Rational enclosure with denominator PRECISION, by integer squaring."""
    scaled = value * PRECISION ** 2
    low = isqrt(scaled.numerator // scaled.denominator)
    high = low if low ** 2 == scaled else low + 1
    return Fraction(low, PRECISION), Fraction(high, PRECISION)


def covering_interval(points):
    # The nearest reachable state is at most half the largest empty arc away.
    ordered = sorted(points, key=cmp_to_key(angular_order))
    cosines = []
    for p, q in zip(ordered, ordered[1:] + ordered[:1]):
        assert p[0] * q[1] - p[1] * q[0] > 0  # each gap is < pi
        cosines.append(p[0] * q[0] + p[1] * q[1])
    gap_cosine = min(cosines)
    inner_lo, inner_hi = sqrt_bounds((1 + gap_cosine) / 2)
    # With visibility sqrt(3)/2, D_max² = (3/8)(1-cos(gap/2)).
    low, _ = sqrt_bounds(Fraction(3, 8) * (1 - inner_hi))
    _, high = sqrt_bounds(Fraction(3, 8) * (1 - inner_lo))
    return [int(low * PRECISION), int(high * PRECISION)]


def read_and_reverse(points, generator):
    # On each |i><j| of the operated qubit, arbitrary partner operators
    # are unaffected. Test the four exact scalar multipliers separately.
    for unit in points:
        entries = (ONE, unit)
        for left in entries:
            for right in entries:
                forward = multiply(left, inverse(right))
                backward = multiply(inverse(left), right)
                assert multiply(backward, forward) == ONE
        # These are ensemble XX expectations divided by the known visibility.
        direct_xx = unit[0]
        shifted_xx = multiply(unit, generator)[0]
        recovered_sine = (generator[0] * direct_xx - shifted_xx) / generator[1]
        assert (direct_xx, recovered_sine) == unit


def run():
    checks = {}
    for name, generator in (("rational_phase", GATE), ("quarter_turn", QUARTER_TURN)):
        samples = {}
        previous_upper = PRECISION
        for budget in (4, 16, 64, 256):
            points = finite_words(generator, budget)
            assert len(points) == (2 * budget + 1 if name == "rational_phase" else 4)
            read_and_reverse(points, generator)
            interval = covering_interval(points)
            assert interval[1] <= previous_upper
            previous_upper = interval[1]
            samples[str(budget)] = {"distinct_words_as_phases": len(points),
                                     "covering_radius_bounds_times_1e12": interval}
        checks[name] = samples
    return {"schmidt_weights": ["3/4", "1/4"], "results": checks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    output = Path(__file__).with_name("relational_gate_orbit_results.json")
    if args.check:
        assert json.loads(output.read_text(encoding="utf-8")) == result
    else:
        with output.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
