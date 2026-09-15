"""Small outward-rounded dyadic intervals for independently checking certificates."""

import math
import unittest
from fractions import Fraction


BITS = 100
SCALE = 1 << BITS


def ceil_div(value, divisor):
    return -((-value) // divisor)


class Interval:
    def __init__(self, lower, upper=None):
        self.lo = int(lower)
        self.hi = int(lower if upper is None else upper)
        if self.lo > self.hi:
            raise ValueError("Invalid interval.")

    @classmethod
    def rational(cls, numerator, denominator=1):
        return cls(numerator * SCALE // denominator, ceil_div(numerator * SCALE, denominator))

    @classmethod
    def exact(cls, value):
        if isinstance(value, cls):
            return value
        if isinstance(value, Fraction):
            return cls.rational(value.numerator, value.denominator)
        if isinstance(value, int):
            return cls(value * SCALE)
        raise TypeError("Use integers or exact fractions, not floats.")

    def __add__(self, other):
        other = self.exact(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.exact(other))

    def __rsub__(self, other):
        return self.exact(other) - self

    def __mul__(self, other):
        other = self.exact(other)
        values = (self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi)
        return Interval(min(values) // SCALE, ceil_div(max(values), SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.exact(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("Interval contains zero.")
        values = [Fraction(x * SCALE, y) for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        lower, upper = min(values), max(values)
        return Interval(lower.numerator // lower.denominator, ceil_div(upper.numerator, upper.denominator))

    def __rtruediv__(self, other):
        return self.exact(other) / self

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Use nonnegative integer powers.")
        result, current = Interval.exact(1), self
        while exponent:
            if exponent & 1:
                result = result * current
            current = current * current
            exponent >>= 1
        return result

    def sqrt(self):
        if self.lo < 0:
            raise ValueError("Negative square-root input.")
        return Interval(math.isqrt(self.lo * SCALE), math.isqrt(self.hi * SCALE) + 1)

    def abs_upper(self):
        return Interval(max(abs(self.lo), abs(self.hi)))

    def widen(self, radius):
        bound = self.exact(radius).abs_upper().hi
        return Interval(self.lo - bound, self.hi + bound)

    def floats(self):
        return [self.lo / SCALE, self.hi / SCALE]


def atan_reciprocal(denominator, terms=40):
    x = Interval.rational(1, denominator)
    value = Interval.exact(0)
    for k in range(terms):
        value += ((-1)**k) * x**(2 * k + 1) / (2 * k + 1)
    # Alternating-series remainder, enlarged symmetrically for convenience.
    return value.widen(x**(2 * terms + 1) / (2 * terms + 1))


def pi_interval():
    # Machin's identity, verified by the tangent addition formula and 0<pi/4<pi/2.
    return 16 * atan_reciprocal(5) - 4 * atan_reciprocal(239)


def sin_interval(x):
    x = Interval.exact(x)
    if x.abs_upper().hi > 4 * SCALE:
        raise ValueError("Reduce the argument to [-4,4] first.")
    total = Interval.exact(0)
    for k in range(32):
        total += (-1)**k * x**(2 * k + 1) / math.factorial(2 * k + 1)
    return total.widen(x.abs_upper()**65 / math.factorial(65))


def cos_interval(x):
    x = Interval.exact(x)
    if x.abs_upper().hi > 4 * SCALE:
        raise ValueError("Reduce the argument to [-4,4] first.")
    total = Interval.exact(0)
    for k in range(33):
        total += (-1)**k * x**(2 * k) / math.factorial(2 * k)
    return total.widen(x.abs_upper()**66 / math.factorial(66))


class CertifiedIntervalTests(unittest.TestCase):
    def test_arithmetic_contains_exact_fraction_answers(self):
        a, b = Fraction(-2, 7), Fraction(5, 11)
        for interval, answer in ((Interval.exact(a) + b, a + b), (Interval.exact(a) * b, a * b),
                                 (Interval.exact(a) / b, a / b), (Interval.exact(a) - b, a - b)):
            self.assertLessEqual(Fraction(interval.lo, SCALE), answer)
            self.assertGreaterEqual(Fraction(interval.hi, SCALE), answer)

    def test_square_root_bounds_by_integer_comparisons(self):
        value = Interval.rational(2)
        answer = value.sqrt()
        self.assertLessEqual(answer.lo**2, 2 * SCALE**2)
        self.assertGreaterEqual(answer.hi**2, 2 * SCALE**2)

    def test_pi_and_small_angle_values_agree_with_independent_floats(self):
        lower, upper = pi_interval().floats()
        self.assertAlmostEqual((lower + upper) / 2, math.pi, places=14)
        for numerator in (-1, 0, 1, 7):
            angle = Interval.rational(numerator, 4)
            for result, expected in ((sin_interval(angle), math.sin(numerator / 4)),
                                     (cos_interval(angle), math.cos(numerator / 4))):
                lower, upper = result.floats()
                self.assertAlmostEqual((lower + upper) / 2, expected, places=14)

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(TypeError):
            Interval.exact(0.1)
        with self.assertRaises(ZeroDivisionError):
            Interval.rational(1) / Interval(-1, 1)
        with self.assertRaises(ValueError):
            sin_interval(Interval.exact(5))
