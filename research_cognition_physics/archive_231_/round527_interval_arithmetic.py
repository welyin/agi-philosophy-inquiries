"""Outward binary intervals for round 527's finite calibration algorithm.

Each cached factory owns its scale and transcendental enclosure.  The frozen
round-526 module and its global SCALE are neither imported nor modified.
"""
from fractions import Fraction
from functools import lru_cache
import math
import operator
from types import SimpleNamespace


class PrecisionNeeded(ArithmeticError):
    """An enclosure has not yet certified an operation's stated domain.

    Increasing precision is appropriate only when the caller has established
    the exact expression's domain independently.  Invalid exact input is not
    repaired by this exception or by arbitrarily increasing precision.
    """


def _ceildiv(numerator, denominator):
    return -((-numerator) // denominator)


@lru_cache(maxsize=None)
def arithmetic(bits=64):
    """Return an independent interval API at any integer precision >= 64.

    Interval endpoints are integers divided by Interval.scale.  Addition,
    subtraction, multiplication, and division by a certified positive interval
    round outward.  ``gradient(q)`` is the original 525 local three-field force
    gradient; ``leading(x, T)`` is c*x*F(x)*T**2/2.

    exp is used only on [-1, 1].  At 64 bits it uses the inherited 24th-order
    Taylor polynomial with remainder 3/25!.  Its order then grows with bits, so
    both arithmetic error and the analytic Taylor remainder vanish as precision
    increases; there is no fixed-error floor in adaptive finite graph sums.
    """
    bits = operator.index(bits)
    if bits < 64:
        raise ValueError('binary interval precision must be at least 64 bits')
    binary_scale = 1 << bits
    exp_order = 24 + (bits - 64 + 7) // 8

    class Interval:
        __slots__ = ('lo', 'hi')
        scale = binary_scale

        def __init__(self, lo, hi=None):
            self.lo = operator.index(lo)
            self.hi = self.lo if hi is None else operator.index(hi)
            if self.lo > self.hi:
                raise ValueError('interval endpoints are reversed')

        @staticmethod
        def of(value):
            if isinstance(value, Interval):
                return value
            value = Fraction(value) * binary_scale
            return Interval(value.numerator // value.denominator,
                            _ceildiv(value.numerator, value.denominator))

        def __add__(self, other):
            other = Interval.of(other)
            return Interval(self.lo + other.lo, self.hi + other.hi)

        __radd__ = __add__

        def __neg__(self):
            return Interval(-self.hi, -self.lo)

        def __sub__(self, other):
            return self + -Interval.of(other)

        def __rsub__(self, other):
            return Interval.of(other) + -self

        def __mul__(self, other):
            other = Interval.of(other)
            products = [a*b for a in (self.lo, self.hi)
                        for b in (other.lo, other.hi)]
            return Interval(min(products) // binary_scale,
                            _ceildiv(max(products), binary_scale))

        __rmul__ = __mul__

        def __truediv__(self, other):
            other = Interval.of(other)
            if other.lo <= 0:
                raise PrecisionNeeded('division requires a certified positive denominator')
            quotients = [Fraction(a * binary_scale, b)
                         for a in (self.lo, self.hi)
                         for b in (other.lo, other.hi)]
            low, high = min(quotients), max(quotients)
            return Interval(low.numerator // low.denominator,
                            _ceildiv(high.numerator, high.denominator))

        def __rtruediv__(self, other):
            return Interval.of(other) / self

        def __repr__(self):
            return f'Interval({self.lo}, {self.hi}; scale=2**{bits})'

    def exp_interval(value):
        value = Interval.of(value)
        if value.lo < -binary_scale or value.hi > binary_scale:
            raise PrecisionNeeded('exp enclosure must be contained in [-1, 1]')
        term = total = Interval.of(1)
        for k in range(1, exp_order + 1):
            term = term * value / k
            total = total + term
        # Lagrange remainder: e*|u|^(m+1)/(m+1)! <= e/(m+1)! < 3/(m+1)!.
        remainder = _ceildiv(3 * binary_scale, math.factorial(exp_order + 1))
        return Interval(total.lo - remainder, total.hi + remainder)

    root = math.isqrt(5 * binary_scale * binary_scale)
    sqrt5 = Interval(root, root + 1)
    vector = (-Interval.of(1) / sqrt5,
              Interval.of(2) / sqrt5,
              Interval.of(-1))

    def gradient(q):
        if len(q) != 3:
            raise ValueError('the local reference gradient has three components')
        q = [Interval.of(value) for value in q]
        gate = 1 / (1 + exp_interval(-2 * (q[0] - 1)))
        derivative = 2 * gate * (1 - gate)
        combination = sum(v * value for v, value in zip(vector, q))
        out = [gate * combination * v for v in vector]
        out[0] = out[0] + derivative * combination * combination / 2
        out[1] = out[1] + Fraction(5, 4) * q[1]
        return out

    def leading(x, time):
        q = [Interval.of(x), Interval.of(0), Interval.of(0)]
        return -gradient(q)[2] * time * time / 2

    return SimpleNamespace(bits=bits, scale=binary_scale, Interval=Interval,
                           gradient=gradient, leading=leading,
                           exp_interval=exp_interval, exp_order=exp_order,
                           sqrt5=sqrt5, vector=vector,
                           PrecisionNeeded=PrecisionNeeded)
