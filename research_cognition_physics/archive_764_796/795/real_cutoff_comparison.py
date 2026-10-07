"""Round 795: exact diagnostics for real repair and causal cutoff comparison.

These finite algebras calibrate the identities in the report, not continuum
anomaly coefficients or a positive state of the original interacting theory.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json
import sys

sys.dont_write_bytecode = True
import numpy as np
import cauchy_state_probe

HERE = Path(__file__).resolve().parent
ORDER = 4


@dataclass(frozen=True)
class C:
    re: Q = Q(0)
    im: Q = Q(0)

    def __post_init__(self):
        object.__setattr__(self, "re", Q(self.re))
        object.__setattr__(self, "im", Q(self.im))

    def __add__(self, other):
        other = cast(other)
        return C(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return C(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-cast(other))

    def __rsub__(self, other):
        return cast(other) + (-self)

    def __mul__(self, other):
        other = cast(other)
        return C(self.re * other.re - self.im * other.im,
                 self.re * other.im + self.im * other.re)

    __rmul__ = __mul__

    def conjugate(self):
        return C(self.re, -self.im)

    def __bool__(self):
        return bool(self.re or self.im)

    def __str__(self):
        return str(self.re) if not self.im else f"({self.re})+({self.im})i"


def cast(value):
    return value if isinstance(value, C) else C(value)


def matrix(rows):
    return np.array([[cast(x) for x in row] for row in rows], dtype=object)


ZERO = matrix([[0, 0], [0, 0]])
IDENTITY = matrix([[1, 0], [0, 1]])


def scale(a, c):
    return np.array([[cast(c) * x for x in row] for row in a], dtype=object)


def dagger(a):
    return np.array([[a[j, i].conjugate() for j in range(2)]
                     for i in range(2)], dtype=object)


def equal(a, b):
    return all(not x for x in (a - b).flat)


def clean(p):
    return {k: v for k, v in p.items() if not equal(v, ZERO)}


def add(p, q):
    out = {k: v.copy() for k, v in p.items()}
    for k, v in q.items():
        out[k] = out.get(k, ZERO) + v
    return clean(out)


def mul(p, q):
    out = {}
    for a, x in p.items():
        for b, y in q.items():
            key = (a[0] + b[0], a[1] + b[1])
            if sum(key) <= ORDER:
                out[key] = out.get(key, ZERO) + x @ y
    return clean(out)


def adj(p):
    return {k: dagger(v) for k, v in p.items()}


def same(p, q):
    return all(equal(p.get(k, ZERO), q.get(k, ZERO)) for k in p.keys() | q.keys())


def const(a):
    return {(0, 0): a}


ONE = const(IDENTITY)


def exp(p):
    assert (0, 0) not in p
    out, power = ONE, ONE
    for n in range(1, ORDER + 1):
        power = mul(power, p)
        out = add(out, {k: scale(v, Q(1, factorial(n))) for k, v in power.items()})
    return out


def real_part(p):
    return clean({k: np.array([[C(x.re) for x in row] for row in a], dtype=object)
                  for k, a in p.items()})


def repair_probe():
    """Ordinary conjugation of a real MC complex; not a Hilbert adjoint."""
    d = matrix([[0, 1], [0, 0]])
    e = matrix([[0, 0], [1, 0]])
    j = e @ d - d @ e
    rows = []
    for g in (Q(0), Q(1, 4), Q(1, 2), Q(1)):
        r = g * (1 - g)
        a = C(2, 3) * r
        transform = add(ONE, {(1, 0): scale(e, a)})
        inverse = add(ONE, {(1, 0): scale(e, -a)})
        assert same(mul(transform, inverse), ONE)
        complex_solution = mul(mul(transform, const(d)), inverse)
        assert not mul(complex_solution, complex_solution)
        averaged = real_part(complex_solution)
        defect = mul(averaged, averaged)
        expected = clean({(2, 0): scale(IDENTITY, 9 * r * r)})
        assert same(defect, expected)
        # Recompute the next equation after selecting the real first coefficient.
        first = scale(j, 2 * r)
        second = scale(e, -4 * r * r)
        assert equal(d @ second + second @ d, -(first @ first))
        repaired = clean({(0, 0): d, (1, 0): first, (2, 0): second})
        assert not mul(repaired, repaired)
        assert same(repaired, real_part(repaired))
        if r == 0:
            assert same(repaired, const(d))
        else:
            assert defect
        rows.append(dict(g=str(g), collar_marker=str(r),
                         whole_average_square_coefficient=str(9 * r * r),
                         sequential_real_repair_square_zero=True))
    return dict(samples=rows, exact_complex_MC_and_real_recursive_solutions=True,
                whole_solution_averaging_rejected=True,
                support_marker_only_not_spacetime_support_proof=True,
                physical_BV_conjugation_not_simulated=True)


def expectation(a, psi):
    return sum((psi[i].conjugate() * a[i, j] * psi[j]
                for i in range(2) for j in range(2)), C())


def cutoff_probe():
    x = matrix([[0, 1], [1, 0]])
    z = matrix([[1, 0], [0, -1]])
    contact = matrix([[2, 1], [1, -1]])
    # Three time ordered source pulses, including a genuine mixed source term.
    source = mul(mul(exp({(1, 0): scale(x, C(0, Q(1, 2)))}),
                     exp({(0, 1): scale(z, C(0, 1)),
                          (1, 1): scale(contact, C(0, Q(2, 7)))})),
                 exp({(1, 0): scale(x, C(0, Q(1, 2)))}))
    assert same(mul(adj(source), source), ONE)
    old_past = matrix([[Q(3, 5), -Q(4, 5)], [Q(4, 5), Q(3, 5)]])
    new_past = matrix([[Q(5, 13), -Q(12, 13)], [Q(12, 13), Q(5, 13)]])
    old_future = matrix([[0, -1], [1, 0]])
    new_future = matrix([[Q(8, 17), -Q(15, 17)], [Q(15, 17), Q(8, 17)]])

    def relative(future, past):
        s0 = future @ past
        sz = mul(mul(const(future), source), const(past))
        return mul(const(dagger(s0)), sz)

    old = relative(old_future, old_past)
    new = relative(new_future, new_past)
    assert same(old, relative(new_future, old_past))
    assert same(new, relative(old_future, new_past))
    u = dagger(old_past) @ new_past
    assert equal(dagger(u) @ u, IDENTITY)

    def alpha(p):
        return mul(mul(const(dagger(u)), p), const(u))

    assert same(alpha(old), new)
    assert same(alpha(adj(old)), adj(new))
    psi = (C(1), C(0))
    transported = tuple(sum((dagger(u)[i, j] * psi[j] for j in range(2)), C())
                        for i in range(2))
    assert sum((v.conjugate() * v for v in transported), C()) == C(1)
    checked = 0
    for key, a in old.items():
        assert expectation(a, psi) == expectation(new.get(key, ZERO), transported)
        checked += 1
    a = scale(old[(1, 0)], C(0, -1))
    b = scale(old[(0, 1)], C(0, -1))
    assert equal(a, dagger(a)) and equal(b, dagger(b))
    assert not equal(a @ b, b @ a)
    aa, bb = dagger(u) @ a @ u, dagger(u) @ b @ u
    for p, q in ((a @ b, aa @ bb), (b @ a, bb @ aa),
                 (a @ b @ a, aa @ bb @ aa)):
        assert equal(dagger(u) @ p @ u, q)
        assert expectation(p, psi) == expectation(q, transported)
    # Same algebra and the same numeric state vector need not be the same state.
    old_mean = expectation(a, psi)
    unmatched = expectation(aa, psi)
    assert old_mean != unmatched
    f = a + scale(b, C(0, Q(2, 3)))
    square_value = expectation(dagger(f) @ f, psi)
    assert not square_value.im and square_value.re >= 0
    transported_f = dagger(u) @ f @ u
    assert expectation(dagger(transported_f) @ transported_f, transported) == square_value
    without_contact = mul(mul(exp({(1, 0): scale(x, C(0, Q(1, 2)))}),
                             exp({(0, 1): scale(z, C(0, 1))})),
                         exp({(1, 0): scale(x, C(0, Q(1, 2)))}))
    assert not equal(source.get((1, 1), ZERO), without_contact.get((1, 1), ZERO))
    return dict(total_source_degree=ORDER, source_coefficients_checked=checked,
                relative_unitarity=True, future_switch_cancels=True,
                one_source_independent_unitary_compares_both_sources=True,
                mixed_source_contact_nonzero=True, adjoints_preserved=True,
                noncommuting_ordered_products_checked=3,
                original_vector_untransported_mean=str(unmatched),
                original_mean=str(old_mean),
                jointly_transported_mean=str(expectation(aa, transported)),
                positive_square_value=str(square_value),
                finite_state_transport_identity=True,
                original_continuum_state_match_proven_by_this_probe=False)


def run():
    cauchy = cauchy_state_probe.run()
    assert cauchy == json.loads((HERE/'cauchy_state_probe_results.json').read_text(encoding='utf-8'))
    return dict(round=795, fresh_test_groups=3, real_repair=repair_probe(),
                causal_cutoff=cutoff_probe(), first_order_cauchy=cauchy,
                all_checks_passed=True, exact_arithmetic="Gaussian rationals / Fraction",
                scope="Finite calibration of recursive reality, common unitary cutoff comparison, and initial-data matching. Continuum conclusions rely on the report's assumptions and proof.",
                original_interacting_positive_state_constructed=False)


if __name__ == '__main__':
    result = run()
    (HERE/'real_cutoff_comparison_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
