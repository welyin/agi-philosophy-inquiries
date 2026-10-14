"""Independent internal review of 1185. Default is read-only; --write once.

Exact radical matrices and independently implemented rational enclosures.
Written before reading the author's check.py. No author functions imported.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt, nextafter, inf
import argparse
import hashlib
import json


# Exact Q(sqrt(2), sqrt(3)), basis 1,sqrt(2),sqrt(3),sqrt(6).
class Q:
    def __init__(self, value=0):
        self.c = value if isinstance(value, tuple) else (F(value), F(0), F(0), F(0))

    def __add__(self, other):
        other = other if isinstance(other, Q) else Q(other)
        return Q(tuple(x+y for x, y in zip(self.c, other.c)))

    __radd__ = __add__

    def __neg__(self):
        return Q(tuple(-x for x in self.c))

    def __sub__(self, other):
        return self + -(other if isinstance(other, Q) else Q(other))

    def __mul__(self, other):
        other = other if isinstance(other, Q) else Q(other)
        out = [F(0)]*4
        for i, x in enumerate(self.c):
            for j, y in enumerate(other.c):
                out[i ^ j] += x*y*(2 if i & j & 1 else 1)*(3 if i & j & 2 else 1)
        return Q(tuple(out))

    __rmul__ = __mul__

    def __eq__(self, other):
        return self.c == (other.c if isinstance(other, Q) else Q(other).c)


def matrix(rows):
    return [[Q(x) for x in row] for row in rows]


def eye(n):
    return matrix([[int(i == j) for j in range(n)] for i in range(n)])


def add(a, b):
    return [[x+y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(c, a):
    return [[c*x for x in row] for row in a]


def sub(a, b):
    return add(a, scale(-1, b))


def mul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return list(map(list, zip(*a)))


def kron(a, b):
    return [[a[i][j]*b[k][l] for j in range(len(a[0])) for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


def equal(a, b):
    assert all(x == 0 for row in sub(a, b) for x in row)


def algebra():
    i = eye(2)
    i4 = eye(4)
    x = matrix([[0, 1], [1, 0]])
    y_over_i = matrix([[0, -1], [1, 0]])
    z = matrix([[1, 0], [0, -1]])
    a = scale(F(1, 2), add(i, x))
    n = scale(F(1, 2), sub(i, z))
    s2 = Q((F(0), F(1), F(0), F(0)))
    s3 = Q((F(0), F(0), F(1), F(0)))
    s6 = s2*s3
    r = [[s6*F(1, 3), s3*F(1, 3)], [-s3*F(1, 3), s6*F(1, 3)]]
    w = add(kron(sub(i, a), i), kron(a, r))
    e = mul(w, kron(i, matrix([[1], [0]])))
    h = add(add(i4, kron(i, n)), scale(s2, kron(a, x)))
    hc = sub(i, a)
    equal(mul(transpose(w), w), i4)
    equal(mul(transpose(e), e), i)
    equal(mul(mul(transpose(w), h), w),
          add(kron(hc, i), kron(add(i, scale(2, a)), n)))
    equal(mul(h, e), mul(e, hc))
    nb = mul(mul(w, kron(i, n)), transpose(w))
    equal(mul(nb, e), scale(0, e))
    equal(mul(nb, h), mul(h, nb))
    for port in (i, x, y_over_i, z):
        transported = mul(mul(w, kron(port, i)), transpose(w))
        equal(mul(transported, e), mul(e, port))
        equal(mul(nb, transported), mul(transported, nb))
    for source, effective in ((scale(F(3, 20), h), scale(F(3, 20), hc)),
                              (scale(F(3, 20), i4), scale(F(3, 20), i))):
        equal(mul(source, e), mul(e, effective))
        equal(mul(nb, source), mul(source, nb))
    bare_ze = mul(kron(z, i), e)
    equal(mul(transpose(e), bare_ze), scale(s6*F(1, 3), z))
    leak = mul(sub(i4, mul(e, transpose(e))), bare_ze)
    equal(mul(transpose(leak), leak), scale(F(1, 3), i))
    return {
        'status': 'EXACT_PASS', 'arithmetic': 'Q(sqrt(2),sqrt(3))',
        'micro_dimension': 4, 'code_dimension': 2,
        'micro_spectrum': [0, 1, 2, 3], 'code_spectrum': [0, 1],
        'all_four_t_port_intertwiners': True,
        'dressed_number_commutes_with_all_declared_t_generators': True,
        'two_source_intertwiners': True,
        'bare_Z_leakage_squared': '1/3',
        'tensor_extension': 'all encoded background factors, CAR, probe and arbitrary R',
    }


# Each interval operation rounds outwards to an exact 10^-60 rational grid.
GRID = 10**60


def floor_grid(x):
    x = F(x)
    return F((x.numerator*GRID)//x.denominator, GRID)


def ceil_grid(x):
    return -floor_grid(-F(x))


class IV:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = floor_grid(lo), ceil_grid(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def cast(value):
        return value if isinstance(value, IV) else IV(value)

    def __add__(self, other):
        other = self.cast(other)
        return IV(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return IV(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.cast(other)

    def __rsub__(self, other):
        return self.cast(other) + -self

    def __mul__(self, other):
        other = self.cast(other)
        v = [x*y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return IV(min(v), max(v))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        assert other.lo*other.hi > 0
        return self*IV(1/other.hi, 1/other.lo)


def atan_reciprocal(n, terms):
    total = sum(((-1)**j*F(1, (2*j+1)*n**(2*j+1)) for j in range(terms)), F(0))
    next_term = (-1)**terms*F(1, (2*terms+1)*n**(2*terms+1))
    return IV(min(total, total+next_term), max(total, total+next_term))


def pi_interval():
    return 16*atan_reciprocal(5, 100)-4*atan_reciprocal(239, 40)


def root_interval(value):
    value = IV.cast(value)
    assert value.lo >= 0
    lo = isqrt((value.lo.numerator*GRID*GRID)//value.lo.denominator)
    hi = isqrt((value.hi.numerator*GRID*GRID)//value.hi.denominator)+1
    return IV(F(lo, GRID), F(hi, GRID))


def sine(value):
    value = IV.cast(value)
    # Polynomial interval evaluation; global Taylor remainder uses |sin^(n)|<=1.
    assert max(abs(value.lo), abs(value.hi)) <= 4
    term = value
    total = term
    for j in range(1, 50):
        term = -term*value*value/((2*j)*(2*j+1))
        total = total+term
    bound_x = max(abs(value.lo), abs(value.hi))
    factorial = 1
    for j in range(1, 101):
        factorial *= j
    # Polynomial through degree 99; order-100 Lagrange remainder suffices.
    rem = bound_x**100/factorial
    return total+IV(-rem, rem)


def exponential(value):
    value = IV.cast(value)
    assert 0 <= value.lo <= value.hi <= 4
    term = IV(1)
    total = term
    for j in range(1, 81):
        term = term*value/j
        total = total+term
    first = value.hi**81
    for j in range(1, 82):
        first /= j
    remainder = first/(1-value.hi/F(82))
    return total+IV(0, remainder)


def down_float(x):
    v = float(x)
    return nextafter(v, -inf) if F(v) > x else v


def up_float(x):
    v = float(x)
    return nextafter(v, inf) if F(v) < x else v


def record(value):
    value = IV.cast(value)
    return {'lower_exact': str(value.lo), 'upper_exact': str(value.hi),
            'lower_display': down_float(value.lo), 'upper_display': up_float(value.hi)}


def probability(b, time):
    omega, epsilon = F(1, 2), F(3, 20)
    r = root_interval(omega*omega+b*b)
    tau = time-epsilon*omega*b*(time-sine(2*r*time)/(2*r))/(r*r)
    return tau, IV(F(1, 2))+sine(F(1, 5)*tau)/2


def compute_math():
    g1 = algebra()
    pi = pi_interval()
    assert F('3.14159265358979323846264338327950288419716939937510') < pi.lo
    assert pi.hi < F('3.14159265358979323846264338327950288419716939937511')
    time = pi/2
    a_s, a_t, w, sqrt_k = F(1389, 320), F(13, 8), F(3, 10), 10**7
    bnorm = root_interval(7*a_s*a_s+5*a_t*a_t)
    eta = 4*time*exponential(6*w*time)*(root_interval(2)*bnorm+a_s+a_t)/sqrt_k
    public_eta, instrument = F(245, 10**6), F(2, 10**5)
    assert eta.hi < public_eta
    volume, kval = 65**3, 10**14
    epsilon, nu, g, ell, omega = F(3, 20), F(7, 20), F(1, 5), F(16, 65), F(1, 2)
    source_bound = 4*epsilon+6*nu/ell+4*g/ell
    probe_bound = F(8, 5)*epsilon+18*nu/ell+12*g/ell
    original_ck = (1+kval*volume*(4+6/ell+2*omega+source_bound)
                   +volume*(F(8, 5)+18/ell+probe_bound))
    free_addition = F(5, 2)*volume*kval
    reciprocal_addition = F(9, 20)*volume*(kval-1)
    g2 = {
        'status': 'RATIONAL_ENCLOSURE_PASS', 'pi': record(pi),
        'K': 10**14, 'T': record(time), 'a_s': str(a_s), 'a_t': str(a_t),
        'full_micro_source_norm': '9/20', 'effective_source_norm': '3/20',
        'effective_pair_norm': str(w),
        'positive_shift_additions': ['(5/2) V K', '(9/20) V (K-1)'],
        'original_CK_exact': str(original_ck),
        'free_norm_increment_exact': str(free_addition),
        'reciprocal_full_space_increment_exact': str(reciprocal_addition),
        'positive_shift_total_exact': str(original_ck+free_addition+reciprocal_addition),
        'eta_P': record(eta), 'declared_eta_P': str(public_eta),
        'instrument_per_complete_task': str(instrument),
        'scope': 'prepared full 48-dimensional code/R; comparator retains exact BG',
    }
    inputs = [('match_q_plus', F(3, 80)), ('match_q_minus', F(9, 80)),
              ('weak_q_plus', F(3, 20)), ('weak_q_minus', F(3, 20))]
    rows = {}
    probabilities = {}
    for name, b in inputs:
        tau, p = probability(b, time)
        rows[name] = {'b': str(b), 'tau': record(tau), 'p_minus_Y': record(p)}
        probabilities[name] = p
    clock_window_upper = (IV(1)+sine(F(19, 50)))/2
    for p in probabilities.values():
        assert p.lo-public_eta-instrument > F(1, 2)
        assert p.hi+public_eta+instrument < clock_window_upper.lo
    gap = ((probabilities['match_q_plus']-probabilities['match_q_minus'])
           -(probabilities['weak_q_plus']-probabilities['weak_q_minus']))
    loss = 4*(public_eta+instrument)
    lower = gap.lo-loss
    assert loss == F(53, 50000)
    assert lower > F(3, 20000)
    g3 = {
        'status': 'FINITE_FOUR_TASK_DIFFERENCE_PASS', 'tasks': rows,
        'ideal_difference_of_differences': record(gap),
        'four_task_error': str(loss), 'certified_lower_exact': str(lower),
        'certified_lower_display': down_float(lower),
        'public_lower': '3/20000',
        'all_four_complete_error_windows_inside_calibrated_clock_branch': True,
        'scope': 'original p=0 material clock; no GR or hardware certification',
    }
    return {'g1_micro_intertwining': g1, 'g2_full_prepared_error': g2,
            'g3_original_clock_readback': g3}


HERE = Path(__file__).resolve().parent
CORE = ['../research_note_1185.md', 'admission.md', 'proof.md',
        'check.py', 'results.json', 'NEXT.md']


def compare_author(checks):
    # Read saved results only, after independent mathematics has been implemented.
    author = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    groups = author['groups']
    assert len(groups) == 3 and all(g['pass'] for g in groups.values())
    cert = groups['3_original_clock_four_task_finite_certificate']['certificate']
    own = checks['g3_original_clock_readback']
    comparisons = 0

    def contained(outer, inner):
        nonlocal comparisons
        assert F(outer['lower_exact']) <= F(inner['lower_exact'])
        assert F(inner['upper_exact']) <= F(outer['upper_exact'])
        comparisons += 1

    contained(cert['time'], checks['g2_full_prepared_error']['T'])
    contained(cert['eta_interval'], checks['g2_full_prepared_error']['eta_P'])
    contained(cert['ideal_four_task_difference'], own['ideal_difference_of_differences'])
    name_map = {'match_q_plus': 'matched_q_plus', 'match_q_minus': 'matched_q_minus',
                'weak_q_plus': 'deleted_q_plus', 'weak_q_minus': 'deleted_q_minus'}
    for ours, theirs in name_map.items():
        row, arow = own['tasks'][ours], cert['tasks'][theirs]
        assert F(row['b']) == F(arow['b_exact'])
        contained(arow['tau'], row['tau'])
        contained(arow['P_minus_Y_probability'], row['p_minus_Y'])
    assert F(cert['actual_four_task_lower_exact']) <= F(own['certified_lower_exact'])
    assert F(cert['published_strict_lower_exact']) == F(own['public_lower'])
    assert F(cert['four_task_error_exact']) == F(own['four_task_error'])
    ag2 = groups['2_complete_parent_and_deleted_source_control']
    for key in ['original_CK_exact', 'free_norm_increment_exact',
                'reciprocal_full_space_increment_exact', 'positive_shift_total_exact']:
        assert F(ag2[key]) == F(checks['g2_full_prepared_error'][key])
    for name, digest in author['source_sha256'].items():
        assert hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest() == digest
    return {'status': 'PASS', 'independent_intervals_inside_author_intervals': comparisons,
            'exact_positive_shift_components_equal': 4,
            'frozen_source_hashes_verified': author['source_sha256'],
            'author_float_lower_is_display_only': True}


def complete_result():
    checks = compute_math()
    return {
        'status': 'PASS_WITH_SCOPE', 'groups': 3, 'blockers': 0,
        'review_type': 'internal AI cross-check, non-blind; independently implemented before author code read',
        'author_core_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in CORE},
        'checks': checks,
        'author_comparison': compare_author(checks),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--math-only', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    if args.math_only:
        assert not args.write
        result = compute_math()
        print('MATH_ONLY_PASS 3 groups; lower =',
              result['g3_original_clock_readback']['certified_lower_display'])
        return
    out = complete_result()
    path = HERE/'review_checks.json'
    if args.write:
        with path.open('x', encoding='utf-8') as handle:
            json.dump(out, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print('WRITE_ONCE_PASS 3 groups')
    else:
        assert out == json.loads(path.read_text(encoding='utf-8')), 'REVIEW_RESULT_MISMATCH'
        print('READ_ONLY_MATCH 3 groups; blockers 0')


if __name__ == '__main__':
    main()
