"""Independent rational/Simpson review of 1184; default never writes.

The author code was read before this implementation. No author functions are
imported. Four groups audit the same four scientific groups, not new work.
"""
from fractions import Fraction as F
from math import isqrt, factorial, comb
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
GRID = 10**44


def down(x):
    x = F(x)
    return F((x.numerator*GRID)//x.denominator, GRID)


def up(x):
    return -down(-F(x))


class R:
    def __init__(self, lo=0, hi=None):
        if isinstance(lo, R) and hi is None:
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo = down(lo)
            self.hi = up(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, b):
        b = R(b)
        return R(self.lo+b.lo, self.hi+b.hi)
    __radd__ = __add__

    def __neg__(self):
        return R(-self.hi, -self.lo)

    def __sub__(self, b):
        return self+-R(b)

    def __rsub__(self, b):
        return R(b)+-self

    def __mul__(self, b):
        b = R(b)
        ends = [x*y for x in (self.lo, self.hi) for y in (b.lo, b.hi)]
        return R(min(ends), max(ends))
    __rmul__ = __mul__

    def __truediv__(self, b):
        b = R(b)
        assert b.lo*b.hi > 0
        return self*R(1/b.hi, 1/b.lo)

    def __rtruediv__(self, b):
        return R(b)/self

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        result = R(1)
        for _ in range(n):
            result = result*self
        return result

    def sqrt(self):
        assert self.lo >= 0
        a = isqrt(self.lo.numerator*GRID**2//self.lo.denominator)
        b = isqrt(self.hi.numerator*GRID**2//self.hi.denominator)
        if b*b*self.hi.denominator < self.hi.numerator*GRID**2:
            b += 1
        return R(F(a, GRID), F(b, GRID))

    def pad(self, radius):
        radius = R(radius).hi
        return R(self.lo-radius, self.hi+radius)

    def trig(self, cosine=False):
        # Evaluate at the rational midpoint. Global Lipschitz constant 1
        # encloses the entire input interval, independent of monotonicity.
        center = (self.lo+self.hi)/2
        radius = (self.hi-self.lo)/2
        assert abs(center) <= 3
        x = R(center)
        term = R(1) if cosine else x
        total = term
        for j in range(1, 36):
            denominator = (2*j-1)*(2*j) if cosine else (2*j)*(2*j+1)
            term = -term*x*x/denominator
            total = total+term
        degree = 72 if cosine else 73
        return total.pad(radius+F(3**degree, factorial(degree)))

    def sin(self):
        return self.trig(False)

    def cos(self):
        return self.trig(True)

    def exp(self):
        assert 0 <= self.lo and self.hi <= 2
        term, total = R(1), R(1)
        for j in range(1, 81):
            term = term*self/j
            total = total+term
        # Lagrange remainder, using e^2 < 8.
        return total.pad(F(8*2**81, factorial(81)))

    def data(self):
        return {'lower': str(self.lo), 'upper': str(self.hi),
                'display_only': [float(self.lo), float(self.hi)]}


def atan_inverse(n):
    value = sum((F((-1)**j, (2*j+1)*n**(2*j+1))
                 for j in range(100)), F(0))
    return R(value).pad(F(1, 201*n**201))


# atan(1/2)+atan(1/3)=pi/4; distinct from the author's Machin formula.
PI = 4*(atan_inverse(2)+atan_inverse(3))
ELL = F(16, 65)
TIME = PI/2
CARRIER = PI/8
KSTAR = (CARRIER*ELL).sin()/ELL
W = R(1609).sqrt()/40
NMAX = F(1)+F(36, 1609)
M1, M2, M3 = F(3, 250), F(9, 800), F(3, 250)
D2 = NMAX**2+M1
D4 = NMAX**4+6*NMAX**2*M1+3*M1**2+4*NMAX*M2+M3


def proper(sign, t):
    # sign labels q=+/-1/2. All constants reduced independently.
    return t+sign*F(18, 1609)*(t-(W*t).sin()/W)


def transport(sign, t, panels=128):
    assert panels % 2 == 0
    values = [proper(sign, t*F(j, panels)).cos() for j in range(panels+1)]
    total = values[0]+values[-1]
    for j in range(1, panels):
        total = total+(4 if j % 2 else 2)*values[j]
    integral = total*t/(3*panels)
    error = D4*t.hi**5/(180*panels**4)
    integral = integral.pad(error)
    tau = proper(sign, t)
    isotropic = t+F(7, 3)*(tau-t)
    anisotropic = R(3).sqrt()*integral/10
    return tau, [isotropic+anisotropic, isotropic-anisotropic, isotropic], error


def mathematics():
    assert PI.lo > F('3.1415926535897932384626433832794')
    assert PI.hi < F('3.1415926535897932384626433832796')
    assert R(F(18, 1609)*W).hi < M1
    assert R(F(9, 800)*W).hi < M3
    assert D2 < F(3, 2) and D4 < F(3, 2)

    # Exact finite Fourier identities, independent of gigantic spatial arrays.
    coeff = {m: F(comb(96, 48+m), 2**96) for m in range(-48, 49)}
    norm = sum((c*c for c in coeff.values()), F(0))
    second = sum((m*m*c*c for m, c in coeff.items()), F(0))/norm
    assert norm == F(comb(192, 96), 4**96) and norm > F(1, 20)
    assert second == F(2304, 191)
    lx, nx, carrier = 2097152, 8519680, 131072
    assert F(lx, nx) == ELL and carrier > 96 and carrier+48 < nx//2
    volume, copies = nx*65**2, 10**16
    b_source = 2*F(3, 20)*2+6*F(7, 20)/ELL+4*F(1, 5)/ELL
    b_probe = 2*F(3, 20)*F(4, 5)+18*F(7, 20)/ELL+12*F(1, 5)/ELL
    c_old = 1+copies*volume*(4+6/ELL+1+b_source)+volume*(F(8, 5)+18/ELL+b_probe)
    baseline = c_old+volume*(copies-1)*F(3, 40)
    assert baseline == 1+volume*(F(6693, 160)*copies+F(88187, 800))
    assert baseline == 15057429504000003967935262721
    # The second platform is displaced by an integer number of carrier periods.
    assert carrier % 2 == 0
    c_squared = (3*PI/16).cos()**2
    assert c_squared.hi < F(173, 250)
    tail = R(20*F(173, 250)**96).sqrt()
    assert tail.hi < F(9454, 10**11)
    rho = (2*PI/lx)*R(second).sqrt()

    a_s, a_t = F(1389, 320), F(13, 8)
    bb = 7*a_s*a_s+5*a_t*a_t
    source = 4*TIME*(F(9, 10)*TIME).exp()/10**8 * (R(2*bb).sqrt()+a_s+a_t)
    distance = nx//32-2
    assert distance == 266238 and F(31, 20)/ELL*TIME.hi < 10
    assert F(30, 16638) < F(1, 500) and 2*10**30 < 500**12
    hop = F(1, 10**30)
    platform = TIME*2*(a_s+a_t)*(2*tail+hop)
    envelope = F(31, 20)*TIME*rho
    assert source.hi < F(6, 10**6)
    assert platform.hi < F(4, 10**6)
    assert envelope.hi < F(26, 10**6)
    budget = F(56, 10**6)
    assert (source+platform+envelope+F(1, 50000)).hi < budget

    tplus = R(F('1.558104078501796'))
    tauplus, plus, errorplus = transport(1, tplus)
    tauminus, minus, errorminus = transport(-1, TIME)
    mismatch = tauplus-tauminus
    assert -F(1, 10**14) < mismatch.lo <= mismatch.hi < F(1, 10**14)
    lipschitz = F(5, 2)*CARRIER*F(31, 20)/(F(17, 20)*R(F(19, 50)).cos())
    assert lipschitz.hi < 2
    clock_upper = (1+R(F(19, 50)).sin())/2
    for tau in (tauplus, tauminus):
        probability = (1+(tau/5).sin())/2
        assert probability.lo-budget > F(1, 2)
        assert probability.hi+budget < clock_upper.lo

    phase_bound = (CARRIER-KSTAR)*F(31, 20)*TIME/4
    width = 96*PI/lx
    full_band = F(31, 20)*TIME*ELL**2*((CARRIER+width)**3+width**3)/6
    assert phase_bound.hi < F(38, 10**5)
    assert full_band.hi < F(16, 10000)
    direct_full = full_band+source+platform
    assert direct_full.hi < F(1501, 10**6)
    assert (direct_full+F(1, 50000)).hi < F(1521, 10**6)
    rows = []
    for axis, xp, xm in zip('xyz', plus, minus):
        gap = ((KSTAR*xp).sin()-(KSTAR*xm).sin())/4
        cg = ((CARRIER*xp).sin()-(CARRIER*xm).sin())/4
        strict = gap-6*budget-F(1, 10**12)
        continuous_record = strict-2*F(38, 10**5)
        assert strict.lo > F(88, 100000)
        assert continuous_record.lo > F(12, 100000)
        rows.append({'axis': axis, 'xi_plus': xp.data(), 'xi_minus': xm.data(),
                     'ideal_lattice_gap': gap.data(), 'ideal_continuous_gap': cg.data(),
                     'common_clock_lattice_lower': strict.data(),
                     'continuum_record_robust_lower': continuous_record.data()})
    result = {
        'status': 'PASS_WITH_SCOPE', 'groups': 4, 'blockers': 0,
        'new_scientific_groups': 0,
        'method': 'rational outward bounds; independent atan identity and composite Simpson',
        'implementation_independence': 'author code read first; no author imports or function reuse',
        'group1_packet': {'normalization': str(norm), 'fourier_second_moment': str(second),
                          'frequency_support': [-48, 48], 'guard_hops': distance,
                          'spatial_sites': volume, 'copies': copies,
                          'complete_Fock_baseline': str(baseline),
                          'per_platform_unknown_dimension': 48, 'arbitrary_R': True},
        'group2_error': {'source_bound': source.data(), 'tail_bound': tail.data(),
                         'platform_bound': platform.data(), 'envelope_bound': envelope.data(),
                         'common_probability_budget': str(budget),
                         'probability_lattice_to_continuum': phase_bound.data(),
                         'general_full_state_lattice_to_continuum': full_band.data()},
        'group3_clock_axes': {'pi': PI.data(), 'second_derivative_bound': str(D2),
                              'fourth_derivative_bound': str(D4), 'simpson_panels': 128,
                              'quadrature_remainders': [str(errorplus), str(errorminus)],
                              'tau_plus': tauplus.data(), 'tau_minus': tauminus.data(),
                              'clock_mismatch': mismatch.data(), 'clock_curve_lipschitz': lipschitz.data(),
                              'rows': rows},
        'group4_scope': {'negative_control': 'Gs=+Z has identical Hartree probe histories under q sign change; not identical exact BG',
                         'direct_actual_to_continuous_wide_packet': direct_full.data(),
                         'with_instrument': (direct_full+F(1, 50000)).data(),
                         'timing': 'predeclared known endpoints; no adaptive single-shot stop',
                         'no_claims': ['asynchronous 96-code handover', 'instantaneous volume rate',
                                       'non-CMC', 'Einstein equations', 'empirical implementation']}}
    return result


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bind(result):
    author = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    comparisons = []

    def compare(name, mine, theirs):
        ml, mh = F(mine['lower']), F(mine['upper'])
        al, ah = F(theirs['lower']), F(theirs['upper'])
        assert max(ml, al) <= min(mh, ah), name
        comparisons.append(name)

    for key, value in result['group2_error'].items():
        if isinstance(value, dict):
            compare(key, value, author[key])
    for key in ('tau_plus', 'tau_minus', 'clock_mismatch', 'clock_curve_lipschitz'):
        compare(key, result['group3_clock_axes'][key], author[key])
    for mine, theirs in zip(result['group3_clock_axes']['rows'], author['rows']):
        assert mine['axis'] == theirs['axis']
        for key in mine:
            if key != 'axis':
                compare(mine['axis']+'/'+key, mine[key], theirs[key])
    core_names = ['../research_note_1184.md', 'admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md']
    source_paths = ['1176/proof.md', '1178/proof.md', '1179/proof.md',
                    '../archive_1120_1158/1128/proof.md',
                    '_admission/after1174_direction/mean_field_bound.md']
    result['author_comparison'] = {'overlapping_certified_intervals': comparisons,
                                   'author_groups': author['groups']}
    result['author_core_sha256'] = {name: hash_file(HERE/name) for name in core_names}
    result['source_sha256'] = {name: hash_file(ROOT/name) for name in source_paths}
    result['reviewer_code_sha256'] = hash_file(Path(__file__))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--math-only', action='store_true')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = mathematics()
    if args.math_only:
        print('MATH_ONLY_PASS 4 groups; exact Fourier and independent Simpson certificate')
        print(json.dumps({r['axis']: r['continuum_record_robust_lower']['display_only']
                          for r in result['group3_clock_axes']['rows']}))
        return
    result = bind(result)
    path = HERE/'review_checks.json'
    if args.write:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('WRITE_ONCE_PASS 4 groups; blockers 0')
    else:
        assert result == json.loads(path.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH 4 groups; blockers 0')


if __name__ == '__main__':
    main()
