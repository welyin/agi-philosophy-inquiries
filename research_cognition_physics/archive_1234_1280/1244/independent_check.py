"""Independent 1244 audit: rational radar budget and convex-integral calibration.

No author code is imported.  The transcendental enclosures use Bernoulli
repeated squaring and midpoint/trapezoid integration, not the author's Taylor
series.  --write creates the initial result exclusively; default is read-only.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
GRID = 10**45


def down(q):
    return F((q.numerator * GRID) // q.denominator, GRID)


def up(q):
    return -down(-q)


def pack(q):
    return {"rational": str(q), "decimal_diagnostic": float(q)}


class Box:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi

    def scale(self, q):
        q = F(q)
        return Box(min(q*self.lo, q*self.hi), max(q*self.lo, q*self.hi))

    def add(self, other):
        if not isinstance(other, Box):
            other = Box(other)
        return Box(self.lo+other.lo, self.hi+other.hi)

    def sub(self, other):
        if not isinstance(other, Box):
            other = Box(other)
        return self.add(other.scale(-1))

    def output(self):
        return {"lower": pack(self.lo), "upper": pack(self.hi),
                "width": pack(self.hi-self.lo)}


def bernoulli_exp(z, squarings=40):
    """(1+z/n)^n <= exp(z) <= (1-z/n)^(-n), n=2**squarings.

    Every squaring is outward rounded. Positivity makes interval squaring
    monotone. The starting bounds follow from log(1+x)<=x<=-log(1-x).
    """
    n = 2**squarings
    assert 0 <= z < n
    lower, upper = down(1+z/n), up(1/(1-z/n))
    for _ in range(squarings):
        lower, upper = down(lower*lower), up(upper*upper)
    return Box(lower, upper)


def tanh_box(z):
    b = bernoulli_exp(2*z)
    return Box(down((b.lo-1)/(b.lo+1)), up((b.hi-1)/(b.hi+1)))


def atanh_by_convex_integral(z, cells=512):
    """Integral of f(t)=1/(1-t*t): midpoint <= integral <= trapezoid.

    f''=(2+6t*t)/(1-t*t)^3 > 0 on [0,1). For interval endpoints,
    monotonicity of atanh permits the lower and upper integrals separately.
    Rounding each summand to a fixed rational grid avoids denominator growth.
    """
    assert 0 <= z.lo <= z.hi < 1
    lower_sum = F(0)
    for j in range(cells):
        t = z.lo * F(2*j+1, 2*cells)
        lower_sum += down(1/(1-t*t))
    upper_sum = up((1+1/(1-z.hi*z.hi))/2)
    for j in range(1, cells):
        t = z.hi * F(j, cells)
        upper_sum += up(1/(1-t*t))
    return Box(down(z.lo*lower_sum/cells), up(z.hi*upper_sum/cells))


def finite_budget():
    eps, G, T, L = F(1, 10000), F(1, 1000), F(1, 2), F(1, 10)
    eta_clock, eta_direction, eta_model = F(1, 100000), F(1, 10000), F(1, 100000)
    domain_guards = []
    for epsilon in (F(0), eps, F(1, 4)):
        cplus2 = (1+epsilon)/(1-epsilon)
        cminus2 = 1/cplus2
        assert cplus2 <= (1+2*epsilon)**2
        assert cminus2 >= (1-2*epsilon)**2
        assert cplus2 < F(9, 4)
        assert 1/(1-epsilon*epsilon) <= 4
        domain_guards.append({"epsilon": pack(epsilon), "cplus_squared": pack(cplus2)})
    assert 8*G*T < 1
    # sqrt(3)<7/4, |v|<3/2 imply ||W||_1<4 and the 56G bound.
    assert F(7, 4)**2 > 3
    assert 1+F(7, 4)*F(3, 2) < 4
    acceleration_over_G = (F(7, 4)+F(3, 2))*16
    assert acceleration_over_G < 56
    frame_bound = 8*G*T/(1-8*G*T)
    actual_exp = bernoulli_exp(4*G*T)
    assert actual_exp.hi-1 < frame_bound
    dmax = L*(1+2*eps)
    geometry_terms = (4*eps*L, dmax*frame_bound, 28*G*dmax*dmax)
    geometric = sum(geometry_terms, F(0))
    instrument = 2*eta_clock+dmax*eta_direction+eta_model
    pair = 3*(geometric+instrument)
    assert geometric == F(37443310581, 51875000000000)
    assert pair == F(118555242993, 51875000000000)
    author = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    assert pair == F(author['groups'][0]['two_observer_error_bound']['exact'])

    # Same flat geometry and exact radar times, but a reversed direction gives
    # a finite position error. The removed contract is the direction identity.
    emission, reception = F(1, 10), F(3, 10)
    radius, radar_time = (reception-emission)/2, (emission+reception)/2
    direction_delete_squared = (2*radius)**2
    assert direction_delete_squared > 0
    # Identity observers and q on the common axis; a twice-running B display
    # gives a nonzero temporal discrepancy despite epsilon=G=0.
    clock_delete_error = radar_time
    assert clock_delete_error > 0
    # Boosts act on the null Euclidean vector (1,-1) by the factor D exactly.
    amplification = []
    for D in (F(2), F(4), F(8)):
        gamma, gamma_v = (D+1/D)/2, (D-1/D)/2
        assert gamma*gamma-gamma_v*gamma_v == 1
        out = (gamma+gamma_v, -gamma_v-gamma)
        assert out == (D, -D)
        amplification.append({"D": pack(D), "squared_norm_ratio": pack(D*D)})
    return {
        'status': 'PASS', 'cone_guards': domain_guards,
        'acceleration_constant_safe_upper': pack(acceleration_over_G),
        'frame_exponential_independent_enclosure': actual_exp.output(),
        'frame_rational_bound': pack(frame_bound), 'dmax': pack(dmax),
        'geometry_terms': [pack(q) for q in geometry_terms],
        'geometric_error': pack(geometric), 'instrument_error': pack(instrument),
        'two_observer_error': pack(pair), 'author_exact_budget_match': True,
        'deletions': {
            'flat_direction_identity_removed_error_squared': pack(direction_delete_squared),
            'flat_clock_calibration_removed_error': pack(clock_delete_error),
            'unbounded_rapidity_amplification': amplification,
        },
        'scope': 'Declared finite tube/readout contracts only, not an implemented two-observer apparatus.'
    }


def de_sitter_case(H, emission, D, compare_author=False):
    gamma, gamma_v = (D+1/D)/2, (D-1/D)/2
    assert gamma*gamma-gamma_v*gamma_v == 1
    assert D*D-2*gamma*D+1 == 0
    u = tanh_box(H*emission/2)
    assert D*D*u.hi < 1
    reception = atanh_by_convex_integral(u.scale(D)).scale(2/H)
    returned = atanh_by_convex_integral(u.scale(D*D)).scale(2/H)
    assert reception.hi < returned.lo
    transported_t = returned.scale(1/(2*D)).add(D*emission/2)
    transported_x = returned.scale(1/(2*D)).sub(D*emission/2)
    dt, dx = transported_t.sub(reception), transported_x
    assert dt.lo > 0 and dx.lo > 0
    assert returned.hi-returned.lo < F(1, 10**9)
    assert reception.hi-reception.lo < F(1, 10**9)
    # Independent analytic positivity: first nonzero positive atanh-series term.
    # Exact spatial mismatch = [atanh(D^2 u)-D^2 atanh(u)]/(H D).
    spatial_positive_lower = (D**5-D)*u.lo**3/(3*H)
    assert 0 < spatial_positive_lower < dx.hi
    # Exact normal-coordinate relation of q on B, checked symbolically in tau_B.
    normal_time_coefficient = gamma*gamma-gamma_v*gamma_v
    normal_space_coefficient = -gamma_v*gamma+gamma*gamma_v
    assert (normal_time_coefficient, normal_space_coefficient) == (1, 0)
    # The H=0 tangent-space calibration has no discrepancy.
    flat_B, flat_R = D*emission, D*D*emission
    assert D*emission/2+flat_R/(2*D)-flat_B == 0
    assert flat_R/(2*D)-D*emission/2 == 0
    author_included = None
    if compare_author:
        author = json.loads((HERE/'results.json').read_text(encoding='utf-8'))['groups'][1]
        author_included = {}
        for name, b in (('B_reception_proper_time', reception), ('A_return_proper_time', returned)):
            lo, hi = (F(author[name][key]['exact']) for key in ('lower', 'upper'))
            assert b.lo <= lo <= hi <= b.hi
            author_included[name] = True
    return {
        'H': pack(H), 'emission': pack(emission), 'D': pack(D),
        'relative_velocity': pack(gamma_v/gamma), 'tanh_enclosure': u.output(),
        'B_reception': reception.output(), 'A_return': returned.output(),
        'transported_radar_minus_B': {'time': dt.output(), 'space': dx.output()},
        'analytic_positive_spatial_lower': pack(spatial_positive_lower),
        'exact_RNC_coefficients': [pack(normal_time_coefficient), pack(normal_space_coefficient)],
        'flat_limit_discrepancy': pack(F(0)), 'author_intervals_contained': author_included,
    }


def curved_calibrations():
    cases = [de_sitter_case(F(1, 10), F(1, 10), F(2), True),
             de_sitter_case(F(1, 5), F(1, 10), F(3, 2))]
    return {'status': 'PASS', 'cases': cases,
            'certification': {'exp_method': 'Bernoulli product enclosure with 40 repeated squarings',
                              'atanh_method': 'convex-integrand midpoint lower / trapezoid upper; 512 cells',
                              'outward_grid': str(GRID)},
            'scope': 'Specified curved geometry and ideal clock/null contracts; neither case certifies group 1 domain bounds or full source/feedback.'}


def compute():
    names = ['admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md']
    return {'round': 1244, 'audit': 'independent', 'status': 'PASS', 'science_groups': 2,
            'author_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
            'independent_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'groups': {'finite_radar_budget': finite_budget(), 'curved_calibrations': curved_calibrations()},
            'CT_adopted': False, 'complete_apparatus_realization': False,
            'full_active_SR_R_proved': False, 'GR_dynamics_derived': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = compute()
    path = HERE/'independent_results.json'
    if args.write:
        with path.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('CREATED_ONCE PASS independent science_groups=2')
    else:
        saved = json.loads(path.read_text(encoding='utf-8'))
        assert saved == result, 'Independent results or bound author files changed; no writes performed'
        print('READ_ONLY_MATCH PASS independent science_groups=2')
    print('two_observer_bound', result['groups']['finite_radar_budget']['two_observer_error']['rational'])


if __name__ == '__main__':
    main()
