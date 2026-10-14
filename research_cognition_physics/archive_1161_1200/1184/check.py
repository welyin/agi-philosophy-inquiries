"""Finite dual-platform certificate. Default is read-only; --write saves results.

All certification uses outward rounded integer intervals, not NumPy floats.
No gigantic many-body state is simulated. See proof.md for the operator bounds.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt, factorial, comb
import argparse
import json

HERE = Path(__file__).resolve().parent
S = 10**36


def ceildiv(a, b):
    return -((-a)//b)


class I:
    def __init__(self, x=0, hi=None, raw=False):
        if raw:
            self.lo, self.hi = x, hi
        elif isinstance(x, I):
            self.lo, self.hi = x.lo, x.hi
        else:
            x = F(x)
            self.lo = x.numerator*S//x.denominator
            self.hi = ceildiv(x.numerator*S, x.denominator)
        assert self.lo <= self.hi

    def __add__(self, other):
        o = I(other)
        return I(self.lo+o.lo, self.hi+o.hi, True)
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo, True)

    def __sub__(self, other):
        return self+-I(other)

    def __rsub__(self, other):
        return I(other)+-self

    def __mul__(self, other):
        o = I(other)
        p = [a*b for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(p)//S, ceildiv(max(p), S), True)
    __rmul__ = __mul__

    def reciprocal(self):
        assert self.lo > 0 or self.hi < 0
        return I(S*S//self.hi, ceildiv(S*S, self.lo), True)

    def __truediv__(self, other):
        return self*I(other).reciprocal()

    def __rtruediv__(self, other):
        return I(other)*self.reciprocal()

    def __pow__(self, n):
        assert n >= 0
        v = I(1)
        for _ in range(n):
            v = v*self
        return v

    def sqrt(self):
        assert self.lo >= 0
        low = isqrt(self.lo*S)
        high = isqrt(self.hi*S)
        if high*high < self.hi*S:
            high += 1
        return I(low, high, True)

    def padded(self, r):
        r = I(r).hi
        return I(self.lo-r, self.hi+r, True)

    def sin(self):
        assert max(abs(self.lo), abs(self.hi)) <= 3*S
        z, term, total = self*self, self, self
        for n in range(1, 30):
            term = -term*z/((2*n)*(2*n+1))
            total = total+term
        return total.padded(F(3**61, factorial(61)))

    def cos(self):
        assert max(abs(self.lo), abs(self.hi)) <= 3*S
        z, term, total = self*self, I(1), I(1)
        for n in range(1, 30):
            term = -term*z/((2*n-1)*(2*n))
            total = total+term
        return total.padded(F(3**60, factorial(60)))

    def exp(self):
        assert 0 <= self.lo <= self.hi <= 2*S
        term, total = I(1), I(1)
        for n in range(1, 81):
            term = term*self/n
            total = total+term
        return total.padded(F(2*2**81, factorial(81)))

    def below(self, x):
        return F(self.hi, S) < F(x)

    def above(self, x):
        return F(self.lo, S) > F(x)

    def dump(self):
        # Exact rational fields are authoritative. Floats are display only.
        return dict(lower=str(F(self.lo, S)), upper=str(F(self.hi, S)),
                    display=[self.lo/S, self.hi/S])


def arctan_recip(d):
    v = sum(((-1)**n*F(1, (2*n+1)*d**(2*n+1)) for n in range(80)), F(0))
    return I(v).padded(F(1, 161*d**161))


PI = 16*arctan_recip(5)-4*arctan_recip(239)
ELL, EPS, NU, G, OM = map(I, (F(16,65), F(3,20), F(7,20), F(1,5), F(1,2)))
T, P = PI/2, PI/8
KSTAR = (P*ELL).sin()/ELL


def tau(q, t):
    b = -EPS*OM*q
    r2 = OM*OM+b*b
    w = 2*r2.sqrt()
    return t-EPS*OM*b/r2*(t-(w*t).sin()/w)


def clocks_and_axes(q, t, n=1024):
    proper = tau(q, t)
    # Composite midpoint: |(cos tau)''| <= N^2+|N'| < 3/2.
    # Error <= (3/2)*t^3/(24*n^2); all evaluations are intervals.
    integral = sum((tau(q, t*F(2*j+1,2*n)).cos() for j in range(n)), I(0))*t/n
    integral = integral.padded(F(3,2)*F(t.hi,S)**3/(24*n*n))
    z = t+NU/EPS*(proper-t)
    shape = G*(1-I(q)*q).sqrt()*integral
    axes = [z+shape, z-shape, z]
    return proper, axes


def compute():
    assert PI.above('3.1415926535897932384626433832794')
    assert PI.below('3.1415926535897932384626433832796')
    Lx, Nx, K = 16*131072, 65*131072, 10**16
    assert F(Lx,Nx) == F(16,65)
    normalization = F(comb(192,96),4**96)
    assert normalization > F(1,20)
    assert ((PI*F(3,16)).cos()**2).below(F(173,250))
    tail = (I(20)*I(F(173,250))**96).sqrt()
    assert tail.below(F(1,10**7))
    # Four carrier envelopes have exact disjoint Fourier support.
    carrier_x = 131072
    assert 2*48 < carrier_x and carrier_x+48 < Nx//2
    guard_hops = Nx//32-2
    assert guard_hops >= 16638
    # JT < 10 and e < 3 imply (e JT/d)^d < (1/500)^d.
    assert (I(F(31,20))/ELL*T).below(10)
    assert F(30,16638) < F(1,500)
    # A coarse rational bound replaces the astronomically tiny exact tail.
    hop_tail = I(F(1,10**30))
    assert 2*10**30 < 500**12 and guard_hops >= 12
    a_s, a_t = EPS*F(1,2)+3*NU/ELL, 2*G/ELL
    B = (7*a_s*a_s+5*a_t*a_t).sqrt()
    source = 4*T*(6*I(F(3,20))*T).exp()/(I(K).sqrt())*((I(2).sqrt())*B+a_s+a_t)
    assert source.below(F(6,10**6))
    # Full code tail <= 2*tail, Hamiltonian difference <= 2(a_s+a_t).
    platform = T*2*(a_s+a_t)*(2*tail+hop_tail)
    assert platform.below(F(4,10**6))
    rho = 96*PI/(Lx*I(191).sqrt())
    envelope = I(F(31,20))*T*rho
    assert envelope.below(F(26,10**6))
    delta = F(56,10**6)  # 6 + 4 + 26 + 20 micro: source/platform/envelope/instrument.
    phase_lattice = (P-KSTAR)*I(F(31,20))*T/4
    assert phase_lattice.below(F(38,10**5))
    # Full-state continuum error is separate from the selected probability bound.
    width = 96*PI/Lx
    full_cont = I(F(31,20))*T*ELL**2*((P+width)**3+width**3)/6
    assert full_cont.below(F(16,10000))
    # Endpoint is a predeclared decimal rational, not an ideal adaptive stop.
    tplus = I('1.558104078501796')
    tminus = T
    properplus, xiplus = clocks_and_axes(F(1,2), tplus)
    properminus, ximinus = clocks_and_axes(F(-1,2), tminus)
    mismatch = properplus-properminus
    assert mismatch.above(F(-1,10**14)) and mismatch.below(F(1,10**14))
    # Matching at a common clock probability: entire enlarged proper-time window.
    lip = I(F(5,2))*P*I(F(31,20))/(I(F(17,20))*I(F(38,100)).cos())
    assert lip.below(2)
    clocks = [(1+(p/5).sin())/2 for p in (properplus,properminus)]
    for clock in clocks:
        assert (clock-delta).above(F(1,2))
        assert (clock+delta).below(F((1+I(F(38,100)).sin()).lo,2*S))
    rows = []
    for axis, xp, xm in zip('xyz', xiplus, ximinus):
        lattice_gap = ((KSTAR*xp).sin()-(KSTAR*xm).sin())/4
        cont_gap = ((P*xp).sin()-(P*xm).sin())/4
        # 1e-12 covers common-clock mismatch, far above the proved endpoint error.
        certified_lattice = lattice_gap-6*delta-F(1,10**12)
        certified_cont = lattice_gap-6*delta-2*F(38,10**5)-F(1,10**12)
        assert certified_lattice.above(F(88,100000))
        assert certified_cont.above(F(12,100000))
        rows.append(dict(axis=axis, xi_plus=xp.dump(), xi_minus=xm.dump(),
                         ideal_lattice_gap=lattice_gap.dump(), ideal_continuous_gap=cont_gap.dump(),
                         common_clock_lattice_lower=certified_lattice.dump(),
                         continuum_record_robust_lower=certified_cont.dump()))
    return dict(schema='round1184_dual_platform_v1', groups=4,
                certification='integer outward intervals; exact rational fields authoritative',
                resources=dict(Lx=Lx, Nx=Nx, transverse_sites=65, copies=K,
                               spatial_sites=Nx*65**2, per_platform_unknown_dimension=48,
                               guard_hops_lower=guard_hops, scalar_carrier_mode=carrier_x),
                packet_normalization=str(normalization),
                source_bound=source.dump(), tail_bound=tail.dump(), platform_bound=platform.dump(),
                envelope_bound=envelope.dump(), instrument_input='1/50000',
                common_probability_budget=str(delta), probability_lattice_to_continuum=phase_lattice.dump(),
                general_full_state_lattice_to_continuum=full_cont.dump(),
                tplus=str(F('1.558104078501796')), tminus=T.dump(),
                tau_plus=properplus.dump(), tau_minus=properminus.dump(),
                clock_mismatch=mismatch.dump(), clock_curve_lipschitz=lip.dump(),
                rows=rows, negative_control='Gs=+Z: z_s even in b; identical Hartree clock/axis histories only',
                limits=['new finite preparation/timing/readout contract', 'no 96-code asynchronous handover',
                        'cumulative transport, not instantaneous volume rate or non-CMC',
                        'no Einstein dynamics or dark-sector identification', 'no empirical experiment'])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    result = compute()
    target = HERE/'results.json'
    if args.write:
        if target.exists():
            raise SystemExit('Refusing to overwrite an existing result')
        target.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        print('WROTE 4 groups; integer interval certificate')
    else:
        assert json.loads(target.read_text(encoding='utf-8')) == result
        print('READ_ONLY_MATCH 4 groups; integer interval certificate')


if __name__ == '__main__':
    main()
