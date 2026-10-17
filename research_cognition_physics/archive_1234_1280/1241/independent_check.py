"""1241 independent finite-state/scale audit, no actual detector claim.

Uses uniform composite Simpson Fourier integration, one-dimensional radial
Maxwell integration (not the author's product Hermite rule), and Fraction
enclosures for the scale certificate. Default only compares frozen results.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'independent_results.json'
WATCH = [HERE/n for n in ['admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md']]


def need(value, msg):
    if not bool(value):
        raise AssertionError(msg)


def rat(x):
    return dict(fraction=str(x), decimal=float(x))


def simpson(values, step):
    return step/3*(values[0]+values[-1]+4*np.sum(values[1:-1:2])+2*np.sum(values[2:-1:2]))


def pi_interval():
    # Machin formula, alternating rational arctangent series.
    def atan(z, n=40):
        total = sum(((-1)**k*z**(2*k+1)/F(2*k+1) for k in range(n)), F(0))
        remainder = z**(2*n+1)/F(2*n+1)
        return (total, total+remainder) if n % 2 == 0 else (total-remainder, total)
    a, b = atan(F(1, 5)), atan(F(1, 239))
    lo, hi = 16*a[0]-4*b[1], 16*a[1]-4*b[0]
    # Outward decimal rounding keeps the saved exact certificate readable.
    scale = 10**40
    lower = F((lo.numerator*scale)//lo.denominator, scale)
    upper = F((hi.numerator*scale+hi.denominator-1)//hi.denominator, scale)
    need(lower <= lo <= hi <= upper, 'outward pi rounding')
    return lower, upper


def sqrt_interval(q, digits=40):
    scale = 10**digits
    n = math.isqrt(q.numerator*scale*scale//q.denominator)
    lo, hi = F(n, scale), F(n+1, scale)
    need(lo*lo <= q < hi*hi, 'exact rational square-root enclosure')
    return lo, hi


def packet_group():
    rows = []
    points = np.array([0., .7, math.pi, -math.pi, 2*math.pi, 7., 13.])
    exact_fourier = []
    for u in points:
        if abs(abs(u)-math.pi) < 1e-14:
            val = -1j*np.sign(u)/2
        else:
            val = math.pi*(1+np.exp(-1j*u))/(math.pi**2-u*u)
        exact_fourier.append(val/math.sqrt(math.pi))
    for n in [8192, 16384]:
        x = np.linspace(0., 1., n+1)
        f = np.sqrt(2)*np.sin(np.pi*x)
        df = np.sqrt(2)*np.pi*np.cos(np.pi*x)
        values = [float(simpson(x**k*f*f, 1/n)) for k in range(3)]
        variance = float(simpson(df*df, 1/n))
        transform = np.array([simpson(f*np.exp(-1j*u*x), 1/n)/np.sqrt(2*np.pi) for u in points])
        err = float(np.max(abs(transform-np.array(exact_fourier))))
        expected = [1., .5, 1/3-1/(2*np.pi**2)]
        need(max(abs(a-b) for a,b in zip(values, expected)) < 1e-13, 'spectral normalization and first two moments')
        need(abs(variance-np.pi**2) < 1e-12 and err < 2e-13, 'H1 time moment and direct Fourier amplitude')
        rows.append(dict(subintervals=n, energy_moments=values, time_second_moment=variance,
                         Fourier_error=err))
    # Direct real-time density from the algebraic Fourier transform; use the
    # removable values only at the two exact singular points, if any occur.
    U, n = 512., 2**18
    t = np.linspace(-U, U, n+1)
    denominator = np.pi**2-t*t
    density = 4*np.pi*np.cos(t/2)**2/denominator**2
    near = abs(abs(t)-np.pi) < 1e-12
    density[near] = 1/(4*np.pi)
    integral = float(simpson(density, 2*U/n))
    second = float(simpson(density*t*t, 2*U/n))
    # For |t|>=U>=2pi, p(t)<=64pi/(9t^4). This explicit omitted
    # tail and Plancherel, rather than finite sampling, normalize the POVM.
    tail0, tail2 = 128*np.pi/(27*U**3), 128*np.pi/(9*U)
    need(integral < 1+1e-12 and integral+tail0 > 1, 'direct probability and omitted tail')
    need(second < np.pi**2 < second+tail2, 'direct time moment and omitted tail')
    need(np.max(density) <= 4/np.pi**3+1e-13, 'finite diagnostic of analytic density supremum')
    pilo, pihi = pi_interval()
    need(pihi-pilo < F(1, 10**38), 'tight rational pi interval')
    # The success upper bound used in deletion follows from the exact L1 norm
    # sqrt(2)*2/pi, not from the numerical density maximum.
    success = 8*F(10**4)*F(1,10**8)/pilo**3
    need(success < F(1, 10000), 'small-energy state-family success upper bound')
    return dict(status='PASS', independent_method='composite Simpson plus algebraic Fourier formula',
                frequency_checks=rows, direct_time_cutoff=U, direct_time_probability=integral,
                omitted_probability_upper=tail0, direct_time_second=second,
                omitted_second_upper=tail2, pi_enclosure=[rat(pilo),rat(pihi)],
                low_energy_success_upper=rat(success),
                canonical_time_POVM_is_not_a_built_arrival_detector=True)


def radial_gaussian(chi, n):
    T = .82
    s = np.linspace(0., 16., n+1)
    pdf = np.sqrt(2/np.pi)*s*s*np.exp(-s*s/2)
    v = chi/(2*T)
    p2 = v*s*s
    # Stable global spectral dispersion difference, dimensionless Mc^2=1.
    gap = p2*p2/(4*(1+p2/2+np.sqrt(1+p2)))
    norm2 = simpson(pdf*4*np.sin(T*gap/(2*chi))**2, 16/n)
    p8 = simpson(pdf*p2**4, 16/n)
    kinetic = simpson(pdf*p2/(np.sqrt(1+p2)+1), 16/n)
    norm = float(np.sqrt(norm2))
    bound = np.sqrt(945)*chi/(32*T)
    need(np.all(gap >= 0) and np.all(gap <= p2*p2/8+1e-18), 'all-spectrum remainder identity')
    need(abs(float(p8)/(945*v**4)-1) < 1e-12, '3D radial eighth moment')
    need(0 < norm < bound and 0 < kinetic < 3*chi/(4*T), 'relativistic norm and total kinetic upper bounds')
    return dict(chi=chi, radial_subintervals=n, p8=float(p8), analytic_p8=945*v**4,
                norm=norm, norm_upper=float(bound), mean_kinetic=float(kinetic),
                mean_kinetic_upper=3*chi/(4*T))


def endpoint_window_group():
    # Maxwell moment E[s^8]=(3)(5)(7)(9)=945: this is not the 1D 105 moment.
    moment = math.prod(range(3, 10, 2))
    need(moment == 945, 'three-dimensional Gaussian moment identity')
    radials = []
    for chi in [.02, .001]:
        coarse, fine = radial_gaussian(chi, 8192), radial_gaussian(chi, 16384)
        need(abs(coarse['norm']-fine['norm']) < 1e-13, 'independent radial resolution diagnostic')
        radials.append(fine)
    # Maxwell omitted probability at s=16 is bounded without assigning zero
    # to the tail: I_2=A exp(-A^2/2)+I_0, I_0<=exp(-A^2/2)/A.
    tail_prob = math.sqrt(2/math.pi)*(16+1/16)*math.exp(-128)
    need(4*tail_prob < 1e-52, 'bounded numerical radial omission')
    T, r, b, zeta = F(41,50), F(1,10**8), F(1,1000), F(1,10**8)
    dt, dx, pt, px = F(1,10**8), F(1,10**9), F(1,10**4), F(1,10**4)
    pilo, pihi = pi_interval()
    rootlo, roothi = sqrt_interval(F(945))
    crlo, crhi = rootlo/(32*T), roothi/(32*T)
    time_lo = r*dt*F(1,100)/pihi
    position_lo, position_hi = px/(3*T/dx**2+crhi), px/(3*T/dx**2+crlo)
    need(position_hi < time_lo, 'localization sets the stated sufficient chi threshold')
    mustar = b*zeta/2
    # Strictly conservative threshold based on the lower admissible chi.
    Lbound = sqrt_interval(1/(position_lo*mustar))[1]
    ell, m, B = F(10**35), F(10**5), F(10**12)
    chi, mu = 1/(m*ell), m/ell
    need(chi*mu == 1/ell**2, 'same fixed hbar,c,G product, not independent tuning')
    masslo, masshi = 1/(ell*position_lo), mustar*ell
    time_fail = pihi**2/(B*dt)**2
    position_fail = (3*T/dx**2+crhi)*chi
    compactness = 2*mu/b
    energy_cap = chi*B
    need(2*b < F(3,100) and ell > Lbound and masslo < m < masshi,
         'same finite scale lies in the conservative mass interval')
    need(energy_cap < r and time_fail < pt and position_fail < px and compactness < zeta,
         'all explicit preparation/screening inequalities')
    need(chi*T/2+chi*T/2 == chi*T, 'maximum uncorrelated free variance includes spreading')
    # A size choice is distinct from the COM variance; no gravitational
    # observable error is inferred from this dimensionless compactness.
    heaviness = 2*(m*10**20/ell)/b
    need(heaviness > zeta, 'dropping baseline mass falsely admits this heavy body')
    small_L = F(10**10)
    small_lo = 1/(small_L*position_hi)
    small_hi = mustar*small_L
    need(small_lo > small_hi, 'even relaxed sufficient interval is empty at selected small L')
    # Physical background identity after the root's unit clarification:
    # H0=2c/(3L), rho_mass=c^2/(6pi G L^2), so Friedmann's ratio is 1.
    need(F(8,6) == F(3)*F(4,9), 'mass-density Friedmann scaling coefficient')
    return dict(status='PASS', independent_method='radial Maxwell quadrature and rational interval arithmetic',
                Gaussian_radial_cases=radials, radial_norm_squared_omission_upper=4*tail_prob,
                C_R_interval=[rat(crlo),rat(crhi)], chi_position_interval=[rat(position_lo),rat(position_hi)],
                conservative_length_threshold=rat(Lbound), chi_witness=rat(chi),mu_witness=rat(mu),
                conservative_mass_interval=[rat(masslo),rat(masshi)], energy_cap_ratio=rat(energy_cap),
                time_window_failure_upper=rat(time_fail), single_time_position_failure_upper=rat(position_fail),
                body_compactness=rat(compactness), free_relativistic_norm_upper=rat(crhi*chi),
                complete_mean_COM_kinetic_upper=rat(3*chi/(4*T)),
                heavy_body_screen_deletion=rat(heaviness),
                small_L_sufficient_certificate_interval=[rat(small_lo),rat(small_hi)],
                neither_empty_certificate_nor_compactness_is_a_physical_no_go=True,
                no_joint_repeated_measurement_or_curved_source_response_certified=True)


def hashes():
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in WATCH}


def run():
    before = hashes()
    result = dict(round=1241,status='PASS',independent_scientific_groups=2,
                  positive_energy_time_preparation=packet_group(),
                  endpoint_and_fixed_constant_window=endpoint_window_group(),source_sha256=before,
                  actual_arrival_instrument=False,actual_multiple_measurement_history=False,
                  gravitational_return_error_bound=False, CT='candidate, not adopted')
    need(before == hashes(), 'frozen author sources unchanged')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--diagnostic',action='store_true')
    args=parser.parse_args()
    need(not(args.write and args.diagnostic),'choose one output mode')
    result=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('CREATED_ONCE: 1241 independent, two groups PASS')
    elif args.diagnostic:
        print(json.dumps(result,ensure_ascii=False,indent=2))
    else:
        need(result==json.loads(OUT.read_text('utf-8')),'read-only data/source match')
        print('READ_ONLY_MATCH: 1241 independent, two groups PASS')
