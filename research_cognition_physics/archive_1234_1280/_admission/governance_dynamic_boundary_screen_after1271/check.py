"""Bounded diagnostics only; no physical apparatus or gravity-generation claim."""
from pathlib import Path
from fractions import Fraction
import argparse
import json
import math
import numpy as np

D = Path(__file__).resolve().parent


def integrate(y, x):
    return np.trapezoid(y, x) if hasattr(np, 'trapezoid') else np.trapz(y, x)


def compute():
    x = np.linspace(0., 1., 65537)
    g = math.sqrt(2.) * np.sin(math.pi*x)
    gp = math.sqrt(2.) * math.pi * np.cos(math.pi*x)
    waves = []
    means, variances, probabilities = [], [], []
    t = 1./(4.*math.pi)
    for sign in [1, -1]:
        k = sign/(2.*math.pi)
        wave = np.exp(1j*k*x)*g
        derivative = np.exp(1j*k*x)*(gp+1j*k*g)
        d_wave = -2j*math.pi*derivative
        norm = float(integrate(abs(wave)**2, x))
        mean = integrate(np.conj(wave)*d_wave, x)
        second = float(integrate(abs(d_wave)**2, x))
        assert abs(norm-1.) < 1e-12
        assert abs(mean-sign) < 1e-12
        variance = second-float(mean.real)**2
        assert abs(variance-4.*math.pi**4) < 1e-9
        # e^{itD} acts by x -> x + 1/2; integrate overlap directly.
        y = np.linspace(0., .5, 65537)
        a = np.exp(1j*k*y)*math.sqrt(2.)*np.sin(math.pi*y)
        b = np.exp(1j*k*(y+.5))*math.sqrt(2.)*np.sin(math.pi*(y+.5))
        overlap = integrate(np.conj(a)*b, y)
        assert abs(overlap-np.exp(1j*sign*t)/math.pi) < 1e-9
        probabilities.append(float((1.+overlap.imag)/2.))
        means.append(float(mean.real))
        variances.append(variance)
        waves.append(wave)
    assert np.max(abs(abs(waves[0])**2-abs(waves[1])**2)) < 1e-14
    gap = probabilities[0]-probabilities[1]
    lower = Fraction(245, 11616)
    assert lower > Fraction(1, 50)
    assert gap > float(lower)
    assert gap-2./200. > 1./100.
    # Triangle contractions: Plancherel weight = n while ||a_n|| = 1.
    weight_checks = []
    for n in [1, 2, 4]:
        length = 2.*math.pi*n
        p = np.linspace(-length, length, 65537)
        value = float(integrate(1.-abs(p)/length, p)/(2.*math.pi))
        assert abs(value-n) < 1e-12
        weight_checks.append({'n': n, 'norm': 1, 'weight': value})
    # Normal finite-window average: the L*epsilon/2 bound is attained by L*|u|.
    L, eps = Fraction(3, 2), Fraction(1, 10)
    average = L*eps/2
    assert average == Fraction(3, 40)
    return {
        'status': 'PASS_BOUNDED_DIAGNOSTICS',
        'scope': 'Formal named-cut states and bounded-query qualification; not an actual process.',
        'cut_means': means,
        'cut_variances': variances,
        'effect_probabilities': probabilities,
        'effect_gap': gap,
        'rigorous_gap_lower': str(lower),
        'physical_actual_instrument': False,
        'plancherel_contractions': weight_checks,
        'finite_window_lipschitz_bound': str(average),
        'new_rounds': 0,
        'new_axioms': 0,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = compute()
    path = D/'results.json'
    if args.write:
        with path.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write('\n')
        print('WRITE_ONCE_PASS')
    else:
        assert json.loads(path.read_text(encoding='utf-8-sig')) == result
        print('READ_ONLY_MATCH; bounded qualification only')


if __name__ == '__main__':
    main()
