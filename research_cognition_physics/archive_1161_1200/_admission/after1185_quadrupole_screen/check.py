"""Exact unnumbered quadrupole kernel screen, not an actual instrument.

Fraction polynomial jets and time integrals. No numerical quadrature, no
source-preparation/quantum-channel/xi-selection claim. Default is read-only.
"""
from fractions import Fraction as F
from math import factorial, comb
from pathlib import Path
import argparse
import hashlib
import json


def clean(poly):
    return {m: F(c) for m, c in poly.items() if c}


def add(a, b):
    result = dict(a)
    for m, c in b.items():
        result[m] = result.get(m, F(0))+c
    return clean(result)


def scale(c, a):
    return clean({m: c*x for m, x in a.items()})


def sub(a, b):
    return add(a, scale(-1, b))


def mul(a, b, degree=None):
    result = {}
    for m, c in a.items():
        for n, d in b.items():
            mn = tuple(x+y for x, y in zip(m, n))
            if degree is None or sum(mn) <= degree:
                result[mn] = result.get(mn, F(0))+c*d
    return clean(result)


def variable(index, count):
    m = [0]*count
    m[index] = 1
    return {tuple(m): F(1)}


def derivative_at_zero(poly, powers):
    return poly.get(tuple(powers), F(0))*prod(factorial(n) for n in powers)


def prod(values):
    result = 1
    for value in values:
        result *= value
    return result


def radial_jet(exponent):
    # R=1, around (0,0,1); restore powers of R by homogeneity.
    x, y, z = [variable(i, 3) for i in range(3)]
    u = add(scale(2, z), add(mul(x, x), add(mul(y, y), mul(z, z))))
    one = {(0, 0, 0): F(1)}
    result, power, coefficient = dict(one), dict(one), F(1)
    for n in range(1, 5):
        power = mul(power, u, degree=4)
        coefficient *= (exponent-(n-1))/n
        result = add(result, scale(coefficient, power))
    return result


def d_squared(poly):
    return (derivative_at_zero(poly, (4, 0, 0))
            +derivative_at_zero(poly, (0, 4, 0))
            -2*derivative_at_zero(poly, (2, 2, 0)))


def conservation_and_contact():
    k = [variable(i, 4) for i in range(3)]
    omega = variable(3, 4)
    e = [F(1), F(-1), F(0)]
    rho = scale(-1, add(mul(k[0], k[0]), scale(-1, mul(k[1], k[1]))))
    v = [scale(-e[i], mul(omega, k[i])) for i in range(3)]
    stress = [[scale(-e[i], mul(omega, omega)) if i == j else {}
               for j in range(3)] for i in range(3)]
    density_divergence = scale(-1, mul(omega, rho))
    for i in range(3):
        density_divergence = add(density_divergence, mul(k[i], v[i]))
    assert not density_divergence
    for i in range(3):
        momentum_divergence = scale(-1, mul(omega, v[i]))
        for j in range(3):
            momentum_divergence = add(momentum_divergence, mul(k[j], stress[j][i]))
        assert not momentum_divergence
    assert not add(add(stress[0][0], stress[1][1]), stress[2][2])
    k2 = {}
    for ki in k:
        k2 = add(k2, mul(ki, ki))
    s = [[sub(k2 if i == j else {}, mul(k[i], k[j]))
          for j in range(3)] for i in range(3)]
    trace = add(add(s[0][0], s[1][1]), s[2][2])
    assert not sub(trace, scale(2, k2))
    for i in range(3):
        divergence = {}
        for j in range(3):
            divergence = add(divergence, mul(k[j], s[i][j]))
        assert not divergence
        for j in range(3):
            # Transverse tensor: TT = S - P trace(S)/2, multiply by k^2.
            assert not sub(mul(k2, s[i][j]), scale(F(1, 2), mul(s[i][j], trace)))
    return {'four_divergences_zero': True, 'quadrupole_stress_trace_zero': True,
            'contact_source_density_and_current_zero': True,
            'contact_source_TT_zero': True,
            'contact_pairing': 'integral dt d3x (d_t^2 chi) rho_B'}


def spatial_coefficients():
    # Common factor 1/(4*pi); G1=r^-1 and G2=-r/2 in these units.
    g1 = radial_jet(F(-1, 2))
    g2 = scale(F(-1, 2), radial_jet(F(1, 2)))
    d4_g1 = d_squared(g1)
    direct_momentum = 2*(derivative_at_zero(g1, (2, 0, 0))
                         +derivative_at_zero(g1, (0, 2, 0)))
    longitudinal_momentum = d_squared(g2)/2
    angular_node = (derivative_at_zero(g1, (2, 0, 0))
                    -derivative_at_zero(g1, (0, 2, 0)))
    assert d4_g1 == 12
    assert direct_momentum == -4
    assert longitudinal_momentum == 1
    assert angular_node == 0
    density = d4_g1/2
    momentum_r2 = direct_momentum+longitudinal_momentum
    momentum_m2 = -d4_g1/12  # G2 average changes by -m2 G1/6.
    assert (density, momentum_r2, momentum_m2) == (6, -3, -1)
    return {'normalization': 'I_C=(4*pi*R^5)^-1 [c0 int(ab)+(cR R^2+cm m2) int(a_prime b_prime)]',
            'c0': str(density), 'cR': str(momentum_r2), 'cm': str(momentum_m2),
            'jet_D_squared_G1_in_1_over_4pi_units': str(d4_g1),
            'direct_momentum_in_1_over_4pi_units': str(direct_momentum),
            'longitudinal_momentum_in_1_over_4pi_units': str(longitudinal_momentum),
            'point_density_clock_angular_node': str(angular_node),
            'exact_biharmonic_average_coefficient': '-m2/6'}


def derivative_1d(coefficients):
    return [F(i)*c for i, c in enumerate(coefficients) if i]


def polynomial_product(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def integral_unit(coefficients):
    return sum((c/F(i+1) for i, c in enumerate(coefficients)), F(0))


def temporal_witness():
    p = list(map(F, [0, 0, 0, 1, -3, 3, -1]))
    endpoint = p
    for order in range(3):
        assert endpoint[0] == 0 and sum(endpoint) == 0
        endpoint = derivative_1d(endpoint)
    dp = derivative_1d(p)
    i0 = integral_unit(polynomial_product(p, p))
    i1 = integral_unit(polynomial_product(dp, dp))
    assert i0 == F(1, 12012) and i1 == F(1, 770)
    # f is normalized (1-r^2/a^2)^4_+; 4*pi and a powers cancel in moment ratio.
    radial_normalization = sum(((-1)**j*F(comb(4, j), 2*j+3) for j in range(5)), F(0))
    radial_moment = sum(((-1)**j*F(comb(4, j), 2*j+5) for j in range(5)), F(0))
    moment_ratio = radial_moment/radial_normalization
    assert moment_ratio == F(3, 13)
    radius, R, time, q0 = F(1, 4), F(3), F(1), F(1)
    m2 = 2*radius*radius*moment_ratio
    assert m2 == F(3, 104)
    assert 0 < time < R-2*radius
    pi_times_pairing = q0*q0/(4*R**5)*(6*time*i0-(3*R*R+m2)*i1/time)
    assert pi_times_pairing < 0
    assert 6*i0 == F(1, 2002)
    # C-infinity sign certificate: pi>3, T<R, m2>=0 imply
    # 6-(3 R^2+m2) pi^2/T^2 < 6-3*9 = -21.
    assert F(6)-3*3**2 == -21
    return {'p_coefficients_ascending': [str(x) for x in p],
            'endpoint_derivatives_0_through_2_zero': True,
            'integral_p_squared': str(i0), 'integral_derivative_squared': str(i1),
            'poincare_ratio_for_polynomial': str(i1/i0),
            'radial_moment_over_radius_squared': str(moment_ratio),
            'fixed_example': {'R': str(R), 'rA': str(radius), 'rB': str(radius),
                              'T': str(time), 'q0': str(q0), 'm2': str(m2),
                              'pi_times_I_C_exact': str(pi_times_pairing),
                              'strictly_negative': True,
                              'whole_supports_spacelike': True},
            'C_infinity_sign': 'Poincare; pi>3 and T<R imply normalized bracket < -21',
            'actual_quantum_source_or_instrument_claimed': False}


HERE = Path(__file__).resolve().parent


def calculate():
    return {'status': 'UNNUMBERED_MATHEMATICAL_KERNEL_SCREEN_PASS',
            'date': '2026-10-13', 'new_research_rounds': 0,
            'conservation_and_contact': conservation_and_contact(),
            'complete_constraint_pairing': spatial_coefficients(),
            'finite_temporal_witness': temporal_witness(),
            'scope': {'complete_Cg_used': True, 'trace_zero_specific_to_chosen_e': True,
                      'GR_retarded_complete_pairing_spacelike': '0',
                      'deformed_pairing': '(xi-1) I_C',
                      'xi_selected': False, 'actual_preparation_certified': False,
                      'unknown_input_R_channel_certified': False,
                      'Einstein_derived': False},
            'proof_sha256': hashlib.sha256((HERE/'proof.md').read_bytes()).hexdigest(),
            'check_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write:
        with path.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print('WRITE_ONCE_PASS: unnumbered mathematical kernel screen')
    else:
        assert result == json.loads(path.read_text(encoding='utf-8')), 'READ_ONLY_MISMATCH'
        print('READ_ONLY_MATCH: unnumbered mathematical kernel screen')
    print('pi * I_C =', result['finite_temporal_witness']['fixed_example']['pi_times_I_C_exact'])


if __name__ == '__main__':
    main()
