"""Exact local screening only. No ODE, new research round or device licence."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


class D:
    def __init__(self, value, derivative=0):
        self.v, self.d = F(value), F(derivative)

    def __add__(self, other):
        other = lift(other)
        return D(self.v + other.v, self.d + other.d)

    __radd__ = __add__

    def __neg__(self):
        return D(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-lift(other))

    def __rsub__(self, other):
        return lift(other) + (-self)

    def __mul__(self, other):
        other = lift(other)
        return D(self.v * other.v, self.d * other.v + self.v * other.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = lift(other)
        return D(self.v / other.v,
                 (self.d * other.v - self.v * other.d) / other.v**2)

    def __rtruediv__(self, other):
        return lift(other) / self

    def __pow__(self, power):
        return D(self.v**power, power * self.v**(power-1) * self.d) if power else D(1)


def lift(value):
    return value if isinstance(value, D) else D(value)


def source_and_lagrangian(N, a, pdot, xdot, p, x, dust_com, c, ed):
    """Per coordinate volume, common n_c,0=1; omit metric-free symplectic terms."""
    n = a**-3
    kinetic = ((pdot / N)**2 + (xdot / N)**2) / 2
    potential = (p*p + x*x) / 2
    rod = F(3, 20) * (1 + p/F(5)) / a**2
    dust = dust_com * n
    bare_rho = kinetic + potential + rod + dust
    bare_pressure = kinetic - potential - rod / 3
    lagrangian = N * (a**3 * (kinetic - potential - rod)
                      - dust_com - ed - c * bare_rho)
    rho = bare_rho + n*ed + n*c*(bare_rho - 2*kinetic)
    pressure = bare_pressure + n*c*(2*rod/3 + dust)
    return lagrangian, rho, pressure, bare_rho, kinetic


def run():
    project = next(p for p in Path(__file__).resolve().parents
                   if (p / 'research_cognition_physics').is_dir())
    stage = project / 'research_cognition_physics/archive_1234_'
    sources = {}
    for relative in ('1247/proof.md', '1249/canonical.md'):
        sources[relative] = hashlib.sha256((stage / relative).read_bytes()).hexdigest()

    c, ed = F(1,1000), F(201,20000)
    # Independent, nonzero-velocity local jet, not claimed to be the saved trajectory.
    a, pdot, xdot, p, x, dust = map(F, (2, F(2,5), F(-1,3), F(1,10), F(1,5), F(239,250)))
    lag_N, rho, pressure, bare_rho, kinetic = source_and_lagrangian(
        D(1,1), D(a), pdot, xdot, p, x, dust, c, ed)
    lag_a = source_and_lagrangian(D(1), D(a,1), pdot, xdot, p, x, dust, c, ed)[0]
    assert -lag_N.d/a**3 == rho.v
    assert lag_a.d/(3*a*a) == pressure.v
    lag_p = source_and_lagrangian(D(1), D(a), D(pdot,1), D(xdot), p, x, dust, c, ed)[0]
    lag_x = source_and_lagrangian(D(1), D(a), D(pdot), D(xdot,1), p, x, dust, c, ed)[0]
    d = 1-c/a**3
    assert lag_p.d == a**3*d*pdot
    assert lag_x.d == a**3*d*xdot
    legendre = lag_p.d*pdot + lag_x.d*xdot - lag_N.v
    assert legendre == a**3*rho.v
    naive = bare_rho.v + ed/a**3 + c/a**3*bare_rho.v
    omission = rho.v-naive
    assert omission == -2*c/a**3*kinetic.v and omission != 0

    # Original 1247 material initial jet, with a newly specified positive pointer.
    _, r0, p0, bare0, k0 = source_and_lagrangian(
        D(1), D(1), F(0), F(0), F(0), F(1,10), F(191,200), c, ed)
    assert bare0.v == F(111,100) and k0.v == 0
    assert r0.v == F(28029,25000)
    assert p0.v == F(-10789,200000)
    y_rate = F(1,100) + F(1,10)*bare0.v
    theta_rate = 1+c
    assert y_rate == F(121,1000)
    lag_J = source_and_lagrangian(
        D(1), D(1), F(0), F(0), F(0), F(1,10), D(F(191,200),1), c, ed)[0]
    kd = D(F(1,100),1)
    lag_kd = source_and_lagrangian(
        D(1), D(1), F(0), F(0), F(0), F(1,10), F(191,200),
        F(1,10)*kd, F(1,100)+kd**2/2)[0]
    assert -lag_J.d == theta_rate
    assert -lag_kd.d == y_rate
    naive_density = (y_rate/theta_rate-F(1,100))/F(1,10)
    assert naive_density != bare0.v
    assert ((y_rate/theta_rate)*theta_rate-F(1,100))/F(1,10) == bare0.v
    cone_ratio = (1+c)/(1-c)
    assert cone_ratio == F(1001,999) and cone_ratio > 1
    assert (1-c)*cone_ratio-(1+c) == 0
    assert 1-cone_ratio < 0

    return {
        'status': 'PASS_LOCAL_SCREEN_ONLY',
        'scope': 'Exact local algebra; no finite-window, quantum, causal licence or new science group.',
        'sources_sha256': sources,
        'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'nonzero_velocity_jet': {
            'rho_from_lapse': str(rho.v), 'pressure_from_scale': str(pressure.v),
            'p_momentum': str(lag_p.d), 'x_momentum': str(lag_x.d),
            'legendre_density': str(legendre/a**3),
            'missing_velocity_term_in_naive_source': str(omission)},
        'specified_positive_pointer_initial_jet': {
            'bare_rho': str(bare0.v), 'new_rho_and_H_squared': str(r0.v),
            'new_pressure': str(p0.v), 'pointer_bare_energy': str(ed),
            'pointer_proper_rate': str(y_rate), 'theta_proper_rate': str(theta_rate),
            'pointer_theta_rate': str(y_rate/theta_rate),
            'wrong_density_if_theta_is_unmodified': str(naive_density),
            'same_history_omitted_baseline_read_rate': str(F(1,10)*F(1,100)),
            'scalar_characteristic_speed_squared': str(cone_ratio),
            'g_covector_norm_for_unit_spatial_wavevector': str(1-cone_ratio)},
        'conclusion': 'The positive-coupling local scalar principal mode exceeds the original g cone; stop candidate.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true', help='Exclusive first save only.')
    args = parser.parse_args()
    result = run()
    destination = Path(__file__).with_name('results.json')
    if args.write:
        with destination.open('x', encoding='utf-8') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
        print('SAVED_LOCAL_SCREEN_ONLY')
    else:
        assert result == json.loads(destination.read_text(encoding='utf-8')), 'Saved result mismatch'
        print('READ_ONLY_MATCH_LOCAL_SCREEN_ONLY')
