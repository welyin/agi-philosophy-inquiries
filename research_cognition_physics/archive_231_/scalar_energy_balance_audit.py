"""Round 303: Hilbert stress, action variation and averaged scalar energy balance."""
import unittest
import numpy as np
from conformal_volume_geometry_audit import de_sitter_jet
from curved_scalar_propagation_audit import exact_solution, field_from_chi
from growing_stream_audit import main


def quadrature(function, left=-4., right=-1., order=64):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x = (right-left)*nodes/2+(right+left)/2
    return float((right-left)/2*np.dot(weights, [function(t) for t in x]))


def averaged_stress(eta, qjet, xi, mass=0., k=2., hubble=0.5):
    """Average Re(q exp(ikx)) over one period; p is spatial trace/3."""
    a, da, dda = de_sitter_jet(eta, hubble)
    h, dh = da/a, dda/a-(da/a)**2
    q, dq, ddq = qjet
    f, df = abs(q)**2/2, float(np.real(q.conjugate()*dq))
    ddf = abs(dq)**2+float(np.real(q.conjugate()*ddq))
    kinetic, dkinetic = abs(dq)**2/2, float(np.real(dq.conjugate()*ddq))
    rho0 = (kinetic+k*k*f)/(2*a*a)+mass*mass*f/2
    p0 = (kinetic-k*k*f/3)/(2*a*a)-mass*mass*f/2
    drho0 = (dkinetic+k*k*df-2*h*(kinetic+k*k*f))/(2*a*a)+mass*mass*df/2
    extra = 3*h*h*f+3*h*df
    dextra = 6*h*dh*f+3*h*h*df+3*dh*df+3*h*ddf
    rho = rho0+xi*extra/a**2
    pressure = p0+xi*(-(2*dh+h*h)*f-ddf-h*df)/a**2
    drho = drho0+xi*(dextra-2*h*extra)/a**2
    scalar = 6*dda/a**3
    residual = ddq+2*h*dq+(k*k+a*a*(mass*mass+xi*scalar))*q
    return {'rho': float(rho), 'p': float(pressure), 'drho': float(drho),
            'balance': float(drho+3*h*(rho+pressure)),
            'off_shell_source': float(np.real(residual*dq.conjugate())/(2*a*a)),
            'minimal_only_balance': float(drho0+3*h*(rho0+p0)),
            'comoving_energy': float(a**3*rho),
            'comoving_work_rate': float(3*h*a**3*pressure),
            'radiation_invariant': float(a**4*rho),
            'trace': float(-rho+3*pressure)}


def solution_stress(eta, xi):
    qjet = field_from_chi(de_sitter_jet(eta), exact_solution(eta, 2., xi))
    return averaged_stress(eta, qjet, xi)


def test_profile(eta):
    return (0.2*eta*eta+0.1*eta-0.3, 0.4*eta+0.1, 0.4)


def compact_variation(eta):
    # v and v' vanish at both endpoints, including curvature boundary terms.
    u, du = (eta+4)*(eta+1), 2*eta+5
    return (u*u/16, u*du/8, (du*du+2*u)/8)


def action(field_shift=0., lapse_shift=0., xi=1/6, mass=0.3, k=2.):
    def density(eta):
        a, da, dda = de_sitter_jet(eta)
        v, dv, ddv = compact_variation(eta)
        q, dq, _ = test_profile(eta)
        q, dq = q+field_shift*v, dq+field_shift*dv
        n, dn = a+lapse_shift*v, da+lapse_shift*dv
        if n <= 0:
            raise ValueError('Positive lapse required.')
        scalar = 6*(dda/(a*n*n)+da*da/(a*a*n*n)-da*dn/(a*n**3))
        return -n*a**3*(-dq*dq/n**2+(k*k/a**2+mass*mass+xi*scalar)*q*q)/4
    return quadrature(density)


def variation_errors(xi=1/6):
    step = 1e-5
    numerical_field = (action(field_shift=step, xi=xi)-action(field_shift=-step, xi=xi))/(2*step)
    numerical_lapse = (action(lapse_shift=step, xi=xi)-action(lapse_shift=-step, xi=xi))/(2*step)

    def field_density(eta):
        a, da, dda = de_sitter_jet(eta)
        q, dq, ddq = test_profile(eta)
        residual = ddq+2*da/a*dq+(4+a*a*(0.3**2+xi*6*dda/a**3))*q
        return -a*a*residual*compact_variation(eta)[0]/2

    field_exact = quadrature(field_density)
    lapse_exact = quadrature(lambda t: -de_sitter_jet(t)[0]**3*
                            averaged_stress(t, test_profile(t), xi, mass=0.3)['rho']*
                            compact_variation(t)[0])
    return {'field_variation_error': abs(numerical_field-field_exact),
            'lapse_variation_error': abs(numerical_lapse-lapse_exact)}


def report():
    branches = []
    grid = np.linspace(-4, -1, 33)
    for xi in (0., 1/6):
        values = [solution_stress(t, xi) for t in grid]
        work = quadrature(lambda t: solution_stress(t, xi)['comoving_work_rate'])
        delta_energy = values[-1]['comoving_energy']-values[0]['comoving_energy']
        branches.append({'xi': xi, 'initial_rho': values[0]['rho'], 'final_rho': values[-1]['rho'],
                         'initial_comoving_energy': values[0]['comoving_energy'],
                         'final_comoving_energy': values[-1]['comoving_energy'],
                         'integrated_pressure_work': work,
                         'integrated_balance_error': abs(delta_energy+work),
                         'max_local_balance_error': max(abs(v['balance']) for v in values),
                         'max_omitted_improvement_residual': max(abs(v['minimal_only_balance']) for v in values),
                         'max_trace_magnitude': max(abs(v['trace']) for v in values),
                         'radiation_invariant_range': [min(v['radiation_invariant'] for v in values),
                                                       max(v['radiation_invariant'] for v in values)]})
    return {'round': 303,
            'scope': 'Classical test scalar, given de Sitter background, Hilbert stress including nonminimal improvement; local covariant energy accounting, not a selected gravitational dynamics or global universe energy.',
            'inputs': {'spatial_average': 'one period of Re(q exp(ikx))', 'H': 0.5, 'k': 2.,
                       'xi_branches': [0., 1/6], 'backreaction_included': False},
            'action_variation_checks': [{'xi': xi, **variation_errors(xi)} for xi in (0., 1/6, 0.5)],
            'branches': branches,
            'conclusion': 'Full stress conserves for either coupling; comoving matter energy changes by pressure work. Conservation does not select xi or Einstein dynamics.'}


class EnergyTests(unittest.TestCase):
    def test_01_field_variation_matches_wave_equation(self):
        for xi in (0., 1/6, 0.5):
            self.assertLess(variation_errors(xi)['field_variation_error'], 1e-8)

    def test_02_lapse_variation_matches_full_energy_density(self):
        for xi in (0., 1/6, 0.5):
            self.assertLess(variation_errors(xi)['lapse_variation_error'], 1e-8)

    def test_03_off_shell_noether_identity(self):
        for xi in (0., 1/6, 0.5, -0.2):
            for eta in (-3.7, -2.1, -1.2):
                qjet = (eta*eta+1j*eta, 2*eta+1j, 2.)
                stress = averaged_stress(eta, qjet, xi, mass=0.3)
                self.assertAlmostEqual(stress['balance'], stress['off_shell_source'], places=9)

    def test_04_on_shell_conservation_both_couplings(self):
        for xi in (0., 1/6):
            for eta in np.linspace(-4, -1, 19):
                self.assertLess(abs(solution_stress(eta, xi)['balance']), 1e-12)

    def test_05_omitting_curvature_stress_breaks_balance(self):
        value = solution_stress(-4., 1/6)
        self.assertGreater(abs(value['minimal_only_balance']), 0.01)
        self.assertLess(abs(value['balance']), 1e-12)

    def test_06_conformal_trace_and_redshift(self):
        for eta in np.linspace(-4, -1, 19):
            stress = solution_stress(eta, 1/6)
            self.assertLess(abs(stress['trace']), 1e-12)
            self.assertAlmostEqual(stress['radiation_invariant'], 0.5)

    def test_07_integrated_pressure_work(self):
        for row in report()['branches']:
            self.assertLess(row['integrated_balance_error'], 1e-11)
            self.assertGreater(abs(row['final_comoving_energy']-row['initial_comoving_energy']), 0.5)

    def test_08_spatial_average_normalization(self):
        angle = np.arange(128)*2*np.pi/128
        q, dq = 0.7+0.3j, -0.2+0.8j
        phi = np.real(q*np.exp(1j*angle))
        derivative = np.real(dq*np.exp(1j*angle))
        self.assertAlmostEqual(float(np.mean(phi*phi)), abs(q)**2/2)
        self.assertAlmostEqual(float(np.mean(derivative**2)), abs(dq)**2/2)
        self.assertAlmostEqual(float(np.mean(phi*derivative)), np.real(q.conjugate()*dq)/2)

    def test_09_zero_field_and_quadratic_resource_scaling(self):
        zero = averaged_stress(-2., (0j, 0j, 0j), 1/6)
        self.assertTrue(all(v == 0 for v in zero.values()))
        jet = field_from_chi(de_sitter_jet(-2.), exact_solution(-2., 2., 1/6))
        normal = averaged_stress(-2., jet, 1/6)
        scaled = averaged_stress(-2., tuple(3*v for v in jet), 1/6)
        self.assertAlmostEqual(scaled['rho'], 9*normal['rho'])
        self.assertAlmostEqual(scaled['p'], 9*normal['p'])


if __name__ == '__main__':
    main(__name__, 'scalar_energy_balance_audit', report)
