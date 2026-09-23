"""Round 319: nonlinear Wald focusing in a supplied classical scalar-tensor theory."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main


def fields(lam, xi=0.8, amplitude=0.6, shear=0.03, kappa=1.):
    s, c = np.sin(np.pi*np.asarray(lam)), np.cos(np.pi*np.asarray(lam))
    phi = amplitude*s*s
    dphi = 2*np.pi*amplitude*s*c
    ddphi = 2*np.pi**2*amplitude*(c*c-s*s)
    beta = shear*s**4
    dbeta = 4*np.pi*shear*s**3*c
    ddbeta = 4*np.pi**2*shear*(3*s*s*c*c-s**4)
    f = 1-kappa*xi*phi*phi
    df = -2*kappa*xi*phi*dphi
    ddf = -2*kappa*xi*(dphi*dphi+phi*ddphi)
    return phi, dphi, ddphi, beta, dbeta, ddbeta, f, df, ddf


def rk4(rhs, terminal, steps):
    grid = np.linspace(1., 0., steps+1)
    values = np.empty((steps+1, len(terminal)))
    values[0] = terminal
    h = -1./steps
    for j, lam in enumerate(grid[:-1]):
        y = values[j]
        k1 = rhs(lam, y); k2 = rhs(lam+h/2, y+h*k1/2)
        k3 = rhs(lam+h/2, y+h*k2/2); k4 = rhs(lam+h, y+h*k3)
        values[j+1] = y+h*(k1+2*k2+2*k3+k4)/6
    return grid[::-1].copy(), values[::-1].copy()


@lru_cache(maxsize=24)
def solve(steps=512, xi=0.8, amplitude=0.6, shear=0.03, kappa=1., final_slope=0.):
    if steps < 8 or kappa <= 0 or max(0., kappa*xi)*amplitude**2 >= 1:
        raise ValueError('Require resolved positive-kappa patch with strictly positive F.')
    def rhs(lam, y):
        _, dp, _, _, db, _, f, _, ddf = fields(lam, xi, amplitude, shear, kappa)
        curvature = (kappa*dp*dp+ddf)/f
        return np.array([y[1], -(curvature/2+db*db)*y[0]])
    lam, y = rk4(rhs, [1., final_slope], steps)
    r, dr = y.T
    if np.any(r <= 0):
        raise ValueError('A transverse caustic invalidates this positive area patch.')
    phi, dp, ddp, beta, db, ddb, f, df, ddf = fields(lam, xi, amplitude, shear, kappa)
    ricci = (kappa*dp*dp+ddf)/f
    ddr = -(ricci/2+db*db)*r
    theta = 2*dr/r; z = df/f
    generalized = theta+z
    focusing = theta*theta/2+2*db*db+kappa*dp*dp/f+z*z
    return {'lambda': lam, 'r': r, 'dr': dr, 'ddr': ddr, 'phi': phi, 'dphi': dp,
            'beta': beta, 'dbeta': db, 'ddbeta': ddb, 'F': f, 'dF': df, 'ddF': ddf,
            'Rkk': ricci, 'area': r*r, 'area_rate': 2*r*dr, 'theta': theta,
            'Theta': generalized, 'W': f*r*r, 'W_rate': f*r*r*generalized,
            'focusing_density': focusing, 'kappa': kappa, 'xi': xi}


def solve_entropy(steps=512, xi=0.8, amplitude=0.6, shear=0.03, kappa=1.):
    def rhs(lam, y):
        _, dp, _, _, db, _, f, df, _ = fields(lam, xi, amplitude, shear, kappa)
        theta = y[0]-df/f
        return np.array([-theta*theta/2-2*db*db-kappa*dp*dp/f-(df/f)**2, y[0]*y[1]])
    return rk4(rhs, [0., 1.], steps)[1]


def derivative(values, h):
    return (values[:-4]-8*values[1:-3]+8*values[3:-1]-values[4:])/(12*h)


def simpson(values, h):
    if (len(values)-1) % 2:
        raise ValueError('Even number of subintervals required.')
    return h/3*(values[0]+values[-1]+4*np.sum(values[1:-1:2])+2*np.sum(values[2:-2:2]))


def report():
    data = solve()
    minimum = int(np.argmin(data['area_rate']))
    fine = solve(2048)
    convergence = []
    for n in (32, 64, 128):
        coarse = solve(n)
        convergence.append({'steps': n, 'initial_r_error': float(abs(coarse['r'][0]-fine['r'][0]))})
    return {'round': 319,
            'scope': 'Classical scalar-tensor action F(phi)R/(2*kappa)-1/2(grad phi)^2, kappa=8*pi*G>0, F=1-kappa*xi*phi^2>0. Exact plane-wave reduction; future equilibrium and a regular null patch are supplied. No quantum entropy or derivation of the field equations.',
            'parameters': {'xi': 0.8, 'amplitude': 0.6, 'shear': 0.03, 'kappa': 1.},
            'minimum_F': float(np.min(data['F'])), 'minimum_transverse_radius': float(np.min(data['r'])),
            'initial_W': float(data['W'][0]), 'final_W': float(data['W'][-1]),
            'minimum_W_rate': float(np.min(data['W_rate'])),
            'area_decrease_control': {'lambda': float(data['lambda'][minimum]),
                                      'area_rate': float(data['area_rate'][minimum]),
                                      'W_rate': float(data['W_rate'][minimum]),
                                      'F_rate': float(data['dF'][minimum])},
            'integrated_focusing_error': float(abs(data['Theta'][0]-data['Theta'][-1]-simpson(data['focusing_density'], 1/512))),
            'independent_entropy_ODE_error': float(np.max(abs(solve_entropy()[:, 1]-data['W']))),
            'RK4_convergence': convergence,
            'wrong_future_boundary_control': {'terminal_Theta': -0.4,
                                              'initial_W': float(solve(amplitude=0., shear=0., final_slope=-0.2)['W'][0]),
                                              'final_W': 1.},
            'next_interface': 'Conformal Einstein-frame area and affine generator must be checked together, including the induced scalar kinetic term.'}


class Checks(unittest.TestCase):
    def test_01_null_metric_equation_independent_difference(self):
        d = solve(1024); h = 1/1024
        second = (-d['r'][4:]+16*d['r'][3:-1]-30*d['r'][2:-2]+16*d['r'][1:-3]-d['r'][:-4])/(12*h*h)
        geometric = -2*(second/d['r'][2:-2]+d['dbeta'][2:-2]**2)
        residual = d['F'][2:-2]*geometric-d['dphi'][2:-2]**2-d['ddF'][2:-2]
        self.assertLess(float(np.max(abs(residual))), 2e-7)

    def test_02_generalized_focusing_by_differentiation(self):
        d = solve(1024)
        residual = derivative(d['Theta'], 1/1024)+d['focusing_density'][2:-2]
        self.assertLess(float(np.max(abs(residual))), 2e-7)

    def test_03_Wald_increases_while_area_can_decrease(self):
        d = solve()
        self.assertGreater(np.min(d['F']), 0.7)
        self.assertGreaterEqual(np.min(d['W_rate']), -1e-10)
        self.assertGreaterEqual(np.min(np.diff(d['W'])), -1e-10)
        self.assertLess(np.min(d['area_rate']), -0.05)

    def test_04_RK4_convergence(self):
        errors = [r['initial_r_error'] for r in report()['RK4_convergence']]
        self.assertTrue(all(12 < a/b < 20 for a, b in zip(errors, errors[1:])))

    def test_05_independent_generalized_expansion_ODE(self):
        d = solve(); alternate = solve_entropy()
        self.assertLess(float(np.max(abs(alternate[:, 1]-d['W']))), 2e-8)
        self.assertLess(float(np.max(abs(alternate[:, 0]-d['Theta']))), 2e-8)

    def test_06_integrated_focusing_with_terminal_condition(self):
        self.assertLess(report()['integrated_focusing_error'], 2e-8)
        self.assertLess(abs(solve()['Theta'][-1]), 1e-14)

    def test_07_minimal_coupling_reduces_to_area(self):
        d = solve(xi=0.)
        np.testing.assert_allclose(d['W'], d['area'], atol=1e-13)
        self.assertGreaterEqual(np.min(d['area_rate']), -1e-12)
        np.testing.assert_allclose(d['Rkk'], d['dphi']**2, atol=1e-13)

    def test_08_zero_field_zero_shear_is_flat(self):
        d = solve(amplitude=0., shear=0.)
        np.testing.assert_allclose(d['r'], 1., atol=1e-14)
        np.testing.assert_allclose(d['Theta'], 0., atol=1e-14)

    def test_09_positive_F_alone_does_not_supply_future_equilibrium(self):
        d = solve(amplitude=0., shear=0., final_slope=-0.2)
        self.assertGreater(d['W'][0], d['W'][-1])
        self.assertLess(np.max(d['Theta']), 0.)
        np.testing.assert_allclose(d['Rkk'], 0.)

    def test_10_zero_or_negative_F_is_not_silently_regularized(self):
        for amplitude in (1., 1.1):
            with self.assertRaises(ValueError):
                solve(xi=1., amplitude=amplitude)


if __name__ == '__main__':
    main(__name__, 'wald_focusing_audit', report)
