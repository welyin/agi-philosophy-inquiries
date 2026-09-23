"""Round 302: scalar propagation on the specified round-301 curved background."""
import unittest
import numpy as np
from conformal_volume_geometry_audit import de_sitter_jet, ricci_from_connection
from growing_stream_audit import main


def box_mode(lapse, scale, field_jet, k):
    n, dn, _ = lapse
    a, da, _ = scale
    phi, dphi, ddphi = field_jet
    return -(ddphi+(3*da/a-dn/n)*dphi)/n**2-k*k*phi/a**2


def field_from_chi(scale, chi_jet):
    a, da, dda = scale
    chi, dchi, ddchi = chi_jet
    return (chi/a, dchi/a-chi*da/a**2,
            ddchi/a-2*dchi*da/a**2+chi*(2*da**2/a**3-dda/a**2))


def omega_squared(eta, k, xi, mass=0., hubble=0.5):
    a, _, dda = de_sitter_jet(eta, hubble)
    return k*k+a*a*mass*mass+(6*xi-1)*dda/a


def fundamental_mode(eta, k, xi):
    """Two exactly known massless branches; return one mode and its conjugate."""
    plane = np.exp(-1j*k*eta)/np.sqrt(2*k)
    if xi == 1/6:
        return np.array([plane, -1j*k*plane, -k*k*plane])
    if xi != 0:
        raise ValueError('Analytic basis implemented only for xi=0 or xi=1/6.')
    f, df, ddf = 1-1j/(k*eta), 1j/(k*eta**2), -2j/(k*eta**3)
    return np.array([f*plane, (df-1j*k*f)*plane,
                     (ddf-2j*k*df-k*k*f)*plane])


def initial_state(k=2., left=-4.):
    # Identical finite-time Cauchy data for both couplings, not a derived vacuum.
    plane = np.exp(-1j*k*left)/np.sqrt(2*k)
    return np.array([plane, -1j*k*plane])


def exact_solution(eta, k, xi, left=-4.):
    initial = fundamental_mode(left, k, xi)[:2]
    coefficients = np.linalg.solve(np.column_stack((initial, initial.conj())), initial_state(k, left))
    basis = fundamental_mode(eta, k, xi)
    return coefficients[0]*basis+coefficients[1]*basis.conj()


def integrate_mode(xi, steps, k=2., mass=0., left=-4., right=-1.):
    times = np.linspace(left, right, steps+1)
    step = (right-left)/steps
    values = np.empty((steps+1, 2), dtype=complex)
    values[0] = initial_state(k, left)

    def rhs(t, y):
        return np.array([y[1], -omega_squared(t, k, xi, mass)*y[0]])

    for i, t in enumerate(times[:-1]):
        y = values[i]
        d1 = rhs(t, y)
        d2 = rhs(t+step/2, y+step*d1/2)
        d3 = rhs(t+step/2, y+step*d2/2)
        d4 = rhs(t+step, y+step*d3)
        values[i+1] = y+step*(d1+2*d2+2*d3+d4)/6
    return times, values


def wronskian(values):
    return values[:, 0]*values[:, 1].conj()-values[:, 0].conj()*values[:, 1]


def conformal_identity_error():
    errors = []
    for eta in np.linspace(-4, -1, 13):
        # Off-shell polynomial test, so both sides are generally nonzero.
        psi = (eta**3+1j*eta, 3*eta**2+1j, 6*eta)
        ajet = de_sitter_jet(eta)
        phi = field_from_chi(ajet, psi)
        _, scalar = ricci_from_connection(ajet, ajet)
        lhs = box_mode(ajet, ajet, phi, 1.3)-scalar*phi[0]/6
        rhs = (-psi[2]-1.3**2*psi[0])/ajet[0]**3
        errors.append(abs(lhs-rhs))
    return float(max(errors))


def report():
    branches = []
    for xi in (0., 1/6):
        exact = exact_solution(-1., 2., xi)
        errors, drifts = [], []
        for steps in (64, 128, 256):
            _, values = integrate_mode(xi, steps)
            errors.append(float(np.max(abs(values[-1]-exact[:2]))))
            drifts.append(float(np.max(abs(wronskian(values)-1j))))
        branches.append({'xi': xi, 'end_phi_abs': float(abs(exact[0])/2),
                         'end_chi': [float(exact[0].real), float(exact[0].imag)],
                         'rk4_endpoint_errors': errors,
                         'error_ratios': [errors[i]/errors[i+1] for i in (0, 1)],
                         'max_wronskian_drifts': drifts})
    return {'round': 302,
            'scope': 'Classical linear test-field equations on a supplied 3+1 de Sitter metric; no quantum vacuum selection, backreaction, discrete convergence or gravity derivation.',
            'inputs': {'H': 0.5, 'k': 2., 'm': 0., 'eta_interval': [-4, -1],
                       'rk4_steps': [64, 128, 256], 'same_finite_time_Cauchy_data': True},
            'conformal_identity_max_error': conformal_identity_error(),
            'branches': branches,
            'end_amplitude_ratio_minimal_to_conformal': branches[0]['end_phi_abs']/branches[1]['end_phi_abs'],
            'same_characteristic_cone': True,
            'BD_continuum': {'raw_xi': 0.5, 'conformal_xi': 1/6,
                             'constant_field_B1': -1.5, 'conformal_plane_residual_B_over_phi': -1.,
                             'R_coefficient_added_for_minimal': 0.5,
                             'R_coefficient_added_for_conformal': 1/3,
                             'finite_causal_set_operator_simulated': False}}


class PropagationTests(unittest.TestCase):
    def test_01_off_shell_conformal_identity(self):
        self.assertLess(conformal_identity_error(), 5e-12)

    def test_02_analytic_bases_solve_both_mode_equations(self):
        for xi in (0., 1/6):
            for eta in (-4., -2., -1.):
                mode = fundamental_mode(eta, 2., xi)
                self.assertLess(abs(mode[2]+omega_squared(eta, 2., xi)*mode[0]), 1e-14)

    def test_03_identical_initial_data_distinct_propagation(self):
        for xi in (0., 1/6):
            np.testing.assert_allclose(exact_solution(-4., 2., xi)[:2], initial_state())
        self.assertGreater(abs(exact_solution(-1., 2., 0.)[0]-exact_solution(-1., 2., 1/6)[0]), 0.1)

    def test_04_fourth_order_convergence_against_exact_solutions(self):
        for branch in report()['branches']:
            for ratio in branch['error_ratios']:
                self.assertTrue(14. < ratio < 18., ratio)
            self.assertLess(branch['rk4_endpoint_errors'][-1], 2e-8)

    def test_05_wronskian_conservation_and_field_normalization(self):
        for xi in (0., 1/6):
            times, values = integrate_mode(xi, 256)
            self.assertLess(np.max(abs(wronskian(values)-1j)), 1e-8)
            for eta in times[::32]:
                chi = exact_solution(eta, 2., xi)
                ajet = de_sitter_jet(eta)
                phi = field_from_chi(ajet, chi)
                conserved = ajet[0]**2*(phi[0]*phi[1].conjugate()-phi[0].conjugate()*phi[1])
                self.assertLess(abs(conserved-1j), 1e-14)

    def test_06_rescaled_and_covariant_equations_agree_off_shell(self):
        for xi in (0., 1/6, 0.5):
            for eta in (-3., -1.2):
                ajet = de_sitter_jet(eta)
                chi = np.array([eta*eta+1j, 2*eta, 2])
                phi = field_from_chi(ajet, chi)
                _, scalar = ricci_from_connection(ajet, ajet)
                lhs = box_mode(ajet, ajet, phi, 2.)-(0.3**2+xi*scalar)*phi[0]
                rhs = -(chi[2]+omega_squared(eta, 2., xi, 0.3)*chi[0])/ajet[0]**3
                self.assertLess(abs(lhs-rhs), 1e-13)

    def test_07_BD_is_not_the_conformal_scalar_operator(self):
        for eta in (-3., -1.2):
            ajet = de_sitter_jet(eta)
            phi = field_from_chi(ajet, fundamental_mode(eta, 2., 1/6))
            box = box_mode(ajet, ajet, phi, 2.)
            raw = box-1.5*phi[0]
            self.assertLess(abs(raw/phi[0]+1), 1e-13)
            self.assertLess(abs(raw+phi[0]), 1e-13)
            constant = box_mode(ajet, ajet, (1., 0., 0.), 0.)-1.5
            self.assertEqual(constant, -1.5)

    def test_08_operator_coordinate_invariance(self):
        for eta in (-3., -1.2):
            ajet = de_sitter_jet(eta)
            a = ajet[0]
            jet = (eta**2+1j*eta, 2*eta+1j, 2.)
            zjet = (jet[0], -eta*jet[1], eta**2*jet[2]+eta*jet[1])
            original = box_mode(ajet, ajet, jet, 1.7)
            changed = box_mode((2., 0., 0.), (a, a, a), zjet, 1.7)
            self.assertLess(abs(original-changed), 1e-13)

    def test_09_flat_limit_and_mass_term(self):
        for xi in (0., 1/6, 0.5):
            for eta in (-3., -1.):
                difference = omega_squared(eta, 2., xi, 0.3)-omega_squared(eta, 2., xi)
                a = de_sitter_jet(eta)[0]
                self.assertAlmostEqual(difference, a*a*0.3**2)
        field = (2+1j, -1+3j, -2j)
        self.assertEqual(box_mode((1, 0, 0), (1, 0, 0), field, 2.), -field[2]-4*field[0])


if __name__ == '__main__':
    main(__name__, 'curved_scalar_propagation_audit', report)
