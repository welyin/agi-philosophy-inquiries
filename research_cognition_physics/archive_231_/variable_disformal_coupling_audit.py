"""Round 340: functional, non-exponential disformal attraction in a slope band."""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main
from conformal_probe_motion_audit import integrate
from disformal_occupied_attractor_audit import initial_charge
from disformal_slope_robustness_audit import phase_rhs as constant_slope_rhs
from disformal_characteristics_audit import analytic_squares, matrices, spectrum

ROOT6 = np.sqrt(6.)


def coupling_data(phi, amplitude=.01, frequency=1., phase=0.):
    """Return log D, beta=d log D/d phi, and beta_phi, without exponent overflow."""
    argument = frequency*phi+phase
    logd = (np.log(.6)+ROOT6*phi
            +amplitude/(ROOT6*frequency)*(np.sin(argument)-np.sin(phase)))
    beta = ROOT6+amplitude/ROOT6*np.cos(argument)
    beta_phi = -amplitude*frequency/ROOT6*np.sin(argument)
    return logd, beta, beta_phi


def phase_rhs(y, amplitude=.01, frequency=1., phase=0.):
    """Autonomous (f,q,phi) flow with log(H/H0) as one passive quadrature."""
    fraction, q, phi, logh = y
    if not (0 <= fraction < 1 and 0 < q < 1):
        raise ValueError('Outside the positive homogeneous branch.')
    _, beta, _ = coupling_data(phi, amplitude, frequency, phase)
    old = constant_slope_rhs([fraction, q, logh], beta)
    return np.array([old[0], old[1], ROOT6*np.sqrt(1-fraction), old[2]])


@lru_cache(maxsize=20)
def solve_phase(fraction=.1, amplitude=.01, frequency=1., phase=0.,
                end=100., steps=8192):
    return integrate(lambda y: phase_rhs(y, amplitude, frequency, phase),
                     np.array([fraction, .4, 0., 0.]), end, steps)


def canonical_quantities(y, charge, amplitude=.01, frequency=1., phase=0.):
    n, phi, momentum, hubble, time, work = y
    a = np.exp(n)
    logd, beta, _ = coupling_data(phi, amplitude, frequency, phase)
    coupling = np.exp(logd)
    low = 0.
    high = min(momentum/a**3, 1/np.sqrt(coupling))
    for _ in range(60):
        v = (low+high)/2
        c = 1-coupling*v*v
        recovered = a**3*v+charge*coupling*v/(2*a**3*np.sqrt(c))
        if recovered < momentum:
            low = v
        else:
            high = v
    v = (low+high)/2
    q = coupling*v*v
    c = 1-q
    matter = charge/(2*a**6*np.sqrt(c))
    rho = .5*v*v+matter
    pressure = .5*v*v+c*matter
    return v, q, matter, rho, pressure, beta


def canonical_rhs(y, charge, amplitude=.01, frequency=1., phase=0.):
    n, phi, momentum, hubble, time, work = y
    v, q, matter, rho, pressure, beta = canonical_quantities(
        y, charge, amplitude, frequency, phase)
    return np.array([1., v/hubble,
                     np.exp(3*n)*matter*beta*q/(2*hubble),
                     -(rho+pressure)/(2*hubble), 1/hubble,
                     3*np.exp(3*n)*pressure])


@lru_cache(maxsize=8)
def solve_canonical(fraction=.1, amplitude=.01, frequency=1., phase=0.,
                    end=1., steps=2048):
    initial, charge = initial_charge(fraction)
    path = integrate(
        lambda y: canonical_rhs(y, charge, amplitude, frequency, phase),
        initial, end, steps)
    return path, charge


def routhian(a, phi, v, charge, amplitude=.01, frequency=1., phase=0.):
    coupling = np.exp(coupling_data(phi, amplitude, frequency, phase)[0])
    return .5*a**3*v*v-charge*np.sqrt(1-coupling*v*v)/(2*a**3)


@lru_cache(maxsize=4)
def independent_audit(frequency=1.):
    path, charge = solve_canonical(frequency=frequency)
    values = np.array([canonical_quantities(y, charge, frequency=frequency)
                       for y in path])
    reduced = solve_phase(frequency=frequency, end=1., steps=2048)
    return {
        'maximum_fraction_difference': float(np.max(abs(
            values[:, 2]/values[:, 3]-reduced[:, 0]))),
        'maximum_q_difference': float(np.max(abs(values[:, 1]-reduced[:, 1]))),
        'maximum_phi_difference': float(np.max(abs(path[:, 1]-reduced[:, 2]))),
        'maximum_log_H_difference': float(np.max(abs(
            np.log(path[:, 3]/path[0, 3])-reduced[:, 3]))),
        'maximum_friedmann_residual': float(np.max(abs(
            3*path[:, 3]**2-values[:, 3]))),
        'maximum_pressure_work_residual': float(np.max(abs(
            np.exp(3*path[:, 0])*values[:, 3]+path[:, 5]-values[0, 3]))),
    }


def trajectory_row(fraction, amplitude, frequency):
    path = solve_phase(fraction, amplitude, frequency)
    speeds = analytic_squares(*path[-1, :2])
    beta = coupling_data(path[:, 2], amplitude, frequency)[1]
    return {'initial_fraction': fraction,
            'slope_band_amplitude': amplitude, 'field_frequency': frequency,
            'end_log_scale_factor': 100.,
            'end_fraction': float(path[-1, 0]), 'end_q': float(path[-1, 1]),
            'end_phi': float(path[-1, 2]),
            'end_phi_mode_speed_squared': float(speeds[0]),
            'end_theta_mode_speed_squared': float(speeds[1]),
            'maximum_sampled_slope_offset': float(np.max(abs(ROOT6*beta-6))),
            'fraction_lower_bound': float(fraction*np.exp(-1.4))}


def report():
    witness = []
    for frequency in (1., 1000., 1e6):
        logd, beta, beta_phi = coupling_data(0., frequency=frequency,
                                            phase=np.pi/2)
        witness.append({'field_frequency': frequency,
                        'D_at_phi_zero': float(np.exp(logd)),
                        'beta_at_phi_zero': float(beta),
                        'beta_phi_at_phi_zero': float(beta_phi),
                        'global_slope_offset_bound': .01})
    fine = solve_phase(frequency=7.)
    coarse = solve_phase(frequency=7., steps=4096)
    return {
        'round': 340,
        'scope': 'Classical Einstein plus a canonical phi and a disformal massless chi, zero potentials, expanding spatially flat homogeneous backgrounds with positive phi velocity. For positive C2 coupling D(phi) on phi>=phi0, a global logarithmic-slope band extends the round-338 bootstrap to a function class. No uniform bound on beta_phi is required for homogeneous attraction or principal hyperbolicity; lower-derivative perturbations, actual EFT frequency windows, nonlinear stability and derivation from cognition remain unproved.',
        'literature_interface': {
            'source': 'https://arxiv.org/html/1510.01650',
            'location': 'General coupling action and section 3.1, equations 43-46; only the principal part is imported, with conformal factor one and zero potentials.'},
        'analytic_class': {
            'initial_fraction_interval': [.09, .1],
            'initial_q_interval': '(0,0.4]',
            'global_beta_sqrt6_minus6_bound': .01,
            'temporary_fraction_lower': .02,
            'proved_fraction_lower': float(.09*np.exp(-1.4)),
            'q_bound': 'q(N) <= q0 exp(-0.0225 N)',
            'phi_prime_interval': [float(np.sqrt(5.4)), ROOT6],
            'uniform_beta_phi_bound_required': False,
            'low_frequency_stability_claimed': False},
        'nonexponential_trajectories': [
            trajectory_row(.09, .01, 1.),
            trajectory_row(.1, -.01, 1.),
            trajectory_row(.1, .01, 7.)],
        'independent_routh_charge_einstein_checks': {
            'frequency_one': independent_audit(1.),
            'frequency_seven': independent_audit(7.)},
        'phase_step_refinement_error_f_q': float(np.max(abs(
            fine[::2, :2]-coarse[:, :2]))),
        'zero_occupation_nonattracting_control': {
            'amplitude': .01, 'frequency': 1.,
            'exact_solution': 'f=0; phi=sqrt(6) N; q=0.4 exp[0.01 sin(sqrt(6) N)/sqrt(6)]',
            'all_future_q_bounds': [float(.4*np.exp(-.01/ROOT6)),
                                   float(.4*np.exp(.01/ROOT6))]},
        'unbounded_curvature_family_at_identical_local_state': witness,
        'next_interface': 'Separate global functional selection from finite-window evidence; examine lower-derivative perturbations or an EFT validity window without treating a slope-band theorem as universal physical coupling.'}


class Checks(unittest.TestCase):
    def test_01_coupling_derivatives_come_from_nonexponential_log_function(self):
        for frequency in (1., 7.):
            for phi in (0., .37, 1.2):
                h = 1e-5/frequency
                center = coupling_data(phi, frequency=frequency)
                left = coupling_data(phi-h, frequency=frequency)
                right = coupling_data(phi+h, frequency=frequency)
                self.assertLess(abs((right[0]-left[0])/(2*h)-center[1]), 4e-10)
                self.assertLess(abs((right[1]-left[1])/(2*h)-center[2]), 3e-10)

    def test_02_zero_modulation_recovers_known_exponential_flow(self):
        for f, q, phi in ((.09, .4, 0.), (.05, .03, 2.), (.1, .1, 19.)):
            got = phase_rhs([f, q, phi, 0.], amplitude=0.)
            expected = constant_slope_rhs([f, q, 0.], ROOT6)
            np.testing.assert_allclose(got[[0, 1, 3]], expected, atol=1e-15)
            self.assertAlmostEqual(got[2], np.sqrt(6*(1-f)), places=14)

    def test_03_variable_slope_bootstrap_is_pointwise_and_closes_first_exit(self):
        for phi in np.linspace(0., 2*np.pi, 19):
            for f in (.02, .03, .06, .1):
                for q in (.001, .1, .4):
                    fp, qp, _, _ = phase_rhs([f, q, phi, 0.])
                    self.assertLessEqual(qp/q, -9*f/8+1e-14)
                    self.assertLess(fp, 0.)
                    self.assertLessEqual((fp/f)/qp, 3.5+1e-12)
        self.assertGreater(.09*np.exp(-1.4), .02)

    def test_04_nonexponential_future_paths_obey_uniform_analytic_envelope(self):
        for f, amp, freq in ((.09, .01, 1.), (.1, -.01, 1.), (.1, .01, 7.)):
            path = solve_phase(f, amp, freq)
            n = np.linspace(0., 100., len(path))
            self.assertTrue(np.all(path[:, 0] >= f*np.exp(-1.4)-1e-12))
            self.assertTrue(np.all(path[:, 1] <= .4*np.exp(-.0225*n)+1e-12))
            self.assertTrue(np.all(np.diff(path[:, 0]) <= 1e-14))
            self.assertTrue(np.all(np.diff(path[:, 1]) <= 1e-14))
            self.assertTrue(np.all(path[:, 2] >= np.sqrt(5.4)*n-1e-10))
            self.assertTrue(np.all(path[:, 2] <= ROOT6*n+1e-10))

    def test_05_true_mixed_principal_modes_remain_stable_and_converge(self):
        path = solve_phase(frequency=7.)
        for f, q, _, _ in path[::64]:
            speeds = analytic_squares(f, q)
            np.testing.assert_allclose(spectrum(f, q), np.sort(speeds), atol=2e-14)
            self.assertTrue(np.all(speeds >= .6-1e-14))
            self.assertTrue(np.all(1-speeds <= q+1e-14))
            k, g = matrices(f, q)
            self.assertGreater(np.linalg.eigvalsh(k)[0], 0.)
            self.assertGreater(np.linalg.eigvalsh(g)[0], 0.)
        self.assertLess(max(1-analytic_squares(*path[-1, :2])), 1e-5)

    def test_06_nonconstant_routh_source_is_independently_varied(self):
        y, charge = initial_charge(.1)
        v = np.sqrt(2/3)
        h = 1e-5
        for phi in (0., .05):
            y = y.copy(); y[1] = phi
            numerical = (routhian(1., phi+h, v, charge, frequency=7.)
                         -routhian(1., phi-h, v, charge, frequency=7.))/(2*h)
            logd, beta, _ = coupling_data(phi, frequency=7.)
            q = np.exp(logd)*v*v
            matter = charge/(2*np.sqrt(1-q))
            expected = matter*beta*q/2
            self.assertLess(abs(numerical-expected), 3e-10)

    def test_07_independent_canonical_charge_and_einstein_flows_agree(self):
        for frequency in (1., 7.):
            result = independent_audit(frequency)
            for key in ('maximum_fraction_difference', 'maximum_q_difference',
                        'maximum_phi_difference', 'maximum_log_H_difference'):
                self.assertLess(result[key], 4e-9)
            self.assertLess(result['maximum_friedmann_residual'], 2e-10)
            self.assertLess(result['maximum_pressure_work_residual'], 2e-10)

    def test_08_refinement_resolves_modulation_without_changing_decay_conclusion(self):
        fine = solve_phase(frequency=7.)
        coarse = solve_phase(frequency=7., steps=4096)
        self.assertLess(np.max(abs(fine[::2, :2]-coarse[:, :2])), 2e-8)
        full, _ = solve_canonical(frequency=7.)
        half, _ = solve_canonical(frequency=7., steps=1024)
        self.assertLess(np.max(abs(full[::2, :4]-half[:, :4])), 2e-9)

    def test_09_bounded_slope_does_not_bound_coupling_curvature(self):
        derivatives = []
        for frequency in (1., 1000., 1e6):
            logd, beta, beta_phi = coupling_data(0., frequency=frequency,
                                                phase=np.pi/2)
            self.assertAlmostEqual(np.exp(logd), .6, places=14)
            self.assertAlmostEqual(beta, ROOT6, places=14)
            derivatives.append(abs(beta_phi))
            rhs = phase_rhs([.1, .4, 0., 0.], frequency=frequency,
                            phase=np.pi/2)
            np.testing.assert_allclose(
                rhs, phase_rhs([.1, .4, 0., 0.], amplitude=0.), atol=1e-14)
        self.assertAlmostEqual(derivatives[-1]/derivatives[0], 1e6, places=7)
        self.assertGreater(derivatives[-1], 4000.)

    def test_10_global_charge_identity_checks_long_phase_trajectory(self):
        path = solve_phase(frequency=7.)
        f, q, phi, logh = path.T
        n = np.linspace(0., 100., len(path))
        # P_chi^2 = 6 a^6 H^2 f sqrt(1-q) is exactly conserved.
        log_charge_ratio = (6*n+2*logh+np.log(f/f[0])
                            +.5*np.log((1-q)/(1-q[0])))
        # q = 6 D(phi) (1-f) H^2, normalized by its initial value.
        logd = coupling_data(phi, frequency=7.)[0]
        log_q_error = (np.log(q/q[0])-(logd-logd[0])
                       -np.log((1-f)/(1-f[0]))-2*logh)
        self.assertLess(np.max(abs(log_charge_ratio)), 2e-8)
        self.assertLess(np.max(abs(log_q_error)), 2e-7)

    def test_11_zero_occupation_has_exact_nonvanishing_oscillatory_q(self):
        path = solve_phase(fraction=0., end=5., steps=2048)
        n = np.linspace(0., 5., len(path))
        exact = .4*np.exp(.01/ROOT6*np.sin(ROOT6*n))
        np.testing.assert_allclose(path[:, 1], exact, atol=1e-13, rtol=1e-13)
        self.assertTrue(np.all(path[:, 0] == 0.))
        self.assertGreater(.4*np.exp(-.01/ROOT6), .39)


if __name__ == '__main__':
    main(__name__, 'variable_disformal_coupling_audit', report)
