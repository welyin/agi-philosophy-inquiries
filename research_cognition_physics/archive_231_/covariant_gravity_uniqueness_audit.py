"""Round 344: common metric and flat-space health do not uniquely select GR.

The specified counter-family is metric f(R)=R+alpha*R**2, alpha>0, with
the same universally minimally coupled matter action in the Jordan metric.
Only classical local/linear health is claimed. This is not a model of the
entire cognitive operation contract and is not a nonlinear plane-wave solution.
"""
import unittest
import numpy as np
from growing_stream_audit import main

ETA = np.diag([-1., 1., 1., 1.])
PAIRS = [(i, j) for i in range(4) for j in range(i, 4)]


def invariant_square(k):
    """k has lower indices; exp(i k_mu x^mu) is the Fourier convention."""
    return float(np.asarray(k) @ ETA @ np.asarray(k))


def linear_ricci(h, k):
    """Direct second-derivative contraction of a Fourier metric amplitude."""
    k = np.asarray(k)
    trace_h = np.trace(ETA @ h)
    contraction = h @ ETA @ k
    ricci = .5*(-np.outer(k, contraction)-np.outer(contraction, k)
                + invariant_square(k)*h + trace_h*np.outer(k, k))
    return ricci, np.trace(ETA @ ricci)


def curvature_from_connection(h, k):
    """Independent connection -> Riemann -> Ricci path, retaining all indices."""
    connection = np.zeros((4, 4, 4), dtype=complex)
    for r in range(4):
        for mu in range(4):
            for nu in range(4):
                connection[r, mu, nu] = .5j*sum(
                    ETA[r, s]*(k[mu]*h[s, nu]+k[nu]*h[s, mu]-k[s]*h[mu, nu])
                    for s in range(4))
    riemann = np.empty((4, 4, 4, 4), dtype=complex)
    for r in range(4):
        for s in range(4):
            for mu in range(4):
                for nu in range(4):
                    riemann[r, s, mu, nu] = 1j*(
                        k[mu]*connection[r, nu, s]-k[nu]*connection[r, mu, s])
    ricci = np.einsum('rsrn->sn', riemann)
    lowered = np.einsum('ar,rsuv->asuv', ETA, riemann)
    return ricci, np.trace(ETA @ ricci), lowered


def field_operator(h, k, alpha):
    """Full linearized f(R) metric equation, not just its trace."""
    ricci, scalar = linear_ricci(h, k)
    einstein = ricci-.5*ETA*scalar
    correction = 2*alpha*(np.outer(k, k)-ETA*invariant_square(k))*scalar
    return einstein+correction


def gauge_amplitude(k, xi):
    # A common factor i may be absorbed in the complex xi amplitude.
    return np.outer(k, xi)+np.outer(xi, k)


def scalar_wave(alpha=.25, spatial=(0., 0., .8), ricci_amplitude=.002):
    mass_squared = 1/(6*alpha)
    spatial = np.asarray(spatial)
    omega = np.sqrt(spatial @ spatial+mass_squared)
    k = np.r_[-omega, spatial]
    h = -2*alpha*ricci_amplitude*ETA
    return h, k


def pack(h):
    return np.array([h[i, j] for i, j in PAIRS])


def symbol_matrix(k, alpha):
    columns = []
    for i, j in PAIRS:
        h = np.zeros((4, 4))
        h[i, j] = h[j, i] = 1.
        columns.append(pack(field_operator(h, k, alpha)))
    return np.column_stack(columns)


def nullity(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    return matrix.shape[1]-int(np.count_nonzero(singular > 1e-10*singular[0]))


def mode_count(k, alpha):
    gauge = np.column_stack([pack(gauge_amplitude(k, row)) for row in np.eye(4)])
    gauge_rank = int(np.linalg.matrix_rank(gauge, tol=1e-10))
    full_nullity = nullity(symbol_matrix(k, alpha))
    return {'fourier_kernel_dimension': full_nullity,
            'pure_gauge_dimension': gauge_rank,
            'physical_polarizations_on_this_shell': full_nullity-gauge_rank}


def exact_metric_equation_at_point(ricci, hessian_r, alpha):
    """Evaluate the exact covariant equation in a local orthonormal frame.

    hessian_r is covariant nabla_mu nabla_nu R. This algebraic trace check
    does not purport to generate arbitrary nonlinear geometric jets.
    """
    scalar = np.trace(ETA @ ricci)
    fprime = 1+2*alpha*scalar
    f = scalar+alpha*scalar*scalar
    return (fprime*ricci-.5*f*ETA
            +2*alpha*(ETA*np.trace(ETA @ hessian_r)-hessian_r))


def potential(varphi, alpha=.25, kappa=1.):
    b = np.sqrt(2*kappa/3)
    return np.expm1(-b*np.asarray(varphi))**2/(8*kappa*alpha)


def averaged_scalar_energy(alpha, spatial, ricci_amplitude, kappa=1.):
    """Positive quadratic canonical energy per volume, averaged over phase."""
    mass_squared = 1/(6*alpha)
    momentum_squared = np.asarray(spatial) @ np.asarray(spatial)
    amplitude = np.sqrt(6/kappa)*alpha*ricci_amplitude
    omega_squared = momentum_squared+mass_squared
    # Directly average 1/2(phi_dot^2 + |grad phi|^2 + m^2 phi^2).
    return .25*amplitude**2*(omega_squared+momentum_squared+mass_squared)


def witness(alpha=.25, spatial=(0., 0., .8), amplitude=.002):
    h, k = scalar_wave(alpha, spatial, amplitude)
    _, scalar, riemann = curvature_from_connection(h, k)
    tidal = np.real(riemann[0, 1:, 0, 1:])
    m2 = 1/(6*alpha)
    return {'alpha': alpha, 'ricci_amplitude': amplitude,
            'frequency': float(-k[0]), 'spatial_wavevector': list(spatial),
            'mass_squared': m2,
            'group_speed': float(np.linalg.norm(spatial)/(-k[0])),
            'reconstructed_R': float(np.real(scalar)),
            'f_R_min_over_phase_to_linear_order': 1-abs(2*alpha*amplitude),
            'f_R_max_over_phase_to_linear_order': 1+abs(2*alpha*amplitude),
            'f_R_equation_max_residual': float(np.max(abs(field_operator(h, k, alpha)))),
            'Einstein_vacuum_equation_max_residual': float(np.max(abs(field_operator(h, k, 0.)))),
            'tidal_amplitude_matrix': tidal.tolist(),
            'longitudinal_over_transverse_tidal_amplitude': float(m2/k[0]**2),
            'canonical_scalar_energy_per_volume_phase_average': float(
                averaged_scalar_energy(alpha, spatial, amplitude)),
            'f_R_modes_at_scalar_shell': mode_count(k, alpha),
            'Einstein_modes_at_same_timelike_shell': mode_count(k, 0.)}


def report():
    alpha = .25
    null_k = np.array([-1., 0., 0., 1.])
    off_shell_k = np.array([-1.7, .2, .1, .8])
    rng = np.random.default_rng(344)
    max_path_difference = 0.
    max_divergence = 0.
    for _ in range(24):
        h = rng.normal(size=(4, 4)); h = (h+h.T)/2
        k = rng.normal(size=4)
        ricci, scalar = linear_ricci(h, k)
        other, other_scalar, _ = curvature_from_connection(h, k)
        max_path_difference = max(max_path_difference, np.max(abs(ricci-other)),
                                  abs(scalar-other_scalar))
        max_divergence = max(max_divergence,
                             np.max(abs((ETA @ k) @ field_operator(h, k, alpha))))
    return {'round': 344, 'baseline': 'Frozen research through round 341.',
            'scope': 'Specified smooth four-dimensional metric R+alpha R^2 theory near flat vacuum, with identical minimally coupled matter. Tests exact geometric uniqueness, not the full cognitive contract or low-energy GR universality.',
            'hypothesis': 'A common matter metric, covariance and positive linear energy near flat space may suffice to force the Einstein equation.',
            'result': 'This sufficiency claim is false for the declared geometric conditions: metric R+alpha R^2 with alpha>0 retains one healthy extra gravitational scalar.',
            'added_inputs': ['Smooth 3+1 Lorentz geometry, Levi-Civita metric connection and a local classical action.',
                             'Identical specified matter action S_m[g,psi], minimally and universally coupled to the physical Jordan metric in both comparisons.',
                             'alpha>0, kappa>0, F=1+2alpha R>0 near the flat vacuum; linearized health only.'],
            'lovelock_interface': {'primary_proof': 'https://arxiv.org/html/1005.2386',
                                  'precise_input': 'Second-order natural, symmetric, identically divergence-free two-tensor constructed solely from the metric.',
                                  'dimension': 4,
                                  'basis': ['metric', 'Einstein tensor'],
                                  'missing_selection': 'Second order in the pure-metric field equations is not supplied by covariance or linearized positivity.',
                                  'constants_not_fixed': ['Newton coupling', 'cosmological constant']},
            'linear_wave_interface': 'https://arxiv.org/html/1104.0819',
            'sign_convention': '(-+++), R^rho_sigma_mu_nu = d_mu Gamma^rho_nu_sigma - d_nu Gamma^rho_mu_sigma + ...',
            'wave_witness': witness(),
            'null_shell_mode_count': {'f_R': mode_count(null_k, alpha),
                                      'Einstein': mode_count(null_k, 0.)},
            'generic_off_shell_mode_count': {'f_R': mode_count(off_shell_k, alpha),
                                             'Einstein': mode_count(off_shell_k, 0.)},
            'independent_connection_vs_direct_Ricci_max_error': float(max_path_difference),
            'off_shell_linear_Bianchi_max_residual': float(max_divergence),
            'derivative_orders': {'Einstein_symbol_scaling_degree': 2,
                                  'R_squared_correction_symbol_scaling_degree': 4},
            'scalar_tensor_interface': {'F': '1+2alpha R',
                                        'canonical_field': 'sqrt(3/(2kappa))*log F',
                                        'potential': '(1-exp(-sqrt(2kappa/3)*varphi))^2/(8kappa alpha)',
                                        'mass_squared_at_flat_vacuum': 1/(6*alpha),
                                        'transformed_matter': 'S_m[F^(-1)g_E,psi]; the scalar cannot be discarded while claiming unchanged measurements.'},
            'limits': ['The displayed metric wave solves the complete linearized equations, not the full nonlinear equations at finite amplitude.',
                       'An infinite monochromatic wave is a Fourier diagnostic; quoted energy is a phase average per unit volume, not a finite total energy.',
                       'No global nonlinear stability, quantum completion or empirical fit is claimed.',
                       'No proof is supplied that this counter-family realizes the full cognitive F+U+C+P operation contract.',
                       'This falsifies only the listed geometric sufficiency claim, not the possibility of a future cognitive derivation of GR.'],
            'next_interface': 'Either derive the Lovelock pure-metric second-order restriction independently from the cognitive contract, or derive a controlled low-energy elimination of the scalar with state and error bounds.'}


class Checks(unittest.TestCase):
    def test_01_connection_and_direct_curvature_agree_for_general_fourier_data(self):
        rng = np.random.default_rng(1344)
        for _ in range(30):
            h = rng.normal(size=(4, 4)); h = (h+h.T)/2
            k = rng.normal(size=4)
            r, scalar = linear_ricci(h, k)
            other, s2, riemann = curvature_from_connection(h, k)
            np.testing.assert_allclose(r, other, atol=3e-14)
            self.assertAlmostEqual(scalar, s2.real, places=12)
            np.testing.assert_allclose(riemann, -riemann.swapaxes(0, 1), atol=3e-14)
            np.testing.assert_allclose(riemann, -riemann.swapaxes(2, 3), atol=3e-14)
            np.testing.assert_allclose(riemann, riemann.transpose(2, 3, 0, 1), atol=3e-14)

    def test_02_curvature_and_field_equations_remove_four_pure_gauges(self):
        rng = np.random.default_rng(2344)
        for _ in range(20):
            k = rng.normal(size=4)
            h = gauge_amplitude(k, rng.normal(size=4))
            r, scalar, riemann = curvature_from_connection(h, k)
            np.testing.assert_allclose(riemann, 0., atol=1e-13)
            np.testing.assert_allclose(field_operator(h, k, .25), 0., atol=1e-12)
            self.assertLess(abs(scalar), 1e-13)
            np.testing.assert_allclose(r, 0., atol=1e-13)

    def test_03_exact_nonlinear_trace_cancels_the_algebraic_R_squared_terms(self):
        rng = np.random.default_rng(3344)
        for alpha in (.04, .25, 2.):
            for _ in range(12):
                r = rng.normal(size=(4, 4)); r = (r+r.T)/2
                d2r = rng.normal(size=(4, 4)); d2r = (d2r+d2r.T)/2
                result = np.trace(ETA @ exact_metric_equation_at_point(r, d2r, alpha))
                expected = 6*alpha*np.trace(ETA @ d2r)-np.trace(ETA @ r)
                self.assertAlmostEqual(result, expected, places=12)

    def test_04_reconstructed_scalar_wave_solves_all_ten_field_equations(self):
        for alpha in (.04, .25, 2.):
            for spatial in ((0., 0., 0.), (0., 0., .8), (.2, -.4, .7)):
                h, k = scalar_wave(alpha, spatial)
                _, scalar = linear_ricci(h, k)
                self.assertAlmostEqual(scalar, .002, places=14)
                np.testing.assert_allclose(field_operator(h, k, alpha), 0., atol=2e-16)
                self.assertAlmostEqual(-invariant_square(k), 1/(6*alpha), places=13)

    def test_05_same_wave_is_not_an_Einstein_vacuum_wave(self):
        h, k = scalar_wave()
        self.assertGreater(np.max(abs(field_operator(h, k, 0.))), 1e-4)
        self.assertAlmostEqual(np.trace(ETA @ field_operator(h, k, 0.)), -.002, places=14)
        # Adding any infinitesimal coordinate change cannot remove this R.
        changed = h+.01*gauge_amplitude(k, np.array([.3, -.2, .4, .7]))
        self.assertAlmostEqual(linear_ricci(changed, k)[1], .002, places=14)

    def test_06_two_transverse_traceless_null_modes_survive_in_both_theories(self):
        k = np.array([-1., 0., 0., 1.])
        plus = np.diag([0., 1., -1., 0.])
        cross = np.zeros((4, 4)); cross[1, 2] = cross[2, 1] = 1.
        for h in (plus, cross):
            for alpha in (0., .04, .25, 2.):
                np.testing.assert_allclose(field_operator(h, k, alpha), 0., atol=1e-14)
            self.assertEqual(linear_ricci(h, k)[1], 0.)

    def test_07_symbol_kernel_distinguishes_physical_modes_from_gauge(self):
        _, scalar_k = scalar_wave()
        null_k = np.array([-1., 0., 0., 1.])
        off_k = np.array([-1.7, .2, .1, .8])
        self.assertEqual(mode_count(scalar_k, .25)['physical_polarizations_on_this_shell'], 1)
        self.assertEqual(mode_count(scalar_k, 0.)['physical_polarizations_on_this_shell'], 0)
        for alpha in (0., .25):
            self.assertEqual(mode_count(null_k, alpha)['physical_polarizations_on_this_shell'], 2)
            self.assertEqual(mode_count(off_k, alpha)['physical_polarizations_on_this_shell'], 0)

    def test_08_linear_Bianchi_identity_holds_off_shell(self):
        rng = np.random.default_rng(8344)
        for _ in range(24):
            h = rng.normal(size=(4, 4)); h = (h+h.T)/2
            k = rng.normal(size=4)
            for alpha in (0., .25, 2.):
                np.testing.assert_allclose((ETA @ k) @ field_operator(h, k, alpha),
                                           0., atol=2e-12)

    def test_09_geodesic_deviation_has_a_gauge_invariant_scalar_tidal_signal(self):
        alpha = .25; amplitude = .002
        for spatial in ((0., 0., .8), (.2, -.4, .7)):
            h, k = scalar_wave(alpha, spatial, amplitude)
            tidal = curvature_from_connection(h, k)[2][0, 1:, 0, 1:]
            expected = -alpha*amplitude*(k[0]**2*np.eye(3)-np.outer(spatial, spatial))
            np.testing.assert_allclose(tidal, expected, atol=1e-16)
            h = h+gauge_amplitude(k, np.array([.13, -.08, .11, -.23]))
            np.testing.assert_allclose(curvature_from_connection(h, k)[2][0, 1:, 0, 1:],
                                       expected, atol=1e-15)
        row = witness()
        self.assertGreater(row['longitudinal_over_transverse_tidal_amplitude'], 0.)
        self.assertLess(row['longitudinal_over_transverse_tidal_amplitude'], 1.)

    def test_10_auxiliary_scalar_elimination_recovers_the_original_action_density(self):
        for alpha in (.04, .25, 2.):
            for scalar in np.linspace(-.1, .1, 11):
                auxiliary = 1+2*alpha*scalar
                jordan = auxiliary*scalar-(auxiliary-1)**2/(4*alpha)
                self.assertAlmostEqual(jordan, scalar+alpha*scalar*scalar, places=14)
                self.assertAlmostEqual(scalar-(auxiliary-1)/(2*alpha), 0., places=14)

    def test_11_canonical_potential_mass_and_linear_energy_are_positive(self):
        for alpha in (.04, .25, 2.):
            eps = 1e-4
            numerical_mass_squared = (potential(eps, alpha)+potential(-eps, alpha))/eps**2
            self.assertAlmostEqual(numerical_mass_squared/(1/(6*alpha)), 1., places=7)
            self.assertTrue(np.all(potential(np.linspace(-2., 2., 41), alpha) >= 0.))
            for spatial in ((0., 0., 0.), (0., 0., .8), (.2, -.4, .7)):
                self.assertGreater(averaged_scalar_energy(alpha, spatial, .002), 0.)
                h, k = scalar_wave(alpha, spatial)
                self.assertGreater(-k[0], 0.)
        # The opposite alpha sign fails the same local stability criterion.
        self.assertLess(1/(6*(-.25)), 0.)

    def test_12_Einstein_frame_cancels_linear_metric_wave_but_not_the_matter_scalar(self):
        alpha = .25; amplitude = .002
        h, _ = scalar_wave(alpha, ricci_amplitude=amplitude)
        delta_f = 2*alpha*amplitude
        h_einstein = h+delta_f*ETA
        np.testing.assert_allclose(h_einstein, 0., atol=1e-16)
        varphi_amplitude = np.sqrt(3/2)*delta_f
        self.assertGreater(varphi_amplitude, 0.)
        # For the same massive matter point action, m_E = m_J / sqrt(F).
        for f in (.9, 1., 1.1):
            mass_e = 2./np.sqrt(f)
            self.assertAlmostEqual(mass_e*np.sqrt(f), 2., places=14)
        # d log(m_E)/d varphi = -1/sqrt(6), not zero.
        eps = 1e-5
        masses = 2*np.exp(-np.array([-eps, eps])/np.sqrt(6))
        derivative = (np.log(masses[1])-np.log(masses[0]))/(2*eps)
        self.assertAlmostEqual(derivative, -1/np.sqrt(6), places=10)

    def test_13_R_squared_metric_equation_has_fourth_not_second_derivative_order(self):
        h = np.diag([.3, .7, -.4, .2])
        k = np.array([-1.3, .2, -.1, .6])
        base = field_operator(h, k, 0.)
        correction = field_operator(h, k, .25)-base
        self.assertGreater(np.max(abs(correction)), .1)
        for scale in (.5, 2., 3.):
            scaled_base = field_operator(h, scale*k, 0.)
            scaled_correction = field_operator(h, scale*k, .25)-scaled_base
            np.testing.assert_allclose(scaled_base, scale**2*base, rtol=2e-14, atol=1e-13)
            np.testing.assert_allclose(scaled_correction, scale**4*correction,
                                       rtol=2e-14, atol=1e-13)


if __name__ == '__main__':
    main(__name__, 'covariant_gravity_uniqueness_audit', report)
