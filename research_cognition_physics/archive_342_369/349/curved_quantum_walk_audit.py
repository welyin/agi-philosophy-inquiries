"""Round 349: a static curved-Dirac interface with a local unitary walk.

Uses a specified optical coordinate, not the general four-component paired
construction. All numerical work uses NumPy; errors are state/H1 bounds.
"""
from functools import lru_cache
from math import factorial
import unittest

import numpy as np

from growing_stream_audit import main


V0 = .7
EPS = .3
MASS = .8
TIME = .5
X = np.array([[0., 1.], [1., 0.]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
INIT_NORM2 = 1+.2**2/2+.35**2


def background(y, epsilon=EPS):
    y = np.asarray(y)
    v = V0*(1+epsilon*np.cos(2*np.pi*y))
    vy = -V0*epsilon*2*np.pi*np.sin(2*np.pi*y)
    vyy = -V0*epsilon*(2*np.pi)**2*np.cos(2*np.pi*y)
    physical_x = V0*(y+epsilon*np.sin(2*np.pi*y)/(2*np.pi))
    return physical_x, v, vy, vyy


def initial(y):
    y = np.asarray(y)
    return np.stack((1+.2*np.sin(2*np.pi*y), .35*np.exp(2j*np.pi*y)),
                    axis=-1)/np.sqrt(INIT_NORM2)


def derivative_initial_norm():
    return float(2*np.pi*np.sqrt((.2**2/2+.35**2)/INIT_NORM2))


def initial_modes(cutoff):
    out = np.zeros((2*cutoff+1, 2), dtype=complex)
    out[cutoff, 0] = 1.
    out[cutoff+1, 0] = -.1j
    out[cutoff-1, 0] = .1j
    out[cutoff+1, 1] = .35
    return out/np.sqrt(INIT_NORM2)


def coin(chi, y, duration, mass=MASS, epsilon=EPS):
    v = background(y, epsilon)[1]
    theta = duration*mass*v
    return np.cos(theta)[:, None]*chi-1j*np.sin(theta)[:, None]*chi[:, ::-1]


def one_step(chi, mass=MASS, epsilon=EPS):
    count = len(chi)
    y = np.arange(count)/count
    rotated = coin(chi, y, 1/count, mass, epsilon)
    out = np.empty_like(rotated)
    out[:, 0] = np.roll(rotated[:, 0], 1)
    out[:, 1] = np.roll(rotated[:, 1], -1)
    return out


def inverse_step(chi, mass=MASS, epsilon=EPS):
    count = len(chi)
    shifted = np.empty_like(chi)
    shifted[:, 0] = np.roll(chi[:, 0], -1)
    shifted[:, 1] = np.roll(chi[:, 1], 1)
    return coin(shifted, np.arange(count)/count, -1/count, mass, epsilon)


@lru_cache(maxsize=None)
def solve_walk(count, time=TIME, mass=MASS, epsilon=EPS):
    steps = int(round(time*count))
    if abs(steps/count-time) > 1e-12:
        raise ValueError("Use an integer number of steps.")
    chi = initial(np.arange(count)/count)
    for _ in range(steps):
        chi = one_step(chi, mass, epsilon)
    return chi


def fourier_hamiltonian(cutoff, mass=MASS, epsilon=EPS):
    modes = np.arange(-cutoff, cutoff+1)
    kinetic = np.kron(np.diag(2*np.pi*modes), Z)
    vhat = np.diag(np.full(len(modes), V0))
    vhat += np.diag(np.full(len(modes)-1, V0*epsilon/2), 1)
    vhat += np.diag(np.full(len(modes)-1, V0*epsilon/2), -1)
    return kinetic+mass*np.kron(vhat, X)


@lru_cache(maxsize=None)
def continuum_coefficients(cutoff=18, time=TIME, mass=MASS, epsilon=EPS):
    h = fourier_hamiltonian(cutoff, mass, epsilon)
    energies, vecs = np.linalg.eigh(h)
    state = initial_modes(cutoff).reshape(-1)
    out = vecs @ (np.exp(-1j*time*energies)*(vecs.conj().T @ state))
    return out.reshape(-1, 2)


def evaluate_modes(coefficients, y):
    cutoff = (len(coefficients)-1)//2
    phases = np.exp(2j*np.pi*np.outer(y, np.arange(-cutoff, cutoff+1)))
    return phases @ coefficients


def derivative_samples(values):
    count = len(values)
    frequencies = 2*np.pi*np.fft.fftfreq(count, 1/count)
    return np.fft.ifft(1j*frequencies[:, None]*np.fft.fft(values, axis=0), axis=0)


def optical_generator(chi, mass=MASS, epsilon=EPS, omit_connection=False):
    count = len(chi)
    y = np.arange(count)/count
    _, v, vy, _ = background(y, epsilon)
    out = -derivative_samples(chi)*np.array([1., -1.])[None, :]
    out -= 1j*mass*v[:, None]*chi[:, ::-1]
    if omit_connection:
        out += .5*(vy/v)[:, None]*chi*np.array([1., -1.])[None, :]
    return out


def physical_generator(psi, y, mass=MASS, epsilon=EPS, omit_connection=False):
    _, v, vy, _ = background(y, epsilon)
    # d/dx = (1/v) d/dy, so v*d_x psi = d_y psi.
    out = -derivative_samples(psi)*np.array([1., -1.])[None, :]
    if not omit_connection:
        out -= .5*(vy/v)[:, None]*psi*np.array([1., -1.])[None, :]
    return out-1j*mass*v[:, None]*psi[:, ::-1]


def error_constants(time=TIME, mass=MASS, epsilon=EPS):
    vmax = V0*(1+abs(epsilon))
    vprime = 2*np.pi*V0*abs(epsilon)
    m1 = derivative_initial_norm()+time*abs(mass)*vprime
    c1 = time*abs(mass)*(vprime+2*vmax*m1)/2
    return {'derivative_bound_M1':m1, 'continuum_L2_coefficient_C1':c1,
            'sampled_L2_coefficient':float(np.sqrt(c1*c1+4*m1*c1)),
            'vmax':vmax, 'v_y_sup':vprime}


def galerkin_tail_bounds(cutoff, time=TIME, mass=MASS, epsilon=EPS):
    # Initial modes |k|<=1; a variable-potential factor shifts at most one mode.
    # Constant mass belongs to the free block-diagonal propagator.
    q = time*abs(mass)*V0*abs(epsilon)
    first = q**cutoff/factorial(cutoff)
    ratio = q/(cutoff+1)
    weighted_ratio = q*(cutoff+2)/(cutoff+1)**2
    if weighted_ratio >= 1:
        raise ValueError("Choose a larger cutoff for the simple geometric tail bound.")
    # Upper bounds on the entire infinite sums, not truncated numerical sums.
    l2 = 2*first/(1-ratio)
    derivative = 4*np.pi*(cutoff+1)*first/(1-weighted_ratio)
    return float(l2), float(derivative)


def ricci_from_connection(v, vx, vxx):
    gamma = np.zeros((2, 2, 2))
    deriv = np.zeros((2, 2, 2, 2))
    gamma[0, 0, 1] = gamma[0, 1, 0] = vx/v
    gamma[1, 0, 0] = v*vx
    deriv[1, 0, 0, 1] = deriv[1, 0, 1, 0] = vxx/v-vx*vx/v**2
    deriv[1, 1, 0, 0] = vx*vx+v*vxx
    ricci = np.zeros((2, 2))
    for mu in range(2):
        for nu in range(2):
            for lam in range(2):
                ricci[mu, nu] += deriv[lam, lam, mu, nu]-deriv[nu, lam, mu, lam]
                for sig in range(2):
                    ricci[mu, nu] += (gamma[lam, mu, nu]*gamma[sig, lam, sig]
                                      -gamma[sig, mu, lam]*gamma[lam, nu, sig])
    return float(-ricci[0, 0]/v**2+ricci[1, 1])


def comparison(count):
    y = np.arange(count)/count
    _, v, _, _ = background(y)
    chi = solve_walk(count)
    reference = evaluate_modes(continuum_coefficients(), y)
    psi = chi/np.sqrt(v)[:, None]
    e = float(np.linalg.norm(chi-reference)/np.sqrt(count))
    correct_norm = float(np.sum(v[:, None]*abs(psi)**2)/count)
    wrong_norm = float(np.sum(abs(psi)**2)/count)
    c = error_constants()
    return {'sites':count, 'optical_spacing_and_step':1/count,
            'steps':int(TIME*count), 'physical_circumference':V0,
            'state_L2_error':e, 'sampled_analytic_error_upper':c['sampled_L2_coefficient']/count,
            'weighted_physical_norm':correct_norm,
            'unweighted_physical_samples_sum':wrong_norm,
            'step_physical_distance_upper':c['vmax']/count}


def report():
    count = 1024
    y = np.arange(count)/count
    _, v, vy, vyy = background(y)
    chi = initial(y)
    right = optical_generator(chi)
    wrong = optical_generator(chi, omit_connection=True)
    right_rate = float(2*np.vdot(chi, right).real/count)
    wrong_rate = float(2*np.vdot(chi, wrong).real/count)
    direct_wrong = float(np.mean(vy/v*(abs(chi[:, 0])**2-abs(chi[:, 1])**2)))
    vx = vy/v
    vxx = (vyy*v-vy*vy)/v**3
    curvature = -2*vxx/v
    massless = solve_walk(128, mass=0.)
    expected = initial(np.arange(128)/128)
    expected[:, 0] = np.roll(expected[:, 0], 64)
    expected[:, 1] = np.roll(expected[:, 1], -64)
    return {
        'round':349,
        'scope':'A specified static positive-lapse 1+1 curved-Dirac subclass. A uniform optical-coordinate conditional shift and position-dependent mass coin are exactly local and unitary; a supplied nonuniform physical sampling map gives B=-v sigma_z and its half-derivative connection. State-specific H1 and finite-time convergence is proved. This is not the general four-component paired walk construction and not dynamical geometry or GR.',
        'background':{'v0':V0, 'epsilon':EPS, 'mass':MASS, 'time':TIME,
                      'x_of_y':'v0*(y+epsilon*sin(2*pi*y)/(2*pi))',
                      'v_of_y':'v0*(1+epsilon*cos(2*pi*y))',
                      'metric':'ds^2=-v(x)^2 dt^2+dx^2=v(y)^2*(-dt^2+dy^2)',
                      'scalar_curvature_min':float(np.min(curvature)),
                      'scalar_curvature_max':float(np.max(curvature)),
                      'curvature_convention':'R_mu_nu=partial_lam Gamma^lam_mu_nu-partial_nu Gamma^lam_mu_lam+Gamma^lam_mu_nu Gamma^sig_lam_sig-Gamma^sig_mu_lam Gamma^lam_nu_sig'},
        'exact_connections':{
            'field_map':'chi(y)=sqrt(v(x(y)))*psi(x(y)); L2(dy) and L2(dx) are unitarily related.',
            'physical_PDE':'partial_t psi=-v sigma_z partial_x psi-(v_x/2) sigma_z psi-i*m*v sigma_x psi.',
            'continuity':'partial_t(psi^dagger psi)+partial_x(v psi^dagger sigma_z psi)=0.',
            'optical_PDE':'partial_t chi=-sigma_z partial_y chi-i*m*v(y) sigma_x chi.'},
        'error_constants':error_constants(),
        'galerkin_reference':{'cutoff_mode':18, 'dimension':74,
                              'L2_tail_upper':galerkin_tail_bounds(18)[0],
                              'derivative_tail_upper':galerkin_tail_bounds(18)[1],
                              'note':'Analytic truncation bound is separate from finite eigensolver floating-point error; no infinite-dimensional operator-norm claim.'},
        'refinement':[comparison(n) for n in (32, 64, 128, 256)],
        'connection_negative_control':{
            'correct_initial_norm_rate':right_rate,
            'omitted_connection_initial_norm_rate':wrong_rate,
            'independent_integral_rate':direct_wrong},
        'massless_exact_shift_error':float(np.max(abs(massless-expected))),
        'geometry_interpretation_degeneracy':{
            'curved':'Specified static conformal metric with constant rest mass m and densitized spinor chi.',
            'flat':'Flat metric in (t,y) with specified external mass profile M(y)=m*v(y).',
            'same_propagation_equation':True,
            'missing_discriminator':'Independent physical rods, clock standards and universal matter coupling, or geometry dynamics.'},
        'added_inputs':['static positive bounded lapse v', 'optical and physical coordinate identification', 'single-particle two-component state and mass parameter', 'periodic boundary and smooth initial data', 'dt=optical spacing and position-dependent unitary coins', 'flat dx probability density convention'],
        'limitations':['No position-dependent geometry is selected by FUCP.', 'The target lapse is prescribed and has no matter backreaction equation.', 'The local norm connection is fixed by the chosen measure, not a proof of the full spin connection in every dimension.', 'Variable coefficients do not preserve a fixed Fourier band; estimates use H1 control instead.', 'The flat variable-mass description has the same selected propagation data.', 'No three-dimensional space, Einstein equation, cosmological constant, universal metric or autonomous infinite controller is derived.'],
        'sources':['https://arxiv.org/html/1505.07023', 'https://arxiv.org/html/1609.00305']}


class Checks(unittest.TestCase):
    def test_01_coordinate_map_positive_lapse_and_independent_derivative(self):
        y = np.linspace(0., 1., 31, endpoint=False)
        x, v, _, _ = background(y)
        h = 1e-6
        numerical = (background(y+h)[0]-background(y-h)[0])/(2*h)
        np.testing.assert_allclose(numerical, v, atol=2e-10)
        self.assertGreater(np.min(v), .48)
        self.assertLess(np.max(v), .92)
        self.assertAlmostEqual(background(1.)[0]-background(0.)[0], V0)

    def test_02_one_step_unitarity_and_reverse_on_arbitrary_vectors(self):
        rng = np.random.default_rng(349002)
        chi = rng.normal(size=(37, 2))+1j*rng.normal(size=(37, 2))
        result = one_step(chi)
        self.assertAlmostEqual(np.linalg.norm(result), np.linalg.norm(chi), places=13)
        np.testing.assert_allclose(inverse_step(result), chi, atol=2e-15)

    def test_03_strict_step_support_and_physical_distance(self):
        count = 64
        chi = np.zeros((count, 2), dtype=complex)
        chi[20, :] = 1/np.sqrt(2)
        for _ in range(5):
            chi = one_step(chi)
        outside = np.r_[0:15, 26:64]
        np.testing.assert_array_equal(chi[outside], 0.)
        x = background(np.arange(count+1)/count)[0]
        self.assertLessEqual(np.max(np.diff(x)), V0*(1+EPS)/count)

    def test_04_probability_measure_transformation_is_exact_on_the_grid(self):
        y = np.arange(128)/128
        v = background(y)[1]
        chi = solve_walk(128)
        psi = chi/np.sqrt(v)[:, None]
        self.assertAlmostEqual(np.mean(np.sum(abs(chi)**2, axis=1)),
                               np.mean(v*np.sum(abs(psi)**2, axis=1)), places=14)

    def test_05_independent_generators_agree_under_the_field_map(self):
        y = np.arange(1024)/1024
        v = background(y)[1]
        chi = initial(y)
        psi = chi/np.sqrt(v)[:, None]
        transformed = np.sqrt(v)[:, None]*physical_generator(psi, y)
        np.testing.assert_allclose(transformed, optical_generator(chi), atol=2e-12)

    def test_06_connection_is_needed_for_norm_and_local_continuity(self):
        count = 1024
        y = np.arange(count)/count
        _, v, vy, _ = background(y)
        chi = initial(y)
        psi = chi/np.sqrt(v)[:, None]
        dpsi = physical_generator(psi, y)
        rho_t = 2*np.sum((psi.conj()*dpsi).real, axis=1)
        current = v*(abs(psi[:, 0])**2-abs(psi[:, 1])**2)
        current_x = derivative_samples(current[:, None])[:, 0].real/v
        np.testing.assert_allclose(rho_t+current_x, 0., atol=3e-12)
        right_rate = 2*np.vdot(chi, optical_generator(chi)).real/count
        wrong_rate = 2*np.vdot(chi, optical_generator(chi, omit_connection=True)).real/count
        exact_wrong = np.mean(vy/v*(abs(chi[:, 0])**2-abs(chi[:, 1])**2))
        self.assertLess(abs(right_rate), 1e-13)
        self.assertAlmostEqual(wrong_rate, exact_wrong, places=13)
        self.assertGreater(abs(wrong_rate), .1)

    def test_07_fourier_H_is_Hermitian_and_matches_spectral_position_generator(self):
        cutoff = 10
        y = np.arange(128)/128
        coefficients = initial_modes(cutoff)
        h_action = (-1j*fourier_hamiltonian(cutoff) @ coefficients.reshape(-1)).reshape(-1, 2)
        np.testing.assert_allclose(evaluate_modes(h_action, y),
                                   optical_generator(initial(y)), atol=5e-13)
        h = fourier_hamiltonian(cutoff)
        np.testing.assert_array_equal(h, h.conj().T)

    def test_08_Galerkin_reference_refines_beyond_working_accuracy(self):
        y = np.arange(256)/256
        a = evaluate_modes(continuum_coefficients(12), y)
        b = evaluate_modes(continuum_coefficients(18), y)
        self.assertLess(np.linalg.norm(a-b)/16, 5e-13)
        self.assertLess(galerkin_tail_bounds(12)[0], 1e-20)

    def test_09_first_order_convergence_obeys_the_analytic_H1_bound(self):
        rows = [comparison(n) for n in (32, 64, 128, 256)]
        for row in rows:
            self.assertLess(row['state_L2_error'], row['sampled_analytic_error_upper'])
            self.assertLess(abs(row['weighted_physical_norm']-1), 3e-14)
        for coarse, fine in zip(rows, rows[1:]):
            ratio = coarse['state_L2_error']/fine['state_L2_error']
            self.assertGreater(ratio, 1.9)
            self.assertLess(ratio, 2.1)

    def test_10_initial_H1_norm_and_later_regular_state_bound(self):
        count = 1024
        y = np.arange(count)/count
        self.assertAlmostEqual(np.linalg.norm(derivative_samples(initial(y)))/np.sqrt(count),
                               derivative_initial_norm(), places=12)
        coefficients = continuum_coefficients()
        modes = np.arange(-18, 19)
        later = np.linalg.norm(2j*np.pi*modes[:, None]*coefficients)
        self.assertLessEqual(later, error_constants()['derivative_bound_M1'])

    def test_11_massless_optical_shift_is_exact_despite_nonconstant_curvature(self):
        count = 128
        y = np.arange(count)/count
        target = initial(y)
        target[:, 0] = np.roll(target[:, 0], 64)
        target[:, 1] = np.roll(target[:, 1], -64)
        np.testing.assert_array_equal(solve_walk(count, mass=0.), target)
        _, v, vy, vyy = background(y)
        curvature = -2*(vyy*v-vy**2)/v**4
        self.assertGreater(np.ptp(curvature), 1.)

    def test_12_connection_contraction_and_scalar_curvature_formula(self):
        for y in (0., .13, .25, .49, .8):
            _, v, vy, vyy = background(y)
            vx = vy/v
            vxx = (vyy*v-vy*vy)/v**3
            self.assertAlmostEqual(ricci_from_connection(v, vx, vxx), -2*vxx/v,
                                   places=11)

    def test_13_nonconstant_mass_generates_modes_outside_the_initial_band(self):
        coefficients = continuum_coefficients()
        modes = np.arange(-18, 19)
        outside = float(np.sum(abs(coefficients[abs(modes) > 1])**2))
        self.assertGreater(outside, 1e-7)
        self.assertLess(outside, .01)

    def test_14_the_physical_unweighted_grid_norm_is_not_the_probability(self):
        row = comparison(128)
        self.assertGreater(abs(row['unweighted_physical_samples_sum']-1), .1)
        self.assertLess(abs(row['weighted_physical_norm']-1), 3e-14)


if __name__ == '__main__':
    main(__name__, 'curved_quantum_walk_audit', report)
