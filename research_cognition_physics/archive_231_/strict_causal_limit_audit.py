"""Round 347: strict finite delay, matrix analyticity and a controlled tail limit."""
from functools import lru_cache
from math import factorial
import unittest

import numpy as np

from growing_stream_audit import main


I2 = np.eye(2, dtype=complex)
X = np.array([[0., 1.], [1., 0.]], dtype=complex)
Y = np.array([[0., -1j], [1j, 0.]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
PAULI = (I2, X, Y, Z)


def comm(a, b):
    return a @ b-b @ a


def evolved(h, a, time):
    """Entire extension uses inverse, not adjoint, at complex time."""
    vals, vecs = np.linalg.eigh(h)
    left = (vecs*np.exp(1j*time*vals)) @ vecs.conj().T
    right = (vecs*np.exp(-1j*time*vals)) @ vecs.conj().T
    return left @ a @ right


def derivative_commutators(h, a, b, orders):
    out = []
    current = a.copy()
    for _ in range(orders):
        out.append(comm(current, b))
        current = 1j*comm(h, current)
    return out


def matrix_time_series(h, a, time, orders=50):
    out = a.copy()
    term = a.copy()
    for n in range(1, orders):
        term = (1j*time/n)*comm(h, term)
        out += term
    return out


def interaction_part_two_qubits(h):
    return sum(np.trace(np.kron(a, b) @ h)/4*np.kron(a, b)
               for a in PAULI[1:] for b in PAULI[1:])


def all_pair_first_derivative_norm(h):
    return max(float(np.linalg.norm(comm(comm(h, np.kron(a, I2)),
                                        np.kron(I2, b)), 2))
               for a in PAULI[1:] for b in PAULI[1:])


def cyclic_hamiltonian(size, spacing, speed=1.):
    shift = np.roll(np.eye(size), 1, axis=0)
    return speed/(2*spacing)*(2*np.eye(size)-shift-shift.T)


def ring_kernel(size, z):
    theta = 2*np.pi*np.arange(size)/size
    return np.fft.ifft(np.exp(-1j*z*(1-np.cos(theta))))


def bessel_real_quadrature(order, z, nodes=8192):
    theta = 2*np.pi*np.arange(nodes)/nodes
    return np.mean(np.exp(1j*(order*theta-z*np.sin(theta))))


def bessel_tail_contour(order, z, nodes=2048):
    """For n>z>=0, shifted integral computes tiny tails without cancellation."""
    if z == 0.:
        return complex(order == 0)
    if order <= z:
        raise ValueError("This saddle contour is for n > z > 0.")
    shift = np.arccosh(order/z)
    root = np.sqrt(order*order-z*z)
    theta = 2*np.pi*np.arange(nodes)/nodes
    scaled = np.exp(root*(np.cos(theta)-1)
                    +1j*order*(theta-np.sin(theta)))
    return np.exp(-order*shift+root)*np.mean(scaled)


def log_endpoint_bound(order, z):
    if not order > z > 0:
        raise ValueError("Require n > z > 0.")
    return float(-order*np.arccosh(order/z)+np.sqrt(order*order-z*z))


def rate(distance, speed_time):
    if not distance > speed_time > 0:
        raise ValueError("Strictly outside the ballistic cone.")
    return float(distance*np.arccosh(distance/speed_time)
                 -np.sqrt(distance*distance-speed_time*speed_time))


def log_total_tail_bound(first_order, z):
    shift = np.arccosh(first_order/z)
    return float(np.log(2)+2*log_endpoint_bound(first_order, z)
                 -np.log(-np.expm1(-2*shift)))


def log_alias_bound(order, z, size):
    distance = size-abs(order)
    if distance <= z:
        raise ValueError("The first periodic image must be outside its cone.")
    shift = np.arccosh(distance/z)
    # Images m!=0 have distances >= size-|n|+(|m|-1)*size.
    return float(np.log(2)-shift*distance+z*np.sinh(shift)
                 -np.log(-np.expm1(-shift*size)))


@lru_cache(maxsize=None)
def continuum_row(spacing, distance=1.5, time=1., speed=1.):
    ratio = distance/spacing
    first = int(round(ratio))
    if abs(first-ratio) > 1e-10:
        raise ValueError("Choose a threshold exactly on a site for this table.")
    z = speed*time/spacing
    terms = 64
    values = [bessel_tail_contour(n, z) for n in range(first, first+terms)]
    tail = float(2*sum(abs(v)**2 for v in values))
    size = 4096
    fft_value = ring_kernel(size, z)[first]
    exact = np.exp(-1j*z)*(1j)**first*values[0]
    return {'spacing':spacing, 'first_outside_site':first, 'z':z,
            'coupling_J':speed/(2*spacing), 'physical_threshold':distance,
            'localized_preparation_mean_energy':speed/spacing,
            'single_endpoint_probability':float(abs(values[0])**2),
            'two_sided_tail_probability_partial_sum':tail,
            'log_omitted_probability_upper':log_total_tail_bound(first+terms, z),
            'two_sided_tail_probability_upper':float(np.exp(log_total_tail_bound(first, z))),
            'rate_I':rate(distance, speed*time),
            'minus_spacing_log_tail':float(-spacing*np.log(tail)),
            'FFT_endpoint_absolute_discrepancy':float(abs(fft_value-exact)),
            'log_periodic_image_error_upper':log_alias_bound(first, z, size),
            'FFT_is_relative_tail_measurement':bool(abs(exact) > 1e-12),
            'warning':'Contour quadrature is numerical, not interval arithmetic; analytic inequalities provide the rigorous bounds.'}


def report():
    h = np.kron(X, X)+.3*np.kron(Z, I2)
    a = np.kron(Z, I2)
    b = np.kron(I2, Z)
    ds = derivative_commutators(h, a, b, 16)
    onsite = .7*np.kron(X, I2)+.4*np.kron(I2, Z)
    rows = [continuum_row(1/n) for n in (4, 8, 16, 32, 64, 128)]
    z = 32.
    kernel = ring_kernel(4096, z)
    sites = np.fft.fftfreq(len(kernel))*len(kernel)
    moment = float(np.sum(abs(kernel)**2*(sites/32)**2))
    return {
        'round':347,
        'scope':'Fixed finite-dimensional autonomous Hamiltonian and fixed local algebras: exact commutator silence on any open time interval implies silence at every time. Separately, a specified one-dimensional tight-binding single-particle response admits an exponentially controlled continuum tail limit. The latter is not a full-algebra quantum-field microcausality proof or a derivation of GR.',
        'analytic_finite_delay_result':{
            'conditions':['finite Hilbert dimension', 'fixed time-independent H', 'fixed local A and B', 'positive interval of exact zero commutator'],
            'conclusion':'[exp(itH) A exp(-itH),B] is identically zero if it vanishes on an open real interval.',
            'all_pair_consequence':'If every pair of distinct tensor sites has a strictly positive exact delay for all local observables, H is a sum of on-site terms plus a scalar.',
            'finite_derivative_certificate':'For Hilbert dimension D, derivatives n=0,...,D^2-1 suffice for exact symbolic zero certification by Cayley-Hamilton; floating thresholds are not that certificate.'},
        'finite_matrix_examples':{
            'interacting_first_nonzero_derivative':next(i for i, d in enumerate(ds) if np.linalg.norm(d) > 1e-12),
            'interaction_part_norm':float(np.linalg.norm(interaction_part_two_qubits(h), 2)),
            'onsite_max_first_cross_derivative':all_pair_first_derivative_norm(onsite),
            'interacting_commutator_at_0_001':float(np.linalg.norm(comm(evolved(h, a, .001), b), 2))},
        'continuum_scaling':{
            'H':'H_a=c(2I-S-S_dagger)/(2a), J_a=c/(2a), hbar=1.',
            'kernel':'K_n(t)=exp(-iz) i^n J_n(z), z=ct/a.',
            'single_endpoint_bound':'|K_n| <= exp[-I(x,ct)/a], x=na>ct.',
            'I':'x*acosh(x/(ct))-sqrt(x*x-c*c*t*t).',
            'total_probability_bound':'P(|n|>=N) <= 2 exp[-2I(Na,ct)/a]/(1-exp[-2acosh(Na/(ct))]).',
            'per_bond_full_spin_norm_if_realized_by_exchange':'2J_a=c/a, not bounded uniformly as a->0.',
            'fixed_physical_low_momentum_energy':'c*(1-cos(a*k))/a -> 0; tail concentration alone is not a nontrivial Lorentz dispersion limit.'},
        'outside_cone_rows':rows,
        'nontrivial_spreading_check':{
            'spacing':1/32, 'time':1., 'speed':1.,
            'physical_second_moment':moment,
            'exact_second_moment':.5,
            'interpretation':'An initially site-localized state changes with the cutoff and has a broad momentum distribution; this spreading does not contradict freezing of fixed low physical momentum.'},
        'added_inputs':['one-dimensional lattice and local tensor identification', 'speed scale c', 'coupling cutoff dependence J_a=c/(2a)', 'localized reference preparation and single-particle readout', 'infinite-volume then small-spacing limiting prescription'],
        'limitations':['No finite numerical sampling can establish exact silence on an interval.', 'The positive tail-limit result is a specified propagation kernel, not all quantum observables or arbitrary many-body states.', 'The site-localized initial state has mean energy c/a; this is not a fixed finite-energy continuum preparation.', 'No unique spatial dimension, reference vacuum, Lorentz field equation, universal coupling or Einstein dynamics is derived.', 'Hegerfeldt positive-observable results are not interchangeable with operational signal differences or commutator conditions.'],
        'sources':['https://arxiv.org/html/quant-ph/9809030', 'https://arxiv.org/abs/2408.02108', 'https://www-fourier.univ-grenoble-alpes.fr/~joye/extail.pdf']}


class Checks(unittest.TestCase):
    def test_01_complex_time_entire_extension_matches_nested_series(self):
        h = np.kron(X, X)+.3*np.kron(Z, I2)
        a = np.kron(Y, I2)
        for t in (.2, -.3, .2+.1j):
            np.testing.assert_allclose(evolved(h, a, t), matrix_time_series(h, a, t),
                                       atol=2e-14)

    def test_02_an_interaction_has_an_immediate_nonzero_cross_derivative(self):
        h = np.kron(X, X)
        a = np.kron(Z, I2)
        b = np.kron(I2, Z)
        ds = derivative_commutators(h, a, b, 16)
        self.assertEqual(np.linalg.norm(ds[0]), 0.)
        self.assertAlmostEqual(np.linalg.norm(ds[1], 2), 4., places=13)
        for t in (.001, .01):
            numerical = np.linalg.norm(comm(evolved(h, a, t), b), 2)
            self.assertAlmostEqual(numerical, 2*abs(np.sin(2*t)), places=13)

    def test_03_silence_for_a_single_readout_does_not_remove_interactions(self):
        h = np.kron(Z, Z)
        a = np.kron(X, I2)
        b = np.kron(I2, Z)
        for derivative in derivative_commutators(h, a, b, 16):
            np.testing.assert_array_equal(derivative, 0.)
        for t in (.01, .8, 5.):
            np.testing.assert_allclose(comm(evolved(h, a, t), b), 0., atol=1e-14)
        self.assertGreater(all_pair_first_derivative_norm(h), 0.)

    def test_04_all_onsite_H_terms_preserve_all_remote_silence(self):
        h = .7*np.kron(X, I2)+.4*np.kron(I2, Z)
        self.assertEqual(all_pair_first_derivative_norm(h), 0.)
        np.testing.assert_array_equal(interaction_part_two_qubits(h), 0.)
        for a in PAULI[1:]:
            for b in PAULI[1:]:
                np.testing.assert_allclose(comm(evolved(h, np.kron(a, I2), .73),
                                                np.kron(I2, b)), 0., atol=2e-14)

    def test_05_every_two_site_Pauli_interaction_is_detected_by_the_all_pair_condition(self):
        for a in PAULI[1:]:
            for b in PAULI[1:]:
                h = np.kron(a, b)
                np.testing.assert_array_equal(interaction_part_two_qubits(h), h)
                self.assertAlmostEqual(all_pair_first_derivative_norm(h), 4., places=13)

    def test_06_periodic_matrix_exponential_equals_the_independent_FFT_kernel(self):
        h = cyclic_hamiltonian(31, .25)
        vals, vecs = np.linalg.eigh(h)
        direct = ((vecs*np.exp(-1j*.7*vals)) @ vecs.conj().T)[:, 0]
        np.testing.assert_allclose(direct, ring_kernel(31, .7/.25), atol=3e-14)

    def test_07_shifted_contour_matches_real_integral_at_resolvable_tails(self):
        for n, z in ((4, 2.), (8, 5.), (16, 10.), (24, 16.)):
            exact = bessel_real_quadrature(n, z)
            shifted = bessel_tail_contour(n, z)
            self.assertLess(abs(exact-shifted), 2e-15)
            self.assertLess(abs(shifted.imag), 1e-15)

    def test_08_saddle_is_the_optimal_positive_contour_shift(self):
        for n, z in ((6, 4.), (12, 8.), (48, 32.), (192, 128.)):
            s = np.arccosh(n/z)
            value = -n*s+z*np.sinh(s)
            self.assertAlmostEqual(value, log_endpoint_bound(n, z), places=12)
            self.assertAlmostEqual(-n+z*np.cosh(s), 0., places=12)
            for delta in (-.1, .1):
                self.assertLessEqual(value, -n*(s+delta)+z*np.sinh(s+delta))

    def test_09_single_point_and_total_tail_bounds_have_the_correct_probability_factor(self):
        for inv in (4, 8, 16, 32, 64, 128):
            row = continuum_row(1/inv)
            self.assertLess(row['single_endpoint_probability'],
                            np.exp(-2*row['rate_I']/row['spacing']))
            self.assertLess(row['two_sided_tail_probability_partial_sum'],
                            row['two_sided_tail_probability_upper'])
            self.assertLess(row['log_omitted_probability_upper'],
                            np.log(row['two_sided_tail_probability_partial_sum'])-100)

    def test_10_total_tail_refinement_really_concentrates_outside_a_fixed_cone(self):
        values = [continuum_row(1/n) for n in (4, 8, 16, 32, 64, 128)]
        self.assertTrue(all(b['two_sided_tail_probability_upper'] <
                            a['two_sided_tail_probability_upper']
                            for a, b in zip(values, values[1:])))
        target = 2*rate(1.5, 1.)
        errors = [abs(row['minus_spacing_log_tail']-target) for row in values]
        self.assertTrue(all(b < a for a, b in zip(errors, errors[1:])))

    def test_11_periodic_image_bound_and_resolvable_fft_agree(self):
        for inv in (4, 8, 16, 32):
            row = continuum_row(1/inv)
            self.assertLess(row['log_periodic_image_error_upper'], -1000.)
            self.assertLess(row['FFT_endpoint_absolute_discrepancy'], 3e-14)
        # Tiny FFT residuals are not used to measure the tiniest tail.
        self.assertFalse(continuum_row(1/128)['FFT_is_relative_tail_measurement'])

    def test_12_contour_quadrature_convergence_and_analytic_recurrence(self):
        for n, z in ((12, 8.), (48, 32.), (192, 128.)):
            fine = bessel_tail_contour(n, z, 4096)
            coarse = bessel_tail_contour(n, z, 2048)
            self.assertLess(abs(fine-coarse)/abs(fine), 1e-12)
            recurrence = (bessel_tail_contour(n-1, z)
                          +bessel_tail_contour(n+1, z)-2*n/z*fine)
            self.assertLess(abs(recurrence)/abs(fine), 1e-11)

    def test_13_normalization_and_nontrivial_physical_second_moment(self):
        for spacing in (.25, .125, .0625, .03125):
            kernel = ring_kernel(4096, 1/spacing)
            sites = np.fft.fftfreq(len(kernel))*len(kernel)
            self.assertAlmostEqual(np.sum(abs(kernel)**2), 1., places=13)
            self.assertAlmostEqual(np.sum(abs(kernel)**2*(spacing*sites)**2), .5,
                                   places=12)

    def test_14_speed_fixes_a_nonuniform_coupling_and_not_Lorentz_dispersion(self):
        previous = None
        for spacing in (.25, .125, .0625, .03125):
            coupling = 1/(2*spacing)
            self.assertAlmostEqual(2*coupling*spacing, 1.)
            energy = 2/spacing*np.sin(spacing*.4/2)**2
            self.assertLessEqual(energy, spacing*.4**2/2+1e-15)
            if previous is not None:
                self.assertLess(energy, previous)
            previous = energy


if __name__ == '__main__':
    main(__name__, 'strict_causal_limit_audit', report)
