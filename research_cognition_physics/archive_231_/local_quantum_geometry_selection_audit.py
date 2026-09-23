"""Round 343: full quantum operations and local drift do not select geometry.

No restriction to the one-excitation sector is used for the locality theorem.
That invariant sector supplies an operationally preparable response diagnostic.
"""
from functools import lru_cache
from itertools import product
from math import comb, factorial
import unittest

import numpy as np

from growing_stream_audit import main


def torus(dimension, side):
    if dimension < 1 or side < 3:
        raise ValueError("Use positive d and L >= 3; L=2 merges the two neighbors.")
    points = list(product(range(side), repeat=dimension))
    index = {p: i for i, p in enumerate(points)}
    edges = set()
    for i, p in enumerate(points):
        for axis in range(dimension):
            q = list(p)
            q[axis] = (q[axis]+1) % side
            edges.add(tuple(sorted((i, index[tuple(q)]))))
    return points, sorted(edges)


def laplacian(count, edges):
    answer = np.zeros((count, count))
    for x, y in edges:
        answer[x, x] += 1.
        answer[y, y] += 1.
        answer[x, y] -= 1.
        answer[y, x] -= 1.
    return answer


def full_exchange(count, edges, coupling=1.):
    """Qubit bit x is site x; positive bond term is J (I-SWAP)."""
    size = 1 << count
    h = np.zeros((size, size))
    for x, y in edges:
        for bits in range(size):
            if ((bits >> x) & 1) != ((bits >> y) & 1):
                swapped = bits ^ (1 << x) ^ (1 << y)
                h[bits, bits] += coupling
                h[swapped, bits] -= coupling
    return h


def local_pauli(count, site, axis):
    size = 1 << count
    out = np.zeros((size, size), dtype=complex)
    for bits in range(size):
        sign = 1-2*((bits >> site) & 1)
        if axis == 'z':
            out[bits, bits] = sign
        elif axis == 'x':
            out[bits ^ (1 << site), bits] = 1.
        elif axis == 'y':
            out[bits ^ (1 << site), bits] = 1j*sign
        else:
            raise ValueError(axis)
    return out


def unitary(h, time):
    energy, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*energy*time)) @ vectors.conj().T


def excitation_indices(count):
    return np.array([1 << x for x in range(count)])


def dispersion(k, coupling):
    # sin form remains accurate in the low-momentum limit.
    return float(4*coupling*np.sum(np.sin(np.asarray(k)/2)**2))


def dispersion_remainder_bound(k, coupling):
    return float(coupling*np.sum(np.asarray(k)**4)/12)


def group_velocity(k, coupling):
    return 2*coupling*np.sin(np.asarray(k))


def ball_count(dimension, radius):
    return sum(2**k*comb(dimension, k)*comb(radius, k)
               for k in range(min(dimension, radius)+1))


def ball_count_enumerated(dimension, radius):
    return sum(sum(abs(x) for x in p) <= radius
               for p in product(range(-radius, radius+1), repeat=dimension))


def poisson_series_tail(rate, distance):
    """Sum_{n>=distance} rate**n/n!, without small-rate cancellation."""
    term = rate**distance/factorial(distance)
    result = term
    for n in range(distance+1, distance+1000):
        term *= rate/n
        result += term
        if abs(term) < 1e-17*max(abs(result), 1e-300):
            return result
    raise RuntimeError("Series did not converge")


def lr_bound(time, distance, max_degree, coupling, support_size=1):
    # From NOS (2006), recurrence (2.9)-(2.12):
    # first edge <= |X| Delta, later edge <= 2 Delta, norm per edge 2J.
    # Loose chain count a_n <= |X| (4 Delta J)^n; Duhamel contributes (2t)^n.
    return 2*support_size*poisson_series_tail(
        8*max_degree*coupling*abs(time), distance)


def lr_exponential_bound(time, distance, max_degree, coupling,
                         support_size=1, mu=1.):
    return float(2*support_size*np.exp(
        -mu*distance+8*max_degree*coupling*np.exp(mu)*abs(time)))


def steering_data(rho, pieces):
    vals, vecs = np.linalg.eigh(rho)
    root = (vecs*np.sqrt(vals)) @ vecs.conj().T
    invroot = (vecs*(1/np.sqrt(vals))) @ vecs.conj().T
    effects = [(invroot @ p @ invroot).T for p in pieces]
    recovered = [root @ e.T @ root for e in effects]
    return effects, recovered


@lru_cache(maxsize=None)
def full_commutator_samples():
    count = 6
    edges = [(i, i+1) for i in range(count-1)]
    h = full_exchange(count, edges)
    a = local_pauli(count, 0, 'x')
    rows = []
    for site in (1, 2, 3, 5):
        b = local_pauli(count, site, 'z')
        for t in (.002, .02):
            u = unitary(h, t)
            evolved = u.conj().T @ a @ u
            value = float(np.linalg.norm(evolved @ b-b @ evolved, 2))
            bound = lr_bound(t, site, 2, 1.)
            rows.append({'distance':site, 'time':t,
                         'full_operator_commutator_norm':value,
                         'path_series_bound':bound})
    return rows


def report():
    rows = []
    for d, side in ((1, 64), (2, 8), (3, 4), (4, 4)):
        pts, edges = torus(d, side)
        coupling = 1./d
        energies = np.linalg.eigvalsh(coupling*laplacian(len(pts), edges))
        positive = energies[energies > 1e-10]
        rows.append({'dimension':d, 'side':side, 'qubits':len(pts),
                     'degree':2*d, 'bond_coupling':coupling,
                     'incident_bond_norm_budget':4.,
                     'total_norm_upper_per_site':2.,
                     'excitation_gap':float(positive[0]),
                     'gap_formula':dispersion([2*np.pi/side], coupling),
                     'gap_multiplicity':int(np.sum(abs(energies-positive[0]) < 1e-9)),
                     'volume_growth_exponent':d,
                     'gapless_response_dynamical_exponent':2,
                     'max_l1_group_velocity_bound':2.,
                     'LR_mu_one_speed_upper':float(16*np.e)})
    small = []
    for radius in (4, 16, 64):
        small.append({'radius':radius,
                      'ball_counts':{str(d):ball_count(d, radius) for d in (1, 2, 3, 4)}})
    krows = []
    for wave in (.2, .1, .05):
        energy = dispersion([wave, 0., 0.], 1/3)
        krows.append({'k':wave, 'energy':energy,
                      'energy_over_abs_k':energy/wave,
                      'energy_over_k_squared':energy/wave**2})
    return {
        'round':343,
        'scope':'All finite tensor-qubit quantum objects, all density matrices and finite CP instruments, with a separately specified nonconstant natural exchange Hamiltonian on arbitrarily large d-dimensional graphs. Full-spin Lieb-Robinson locality is proved; a preparable invariant excitation response diagnoses the selected locality. This is not a no-go theorem for every encoded or emergent GR sector.',
        'full_contract_witness':{
            'F':'Ordinary finite quantum states, effects, instruments, tensor products, discard.',
            'U':'Same-size purifying partner; E_i=(rho^(-1/2) tau_i rho^(-1/2))^T on the support, completed on the kernel.',
            'C':'Full projective unitary groups, including every finite composite.',
            'P':'Local extensions of the same CP map coincide on all joint states.',
            'Time':'exp(-itH) is a nonconstant continuous one-parameter unitary group.'},
        'uniform_resource_class':{
            'choice':'J_d=J_star/d, J_star=1',
            'per_site_incident_interaction_norm':4.,
            'per_site_total_H_norm_upper':2.,
            'LR_speed_bound_at_mu_1':float(16*np.e),
            'note':'The bounds do not depend on volume or dimension. Site degree itself is 2d and is not fixed by the operation contract.'},
        'finite_volume_spectral_witnesses':rows,
        'lattice_ball_counts':small,
        'low_momentum_response':krows,
        'full_algebra_commutator_checks':full_commutator_samples(),
        'logical_conclusion':{
            'disproved':'F+U+C+P+nontrivial Time plus finite-range bounded local drift with a volume-uniform LR estimate uniquely selects d=3 or linear gapless dispersion in its directly measured spatial response.',
            'not_disproved':'There exists a further operationally justified sector, state, coarse-graining and scaling prescription yielding GR; an enlarged cognitive contract may select it.',
            'strict_causality':'Lieb-Robinson exponential suppression is not exact vanishing outside a Lorentz light cone.'},
        'added_inputs':['tensor-site identification and graph family', 'natural Hamiltonian and coupling normalization', 'polarized reference and local excitation readout', 'the identification of physical distance with the recovered natural-interaction locality', 'thermodynamic and low-momentum limiting prescriptions'],
        'primary_source':'Nachtergaele, Ogata and Sims (2006), https://arxiv.org/html/math-ph/0603064, Theorems 2.1/2.2 and recurrence (2.9)-(2.12).'}


class Checks(unittest.TestCase):
    def test_01_torus_graphs_are_regular_and_closed(self):
        for d, side in ((1, 7), (2, 4), (3, 3), (4, 3)):
            pts, edges = torus(d, side)
            lap = laplacian(len(pts), edges)
            np.testing.assert_array_equal(lap.sum(axis=1), 0.)
            np.testing.assert_array_equal(np.diag(lap), 2*d)
            self.assertEqual(len(edges), d*side**d)

    def test_02_full_bond_is_positive_bounded_entangling(self):
        h = full_exchange(2, [(0, 1)])
        np.testing.assert_allclose(np.linalg.eigvalsh(h), [0., 0., 0., 2.], atol=1e-14)
        vector = unitary(h, np.pi/4)[:, 1].reshape(2, 2)
        reduced = vector @ vector.conj().T
        np.testing.assert_allclose(reduced, np.eye(2)/2, atol=1e-14)

    def test_03_full_H_preserves_number_but_operations_need_not(self):
        count = 4
        h = full_exchange(count, [(0, 1), (1, 2), (2, 3), (3, 0)])
        n = np.diag([i.bit_count() for i in range(1 << count)])
        np.testing.assert_array_equal(h @ n-n @ h, 0.)
        x = local_pauli(count, 0, 'x')
        self.assertGreater(np.linalg.norm(x @ n-n @ x), 1.)

    def test_04_excitation_sector_is_extracted_from_the_full_spin_model(self):
        pts, edges = torus(1, 5)
        h = full_exchange(len(pts), edges, .7)
        inds = excitation_indices(len(pts))
        np.testing.assert_allclose(h[np.ix_(inds, inds)],
                                   .7*laplacian(len(pts), edges), atol=1e-14)
        complement = np.setdiff1d(np.arange(len(h)), inds)
        np.testing.assert_array_equal(h[np.ix_(complement, inds)], 0.)
        self.assertGreaterEqual(np.linalg.eigvalsh(h)[0], -1e-12)

    def test_05_full_unitary_response_equals_sector_propagation(self):
        pts, edges = torus(1, 5)
        h = full_exchange(len(pts), edges)
        inds = excitation_indices(len(pts))
        for t in (.03, .4):
            full = unitary(h, t)[np.ix_(inds, inds)]
            projected = unitary(laplacian(len(pts), edges), t)
            np.testing.assert_allclose(full, projected, atol=3e-14)

    def test_06_independent_fourier_spectrum_all_dimensions(self):
        for d, side in ((1, 13), (2, 5), (3, 4), (4, 3)):
            pts, edges = torus(d, side)
            matrix_spectrum = np.linalg.eigvalsh(laplacian(len(pts), edges)/d)
            fourier = sorted(dispersion(2*np.pi*np.asarray(k)/side, 1/d) for k in pts)
            np.testing.assert_allclose(matrix_spectrum, fourier, atol=3e-14)

    def test_07_short_time_born_probabilities_recover_adjacency(self):
        pts, edges = torus(2, 5)
        h = laplacian(len(pts), edges)/2
        t = 1e-3
        u = unitary(h, t)
        inferred = abs(u[:, 0])**2/t**2
        target = (h[:, 0] < 0).astype(float)/4
        inferred[0] = target[0] = 0.
        np.testing.assert_allclose(inferred, target, atol=3e-7)
        self.assertEqual(np.count_nonzero(inferred > .1), 4)

    def test_08_operational_ball_count_exact_and_large_scale_exponent(self):
        for d in (1, 2, 3, 4):
            for radius in (0, 1, 2, 3):
                self.assertEqual(ball_count(d, radius), ball_count_enumerated(d, radius))
            r = 100000
            ratio = ball_count(d, 2*r)/ball_count(d, r)
            self.assertLess(abs(ratio/2**d-1), 3e-5)

    def test_09_quadratic_gapless_response_has_controlled_remainder(self):
        rng = np.random.default_rng(343009)
        for d in (1, 2, 3, 4):
            for _ in range(20):
                k = rng.uniform(-.2, .2, d)
                quadratic = np.dot(k, k)/d
                error = quadratic-dispersion(k, 1/d)
                self.assertGreaterEqual(error, -1e-15)
                self.assertLessEqual(error, dispersion_remainder_bound(k, 1/d)+1e-15)
            self.assertAlmostEqual(dispersion([1e-3], 1/d)/dispersion([5e-4], 1/d),
                                   4., places=6)

    def test_10_resource_and_speed_bounds_are_dimension_independent(self):
        for d in (1, 2, 3, 4, 10, 100):
            j = 1/d
            self.assertAlmostEqual(2*d*(2*j), 4.)
            self.assertAlmostEqual(d*(2*j), 2.)
            self.assertAlmostEqual(8*(2*d)*j*np.e, 16*np.e)
            self.assertLessEqual(np.sum(abs(group_velocity(np.full(d, np.pi/2), j))), 2+1e-14)

    def test_11_full_algebra_commutators_obey_the_path_bound(self):
        for row in full_commutator_samples():
            self.assertLessEqual(row['full_operator_commutator_norm'],
                                 row['path_series_bound']+3e-14)
            self.assertLessEqual(row['path_series_bound'],
                                 lr_exponential_bound(row['time'], row['distance'], 2, 1.)+1e-14)

    def test_12_finite_range_time_evolution_has_tails_not_a_strict_cone(self):
        h = laplacian(7, [(x, x+1) for x in range(6)])
        distance = 3
        for n in range(distance):
            self.assertEqual(np.linalg.matrix_power(h, n)[distance, 0], 0.)
        self.assertEqual(np.linalg.matrix_power(h, distance)[distance, 0], -1.)
        for t in (.002, .001):
            amplitude = unitary(h, t)[distance, 0]
            leading = (1j*t)**distance/factorial(distance)
            self.assertGreater(abs(amplitude), 0.)
            self.assertLess(abs(amplitude/leading-1), .006)

    def test_13_same_size_steering_uses_the_whole_quantum_state_space(self):
        rng = np.random.default_rng(343013)
        vectors = rng.normal(size=(5, 4))+1j*rng.normal(size=(5, 4))
        pieces = [np.outer(v, v.conj()) for v in vectors]
        normalization = sum(np.trace(p).real for p in pieces)
        pieces = [p/normalization for p in pieces]
        rho = sum(pieces)
        effects, recovered = steering_data(rho, pieces)
        np.testing.assert_allclose(sum(effects), np.eye(4), atol=3e-14)
        for effect, target, actual in zip(effects, pieces, recovered):
            self.assertGreaterEqual(np.linalg.eigvalsh(effect)[0], -2e-14)
            np.testing.assert_allclose(actual, target, atol=3e-14)

    def test_14_full_drift_has_nontrivial_Time_group_law(self):
        h = full_exchange(4, [(0, 1), (1, 2), (2, 3)])
        np.testing.assert_allclose(unitary(h, .37) @ unitary(h, -.12),
                                   unitary(h, .25), atol=3e-14)
        self.assertGreater(np.linalg.norm(unitary(h, .2)-np.eye(len(h))), .5)


if __name__ == '__main__':
    main(__name__, 'local_quantum_geometry_selection_audit', report)
