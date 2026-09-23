"""Round 308: exact commuting kernel, restricted modular flow and continuum cancellation."""
import math
import unittest
import numpy as np
from growing_stream_audit import main
from ground_state_modular_audit import interval_correlation, modular_matrix


def commuting_kernel(n):
    i = np.arange(1, n)
    bonds = i*(n-i)/n**2
    return np.diag(bonds, 1)+np.diag(bonds, -1)


def window_data(n, cutoff=4., parameter=1.):
    c, t = interval_correlation(n), commuting_kernel(n)
    lam, u = np.linalg.eigh(t)
    transformed = u.T@c@u
    z = np.diag(transformed)
    lower = 1/(1+np.exp(cutoff))
    keep = (z >= lower) & (z <= 1-lower)
    exact = np.log((1-z[keep])/z[keep])
    geometric = -np.pi*n*lam[keep]
    delta = exact-geometric
    return {'sites': n, 'cutoff_modular_frequency': cutoff, 'retained_modes': int(keep.sum()),
            'endpoint_modes_unresolved_at_1e_12': int(np.count_nonzero((z <= 1e-12) | (z >= 1-1e-12))),
            'commutator_fro': float(np.linalg.norm(t@c-c@t)),
            'common_basis_c_offdiagonal_fro': float(np.linalg.norm(transformed-np.diag(z))),
            'exact_window_frequencies': exact.tolist(),
            'geometric_window_frequencies': geometric.tolist(),
            'max_generator_error': float(np.max(abs(delta))),
            'flow_parameter': parameter,
            'flow_error_op_on_single_particle_window': float(np.max(2*abs(np.sin(parameter*delta/2)))),
            'flow_error_linear_bound': float(abs(parameter)*np.max(abs(delta)))}


def cancellation(m):
    return sum((-1)**p*(2*p+1)*math.comb(2*m+1, m-p) for p in range(m+1))


def series_coefficient(m):
    # alpha_m * beta_m in the known large-N expansion; never treated as exact finite N.
    return math.exp(math.lgamma(m+0.5)-math.lgamma(m+1)+2*m*math.log(2)
                    +math.lgamma(2*m+0.5)-math.lgamma(2*m+2))


def hopping_sum_velocity(z, degree, nearest_only=False):
    total = 0.
    for m in range(degree+1):
        weight = math.comb(2*m+1, m) if nearest_only else cancellation(m)
        total += 2*series_coefficient(m)*weight*z**(2*m+1)
    return total


def report():
    full = []
    for n in (8, 12, 16):
        exact = modular_matrix(interval_correlation(n))
        geometric = -np.pi*n*commuting_kernel(n)
        full.append({'sites': n, 'full_matrix_error_op': float(np.linalg.norm(exact-geometric, 2)),
                     'center_bond_relative_excess': float(exact[n//2-1, n//2]/geometric[n//2-1, n//2]-1)})
    series = []
    for z in (0.1, 0.2, 0.25):
        target = 2*np.pi*z
        for degree in (1, 3, 8):
            all_hoppings = hopping_sum_velocity(z, degree)
            nearest = hopping_sum_velocity(z, degree, True)
            series.append({'z': z, 'maximum_power': 2*degree+1,
                           'all_hoppings_velocity_per_length': all_hoppings,
                           'nearest_only_velocity_per_length': nearest,
                           'geometric_target_per_length': target,
                           'all_hoppings_error': abs(all_hoppings-target),
                           'nearest_only_relative_excess': nearest/target-1})
    return {'round': 308,
            'scope': 'Half-filled 1D free-fermion interval: exact lattice commutation, numerical bounded modular-frequency flow, and a known continuum series cancellation. Not a 3+1D boost or gravity derivation.',
            'full_lattice_counterexamples': full,
            'restricted_flow': [window_data(n) for n in (16, 32, 64, 128, 256)],
            'integer_cancellation_m_0_through_40': [cancellation(m) for m in range(41)],
            'finite_polynomial_continuum_comparison': series,
            'limitations': ['Fixed modular frequency is not a physical energy cutoff.',
                            'The flow norm is on single-particle mode functions, not unrestricted Fock space.',
                            'Only resolved central occupations are logged; endpoint modes are not clipped.',
                            'The large-N power series is an external analytical result; finite-polynomial cancellation alone does not justify exchanging infinite limits.',
                            'No area entropy coefficient, vacuum uniqueness, or entanglement stationarity is derived.']}


class Checks(unittest.TestCase):
    def test_01_exact_commutation_calibration(self):
        for n in (7, 16, 31, 64):
            c, t = interval_correlation(n), commuting_kernel(n)
            np.testing.assert_allclose(c@t, t@c, atol=1e-14)

    def test_02_finite_matrix_not_equal(self):
        rows = report()['full_lattice_counterexamples']
        self.assertGreater(rows[0]['full_matrix_error_op'], 0.5)
        self.assertTrue(all(b['full_matrix_error_op'] > a['full_matrix_error_op'] for a, b in zip(rows, rows[1:])))

    def test_03_common_eigenvectors_resolved(self):
        data = window_data(128)
        self.assertLess(data['common_basis_c_offdiagonal_fro'], 2e-12)
        self.assertGreater(data['endpoint_modes_unresolved_at_1e_12'], 0)

    def test_04_central_frequencies_direct_log(self):
        n = 12
        k = modular_matrix(interval_correlation(n))
        lam, u = np.linalg.eigh(commuting_kernel(n))
        exact = np.diag(u.T@k@u)
        keep = abs(exact) <= 4
        np.testing.assert_allclose(exact[keep], window_data(n)['exact_window_frequencies'], atol=1e-9)

    def test_05_fixed_window_refinement(self):
        rows = [window_data(n) for n in (16, 32, 64, 128, 256)]
        errors = [r['max_generator_error'] for r in rows]
        self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])))
        self.assertLess(errors[-1], 5e-6)
        self.assertGreater(errors[0]/errors[-1], 100)

    def test_06_projected_flow_bound(self):
        for s in (-2., 0., 0.5, 3.):
            row = window_data(32, parameter=s)
            self.assertLessEqual(row['flow_error_op_on_single_particle_window'], row['flow_error_linear_bound']+1e-14)

    def test_07_integer_cancellation(self):
        self.assertEqual([cancellation(m) for m in range(41)], [1]+[0]*40)

    def test_08_all_hoppings_recover_leading_velocity(self):
        for z in (0.1, 0.2, 0.25):
            for degree in (0, 1, 3, 8):
                self.assertAlmostEqual(hopping_sum_velocity(z, degree), 2*np.pi*z, places=13)

    def test_09_nearest_only_bias(self):
        target = 2*np.pi*0.25
        self.assertGreater(hopping_sum_velocity(0.25, 8, True)/target-1, 0.05)
        self.assertAlmostEqual(series_coefficient(0), np.pi)

    def test_10_independent_matrix_flow(self):
        n, s = 12, 3.
        k = modular_matrix(interval_correlation(n))
        t = commuting_kernel(n)
        _, u = np.linalg.eigh(t)
        exact = np.diag(u.T@k@u)
        selected = u[:, abs(exact) <= 4]
        projector = selected@selected.T
        def exponential(h):
            values, vectors = np.linalg.eigh(h)
            return (vectors*np.exp(-1j*s*values))@vectors.T
        direct = np.linalg.norm((exponential(k)-exponential(-np.pi*n*t))@projector, 2)
        self.assertAlmostEqual(direct, window_data(n, parameter=s)['flow_error_op_on_single_particle_window'], places=9)


if __name__ == '__main__':
    main(__name__, 'geometric_modular_window_audit', report)
