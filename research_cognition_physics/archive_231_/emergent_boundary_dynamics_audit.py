"""Round 363: boundary Hamiltonians, physical algebras, and compression.

Finite matrix audit of the interfaces in Marolf's conditional obstruction.
The examples are not gravity models. Low-energy compression is not an
algebra homomorphism unless the required operators preserve the code.
"""
import itertools
import unittest
import numpy as np
from growing_stream_audit import main


I2 = np.eye(2, dtype=complex)
X = np.array([[0., 1.], [1., 0.]], dtype=complex)
Y = np.array([[0., -1j], [1j, 0.]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
PAULI = (I2, X, Y, Z)


def norm(matrix):
    return float(np.linalg.norm(matrix, ord=2))


def comm(a, b):
    return a @ b-b @ a


def kron_all(factors):
    answer = np.ones((1, 1), dtype=complex)
    for factor in factors:
        answer = np.kron(answer, factor)
    return answer


def at_site(operator, site, count):
    return kron_all([operator if i == site else I2 for i in range(count)])


def unitary(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*time*values)) @ vectors.conj().T


def evolve_operator(h, observable, time):
    u = unitary(h, time)
    return u.conj().T @ observable @ u


def compress(v, operator):
    return v.conj().T @ operator @ v


def defects(v, h, boundary, observable):
    p = v @ v.conj().T
    q = np.eye(len(v))-p
    hl, ol = compress(v, h), compress(v, observable)
    eps = norm(compress(v, h-boundary))
    boundary_leak = norm(q @ boundary @ v)
    observable_leak = norm(q @ observable @ v)
    h_leak = norm(q @ h @ v)
    locality = norm(compress(v, comm(boundary, observable)))
    rate = (locality+2*boundary_leak*observable_leak+2*eps*norm(ol))
    return {"compressed_energy_error": eps,
            "boundary_code_leak": boundary_leak,
            "observable_code_leak": observable_leak,
            "hamiltonian_code_leak": h_leak,
            "compressed_locality_defect": locality,
            "effective_commutator": norm(comm(hl, ol)),
            "effective_rate_bound": rate,
            "full_compressed_rate_bound": rate+2*norm(observable)*h_leak}


def compression_family(theta, gap=4.):
    """A commuting physical A/B pair whose arbitrary compression need not commute."""
    basis = np.eye(4, dtype=complex)
    v = np.column_stack(((basis[:, 0]-basis[:, 1])/np.sqrt(2),
                         np.cos(theta)*(basis[:, 0]+basis[:, 1])/np.sqrt(2) -
                         np.sin(theta)*basis[:, 2]))
    p = v @ v.conj().T
    q = np.eye(4)-p
    a = np.kron(Z, I2)  # Strict boundary support.
    b = np.kron(I2, Z)  # Strict interior support.
    h = p @ a @ p+gap*q  # Invariant low-energy code, but not a boundary operator.
    return v, h, a, b


def repetition_code(count):
    v = np.zeros((2**count, 2), dtype=complex)
    v[0, 0] = v[-1, 1] = 1.
    return v


def repetition_model(count, gap=3.):
    v = repetition_code(count)
    boundary = .5*at_site(Z, 0, count)
    penalty = np.zeros_like(boundary)
    for i in range(count-1):
        penalty += .5*gap*(np.eye(2**count) -
                          at_site(Z, i, count) @ at_site(Z, i+1, count))
    h = boundary+penalty
    dressed_x = kron_all([X]*count)
    return v, h, boundary, dressed_x


def span_rank(matrices, tolerance=1e-10):
    columns = np.column_stack([matrix.reshape(-1) for matrix in matrices])
    return int(np.linalg.matrix_rank(columns, tol=tolerance))


def interior_compression_rank(count):
    v = repetition_code(count)
    matrices = [compress(v, kron_all((I2,)+tuple(PAULI[j] for j in labels)))
                for labels in itertools.product(range(4), repeat=count-1)]
    return span_rank(matrices)


def random_hermitian(rng, dimension):
    raw = rng.normal(size=(dimension, dimension))+1j*rng.normal(size=(dimension, dimension))
    return (raw+raw.conj().T)/(2*dimension)


def report():
    theta_star = float(np.arccos(1/np.sqrt(3)))
    family = []
    for theta in (.15, .5, theta_star, 1.4):
        v, h, a, b = compression_family(theta)
        p, q = v @ v.conj().T, np.eye(4)-v @ v.conj().T
        al, bl = compress(v, a), compress(v, b)
        identity = (comm(al, bl)-compress(v, comm(a, b)) +
                    v.conj().T @ a @ q @ b @ v-v.conj().T @ b @ q @ a @ v)
        time = .9
        effective_change = norm(evolve_operator(al, bl, time)-bl)
        # Full boundary evolution leaves b exactly fixed. Compression alone
        # produces a different, unjustified Heisenberg evolution.
        false_prediction_error = norm(evolve_operator(al, bl, time) -
                                      compress(v, evolve_operator(a, b, time)))
        data = defects(v, h, a, b)
        family.append({
            "theta": theta,
            "underlying_boundary_bulk_commutator": norm(comm(a, b)),
            "compression_identity_residual": norm(identity),
            "analytic_effective_commutator": float(2*np.sin(theta)**2*np.cos(theta)),
            "analytic_boundary_leak": float(abs(np.sin(2*theta))),
            "analytic_observable_leak": float(np.sin(theta)*np.sqrt(1+np.cos(theta)**2)),
            "energy_band": np.linalg.eigvalsh(h).tolist(),
            "effective_operator_change_t_0_9": effective_change,
            "wrong_compression_prediction_error_t_0_9": false_prediction_error,
            "observable_square_compression_defect": norm(compress(v, b @ b)-bl @ bl),
            "strong_on_code_energy_difference": norm((h-a) @ v),
            **data})
    encoded = []
    for count in (2, 3, 5):
        v, h, boundary, dressed = repetition_model(count)
        time = .8
        hl = compress(v, h)
        encoded.append({
            "sites": count,
            "boundary_h_preserves_code": norm((np.eye(2**count)-v @ v.conj().T) @ boundary @ v),
            "energy_representative_on_code": norm((h-boundary) @ v),
            "all_strict_interior_compressions_rank": interior_compression_rank(count),
            "full_logical_algebra_rank": span_rank([I2, compress(v, dressed),
                                                   compress(v, boundary)*2,
                                                   compress(v, -2j*boundary @ dressed)]),
            "logical_intertwining_error_t_0_8": norm(unitary(h, time) @ v-v @ unitary(hl, time)),
            "dressed_X_operator_change_t_0_8": norm(evolve_operator(hl, X, time)-X),
            "dressed_X_boundary_commutator": norm(comm(boundary, dressed)),
            "minimum_interior_reconstruction_error_for_X": 1.,
            "excited_sector_gap_above_logical_band": 2.})
    slow = []
    boundary, bulk = np.kron(.5*Z, I2), np.kron(I2, Z)
    for epsilon in (.2, .04, .008):
        h = boundary+.5*epsilon*np.kron(I2, X)
        short = norm(evolve_operator(h, bulk, 1.)-bulk)
        long = norm(evolve_operator(h, bulk, 1/epsilon)-bulk)
        slow.append({"epsilon": epsilon,
                     "global_energy_remainder": norm(h-boundary),
                     "fixed_time_operator_change": short,
                     "fixed_time_bound": epsilon,
                     "time_1_over_epsilon_operator_change": long,
                     "analytic_long_time_value": float(2*np.sin(.5))})
    return {
        "round": 363,
        "scope": ("Finite-dimensional audit of Marolf's conditional boundary-energy "
                  "obstruction: common time, physical observable algebra, and a "
                  "compatible boundary/interior dictionary are required. The "
                  "models are algebraic witnesses, not gravitational theories or "
                  "a no-go theorem for the full cognitive program. Projection "
                  "compression need not preserve products or commutators. "
                  "Approximate bounds are for bounded operators and fixed "
                  "isometric codes, and are independently derived here."),
        "sources": ["https://arxiv.org/abs/1409.2509",
                    "https://link.aps.org/accepted/10.1103/PhysRevLett.114.031104"],
        "compression_counterexamples": family,
        "encoded_boundary_dynamics": encoded,
        "nonuniform_long_time_limit": slow,
        "general_bound": ("||[h,o]|| <= c+2 beta_HB beta_O+2 epsilon ||o||; "
                          "compressed true evolution differs from o by at most "
                          "|t| times this rate plus 2|t| ||O|| leakage_H. "
                          "All constants refer to explicitly defined operator norms."),
        "not_derived": ["choice of constraints or low-energy code",
                        "universal gravitational boundary energy law",
                        "physical locality dictionary and common clock",
                        "Einstein equations, Lorentz geometry, or holography"]
    }


class Checks(unittest.TestCase):
    def test_01_boundary_h_freezes_entire_interior_algebra(self):
        boundary_h = .7*np.kron(X, Z)+.2*np.kron(Z, X)+.3*np.kron(Y, Y)
        h = np.kron(boundary_h, np.eye(4))
        for labels in itertools.product(range(4), repeat=2):
            observable = np.kron(np.eye(4), kron_all([PAULI[j] for j in labels]))
            self.assertLess(norm(comm(h, observable)), 1e-14)
            for time in (.1, .7, 2.):
                self.assertLess(norm(evolve_operator(h, observable, time)-observable), 5e-14)

    def test_02_nonlocal_h_does_not_destroy_equal_time_commutativity(self):
        h = np.kron(X, X)+.6*np.kron(Z, Y)
        a, b = np.kron(Z, I2), np.kron(I2, Z)
        self.assertGreater(norm(comm(h, a)), 1.)
        for time in (.2, .8, 1.3):
            evolved_a, evolved_b = evolve_operator(h, a, time), evolve_operator(h, b, time)
            self.assertLess(norm(comm(evolved_a, evolved_b)), 4e-14)
        self.assertGreater(norm(comm(evolve_operator(h, a, .2), b)), .1)

    def test_03_compression_defect_identity_for_generic_noncommuting_operators(self):
        rng = np.random.default_rng(363)
        for full_dim, code_dim in ((4, 2), (7, 3), (9, 4)):
            raw = rng.normal(size=(full_dim, code_dim))+1j*rng.normal(size=(full_dim, code_dim))
            v = np.linalg.qr(raw)[0]
            q = np.eye(full_dim)-v @ v.conj().T
            a, b = random_hermitian(rng, full_dim), random_hermitian(rng, full_dim)
            lhs = comm(compress(v, a), compress(v, b))
            rhs = (compress(v, comm(a, b))-v.conj().T @ a @ q @ b @ v +
                   v.conj().T @ b @ q @ a @ v)
            np.testing.assert_allclose(lhs, rhs, atol=4e-16)
            self.assertLessEqual(norm(lhs-compress(v, comm(a, b))),
                                 2*norm(q @ a @ v)*norm(q @ b @ v)+1e-14)

    def test_04_commuting_boundary_bulk_compressions_need_not_commute(self):
        for theta in (.15, .5, np.arccos(1/np.sqrt(3)), 1.4):
            v, _, a, b = compression_family(theta)
            self.assertEqual(norm(comm(a, b)), 0.)
            expected_a = np.diag([1., np.cos(2*theta)])
            expected_b = np.array([[0., np.cos(theta)], [np.cos(theta), np.sin(theta)**2]])
            np.testing.assert_allclose(compress(v, a), expected_a, atol=5e-16)
            np.testing.assert_allclose(compress(v, b), expected_b, atol=5e-16)
            self.assertAlmostEqual(norm(comm(compress(v, a), compress(v, b))),
                                   2*np.sin(theta)**2*np.cos(theta), places=14)

    def test_05_exact_weak_energy_match_is_not_physical_boundary_operator(self):
        for theta in (.2, .8, 1.2):
            v, h, a, b = compression_family(theta)
            d = defects(v, h, a, b)
            self.assertLess(d["hamiltonian_code_leak"], 3e-15)
            self.assertLess(d["compressed_energy_error"], 3e-15)
            self.assertGreater(d["boundary_code_leak"], .3)
            self.assertGreater(d["observable_code_leak"], .2)
            self.assertAlmostEqual(d["boundary_code_leak"], abs(np.sin(2*theta)), places=14)
            self.assertAlmostEqual(d["observable_code_leak"],
                                   np.sin(theta)*np.sqrt(1+np.cos(theta)**2), places=14)
            self.assertAlmostEqual(norm((h-a) @ v), d["boundary_code_leak"], places=14)
            np.testing.assert_allclose(np.linalg.eigvalsh(h),
                                       [np.cos(2*theta), 1., 4., 4.], atol=3e-15)

    def test_06_compressed_fictitious_dynamics_disagrees_with_true_boundary_evolution(self):
        theta = float(np.arccos(1/np.sqrt(3)))
        v, h, a, b = compression_family(theta)
        for time in (.2, .9, 3*np.pi/4):
            compressed_true_a = compress(v, evolve_operator(a, b, time))
            logical_prediction = evolve_operator(compress(v, a), compress(v, b), time)
            np.testing.assert_allclose(compressed_true_a, compress(v, b), atol=5e-16)
            expected = 2*np.cos(theta)*abs(np.sin(np.sin(theta)**2*time))
            self.assertAlmostEqual(norm(logical_prediction-compressed_true_a), expected, places=14)
            # It IS valid for the invariant-code H, whose full-space operator is not A.
            np.testing.assert_allclose(compress(v, evolve_operator(h, b, time)),
                                       logical_prediction, atol=5e-15)

    def test_07_general_effective_and_full_time_bounds(self):
        rng = np.random.default_rng(1363)
        for _ in range(7):
            dimension = 6
            raw = rng.normal(size=(dimension, 3))+1j*rng.normal(size=(dimension, 3))
            v = np.linalg.qr(raw)[0]
            h, boundary, observable = [random_hermitian(rng, dimension) for _ in range(3)]
            d = defects(v, h, boundary, observable)
            self.assertLessEqual(d["effective_commutator"], d["effective_rate_bound"]+2e-14)
            for time in (.05, .3, 1.1):
                actual = norm(compress(v, evolve_operator(h, observable, time)) -
                              compress(v, observable))
                self.assertLessEqual(actual, time*d["full_compressed_rate_bound"]+2e-14)
                logical = evolve_operator(compress(v, h), compress(v, observable), time)
                self.assertLessEqual(norm(logical-compress(v, observable)),
                                     time*d["effective_rate_bound"]+2e-14)

    def test_08_code_dynamics_leakage_bound_and_reference(self):
        v, _, h, observable = compression_family(.8)
        q = np.eye(4)-v @ v.conj().T
        leakage = norm(q @ h @ v)
        for time in (.03, .3, 1.2):
            exact = unitary(h, time) @ v
            approximate = v @ unitary(compress(v, h), time)
            error = norm(exact-approximate)
            self.assertLessEqual(error, time*leakage+2e-15)
            self.assertAlmostEqual(norm(np.kron(exact-approximate, np.eye(3))), error, places=14)
            self.assertLessEqual(norm(q @ exact)**2, time*time*leakage*leakage+2e-15)
            heisenberg_error = norm(compress(v, evolve_operator(h, observable, time)) -
                                    evolve_operator(compress(v, h), compress(v, observable), time))
            self.assertLessEqual(heisenberg_error, 2*time*norm(observable)*leakage+2e-14)

    def test_09_code_preservation_prevents_commutator_artifact(self):
        rng = np.random.default_rng(2363)
        v = repetition_code(3)
        p = v @ v.conj().T
        q = np.eye(8)-p
        a0 = random_hermitian(rng, 8)
        a = p @ a0 @ p+q @ a0 @ q
        b = random_hermitian(rng, 8)
        np.testing.assert_allclose(comm(compress(v, a), compress(v, b)),
                                   compress(v, comm(a, b)), atol=2e-16)

    def test_10_true_boundary_logical_dynamics_preserves_unknown_reference(self):
        bell = np.array([1., 0., 0., 1.], dtype=complex)/np.sqrt(2)
        for count in (2, 3, 5):
            v, h, boundary, dressed = repetition_model(count)
            self.assertLess(norm((h-boundary) @ v), 1e-15)
            self.assertLess(norm((np.eye(2**count)-v @ v.conj().T) @ boundary @ v), 1e-15)
            for time in (.1, .8, 1.7):
                lhs = np.kron(unitary(h, time) @ v, I2) @ bell
                rhs = np.kron(v @ unitary(.5*Z, time), I2) @ bell
                np.testing.assert_allclose(lhs, rhs, atol=2e-15)
                np.testing.assert_allclose(compress(v, evolve_operator(h, dressed, time)),
                                           np.cos(time)*X-np.sin(time)*Y, atol=3e-15)

    def test_11_all_strict_interior_compressions_are_diagonal(self):
        for count in (2, 3, 5):
            self.assertEqual(interior_compression_rank(count), 2)
            v = repetition_code(count)
            rng = np.random.default_rng(363+count)
            interior = random_hermitian(rng, 2**(count-1))
            compressed = compress(v, np.kron(I2, interior))
            self.assertEqual(compressed[0, 1], 0.)
            self.assertEqual(compressed[1, 0], 0.)
            # Column action on |0> proves distance >= 1 for any diagonal.
            self.assertGreaterEqual(norm(X-compressed), 1.-1e-15)
        self.assertEqual(norm(X), 1.)  # Zero diagonal attains the lower bound.

    def test_12_dynamical_logical_phase_requires_boundary_support(self):
        for count in (2, 4, 5):
            v, h, boundary, dressed = repetition_model(count)
            self.assertAlmostEqual(norm(comm(boundary, dressed)), 1.)
            for site in range(1, count):
                local_z = at_site(Z, site, count)
                self.assertEqual(norm(comm(boundary, local_z)), 0.)
                self.assertEqual(norm(comm(compress(v, h), compress(v, local_z))), 0.)
            plus = np.array([1., 1.])/np.sqrt(2)
            time = .8
            expectation = np.vdot(plus, evolve_operator(compress(v, h), X, time) @ plus).real
            self.assertAlmostEqual(float(expectation), float(np.cos(time)), places=14)

    def test_13_repetition_code_is_an_isolated_low_energy_sector(self):
        for count in (2, 3, 5):
            v, h, _, _ = repetition_model(count)
            spectrum = np.linalg.eigvalsh(h)
            np.testing.assert_allclose(spectrum[:2], [-.5, .5], atol=1e-15)
            self.assertAlmostEqual(float(spectrum[2]-spectrum[1]), 2.)
            p = v @ v.conj().T
            self.assertLess(norm(comm(h, p)), 1e-15)

    def test_14_small_nonboundary_remainder_only_freezes_fixed_time_windows(self):
        boundary, observable = np.kron(.5*Z, I2), np.kron(I2, Z)
        for epsilon in (.2, .04, .008):
            h = boundary+.5*epsilon*np.kron(I2, X)
            self.assertAlmostEqual(norm(h-boundary), epsilon/2)
            for time in (.2, 1., 1/epsilon):
                deviation = norm(evolve_operator(h, observable, time)-observable)
                self.assertAlmostEqual(deviation, 2*abs(np.sin(epsilon*time/2)), places=13)
                self.assertLessEqual(deviation, epsilon*time+2e-14)
            self.assertGreater(norm(evolve_operator(h, observable, 1/epsilon)-observable), .95)

    def test_15_compression_preserves_first_moments_but_not_measurement_products(self):
        for theta in (.2, .8, np.arccos(1/np.sqrt(3))):
            v, _, _, b = compression_family(theta)
            p = v @ v.conj().T
            q = np.eye(4)-p
            bl = compress(v, b)
            difference = compress(v, b @ b)-bl @ bl
            np.testing.assert_allclose(difference, v.conj().T @ b @ q @ b @ v, atol=7e-16)
            self.assertGreaterEqual(float(np.min(np.linalg.eigvalsh(difference))), -1e-15)
            self.assertAlmostEqual(norm(difference), 1-np.cos(theta)**4, places=14)
            values, vectors = np.linalg.eigh(bl)
            state = vectors[:, 0]
            self.assertAlmostEqual(float(values[0]), -np.cos(theta)**2, places=14)
            first_full = np.vdot(v @ state, b @ v @ state).real
            self.assertAlmostEqual(float(first_full), float(np.vdot(state, bl @ state).real))
            self.assertAlmostEqual(float(np.vdot(v @ state, b @ b @ v @ state).real), 1.)
            self.assertAlmostEqual(float(np.vdot(state, bl @ bl @ state).real),
                                   np.cos(theta)**4, places=14)
            physical_extension = p @ b @ p+q @ b @ q
            self.assertLess(norm(comm(physical_extension, p)), 2e-15)
            self.assertGreater(norm(comm(physical_extension, np.kron(X, I2))), .05)


if __name__ == "__main__":
    main(__name__, "emergent_boundary_dynamics_audit", report)
