"""Round 359: a finite representation-theoretic Weinberg-Witten audit.

The rotation selection lemma is exact. Approximate-covariance and TT-gauge
examples are kinematic tensors, not an interacting graviton construction.
"""
import unittest
import numpy as np
from growing_stream_audit import main


ETA = np.diag([-1., 1., 1., 1.])


def rotation(theta):
    matrix = np.eye(4)
    c, s = np.cos(theta), np.sin(theta)
    matrix[1:3, 1:3] = [[c, -s], [s, c]]
    return matrix


def rotation_generator():
    generator = np.zeros((4, 4))
    generator[1, 2], generator[2, 1] = -1., 1.
    return generator


def symmetric_basis():
    basis = []
    for i in range(4):
        for j in range(i, 4):
            element = np.zeros((4, 4))
            element[i, j] = 1.
            if i != j:
                element[j, i] = 1.
                element /= np.sqrt(2)
            basis.append(element)
    return np.array(basis)


BASIS = symmetric_basis()


def representation(theta, rank=2, symmetric=True):
    r = rotation(theta)
    if rank == 1:
        return r
    if not symmetric:
        return np.kron(r, r)
    transformed = np.einsum("ij,bjk,lk->bil", r, BASIS, r)
    return np.einsum("aij,bij->ab", BASIS, transformed)


def hermitian_generator(rank=2, symmetric=True):
    generator = rotation_generator()
    if rank == 1:
        return -1j*generator
    if not symmetric:
        return -1j*(np.kron(generator, np.eye(4))+
                    np.kron(np.eye(4), generator))
    transformed = np.array([generator @ b+b @ generator.T for b in BASIS])
    return -1j*np.einsum("aij,bij->ab", BASIS, transformed)


def rotation_projector(weight, rank=2, symmetric=True, samples=32):
    """Exact trapezoidal integral for these bounded Fourier polynomials."""
    dimension = 4 if rank == 1 else (10 if symmetric else 16)
    out = np.zeros((dimension, dimension), dtype=complex)
    for theta in 2*np.pi*np.arange(samples)/samples:
        out += np.exp(-1j*weight*theta)*representation(theta, rank, symmetric)
    return out/samples


def kernel_dimension(helicity, rank=2, symmetric=True):
    generator = hermitian_generator(rank, symmetric)
    singular = np.linalg.svd(generator-2*helicity*np.eye(len(generator)),
                             compute_uv=False)
    return int(np.count_nonzero(singular < 1e-10))


def tensor_coordinates(tensor):
    return np.einsum("aij,ij->a", BASIS, tensor)


def covariance_rms(tensor, helicity=2., samples=32):
    values = []
    for theta in 2*np.pi*np.arange(samples)/samples:
        r = rotation(theta)
        defect = r @ tensor @ r.T-np.exp(2j*helicity*theta)*tensor
        values.append(np.linalg.norm(defect)**2)
    return float(np.sqrt(np.mean(values)))


def discrete_defects(tensor, helicity=2., order=4):
    return [float(np.linalg.norm(
        rotation(theta) @ tensor @ rotation(theta).T -
        np.exp(2j*helicity*theta)*tensor))
        for theta in 2*np.pi*np.arange(order)/order]


def boost(axis, beta):
    gamma = 1/np.sqrt(1-beta*beta)
    matrix = np.eye(4)
    matrix[0, 0] = matrix[axis, axis] = gamma
    matrix[0, axis] = matrix[axis, 0] = -gamma*beta
    return matrix


def brick_wall_pair(alpha, energy=1.):
    incoming = energy*np.array([1., np.cos(alpha), 0., np.sin(alpha)])
    outgoing = energy*np.array([1., np.cos(alpha), 0., -np.sin(alpha)])
    transform = boost(1, np.cos(alpha))
    return incoming, outgoing, transform


def nonuniform_example(alpha):
    """A kinematic sequence, not a physical stress form factor."""
    _, _, transform = brick_wall_pair(alpha)
    null_direction = np.array([1., 1., 0., 0.])/np.sqrt(2)
    laboratory = np.outer(null_direction, null_direction)
    brick = transform @ laboratory @ transform.T
    recovered = np.linalg.inv(transform) @ brick @ np.linalg.inv(transform).T
    return {
        "alpha": alpha,
        "gamma": 1/np.sin(alpha),
        "rank2_boost_amplification": 1/np.tan(alpha/2)**2,
        "laboratory_Frobenius_norm": float(np.linalg.norm(laboratory)),
        "brick_Frobenius_norm": float(np.linalg.norm(brick)),
        "predicted_brick_norm": float(np.tan(alpha/2)**2),
        "brick_rotation_RMS_defect": covariance_rms(brick),
        "back_transformed_norm": float(np.linalg.norm(recovered)),
    }


def polarization():
    e = np.array([0., 1., 1j, 0.])/np.sqrt(2)
    return np.outer(e, e)


def restore_temporal_gauge(momentum, tensor):
    """Contravariant h and p; h -> h+p*xi+xi*p sets h^(0,mu)=0."""
    xi = np.zeros(4, dtype=complex)
    xi[0] = -tensor[0, 0]/(2*momentum[0])
    xi[1:] = -(tensor[0, 1:]+momentum[1:]*xi[0])/momentum[0]
    gauge = np.outer(momentum, xi)+np.outer(xi, momentum)
    return tensor+gauge, xi, gauge


def linear_curvature(momentum, tensor):
    """Fourier linearized Riemann, overall Fourier sign is immaterial here."""
    p = ETA @ momentum
    h = ETA @ tensor @ ETA
    curvature = np.zeros((4, 4, 4, 4), dtype=complex)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for ell in range(4):
                    curvature[i, j, k, ell] = .5*(
                        p[j]*p[k]*h[i, ell]+p[i]*p[ell]*h[j, k] -
                        p[i]*p[k]*h[j, ell]-p[j]*p[ell]*h[i, k])
    return curvature


def report():
    generator = hermitian_generator()
    rows = [{"helicity": h, "required_rotation_weight": 2*h,
             "current_rank1_kernel": kernel_dimension(h, 1),
             "symmetric_rank2_kernel": kernel_dimension(h),
             "full_rank2_kernel": kernel_dimension(h, symmetric=False)}
            for h in (0., .5, 1., 1.5, 2., -2.)]
    candidate = np.diag([1., 0., 0., 0.])
    r = rotation(np.pi/4)
    p = np.array([1., 0., 0., 1.])
    transform = boost(1, .6)
    p_new = transform @ p
    raw = transform @ polarization() @ transform.T
    restored, xi, gauge = restore_temporal_gauge(p_new, raw)
    return {
        "round": 359,
        "scope": ("Exact finite rotation-weight lemma underlying the Weinberg-Witten "
                  "obstruction, discrete-rotation and nonuniform-boost kinematic "
                  "counterexamples to weakened arguments, and TT gauge restoration. "
                  "No finite calculation proves all QFT assumptions, constructs a "
                  "stress operator, or excludes all emergent/curved-spacetime gravity."),
        "selection_table": rows,
        "symmetric_rank2_generator_spectrum": np.linalg.eigvalsh(generator).tolist(),
        "h2_generator_smallest_singular_value": float(np.linalg.svd(
            generator-4*np.eye(10), compute_uv=False)[-1]),
        "h2_continuous_projector_norm": float(np.linalg.norm(rotation_projector(4))),
        "C4_aliasing": {
            "candidate": "F^(00)=1; all other components zero",
            "four_rotation_defects": discrete_defects(candidate),
            "pi_over_4_defect": float(np.linalg.norm(r @ candidate @ r.T+candidate)),
            "continuous_rotation_RMS_defect": covariance_rms(candidate),
            "predicted_RMS": float(np.sqrt(2)),
        },
        "nonuniform_boost_sequence": [nonuniform_example(a) for a in (.4, .2, .1, .05)],
        "TT_gauge_example": {
            "boost_beta": .6, "boosted_momentum": p_new.tolist(),
            "raw_time_row_norm": float(np.linalg.norm(raw[0])),
            "restored_time_row_norm": float(np.linalg.norm(restored[0])),
            "gauge_change_norm": float(np.linalg.norm(gauge)),
            "curvature_change_norm": float(np.linalg.norm(
                linear_curvature(p_new, restored)-linear_curvature(p_new, raw))),
            "curvature_norm": float(np.linalg.norm(linear_curvature(p_new, raw))),
            "p_dot_gauge_parameter_abs": float(abs(p_new @ ETA @ xi)),
        },
        "theorem_inputs": [
            "Exact appropriate 3+1 Poincare covariance on physical finite-helicity states",
            "An exactly massless one-particle helicity |h|>1 sector",
            "A conserved rank2 stress operator with Lorentz-covariant physical matrix elements",
            "Its space integral generates that particle's nonzero four-momentum",
            "The near-forward operational charge limit used in original footnote 3",
        ],
        "not_sufficient_inputs": [
            "A tensor-shaped response coefficient or curved acoustic metric",
            "Only C4 rotational symmetry or only low-energy approximate Lorentz dispersion",
            "Posthoc G/kappa definition of stress", "Finite CP instruments and Time alone",
        ],
        "sources": [
            "https://doi.org/10.1016/0370-2693(80)90212-9",
            "https://hep.physics.uoc.gr/Erasmus-IP-2013/files/lectures/chiral-symmetry/Weinberg-Witten.pdf",
        ],
    }


class Checks(unittest.TestCase):
    def test_01_rotations_and_boosts_preserve_the_Minkowski_metric(self):
        for angle in (.1, .7, 2.):
            r = rotation(angle)
            np.testing.assert_allclose(r.T @ ETA @ r, ETA, atol=4e-16)
            np.testing.assert_allclose(r @ rotation(-angle), np.eye(4), atol=4e-16)
        for axis in (1, 2, 3):
            b = boost(axis, .6)
            np.testing.assert_allclose(b.T @ ETA @ b, ETA, atol=4e-16)

    def test_02_symmetric_tensor_basis_and_generator_spectrum(self):
        np.testing.assert_allclose(np.einsum("aij,bij->ab", BASIS, BASIS), np.eye(10),
                                   atol=3e-16)
        j = hermitian_generator()
        np.testing.assert_allclose(j, j.conjugate().T, atol=1e-15)
        np.testing.assert_allclose(np.linalg.eigvalsh(j),
                                   [-2., -1., -1., 0., 0., 0., 0., 1., 1., 2.],
                                   atol=1e-14)

    def test_03_vector_current_and_rank2_thresholds(self):
        self.assertEqual([kernel_dimension(h, 1) for h in (0., .5, 1.)], [2, 1, 0])
        self.assertEqual([kernel_dimension(h) for h in (0., .5, 1., 1.5, 2.)],
                         [4, 2, 1, 0, 0])
        self.assertEqual(kernel_dimension(-2.), 0)

    def test_04_not_assuming_a_symmetric_stress_is_needed_for_the_obstruction(self):
        self.assertEqual(kernel_dimension(2., symmetric=False), 0)
        self.assertEqual(kernel_dimension(-2., symmetric=False), 0)
        self.assertAlmostEqual(np.max(np.linalg.eigvalsh(
            hermitian_generator(symmetric=False))), 2.)

    def test_05_all_allowed_Fourier_projectors_are_complete(self):
        projectors = [rotation_projector(m) for m in range(-2, 3)]
        np.testing.assert_allclose(sum(projectors), np.eye(10), atol=2e-15)
        for p in projectors:
            np.testing.assert_allclose(p @ p, p, atol=2e-15)
            np.testing.assert_allclose(p.conjugate().T, p, atol=2e-15)
        self.assertLess(np.linalg.norm(rotation_projector(4)), 2e-15)

    def test_06_compatible_helicity_one_rank2_is_a_nonzero_control(self):
        e = np.array([0., 1., -1j, 0.])/np.sqrt(2)
        tensor = np.outer(e, e)
        self.assertLess(covariance_rms(tensor, helicity=1.), 2e-15)
        self.assertAlmostEqual(np.linalg.norm(tensor), 1.)

    def test_07_every_nonzero_rank2_tensor_has_h2_RMS_defect(self):
        rng = np.random.default_rng(359)
        for _ in range(12):
            f = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
            self.assertAlmostEqual(covariance_rms(f), np.sqrt(2)*np.linalg.norm(f),
                                   places=12)

    def test_08_generator_defect_has_a_sharp_gap_two(self):
        j = hermitian_generator()
        values, vectors = np.linalg.eigh(j)
        f = vectors[:, np.argmax(values)]
        self.assertAlmostEqual(np.linalg.norm((j-4*np.eye(10)) @ f), 2.)
        rng = np.random.default_rng(360)
        for _ in range(8):
            v = rng.normal(size=10)+1j*rng.normal(size=10)
            self.assertGreaterEqual(np.linalg.norm((j-4*np.eye(10)) @ v)+1e-14,
                                    2*np.linalg.norm(v))

    def test_09_C4_covariance_does_not_imply_continuous_covariance(self):
        f = np.diag([1., 0., 0., 0.])
        self.assertLess(max(discrete_defects(f)), 1e-14)
        r = rotation(np.pi/4)
        self.assertAlmostEqual(np.linalg.norm(r @ f @ r.T+f), 2.)
        self.assertAlmostEqual(covariance_rms(f), np.sqrt(2))

    def test_10_noncollinear_massless_pair_has_a_brick_wall_frame(self):
        for alpha in (.4, .2, .1, .05):
            p, q, b = brick_wall_pair(alpha)
            self.assertAlmostEqual(p @ ETA @ p, 0., places=14)
            self.assertAlmostEqual(q @ ETA @ q, 0., places=14)
            self.assertGreater((p-q) @ ETA @ (p-q), 0.)
            np.testing.assert_allclose(b @ p,
                                       [np.sin(alpha), 0., 0., np.sin(alpha)], atol=3e-14)
            np.testing.assert_allclose(b @ q,
                                       [np.sin(alpha), 0., 0., -np.sin(alpha)], atol=3e-14)

    def test_11_rank2_boost_amplification_is_unbounded(self):
        for alpha in (.4, .2, .1, .05):
            _, _, b = brick_wall_pair(alpha)
            norm = np.linalg.norm(np.kron(b, b), 2)
            self.assertAlmostEqual(norm*np.tan(alpha/2)**2, 1., places=11)

    def test_12_small_brick_frame_errors_can_have_nonzero_laboratory_limit(self):
        for alpha in (.4, .2, .1, .05):
            row = nonuniform_example(alpha)
            self.assertAlmostEqual(row["brick_Frobenius_norm"],
                                   row["predicted_brick_norm"], places=12)
            self.assertAlmostEqual(row["laboratory_Frobenius_norm"], 1., places=14)
            self.assertAlmostEqual(row["back_transformed_norm"], 1., places=9)

    def test_13_physical_TT_representative_requires_gauge_restoration_after_boost(self):
        p = np.array([1., 0., 0., 1.])
        b = boost(1, .6)
        p_new = b @ p
        raw = b @ polarization() @ b.T
        restored, xi, _ = restore_temporal_gauge(p_new, raw)
        self.assertGreater(np.linalg.norm(raw[0]), .1)
        self.assertLess(np.linalg.norm(restored[0]), 2e-15)
        self.assertLess(np.linalg.norm(p_new @ ETA @ restored), 2e-15)
        self.assertLess(abs(np.einsum("ij,ij", ETA, restored)), 2e-15)
        self.assertLess(abs(p_new @ ETA @ xi), 2e-15)

    def test_14_gauge_restoration_changes_the_representative_not_its_curvature(self):
        b = boost(1, .6)
        p = b @ np.array([1., 0., 0., 1.])
        raw = b @ polarization() @ b.T
        restored, _, gauge = restore_temporal_gauge(p, raw)
        self.assertGreater(np.linalg.norm(gauge), .1)
        self.assertLess(np.linalg.norm(linear_curvature(p, gauge)), 2e-15)
        np.testing.assert_allclose(linear_curvature(p, raw),
                                   linear_curvature(p, restored), atol=2e-15)
        self.assertGreater(np.linalg.norm(linear_curvature(p, restored)), 1.)

    def test_15_forward_stress_normalization_carries_nonzero_momentum(self):
        for energy in (.3, 1., 2.):
            p = energy*np.array([1., .6, 0., .8])
            forward = 2*np.outer(p, p)
            np.testing.assert_allclose(forward[0]/(2*p[0]), p, atol=2e-15)
            self.assertGreater(forward[0, 0], 0.)


if __name__ == "__main__":
    main(__name__, "emergent_spin2_obstruction_audit", report)
