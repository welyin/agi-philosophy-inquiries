"""Round 383: free reversal and the sharp obstruction to a Hopf lift.

All continuous lifts R of Bloch antipodes on S^3 have maximal worst-case
two-step label error 2. The theorem is analytic, not inferred from sampling.
J acts on specified source labels, never as a channel on an unknown qubit.
"""
import unittest

import numpy as np

from growing_stream_audit import main
from direction_hopf_faithfulness_audit import (
    I2, PAULI, normalize, projector, bloch, su2, controlled,
    source, trace_distance, helstrom_instrument,
)


def j_reverse(z):
    z = np.asarray(z, complex)
    return np.array([-z[1].conjugate(), z[0].conjugate()])


def real_coordinates(z):
    return np.array([z[0].real, z[0].imag, z[1].real, z[1].imag])


def complex_coordinates(x):
    return np.array([x[0] + 1j * x[1], x[2] + 1j * x[3]])


def phase_reverse(z, phi):
    return np.exp(1j * phi(z)) * j_reverse(z)


def square_factor(z, phi):
    return -np.exp(1j * (phi(phase_reverse(z, phi)) - phi(z)))


def phase_on_bloch(n, kappa=0.9, nu=0.4):
    return float(kappa * n[2] + nu * n[0] * n[1])


def reversible_reverse(z, kappa=0.9, nu=0.4):
    return np.exp(1j * phase_on_bloch(bloch(z), kappa, nu)) * j_reverse(z)


def reversible_inverse(w, kappa=0.9, nu=0.4):
    return -np.exp(1j * phase_on_bloch(-bloch(w), kappa, nu)) * j_reverse(w)


def chord_involution(z, b):
    """Swap the endpoints of the S^3 chord through an interior point -b."""
    x, b = real_coordinates(z), np.asarray(b, float)
    if np.linalg.norm(b) >= 1:
        raise ValueError("The chord centre must be strictly inside the ball.")
    u = x + b
    return complex_coordinates(-b - (1 - np.dot(b, b)) * u / np.dot(u, u))


def linear_phase_witness(coefficients):
    """A computable IVT witness for a non-fibre-invariant continuous phase.

    The extrema of the real linear phase are known exactly. Bisection locates
    a zero of phi(Rz)-phi(z) on a great semicircle joining those extrema.
    """
    a = np.asarray(coefficients, float)
    start = -a / np.linalg.norm(a)
    basis = np.eye(4)[int(np.argmin(abs(start)))]
    tangent = basis - np.dot(basis, start) * start
    tangent /= np.linalg.norm(tangent)
    phi = lambda z: float(np.dot(a, real_coordinates(z)))

    def path(t):
        return complex_coordinates(np.cos(np.pi * t) * start
                                   + np.sin(np.pi * t) * tangent)

    def residual(t):
        z = path(t)
        return phi(phase_reverse(z, phi)) - phi(z)

    lo, hi = 0., 1.
    brackets = (residual(lo), residual(hi))
    for _ in range(60):
        mid = (lo + hi) / 2
        if residual(mid) >= 0:
            lo = mid
        else:
            hi = mid
    z = path((lo + hi) / 2)
    return z, phi, brackets, residual((lo + hi) / 2)


def discontinuous_section(n):
    """Northern section with numerical south-pole protection.

    The analytic cut assigns only the exact south pole separately. Floating-point
    code uses a 1e-14 cap; tests of the section identity avoid its punctured part.
    No exact identity is claimed for all near-pole inputs inside that cap.
    """
    n = np.asarray(n, float)
    if n[2] <= -1 + 1e-14:
        return np.array([0., 1.], complex)
    a = np.sqrt((1 + n[2]) / 2)
    return normalize([a, (n[0] + 1j * n[1]) / (2 * a)])


def discontinuous_involution(z):
    n = bloch(z)
    section = discontinuous_section(n)
    phase = np.vdot(section, z)
    phase /= abs(phase)
    return phase * discontinuous_section(-n)


def finite_program(z):
    """Four explicit orthogonal labels, each with its own known source gate.

    T=U(S x I)U^dagger is a single fixed finite unitary on program,C,Q.
    No universally programmable continuum gate or unknown-state spin flip.
    """
    labels = [np.asarray(z, complex)]
    for _ in range(3):
        labels.append(j_reverse(labels[-1]))
    shift = np.roll(np.eye(4), 1, axis=0)
    u = np.zeros((16, 16), complex)
    for index, label in enumerate(labels):
        u[4 * index:4 * index + 4, 4 * index:4 * index + 4] = controlled(su2(label))
    translated = np.kron(shift, np.eye(4))
    t = u @ translated @ u.conj().T
    return labels, shift, u, t


def random_labels(rng, count):
    a = rng.normal(size=(count, 2)) + 1j * rng.normal(size=(count, 2))
    return a / np.linalg.norm(a, axis=1)[:, None]


def report():
    coefficients = np.array([0.8, -0.3, 0.6, 0.2])
    witness, phi, bracket, residual = linear_phase_witness(coefficients)
    twice = phase_reverse(phase_reverse(witness, phi), phi)
    source_errors = {str(v): trace_distance(source(twice, v), source(witness, v))
                     for v in (0., 0.25, 0.6, 1.)}
    north = np.array([1., 0.], complex)
    equator = normalize([1., 1.])
    kappa, nu = np.pi / 2, 0.4
    labels, _, _, t = finite_program(equator)
    angle = 0.12
    gate = np.diag([np.exp(-0.5j * angle), np.exp(0.5j * angle)])
    eps = float(np.linalg.norm(gate - I2, 2))
    return {
        "round": 383,
        "scope": {
            "positive_interface": "On a complete sphere boundary a continuous fixed-point-free involution suffices for the fixed-trace qubit d<=3 upper bound; no chosen conjugacy to the standard antipodal map or smooth/cellular action is added.",
            "sharp_new_result": "Every continuous R:S^3->S^3 with P_(Rz)=I-P_z satisfies sup_z ||R(Rz)-z||=2, even without assuming R is bijective. Nevertheless P_(R(Rz))=P_z everywhere.",
            "dual_sharp_result": "Every continuous involution on S^3 has worst-case operator-norm error 1 relative to exact projector complementation, whether or not it is free.",
            "source_application": "For the already specified round-381 coherent source of visibility v, the worst two-step source trace distance is exactly v. This is not a new recovery optimum v/2.",
            "inputs": ["A complete sphere boundary or the explicit S^3 source-label model", "Actual endpoint identity stronger than one family of effect readings", "Continuity on the entire label space", "Calibrated source reference for operationally reading fibre phase"],
            "not_proved": ["Sphere boundaries or endpoint reversal from FUCP", "Extra spatial dimensions from an observable source phase", "A CPTP spin flip of an arbitrary unknown qubit", "Uniform behaviour of all reversers from finite sampling", "GR or a complete physical three-dimensionality theorem"],
        },
        "general_continuous_phase_witness": {
            "linear_phase_coefficients": coefficients.tolist(),
            "phase_difference_at_extrema": list(bracket),
            "bisected_phase_difference": residual,
            "two_step_label_error": float(np.linalg.norm(twice - witness)),
            "two_step_effect_error": float(np.linalg.norm(projector(twice) - projector(witness), 2)),
            "source_errors_by_visibility": source_errors,
        },
        "reversible_nonlinear_phase_family": {
            "phase": "kappa*n_z+nu*n_x*n_y",
            "kappa": float(kappa), "nu": nu,
            "square_factor": "-exp(-2i*kappa*n_z)",
            "north_pole_two_step_error": float(np.linalg.norm(reversible_reverse(reversible_reverse(north, kappa, nu), kappa, nu) - north)),
            "equator_two_step_error": float(np.linalg.norm(reversible_reverse(reversible_reverse(equator, kappa, nu), kappa, nu) - equator)),
            "analytic_global_worst_error": 2.,
            "inverse_is_explicit": True,
        },
        "finite_legal_program": {
            "orthogonal_program_labels": 4,
            "program_C_Q_dimension": 16,
            "twice_returns_effect_not_program": True,
            "square_distance_from_identity": float(np.linalg.norm(t @ t - np.eye(16), 2)),
            "fourth_power_residual": float(np.linalg.norm(np.linalg.matrix_power(t, 4) - np.eye(16), 2)),
            "reference_scope": "U^dagger T^k U=S^k x I, tensored with identity on any reference; logical payload is unchanged after decoding while the program shifts.",
        },
        "nonstandard_free_involution": {
            "construction": "Chord endpoint exchange through -b with ||b||<1; a continuous free involution usually different from the coordinate antipode.",
            "b": [0.2, -0.1, 0.15, 0.05],
            "complement_error_at_a_label_parallel_to_b": 1.,
            "scope": "An explicit S^3 label model, not an inferred physical position space.",
        },
        "perturbed_reversal_stability": {
            "implemented_map": "z -> exp(-i*angle*Z/2) Jz",
            "angle": angle,
            "uniform_one_step_label_error": eps,
            "lower_bound_at_sharp_witness_2_minus_2epsilon": 2 - 2 * eps,
            "exact_two_step_error": float(2 * np.cos(angle / 2)),
            "scope": "This perturbation need not exactly reverse the effect. The bound is conditional on a uniform implementation error relative to an exact lift.",
        },
        "primary_sources": [
            "https://arxiv.org/pdf/1008.1134",
            "https://pi.math.cornell.edu/~hatcher/AT/ATch1.pdf",
            "https://arxiv.org/pdf/quant-ph/0108137",
            "https://arxiv.org/pdf/quant-ph/9901053",
        ],
    }


class Checks(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(383)

    def test_01_standard_lift_is_fourth_order_and_reverses_the_projector(self):
        for z in random_labels(self.rng, 20):
            np.testing.assert_allclose(j_reverse(j_reverse(z)), -z)
            np.testing.assert_allclose(projector(j_reverse(z)), I2 - projector(z), atol=4e-16)
            self.assertAlmostEqual(np.vdot(z, j_reverse(z)), 0)
            np.testing.assert_allclose(j_reverse(1j * z), -1j * j_reverse(z))

    def test_02_phase_of_an_exact_lift_is_unique_and_recovers_the_whole_vector(self):
        phi = lambda z: float(0.7 * z[0].real - 0.3 * z[1].imag + bloch(z)[2] ** 2)
        for z in random_labels(self.rng, 35):
            w = phase_reverse(z, phi)
            coefficient = np.vdot(j_reverse(z), w)
            self.assertAlmostEqual(abs(coefficient), 1)
            np.testing.assert_allclose(w, coefficient * j_reverse(z), atol=3e-16)
            np.testing.assert_allclose(projector(w), I2 - projector(z), atol=5e-16)

    def test_03_arbitrary_continuous_phase_square_matches_cocycle_formula(self):
        phi = lambda z: float(1.4 * z[0].real + 0.6 * z[1].imag + 0.3 * bloch(z)[0])
        for z in random_labels(self.rng, 40):
            twice = phase_reverse(phase_reverse(z, phi), phi)
            np.testing.assert_allclose(twice, square_factor(z, phi) * z, atol=1e-15)
            np.testing.assert_allclose(projector(twice), projector(z), atol=1e-15)

    def test_04_extremal_phase_argument_has_a_computable_sharp_witness(self):
        for coefficients in ([0.8, -0.3, 0.6, 0.2], [1., 0.2, -0.1, 0.7],
                             [-0.2, 0.9, 0.6, -0.4]):
            z, phi, bracket, residual = linear_phase_witness(coefficients)
            self.assertGreaterEqual(bracket[0], 0)
            self.assertLessEqual(bracket[1], 0)
            self.assertLess(abs(residual), 2e-15)
            twice = phase_reverse(phase_reverse(z, phi), phi)
            np.testing.assert_allclose(twice, -z, atol=2e-15)
            self.assertAlmostEqual(np.linalg.norm(twice - z), 2)

    def test_05_nonlinear_reversible_family_has_a_two_sided_inverse(self):
        for kappa, nu in ((0.9, 0.4), (np.pi / 2, 0.4), (3.2, -1.7)):
            for z in random_labels(self.rng, 25):
                np.testing.assert_allclose(reversible_inverse(reversible_reverse(z, kappa, nu), kappa, nu), z, atol=4e-15)
                np.testing.assert_allclose(reversible_reverse(reversible_inverse(z, kappa, nu), kappa, nu), z, atol=4e-15)

    def test_06_exact_return_at_poles_does_not_extend_to_the_whole_boundary(self):
        kappa, nu = np.pi / 2, 0.4
        for z in random_labels(self.rng, 30):
            twice = reversible_reverse(reversible_reverse(z, kappa, nu), kappa, nu)
            np.testing.assert_allclose(twice, -np.exp(-2j * kappa * bloch(z)[2]) * z, atol=2e-15)
            self.assertAlmostEqual(np.linalg.norm(twice - z), 2 * abs(np.cos(kappa * bloch(z)[2])))
        for z in (np.array([1., 0.]), np.array([0., 1.])):
            np.testing.assert_allclose(reversible_reverse(reversible_reverse(z, kappa, nu), kappa, nu), z, atol=2e-15)
        for phase in (0., 0.7, 2.1):
            z = normalize([1., np.exp(1j * phase)])
            self.assertAlmostEqual(np.linalg.norm(reversible_reverse(reversible_reverse(z, kappa, nu), kappa, nu) - z), 2)

    def test_07_source_visibility_translates_the_sharp_cycle_obstruction(self):
        z, phi, _, _ = linear_phase_witness([0.8, -0.3, 0.6, 0.2])
        twice = phase_reverse(phase_reverse(z, phi), phi)
        for visibility in (0., 0.25, 0.6, 1.):
            self.assertAlmostEqual(trace_distance(source(twice, visibility), source(z, visibility)), visibility)
        for z in random_labels(self.rng, 25):
            twice = reversible_reverse(reversible_reverse(z))
            self.assertAlmostEqual(trace_distance(source(twice, 0.6), source(z, 0.6)),
                                   0.3 * np.linalg.norm(twice - z))

    def test_08_full_effect_return_coexists_with_a_deterministic_source_instrument_witness(self):
        z = normalize([0.4 + 0.2j, 0.7 - 0.1j])
        twice = j_reverse(j_reverse(z))
        for state in [projector(w) for w in random_labels(self.rng, 15)]:
            self.assertAlmostEqual(np.trace(projector(twice) @ state), np.trace(projector(z) @ state))
        a, b = source(twice, 0.6), source(z, 0.6)
        e, kraus = helstrom_instrument(a, b)
        np.testing.assert_allclose(sum(k.conj().T @ k for k in kraus), np.eye(4), atol=2e-15)
        self.assertAlmostEqual(abs(np.trace(e @ (a - b))), 0.6)

    def test_09_finite_program_is_a_legal_unitary_with_order_four_not_two(self):
        labels, shift, u, t = finite_program(normalize([0.7, 0.2 + 0.3j]))
        np.testing.assert_allclose(u.conj().T @ u, np.eye(16), atol=4e-16)
        np.testing.assert_allclose(t.conj().T @ t, np.eye(16), atol=7e-16)
        np.testing.assert_allclose(np.linalg.matrix_power(t, 4), np.eye(16), atol=2e-15)
        self.assertAlmostEqual(np.linalg.norm(t @ t - np.eye(16), 2), 2)
        np.testing.assert_allclose(u.conj().T @ t @ u, np.kron(shift, np.eye(4)), atol=7e-16)
        np.testing.assert_allclose(labels[2], -labels[0])

    def test_10_finite_program_preserves_an_unknown_payload_with_reference_after_decoding(self):
        labels, shift, u, t = finite_program(normalize([0.7, 0.2 + 0.3j]))
        # A generic CQR pure input, not a known source-state assumption.
        payload = normalize(self.rng.normal(size=8) + 1j * self.rng.normal(size=8))
        initial = np.kron(np.eye(4)[0], payload)
        ur, tr = np.kron(u, I2), np.kron(t, I2)
        encoded = ur @ initial
        decoded = ur.conj().T @ tr @ tr @ encoded
        np.testing.assert_allclose(decoded, np.kron(np.eye(4)[2], payload), atol=1e-15)
        # On a new known |+>|0> source input the same unitary produces label -z.
        fresh = np.array([1., 0., 1., 0.], complex) / np.sqrt(2)
        output = t @ t @ u @ np.kron(np.eye(4)[0], fresh)
        source_vector = output.reshape(4, 4)[2]
        np.testing.assert_allclose(projector(source_vector), source(labels[2]), atol=6e-16)

    def test_11_uniform_implementation_error_gives_a_two_step_stability_bound(self):
        for angle in (0.03, 0.12, 0.4):
            gate = np.diag([np.exp(-0.5j * angle), np.exp(0.5j * angle)])
            eps = np.linalg.norm(gate - I2, 2)
            implemented = lambda z: gate @ j_reverse(z)
            for z in random_labels(self.rng, 25):
                self.assertLessEqual(np.linalg.norm(implemented(z) - j_reverse(z)), eps + 3e-16)
                error = np.linalg.norm(implemented(implemented(z)) - z)
                self.assertAlmostEqual(error, 2 * np.cos(angle / 2))
                self.assertGreaterEqual(error, 2 - 2 * eps - 1e-15)

    def test_12_a_cut_based_involution_pays_for_exact_return_with_discontinuity(self):
        for z in random_labels(self.rng, 30):
            w = discontinuous_involution(z)
            np.testing.assert_allclose(projector(w), I2 - projector(z), atol=3e-14)
            np.testing.assert_allclose(discontinuous_involution(w), z, atol=2e-14)
        north = np.array([1., 0.], complex)
        limit = discontinuous_involution(north)
        distances = []
        for eps in (0.03, 0.01, 0.003):
            z = normalize([1., eps])
            self.assertLess(np.linalg.norm(z - north), 2 * eps)
            distances.append(float(np.linalg.norm(discontinuous_involution(z) - limit)))
        self.assertGreater(min(distances), 1.99)

    def test_13_a_nonstandard_free_involution_saturates_the_dual_error_bound(self):
        b = np.array([0.2, -0.1, 0.15, 0.05])
        difference_from_antipode = []
        for z in random_labels(self.rng, 40):
            w = chord_involution(z, b)
            self.assertAlmostEqual(np.linalg.norm(w), 1)
            np.testing.assert_allclose(chord_involution(w, b), z, atol=8e-16)
            self.assertGreater(np.linalg.norm(w - z), 1)
            difference_from_antipode.append(np.linalg.norm(w + z))
            self.assertLessEqual(np.linalg.norm(projector(w) - (I2 - projector(z)), 2), 1 + 1e-15)
        self.assertGreater(max(difference_from_antipode), 0.4)
        z = complex_coordinates(b / np.linalg.norm(b))
        np.testing.assert_allclose(chord_involution(z, b), -z, atol=4e-16)
        self.assertAlmostEqual(np.linalg.norm(projector(chord_involution(z, b)) - (I2 - projector(z)), 2), 1)


if __name__ == "__main__":
    main(__name__, "free_reversal_dimension_audit", report)
