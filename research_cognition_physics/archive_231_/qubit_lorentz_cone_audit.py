"""Round 372: qubit cone/spinor correspondence and physical filter accounting.

The Pauli identification is an established mathematical interface. Assigning
event displacements to this cone, and selecting a qubit as the geometric
object, remain additional inputs. No spacetime or Einstein dynamics derived.
"""
import unittest
import numpy as np
from growing_stream_audit import main


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
SIGMA = np.array([I, X, Y, Z])
ETA = np.diag([1., -1., -1., -1.])


def hermitian(vector):
    return np.einsum("a,aij->ij", vector, SIGMA)


def coordinates(matrix):
    values = np.array([np.trace(s @ matrix)/2 for s in SIGMA])
    if np.max(abs(values.imag)) > 1e-10:
        raise ValueError("Hermitian input required.")
    return values.real


def density(bloch):
    return hermitian(np.r_[1., bloch])/2


def direction(vector):
    vector = np.asarray(vector, float)
    return vector/np.linalg.norm(vector)


def boost(rapidity, axis=(0., 0., 1.)):
    n_sigma = hermitian(np.r_[0., direction(axis)])
    return np.cosh(rapidity/2)*I+np.sinh(rapidity/2)*n_sigma


def rotation(angle, axis=(0., 0., 1.)):
    n_sigma = hermitian(np.r_[0., direction(axis)])
    return np.cos(angle/2)*I-1j*np.sin(angle/2)*n_sigma


def lorentz(matrix):
    return np.array([[np.trace(sm @ matrix @ sn @ matrix.conj().T).real/2
                      for sn in SIGMA] for sm in SIGMA])


def expected_boost(rapidity, axis):
    n = direction(axis)
    c, s = np.cosh(rapidity), np.sinh(rapidity)
    result = np.eye(4)
    result[0, 0] = c
    result[0, 1:] = s*n
    result[1:, 0] = s*n
    result[1:, 1:] += (c-1)*np.outer(n, n)
    return result


def sqrt_psd(matrix):
    values, vectors = np.linalg.eigh((matrix+matrix.conj().T)/2)
    if values.min() < -2e-13:
        raise ValueError("Positive matrix required.")
    return (vectors*np.sqrt(np.maximum(values, 0))) @ vectors.conj().T


def filters(matrix):
    """Largest scalar multiple K=cA that is trace nonincreasing."""
    scale = 1/np.linalg.norm(matrix, 2)
    success = scale*matrix
    failure = sqrt_psd(I-success.conj().T @ success)
    return success, failure, scale


def dilation(success):
    """Halmos unitary; input flag |0> gives success and failure Kraus blocks."""
    right = sqrt_psd(I-success.conj().T @ success)
    left = sqrt_psd(I-success @ success.conj().T)
    return np.block([[success, -left], [right, success.conj().T]])


def branch(incoming, kraus, reference_dim=1):
    total = np.kron(kraus, np.eye(reference_dim))
    return total @ incoming @ total.conj().T


def reference_marginal(matrix, reference_dim):
    return np.trace(matrix.reshape(2, reference_dim, 2, reference_dim),
                    axis1=0, axis2=2)


def instrument(incoming, kraus, reference_dim=1):
    values = [branch(incoming, k, reference_dim) for k in kraus]
    size = len(incoming)
    recorded = np.zeros((2*size, 2*size), complex)
    recorded[:size, :size], recorded[size:, size:] = values
    return values, recorded


def normalized_filter(matrix, incoming):
    output = matrix @ incoming @ matrix.conj().T
    return output/np.trace(output)


def trace_distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a-b)).sum()/2)


def random_density(size, seed):
    rng = np.random.default_rng(seed)
    matrix = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
    value = matrix @ matrix.conj().T
    return value/np.trace(value)


def spinor_null(spinor):
    return coordinates(np.outer(spinor, spinor.conj()))


def spectral_pair(weight_plus, weight_minus, axis):
    n_sigma = hermitian(np.r_[0., direction(axis)])
    return weight_plus*(I+n_sigma)/2+weight_minus*(I-n_sigma)/2


def affine_witness(rapidity=np.log(2.)):
    matrix = boost(rapidity)
    up, down = density((0, 0, 1)), density((0, 0, -1))
    actual = normalized_filter(matrix, (up+down)/2)
    mixture = (normalized_filter(matrix, up)+normalized_filter(matrix, down))/2
    return trace_distance(actual, mixture)


def cone_case(seed=372):
    rng = np.random.default_rng(seed)
    det_errors, eigen_errors, covariance_errors = [], [], []
    for _ in range(24):
        vector = rng.normal(size=4)
        matrix = hermitian(vector)
        det_errors.append(abs(np.linalg.det(matrix)-vector @ ETA @ vector))
        eigen_errors.append(np.max(abs(np.linalg.eigvalsh(matrix)-
                            [vector[0]-np.linalg.norm(vector[1:]),
                             vector[0]+np.linalg.norm(vector[1:])])))
        a = boost(rng.uniform(-1, 1), rng.normal(size=3)) @ rotation(
            rng.uniform(-2, 2), rng.normal(size=3))
        covariance_errors.append(np.max(abs(
            coordinates(a @ matrix @ a.conj().T)-lorentz(a) @ vector)))
    return {"determinant_max_error": float(max(det_errors)),
            "eigenvalue_max_error": float(max(eigen_errors)),
            "spinor_vector_action_max_error": float(max(covariance_errors))}


def filter_case(rapidity=np.log(2.)):
    matrix = boost(rapidity)
    k, failure, scale = filters(matrix)
    inputs = {"north": density((0, 0, 1)),
              "south": density((0, 0, -1)),
              "maximally_mixed": I/2,
              "equator_plus": density((1, 0, 0)),
              "interior": density((.2, -.3, .4))}
    probabilities = {name: float(np.trace(branch(rho, k)).real)
                     for name, rho in inputs.items()}
    bell = np.eye(2).ravel()/np.sqrt(2)
    bell_rho = np.outer(bell, bell.conj())
    branches, record = instrument(bell_rho, (k, failure), 2)
    probability = float(np.trace(branches[0]).real)
    conditional_ref = reference_marginal(branches[0]/probability, 2)
    average_ref = reference_marginal(sum(branches), 2)
    average_on_plus = sum(
        branch(inputs["equator_plus"], kraus) for kraus in (k, failure))
    return {"rapidity": float(rapidity),
            "maximal_filter_scale": float(scale),
            "success_K_diagonal": np.diag(k).real.tolist(),
            "failure_K_diagonal": np.diag(failure).real.tolist(),
            "success_probabilities": probabilities,
            "uniform_success_minimum": float(np.exp(-2*rapidity)),
            "normalized_map_nonaffinity_trace_distance": affine_witness(rapidity),
            "bell_success_probability": probability,
            "bell_conditional_reference_eigenvalues": np.linalg.eigvalsh(conditional_ref).tolist(),
            "bell_conditional_reference_change": trace_distance(conditional_ref, I/2),
            "bell_average_reference_change": trace_distance(average_ref, I/2),
            "recorded_instrument_trace": float(np.trace(record).real),
            "recorded_instrument_min_eigenvalue": float(np.linalg.eigvalsh(record).min()),
            "averaged_plus_state_bloch": (2*coordinates(average_on_plus)[1:]).tolist(),
            "unitary_dilation_error": float(np.linalg.norm(
                dilation(k).conj().T @ dilation(k)-np.eye(4), 2))}


def passive_pairing_case():
    a = boost(.7, (1, -.2, .4)) @ rotation(.5, (-.3, .6, 1))
    incoming = density((.2, -.1, .4))
    effect = .6*density(direction((1, 2, -.5)))
    inverse = np.linalg.inv(a)
    transformed_state = a @ incoming @ a.conj().T
    transformed_effect = inverse.conj().T @ effect @ inverse
    transformed_unit = inverse.conj().T @ inverse
    return {"pairing_error": float(abs(np.trace(effect @ incoming)-
                         np.trace(transformed_effect @ transformed_state))),
            "transformed_unit_distance_from_identity": float(
                np.linalg.norm(transformed_unit-I, 2)),
            "new_unit_state_pairing": float(np.trace(
                transformed_unit @ transformed_state).real),
            "unchanged_trace_state_value": float(np.trace(transformed_state).real)}


def report():
    a = boost(.8, (1, 2, -1)) @ rotation(.6, (2, -.4, 1))
    la = lorentz(a)
    gamma = .4
    damp0 = np.diag([1., np.sqrt(1-gamma)])
    damp1 = np.array([[0., np.sqrt(gamma)], [0., 0.]])
    damped = sum(branch(I/2, k) for k in (damp0, damp1))
    return {
        "round": 372,
        "scope": ("Established Herm_2(C) future-cone/spinor correspondence with "
                  "a complete CP filter instrument and arbitrary-reference accounting; "
                  "not a derivation of event coordinates, physical 3+1 dimension, "
                  "deterministic quantum boosts, spacetime dynamics or Einstein equations."),
        "prior_difference": ("Round 183 checked generic spin-factor homogeneity; "
                             "this round audits the complex-qubit determinant, "
                             "SL(2,C) action and physical filtering versus passive coordinates."),
        "cone_identity": cone_case(),
        "Lorentz_example": {"metric_residual": float(np.max(abs(la.T @ ETA @ la-ETA))),
                            "determinant": float(np.linalg.det(la)),
                            "future_time_component": float(la[0, 0]),
                            "double_cover_residual": float(np.max(abs(lorentz(-a)-la)))},
        "boost_filter": filter_case(),
        "passive_coordinate_pairing": passive_pairing_case(),
        "deterministic_channel_counterexample": {
            "amplitude_damping_gamma": gamma,
            "input_determinant": .25,
            "output_determinant": float(np.linalg.det(damped).real),
            "analytic_output_determinant": (1-gamma**2)/4,
            "output_bloch": (2*coordinates(damped)[1:]).tolist()},
        "higher_matrix_cone_obstruction": {
            "Herm_3_nonproportional_rank_one_sum_rank": 2,
            "Herm_3_sum_is_boundary": True,
            "Lorentz_nonproportional_extreme_ray_sum_is_interior": True},
        "normalized_state_order_obstruction": {
            "claim": "For normalized rho,sigma, sigma-rho >= 0 implies sigma=rho.",
            "distinct_state_difference_eigenvalues": np.linalg.eigvalsh(
                density((.2, 0, 0))-density((0, 0, .6))).tolist()},
        "added_inputs": ["specified operational qubit and Pauli frame",
                         "chosen rapidity, filter coupling and one blank two-level flag",
                         "success/failure readout retained with state-dependent weights",
                         "event-displacement identification remains unprovided",
                         "two nonnegative spectral weights and a unit effect direction "
                         "do not by themselves define physical echo coordinates"],
        "sources": ["https://arxiv.org/html/quant-ph/0212135",
                    "https://arxiv.org/pdf/quant-ph/0011111"]}


class Checks(unittest.TestCase):
    def test_01_hermitian_eigenvalues_determinant_and_coordinates(self):
        row = cone_case()
        self.assertLess(max(row.values()), 2e-14)
        rng = np.random.default_rng(1372)
        for _ in range(15):
            value = rng.normal(size=4)
            np.testing.assert_allclose(coordinates(hermitian(value)), value, atol=1e-15)

    def test_02_positive_cone_including_boundary_and_past(self):
        n = direction((1, 2, 3))
        for time, radius, allowed in ((1., .7, True), (1., 1., True),
                                      (1., 1.1, False), (-1., .7, False)):
            eigen = np.linalg.eigvalsh(hermitian(np.r_[time, radius*n]))
            self.assertEqual(bool(eigen.min() >= -1e-14), allowed)

    def test_03_spinors_and_rank_one_effect_rays(self):
        rng = np.random.default_rng(2372)
        for _ in range(12):
            psi = rng.normal(size=2)+1j*rng.normal(size=2)
            value = spinor_null(psi)
            self.assertAlmostEqual(float(value @ ETA @ value), 0., places=12)
            self.assertGreater(value[0], 0)
            np.testing.assert_allclose(spinor_null(np.exp(.7j)*psi), value, atol=1e-14)
            n = value[1:]/value[0]
            effect = .7*density(n)
            np.testing.assert_allclose(np.linalg.eigvalsh(effect), [0, .7], atol=1e-14)

    def test_04_sl2_action_preserves_metric_and_future(self):
        for rapidity in (-1.2, 0., .7):
            a = boost(rapidity, (1, 2, 3)) @ rotation(.53, (-2, 1, .3))
            la = lorentz(a)
            np.testing.assert_allclose(la.T @ ETA @ la, ETA, atol=3e-14)
            self.assertAlmostEqual(np.linalg.det(la), 1., places=13)
            self.assertGreaterEqual(la[0, 0], 1.-1e-14)

    def test_05_double_cover_and_composition(self):
        a, b = boost(.9, (1, -.2, .4)), rotation(.4, (.3, -.7, 1))
        np.testing.assert_allclose(lorentz(a @ b), lorentz(a) @ lorentz(b), atol=2e-15)
        np.testing.assert_allclose(lorentz(-a), lorentz(a), atol=1e-15)

    def test_06_boost_and_rotation_match_independent_formulas(self):
        np.testing.assert_allclose(lorentz(boost(.8, (1, 2, 3))),
                                   expected_boost(.8, (1, 2, 3)), atol=2e-15)
        angle = .7
        expected = np.eye(4)
        expected[1:3, 1:3] = [[np.cos(angle), -np.sin(angle)],
                              [np.sin(angle), np.cos(angle)]]
        np.testing.assert_allclose(lorentz(rotation(angle)), expected, atol=1e-15)

    def test_07_boost_is_not_a_trace_preserving_channel(self):
        a = boost(np.log(2.))
        np.testing.assert_allclose(np.linalg.eigvalsh(a.conj().T @ a), [.5, 2.])
        self.assertAlmostEqual(np.trace(a @ density((0, 0, 1)) @ a.conj().T).real, 2.)
        self.assertAlmostEqual(np.trace(a @ density((0, 0, -1)) @ a.conj().T).real, .5)

    def test_08_scaled_filter_complete_and_probability_formula(self):
        n = direction((1, -.7, .2))
        for chi in (0., .4, 1.1):
            a = boost(chi, n)
            k, failure, _ = filters(a)
            np.testing.assert_allclose(k.conj().T @ k+failure.conj().T @ failure, I,
                                       atol=1e-14)
            for r in (np.zeros(3), n, -n, np.array([.2, -.3, .4])):
                p = np.trace(branch(density(r), k)).real
                self.assertAlmostEqual(p, np.exp(-chi)*(np.cosh(chi)+np.sinh(chi)*(n@r)))

    def test_09_full_unitary_dilation_matches_both_branches(self):
        for a in (boost(np.log(2.)), boost(.6, (1, 2, -1)) @ rotation(.3)):
            # Submaximal scale avoids square-root rounding at a saturated singular value.
            k = .9*filters(a)[0]
            u = dilation(k)
            np.testing.assert_allclose(u.conj().T @ u, np.eye(4), atol=2e-14)
            incoming = np.kron(np.array([[1.], [0.]]), I)
            expected = np.vstack([k, sqrt_psd(I-k.conj().T @ k)])
            np.testing.assert_allclose(u @ incoming, expected, atol=2e-14)

    def test_10_complete_instrument_preserves_unknown_reference_on_average(self):
        a = boost(.6, (1, 2, -1)) @ rotation(.3)
        k, failure, _ = filters(a)
        for rdim in (1, 2, 3, 5):
            rho = random_density(2*rdim, 3372+rdim)
            branches, recorded = instrument(rho, (k, failure), rdim)
            self.assertAlmostEqual(np.trace(recorded).real, 1.)
            self.assertGreater(np.linalg.eigvalsh(recorded).min(), -2e-14)
            np.testing.assert_allclose(reference_marginal(sum(branches), rdim),
                                       reference_marginal(rho, rdim), atol=2e-14)
            np.testing.assert_allclose(sum(branches), branch(rho, k, rdim)+
                                       branch(rho, failure, rdim), atol=1e-15)

    def test_11_selected_branch_changes_reference_with_probability_accounted(self):
        row = filter_case()
        self.assertAlmostEqual(row["bell_success_probability"], .625)
        np.testing.assert_allclose(row["bell_conditional_reference_eigenvalues"], [.2, .8])
        self.assertAlmostEqual(row["bell_conditional_reference_change"], .3)
        self.assertLess(row["bell_average_reference_change"], 1e-14)

    def test_12_normalized_boost_is_nonlinear(self):
        self.assertAlmostEqual(affine_witness(), .3)
        self.assertAlmostEqual(affine_witness(0.), 0.)

    def test_13_averaged_filter_is_dephasing_not_boost(self):
        k, failure, _ = filters(boost(np.log(2.)))
        rho = density((.2, -.3, .4))
        actual = branch(rho, k)+branch(rho, failure)
        np.testing.assert_allclose(actual, density((.1, -.15, .4)), atol=1e-15)
        self.assertGreater(trace_distance(actual, normalized_filter(boost(np.log(2.)), rho)), .1)

    def test_14_passive_pairing_changes_the_order_unit(self):
        row = passive_pairing_case()
        self.assertLess(row["pairing_error"], 1e-15)
        self.assertAlmostEqual(row["new_unit_state_pairing"], 1.)
        self.assertGreater(row["transformed_unit_distance_from_identity"], .5)
        self.assertGreater(abs(row["unchanged_trace_state_value"]-1), .1)

    def test_15_two_spectral_weights_plus_direction(self):
        n = direction((1, -.3, .8))
        for plus, minus in ((3., 1.), (2., 0.), (.7, .7)):
            matrix = spectral_pair(plus, minus, n)
            np.testing.assert_allclose(coordinates(matrix),
                np.r_[(plus+minus)/2, (plus-minus)*n/2], atol=1e-15)
            self.assertAlmostEqual(np.linalg.det(matrix).real, plus*minus)

    def test_16_higher_rank_cone_is_not_a_lorentz_cone(self):
        a, b = np.diag([1., 0., 0.]), np.diag([0., 1., 0.])
        self.assertEqual(np.linalg.matrix_rank(a+b), 2)
        self.assertEqual(np.linalg.eigvalsh(a+b).min(), 0.)
        null1, null2 = np.array([1., 1., 0., 0.]), np.array([1., 0., 1., 0.])
        self.assertGreater((null1+null2) @ ETA @ (null1+null2), 0.)

    def test_17_deterministic_cp_map_need_not_preserve_the_determinant(self):
        gamma = .4
        a = np.diag([1., np.sqrt(1-gamma)])
        b = np.array([[0., np.sqrt(gamma)], [0., 0.]])
        np.testing.assert_allclose(a.conj().T @ a+b.conj().T @ b, I)
        output = branch(I/2, a)+branch(I/2, b)
        self.assertAlmostEqual(np.linalg.det(output).real, (1-gamma**2)/4)
        self.assertGreater(abs(np.linalg.det(output).real-.25), .03)

    def test_18_distinct_normalized_states_have_no_positive_order_difference(self):
        delta = density((.2, 0, 0))-density((0, 0, .6))
        self.assertAlmostEqual(np.trace(delta).real, 0.)
        np.testing.assert_allclose(np.linalg.eigvalsh(delta),
                                   [-np.sqrt(.1), np.sqrt(.1)], atol=1e-15)


if __name__ == "__main__":
    main(__name__, "qubit_lorentz_cone_audit", report)
