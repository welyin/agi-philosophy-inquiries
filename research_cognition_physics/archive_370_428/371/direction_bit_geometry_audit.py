"""Round 371: direction-bit theorem interfaces and operational angle recovery."""
import unittest
import numpy as np
from growing_stream_audit import main


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
PAULI = np.array([X, Y, Z])


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def rho(r):
    return (I + np.einsum("i,ijk->jk", np.asarray(r), PAULI)) / 2


def effect(y, visibility=1., offset=.5):
    if not (0 < visibility <= 1 and
            visibility/2 <= offset <= 1-visibility/2):
        raise ValueError("Invalid visibility/offset")
    if not np.isclose(np.linalg.norm(y), 1):
        raise ValueError("Unit measurement direction required")
    return offset*I + (visibility/2)*np.einsum("i,ijk->jk", y, PAULI)


def probability(r, y, visibility=1., offset=.5):
    return float(np.trace(rho(r) @ effect(y, visibility, offset)).real)


def unitary(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values)) @ vectors.conjugate().T


def rotation(axis, angle):
    n = unit(axis)
    cross = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]],
                      [-n[1], n[0], 0]])
    r = np.cos(angle)*np.eye(3)+(1-np.cos(angle))*np.outer(n, n)
    return r+np.sin(angle)*cross


def trace_distance(a, b):
    d = a-b
    return float(np.abs(np.linalg.eigvalsh((d+d.conjugate().T)/2)).sum()/2)


def reduced_b(state):
    return np.trace(state.reshape(2, 2, 2, 2), axis1=0, axis2=2)


def reconstruct_gram(maxima, at_maxima, visibility, offset):
    # Rows are preparation devices; at_maxima[i,j] tests preparation j
    # with the measurement aligned to preparation i's maximizing direction.
    radii = 2*(np.asarray(maxima)-offset)/visibility
    if np.any(radii <= 0):
        raise ValueError("Each calibration preparation must be noncentral")
    gram = radii[:, None]*2*(np.asarray(at_maxima)-offset)/visibility
    return (gram+gram.T)/2


def recover_cosine(gram, at_y, at_z, visibility, offset):
    if np.linalg.eigvalsh(gram)[0] <= 1e-12:
        raise ValueError("Calibration states do not span the direction space")
    vy = 2*(np.asarray(at_y)-offset)/visibility
    vz = 2*(np.asarray(at_z)-offset)/visibility
    return float(vy @ np.linalg.solve(gram, vz))


def calibration_data():
    # Nonorthogonal mixed states; no axes are supplied to the reconstruction.
    r = np.array([[.71, .10, -.05], [.12, .64, .17], [-.14, .19, .58]])
    y, z = unit([.3, -.6, .7]), unit([-.4, .1, .8])
    a, c = .64, .43
    norms = np.linalg.norm(r, axis=1)
    maxima = np.array([probability(v, unit(v), a, c) for v in r])
    at_max = np.array([[probability(vj, unit(vi), a, c) for vj in r]
                       for vi in r])
    at_y = np.array([probability(v, y, a, c) for v in r])
    at_z = np.array([probability(v, z, a, c) for v in r])
    return r, y, z, a, c, norms, maxima, at_max, at_y, at_z


def probability_geometry():
    dirs = np.array([[1., 0, 0], [0, 1., 0], [0, 0, 1.],
                     unit([1, 2, 3]), unit([2, -.5, 1])])
    p = np.array([[probability(x, y) for y in dirs] for x in dirs])
    gram = 2*p-1
    reflection = np.diag([-1., 1., 1.])
    reflected = dirs @ reflection.T
    return {"probabilities": p.tolist(), "gram_eigenvalues": np.linalg.eigvalsh(gram).tolist(),
            "gram_rank": int(np.linalg.matrix_rank(gram, tol=1e-12)),
            "gram_error": float(np.max(abs(gram-dirs @ dirs.T))),
            "reflection_gram_error": float(np.max(abs(gram-reflected @ reflected.T))),
            "oriented_volume_before": float(np.linalg.det(dirs[:3])),
            "oriented_volume_after": float(np.linalg.det(reflected[:3]))}


def bootstrap_geometry():
    r, y, z, a, c, norms, maxima, at_max, py, pz = calibration_data()
    gram = reconstruct_gram(maxima, at_max, a, c)
    estimate = recover_cosine(gram, py, pz, a, c)
    # This bound treats the calibration Gram as exact. Its estimation error
    # is a separate resource and is not silently included in eta.
    eta = 1e-4
    dy = eta*np.array([1, -1, 1])
    dz = eta*np.array([-1, 1, 1])
    noisy = recover_cosine(gram, py+dy, pz+dz, a, c)
    vy, vz = 2*(py-c)/a, 2*(pz-c)/a
    e = 2*np.sqrt(3)*eta/a
    bound = np.linalg.norm(np.linalg.inv(gram), 2)*(e*(np.linalg.norm(vy)+np.linalg.norm(vz))+e*e)
    return {"visibility": a, "offset": c, "radii": norms.tolist(),
            "gram": gram.tolist(), "gram_condition_number": float(np.linalg.cond(gram)),
            "gram_error": float(np.max(abs(gram-r @ r.T))),
            "cosine_recovered": estimate, "cosine_direct": float(y @ z),
            "angle_radians_recovered": float(np.arccos(np.clip(estimate, -1, 1))),
            "cosine_error": abs(estimate-float(y @ z)),
            "measurement_probability_error": eta, "perturbed_cosine_error": abs(noisy-estimate),
            "conditional_perturbation_bound": float(bound)}


def calibration_perturbation():
    r, y, z, a, c, _, maxima, at_max, py, pz = calibration_data()
    gram = reconstruct_gram(maxima, at_max, a, c)
    eta_q, eta_p, eta_read = 1e-4, 1e-4, 1e-4
    q_hat = maxima+eta_q*np.array([1, -1, 1])
    p_hat = at_max+eta_p*np.array([[1, -1, 1], [-1, 1, 1], [1, 1, -1]])
    gram_hat = reconstruct_gram(q_hat, p_hat, a, c)
    py_hat = py+eta_read*np.array([1, -1, 1])
    pz_hat = pz+eta_read*np.array([-1, 1, 1])
    # |r_i| <= 1 bounds each unsymmetrized Gram entry's error.
    entry_bound = 2*(eta_q+eta_p)/a+4*eta_q*eta_p/a**2
    gamma = 3*entry_bound
    lam = float(np.linalg.eigvalsh(gram)[0])
    if gamma >= lam:
        raise ValueError("Calibration bound does not guarantee invertibility")
    vy, vz = 2*(py-c)/a, 2*(pz-c)/a
    e = 2*np.sqrt(3)*eta_read/a
    inverse_part = gamma*np.linalg.norm(vy)*np.linalg.norm(vz)/(lam*(lam-gamma))
    readout_part = (e*(np.linalg.norm(vy)+np.linalg.norm(vz))+e*e)/(lam-gamma)
    recovered = recover_cosine(gram_hat, py_hat, pz_hat, a, c)
    return {"maxima_error_bound": eta_q, "cross_calibration_probability_error_bound": eta_p,
            "final_readout_probability_error_bound": eta_read,
            "gram_operator_error": float(np.linalg.norm(gram_hat-gram, 2)),
            "gram_operator_error_bound": gamma, "true_smallest_gram_eigenvalue": lam,
            "perturbed_smallest_gram_eigenvalue": float(np.linalg.eigvalsh(gram_hat)[0]),
            "recovered_cosine": recovered, "cosine_error": abs(recovered-float(y @ z)),
            "full_calibration_and_readout_error_bound": float(inverse_part+readout_part),
            "visibility_and_offset_assumed_known": True}


def labelled_counterexample():
    x = unit([.4, -.3, .7])
    labels = [np.diag([1., 0.]), np.diag([0., 1.])]
    states = [np.kron(rho(x), label) for label in labels]
    dirs = [x, -x, unit([.1, .2, .9]), np.array([1., 0, 0])]
    gaps = [abs(np.trace((states[0]-states[1]) @ np.kron(effect(y), I)))
            for y in dirs]
    label_effect = np.kron(I, labels[0])
    swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                     [0, 1, 0, 0], [0, 0, 0, 1]], complex)
    after = [swap @ s @ swap.conjugate().T for s in states]
    return {"direction_effect_gap": float(max(gaps)),
            "same_unique_direction_maximum": 1.,
            "full_state_trace_distance": trace_distance(*states),
            "readable_label_probability_gap": float(abs(np.trace((states[0]-states[1]) @ label_effect))),
            "label_gap_after_allowed_swap_then_direction_readout":
                float(abs(np.trace((after[0]-after[1]) @ np.kron(effect([0, 0, 1]), I))))}


def axial_state(n):
    return (np.kron(rho(n), rho(n))+np.kron(rho(-n), rho(-n)))/2


def axial_counterexample():
    n = unit([.2, .5, .7])
    y = unit([.3, -.1, .8])
    state = axial_state(n)
    read = lambda v: float(np.trace(state @ np.kron(rho(v), rho(v))).real)
    return {"state_sign_ambiguity": trace_distance(state, axial_state(-n)),
            "probability_at_y": read(y), "predicted_even_probability": float((1+(n @ y)**2)/4),
            "antipodal_probability_gap": abs(read(y)-read(-y)),
            "L_max": 0., "directed_assumption_1_satisfied": False}


def report():
    return {"round": 371,
            "scope": "Conditional direction-bit theorem audit, reference-free relative angles, "
                     "and extra-label/axis/propagation limitations. No unconditional spatial "
                     "dimension, manifold, natural propagation, or GR derivation.",
            "source": {"url": "https://arxiv.org/html/1206.0630v4",
                       "direction_definition": 14, "dimension_theorem": 26,
                       "two_qubit_theorem": 27, "geometry_protocol": 38,
                       "manifold_example": 39,
                       "numbering_correction": "Item 9 is a lemma, not the dimension theorem."},
            "pure_probability_geometry": probability_geometry(),
            "mixed_noisy_frame_free_bootstrap": bootstrap_geometry(),
            "finite_gram_calibration": calibration_perturbation(),
            "readable_extra_label": labelled_counterexample(),
            "axis_instead_of_oriented_direction": axial_counterexample(),
            "resources_and_limits": [
                "Known nonzero visibility, offset and outcome calibration.",
                "Independently repeatable preparations; no unknown-state cloning.",
                "Three linearly independent calibration Bloch vectors in the tested 3D case.",
                "Angles are invariant under a joint O(3) choice of coordinates.",
                "A single direction leaves an SO(2) stabilizer; Gram data alone do not fix chirality.",
                "Physical device/rotation/direction identification is extra input.",
                "No distance scale, gluing between laboratories, connection, or spacetime metric."]}


class Checks(unittest.TestCase):
    def test_01_trace_rule_and_effect_range(self):
        rng = np.random.default_rng(371)
        for _ in range(30):
            r = unit(rng.normal(size=3))*rng.uniform(0, 1)
            y = unit(rng.normal(size=3))
            for a, c in [(1., .5), (.64, .43)]:
                self.assertAlmostEqual(probability(r, y, a, c), c+a*(r @ y)/2, places=14)
                self.assertGreaterEqual(np.linalg.eigvalsh(effect(y, a, c))[0], -1e-14)
                self.assertLessEqual(np.linalg.eigvalsh(effect(y, a, c))[-1], 1+1e-14)

    def test_02_rotation_covariance_and_double_cover(self):
        axis, angle = unit([.2, .5, -.4]), 1.13
        r, y = np.array([.2, -.1, .6]), unit([.3, .4, .8])
        h = np.einsum("i,ijk->jk", axis, PAULI)/2
        u, rotated = unitary(h, angle), rotation(axis, angle)
        np.testing.assert_allclose(u @ rho(r) @ u.conjugate().T, rho(rotated @ r), atol=1e-14)
        self.assertAlmostEqual(probability(r, rotated @ y), probability(rotated.T @ r, y), places=14)
        np.testing.assert_allclose(unitary(h, angle+2*np.pi), -u, atol=1e-14)

    def test_03_assumption_2_for_qubit_is_radius_and_direction(self):
        for radius in [.12, .57, 1.]:
            n = unit([.4, -.7, .2])
            r = radius*n
            maximum = probability(r, n)-probability(r, -n)
            np.testing.assert_allclose(maximum*n, r, atol=1e-14)
            for y in [-n, unit([.1, .7, .9])]:
                self.assertLess(probability(r, y)-probability(r, -y), maximum)

    def test_04_gram_geometry_and_reflection(self):
        p = probability_geometry()
        self.assertEqual(p["gram_rank"], 3)
        self.assertLess(p["gram_error"], 1e-14)
        self.assertLess(p["reflection_gram_error"], 1e-14)
        self.assertEqual(p["oriented_volume_before"], -p["oriented_volume_after"])

    def test_05_noisy_mixed_bootstrap(self):
        p = bootstrap_geometry()
        self.assertLess(p["gram_error"], 1e-14)
        self.assertLess(p["cosine_error"], 1e-14)
        self.assertLessEqual(p["perturbed_cosine_error"], p["conditional_perturbation_bound"])

    def test_06_no_shared_frame_needed_for_angles(self):
        r, y, z, a, c, _, maxima, at_max, py, pz = calibration_data()
        rot = rotation([.2, -.3, .9], 1.83)
        r2, y2, z2 = r @ rot.T, rot @ y, rot @ z
        py2 = [probability(v, y2, a, c) for v in r2]
        pz2 = [probability(v, z2, a, c) for v in r2]
        gram = reconstruct_gram(maxima, at_max, a, c)
        np.testing.assert_allclose(py2, py, atol=1e-14)
        np.testing.assert_allclose(pz2, pz, atol=1e-14)
        self.assertAlmostEqual(recover_cosine(gram, py2, pz2, a, c), y @ z, places=14)

    def test_07_single_direction_does_not_establish_frame(self):
        n, m = unit([.3, .4, .5]), unit([-.2, .7, .1])
        r = rotation(n, .87)
        np.testing.assert_allclose(r @ n, n, atol=1e-14)
        self.assertGreater(np.linalg.norm(r @ m-m), .1)

    def test_08_singular_calibration_is_not_silently_inverted(self):
        with self.assertRaises(ValueError):
            recover_cosine(np.ones((3, 3)), [.5]*3, [.5]*3, 1, .5)
        with self.assertRaises(ValueError):
            reconstruct_gram([.5]*3, np.full((3, 3), .5), 1, .5)

    def test_09_readable_label_violates_minimality(self):
        p = labelled_counterexample()
        self.assertLess(p["direction_effect_gap"], 1e-14)
        self.assertAlmostEqual(p["full_state_trace_distance"], 1., places=14)
        self.assertAlmostEqual(p["readable_label_probability_gap"], 1., places=14)

    def test_10_hidden_label_requires_restricting_future_operations(self):
        self.assertAlmostEqual(labelled_counterexample()[
            "label_gap_after_allowed_swap_then_direction_readout"], 1., places=14)

    def test_11_axis_information_fails_oriented_direction_assumption(self):
        p = axial_counterexample()
        self.assertLess(p["state_sign_ambiguity"], 1e-14)
        self.assertLess(p["antipodal_probability_gap"], 1e-14)
        self.assertAlmostEqual(p["probability_at_y"], p["predicted_even_probability"], places=14)

    def test_12_full_two_qubit_local_tomography(self):
        rng = np.random.default_rng(37112)
        m = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        s = m @ m.conjugate().T
        s /= np.trace(s)
        basis = [np.kron(a, b) for a in [I, X, Y, Z] for b in [I, X, Y, Z]]
        rebuilt = sum(np.trace(s @ b)*b/4 for b in basis)
        np.testing.assert_allclose(rebuilt, s, atol=1e-14)

    def test_13_allowed_continuous_interaction_is_nonproduct(self):
        u = unitary(np.kron(Z, Z), np.pi/4)
        s = np.kron(rho([1, 0, 0]), rho([1, 0, 0]))
        b = reduced_b(u @ s @ u.conjugate().T)
        self.assertAlmostEqual(float(np.trace(b @ b).real), .5, places=14)

    def test_14_finite_gram_calibration_and_inverse_perturbation(self):
        p = calibration_perturbation()
        self.assertLessEqual(p["gram_operator_error"], p["gram_operator_error_bound"])
        self.assertLess(p["gram_operator_error_bound"], p["true_smallest_gram_eigenvalue"])
        self.assertGreater(p["perturbed_smallest_gram_eigenvalue"], 0)
        self.assertLessEqual(p["cosine_error"], p["full_calibration_and_readout_error_bound"])


if __name__ == "__main__":
    main(__name__, "direction_bit_geometry_audit", report)
