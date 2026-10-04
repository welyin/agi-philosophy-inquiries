"""Round 351: conditional pure-metric hypersurface-deformation selection.

Canonical lambda multiplies pi**2 in H; it is NOT the lambda multiplying K**2
in the Lagrangian. The code tests strong, off-shell bracket identities, not
whether every alternative constraint algorithm is inconsistent. Analytic
proofs in the note cover smooth arbitrary fields. Numerics use a periodic
one-coordinate sector of d-dimensional metrics, without reducing tensor rank.
"""
import unittest
import numpy as np
from growing_stream_audit import main


def grid(size=96):
    return 2*np.pi*np.arange(size)/size


def derivative(values, order=1):
    """Spectral derivative along a periodic coordinate of period 2*pi."""
    n = len(values)
    frequencies = np.fft.fftfreq(n, d=1/n)
    multiplier = (1j*frequencies)**order
    shape = (n,)+(1,)*(values.ndim-1)
    return np.fft.ifft(np.fft.fft(values, axis=0) *
                       multiplier.reshape(shape), axis=0).real


def geometry(q):
    """Compute Ricci from Christoffel symbols, with only x^0 dependence.

    q has shape (grid, d, d). All tensor indices run over d coordinates.
    This path does not use the derived constraint bracket.
    """
    size, d, _ = q.shape
    inverse = np.linalg.inv(q)
    determinant = np.linalg.det(q)
    if np.any(determinant <= 0) or np.min(np.linalg.eigvalsh(q)) <= 0:
        raise ValueError("A positive definite spatial metric is required.")
    dq = derivative(q)
    connection = np.zeros((size, d, d, d))  # Gamma^k_ij
    for k in range(d):
        for i in range(d):
            for j in range(d):
                for ell in range(d):
                    term = ((dq[:, ell, j] if i == 0 else 0) +
                            (dq[:, ell, i] if j == 0 else 0) -
                            (dq[:, i, j] if ell == 0 else 0))
                    connection[:, k, i, j] += .5*inverse[:, k, ell]*term
    dc = derivative(connection)
    ricci = np.zeros_like(q)
    for i in range(d):
        for j in range(d):
            ricci[:, i, j] = dc[:, 0, i, j]
            if j == 0:
                for k in range(d):
                    ricci[:, i, j] -= dc[:, k, i, k]
            for k in range(d):
                for ell in range(d):
                    ricci[:, i, j] += (
                        connection[:, k, k, ell]*connection[:, ell, i, j] -
                        connection[:, k, j, ell]*connection[:, ell, i, k])
    scalar = np.einsum("xij,xij->x", inverse, ricci)
    return np.sqrt(determinant), scalar


def hamiltonian(q, momentum, lapse, lam, cosmological=0.,
                kinetic_scale=1., curvature_scale=1.):
    """Normalized torus integral; momentum is a contravariant density."""
    root, scalar = geometry(q)
    qp = q @ momentum
    trace = np.trace(qp, axis1=1, axis2=2)
    numerator = np.einsum("xij,xji->x", qp, qp)-lam*trace**2
    return float(np.mean(lapse*(kinetic_scale*numerator/root -
                               curvature_scale*root*scalar +
                               2*cosmological*root)))


def flat_gradients(momentum, lapse, lam, cosmological=0.,
                   kinetic_scale=1., curvature_scale=1.):
    """Full functional gradients evaluated at q=I, derived before smearing."""
    _, d, _ = momentum.shape
    identity = np.eye(d)
    trace = np.trace(momentum, axis1=1, axis2=2)
    numerator = np.einsum("xij,xij->x", momentum, momentum)-lam*trace**2
    dp = 2*kinetic_scale*lapse[:, None, None]*(
        momentum-lam*trace[:, None, None]*identity)
    dq = kinetic_scale*lapse[:, None, None]*(
        2*(momentum @ momentum)-2*lam*trace[:, None, None]*momentum -
        .5*numerator[:, None, None]*identity)
    lapse_second = derivative(lapse, 2)
    curvature_gradient = lapse_second[:, None, None]*np.broadcast_to(
        identity, momentum.shape).copy()
    curvature_gradient[:, 0, 0] -= lapse_second
    dq += curvature_scale*curvature_gradient
    dq += cosmological*lapse[:, None, None]*identity
    return dq, dp


def gradient_bracket(momentum, lapse_n, lapse_m, lam, cosmological=0.,
                     kinetic_scale=1., curvature_scale=1.):
    qn, pn = flat_gradients(momentum, lapse_n, lam, cosmological,
                            kinetic_scale, curvature_scale)
    qm, pm = flat_gradients(momentum, lapse_m, lam, cosmological,
                            kinetic_scale, curvature_scale)
    return float(np.mean(np.einsum("xij,xij->x", qn, pm) -
                         np.einsum("xij,xij->x", pn, qm)))


def predicted_bracket(momentum, lapse_n, lapse_m, lam,
                      kinetic_scale=1., curvature_scale=1.):
    """Integrated momentum term plus density-trace anomaly, at flat q."""
    _, d, _ = momentum.shape
    trace = np.trace(momentum, axis1=1, axis2=2)
    v = lapse_n*derivative(lapse_m)-lapse_m*derivative(lapse_n)
    momentum_term = float(-2*np.mean(v*derivative(momentum[:, 0, 0])))
    anomaly = float(2*(1-(d-1)*lam)*np.mean(v*derivative(trace)))
    product = kinetic_scale*curvature_scale
    return {"momentum_term": momentum_term,
            "anomaly": anomaly,
            "bracket": product*(momentum_term+anomaly),
            "scale_product": product}


def action_difference_bracket(momentum, lapse_n, lapse_m, lam, step=1e-3,
                              cosmological=0., kinetic_scale=1.,
                              curvature_scale=1.):
    """Differentiate the actual nonlinear curvature functional, not its Hessian.

    {H[N],H[M]} = d_q H[N](d_pi H[M])-d_q H[M](d_pi H[N]).
    The momentum argument stays fixed; off-diagonal symmetric entries are
    contracted in the full matrix convention, matching the canonical bracket.
    """
    size, d, _ = momentum.shape
    q = np.broadcast_to(np.eye(d), (size, d, d)).copy()
    _, pn = flat_gradients(momentum, lapse_n, lam, cosmological,
                           kinetic_scale, curvature_scale)
    _, pm = flat_gradients(momentum, lapse_m, lam, cosmological,
                           kinetic_scale, curvature_scale)

    def variation(lapse, direction):
        plus = hamiltonian(q+step*direction, momentum, lapse, lam,
                           cosmological, kinetic_scale, curvature_scale)
        minus = hamiltonian(q-step*direction, momentum, lapse, lam,
                            cosmological, kinetic_scale, curvature_scale)
        return (plus-minus)/(2*step)

    return variation(lapse_n, pm)-variation(lapse_m, pn)


def witness(d=3, size=96):
    x = grid(size)
    p = np.zeros((size, d, d))
    p[:, 1, 1] = np.sin(x)
    return p, np.ones(size), 2+np.sin(x)


def general_data(d=3, size=96):
    x = grid(size)
    p = np.zeros((size, d, d))
    for i in range(d):
        p[:, i, i] = .21*np.sin((i+1)*x)+.14*np.cos((i+2)*x)
        for j in range(i):
            p[:, i, j] = p[:, j, i] = .09*np.cos((i+j+1)*x)
    return p, 1+.2*np.cos(x)+.1*np.sin(2*x), (
        1.2+.35*np.sin(x)-.11*np.cos(3*x))


def canonical_lambda(lagrangian_lambda, d):
    denominator = d*lagrangian_lambda-1
    if abs(denominator) < 1e-14:
        raise ValueError("Singular Lagrangian trace Legendre map.")
    return lagrangian_lambda/denominator


def report():
    rows = []
    for d in (2, 3, 4):
        p, n, m = witness(d)
        for lam in (0., 1/d, 1/(d-1), 1.):
            predicted = predicted_bracket(p, n, m, lam)
            rows.append({"d": d, "canonical_lambda": lam,
                         "expected_trace_anomaly": 1-(d-1)*lam,
                         **predicted,
                         "full_gradient_bracket":
                         gradient_bracket(p, n, m, lam)})
    data = general_data()
    convergence = []
    target = predicted_bracket(*data, .17)["bracket"]
    for step in (2e-3, 1e-3, 5e-4):
        raw = action_difference_bracket(*data, .17, step=step)
        half = action_difference_bracket(*data, .17, step=step/2)
        richardson = (4*half-raw)/3
        convergence.append({"step": step, "action_bracket": raw,
                            "analytic_bracket": target,
                            "absolute_error": abs(raw-target),
                            "richardson_bracket": richardson,
                            "richardson_error": abs(richardson-target)})
    scales = []
    for a, b in ((1., 1.), (2., .5), (1., -1.), (1., 0.), (0., 1.)):
        row = predicted_bracket(*data, .5, kinetic_scale=a, curvature_scale=b)
        scales.append({"kinetic_scale": a, "curvature_scale": b, **row,
                       "full_gradient_bracket":
                       gradient_bracket(*data, .5, kinetic_scale=a,
                                        curvature_scale=b)})
    cosmological = [{"Lambda": value,
                     "bracket": gradient_bracket(*data, .5, value)}
                    for value in (-3., 0., .7, 12.)]
    return {
        "round": 351,
        "scope": {
            "kind": "conditional pure-metric classical constraint theorem",
            "selection": "Strong off-shell arbitrary-local-lapse Lorentzian "
                         "HDA selects canonical lambda=1/(d-1) in the "
                         "specified quadratic-momentum, scalar-R ansatz.",
            "added_inputs": ["smooth positive spatial metric in given d>=2",
                             "sole canonical variables q_ij and pi^ij",
                             "local quadratic ultralocal kinetic ansatz",
                             "scalar-curvature potential and spatial constant Lambda",
                             "given Lorentzian deformation algebra, arbitrary lapse",
                             "closed slice or vanishing boundary terms"],
            "not_proved": ["FUCP implies those inputs", "all alternatives inconsistent",
                           "selection of d=3 or Lambda or Newton constant",
                           "general full HKT uniqueness theorem",
                           "quantum constraint closure", "IR emergence of GR"],
            "exceptions": ["added CMC/trace constraints and lapse fixing",
                           "projectable lapse", "ultralocal/zero-product branch",
                           "opposite-sign Euclidean algebra",
                           "singular Legendre maps and d=1",
                           "extra canonical fields such as the f(R) scalar"]},
        "normalization": "mean over a 2*pi periodic coordinate; transverse "
                         "torus volume divided out in every functional",
        "witness_is_off_shell": True,
        "witness_rows": rows,
        "nonlinear_action_directional_difference": convergence,
        "coefficient_sign_and_ultralocal_rows": scales,
        "cosmological_constant_rows": cosmological,
        "lambda_convention": [{"d": d, "Lagrangian_lambda": 1.,
                              "canonical_lambda": canonical_lambda(1., d)}
                             for d in (2, 3, 4)],
        "sources": [
            "https://doi.org/10.1016/0003-4916(73)90096-1",
            "https://doi.org/10.1016/0003-4916(76)90112-3",
            "https://arxiv.org/html/1407.1259",
            "https://arxiv.org/html/1004.0055"]}


class Checks(unittest.TestCase):
    def test_01_spectral_derivatives_and_periodic_integration(self):
        x = grid()
        np.testing.assert_allclose(derivative(np.sin(3*x)), 3*np.cos(3*x),
                                   atol=2e-13)
        f, g = np.cos(2*x)+np.sin(x), np.sin(2*x)
        self.assertLess(abs(np.mean(f*derivative(g)+derivative(f)*g)), 1e-14)

    def test_02_ricci_from_connection_matches_conformal_identity(self):
        x = grid()
        omega = .13*np.sin(x)+.03*np.cos(2*x)
        for d in (2, 3, 4):
            q = np.exp(2*omega)[:, None, None]*np.eye(d)
            root, scalar = geometry(q)
            expected = np.exp(-2*omega)*(
                -2*(d-1)*derivative(omega, 2) -
                (d-1)*(d-2)*derivative(omega)**2)
            np.testing.assert_allclose(root, np.exp(d*omega), atol=2e-14)
            np.testing.assert_allclose(scalar, expected, atol=2e-11)

    def test_03_actual_action_variation_matches_full_gradient(self):
        p, n, _ = general_data()
        q = np.broadcast_to(np.eye(3), p.shape).copy()
        direction = .23*p+.07*np.eye(3)
        grad, _ = flat_gradients(p, n, .17, cosmological=.7)
        expected = float(np.mean(np.einsum("xij,xij->x", grad, direction)))
        step = 2e-4
        measured = (hamiltonian(q+step*direction, p, n, .17, .7) -
                    hamiltonian(q-step*direction, p, n, .17, .7))/(2*step)
        self.assertLess(abs(measured-expected), 3e-9)

    def test_04_full_gradients_match_integrated_bracket(self):
        for d in (2, 3, 4):
            data = general_data(d)
            for lam in (-.3, 0., 1/d, 1/(d-1), .91):
                self.assertAlmostEqual(
                    gradient_bracket(*data, lam),
                    predicted_bracket(*data, lam)["bracket"], places=12)

    def test_05_divergence_free_witness_has_nonzero_trace_gradient(self):
        for d in (2, 3, 4):
            p, n, m = witness(d)
            np.testing.assert_array_equal(p[:, 0, :], 0.)
            for lam in (0., 1/d, 1/(d-1), 1.):
                row = predicted_bracket(p, n, m, lam)
                self.assertEqual(row["momentum_term"], 0.)
                self.assertAlmostEqual(row["anomaly"], 1-(d-1)*lam, places=13)
                self.assertAlmostEqual(gradient_bracket(p, n, m, lam),
                                       row["anomaly"], places=13)

    def test_06_nonlinear_action_bracket_converges_independently(self):
        data = general_data()
        expected = predicted_bracket(*data, .17)["bracket"]
        coarse = action_difference_bracket(*data, .17, step=2e-3)
        fine = action_difference_bracket(*data, .17, step=1e-3)
        self.assertLess(abs(fine-expected), abs(coarse-expected)/3.5)
        self.assertLess(abs((4*fine-coarse)/3-expected), 3e-10)

    def test_07_selected_lambda_closes_with_nonzero_momentum_constraint(self):
        for d in (2, 3, 4):
            data = general_data(d)
            row = predicted_bracket(*data, 1/(d-1))
            self.assertGreater(abs(row["momentum_term"]), 1e-4)
            self.assertEqual(row["anomaly"], 0.)
            self.assertAlmostEqual(gradient_bracket(*data, 1/(d-1)),
                                   row["momentum_term"], places=13)

    def test_08_cosmological_constant_is_not_selected(self):
        data = general_data()
        target = gradient_bracket(*data, .5, 0.)
        for value in (-3., .7, 12.):
            self.assertAlmostEqual(gradient_bracket(*data, .5, value),
                                   target, places=12)

    def test_09_kinetic_curvature_product_controls_algebra_sign(self):
        data = general_data()
        reference = predicted_bracket(*data, .5)["bracket"]
        for a, b in ((1., 1.), (2., .5), (1., -1.), (-2., -.5)):
            self.assertAlmostEqual(
                gradient_bracket(*data, .5, kinetic_scale=a, curvature_scale=b),
                a*b*reference, places=13)

    def test_10_ultralocal_and_projectable_branches_do_not_select_lambda(self):
        p, n, m = general_data()
        for lam in (-.3, .17, .5, 1.):
            self.assertAlmostEqual(gradient_bracket(
                p, n, m, lam, curvature_scale=0.), 0., places=14)
            self.assertAlmostEqual(gradient_bracket(
                p, n, m, lam, kinetic_scale=0.), 0., places=14)
            self.assertAlmostEqual(gradient_bracket(
                p, np.ones(len(p)), 2*np.ones(len(p)), lam), 0., places=14)

    def test_11_cmc_subset_removes_anomaly_without_global_selection(self):
        p, n, m = general_data()
        trace = np.trace(p, axis1=1, axis2=2)
        p = p-(trace-.41)[:, None, None]*np.eye(3)/3
        for lam in (-.3, .17, .5, 1.):
            row = predicted_bracket(p, n, m, lam)
            self.assertLess(abs(row["anomaly"]), 1e-14)
            self.assertAlmostEqual(gradient_bracket(p, n, m, lam),
                                   row["momentum_term"], places=13)

    def test_12_legendre_convention_and_degeneracies(self):
        rng = np.random.default_rng(351)
        for d in (2, 3, 4):
            k = rng.normal(size=(d, d))
            k = (k+k.T)/2
            for lam_l in (-.2, .7, 1., 1.4):
                p = k-lam_l*np.trace(k)*np.eye(d)
                lam_h = canonical_lambda(lam_l, d)
                restored = p-lam_h*np.trace(p)*np.eye(d)
                np.testing.assert_allclose(restored, k, atol=2e-14)
            self.assertAlmostEqual(canonical_lambda(1., d), 1/(d-1))
            with self.assertRaises(ValueError):
                canonical_lambda(1/d, d)
            # A different degeneracy: canonical trace Hessian at lambda_H=1/d.
            np.testing.assert_allclose(np.eye(d)-np.trace(np.eye(d))*np.eye(d)/d,
                                       0., atol=1e-14)

    def test_13_independent_action_recovers_witness_and_mesh_stability(self):
        for size in (64, 96):
            data = witness(3, size)
            coarse = action_difference_bracket(*data, 0., step=1e-3)
            fine = action_difference_bracket(*data, 0., step=5e-4)
            self.assertLess(abs((4*fine-coarse)/3-1.), 3e-10)


if __name__ == "__main__":
    main(__name__, "hypersurface_constraint_selection_audit", report)
