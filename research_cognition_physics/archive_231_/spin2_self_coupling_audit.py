"""Round 358: Deser's conditional first-order spin-two bootstrap.

Independent auxiliary-background variation, connection elimination, curvature,
and boundary checks accompany the analytic derivation. No cognition-to-spin-two
selection or unconstrained uniqueness among higher-derivative theories is made.
"""
from functools import lru_cache
import unittest
import numpy as np
from growing_stream_audit import main


ETA = np.diag([-1., 1., 1., 1.])
SLOTS = [(a, m, n) for a in range(4) for m in range(4) for n in range(m, 4)]
WAVE = np.array([.4, 1., -.35, .2])


def derivative(values):
    size = len(values)
    k = np.fft.fftfreq(size, d=1/size)
    return np.fft.ifft(np.fft.fft(values, axis=0) *
                       (1j*k).reshape((size,)+(1,)*(values.ndim-1)), axis=0).real


def pack(tensor):
    return np.stack([tensor[..., a, m, n] for a, m, n in SLOTS], axis=-1)


def unpack(vector):
    result = np.zeros(vector.shape[:-1]+(4, 4, 4))
    for slot, (a, m, n) in enumerate(SLOTS):
        result[..., a, m, n] = result[..., a, n, m] = vector[..., slot]
    return result


def trace_connection(connection):
    return np.einsum("xrar->xa", connection)


def compatibility(density, connection):
    """M_G(C)_a^mn in covariant differentiation of a weight-one density."""
    trace = trace_connection(connection)
    return (np.einsum("xmar,xrn->xamn", connection, density) +
            np.einsum("xnar,xmr->xamn", connection, density) -
            np.einsum("xa,xmn->xamn", trace, density))


@lru_cache(maxsize=None)
def compatibility_matrix():
    basis = unpack(np.eye(40))
    density = np.broadcast_to(ETA, (40, 4, 4))
    return pack(compatibility(density, basis)).T


def linear_inverse(rhs):
    return unpack(np.linalg.solve(compatibility_matrix(), pack(rhs).T).T)


def density_gradient(density):
    return np.einsum("a,xmn->xamn", WAVE, derivative(density))


def linear_connection(h):
    return -linear_inverse(density_gradient(h))


def nonlinear_connection_correction(h, connection):
    """B(h,C)=-M_eta^{-1} M_h(C), so C=C_L(h)+kappa B on shell."""
    return -linear_inverse(compatibility(h, connection))


def derivative_ricci(connection):
    dc = derivative(connection)
    trace_derivative = trace_connection(dc)
    return (np.einsum("a,xamn->xmn", WAVE, dc) -
            .5*np.einsum("m,xn->xmn", WAVE, trace_derivative) -
            .5*np.einsum("n,xm->xmn", WAVE, trace_derivative))


def quadratic_ricci(connection):
    trace = trace_connection(connection)
    return (np.einsum("xamn,xa->xmn", connection, trace) -
            np.einsum("xamb,xban->xmn", connection, connection))


def metric_from_density(density):
    determinant = np.linalg.det(density)
    if (np.any(determinant >= 0) or
            np.any(np.sum(np.linalg.eigvalsh(density) < 0, axis=-1) != 1)):
        raise ValueError("Use a nonsingular Lorentz-signature density branch.")
    return np.sqrt(-determinant)[:, None, None]*np.linalg.inv(density)


def density_from_metric(metric):
    return np.sqrt(-np.linalg.det(metric))[:, None, None]*np.linalg.inv(metric)


def christoffel(metric):
    inverse = np.linalg.inv(metric)
    dg = derivative(metric)
    result = np.zeros((len(metric), 4, 4, 4))
    for a in range(4):
        for m in range(4):
            for n in range(4):
                for r in range(4):
                    terms = (WAVE[m]*dg[:, r, n] + WAVE[n]*dg[:, r, m] -
                             WAVE[r]*dg[:, m, n])
                    result[:, a, m, n] += .5*inverse[:, a, r]*terms
    return result


def action_parts(h, connection, kappa):
    p, q = derivative_ricci(connection), quadratic_ricci(connection)
    density = ETA+kappa*h
    free = np.einsum("xmn,xmn->x", h, p)+np.einsum("mn,xmn->x", ETA, q)
    cubic = kappa*np.einsum("xmn,xmn->x", h, q)
    surface = np.einsum("mn,xmn->x", ETA, p)/kappa
    palatini = np.einsum("xmn,xmn->x", density, p/kappa+q)
    return free, cubic, surface, palatini


def integrated_palatini(density, connection):
    """Integration-by-parts density: no derivatives on the connection."""
    dg = derivative(density)
    trace = trace_connection(connection)
    kinetic = (-np.einsum("a,xmn,xamn->x", WAVE, dg, connection) +
               np.einsum("n,xmn,xm->x", WAVE, dg, trace))
    return float(np.mean(kinetic+
                         np.einsum("xmn,xmn->x", density,
                                   quadratic_ricci(connection))))


def auxiliary_action(h, connection, auxiliary_density):
    """Covariantized free action with C a tensor and h a weight-one density.

    This computes the Rosenfeld derivative independently of the bilinear B.
    """
    auxiliary_connection = christoffel(metric_from_density(auxiliary_density))
    trace_a = trace_connection(auxiliary_connection)
    trace_c = trace_connection(connection)
    correction = (
        np.einsum("xr,xrmn->xmn", trace_a, connection) -
        np.einsum("xram,xarn->xmn", auxiliary_connection, connection) -
        np.einsum("xran,xamr->xmn", auxiliary_connection, connection) +
        np.einsum("xrmn,xr->xmn", auxiliary_connection, trace_c))
    kinetic = np.einsum("xmn,xmn->x", h,
                        derivative_ricci(connection)+correction)
    potential = np.einsum("xmn,xmn->x", auxiliary_density,
                          quadratic_ricci(connection))
    return float(np.mean(kinetic+potential))


def auxiliary_directional_difference(h, connection, direction, step=1e-3):
    return (auxiliary_action(h, connection, ETA+step*direction) -
            auxiliary_action(h, connection, ETA-step*direction))/(2*step)


def fields(size=96):
    x = 2*np.pi*np.arange(size)/size
    rng = np.random.default_rng(358)
    matrices = []
    for _ in range(4):
        raw = rng.normal(size=(4, 4))
        matrices.append(.035*(raw+raw.T))
    h = (np.sin(x)[:, None, None]*matrices[0] +
         np.cos(2*x)[:, None, None]*matrices[1])
    direction = (np.cos(x)[:, None, None]*matrices[2] +
                 np.sin(3*x)[:, None, None]*matrices[3])
    c = unpack(.07*(np.sin(x)[:, None]*rng.normal(size=40) +
                     np.cos(2*x)[:, None]*rng.normal(size=40)))
    return h, c, direction


def linear_metric(h):
    trace = np.einsum("mn,xmn->x", ETA, h)
    return -np.einsum("ma,xab,bn->xmn", ETA, h, ETA)+.5*trace[:, None, None]*ETA


def linear_christoffel(metric_perturbation):
    """Direct differential formula; independent of the 40x40 metricity inverse."""
    dg = derivative(metric_perturbation)
    result = np.zeros((len(dg), 4, 4, 4))
    for a in range(4):
        for m in range(4):
            for n in range(4):
                for r in range(4):
                    result[:, a, m, n] += .5*ETA[a, r]*(
                        WAVE[m]*dg[:, r, n]+WAVE[n]*dg[:, r, m] -
                        WAVE[r]*dg[:, m, n])
    return result


def linear_einstein(metric_perturbation):
    ricci = derivative_ricci(linear_christoffel(metric_perturbation))
    trace = np.einsum("mn,xmn->x", ETA, ricci)
    return ricci-.5*trace[:, None, None]*ETA


def algebraic_weyl():
    """A Ricci-flat algebraic curvature tensor, not an assumed spacetime solution."""
    electric = np.diag([-2., 1., 1.])
    epsilon = np.zeros((3, 3, 3))
    for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        epsilon[a, b, c] = 1.
        epsilon[a, c, b] = -1.
    riemann = np.zeros((4, 4, 4, 4))
    for i in range(3):
        for j in range(3):
            riemann[0, i+1, 0, j+1] = electric[i, j]
            riemann[i+1, 0, 0, j+1] = -electric[i, j]
            riemann[0, i+1, j+1, 0] = -electric[i, j]
            riemann[i+1, 0, j+1, 0] = electric[i, j]
    riemann[1:, 1:, 1:, 1:] = -np.einsum(
        "ijm,kln,mn->ijkl", epsilon, epsilon, electric)
    return riemann


def curvature_cubic(riemann):
    raised = np.einsum("abmn,mc,nd->abcd", riemann, ETA, ETA)
    return float(np.einsum("abcd,cdef,efab", raised, raised, raised))


def maximum(tensor):
    return float(np.max(np.abs(tensor)))


def report():
    h, c, direction = fields()
    kappa = .7
    q = quadratic_ricci(c)
    b = nonlinear_connection_correction(h, c)
    candidate = q+derivative_ricci(b)
    expectation = float(np.mean(np.einsum("xmn,xmn->x", direction, candidate)))
    naive = float(np.mean(np.einsum("xmn,xmn->x", direction, q)))
    free, cubic, surface, palatini = action_parts(h, c, kappa)
    density = ETA+kappa*h
    gamma = christoffel(metric_from_density(density))
    eliminated_c = gamma/kappa
    eliminated_b = nonlinear_connection_correction(h, eliminated_c)
    eliminated_source = (quadratic_ricci(eliminated_c) +
                         derivative_ricci(eliminated_b))
    exact_ricci = derivative_ricci(gamma)+quadratic_ricci(gamma)
    source_identity = (derivative_ricci(linear_connection(h)) +
                       kappa*eliminated_source-exact_ricci/kappa)
    variations = []
    for eps in (.01, .005, .0025):
        observed = auxiliary_directional_difference(h, c, direction, eps)
        variations.append({"step": eps, "observed": observed,
                           "error": abs(observed-expectation)})
    step = 1e-4
    connection_derivative = (integrated_palatini(density, gamma+step*c) -
                             integrated_palatini(density, gamma-step*c))/(2*step)
    wrong_connection_derivative = (
        integrated_palatini(density, gamma+(.3+step)*c) -
        integrated_palatini(density, gamma+(.3-step)*c))/(2*step)
    riemann = algebraic_weyl()
    return {
        "round": 358,
        "scope": ("Conditional local two-derivative single massless spin-two "
                  "bootstrap. It supplies the nonlinear Einstein action after "
                  "the free Pauli-Fierz field and self-coupling requirements are "
                  "given; it does not derive this field, Lorentz geometry, or "
                  "the derivative restriction from cognitive principles. "
                  "Numerics sample smooth one-phase fields with all four "
                  "derivative components nonzero; analytic identities provide "
                  "the general quantifiers. No unrestricted uniqueness claim."),
        "sources": [
            "https://arxiv.org/abs/gr-qc/0411023",
            "https://arxiv.org/abs/0910.2975",
            "https://arxiv.org/abs/hep-th/0007220"],
        "phase_covector": WAVE.tolist(),
        "grid_points": len(h),
        "kappa": kappa,
        "compatibility_matrix_condition": float(np.linalg.cond(compatibility_matrix())),
        "pointwise_action_identity_residual": maximum(free+cubic+surface-palatini),
        "boundary_density_maximum": maximum(surface),
        "boundary_integral": float(np.mean(surface)),
        "auxiliary_stress": {
            "analytic_directional_variation": expectation,
            "naive_quadratic_Q_only": naive,
            "missing_derivative_term": expectation-naive,
            "finite_difference_sequence": variations,
            "richardson_error": abs((4*variations[1]["observed"] -
                                     variations[0]["observed"])/3-expectation)},
        "connection_elimination": {
            "metricity_residual": maximum(density_gradient(density) +
                                           compatibility(density, gamma)),
            "implicit_source_relation_residual": maximum(
                eliminated_c-linear_connection(h)-kappa*eliminated_b),
            "nonlinear_minus_linear_connection": maximum(
                eliminated_c-linear_connection(h)),
            "palatini_connection_variation_at_levi_civita": connection_derivative,
            "variation_at_wrong_connection": wrong_connection_derivative,
            "linear_ricci_plus_source_equals_full_ricci": maximum(source_identity),
            "chosen_metric_is_not_a_vacuum_solution": maximum(exact_ricci)},
        "higher_derivative_boundary": {
            "algebraic_ricci_residual": maximum(np.einsum("ac,abcd->bd", ETA, riemann)),
            "curvature_cubic": curvature_cubic(riemann),
            "cubic_at_amplitudes": [{"amplitude": eps,
                                     "value": curvature_cubic(eps*riemann)}
                                    for eps in (.1, .2, .4)],
            "interpretation": ("An independent six-derivative cubic curvature "
                               "operator has no quadratic flat-space term. "
                               "This is outside the two-derivative theorem, "
                               "not a claim of a healthy ultraviolet completion.")}
    }


class Checks(unittest.TestCase):
    def test_01_metric_density_roundtrip_and_signature(self):
        h, _, _ = fields()
        density = ETA+.7*h
        np.testing.assert_allclose(density_from_metric(metric_from_density(density)),
                                   density, atol=1e-15)
        with self.assertRaises(ValueError):
            metric_from_density(np.broadcast_to(np.eye(4), (4, 4, 4)))
        with self.assertRaises(ValueError):
            metric_from_density(np.broadcast_to(-ETA, (4, 4, 4)))

    def test_02_auxiliary_connection_inversion_is_full_rank(self):
        matrix = compatibility_matrix()
        self.assertEqual(np.linalg.matrix_rank(matrix), 40)
        self.assertLess(np.linalg.cond(matrix), 3.)
        _, connection, _ = fields()
        np.testing.assert_allclose(
            linear_inverse(compatibility(np.broadcast_to(ETA, (96, 4, 4)), connection)),
            connection, atol=3e-16)

    def test_03_free_connection_matches_direct_linear_christoffel(self):
        h, _, _ = fields()
        np.testing.assert_allclose(linear_connection(h),
                                   linear_christoffel(linear_metric(h)), atol=3e-15)

    def test_04_pure_gauge_has_zero_linear_curvature(self):
        x = 2*np.pi*np.arange(96)/96
        xi = np.sin(2*x)[:, None]*np.array([.2, -.1, .3, .5])
        gradient = np.einsum("m,xn->xmn", WAVE, derivative(xi))
        pure_gauge = gradient+gradient.transpose(0, 2, 1)
        self.assertLess(maximum(linear_einstein(pure_gauge)), 3e-13)
        h = linear_metric(pure_gauge)  # Trace reversal is an involution in 4D.
        self.assertLess(maximum(derivative_ricci(linear_connection(h))), 3e-13)

    def test_05_cubic_bootstrap_matches_full_palatini_pointwise(self):
        h, c, _ = fields()
        for kappa in (.2, .7, 1.1):
            free, cubic, surface, full = action_parts(h, c, kappa)
            np.testing.assert_allclose(free+cubic+surface, full, atol=8e-16)

    def test_06_surface_term_is_nonzero_but_integrates_to_zero(self):
        h, c, _ = fields()
        free, cubic, surface, full = action_parts(h, c, .7)
        self.assertGreater(maximum(surface), .1)
        self.assertLess(abs(float(np.mean(surface))), 5e-17)
        self.assertAlmostEqual(float(np.mean(free+cubic)), float(np.mean(full)), places=14)
        self.assertAlmostEqual(integrated_palatini(ETA+.7*h, .7*c)/.7**2,
                               float(np.mean(full)), places=14)

    def test_07_independent_auxiliary_variation_recovers_entire_source(self):
        data = report()["auxiliary_stress"]
        errors = [row["error"] for row in data["finite_difference_sequence"]]
        self.assertLess(data["richardson_error"], 2e-12)
        self.assertLess(errors[-1], 3e-10)
        self.assertGreater(errors[0]/errors[1], 3.5)
        self.assertGreater(errors[1]/errors[2], 3.5)

    def test_08_quadratic_Q_alone_is_not_the_full_stress_source(self):
        data = report()["auxiliary_stress"]
        self.assertGreater(abs(data["missing_derivative_term"]), 1e-4)

    def test_09_palatini_connection_variation_vanishes_only_on_correct_connection(self):
        data = report()["connection_elimination"]
        self.assertLess(abs(data["palatini_connection_variation_at_levi_civita"]), 1e-11)
        self.assertGreater(abs(data["variation_at_wrong_connection"]), 1e-4)

    def test_10_full_connection_has_exact_bilinear_implicit_relation(self):
        data = report()["connection_elimination"]
        self.assertLess(data["metricity_residual"], 2e-13)
        self.assertLess(data["implicit_source_relation_residual"], 2e-13)
        self.assertGreater(data["nonlinear_minus_linear_connection"], .001)

    def test_11_h_variation_gives_ricci_equation(self):
        h, c, direction = fields()
        kappa, eps = .7, 1e-5
        def action(perturbation):
            free, cubic, _, _ = action_parts(perturbation, c, kappa)
            return np.mean(free+cubic)
        actual = (action(h+eps*direction)-action(h-eps*direction))/(2*eps)
        expected = np.mean(np.einsum("xmn,xmn->x", direction,
                             derivative_ricci(c)+kappa*quadratic_ricci(c)))
        self.assertAlmostEqual(float(actual), float(expected), places=11)

    def test_12_elimination_yields_full_self_source_off_vacuum(self):
        data = report()["connection_elimination"]
        self.assertLess(data["linear_ricci_plus_source_equals_full_ricci"], 2e-12)
        self.assertGreater(data["chosen_metric_is_not_a_vacuum_solution"], .01)

    def test_13_explicit_local_improvement_is_linear_einstein_exact(self):
        x = 2*np.pi*np.arange(96)/96
        scalar = np.sin(x)+.3*np.cos(3*x)
        kk = float(WAVE @ ETA @ WAVE)
        improvement = derivative(derivative(scalar))[:, None, None]*(
            np.outer(WAVE, WAVE)-ETA*kk)
        local_shift = -scalar[:, None, None]*ETA
        np.testing.assert_allclose(linear_einstein(local_shift), improvement, atol=2e-13)
        divergence = np.einsum("m,xmn->xn", ETA @ WAVE, derivative(improvement))
        self.assertLess(maximum(divergence), 3e-13)

    def test_14_higher_derivative_curvature_invariant_is_not_a_ricci_polynomial(self):
        curvature = algebraic_weyl()
        np.testing.assert_allclose(curvature, -curvature.swapaxes(0, 1))
        np.testing.assert_allclose(curvature, -curvature.swapaxes(2, 3))
        np.testing.assert_allclose(curvature, curvature.transpose(2, 3, 0, 1))
        np.testing.assert_allclose(curvature+curvature.transpose(0, 2, 3, 1) +
                                   curvature.transpose(0, 3, 1, 2), 0.)
        np.testing.assert_allclose(np.einsum("ac,abcd->bd", ETA, curvature), 0.)
        self.assertAlmostEqual(curvature_cubic(curvature), 96.)

    def test_15_curvature_cube_leaves_the_quadratic_free_action_unchanged(self):
        curvature = algebraic_weyl()
        for amplitude in (.01, .1, .2, .4):
            self.assertAlmostEqual(curvature_cubic(amplitude*curvature),
                                   96*amplitude**3, places=12)
        step = 1e-3
        second = (curvature_cubic(step*curvature) +
                  curvature_cubic(-step*curvature))/step**2
        self.assertEqual(second, 0.)


if __name__ == "__main__":
    main(__name__, "spin2_self_coupling_audit", report)
