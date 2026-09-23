"""Round 366: classical parametrized scalar field and genuine slice changes.

Specified Minkowski cylinder, smooth spacelike embeddings, canonical continuum
Poisson brackets. FFT grids check selected smooth data; they do not define an
exact finite-grid constraint algebra or a quantum anomaly-free representation.
"""
import unittest
import numpy as np
from growing_stream_audit import main


def grid(size=128):
    return 2*np.pi*np.arange(size)/size


def derivative(values):
    size = values.shape[-1]
    modes = np.fft.fftfreq(size, 1/size)
    return np.fft.ifft(1j*modes*np.fft.fft(values, axis=-1), axis=-1).real


def integral(values):
    return float(2*np.pi*np.mean(values))


def geometry(state):
    # Store X=x+xi, never take an FFT derivative of the winding coordinate x.
    tp, xp, fp = derivative(state[:3])
    xp = 1+xp
    q = xp*xp-tp*tp
    if np.min(q) <= 0 or np.min(xp) <= 0:
        raise ValueError("An oriented spacelike embedding is required.")
    return tp, xp, fp, q, np.sqrt(q)


def constraints(state):
    tp, xp, fp, q, radius = geometry(state)
    pt, px, momentum = state[3:]
    numerator = pt*xp+px*tp+(momentum**2+fp**2)/2
    normal = numerator/radius
    tangent = pt*tp+px*xp+momentum*fp
    plus = (pt+px)*(tp+xp)/2+(momentum+fp)**2/4
    minus = (pt-px)*(tp-xp)/2-(momentum-fp)**2/4
    return normal, tangent, plus, minus


def functional(state, lapse, shift):
    normal, tangent, _, _ = constraints(state)
    return integral(lapse*normal+shift*tangent)


def gradient(state, lapse, shift):
    """Exact continuum functional gradients evaluated on a periodic grid."""
    tp, xp, fp, q, radius = geometry(state)
    pt, px, momentum = state[3:]
    numerator = pt*xp+px*tp+(momentum**2+fp**2)/2
    grad_q = -derivative(np.array([
        lapse*(px/radius+numerator*tp/radius**3)+shift*pt,
        lapse*(pt/radius-numerator*xp/radius**3)+shift*px,
        lapse*fp/radius+shift*momentum]))
    grad_p = np.array([lapse*xp/radius+shift*tp,
                       lapse*tp/radius+shift*xp,
                       lapse*momentum/radius+shift*fp])
    return grad_q, grad_p


def bracket(grad_f, grad_g):
    fq, fp = grad_f
    gq, gp = grad_g
    return integral(np.sum(fq*gp-fp*gq, axis=0))


def chiral_gradient(state, smear, sign):
    tp, xp, fp, _, _ = geometry(state)
    pt, px, momentum = state[3:]
    p_chiral = (pt+sign*px)/2
    xp_chiral = tp+sign*xp
    current = momentum+sign*fp
    gx = -derivative(smear*p_chiral)
    grad_q = np.array([gx, sign*gx, -derivative(smear*current/2)])
    grad_p = np.array([smear*xp_chiral/2,
                       sign*smear*xp_chiral/2,
                       sign*smear*current/2])
    return grad_q, grad_p


def off_shell(size=128):
    x = grid(size)
    state = np.array([
        .16*np.sin(x)+.035*np.cos(3*x),
        .07*np.sin(2*x)+.01*np.cos(x),
        .3*np.cos(x)+.12*np.sin(3*x),
        .21+.17*np.cos(2*x)+.05*np.sin(x),
        -.11+.13*np.sin(3*x)+.04*np.cos(x),
        .19+.23*np.cos(2*x)-.08*np.sin(x)])
    f = .8+.21*np.cos(x)+.11*np.sin(2*x)
    g = -.13+.17*np.sin(x)-.09*np.cos(3*x)
    return state, f, g


def closure_case(size=128):
    state, f, g = off_shell(size)
    zero = np.zeros_like(f)
    h, d, plus, minus = constraints(state)
    _, _, _, q, _ = geometry(state)
    lie = f*derivative(g)-g*derivative(f)
    hf = gradient(state, f, zero)
    hg = gradient(state, g, zero)
    df = gradient(state, zero, f)
    dg = gradient(state, zero, g)
    cases = {
        "DD": (bracket(df, dg), integral(lie*d)),
        "DH": (bracket(df, hg), integral(f*derivative(g)*h)),
        "HH": (bracket(hf, hg), integral(lie*d/q)),
        "plus": (bracket(chiral_gradient(state, f, 1),
                         chiral_gradient(state, g, 1)), integral(lie*plus)),
        "minus": (bracket(chiral_gradient(state, f, -1),
                          chiral_gradient(state, g, -1)), integral(lie*minus)),
        "cross": (bracket(chiral_gradient(state, f, 1),
                          chiral_gradient(state, g, -1)), 0.)}
    return {key: {"poisson": value[0], "expected": value[1],
                  "absolute_error": abs(value[0]-value[1])}
            for key, value in cases.items()}


def finite_variation_case(step):
    state, f, g = off_shell()
    zero = np.zeros_like(f)
    gq, gp = gradient(state, g, zero)
    direction = np.concatenate([gp, -gq])
    measured = (functional(state+step*direction, f, zero)
                -functional(state-step*direction, f, zero))/(2*step)
    expected = bracket(gradient(state, f, zero), (gq, gp))
    return {"step": step, "directional_difference": measured,
            "poisson": expected, "absolute_error": abs(measured-expected)}


def path_coefficients(time, path):
    if path not in ("A", "B"):
        raise ValueError("Choose path A or B.")
    bend = int(path == "B")
    a = .15*time+bend*.025*np.sin(np.pi*time)
    b = .08*time-bend*.02*np.sin(np.pi*time)
    adot = .15+bend*.025*np.pi*np.cos(np.pi*time)
    bdot = .08-bend*.02*np.pi*np.cos(np.pi*time)
    addot = -bend*.025*np.pi**2*np.sin(np.pi*time)
    bddot = bend*.02*np.pi**2*np.sin(np.pi*time)
    return a, b, adot, bdot, addot, bddot


def embedding(time, x, path):
    a, b, adot, bdot, _, _ = path_coefficients(time, path)
    t = .4*time+a*np.sin(x)
    xi = b*np.sin(2*x)
    td = .4+adot*np.sin(x)
    xd = bdot*np.sin(2*x)
    tp = a*np.cos(x)
    xp = 1+2*b*np.cos(2*x)
    return t, xi, td, xd, tp, xp


def smears(time, x, path):
    _, _, td, xd, tp, xp = embedding(time, x, path)
    q = xp*xp-tp*tp
    lapse = (xp*td-tp*xd)/np.sqrt(q)
    shift = (-tp*td+xp*xd)/q
    return lapse, shift


def ambient(u, v):
    f = .3*np.sin(2*u)+.1*np.cos(3*u)
    g = .2*np.cos(v)-.07*np.sin(4*v)
    fu = .6*np.cos(2*u)-.3*np.sin(3*u)
    gv = -.2*np.sin(v)-.28*np.cos(4*v)
    return f+g, fu, gv


def pulled_state(time, x, path):
    t, xi, _, _, tp, xp = embedding(time, x, path)
    field, fu, gv = ambient(t+x+xi, t-x-xi)
    plus_derivative, minus_derivative = tp+xp, tp-xp
    momentum = plus_derivative*fu-minus_derivative*gv
    pt = -plus_derivative*fu**2+minus_derivative*gv**2
    px = -plus_derivative*fu**2-minus_derivative*gv**2
    return np.array([t, xi, field, pt, px, momentum])


def flow(time, state, x, path):
    # Prescribed lapse/shift are held fixed while taking functional derivatives.
    lapse, shift = smears(time, x, path)
    gq, gp = gradient(state, lapse, shift)
    return np.concatenate([gp, -gq])


def evolve(path, steps, size=128):
    x = grid(size)
    state = pulled_state(0., x, path)
    dt = 1/steps
    for step in range(steps):
        time = step*dt
        k1 = flow(time, state, x, path)
        k2 = flow(time+dt/2, state+dt*k1/2, x, path)
        k3 = flow(time+dt/2, state+dt*k2/2, x, path)
        k4 = flow(time+dt, state+dt*k3, x, path)
        state += dt*(k1+2*k2+2*k3+k4)/6
    return state


def evolution_case(path, steps, size=128):
    x = grid(size)
    state = evolve(path, steps, size)
    target = pulled_state(1., x, path)
    h, d, _, _ = constraints(state)
    return {"path": path, "steps": steps, "grid_size": size,
            "six_field_max_error": float(np.max(np.abs(state-target))),
            "normal_constraint_max": float(np.max(np.abs(h))),
            "tangent_constraint_max": float(np.max(np.abs(d))),
            "killing_energy": integral(-state[3])}


def finite_difference(function, point, axis, step=1e-4):
    point = np.array(point, dtype=float)
    values = []
    for multiple in (-2, -1, 1, 2):
        offset = point.copy()
        offset[axis] += multiple*step
        values.append(function(*offset))
    return (values[0]-8*values[1]+8*values[2]-values[3])/(12*step)


def embedding_jets(time, x, path="B"):
    a, b, adot, bdot, addot, bddot = path_coefficients(time, path)
    jac = np.array([[.4+adot*np.sin(x), a*np.cos(x)],
                    [bdot*np.sin(2*x), 1+2*b*np.cos(2*x)]])
    second = np.array([
        [[addot*np.sin(x), adot*np.cos(x)],
         [adot*np.cos(x), -a*np.sin(x)]],
        [[bddot*np.sin(2*x), 2*bdot*np.cos(2*x)],
         [2*bdot*np.cos(2*x), -4*b*np.sin(2*x)]]])
    return jac, second


def metric(time, x):
    jac, _ = embedding_jets(time, x)
    return jac.T@np.diag([-1., 1.])@jac


def connection(time, x):
    jac, second = embedding_jets(time, x)
    return np.einsum("am,mbc->abc", np.linalg.inv(jac), second)


def curvature_case(time=.43, x=1.13):
    gamma = connection(time, x)
    dg = np.array([finite_difference(metric, (time, x), axis)
                   for axis in (0, 1)])
    invg = np.linalg.inv(metric(time, x))
    metric_gamma = np.zeros((2, 2, 2))
    for a, b, c, e in np.ndindex(2, 2, 2, 2):
        metric_gamma[a, b, c] += .5*invg[a, e]*(
            dg[b, e, c]+dg[c, e, b]-dg[e, b, c])
    dgamma = np.array([finite_difference(connection, (time, x), axis)
                       for axis in (0, 1)])
    riemann = np.zeros((2, 2, 2, 2))
    for a, b, c, d in np.ndindex(2, 2, 2, 2):
        riemann[a, b, c, d] = dgamma[c, a, d, b]-dgamma[d, a, c, b]
        for e in range(2):
            riemann[a, b, c, d] += (
                gamma[a, c, e]*gamma[e, d, b]
                -gamma[a, d, e]*gamma[e, c, b])
    return {"metric_connection_agreement": float(np.max(abs(metric_gamma-gamma))),
            "connection_nonzero": float(np.max(abs(gamma))),
            "riemann_max": float(np.max(abs(riemann))),
            "metric_determinant": float(np.linalg.det(metric(time, x)))}


def report():
    x = grid()
    initial_energy = integral(-pulled_state(0., x, "A")[3])
    exact_energy = 2*np.pi*(.6**2+.3**2+.2**2+.28**2)/2
    paths = [evolution_case(path, steps) for path in ("A", "B")
             for steps in (40, 80, 160)]
    endpoint_difference = float(np.max(abs(evolve("A", 160)-evolve("B", 160))))
    minimum_q, minimum_lapse = 1., 1.
    energy_errors, constraint_errors, stress = [], [], []
    for path in ("A", "B"):
        for time in np.linspace(0., 1., 21):
            state = pulled_state(time, x, path)
            minimum_q = min(minimum_q, float(np.min(geometry(state)[3])))
            minimum_lapse = min(minimum_lapse, float(np.min(smears(time, x, path)[0])))
            energy_errors.append(abs(integral(-state[3])-exact_energy))
            normal, tangent, _, _ = constraints(state)
            constraint_errors.append(float(np.max(abs(np.array([normal, tangent])))))
            t, xi = state[:2]
            _, fu, gv = ambient(t+x+xi, t-x-xi)
            stress.append(float(np.max(fu**2+gv**2)))
    return {
        "round": 366,
        "scope": ("Specified classical scalar on a fixed Minkowski cylinder; "
                  "continuum analytic HDA, sampled Fourier checks and full canonical "
                  "slice evolution. No independent metric dynamics, no quantum "
                  "constraint/anomaly result and no cognition-to-GR derivation."),
        "grid_size": len(x),
        "closure_off_shell": closure_case(),
        "functional_variation": [finite_variation_case(step)
                                 for step in (.01, .005, .0025)],
        "slice_evolution": paths,
        "two_path_endpoint_difference": endpoint_difference,
        "exact_slices": {"minimum_q": minimum_q, "minimum_lapse": minimum_lapse,
                         "maximum_constraint_error": max(constraint_errors),
                         "killing_energy": initial_energy,
                         "analytic_killing_energy": exact_energy,
                         "maximum_energy_error": max(energy_errors),
                         "maximum_ambient_energy_density": max(stress)},
        "flat_pullback_geometry": curvature_case(),
        "added_inputs": ["fixed smooth 1+1 Minkowski cylinder and massless scalar action",
                         "oriented invertible spacelike embeddings",
                         "canonical Poisson fields and periodic boundary conditions",
                         "prescribed two lapse/shift histories and smooth initial solution",
                         "FFT collocation and finite-step RK4 for numerical checks"],
        "not_selected": ["independent metric phase space or gravitational degrees of freedom",
                         "Einstein field equations or metric backreaction",
                         "dimension, matter action, quantum constraint representation"]}


class Checks(unittest.TestCase):
    def test_01_periodic_derivative_and_winding(self):
        x = grid()
        np.testing.assert_allclose(derivative(np.sin(3*x)), 3*np.cos(3*x), atol=8e-14)
        state, _, _ = off_shell()
        self.assertAlmostEqual(float(np.mean(geometry(state)[1])), 1., places=14)

    def test_02_chiral_and_normal_constraints_same_normalization(self):
        state, _, _ = off_shell()
        h, d, plus, minus = constraints(state)
        np.testing.assert_allclose(d, plus+minus, atol=1e-15)
        np.testing.assert_allclose(h, (plus-minus)/geometry(state)[4], atol=1e-15)

    def test_03_all_six_functional_derivatives(self):
        state, f, g = off_shell()
        analytic = np.concatenate(gradient(state, f, g))
        x = grid()
        direction = .11*np.cos(2*x)+.07*np.sin(3*x)
        eps = 2e-6
        for index in range(6):
            displacement = np.zeros_like(state)
            displacement[index] = direction
            difference = (functional(state+eps*displacement, f, g)
                          -functional(state-eps*displacement, f, g))/(2*eps)
            self.assertLess(abs(difference-integral(analytic[index]*direction)), 3e-9)

    def test_04_both_chiral_witt_algebras_off_shell(self):
        cases = closure_case()
        for sign in ("plus", "minus"):
            self.assertLess(cases[sign]["absolute_error"], 2e-13)
            self.assertGreater(abs(cases[sign]["expected"]), 1e-3)

    def test_05_opposite_chiral_sectors_commute(self):
        self.assertLess(closure_case()["cross"]["absolute_error"], 2e-13)

    def test_06_all_hypersurface_deformation_brackets_off_shell(self):
        for size in (64, 128, 256):
            cases = closure_case(size)
            for name in ("DD", "DH", "HH"):
                self.assertLess(cases[name]["absolute_error"], 2e-12)
                self.assertGreater(abs(cases[name]["expected"]), 1e-4)

    def test_07_actual_nonlinear_functional_variation_converges(self):
        cases = [finite_variation_case(step) for step in (.01, .005, .0025)]
        self.assertLess(cases[-1]["absolute_error"], 1e-6)
        for coarse, fine in zip(cases, cases[1:]):
            self.assertGreater(coarse["absolute_error"]/fine["absolute_error"], 3.95)
            self.assertLess(coarse["absolute_error"]/fine["absolute_error"], 4.05)

    def test_08_metric_structure_function_cannot_be_dropped(self):
        state, f, g = off_shell()
        _, d, _, _ = constraints(state)
        wrong = integral((f*derivative(g)-g*derivative(f))*d)
        self.assertGreater(abs(wrong-closure_case()["HH"]["poisson"]), .001)

    def test_09_exact_pulled_solution_satisfies_all_constraints(self):
        for path in ("A", "B"):
            for time in np.linspace(0., 1., 7):
                values = constraints(pulled_state(time, grid(), path))
                self.assertLess(float(np.max(abs(np.array(values)))), 2e-13)

    def test_10_full_canonical_flow_matches_slice_derivative(self):
        for path in ("A", "B"):
            x, time = grid(), .37
            expected = finite_difference(lambda t, dummy: pulled_state(t, x, path),
                                         (time, 0.), 0)
            actual = flow(time, pulled_state(time, x, path), x, path)
            self.assertLess(float(np.max(abs(expected-actual))), 1e-10)

    def test_11_two_paths_are_spacelike_future_foliations(self):
        for path in ("A", "B"):
            for time in np.linspace(0., 1., 21):
                state = pulled_state(time, grid(), path)
                self.assertGreater(float(np.min(geometry(state)[3])), .65)
                self.assertGreater(float(np.min(smears(time, grid(), path)[0])), .15)

    def test_12_two_independent_full_canonical_integrations(self):
        target = pulled_state(1., grid(), "A")
        final = []
        for path in ("A", "B"):
            coarse = evolve(path, 40)
            fine = evolve(path, 80)
            err0 = float(np.max(abs(coarse-target)))
            err1 = float(np.max(abs(fine-target)))
            self.assertGreater(err0/err1, 14.)
            self.assertLess(err0/err1, 18.)
            self.assertLess(err1, 2e-7)
            final.append(fine)
        self.assertLess(float(np.max(abs(final[0]-final[1]))), 3e-7)

    def test_13_nonzero_killing_energy_survives_vanishing_constraints(self):
        energy = 2*np.pi*(.6**2+.3**2+.2**2+.28**2)/2
        self.assertGreater(energy, 1.)
        for path in ("A", "B"):
            for time in (0., .27, .69, 1.):
                state = pulled_state(time, grid(), path)
                self.assertLess(abs(integral(-state[3])-energy), 3e-14)
            evolved = evolve(path, 80)
            self.assertLess(abs(integral(-evolved[3])-energy), 3e-14)

    def test_14_nontrivial_metric_components_still_flat(self):
        row = curvature_case()
        self.assertLess(row["metric_determinant"], -.05)
        self.assertGreater(row["connection_nonzero"], .1)
        self.assertLess(row["metric_connection_agreement"], 5e-11)
        self.assertLess(row["riemann_max"], 5e-11)

    def test_15_pulled_slice_changes_while_ambient_solution_fixed(self):
        x = grid()
        start = pulled_state(0., x, "A")
        end = pulled_state(1., x, "A")
        self.assertGreater(float(np.ptp(geometry(end)[3])), .5)
        self.assertGreater(float(np.max(abs(start[2]-end[2]))), .1)
        np.testing.assert_allclose(end, pulled_state(1., x, "B"), atol=1e-15)


if __name__ == "__main__":
    main(__name__, "parametrized_field_constraint_audit", report)
