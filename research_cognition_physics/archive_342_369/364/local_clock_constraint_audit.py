"""Round 364: ordinary interaction, local clock consistency, and constraint closure."""
import unittest
import numpy as np
from growing_stream_audit import main
from collective_poisson_limit_audit import collective_full, tensor_qubits

X = np.array([[0., 1.], [1., 0.]], complex)
Y = np.array([[0., -1j], [1j, 0.]], complex)
Z = np.diag([1., -1.]).astype(complex)


def comm(a, b):
    return a @ b-b @ a


def unitary(h, time):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*time*values)) @ vectors.conjugate().T


def quantum_densities(n=1, omega=(.7, -.45), coupling=.35):
    spin = [2*op/n for op in collective_full(n)]
    identity = np.eye(2**n)
    a = [np.kron(op, identity) for op in spin]
    b = [np.kron(identity, op) for op in spin]
    h1 = omega[0]*a[0]+coupling*a[2] @ b[2]/2
    h2 = omega[1]*b[0]+coupling*a[2] @ b[2]/2
    d = coupling*(omega[1]*a[2] @ b[1]-omega[0]*a[1] @ b[2])
    return h1, h2, d


def densities(m, omega=(1., 1.), coupling=1.):
    return np.array([omega[0]*m[0, 0]+coupling*m[0, 2]*m[1, 2]/2,
                     omega[1]*m[1, 0]+coupling*m[0, 2]*m[1, 2]/2])


def density_gradient(m, site, omega=(1., 1.), coupling=1.):
    result = np.zeros((2, 3))
    result[site, 0] = omega[site]
    result[0, 2] = coupling*m[1, 2]/2
    result[1, 2] = coupling*m[0, 2]/2
    return result


def d_value(m, omega=(1., 1.), coupling=1.):
    return coupling*(omega[1]*m[0, 2]*m[1, 1]-omega[0]*m[0, 1]*m[1, 2])


def d_gradient(m, omega=(1., 1.), coupling=1.):
    result = np.zeros((2, 3))
    result[0, 1] = -coupling*omega[0]*m[1, 2]
    result[0, 2] = coupling*omega[1]*m[1, 1]
    result[1, 1] = coupling*omega[1]*m[0, 2]
    result[1, 2] = -coupling*omega[0]*m[0, 1]
    return result


def poisson(m, grad_f, grad_g):
    return float(2*np.sum(m*np.cross(grad_f, grad_g)))


def constraint_witnesses():
    z = .5
    x = -.5*z*z
    y = np.sqrt(1-z*z-x*x)
    first = np.array([[x, y, z], [x, -y, z]])
    z2 = 2*(np.sqrt(2)-1)
    second = np.array([[1-np.sqrt(2), 0., np.sqrt(z2)]]*2)
    output = []
    for m in (first, second):
        output.append({"spins": m.tolist(), "h1_h2": densities(m).tolist(),
                       "D_equals_bracket_h1_h2": d_value(m),
                       "bracket_h1_D": poisson(m, density_gradient(m, 0), d_gradient(m))})
    return output


def process_distance(u, v):
    """Trace distance of outputs on a normalized maximally entangled input."""
    overlap = np.trace(u.conjugate().T @ v)/len(u)
    return float(np.sqrt(max(0., 1-abs(overlap)**2)))


def fixed_clock_loop(step, h1, h2):
    u1, u2 = unitary(h1, step), unitary(h2, step)
    return {"step": step, "unknown_reference_witness_distance": process_distance(u2 @ u1, u1 @ u2),
            "order_difference_norm": float(np.linalg.norm(u2 @ u1-u1 @ u2, 2))}


def flat_unitary(t1, t2, h1, h2):
    """A deliberately added clock connection, preserving the original diagonal."""
    return unitary(h1+h2, (t1+t2)/2) @ unitary(h1-h2, (t1-t2)/2)


def flat_generators(t1, t2, h1, h2):
    total, difference = h1+h2, h1-h2
    rotation = unitary(total, (t1+t2)/2)
    dressed = rotation @ difference @ rotation.conjugate().T
    return (total+dressed)/2, (total-dressed)/2


def flat_generator_derivative(t1, t2, h1, h2):
    total = h1+h2
    k1, k2 = flat_generators(t1, t2, h1, h2)
    rdot = -.5j*comm(total, k1-k2)  # either clock derivative of dressed K
    return rdot/2, -rdot/2


def path_propagation(end, order, steps, h1, h2):
    position = np.zeros(2)
    result = np.eye(len(h1), dtype=complex)
    for axis in order:
        dt = end[axis]/steps
        for _ in range(steps):
            midpoint = position.copy()
            midpoint[axis] += dt/2
            generator = flat_generators(*midpoint, h1, h2)[axis]
            result = unitary(generator, dt) @ result
            position[axis] += dt
    return result


def report():
    h1, h2, _ = quantum_densities()
    end = (.6, .4)
    target = flat_unitary(*end, h1, h2)
    paths = []
    for steps in (200, 400):
        first = path_propagation(end, (0, 1), steps, h1, h2)
        second = path_propagation(end, (1, 0), steps, h1, h2)
        paths.append({"steps_per_axis": steps,
                      "path_difference_norm": float(np.linalg.norm(first-second, 2)),
                      "first_vs_exact_endpoint_norm": float(np.linalg.norm(first-target, 2)),
                      "second_vs_exact_endpoint_norm": float(np.linalg.norm(second-target, 2))})
    k1, k2 = flat_generators(*end, h1, h2)
    dk1, dk2 = flat_generator_derivative(*end, h1, h2)
    return {
        "round": 364,
        "scope": "Specified interacting two-spin family: ordinary common-time consistency does not make local energy densities first-class constraints or yield the gravitational hypersurface-deformation algebra.",
        "hypothesis_tested": "Do a shared history, local clocks, and interaction automatically turn ordinary evolution into gravitational gauge constraints?",
        "classical_constraint_witnesses": constraint_witnesses(),
        "quantum_fixed_clock_curvature_norm": float(np.linalg.norm(1j*comm(h1, h2), 2)),
        "fixed_clock_order_witnesses": [fixed_clock_loop(step, h1, h2) for step in (.2, .1, .05, .025)],
        "maximally_entangled_small_step_coefficient": float(.35*np.sqrt(.7**2+.45**2)),
        "added_flat_connection": {
            "endpoint": list(end),
            "curvature_norm": float(np.linalg.norm(dk2-dk1+1j*comm(k1, k2), 2)),
            "sum_generator_error": float(np.linalg.norm(k1+k2-h1-h2, 2)),
            "local_generator_change": float(np.linalg.norm(k1-h1, 2)),
            "diagonal_original_dynamics_error": float(np.linalg.norm(
                flat_unitary(.7, .7, h1, h2)-unitary(h1+h2, .7), 2)),
            "independent_midpoint_path_integration": paths},
        "explicit_inputs": ["Two unit classical spins with the inherited Pauli Poisson bracket",
                            "Transverse terms and ZZ coupling; symmetric assignment of interaction energy",
                            "For the constraint test only: proposed local clock momenta or h_v=0",
                            "For the positive extension only: an explicitly chosen nonlocal clock-dependent connection"],
        "not_proved": ["All possible constraint completions are inconsistent",
                       "Local clock readings are unphysical", "Einstein dynamics or an emergent spacetime",
                       "Physical construction of ideal canonical clocks or their resource-free use"],
        "sources": ["https://arxiv.org/html/1309.1103", "https://doi.org/10.1016/0003-4916(73)90096-1"],
    }


class Checks(unittest.TestCase):
    def test_01_poisson_density_bracket_and_sign_from_gradients(self):
        rng = np.random.default_rng(364)
        for _ in range(30):
            m = rng.normal(size=(2, 3)); m /= np.linalg.norm(m, axis=1)[:, None]
            omega = tuple(rng.normal(size=2)); coupling = float(rng.normal())
            actual = poisson(m, density_gradient(m, 0, omega, coupling),
                             density_gradient(m, 1, omega, coupling))
            self.assertAlmostEqual(actual, d_value(m, omega, coupling), places=13)

    def test_02_primary_constraint_surface_is_not_first_class(self):
        row = constraint_witnesses()[0]
        np.testing.assert_allclose(row["h1_h2"], 0., atol=1e-15)
        self.assertGreater(abs(row["D_equals_bracket_h1_h2"]), .8)
        np.testing.assert_allclose(np.linalg.norm(row["spins"], axis=1), 1., atol=1e-15)

    def test_03_adding_first_bracket_still_does_not_close_constraints(self):
        row = constraint_witnesses()[1]
        np.testing.assert_allclose(row["h1_h2"], 0., atol=2e-16)
        self.assertEqual(row["D_equals_bracket_h1_h2"], 0.)
        self.assertAlmostEqual(row["bracket_h1_D"], -4*(np.sqrt(2)-1), places=13)
        np.testing.assert_allclose(np.linalg.norm(row["spins"], axis=1), 1., atol=1e-15)

    def test_04_D_gradient_against_directional_difference(self):
        m = np.array([[.2, -.3, .4], [-.6, .1, .8]])
        direction = np.array([[.7, -.1, .2], [.4, .2, -.3]])
        eps = 1e-5
        derivative = (d_value(m+eps*direction)-d_value(m-eps*direction))/(2*eps)
        self.assertAlmostEqual(derivative, np.sum(d_gradient(m)*direction), places=11)

    def test_05_full_collective_commutator_at_three_N_values(self):
        for n in (1, 2, 3):
            h1, h2, d = quantum_densities(n)
            np.testing.assert_allclose(-1j*n*comm(h1, h2), d, atol=2e-15)

    def test_06_classical_curvature_is_product_quantum_expectation(self):
        m = np.array([[.3, .4, np.sqrt(.75)], [-.4, .2, np.sqrt(.8)]])
        omega, coupling = (.7, -.45), .35
        for n in (1, 2):
            states = [tensor_qubits(n, row[2], np.arctan2(row[1], row[0])) for row in m]
            psi = np.kron(*states)
            _, _, d = quantum_densities(n, omega, coupling)
            self.assertAlmostEqual(np.vdot(psi, d @ psi).real, d_value(m, omega, coupling), places=13)

    def test_07_fixed_clock_order_is_observable_with_a_reference(self):
        h1, h2, _ = quantum_densities()
        coefficients = []
        for step in (.1, .05, .025):
            row = fixed_clock_loop(step, h1, h2)
            coefficients.append(row["unknown_reference_witness_distance"]/step**2)
        expected = .35*np.sqrt(.7**2+.45**2)
        self.assertLess(abs(coefficients[-1]-expected), 3e-5)
        self.assertLess(abs(coefficients[-1]-expected), abs(coefficients[0]-expected))

    def test_08_noninteracting_and_no_transverse_controls(self):
        for omega, coupling in (((.7, -.45), 0.), ((0., 0.), .35)):
            h1, h2, _ = quantum_densities(1, omega, coupling)
            np.testing.assert_allclose(comm(h1, h2), 0., atol=1e-15)
            row = fixed_clock_loop(.7, h1, h2)
            self.assertLess(row["order_difference_norm"], 2e-15)

    def test_09_flat_connection_curvature_with_independent_differences(self):
        h1, h2, _ = quantum_densities()
        t1, t2, eps = .6, .4, 1e-4
        k1, k2 = flat_generators(t1, t2, h1, h2)
        def derivative(axis, component):
            values = []
            for multiple in (-2, -1, 1, 2):
                point = np.array([t1, t2]); point[axis] += multiple*eps
                values.append(flat_generators(*point, h1, h2)[component])
            return (values[0]-8*values[1]+8*values[2]-values[3])/(12*eps)
        f = derivative(0, 1)-derivative(1, 0)+1j*comm(k1, k2)
        self.assertLess(np.linalg.norm(f, 2), 8e-12)

    def test_10_flat_extension_preserves_diagonal_and_origin_generators(self):
        h1, h2, _ = quantum_densities()
        k1, k2 = flat_generators(0., 0., h1, h2)
        np.testing.assert_allclose(k1, h1, atol=2e-15)
        np.testing.assert_allclose(k2, h2, atol=2e-15)
        for time in (-.5, 0., .3, 1.2):
            np.testing.assert_allclose(flat_unitary(time, time, h1, h2),
                                       unitary(h1+h2, time), atol=3e-15)

    def test_11_two_independent_path_integrations_converge(self):
        h1, h2, _ = quantum_densities()
        end = (.6, .4); target = flat_unitary(*end, h1, h2)
        for order in ((0, 1), (1, 0)):
            coarse = np.linalg.norm(path_propagation(end, order, 200, h1, h2)-target, 2)
            fine = np.linalg.norm(path_propagation(end, order, 400, h1, h2)-target, 2)
            self.assertLess(fine, 1e-7)
            self.assertGreater(coarse/fine, 3.9)
            self.assertLess(coarse/fine, 4.1)

    def test_12_flat_extension_does_not_require_special_interaction(self):
        rng = np.random.default_rng(1364)
        for size in (2, 3, 5):
            mats = []
            for _ in range(2):
                a = rng.normal(size=(size, size))+1j*rng.normal(size=(size, size))
                mats.append((a+a.conjugate().T)/2)
            a, b = mats
            k1, k2 = flat_generators(.2, -.1, a, b)
            dk1, dk2 = flat_generator_derivative(.2, -.1, a, b)
            np.testing.assert_allclose(dk2-dk1+1j*comm(k1, k2), 0., atol=2e-14)
            np.testing.assert_allclose(k1+k2, a+b, atol=2e-15)


if __name__ == "__main__":
    main(__name__, "local_clock_constraint_audit", report)
