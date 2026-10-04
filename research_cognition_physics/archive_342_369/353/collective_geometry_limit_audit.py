"""Round 353: collective quantum mean limit, its time window and feedback cost.

The G register is an interface placeholder; this does not derive a metric or GR.
Finite full-joint matrices check closed-form central-spin formulae independently.
"""
import unittest
import numpy as np
from growing_stream_audit import main


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
PLUS = np.ones(2, complex)/np.sqrt(2)


def density(vector):
    return np.outer(vector, vector.conj())


def tensor_power(matrix, n):
    result = np.array(1., complex)
    for _ in range(n):
        result = np.kron(result, matrix)
    return result


def bloch(x=0., y=0., z=0.):
    return (I+x*X+y*Y+z*Z)/2


def product_geometry(n, mean, pure=True):
    x = np.sqrt(max(0., 1-mean*mean))*(1 if pure else .37)
    return tensor_power(bloch(x=x, z=mean), n)


def q_eigenvalues(n):
    return np.array([(n-2*index.bit_count())/n for index in range(2**n)])


def joint_phases(n, theta, reference_dimension=1):
    generator = np.kron(q_eigenvalues(n), [1., -1.])
    return np.repeat(np.exp(-1j*theta*generator), reference_dimension)


def joint_evolve(rho, n, theta, reference_dimension=1):
    phases = joint_phases(n, theta, reference_dimension)
    return phases[:, None]*rho*phases.conj()[None, :]


def matter_reference(rho, n, reference_dimension=1):
    d, m = 2**n, 2*reference_dimension
    return np.trace(rho.reshape(d, m, d, m), axis1=0, axis2=2)


def distance(first, second):
    return float(np.sum(abs(np.linalg.eigvalsh(first-second)))/2)


def coherence(n, mean, theta):
    return complex((np.cos(2*theta/n)-1j*mean*np.sin(2*theta/n))**n)


def centered_coherence(n, mean, theta):
    # Equivalent integer power, with the extensive mean phase removed first.
    z = np.cos(2*theta/n)-1j*mean*np.sin(2*theta/n)
    return complex(np.exp(n*(np.log(z)+2j*theta*mean/n)))


def channel(rho_mr, factor):
    r = rho_mr.shape[0]//2
    result = rho_mr.copy().reshape(2, r, 2, r)
    result[0, :, 1, :] *= factor
    result[1, :, 0, :] *= factor.conjugate()
    return result.reshape(2*r, 2*r)


def channel_error(n, mean, theta):
    return float(abs(centered_coherence(n, mean, theta)-1)/2)


def geometry_operator(n, slot, operator):
    result = np.array(1., complex)
    for j in range(n):
        result = np.kron(result, operator if j == slot else I)
    return result


def expectation(rho, operator):
    return float(np.trace(rho@operator).real)


def negativity(rho, n):
    d = 2**n
    transpose = rho.reshape(d, 2, d, 2).transpose(0, 3, 2, 1)
    return float((np.sum(abs(np.linalg.eigvalsh(transpose.reshape(2*d, 2*d))))-1)/2)


def cat_geometry(n, mean, coherent=True):
    vector = np.zeros(2**n, complex)
    vector[0], vector[-1] = np.sqrt((1+mean)/2), np.sqrt((1-mean)/2)
    rho = density(vector)
    return rho if coherent else np.diag(np.diag(rho))


def covariance_identity(rho, n):
    q = q_eigenvalues(n)
    probabilities = np.diag(rho).real
    z = np.array([[1-2*((index >> (n-1-j)) & 1)
                   for index in range(2**n)] for j in range(n)], float)
    means = z@probabilities
    covariance = (z*probabilities)@z.T-np.outer(means, means)
    variance = float(probabilities@(q*q)-(probabilities@q)**2)
    return variance, float(covariance.sum()/n**2), covariance


def fixed_time_rows():
    mean, theta = .3, .7
    return [{"N": n, "theta": theta, "mean": mean,
             "half_diamond": channel_error(n, mean, theta),
             "variance_bound": theta**2*(1-mean**2)/n,
             "N_times_error": n*channel_error(n, mean, theta)}
            for n in (10, 100, 1000, 10000)]


def scaled_time_rows():
    mean, c = .3, .6
    limit = np.exp(-2*c*c*(1-mean*mean))
    return [{"N": n, "theta": float(c*np.sqrt(n)),
             "centered_coherence_real": centered_coherence(n, mean, c*np.sqrt(n)).real,
             "centered_coherence_imag": centered_coherence(n, mean, c*np.sqrt(n)).imag,
             "half_diamond": channel_error(n, mean, c*np.sqrt(n)),
             "limiting_half_diamond": float((1-limit)/2)}
            for n in (100, 1000, 10000, 100000)]


def response_rows():
    # Pure geometry |+>^N and matter <Z>=0.6. Average Y is a norm-one quantity.
    result = []
    for n in (100, 1000, 10000):
        for label, theta in (("fixed", .7), ("sqrt_N", .6*np.sqrt(n)),
                             ("linear_N", .3*n)):
            result.append({"N": n, "scaling": label, "theta": float(theta),
                           "average_Y_response": float(.6*np.sin(2*theta/n)),
                           "sum_Y_response": float(.6*n*np.sin(2*theta/n)),
                           "matter_half_diamond_mean0": channel_error(n, 0., theta)})
    return result


def report():
    n, theta = 4, np.pi/4
    matter = density(PLUS)
    cat = joint_evolve(np.kron(cat_geometry(n, 0.), matter), n, theta)
    mixture = joint_evolve(np.kron(cat_geometry(n, 0., False), matter), n, theta)
    return {
        "round": 353,
        "scope": {
            "proved": [
                "For independent product geometry and arbitrary MR reference, exact half-diamond error is |Gamma_N-exp(-2i theta mu)|/2 <= min(1,theta^2(1-mu^2)/N).",
                "Fixed-theta convergence is not uniform: theta=c sqrt(N) leaves finite Gaussian dephasing in the mean rotating frame.",
                "Same-mean correlated cat geometry keeps Var(Q_N)=1-mu^2 and can prevent the mean-unitary limit.",
                "The same closed joint unitary changes noncommuting geometry observables and creates geometry-matter entanglement.",
                "One probe with fixed total coupling gives average geometry response O(theta/N); finite fixed-time macroscopic feedback is not obtained."],
            "inputs": ["N geometry qubits and one matter qubit", "Q_N=N^-1 sum_j Z_j and H=g Q_N tensor Z_M", "Initially fixed geometry independent of arbitrary MR", "Product geometry or separately specified correlated geometry", "No geometry or matter self-Hamiltonian"],
            "not_proved": ["Metric or Einstein dynamics", "Spatial dimension or locality", "A finite classical backreaction limit for macroscopic matter", "Universal validity of nonlinear semiclassical equations", "Preparation, reset or repeated-use resource closure"]},
        "fixed_time": fixed_time_rows(),
        "sqrt_N_time": scaled_time_rows(),
        "backreaction_scaling": response_rows(),
        "cat_versus_mixture": {
            "N": n, "theta": float(theta), "mean": 0., "variance": 1.,
            "cat_matter_error": distance(matter_reference(cat, n), matter),
            "mixture_matter_error": distance(matter_reference(mixture, n), matter),
            "same_reduced_channel_distance": distance(matter_reference(cat, n), matter_reference(mixture, n)),
            "cat_joint_negativity": negativity(cat, n),
            "mixture_joint_negativity": negativity(mixture, n)},
        "joint_account_example": {
            "N": 4, "theta": .7,
            "global_trace_distance_from_mean_product_mu0": float(np.sqrt(1-np.cos(.7/4)**8)),
            "matter_trace_distance_mu0": channel_error(4, 0., .7),
            "joint_negativity_mu0": float(np.sqrt(1-abs(coherence(4, 0., .7))**2)/2),
            "sum_Y_Z_covariance_mu0": float(4*np.sin(2*.7/4)),
            "average_Y_Z_covariance_mu0": float(np.sin(2*.7/4))},
        "resource_account": {"geometry_qubits": "N", "GM_Hilbert_dimension": "2^(N+1)", "interaction_operator_norm": "abs(g)", "individual_pair_coupling": "g/N", "time": "theta/g", "recurrence_at_theta": "pi*N gives (-1)^N I on GM"},
        "primary_source": "Cucchietti, Paz, Zurek, Phys. Rev. A 72, 052113 (2005), https://arxiv.org/html/quant-ph/0508184, section II equations (1), (8), (12)-(16), (22). The product coherence and Gaussian mechanism are known; this round adds a scoped internal interface audit, not an originality claim."}


class Checks(unittest.TestCase):
    def test_01_spectrum_norm_and_product_variance(self):
        for n in range(1, 6):
            q = q_eigenvalues(n)
            self.assertEqual(float(max(abs(q))), 1.)
            for mu in (-.7, 0., .3, 1.):
                probabilities = np.diag(product_geometry(n, mu)).real
                self.assertAlmostEqual(float(probabilities@q), mu)
                self.assertAlmostEqual(float(probabilities@(q*q)-mu*mu), (1-mu*mu)/n)

    def test_02_full_joint_matches_product_coherence(self):
        matter = bloch(.3, -.4, .2)
        for n in range(1, 6):
            for mu in (-.6, .2):
                for pure in (False, True):
                    output = joint_evolve(np.kron(product_geometry(n, mu, pure), matter), n, .73)
                    np.testing.assert_allclose(matter_reference(output, n), channel(matter, coherence(n, mu, .73)), atol=2e-15)

    def test_03_bell_reference_attains_half_diamond(self):
        bell = density(np.array([1, 0, 0, 1], complex)/np.sqrt(2))
        for n in (1, 2, 4):
            mu, theta = .41, .82
            output = joint_evolve(np.kron(product_geometry(n, mu), bell), n, theta, 2)
            reduced = matter_reference(output, n, 2)
            mean = channel(bell, np.exp(-2j*theta*mu))
            self.assertAlmostEqual(distance(reduced, mean), channel_error(n, mu, theta))

    def test_04_variance_bound_and_endpoint_eigenstates(self):
        for n in (1, 2, 5, 17, 1000):
            for mu in np.linspace(-1, 1, 11):
                for theta in (.001, .3, 1., 3., 8.):
                    self.assertLessEqual(channel_error(n, mu, theta), min(1., theta**2*(1-mu*mu)/n)+3e-13)
        for mu in (-1., 1.):
            self.assertLess(channel_error(103, mu, 1.3), 1e-13)

    def test_05_fixed_time_leading_coefficient(self):
        rows = fixed_time_rows()
        expected = .7**2*(1-.3**2)
        self.assertLess(abs(rows[-1]["N_times_error"]-expected), 3e-5)
        self.assertTrue(all(rows[j]["half_diamond"] > rows[j+1]["half_diamond"] for j in range(len(rows)-1)))

    def test_06_sqrt_N_gaussian_nonuniform_limit(self):
        for mu in (0., .3, -.6):
            c, n = .6, 100000
            target = np.exp(-2*c*c*(1-mu*mu))
            self.assertLess(abs(centered_coherence(n, mu, c*np.sqrt(n))-target), .001)
            self.assertGreater(channel_error(n, mu, c*np.sqrt(n)), .18)

    def test_07_exact_heisenberg_conjugate_response(self):
        n, theta = 3, .91
        phases = joint_phases(n, theta)
        for j in range(n):
            y = np.kron(geometry_operator(n, j, Y), I)
            xz = np.kron(geometry_operator(n, j, X), Z)
            direct = phases.conj()[:, None]*y*phases[None, :]
            expected = np.cos(2*theta/n)*y+np.sin(2*theta/n)*xz
            np.testing.assert_allclose(direct, expected, atol=1e-15)

    def test_08_geometry_mean_and_joint_covariance(self):
        n, mu, theta = 4, .3, .7
        x = np.sqrt(1-mu*mu)
        ybar = sum(geometry_operator(n, j, Y) for j in range(n))/n
        for matter_z in (0., .6):
            output = joint_evolve(np.kron(product_geometry(n, mu), bloch(x=np.sqrt(1-matter_z**2), z=matter_z)), n, theta)
            self.assertAlmostEqual(expectation(output, np.kron(ybar, I)), x*matter_z*np.sin(2*theta/n))
            covariance = expectation(output, np.kron(ybar, Z))-expectation(output, np.kron(ybar, I))*matter_z
            self.assertAlmostEqual(covariance, x*(1-matter_z**2)*np.sin(2*theta/n))

    def test_09_global_pure_distance_and_entanglement(self):
        for n in range(1, 5):
            theta = .7
            initial = tensor_power(PLUS, n+1)
            final = joint_phases(n, theta)*initial
            rho, target = density(final), density(initial)
            self.assertAlmostEqual(distance(rho, target), np.sqrt(1-np.cos(theta/n)**(2*n)))
            self.assertAlmostEqual(negativity(rho, n), np.sqrt(1-abs(coherence(n, 0., theta))**2)/2)
            self.assertLessEqual(distance(rho, target), theta/np.sqrt(n)+1e-14)

    def test_10_same_mean_cat_channel_does_not_converge(self):
        mu, theta, matter = .3, .61, density(PLUS)
        factor = np.cos(2*theta)-1j*mu*np.sin(2*theta)
        for n in range(1, 6):
            geometry = cat_geometry(n, mu)
            variance, _, _ = covariance_identity(geometry, n)
            self.assertAlmostEqual(variance, 1-mu*mu)
            output = joint_evolve(np.kron(geometry, matter), n, theta)
            np.testing.assert_allclose(matter_reference(output, n), channel(matter, factor), atol=1e-15)

    def test_11_cat_and_classical_mixture_have_different_joint_account(self):
        n, theta = 4, np.pi/4
        outputs = [joint_evolve(np.kron(cat_geometry(n, 0., coherent), density(PLUS)), n, theta) for coherent in (True, False)]
        np.testing.assert_allclose(matter_reference(outputs[0], n), matter_reference(outputs[1], n), atol=1e-15)
        self.assertAlmostEqual(negativity(outputs[0], n), .5)
        self.assertAlmostEqual(negativity(outputs[1], n), 0.)

    def test_12_cat_logical_conjugate_observable(self):
        n, mu, mz, theta = 4, .3, .6, .71
        logical_y = np.zeros((2**n, 2**n), complex)
        logical_y[0, -1], logical_y[-1, 0] = -1j, 1j
        output = joint_evolve(np.kron(cat_geometry(n, mu), bloch(z=mz)), n, theta)
        self.assertAlmostEqual(expectation(output, np.kron(logical_y, I)), np.sqrt(1-mu*mu)*mz*np.sin(2*theta))

    def test_13_covariance_identity_and_nonshrinking_correlations(self):
        n, mu = 5, .3
        for geometry in (product_geometry(n, mu), cat_geometry(n, mu), cat_geometry(n, mu, False)):
            direct, covariance_sum, _ = covariance_identity(geometry, n)
            self.assertAlmostEqual(direct, covariance_sum)
        _, _, covariance = covariance_identity(cat_geometry(n, mu), n)
        np.testing.assert_allclose(covariance, np.full((n, n), 1-mu*mu), atol=2e-15)

    def test_14_response_time_scales_are_distinct(self):
        n = 100000
        self.assertAlmostEqual(n*.6*np.sin(2*.7/n), .84, places=9)
        self.assertLess(.6*np.sin(2*.6/np.sqrt(n)), .003)
        self.assertAlmostEqual(.6*np.sin(2*(.3*n)/n), .6*np.sin(.6))
        self.assertGreater(channel_error(n, 0., .3*n), .499)

    def test_15_closed_joint_energy_and_finite_recurrence(self):
        for n in range(1, 5):
            phases = joint_phases(n, np.pi*n)
            np.testing.assert_allclose(phases, np.full(2**(n+1), (-1)**n), atol=4e-15)
            initial = np.kron(product_geometry(n, .3), bloch(.3, .2, .4))
            h = np.diag(np.kron(q_eigenvalues(n), [1., -1.]))
            self.assertAlmostEqual(expectation(initial, h), expectation(joint_evolve(initial, n, .7), h))

    def test_16_channel_is_affine_for_fixed_geometry(self):
        rho, sigma, weight = bloch(.2, .3, -.4), bloch(-.5, .1, .2), .37
        factor = coherence(31, .3, .7)
        np.testing.assert_allclose(channel(weight*rho+(1-weight)*sigma, factor), weight*channel(rho, factor)+(1-weight)*channel(sigma, factor), atol=1e-15)

    def test_17_arbitrary_reference_never_exceeds_exact_channel_bound(self):
        rng = np.random.default_rng(353)
        n, mu, theta = 3, .3, .7
        for r in (1, 2, 3):
            matrix = rng.normal(size=(2*r, 2*r))+1j*rng.normal(size=(2*r, 2*r))
            rho = matrix@matrix.conj().T
            rho /= np.trace(rho)
            exact = channel(rho, coherence(n, mu, theta))
            mean = channel(rho, np.exp(-2j*theta*mu))
            self.assertLessEqual(distance(exact, mean), channel_error(n, mu, theta)+1e-14)
            for geometry in (product_geometry(n, mu, False), cat_geometry(n, mu)):
                initial = np.kron(geometry, rho)
                joint = joint_evolve(initial, n, theta, r)
                target = np.kron(geometry, mean)
                variance, _, _ = covariance_identity(geometry, n)
                self.assertLessEqual(distance(joint, target), abs(theta)*np.sqrt(variance)+2e-14)


if __name__ == "__main__":
    main(__name__, "collective_geometry_limit_audit", report)
