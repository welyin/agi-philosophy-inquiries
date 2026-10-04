"""Round 348: finite-range unitary walk and a controlled Dirac band limit.

This is a one-particle walk on a specified lattice, not a full field algebra
or a deduction of spacetime dimension, autonomous timing, or gravitation.
"""
import math
import unittest
import numpy as np
from growing_stream_audit import main


IDENTITY = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)


def coin(angle):
    return np.cos(angle)*IDENTITY-1j*np.sin(angle)*X


def walk_step(momentum, spacing, mass, hop=1):
    phases = np.exp(-1j*hop*spacing*momentum*np.array([1., -1.]))
    return np.diag(phases) @ coin(spacing*mass)


def dirac(momentum, mass, time, speed=1.):
    energy = math.hypot(speed*momentum, mass)
    if energy == 0:
        return IDENTITY.copy()
    return np.cos(time*energy)*IDENTITY - 1j*np.sin(time*energy)/energy*(
        speed*momentum*Z+mass*X)


def quasienergy(momentum, spacing, mass, hop=1):
    phase_k, phase_m = hop*spacing*momentum, spacing*mass
    cosine = np.cos(phase_k)*np.cos(phase_m)
    vector = np.array([np.cos(phase_k)*np.sin(phase_m),
                       np.sin(phase_k)*np.sin(phase_m),
                       np.sin(phase_k)*np.cos(phase_m)])
    return float(np.arctan2(np.linalg.norm(vector), cosine)/spacing)


def group_velocity(momentum, spacing, mass):
    energy = quasienergy(momentum, spacing, mass)
    denominator = np.sin(spacing*energy)
    return float(np.sin(spacing*momentum)*np.cos(spacing*mass)/denominator)


def block_diagonal(blocks):
    width = len(blocks[0])
    result = np.zeros((len(blocks)*width, len(blocks)*width), complex)
    for index, block in enumerate(blocks):
        sl = slice(index*width, (index+1)*width)
        result[sl, sl] = block
    return result


def position_walk(sites, spacing, mass, hop=1):
    shift = np.zeros((sites, sites), complex)
    for source in range(sites):
        shift[(source+hop) % sites, source] = 1
    conditional = np.kron(shift, np.diag([1., 0.])) + np.kron(
        shift.conj().T, np.diag([0., 1.]))
    return conditional @ np.kron(np.eye(sites), coin(spacing*mass))


def momentum_grid(sites, spacing):
    return 2*np.pi*np.fft.fftfreq(sites, d=spacing)


def fourier(sites, spacing):
    momenta = momentum_grid(sites, spacing)
    locations = spacing*np.arange(sites)
    return np.exp(1j*np.outer(locations, momenta))/np.sqrt(sites)


def position_from_blocks(blocks, spacing):
    transform = np.kron(fourier(len(blocks), spacing), IDENTITY)
    return transform @ block_diagonal(blocks) @ transform.conj().T


def operator_error(momentum, spacing, mass, time):
    steps = round(time/spacing)
    if not np.isclose(steps*spacing, time, atol=1e-13):
        raise ValueError("Physical comparison time must be an integer step count.")
    actual = np.linalg.matrix_power(walk_step(momentum, spacing, mass), steps)
    expected = dirac(momentum, mass, time)
    return float(np.linalg.norm(actual-expected, 2))


def pure_density(vector):
    vector = vector/np.linalg.norm(vector)
    return np.outer(vector, vector.conj())


def trace_distance(left, right):
    return float(np.sum(abs(np.linalg.eigvalsh(left-right)))/2)


def measurable_case(momentum=.8, spacing=.2, mass=.7, time=2.):
    actual = np.linalg.matrix_power(walk_step(momentum, spacing, mass),
                                    round(time/spacing))
    expected = dirac(momentum, mass, time)
    initial = np.array([1., 1j])/np.sqrt(2)
    rho, sigma = pure_density(actual@initial), pure_density(expected@initial)
    values, vectors = np.linalg.eigh(rho-sigma)
    effect = vectors[:, values > 0] @ vectors[:, values > 0].conj().T
    probabilities = [float(np.trace(effect@state).real) for state in (rho, sigma)]
    return {"momentum": momentum, "spacing": spacing, "mass": mass, "time": time,
            "operator_error": float(np.linalg.norm(actual-expected, 2)),
            "analytic_trace_distance_upper_bound": time*spacing*abs(mass*momentum),
            "chosen_state_trace_distance": trace_distance(rho, sigma),
            "same_effect_probabilities": probabilities,
            "same_effect_probability_difference": abs(probabilities[0]-probabilities[1])}


def band_case(spacing, mass=.7, cutoff=1., time=2.):
    momenta = np.linspace(-cutoff, cutoff, 401)
    errors = [operator_error(k, spacing, mass, time) for k in momenta]
    energies = [abs(quasienergy(k, spacing, mass)-math.hypot(k, mass))
                for k in momenta]
    return {"spacing": spacing, "steps": round(time/spacing),
            "mass": mass, "cutoff": cutoff, "time": time,
            "sampled_band_operator_error": max(errors),
            "analytic_uniform_operator_error_bound": time*spacing*abs(mass)*cutoff,
            "sampled_band_energy_error": max(energies)}


class Checks(unittest.TestCase):
    def test_01_exact_unitarity_and_determinant(self):
        for spacing in (.03, .2, 1.):
            for mass in (0., .7, 2.):
                for k in np.linspace(-np.pi/spacing, np.pi/spacing, 31):
                    unitary = walk_step(k, spacing, mass)
                    self.assertLess(np.linalg.norm(unitary.conj().T@unitary-
                                                   IDENTITY), 1e-13)
                    self.assertLess(abs(np.linalg.det(unitary)-1), 1e-13)

    def test_02_independent_position_and_momentum_implementations(self):
        for sites in (9, 16, 31):
            spacing, mass = .2, .7
            direct = position_walk(sites, spacing, mass)
            blocks = [walk_step(k, spacing, mass)
                      for k in momentum_grid(sites, spacing)]
            self.assertLess(np.linalg.norm(direct-position_from_blocks(
                blocks, spacing)), 2e-12)

    def test_03_strict_finite_step_support_not_only_small_tails(self):
        sites, origin, steps = 41, 20, 6
        state = np.zeros(2*sites, complex)
        state[2*origin:2*origin+2] = np.array([1., 1j])/np.sqrt(2)
        evolved = np.linalg.matrix_power(position_walk(sites, .1, 1.3),
                                         steps)@state
        probabilities = np.sum(abs(evolved.reshape(sites, 2))**2, axis=1)
        outside = abs(np.arange(sites)-origin) > steps
        self.assertEqual(float(np.sum(probabilities[outside])), 0.)
        self.assertAlmostEqual(float(sum(probabilities)), 1., places=13)

    def test_04_exact_dispersion_from_eigenvalues(self):
        for spacing, mass in ((.2, .7), (.5, 1.4)):
            for k in np.linspace(-np.pi/spacing, np.pi/spacing, 37):
                unitary = walk_step(k, spacing, mass)
                energy = quasienergy(k, spacing, mass)
                self.assertLess(abs(np.trace(unitary)/2-
                    np.cos(spacing*k)*np.cos(spacing*mass)), 1e-13)
                phases = np.sort(np.angle(np.linalg.eigvals(unitary)))
                self.assertLess(np.linalg.norm(phases -
                    np.array([-spacing*energy, spacing*energy])), 2e-13)

    def test_05_energy_squared_leading_correction(self):
        k, mass = .8, .7
        target = -(k*mass)**2/3
        residuals = []
        for spacing in (.2, .1, .05, .025):
            coefficient = (quasienergy(k, spacing, mass)**2-k*k-mass*mass)/spacing**2
            residuals.append(abs(coefficient-target))
        self.assertLess(residuals[-1], 1e-5)
        self.assertTrue(all(a > 3.8*b for a, b in zip(residuals, residuals[1:])))

    def test_06_single_step_commutator_bound(self):
        for spacing in (.03, .2, .7):
            for k, mass in ((0., .7), (.4, 0.), (.4, .7), (-1.3, 1.2)):
                error = np.linalg.norm(walk_step(k, spacing, mass) -
                                       dirac(k, mass, spacing), 2)
                self.assertLessEqual(error, spacing**2*abs(k*mass)+2e-14)

    def test_07_finite_time_uniform_band_bound(self):
        for spacing in (.2, .1, .05, .025):
            case = band_case(spacing)
            self.assertLessEqual(case["sampled_band_operator_error"],
                case["analytic_uniform_operator_error_bound"]+2e-13)
        rows = [band_case(a) for a in (.1, .05, .025)]
        for first, second in zip(rows, rows[1:]):
            self.assertGreater(first["sampled_band_operator_error"]/
                               second["sampled_band_operator_error"], 1.9)
            self.assertGreater(first["sampled_band_energy_error"]/
                               second["sampled_band_energy_error"], 3.8)

    def test_08_arbitrary_reference_on_band_limited_input(self):
        momenta, spacing, mass, time = (-1., -.3, .4, 1.), .1, .7, 2.
        actual = block_diagonal([np.linalg.matrix_power(walk_step(k, spacing, mass),
                                                       20) for k in momenta])
        expected = block_diagonal([dirac(k, mass, time) for k in momenta])
        rng = np.random.default_rng(34808)
        state = rng.normal(size=24)+1j*rng.normal(size=24)
        state /= np.linalg.norm(state)
        left = np.kron(actual, np.eye(3))@state
        right = np.kron(expected, np.eye(3))@state
        distance = trace_distance(pure_density(left), pure_density(right))
        self.assertLessEqual(distance, time*spacing*abs(mass)+2e-13)
        self.assertLess(np.linalg.norm(np.kron(actual-expected, np.eye(3)), 2)-
                        np.linalg.norm(actual-expected, 2), 2e-13)

    def test_09_finite_measurement_difference_is_operational(self):
        case = measurable_case()
        self.assertGreater(case["chosen_state_trace_distance"], .01)
        self.assertAlmostEqual(case["chosen_state_trace_distance"],
                               case["same_effect_probability_difference"], places=13)
        self.assertLessEqual(case["chosen_state_trace_distance"],
                             case["analytic_trace_distance_upper_bound"])

    def test_10_massless_limit_is_exact_at_integer_steps(self):
        for spacing in (.2, .1, .05):
            for k in np.linspace(-np.pi/spacing, np.pi/spacing, 31):
                self.assertLess(operator_error(k, spacing, 0., 2.), 2e-13)
        sites = 21
        state = np.zeros(2*sites, complex)
        state[2*10] = 1
        final = np.linalg.matrix_power(position_walk(sites, .1, 0.), 4)@state
        self.assertAlmostEqual(float(abs(final[2*14])**2), 1., places=13)

    def test_11_group_velocity_respects_exact_step_front(self):
        for spacing, mass in ((.2, .7), (.5, 1.4)):
            for k in np.linspace(-np.pi/spacing+.01, np.pi/spacing-.01, 101):
                self.assertLessEqual(abs(group_velocity(k, spacing, mass)),
                                     1.+1e-13)
                h = 1e-5
                numerical = (quasienergy(k+h, spacing, mass)-
                             quasienergy(k-h, spacing, mass))/(2*h)
                self.assertAlmostEqual(group_velocity(k, spacing, mass),
                                       numerical, places=7)

    def test_12_local_unitarity_allows_exact_two_step_oscillation(self):
        spacing = .2
        mass = np.pi/(2*spacing)  # Fixed coin angle, not the finite-mass scaling.
        for k in np.linspace(-10., 10., 31):
            unitary = walk_step(k, spacing, mass)
            self.assertLess(np.linalg.norm(unitary@unitary+IDENTITY), 1e-13)
            self.assertAlmostEqual(quasienergy(k, spacing, mass),
                                   np.pi/(2*spacing), places=13)

    def test_13_two_local_species_need_not_have_one_limiting_speed(self):
        spacing, mass, time, momentum = .02, .7, 1., .8
        for hop in (1, 2):
            actual = np.linalg.matrix_power(walk_step(momentum, spacing, mass, hop),
                                             round(time/spacing))
            expected = dirac(momentum, mass, time, speed=hop)
            self.assertLessEqual(np.linalg.norm(actual-expected, 2),
                                 time*spacing*abs(mass*hop*momentum)+2e-13)
        self.assertAlmostEqual(quasienergy(.1, spacing, 0., hop=2)/
                               quasienergy(.1, spacing, 0., hop=1), 2., places=13)
        direct_sum = block_diagonal([walk_step(momentum, spacing, mass, r)
                                    for r in (1, 2)])
        self.assertLess(np.linalg.norm(direct_sum.conj().T@direct_sum-
                                       np.eye(4)), 1e-13)

    def test_14_principal_continuous_generator_has_long_range(self):
        sites, spacing = 31, 1.
        generator = position_from_blocks([k*Z for k in momentum_grid(sites, spacing)],
                                         spacing)
        for separation in (2, 4, 9, 15):
            block = generator[2*separation:2*separation+2, :2]
            exact = np.pi/(sites*spacing*abs(np.sin(np.pi*separation/sites)))
            self.assertAlmostEqual(float(abs(block[0, 0])), exact, places=12)
        step = position_walk(sites, spacing, 0.)
        self.assertEqual(float(np.linalg.norm(step[8:10, :2])), 0.)

    def test_15_fractional_log_interpolation_loses_step_locality(self):
        sites, spacing = 31, 1.
        fractional = position_from_blocks(
            [dirac(k, 0., .5) for k in momentum_grid(sites, spacing)], spacing)
        state = np.zeros(2*sites, complex)
        state[0] = 1
        probability = abs(fractional@state).reshape(sites, 2)**2
        outside = np.minimum(np.arange(sites), sites-np.arange(sites)) > 1
        self.assertGreater(float(np.sum(probability[outside])), .05)
        integer_step = position_from_blocks(
            [dirac(k, 0., 1.) for k in momentum_grid(sites, spacing)], spacing)
        self.assertLess(np.linalg.norm(integer_step-position_walk(sites, 1., 0.)),
                        2e-12)

    def test_16_cutoff_cannot_expand_to_the_whole_brillouin_zone(self):
        errors = [operator_error(np.pi/a, a, .7, 2.) for a in (.1, .05, .025)]
        self.assertTrue(all(value > 1.2 for value in errors))
        self.assertLess(abs(errors[-1]-2*np.sin(.7)), .01)

    def test_17_same_dispersion_does_not_mean_same_finite_time_channel(self):
        k, spacing, mass = .8, .2, .7
        usual = walk_step(k, spacing, mass)
        shift = np.diag(np.exp(-1j*spacing*k*np.array([1., -1.])))
        reversed_order = coin(spacing*mass)@shift
        self.assertLess(abs(np.trace(usual)-np.trace(reversed_order)), 1e-13)
        self.assertGreater(np.linalg.norm(np.linalg.matrix_power(usual, 10)-
                                         np.linalg.matrix_power(reversed_order, 10)),
                           .05)


def report():
    sites = 31
    generator = position_from_blocks([k*Z for k in momentum_grid(sites, 1.)], 1.)
    fractional = position_from_blocks(
        [dirac(k, 0., .5) for k in momentum_grid(sites, 1.)], 1.)
    outside = np.minimum(np.arange(sites), sites-np.arange(sites)) > 1
    leaked = float(np.sum(abs(fractional[:, 0]).reshape(sites, 2)[outside]**2))
    return {"round": 348,
            "scope": "Specified one-dimensional single-particle unitary walk, "
                     "dt=a and c=hbar=1; controlled Dirac limit on a fixed "
                     "momentum band. Not a derivation of dimension, QFT or GR.",
            "literature_scope": {"classification_coin_dimension": 2,
                                 "classification_dimensions_considered": [1, 2, 3],
                                 "selects_dimension_three": False},
            "band_comparisons": [band_case(a) for a in (.2, .1, .05, .025)],
            "finite_measurement_witness": measurable_case(),
            "nonuniform_high_momentum_limit": [
                {"spacing": a, "momentum": float(np.pi/a),
                 "operator_error": operator_error(np.pi/a, a, .7, 2.)}
                for a in (.1, .05, .025)],
            "continuous_log_counterexample": {
                "sites": sites, "spacing": 1.,
                "generator_distance_four_entry_magnitude":
                    float(abs(generator[8, 0])),
                "half_step_probability_outside_one_site": leaked,
                "exact_integer_step_range": 1},
            "alternate_rules": {
                "fixed_coin_pi_over_two": "U(k)^2=-I; flat quasienergy, no net spreading",
                "two_species_integer_hops": [1, 2],
                "limiting_speed_ratio": 2},
            "added_inputs": ["integer line and cell identities",
                            "two-dimensional coin and selected conditional shifts",
                            "homogeneous fixed update rule",
                            "physical lattice spacing and dt=a calibration",
                            "coin angle a*m with finite m",
                            "fixed low-momentum band and finite comparison time"],
            "open_bridges": ["all-site field algebra and many-particle structure",
                             "autonomous implementation of step updates",
                             "selection of dimension and common physical metric",
                             "curved geometry, its dynamics and Einstein equations"]}


if __name__ == "__main__":
    main(__name__, "quantum_walk_relativistic_limit_audit", report)
