"""Round 350: joint quantum feedback versus non-affine mean-field closure.

G is a finite register called geometry only as an interface placeholder.
No metric, Einstein equation, or physical gravitational coupling is derived.
"""
import unittest
import numpy as np
from growing_stream_audit import main


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
ZZ = np.kron(Z, Z)


def bloch(x=0., y=0., z=0.):
    return (I+x*X+y*Y+z*Z)/2


def rotation(theta):
    return np.diag(np.exp(-1j*theta*np.array([1., -1.])))


def joint_unitary(theta):
    return np.cos(theta)*np.eye(4)-1j*np.sin(theta)*ZZ


def evolve(rho, theta):
    unitary = joint_unitary(theta)
    return unitary@rho@unitary.conj().T


def partial(rho, dimensions, keep):
    tensor = rho.reshape(tuple(dimensions)*2)
    remaining = list(range(len(dimensions)))
    for label in reversed(range(len(dimensions))):
        if label not in keep:
            index = remaining.index(label)
            tensor = np.trace(tensor, axis1=index, axis2=index+len(remaining))
            remaining.remove(label)
    dimension = int(np.prod([dimensions[i] for i in remaining]))
    return tensor.reshape(dimension, dimension)


def distance(first, second):
    return float(np.sum(abs(np.linalg.eigvalsh(first-second)))/2)


def reduced_channel(rho, mean, theta):
    plus, minus = rotation(theta), rotation(-theta)
    return ((1+mean)/2*(plus@rho@plus.conj().T) +
            (1-mean)/2*(minus@rho@minus.conj().T))


def mean_channel(rho, mean, theta):
    unitary = rotation(theta*mean)
    return unitary@rho@unitary.conj().T


def coherence_factor(mean, theta):
    return np.cos(2*theta)-1j*mean*np.sin(2*theta)


def exact_channel_distance(mean, theta):
    return float(abs(coherence_factor(mean, theta)-np.exp(-2j*theta*mean))/2)


def nonlinear_mean_map(rho, theta):
    mean = float(np.trace(Z@rho).real)
    return mean_channel(rho, mean, theta)


def partial_transpose_second(rho):
    return rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)


def steering_branches(projectors):
    bell_vector = np.array([1, 0, 0, 1])/np.sqrt(2)
    bell = np.outer(bell_vector, bell_vector)
    return [np.einsum("abcd,ca->bd", bell.reshape(2, 2, 2, 2), effect)
            for effect in projectors]


def ensemble_case(s=.6, theta=None):
    if theta is None:
        theta = np.pi/(4*s)
    x = np.sqrt(1-s*s)
    tilted = [bloch(x=x, z=s), bloch(x=-x, z=-s)]
    vertical = [bloch(z=1), bloch(z=-1)]
    out_tilted = sum(nonlinear_mean_map(rho, theta) for rho in tilted)/2
    out_vertical = sum(nonlinear_mean_map(rho, theta) for rho in vertical)/2
    unconditioned = nonlinear_mean_map(I/2, theta)
    effect = (I+Y)/2
    return {"s": s, "theta": float(theta),
            "initial_mixture_distance": distance(sum(tilted)/2, sum(vertical)/2),
            "branch_output_distance": distance(out_tilted, out_vertical),
            "analytic_branch_output_distance": float(abs(x*np.sin(2*theta*s))/2),
            "tilted_Y_plus_probability": float(np.trace(effect@out_tilted).real),
            "vertical_Y_plus_probability": float(np.trace(effect@out_vertical).real),
            "whole_mixture_Y_plus_probability":
                float(np.trace(effect@unconditioned).real)}


class Checks(unittest.TestCase):
    def test_01_fixed_joint_unitarity_group_and_energy(self):
        for theta in (-.3, .2, 1.1):
            unitary = joint_unitary(theta)
            self.assertLess(np.linalg.norm(unitary.conj().T@unitary-np.eye(4)),
                            1e-13)
            self.assertLess(np.linalg.norm(unitary.conj().T@ZZ@unitary-ZZ), 1e-13)
        self.assertLess(np.linalg.norm(joint_unitary(.3)@joint_unitary(.7)-
                                       joint_unitary(1.)), 1e-13)

    def test_02_geometry_conjugate_observables_have_backreaction(self):
        theta = .37
        unitary = joint_unitary(theta)
        actual_x = unitary.conj().T@np.kron(X, I)@unitary
        actual_y = unitary.conj().T@np.kron(Y, I)@unitary
        self.assertLess(np.linalg.norm(actual_x - np.cos(2*theta)*np.kron(X, I) +
            np.sin(2*theta)*np.kron(Y, Z)), 1e-13)
        self.assertLess(np.linalg.norm(actual_y - np.cos(2*theta)*np.kron(Y, I) -
            np.sin(2*theta)*np.kron(X, Z)), 1e-13)
        self.assertLess(np.linalg.norm(unitary.conj().T@np.kron(Z, I)@unitary-
                                       np.kron(Z, I)), 1e-13)

    def test_03_matter_can_change_geometry_Y_reading(self):
        for theta in (.1, .4, np.pi/4):
            outputs = [partial(evolve(np.kron(bloch(x=1), bloch(z=z)), theta),
                               (2, 2), (0,)) for z in (1, -1)]
            self.assertAlmostEqual(distance(*outputs), abs(np.sin(2*theta)),
                                   places=13)
            self.assertAlmostEqual(float(np.trace(Y@outputs[0]).real),
                                   float(np.sin(2*theta)), places=13)

    def test_04_reduced_channel_independent_of_geometry_coherences(self):
        matter = bloch(.4, -.3, .2)
        for mean in (-.8, 0., .7):
            for x in (0., np.sqrt(1-mean*mean)):
                geometry = bloch(x=x, z=mean)
                output = partial(evolve(np.kron(geometry, matter), .4), (2, 2), (1,))
                self.assertLess(np.linalg.norm(output-reduced_channel(matter, mean, .4)),
                                1e-13)

    def test_05_CPTP_kraus_completeness_and_positive_Choi(self):
        bell = np.outer([1, 0, 0, 1], [1, 0, 0, 1])/2
        for mean in (-1., -.4, 0., .8, 1.):
            operators = [np.sqrt((1+mean)/2)*rotation(.7),
                         np.sqrt((1-mean)/2)*rotation(-.7)]
            self.assertLess(np.linalg.norm(sum(k.conj().T@k for k in operators)-I),
                            1e-13)
            choi = sum(np.kron(k, I)@bell@np.kron(k.conj().T, I) for k in operators)
            self.assertGreater(np.linalg.eigvalsh(choi)[0], -1e-13)
            self.assertLess(np.linalg.norm(partial(choi, (2, 2), (1,))-I/2), 1e-13)

    def test_06_uniform_variance_bound_for_mean_approximation(self):
        for mean in np.linspace(-1., 1., 31):
            for theta in np.linspace(-2., 2., 41):
                exact = exact_channel_distance(mean, theta)
                self.assertLessEqual(exact, min(1., theta**2*(1-mean*mean))+2e-14)

    def test_07_reference_input_saturates_exact_half_diamond_distance(self):
        bell = np.outer([1, 0, 0, 1], [1, 0, 0, 1])/2
        for mean, theta in ((0., .3), (.4, .7), (-.8, 1.2)):
            initial = np.kron(bloch(z=mean), bell)
            full_u = np.kron(joint_unitary(theta), I)
            exact = partial(full_u@initial@full_u.conj().T, (2, 2, 2), (1, 2))
            mean_u = np.kron(rotation(theta*mean), I)
            approximate = mean_u@bell@mean_u.conj().T
            self.assertAlmostEqual(distance(exact, approximate),
                                   exact_channel_distance(mean, theta), places=13)

    def test_08_global_correlations_are_first_order_local_error_second_order(self):
        initial = np.kron(bloch(x=1), bloch(x=1))
        for theta in (.03, .2, .5):
            exact = evolve(initial, theta)
            self.assertAlmostEqual(distance(exact, initial), abs(np.sin(theta)),
                                   places=13)
            self.assertAlmostEqual(distance(partial(exact, (2, 2), (1,)),
                                           bloch(x=1)), np.sin(theta)**2, places=13)
            self.assertAlmostEqual(float(np.trace(exact@np.kron(Y, Z)).real),
                                   np.sin(2*theta), places=13)

    def test_09_zero_geometry_variance_makes_mean_replacement_exact(self):
        matter = bloch(.3, .4, -.2)
        for mean in (-1., 1.):
            initial = np.kron(bloch(z=mean), matter)
            for theta in (.2, .7, 1.4):
                self.assertLess(np.linalg.norm(evolve(initial, theta) -
                    np.kron(bloch(z=mean), mean_channel(matter, mean, theta))), 1e-13)
                self.assertLess(exact_channel_distance(mean, theta), 1e-13)

    def test_10_same_matter_channel_does_not_fix_joint_entanglement(self):
        outputs = [evolve(np.kron(geometry, bloch(x=1)), np.pi/4)
                   for geometry in (bloch(x=1), I/2)]
        self.assertLess(np.linalg.norm(partial(outputs[0], (2, 2), (1,)) -
                                       partial(outputs[1], (2, 2), (1,))), 1e-13)
        self.assertLess(np.linalg.eigvalsh(partial_transpose_second(outputs[0]))[0],
                        -.49)
        self.assertGreater(np.linalg.eigvalsh(partial_transpose_second(outputs[1]))[0],
                           -1e-13)

    def test_11_maximal_entanglement_and_reduced_dephasing(self):
        output = evolve(np.kron(bloch(x=1), bloch(x=1)), np.pi/4)
        self.assertAlmostEqual(float(np.trace(output@output).real), 1., places=13)
        for keep in ((0,), (1,)):
            self.assertLess(np.linalg.norm(partial(output, (2, 2), keep)-I/2), 1e-13)

    def test_12_nonlinear_branch_rule_fails_affinity(self):
        case = ensemble_case()
        self.assertLess(case["initial_mixture_distance"], 1e-13)
        self.assertAlmostEqual(case["branch_output_distance"], .4, places=13)
        self.assertAlmostEqual(case["tilted_Y_plus_probability"], .9, places=13)
        self.assertAlmostEqual(case["vertical_Y_plus_probability"], .5, places=13)
        self.assertAlmostEqual(case["whole_mixture_Y_plus_probability"], .5, places=13)

    def test_13_pointwise_reversible_group_is_not_linear_channel(self):
        rho = bloch(.5, .3, .4)
        transformed = nonlinear_mean_map(rho, .7)
        self.assertLess(np.linalg.norm(nonlinear_mean_map(transformed, -.7)-rho), 1e-13)
        self.assertLess(np.linalg.norm(nonlinear_mean_map(
            nonlinear_mean_map(rho, .3), .4)-transformed), 1e-13)
        self.assertLess(np.linalg.norm(np.linalg.eigvalsh(rho)-
                                       np.linalg.eigvalsh(transformed)), 1e-13)
        self.assertGreater(ensemble_case()["branch_output_distance"], .3)

    def test_14_remote_steering_is_checked_without_communicating_outcome(self):
        s, x = .6, .8
        for ensemble in ([bloch(x=x, z=s), bloch(x=-x, z=-s)],
                         [bloch(z=1), bloch(z=-1)]):
            branches = steering_branches([rho.T for rho in ensemble])
            self.assertLess(np.linalg.norm(sum(branches)-I/2), 1e-13)
            for branch, target in zip(branches, ensemble):
                self.assertLess(np.linalg.norm(branch-target/2), 1e-13)

    def test_15_fixed_quantum_channel_preserves_ensemble_equivalence(self):
        theta = np.pi/(4*.6)
        ensembles = ([bloch(x=.8, z=.6), bloch(x=-.8, z=-.6)],
                     [bloch(z=1), bloch(z=-1)])
        for mean in (-.8, 0., .6):
            outputs = [sum(reduced_channel(rho, mean, theta) for rho in states)/2
                       for states in ensembles]
            self.assertLess(distance(*outputs), 1e-13)

    def test_16_initial_correlations_preclude_a_channel_of_matter_marginal_alone(self):
        plus = (np.eye(4)+np.kron(Z, X))/4
        minus = (np.eye(4)-np.kron(Z, X))/4
        for state in (plus, minus):
            self.assertGreater(np.linalg.eigvalsh(state)[0], -1e-13)
            for keep in ((0,), (1,)):
                self.assertLess(np.linalg.norm(partial(state, (2, 2), keep)-I/2), 1e-13)
        outputs = [partial(evolve(state, .4), (2, 2), (1,)) for state in (plus, minus)]
        self.assertAlmostEqual(distance(*outputs), abs(np.sin(.8)), places=13)


def report():
    mean_rows = [{"geometry_mean": mean, "geometry_variance": 1-mean*mean,
                  "theta": theta,
                  "exact_half_diamond_distance": exact_channel_distance(mean, theta),
                  "variance_upper_bound": min(1., theta*theta*(1-mean*mean))}
                 for mean, theta in ((0., .1), (0., .2), (.8, .2), (.8, .7), (1., .7))]
    return {"round": 350,
            "scope": "Finite G tensor M register with a specified ZZ coupling; "
                     "G is not a derived spacetime metric. Fixed independent "
                     "geometry gives a CPTP matter channel; mean closure and "
                     "pre-existing correlations require separate treatment.",
            "mean_approximation": mean_rows,
            "correlation_comparison_theta_0_2": {
                "joint_trace_distance_from_zero_mean_product": float(np.sin(.2)),
                "matter_trace_distance_from_zero_mean_unitary": float(np.sin(.2)**2),
                "Y_geometry_Z_matter_covariance": float(np.sin(.4))},
            "maximally_entangling_theta": float(np.pi/4),
            "maximal_entanglement_negativity": .5,
            "ensemble_affinity_counterexample": ensemble_case(),
            "precorrelated_equal_marginals_theta_0_4": {
                "matter_output_distance": float(abs(np.sin(.8))),
                "initial_geometry_and_matter_marginals": "I/2 on both subsystems"},
            "interpretation_limits": [
                "Z_G itself is conserved; its conjugate X_G/Y_G responds",
                "CPTP reduction assumes fixed geometry initially independent of MR",
                "Nonlinear deterministic branch law is incompatible with affine mixtures",
                "Applying the nonlinear formula once to an unlabeled density "
                "matrix is a different rule from evolving its pure-state branches",
                "No theorem excluding all semiclassical approximations or "
                "stochastic models with linear averaged channels",
                "No Einstein dynamics, geometric constraints, or gravitational "
                "coupling constant derived"]}


if __name__ == "__main__":
    main(__name__, "quantum_geometry_backreaction_audit", report)
