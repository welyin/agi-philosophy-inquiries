"""Round 204: a classical competition model with genuine Kahler geometry.

The geometry describes explicitly stored candidate weights in the open
three-category simplex. An ensemble of such registers is not its mean weight.
Fitness, metric and orientation are declared choices, not cognitive axioms.
No quantum theory, all-SoCA completeness, or physical time is derived here.
"""

import argparse
from fractions import Fraction as F
import json
from math import sqrt
from pathlib import Path
import unittest

import numpy as np
import predictive_geometry_audit as old_geometry


E = np.array([[F(0), F(1)], [F(-1), F(0)]], dtype=object)
I = np.array([[F(1), F(0)], [F(0), F(1)]], dtype=object)
P = (F(1, 2), F(1, 3), F(1, 6))
Q = (F(1, 6), F(1, 2), F(1, 3))


def interior(values):
    p = tuple(F(v) for v in values)
    if len(p) < 2 or any(v <= 0 for v in p) or sum(p) != 1:
        raise ValueError("Use an interior normalized finite probability vector.")
    return p


def select(p, weights):
    p, weights = interior(p), tuple(F(v) for v in weights)
    if len(p) != len(weights) or any(v <= 0 for v in weights):
        raise ValueError("One positive selection weight per candidate required.")
    normalizer = sum(a*b for a, b in zip(p, weights))
    return tuple(a*b/normalizer for a, b in zip(p, weights))


def replicator(p, fitness):
    p, fitness = interior(p), tuple(F(v) for v in fitness)
    if len(p) != len(fitness):
        raise ValueError("One fitness per candidate required.")
    mean = sum(a*b for a, b in zip(p, fitness))
    return tuple(a*(b-mean) for a, b in zip(p, fitness))


def cyclic_fitness(p):
    x, y, z = interior(p)
    return z-y, x-z, y-x


def cyclic_field(p):
    return replicator(p, cyclic_fitness(p))


def fisher(p):
    x, y, z = interior(p)
    return np.array([[1/x+1/z, 1/z], [1/z, 1/y+1/z]], dtype=object)


def inverse_fisher(p):
    x, y, z = interior(p)
    return np.array([[x*(1-x), -x*y], [-x*y, y*(1-y)]], dtype=object)


def complex_structure_divided_by_sqrt_product(p):
    """J/sqrt(x*y*z), rational at rational p; omega(u,v)=g(Ju,v)."""
    x, y, z = interior(p)
    return -(inverse_fisher(p) @ E)/(x*y*z)


def lie_metric(p):
    """Exact Lie derivative of Fisher metric along the cyclic field in x,y."""
    x, y, z = interior(p)
    vx, vy, _ = cyclic_field(p)
    dg_dx = np.array([[-1/x**2+1/z**2, 1/z**2],
                      [1/z**2, 1/z**2]], dtype=object)
    dg_dy = np.array([[1/z**2, 1/z**2],
                      [1/z**2, -1/y**2+1/z**2]], dtype=object)
    derivative = np.array([[2*z-1, -2*x], [2*y, 1-2*z]], dtype=object)
    g = fisher(p)
    return vx*dg_dx+vy*dg_dy+derivative.T @ g+g @ derivative


def stereographic_data(p):
    """Explicit complex chart w=(sqrt(x)+i sqrt(y))/(1+sqrt(z))."""
    x, y, z = (float(v) for v in interior(p))
    a, b, c = sqrt(x), sqrt(y), sqrt(z)
    d = 1+c
    w = complex(a, b)/d
    jac = np.array([[1/(2*a*d)+a/(2*c*d*d), a/(2*c*d*d)],
                    [b/(2*c*d*d), 1/(2*b*d)+b/(2*c*d*d)]])
    metric = 16/(1+abs(w)**2)**2*(jac.T @ jac)
    return w, jac, metric


def rational_grid(denominator=9):
    return [(F(i, denominator), F(j, denominator), F(denominator-i-j, denominator))
            for i in range(1, denominator-1) for j in range(1, denominator-i)]


class CompetitionGeometryAuditTests(unittest.TestCase):
    def test_unrestricted_positive_selection_can_reach_any_interior_target(self):
        for p in rational_grid(7):
            for target in (P, Q, (F(1, 3),)*3):
                weights = tuple(b/a for a, b in zip(p, target))
                self.assertEqual(select(p, weights), target)
                # Bounding all scores by 1 does not change the normalized map.
                self.assertEqual(select(p, tuple(w/max(weights) for w in weights)), target)

    def test_unrestricted_fitness_represents_any_given_tangent(self):
        for p in rational_grid():
            for v in ((F(1), F(-2), F(1)), (p[1], -p[0], p[0]-p[1])):
                self.assertEqual(replicator(p, tuple(a/b for a, b in zip(v, p))), v)

    def test_finite_selection_increment_has_the_declared_continuous_generator(self):
        fitness, h = (F(-1), F(2), F(3)), F(1, 100)
        following = select(P, tuple(1+h*f for f in fitness))
        mean = sum(p*f for p, f in zip(P, fitness))
        self.assertEqual(tuple((b-a)/h for a, b in zip(P, following)),
                         tuple(v/(1+h*mean) for v in replicator(P, fitness)))

    def test_fixed_fitness_is_fisher_gradient_with_variance_growth(self):
        fitness = (F(-1), F(0), F(2))
        for p in rational_grid():
            v = replicator(p, fitness)
            covector = np.array([fitness[0]-fitness[2], fitness[1]-fitness[2]], dtype=object)
            np.testing.assert_array_equal(fisher(p) @ np.array(v[:2], dtype=object), covector)
            mean = sum(a*b for a, b in zip(p, fitness))
            variance = sum(a*(b-mean)**2 for a, b in zip(p, fitness))
            self.assertEqual(sum(a*b for a, b in zip(v, fitness)), variance)

    def test_cyclic_field_preserves_mass_and_product_infinitesimally(self):
        for p in rational_grid():
            v = cyclic_field(p)
            self.assertEqual(sum(v), 0)
            self.assertEqual(sum(a/b for a, b in zip(v, p)), 0)
        self.assertNotEqual(cyclic_field(P), (F(0),)*3)
        self.assertEqual(cyclic_field((F(1, 3),)*3), (F(0),)*3)

    def test_fisher_determinant_and_compatible_complex_structure_exactly(self):
        for p in rational_grid():
            x, y, z = p
            product = x*y*z
            g = fisher(p)
            scaled_j = complex_structure_divided_by_sqrt_product(p)
            self.assertEqual(g[0, 0]*g[1, 1]-g[0, 1]*g[1, 0], 1/product)
            np.testing.assert_array_equal(g @ inverse_fisher(p), I)
            np.testing.assert_array_equal(product*(scaled_j @ scaled_j), -I)
            np.testing.assert_array_equal(product*(scaled_j.T @ g @ scaled_j), g)
            np.testing.assert_array_equal(g @ scaled_j, -E/product)

    def test_explicit_complex_chart_matches_fisher_metric(self):
        for p in rational_grid():
            w, jac, metric = stereographic_data(p)
            self.assertTrue(w.real > 0 and w.imag > 0 and abs(w) < 1)
            self.assertGreater(np.linalg.det(jac), 0)
            np.testing.assert_allclose(metric, np.asarray(fisher(p), dtype=float), rtol=3e-15)

    def test_hamiltonian_identity_for_fisher_area_form(self):
        # i_X omega=d(2 sqrt(x*y*z)); multiply by sqrt(product) to test rationally.
        for p in rational_grid():
            x, y, z = p
            differential_log_product = np.array([1/x-1/z, 1/y-1/z], dtype=object)
            v = np.array(cyclic_field(p)[:2], dtype=object)
            np.testing.assert_array_equal(-E @ v, x*y*z*differential_log_product)

    def test_intrinsic_area_is_not_the_old_projective_pullback(self):
        _, _, metric, symplectic = old_geometry.ray_data(np.asarray(P, dtype=float))
        np.testing.assert_allclose(metric, np.asarray(fisher(P), dtype=float)/4)
        np.testing.assert_array_equal(symplectic, np.zeros((2, 2)))
        self.assertNotEqual(1/sqrt(float(P[0]*P[1]*P[2])), 0)

    def test_hamiltonian_flow_is_not_fisher_isometric(self):
        np.testing.assert_array_equal(lie_metric(P),
                                      np.array([[F(-4, 3), F(-1)], [F(-1), F(0)]], dtype=object))

    def test_average_register_cannot_replace_a_randomized_register_preparation(self):
        midpoint = tuple((a+b)/2 for a, b in zip(P, Q))
        h = F(1, 10)
        def step(p):
            return select(p, tuple(1+h*f for f in cyclic_fitness(p)))
        averaged_outputs = tuple((a+b)/2 for a, b in zip(step(P), step(Q)))
        midpoint_output = step(midpoint)
        self.assertEqual(tuple(a-b for a, b in zip(midpoint_output, averaged_outputs)),
                         (F(0), F(1, 480), F(-1, 480)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=1).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(CompetitionGeometryAuditTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    results = {
        "round": 204,
        "scope": "Declared candidate-weight selection submodule; not all SoCA features",
        "antecedent_rounds": [2, 181, 189, 192, 193, 195, 203],
        "new_question": "Does explicit competition dynamics plus genuine compatible Kahler geometry suffice to infer quantum operational structure?",
        "answer": False,
        "classical_register_values": "interior triples p with sum p=1; a mixture of register values is not their mean",
        "replicator_generator": "p_i*(f_i-sum_j p_j*f_j)",
        "unrestricted_interior_fitness_represents_any_tangent_field": True,
        "fisher_metric_cognitively_unique": False,
        "declared_orientation_and_fitness_are_inputs": True,
        "cyclic_field": ["x*(z-y)", "y*(x-z)", "z*(y-x)"],
        "metric_xy": [["1/x+1/z", "1/z"], ["1/z", "1/y+1/z"]],
        "metric_determinant": "1/(x*y*z)",
        "symplectic_form": "dx wedge dy / sqrt(x*y*z)",
        "complex_coordinate": "(sqrt(x)+i*sqrt(y))/(1+sqrt(z))",
        "metric_in_complex_chart": "16*abs(dw)^2/(1+abs(w)^2)^2",
        "hamiltonian_for_same_field_and_fisher_area": "2*sqrt(x*y*z)",
        "integrability_argument": "Explicit complex chart on the open spherical octant; not inferred from J^2=-I samples",
        "point": [str(v) for v in P],
        "lie_derivative_fisher_at_point": [[str(v) for v in row] for row in lie_metric(P)],
        "flow_is_fisher_isometric": False,
        "old_real_projective_pullback_symplectic_form": "zero; different from the newly chosen intrinsic area form",
        "poisson_bracket_x_y": "sqrt(x*y*z); not affine in x,y",
        "midpoint_vs_mixture_one_step_gap_h_one_tenth": ["0", "1/480", "-1/480"],
        "physical_time_derived": False,
        "quantum_structure_derived": False,
        "full_soca_completeness_refuted_or_proved": False,
        "counterexample_to_U_C_L_reconstruction": False,
        "new_tests": checked.testsRun,
    }
    if args.write_results:
        Path(__file__).with_name("competition_geometry_audit_results.json").write_text(
            json.dumps(results, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
