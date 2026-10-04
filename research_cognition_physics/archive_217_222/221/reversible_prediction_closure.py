"""Round221: reversible predictive closure as a conditional bridge to L.

Finite matrix certificates for a general quotient/observability argument.
Complex matrices here are declared benchmarks, not a reconstruction premise.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np


def hermitian_basis(d):
    basis = []
    for i in range(d):
        h = np.zeros((d, d), complex)
        h[i, i] = 1
        basis.append(h)
    for i in range(d):
        for j in range(i+1, d):
            h = np.zeros((d, d), complex)
            h[i, j] = h[j, i] = 1/np.sqrt(2)
            basis.append(h)
            h = np.zeros((d, d), complex)
            h[i, j], h[j, i] = -1j/np.sqrt(2), 1j/np.sqrt(2)
            basis.append(h)
    return basis


def coordinates(h, basis):
    return np.array([np.trace(b @ h).real for b in basis])


def report_matrix(effects, basis):
    return np.array([coordinates(e, basis) for e in effects])


def projector(v):
    return np.outer(v, v.conj())


def tomography_vectors(d):
    eye = np.eye(d, dtype=complex)
    vectors = list(eye)
    for i in range(d):
        for j in range(i+1, d):
            vectors.extend([(eye[i]+eye[j])/np.sqrt(2),
                            (eye[i]+1j*eye[j])/np.sqrt(2)])
    return vectors


def gate_to_vector(v):
    # All chosen vectors have real, nonnegative first component.
    e0 = np.eye(len(v), dtype=complex)[0]
    delta = e0-v
    norm_sq = np.vdot(delta, delta).real
    if norm_sq < 1e-14:
        return np.eye(len(v), dtype=complex)
    return np.eye(len(v), dtype=complex)-2*projector(delta)/norm_sq


def pullback_matrix(unitary, basis):
    # State-coordinate evolution x' = R x.
    return np.column_stack([coordinates(unitary @ b @ unitary.conj().T, basis)
                            for b in basis])


def quotient_lift(t, r):
    # A lift exists on the linear span iff t r vanishes on ker(t).
    lift = t @ r @ np.linalg.pinv(t)
    residual = t @ r-lift @ t
    return lift, float(np.linalg.norm(residual, ord=2))


def finite_orbit(d):
    basis = hermitian_basis(d)
    vectors = tomography_vectors(d)
    seed = projector(np.eye(d, dtype=complex)[0])
    gates = [gate_to_vector(v) for v in vectors]
    effects = [u.conj().T @ seed @ u for u in gates]
    return basis, gates, effects, report_matrix(effects, basis)


def retained_ability_witness():
    # Exact arithmetic: the old Z report is retained, but adding H exposes phase.
    old_report_plus = (Fraction(1, 2), Fraction(1, 2))
    old_report_minus = old_report_plus
    return {'old_reports': [list(map(str, old_report_plus)), list(map(str, old_report_minus))],
            'future_yes_probabilities': ['1', '0'],
            'minimum_uniform_scalar_forecast_error': '1/2',
            'old_Z_measurement_still_available': True}


class ReversibleClosureTests(unittest.TestCase):
    def test_orthonormal_coordinate_system(self):
        for d in (2, 3, 4):
            basis = hermitian_basis(d)
            gram = report_matrix(basis, basis)
            np.testing.assert_allclose(gram, np.eye(d*d), atol=1e-13)

    def test_one_pure_probe_has_a_finite_separating_orbit(self):
        for d in (2, 3, 4):
            _, gates, effects, t = finite_orbit(d)
            self.assertEqual(np.linalg.matrix_rank(t, tol=1e-11), d*d)
            for v, u, e in zip(tomography_vectors(d), gates, effects):
                np.testing.assert_allclose(u.conj().T @ u, np.eye(d), atol=1e-13)
                np.testing.assert_allclose(e, projector(v), atol=1e-13)
                self.assertAlmostEqual(np.trace(e).real, 1, places=12)

    def test_complete_report_has_correct_closed_updates(self):
        rng = np.random.default_rng(221)
        for d in (2, 3, 4):
            basis, gates, _, t = finite_orbit(d)
            for u in gates:
                r = pullback_matrix(u, basis)
                lift, residual = quotient_lift(t, r)
                self.assertLess(residual, 1e-11)
                g = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
                rho = g @ g.conj().T
                rho /= np.trace(rho)
                x = coordinates(rho, basis)
                np.testing.assert_allclose(lift @ (t @ x), t @ (r @ x), atol=1e-11)

    def test_old_ability_preservation_does_not_imply_report_closure(self):
        basis = hermitian_basis(2)
        z_reports = report_matrix([np.diag([1, 0]), np.diag([0, 1])], basis)
        hadamard = np.array([[1, 1], [1, -1]], complex)/np.sqrt(2)
        _, residual = quotient_lift(z_reports, pullback_matrix(hadamard, basis))
        self.assertAlmostEqual(residual, 1, places=12)
        plus, minus = np.array([1, 1])/np.sqrt(2), np.array([1, -1])/np.sqrt(2)
        xp, xm = coordinates(projector(plus), basis), coordinates(projector(minus), basis)
        np.testing.assert_allclose(z_reports @ xp, z_reports @ xm, atol=1e-13)
        future = z_reports @ pullback_matrix(hadamard, basis)
        np.testing.assert_allclose(future @ xp, [1, 0], atol=1e-13)
        np.testing.assert_allclose(future @ xm, [0, 1], atol=1e-13)
        # The allowed family still contains identity and the old readout.
        self.assertLess(quotient_lift(z_reports, np.eye(4))[1], 1e-13)

    def test_closure_alone_does_not_imply_identification(self):
        basis, gates, _, _ = finite_orbit(3)
        t = report_matrix([np.eye(3)], basis)
        for u in gates:
            self.assertLess(quotient_lift(t, pullback_matrix(u, basis))[1], 1e-12)
        self.assertEqual(np.linalg.matrix_rank(t), 1)

    def test_a_sharp_probe_without_sufficient_controls_does_not_suffice(self):
        basis = hermitian_basis(3)
        diagonal_effects = [projector(v) for v in np.eye(3, dtype=complex)]
        t = report_matrix(diagonal_effects, basis)
        u = np.diag(np.exp(1j*np.array([0.0, 0.7, 1.3])))
        self.assertLess(quotient_lift(t, pullback_matrix(u, basis))[1], 1e-12)
        self.assertEqual(np.linalg.matrix_rank(t), 3)
        seed = diagonal_effects[0]
        np.testing.assert_allclose(u.conj().T @ seed @ u, seed, atol=1e-13)

    def test_closed_statistical_update_need_not_be_stochastic(self):
        basis, _, _, t = finite_orbit(2)
        h = np.array([[1, 1], [1, -1]], complex)/np.sqrt(2)
        lift, residual = quotient_lift(t, pullback_matrix(h, basis))
        self.assertLess(residual, 1e-12)
        self.assertLess(lift.min(), -0.9)
        # These are probabilities of different tests, not exclusive outcomes.
        # The lift is valid on t(C), not on an arbitrary probability simplex.

    def test_pure_probe_is_a_sufficient_certificate_not_necessary(self):
        basis = hermitian_basis(2)
        identity = np.eye(2)
        x = np.array([[0, 1], [1, 0]], complex)
        y = np.array([[0, -1j], [1j, 0]], complex)
        z = np.diag([1, -1])
        eta = 1/3
        effects = [(identity+eta*a)/2 for a in (z, -z, x, y)]
        for e in effects:
            self.assertGreater(np.linalg.eigvalsh(e).min(), 0)
        self.assertEqual(np.linalg.matrix_rank(report_matrix(effects, basis)), 4)


def report():
    dimensions = []
    for d in (2, 3, 4):
        basis, gates, _, t = finite_orbit(d)
        initial = report_matrix([projector(v) for v in np.eye(d, dtype=complex)], basis)
        dimensions.append({'benchmark_hilbert_dimension': d,
                           'initial_diagonal_report_rank': int(np.linalg.matrix_rank(initial)),
                           'separating_probe_orbit_rank': int(np.linalg.matrix_rank(t)),
                           'minimum_extra_linear_coordinates_for_this_benchmark': d*d-d,
                           'max_completed_update_residual': max(quotient_lift(t, pullback_matrix(u, basis))[1]
                                                                for u in gates)})
    return {'round': 221, 'date': '2026-09-20',
            'scope': 'Conditional bridge to L; not a new comparison that re-admits ordinary real QM.',
            'inherited': 'Round189 full framework+U+C+L selects complex state cones; round220 distinguishes prediction from single-copy simulation.',
            'analytic_bridge': {
                'D': 'T(g omega) depends only on T(omega), for every allowed reversible g on the declared whole.',
                'R': 'Reversible pullbacks of available effects separate all states.',
                'equivalence': 'Under the finite-dimensional interior-state framework, D iff g*E=E for every g; R is span(g*E)=V*; D+R implies L.',
                'Jordan_sufficient_certificate_for_R': 'A globally pure allowed effect in E, plus the inherited pure-state-transitive action and compatible invariant self-duality.',
                'pure_probe_claimed_necessary': False,
                'minimum_linear_completion': 'dim(span(g*E))-dim(E)',
                'old_ability_preservation_implies_D': False},
            'benchmark_inputs': ['standard complex matrix states and explicit allowed unitaries',
                                 'chosen restricted readout families are not complete subsystem theories',
                                 'no new accepted cognitive axiom; D and the certificate for R remain candidates'],
            'finite_matrix_certificates': dimensions,
            'old_ability_witness': retained_ability_witness(),
            'sources': [{'url': 'https://arxiv.org/pdf/quant-ph/0307127', 'locations': ['Theorems 1-2, section 2.1'],
                         'used_for': 'Known controlled observability and orbit-span method'},
                        {'url': 'https://arxiv.org/html/2306.00362', 'locations': ['Corollary 6', 'Theorem 7', 'paragraph after Theorem 8'],
                         'used_for': 'Inherited reconstruction and pure-effect/pure-state correspondence under reversible action'}],
            'original_general_observability_theorem_claimed': False,
            'cognitive_necessity_complete': False,
            'next': 'Seek an independent cognitive basis for using the member-product statistics themselves as a closed predictive description under integration; audit access to a separating probe orbit without assuming L or complex controls.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(ReversibleClosureTests))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checked.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        Path(__file__).with_name('reversible_prediction_closure_results.json').write_text(
            json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
