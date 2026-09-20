"""Round220: statistical identification versus single-copy LOCC simulation.

Standard complex QM is an explicit benchmark, not a derived theory. The note
contains universal proofs; finite NumPy checks verify the stated constructions.
Run without --write-results to preserve the stored result file.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
PAULI = (X, Y, Z)


def phi_projector(d):
    if d < 2:
        raise ValueError('The nontrivial bipartite benchmark requires d >= 2')
    v = np.eye(d).reshape(-1) / np.sqrt(d)
    return np.outer(v, v).astype(complex)


def partial_transpose_b(matrix, d):
    return matrix.reshape(d, d, d, d).transpose(0, 3, 2, 1).reshape(d*d, d*d)


def product_ic_effects():
    # Bob uses conjugate tetrahedron vectors, so equal labels have weight 5/2.
    vectors = np.array([[1, 1, 1], [1, -1, -1],
                        [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
    local = [(I2 + sum(v[k]*PAULI[k] for k in range(3)))/4 for v in vectors]
    return [np.kron(a, b.conj()) for a in local for b in local]


def hermitian_coordinates(matrix):
    basis = [np.kron(a, b) for a in (I2, *PAULI) for b in (I2, *PAULI)]
    return np.array([np.trace(h @ matrix).real/4 for h in basis])


def ic_reconstruction():
    effects = product_ic_effects()
    coefficients = np.column_stack([hermitian_coordinates(e) for e in effects])
    weights = np.linalg.solve(coefficients, hermitian_coordinates(phi_projector(2)))
    recovered = sum(w*e for w, e in zip(weights, effects))
    return coefficients, weights, recovered


def pauli_records():
    # Shared classical choice of axis, then both parties measure locally.
    records = []
    for axis, correlation_sign in zip(PAULI, (1, -1, 1)):
        for a in (-1, 1):
            for b in (-1, 1):
                effect = np.kron((I2+a*axis)/2, (I2+b*axis)/2)/3
                match = correlation_sign*a*b == 1
                records.append((effect, match))
    return records


def pauli_binary_effect():
    # Keep a matching outcome with probability 3/4, and reject every other one.
    return sum((0.75 if match else 0)*e for e, match in pauli_records())


def minimax_effect(d):
    return (np.eye(d*d) + d*phi_projector(d))/(d+2)


def op_norm(h):
    return float(np.max(np.abs(np.linalg.eigvalsh(h))))


def random_density(rng, d):
    g = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    rho = g @ g.conj().T
    return rho / np.trace(rho)


def probability(effect, rho):
    return float(np.trace(effect @ rho).real)


def finite_report_decision_success(effects, states, prior):
    # Independent enumeration of the optimal binary decision for each report.
    return sum(max(prior[i]*probability(e, states[i]) for i in range(2))
               for e in effects)


class PredictionSimulationTests(unittest.TestCase):
    def test_product_povm_is_informationally_complete(self):
        effects = product_ic_effects()
        np.testing.assert_allclose(sum(effects), np.eye(4), atol=1e-13)
        for e in effects:
            self.assertGreaterEqual(np.linalg.eigvalsh(e).min(), -1e-13)
        matrix, _, _ = ic_reconstruction()
        self.assertEqual(np.linalg.matrix_rank(matrix, tol=1e-12), 16)

    def test_unique_statistical_weights_are_not_stochastic(self):
        _, weights, recovered = ic_reconstruction()
        expected = np.array([2.5 if i == j else -0.5 for i in range(4) for j in range(4)])
        np.testing.assert_allclose(weights, expected, atol=1e-12)
        np.testing.assert_allclose(recovered, phi_projector(2), atol=1e-13)
        self.assertLess(weights.min(), 0)
        self.assertGreater(weights.max(), 1)

    def test_signed_prediction_for_general_complex_states(self):
        rng = np.random.default_rng(220)
        _, weights, _ = ic_reconstruction()
        effects = product_ic_effects()
        for _ in range(30):
            rho = random_density(rng, 4)
            estimate = sum(w*probability(e, rho) for w, e in zip(weights, effects))
            self.assertAlmostEqual(estimate, probability(phi_projector(2), rho), places=12)

    def test_explicit_local_protocol_attains_qubit_bound(self):
        records = pauli_records()
        np.testing.assert_allclose(sum(e for e, _ in records), np.eye(4), atol=1e-13)
        m = pauli_binary_effect()
        np.testing.assert_allclose(m, minimax_effect(2), atol=1e-13)
        self.assertAlmostEqual(op_norm(m-phi_projector(2)), 0.25, places=12)
        for accept in (m, np.eye(4)-m):
            self.assertGreaterEqual(np.linalg.eigvalsh(partial_transpose_b(accept, 2)).min(), -1e-13)

    def test_general_dimension_primal_and_witness_certificates(self):
        for d in (2, 3, 4, 5):
            p = phi_projector(d)
            m = minimax_effect(d)
            eps = float(Fraction(1, d+2))
            # The analytic lower bound uses W=I/d-P, nonnegative on product PSDs.
            witness = np.eye(d*d)/d-p
            self.assertAlmostEqual(np.trace(witness @ m).real, 0, places=12)
            self.assertAlmostEqual(op_norm(m-p), eps, places=12)
            self.assertAlmostEqual(probability(m, p), 1-eps, places=12)
            q = (np.eye(d*d)-p)/(d*d-1)
            self.assertAlmostEqual(probability(m, q), eps, places=12)
            self.assertGreaterEqual(np.linalg.eigvalsh(partial_transpose_b(m, d)).min(), -1e-12)
            self.assertGreaterEqual(np.linalg.eigvalsh(partial_transpose_b(np.eye(d*d)-m, d)).min(), -1e-12)

    def test_witness_on_independently_generated_product_positive_matrices(self):
        rng = np.random.default_rng(2201)
        for d in (2, 3, 4):
            w = np.eye(d*d)/d-phi_projector(d)
            for _ in range(15):
                a, b = random_density(rng, d), random_density(rng, d)
                self.assertGreaterEqual(np.trace(w @ np.kron(a, b)).real, -1e-12)
        # Samples check the implementation, not the universal separability proof.

    def test_blackwell_task_witness_even_with_complete_statistics(self):
        p = phi_projector(2)
        q = (np.eye(4)-p)/3
        prior = (0.25, 0.75)
        for effects in (product_ic_effects(), [e for e, _ in pauli_records()],
                        [pauli_binary_effect(), np.eye(4)-pauli_binary_effect()]):
            self.assertAlmostEqual(finite_report_decision_success(effects, (p, q), prior), 0.75, places=12)
        self.assertAlmostEqual(finite_report_decision_success([p, np.eye(4)-p], (p, q), prior), 1, places=12)

    def test_unbiased_estimator_and_finite_sample_variance(self):
        rng = np.random.default_rng(2202)
        p = phi_projector(2)
        for rho in [p, (np.eye(4)-p)/3, np.eye(4)/4] + [random_density(rng, 4) for _ in range(20)]:
            mean = sum((1 if match else -0.5)*probability(e, rho) for e, match in pauli_records())
            second = sum((1 if match else 0.25)*probability(e, rho) for e, match in pauli_records())
            true_p = probability(p, rho)
            self.assertAlmostEqual(mean, true_p, places=12)
            self.assertAlmostEqual(second-mean*mean, (1+true_p)/2-true_p**2, places=12)
            self.assertLessEqual(second-mean*mean, 9/16+1e-12)


def report():
    matrix, weights, reconstructed = ic_reconstruction()
    p = phi_projector(2)
    q = (np.eye(4)-p)/3
    return {
        'round': 220,
        'date': '2026-09-20',
        'purpose': 'Audit a proposed bridge from round219 report translation to round189 L without discarding established reconstruction constraints.',
        'inputs': ['standard complex quantum benchmark', 'one unknown bipartite input per simulation',
                   'local instruments and classical communication, no extra shared entanglement or quantum communication',
                   'independent repeated preparation only for the separate estimation task'],
        'inherited_result': 'Local tomography and implementability of all joint measurements differ (round83). Full framework+U+C+L reconstruction remains valid.',
        'analytic_results': {
            'uniform_single_copy_probability_error': '1/(d+2), d>=2',
            'attaining_effect': '(I+d*P_Phi)/(d+2)',
            'lower_bound': 'Separable M implies <Phi|M|Phi> <= Tr(M)/d',
            'finite_task_prior': '[1/(d+2), (d+1)/(d+2)] on [P_Phi,(I-P_Phi)/(d^2-1)]',
            'finite_task_optimal_local_success': '(d+1)/(d+2)',
            'finite_task_global_success': '1',
            'qubit_estimator_values': ['1', '-1/2'],
            'qubit_estimator_variance': '(1+p)/2-p^2 <= 9/16',
            'iid_mean_squared_error_bound': '9/(16*N); not a single-copy measurement simulation'
        },
        'matrix_checks': {
            'informationally_complete_product_povm_rank': int(np.linalg.matrix_rank(matrix, tol=1e-12)),
            'unique_dual_weights_rounded': np.round(weights, 12).tolist(),
            'reconstruction_residual': float(np.linalg.norm(reconstructed-p)),
            'qubit_protocol_uniform_error': op_norm(pauli_binary_effect()-p),
            'ic_local_task_success': finite_report_decision_success(product_ic_effects(), (p, q), (0.25, 0.75)),
            'global_task_success': finite_report_decision_success([p, np.eye(4)-p], (p, q), (0.25, 0.75))
        },
        'source': {'url': 'https://arxiv.org/html/1901.09772', 'locations': ['section III.3, equations (35)-(37)', 'section IV.1'],
                   'used_for': 'Known local verification operator (I+dP)/(d+1); balanced rescaling and task witness derived explicitly in the note.'},
        'original_quantum_theorem_claimed': False,
        'cognitive_derivation_of_L_complete': False,
        'next': 'Motivate distributed statistical predictive completeness L at the level of allowed experiments, without requiring stochastic realization of every global measurement or counting repeated preparations as free.'
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(PredictionSimulationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        Path(__file__).with_name('local_prediction_vs_simulation_results.json').write_text(
            json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
