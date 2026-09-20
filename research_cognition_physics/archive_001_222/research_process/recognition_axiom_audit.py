"""Round 218: quantify the four proposed recognition axioms.

Finite examples check the independent analytic arguments in the note.
Overlap and pure-state entanglement are declared example scores, not derived
definitions of cognition. No classification theorem is proved by a scan.
"""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import unittest

import numpy as np


def state_dimension(d, field):
    if type(d) is not int or d < 1 or field not in ('classical', 'real', 'complex'):
        raise ValueError('Positive integer dimension and a declared field required.')
    return {'classical': d - 1, 'real': d * (d + 1) // 2 - 1,
            'complex': d * d - 1}[field]


def normalized_gaps(energies):
    values = sorted(Fraction(x) for x in energies)
    width = values[-1] - values[0]
    if width <= 0:
        raise ValueError('A nonconstant spectrum is required.')
    return [(values[i + 1] - values[i]) / width for i in range(len(values) - 1)]


def orthogonal_flow(generator, time):
    generator = np.asarray(generator, dtype=float)
    if not np.allclose(generator + generator.T, 0, atol=1e-13, rtol=0):
        raise ValueError('A real antisymmetric generator is required.')
    values, vectors = np.linalg.eigh(1j * generator)
    result = (vectors * np.exp(-1j * time * values)) @ vectors.conj().T
    if np.max(np.abs(result.imag)) > 1e-12:
        raise ArithmeticError('Unexpected imaginary residual in the real flow.')
    return result.real


def marginal(rho, m, n, keep):
    tensor = rho.reshape(m, n, m, n)
    if keep == 'A':
        return np.trace(tensor, axis1=1, axis2=3)
    if keep == 'B':
        return np.trace(tensor, axis1=0, axis2=2)
    raise ValueError('Keep A or B.')


def schmidt_path(m, n, theta):
    """Real normalized path; theta=acos(1/sqrt(r)) gives r equal coefficients."""
    r = min(m, n)
    if r < 2:
        raise ValueError('Both factors must have at least two levels.')
    coefficients = np.zeros((m, n))
    coefficients[0, 0] = np.cos(theta)
    for i in range(1, r):
        coefficients[i, i] = np.sin(theta) / np.sqrt(r - 1)
    return coefficients.reshape(-1)


def pure_entanglement_score(psi, m, n):
    """Normalized linear entropy, only a declared score for bipartite pure states."""
    if not np.isclose(np.vdot(psi, psi), 1):
        raise ValueError('Normalized pure vector required.')
    rho = np.outer(psi, psi.conj())
    reduced = marginal(rho, m, n, 'A')
    return float((1 - np.trace(reduced @ reduced).real) / (1 - 1 / min(m, n)))


def squared_overlap(psi, phi):
    return float(abs(np.vdot(psi, phi)) ** 2)


class RecognitionAxiomAuditTests(unittest.TestCase):
    def test_finite_embedding_dimension_obstruction(self):
        # Examples illustrate, rather than prove, the general monotonicity argument.
        for field in ('classical', 'real', 'complex'):
            for d in (2, 3, 5):
                self.assertGreater(state_dimension(d * d, field), state_dimension(d, field))

    def test_spectral_obstruction_survives_clock_rescaling(self):
        a, b = normalized_gaps([0, 1, 3]), normalized_gaps([0, 1, 4])
        self.assertNotEqual(a, b)
        self.assertNotEqual(a, list(reversed(b)))
        self.assertEqual(normalized_gaps([5, 7, 11]), a)

    def test_enlarged_real_register_preserves_dynamics_and_readout(self):
        rng = np.random.default_rng(21801)
        for m, n in ((2, 2), (2, 3), (3, 4)):
            raw_a, raw_b = rng.normal(size=(m, m)), rng.normal(size=(n, n))
            ka, kb = raw_a - raw_a.T, raw_b - raw_b.T
            raw = rng.normal(size=(n, n))
            rho = raw @ raw.T
            rho /= np.trace(rho)
            oa, ob = orthogonal_flow(ka, .37), orthogonal_flow(kb, .37)
            encoded = np.kron(np.eye(m) / m, rho)
            oab = np.kron(oa, ob)
            evolved = oab @ encoded @ oab.T
            expected = np.kron(np.eye(m) / m, ob @ rho @ ob.T)
            np.testing.assert_allclose(evolved, expected, atol=1e-12, rtol=0)
            np.testing.assert_allclose(marginal(encoded, m, n, 'B'), rho, atol=1e-12)
            np.testing.assert_allclose(marginal(evolved, m, n, 'B'), ob @ rho @ ob.T, atol=1e-12)

    def test_continuous_real_entangled_paths(self):
        for m, n in ((2, 2), (2, 3), (3, 4)):
            endpoint = np.arccos(1 / np.sqrt(min(m, n)))
            scores = [pure_entanglement_score(schmidt_path(m, n, t), m, n)
                      for t in np.linspace(0, endpoint, 19)]
            self.assertAlmostEqual(scores[0], 0)
            self.assertAlmostEqual(scores[-1], 1)
            self.assertGreaterEqual(min(np.diff(scores)), -1e-12)

    def test_real_entanglement_stable_under_local_controls(self):
        rng = np.random.default_rng(21802)
        for m, n in ((2, 2), (2, 3), (3, 4)):
            psi = schmidt_path(m, n, .31)
            ka, kb = rng.normal(size=(m, m)), rng.normal(size=(n, n))
            evolved = np.kron(orthogonal_flow(ka - ka.T, .4),
                              orthogonal_flow(kb - kb.T, -.7)) @ psi
            self.assertAlmostEqual(pure_entanglement_score(psi, m, n),
                                   pure_entanglement_score(evolved, m, n))

    def test_interaction_can_both_create_and_remove_entanglement(self):
        k = np.zeros((4, 4))
        k[3, 0], k[0, 3] = 1, -1
        start = np.eye(4)[0]
        for t in (0, .17, np.pi / 4, 1.1, np.pi / 2):
            psi = orthogonal_flow(k, t) @ start
            self.assertAlmostEqual(pure_entanglement_score(psi, 2, 2), np.sin(2 * t) ** 2)
            np.testing.assert_allclose(orthogonal_flow(-k, t) @ psi, start, atol=1e-12)

    def test_overlap_common_motion_and_independent_drift(self):
        k = np.array([[0., -1.], [1., 0.]])
        psi, phi = np.array([1., 0.]), np.array([np.cos(.3), np.sin(.3)])
        for t in (.2, .9, 1.4):
            o = orthogonal_flow(k, t)
            self.assertAlmostEqual(squared_overlap(o @ psi, o @ phi), squared_overlap(psi, phi))
        self.assertAlmostEqual(squared_overlap(psi, orthogonal_flow(k, np.pi / 2) @ psi), 0)

    def test_real_ghz_entanglement_for_every_bipartition(self):
        # The analytic Schmidt decomposition covers all n; this checks n <= 6.
        for n in range(2, 7):
            psi = np.zeros(2 ** n)
            psi[0] = psi[-1] = 1 / np.sqrt(2)
            for mask in range(1, 2 ** n - 1):
                left = [i for i in range(n) if mask & (1 << i)]
                right = [i for i in range(n) if i not in left]
                matrix = psi.reshape((2,) * n).transpose(left + right).reshape(2 ** len(left), -1)
                singular = np.linalg.svd(matrix, compute_uv=False)
                np.testing.assert_allclose(singular[:2], [1/np.sqrt(2)] * 2, atol=1e-12)
                np.testing.assert_allclose(singular[2:], 0, atol=1e-12)


def report():
    return {
        'round': 218, 'date': '2026-09-20',
        'motivation': 'Evaluate four user-proposed recognition axioms before adding further axioms.',
        'assumptions': ['finite affine state spaces for the dimension no-go',
                        'injective affine embeddings preserve operational distinctions',
                        'fixed subsystem identification and declared evolution class',
                        'squared overlap and normalized pure-state linear entropy are example scores only'],
        'analytic_results': {
            'all_pairs_fixed_host_embeddings_force_equal_finite_dimensions': True,
            'strict_dimension_growth_then_contradicts_all_pairs_embedding': True,
            'monotonicity_under_every_group_element_and_its_inverse_implies_invariance': True,
            'scalar_continuity_implies_full_pure_state_reversible_transitivity': False,
            'more_joint_coordinates_alone_implies_entanglement': False,
            'original_axioms_proved_to_select_complex_quantum_theory': False},
        'illustrations': {
            'qubit_affine_dimension': 3, 'qutrit_affine_dimension': 8,
            'two_qubits_affine_dimension': 15,
            'classical_two_bits_joint_affine_dimension': 3,
            'classical_two_bits_sum_of_marginal_affine_dimensions': 2,
            'incompatible_normalized_spectral_gaps': [
                [str(x) for x in normalized_gaps([0, 1, 3])],
                [str(x) for x in normalized_gaps([0, 1, 4])]],
            'times': [0, float(np.pi/4), float(np.pi/2)],
            'overlap_under_independent_real_rotation': [1., .5, 0.],
            'pure_entanglement_under_joint_real_rotation': [0., 1., 0.]},
        'repaired_contract': {
            'encoding': 'E(rho_B)=I_A/d_A tensor rho_B; D=partial_trace_A',
            'intertwining': 'T_AB(t) E = E T_B(t) for product orthogonal dynamics',
            'resource_cost': 'A register of dimension d_A*d_B and the declared independent preparation',
            'real_pure_entanglement_paths_and_local_stability_exist': True,
            'unbounded_party_real_GHZ_family_exists': True,
            'arbitrary_complex_collective_cognition_formalized_or_proved': False,
            'user_has_accepted_revisions': False},
        'reconstruction_dependencies_remaining': ['operational states, effects and transformations',
                                                'Jordan/HSD and factorization or a justified alternative',
                                                'non-signaling compatible associative composition',
                                                'local tomography',
                                                'complex-qubit premise for the 2014 route'],
        'next': 'Define recognition by declared prediction/control tasks and separate passive retention, active learning and noise recovery; then assess whether that task entails local tomography.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = unittest.TextTestRunner(verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(RecognitionAxiomAuditTests))
    if not result.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['numerical_checks'] = {'run': result.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        Path(__file__).with_name('recognition_axiom_audit_results.json').write_text(
            json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
