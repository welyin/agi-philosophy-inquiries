"""Round 317: nested quantum algebras, relative-entropy loss and response calibration."""
import unittest
import numpy as np
from growing_stream_audit import main
from modular_energy_interface_audit import budget, centered, entropy, gibbs, modular_hamiltonian


I = np.eye(2)
X = np.array([[0., 1.], [1., 0.]])
Y = np.array([[0., -1j], [1j, 0.]])
Z = np.diag([1., -1.])


def tensor(a, b, c):
    return np.kron(np.kron(a, b), c)


def discard_first(rho):
    n = len(rho)
    if n % 2 or n < 4:
        raise ValueError('Require a two-dimensional first factor and a nontrivial remainder.')
    return np.trace(np.asarray(rho).reshape(2, n//2, 2, n//2), axis1=0, axis2=2)


def reference():
    # A supplied noncommuting finite reference generator; not a spacetime boost.
    h = (0.7*tensor(Z, I, I)+0.4*tensor(I, X, I)+0.2*tensor(I, I, Z)
         +0.3*tensor(X, X, I)+0.2*tensor(I, Z, Z)+0.1*tensor(Y, I, X))
    return gibbs(h)


def contrast_tangent(sigma):
    reduced = discard_first(sigma)
    contrast = centered(modular_hamiltonian(sigma)-np.kron(I, modular_hamiltonian(reduced)))
    tangent = contrast*(0.3*np.linalg.eigvalsh(sigma)[0]/np.linalg.norm(contrast, 2))
    return contrast, tangent


def step_balance(rho, sigma, alpha=1.):
    before = budget(rho, sigma)
    after = budget(discard_first(rho), discard_first(sigma))
    entropy_change = after['delta_S']-before['delta_S']
    modular_loss = before['delta_K']-after['delta_K']
    relative_loss = before['relative_entropy']-after['relative_entropy']
    increment = entropy_change+alpha*modular_loss
    return {'alpha': alpha, 'entropy_difference_change': entropy_change,
            'modular_contrast': modular_loss, 'relative_entropy_loss': relative_loss,
            'candidate_generalized_increment': increment,
            'decomposition_error': abs(increment-relative_loss-(alpha-1)*modular_loss)}


def nested_rows(rho, sigma):
    rows = []
    for label in ('ABC', 'BC', 'C'):
        row = budget(rho, sigma)
        rows.append({'algebra': label, **row, 'matched_entropy_potential': row['delta_S']-row['delta_K']})
        if label != 'C':
            rho, sigma = discard_first(rho), discard_first(sigma)
    return rows


def admissible_interval(epsilon, sigma, tangent):
    positive = step_balance(sigma+epsilon*tangent, sigma)
    negative = step_balance(sigma-epsilon*tangent, sigma)
    kappa = positive['modular_contrast']
    if kappa <= 0:
        raise ValueError('This control uses a positively oriented modular contrast.')
    return {'epsilon': epsilon,
            'minimum_alpha': 1-positive['relative_entropy_loss']/kappa,
            'maximum_alpha': 1+negative['relative_entropy_loss']/kappa,
            'positive_loss': positive['relative_entropy_loss'],
            'negative_loss': negative['relative_entropy_loss'], 'modular_contrast': kappa}


def report():
    sigma = reference(); contrast, tangent = contrast_tangent(sigma)
    rows = nested_rows(sigma+0.8*tangent, sigma)
    bad = [step_balance(sigma+eps*tangent, sigma, alpha) | {'epsilon': eps}
           for alpha, eps in ((0.8, 0.1), (1.2, -0.1))]
    return {'round': 317,
            'scope': 'Finite full-rank three-qubit reference and nested partial traces. Data processing yields a conditional generalized-entropy balance after explicitly assigning the area-response candidate. No physical horizon, Einstein equation or Newton constant is derived.',
            'reference_min_eigenvalue': float(np.linalg.eigvalsh(sigma)[0]),
            'tangent_noncommutativity': float(np.linalg.norm(sigma@tangent-tangent@sigma)),
            'contrast_norm': float(np.linalg.norm(contrast)),
            'nested_algebra_rows': rows,
            'matched_finite_increment': step_balance(sigma+0.8*tangent, sigma),
            'wrong_coefficient_counterexamples': bad,
            'two_sided_coefficient_intervals': [admissible_interval(eps, sigma, tangent) for eps in (0.8, 0.4, 0.2, 0.1)],
            'logical_boundary': 'Wall 1105.3445 assumes semiclassical Einstein gravity. Its GSL is a forward consistency result here, not an independent derivation of the same field equation. Two-sided near-reference monotonicity only constrains a response contrast.',
            'next_interface': 'Compare canonical boost energy with improved Hilbert stress and the nonminimal Wald surface term on the same null cut.'}


class Checks(unittest.TestCase):
    def test_01_partial_trace_preserves_dual_expectations(self):
        sigma = reference(); rng = np.random.default_rng(317)
        raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
        observable = raw+raw.conj().T
        self.assertAlmostEqual(float(np.trace(sigma@np.kron(I, observable)).real),
                               float(np.trace(discard_first(sigma)@observable).real), places=12)

    def test_02_nested_relative_entropy_is_monotone(self):
        sigma = reference(); rng = np.random.default_rng(31702)
        for _ in range(12):
            raw = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
            rho = raw@raw.conj().T; rho /= np.trace(rho)
            rows = nested_rows(rho, sigma)
            self.assertTrue(all(a['relative_entropy'] >= b['relative_entropy']-1e-12 for a, b in zip(rows, rows[1:])))

    def test_03_matched_increment_is_relative_entropy_loss(self):
        sigma = reference(); _, tangent = contrast_tangent(sigma)
        for eps in (-0.8, 0.1, 0.8):
            row = step_balance(sigma+eps*tangent, sigma)
            self.assertLess(row['decomposition_error'], 2e-14)
            self.assertGreater(row['candidate_generalized_increment'], 0.)

    def test_04_saturation_when_no_distinguishability_is_discarded(self):
        fixed = np.diag([0.4, 0.6]); sigmab = np.diag([0.7, 0.3]); rhob = np.diag([0.3, 0.7])
        row = step_balance(np.kron(fixed, rhob), np.kron(fixed, sigmab))
        self.assertLess(abs(row['relative_entropy_loss']), 1e-13)

    def test_05_different_reference_breaks_the_claim(self):
        sigma = reference(); marginal = discard_first(sigma)
        wrong = np.eye(4)/4
        # Input rho=sigma has D=0; replacing the output reference gives D>0.
        self.assertGreater(budget(marginal, wrong)['relative_entropy'], 0.01)

    def test_06_wrong_area_coefficient_fails_for_one_perturbation_sign(self):
        for row in report()['wrong_coefficient_counterexamples']:
            self.assertLess(row['candidate_generalized_increment'], -1e-5)
            self.assertLess(row['decomposition_error'], 2e-14)

    def test_07_two_sided_band_shrinks_to_one(self):
        rows = report()['two_sided_coefficient_intervals']
        widths = [r['maximum_alpha']-r['minimum_alpha'] for r in rows]
        for r in rows:
            self.assertLess(r['minimum_alpha'], 1.)
            self.assertGreater(r['maximum_alpha'], 1.)
        for first, second in zip(widths, widths[1:]):
            self.assertTrue(1.9 < first/second < 2.1)

    def test_08_operator_response_condition_is_independent_of_basis(self):
        sigma = reference(); contrast, tangent = contrast_tangent(sigma)
        h = 1e-3; alpha = 1.2
        derivative = (step_balance(sigma+h*tangent, sigma, alpha)['candidate_generalized_increment']
                      -step_balance(sigma-h*tangent, sigma, alpha)['candidate_generalized_increment'])/(2*h)
        expected = float((alpha-1)*np.trace(tangent@contrast).real)
        self.assertLess(abs(derivative-expected), 1e-9)
        self.assertGreater(np.linalg.norm(sigma@tangent-tangent@sigma), 1e-5)

    def test_09_entropy_alone_has_no_monotonicity_under_discard(self):
        bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
        pure = np.outer(bell, bell)
        self.assertGreater(entropy(discard_first(pure))-entropy(pure), 0.6)
        mixed = np.eye(4)/4
        self.assertLess(entropy(discard_first(mixed))-entropy(mixed), -0.6)

    def test_10_coefficient_is_unidentifiable_for_uniform_references(self):
        sigma = np.eye(4)/4; rho = np.diag([0.4, 0.1, 0.2, 0.3])
        values = [step_balance(rho, sigma, alpha) for alpha in (-2., 1., 5.)]
        self.assertTrue(all(abs(v['modular_contrast']) < 1e-13 for v in values))
        self.assertLess(max(v['candidate_generalized_increment'] for v in values)-min(v['candidate_generalized_increment'] for v in values), 1e-13)


if __name__ == '__main__':
    main(__name__, 'nested_entropy_balance_audit', report)
