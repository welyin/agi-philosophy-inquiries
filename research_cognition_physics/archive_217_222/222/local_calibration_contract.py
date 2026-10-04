"""Round 222: finite witnesses for local process calibration and tomography.

The general implication is proved in the note. Matrices are declared benchmarks.
Reuse the existing Python/NumPy runtime and round 221 Hermitian basis.
"""
import argparse
import json
from pathlib import Path
import unittest

import numpy as np

from reversible_prediction_closure import hermitian_basis, coordinates


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
J = np.array([[0, -1], [1, 0]], dtype=complex)


def finite_contract(states, effects, basis):
    """Construct two physical channels with the same isolated state map.

    Columns of S are normalized spanning states; rows of E are an IC POVM.
    S a E = identity. Negative a entries are algebraic coefficients only.
    All actual measurement-to-preparation probabilities w are nonnegative.
    """
    k = len(basis)
    s = np.column_stack([coordinates(x, basis) for x in states])
    e = np.array([coordinates(x, basis) for x in effects])
    if s.shape != (k, k) or e.shape != (k, k):
        raise ValueError('Need square, linearly independent state and effect frames')
    a = np.linalg.solve(s, np.linalg.inv(e))
    q = np.full(k, 1/k)
    thresholds = [1.0]
    for j in range(k):
        for i in range(k):
            if a[j, i] < q[j]:
                thresholds.append(q[j]/(q[j]-a[j, i]))
    p = min(thresholds)/2
    w = p*a+(1-p)*q[:, None]
    outputs = [sum(w[j, i]*states[j] for j in range(k)) for i in range(k)]
    reset_state = sum(q[j]*states[j] for j in range(k))
    return dict(p=p, a=a, w=w, outputs=outputs, reset=reset_state,
                algebraic_residual=float(np.linalg.norm(s @ a @ e-np.eye(k))))


def interior_frames(d):
    """An explicit interior state basis and an implementable IC POVM."""
    basis = hermitian_basis(d)
    traceless = []
    for j in range(d-1):
        h = np.zeros((d, d), complex)
        h[j, j], h[-1, -1] = 1, -1
        traceless.append(h)
    traceless += basis[d:]
    k = d*d
    identity = np.eye(d)
    states = [identity/d]+[identity/d+h/(4*d) for h in traceless]
    tail = [identity/k+h/(2*k*k) for h in traceless]
    effects = [identity-sum(tail)]+tail
    return states, effects, basis


def measure_prepare(value, effects, outputs):
    return sum(np.trace(e @ value)*s for e, s in zip(effects, outputs))


def conditional_b(value, effect, da, db):
    return np.einsum('ca,abcd->bd', effect, value.reshape(da, db, da, db))


def extended_mp(value, effects, outputs, da, db):
    return sum(np.kron(s, conditional_b(value, e, da, db))
               for e, s in zip(effects, outputs))


def extended_noisy_identity(value, reset, p, db):
    da = len(reset)
    marginal = conditional_b(value, np.eye(da), da, db)
    return p*value+(1-p)*np.kron(reset, marginal)


def choi(action, d):
    out = np.zeros((d*d, d*d), complex)
    for j in range(d):
        for k in range(d):
            unit = np.zeros((d, d), complex)
            unit[j, k] = 1
            out += np.kron(unit, action(unit))
    return out


def axis_contract(axes):
    """Measure a uniformly random Pauli axis, then prepare its eigenstate."""
    states = [(I+sign*a)/2 for a in axes for sign in (1, -1)]
    effects = [s/len(axes) for s in states]
    return effects, states, 1/len(axes)


def bell_state():
    v = np.array([1, 0, 0, 1], complex)/np.sqrt(2)
    return np.outer(v, v.conj())


def real_axis_witness():
    effects, states, p = axis_contract([X, Z])
    isolated_gap = max(np.linalg.norm(measure_prepare(h, effects, states)-
                       (p*h+(1-p)*np.trace(h)*I/2)) for h in (I, X, Z))
    rho = bell_state()
    noisy = extended_noisy_identity(rho, I/2, p, 2)
    measured = extended_mp(rho, effects, states, 2, 2)
    readout = (np.eye(4)+np.kron(J, J))/2
    gap = noisy-measured
    return dict(p=p, isolated_max_residual=float(isolated_gap),
                noisy_yes=float(np.trace(readout @ noisy).real),
                measured_yes=float(np.trace(readout @ measured).real),
                trace_distance=float(np.linalg.norm(gap, ord='nuc')/2),
                hidden_gap_formula_residual=float(np.linalg.norm(gap-p*np.kron(J, J)/4)),
                output_min_eigenvalues=[float(np.linalg.eigvalsh(v).min()) for v in (noisy, measured)])


def pure_preservation_constraint_rank(d):
    """K v must be parallel to v for basis vectors and pairwise sums.

    Positivity of Kraus terms makes these necessary for a channel fixing all
    pure inputs. Their solution is K = scalar * identity (proved for all d).
    """
    eye = np.eye(d)
    vectors = list(eye)
    vectors += [(eye[j]+eye[k])/np.sqrt(2) for j in range(d) for k in range(j+1, d)]
    constraints = np.vstack([np.kron(v.reshape(1, d), eye-np.outer(v, v)) for v in vectors])
    return int(np.linalg.matrix_rank(constraints)), float(np.linalg.norm(constraints @ eye.flatten(order='F')))


class CalibrationTests(unittest.TestCase):
    def test_general_constructor_positive_frames_and_probabilities(self):
        for d in (2, 3):
            states, effects, basis = interior_frames(d)
            c = finite_contract(states, effects, basis)
            self.assertGreater(min(np.linalg.eigvalsh(v).min() for v in states+effects+c['outputs']), 0)
            self.assertGreater(c['p'], 0)
            self.assertGreaterEqual(c['w'].min(), 0)
            np.testing.assert_allclose(sum(effects), np.eye(d), atol=2e-14)
            np.testing.assert_allclose([np.trace(v) for v in states+c['outputs']], 1, atol=2e-13)
            np.testing.assert_allclose(c['a'].sum(axis=0), 1, atol=3e-11)
            np.testing.assert_allclose(c['w'].sum(axis=0), 1, atol=2e-13)

    def test_general_constructor_agrees_on_entire_complex_space(self):
        for d in (2, 3):
            states, effects, basis = interior_frames(d)
            c = finite_contract(states, effects, basis)
            for h in basis:
                np.testing.assert_allclose(measure_prepare(h, effects, c['outputs']),
                    c['p']*h+(1-c['p'])*np.trace(h)*c['reset'], atol=2e-13)

    def test_both_constructed_channels_are_completely_positive(self):
        for d in (2, 3):
            states, effects, basis = interior_frames(d)
            c = finite_contract(states, effects, basis)
            m = choi(lambda h: measure_prepare(h, effects, c['outputs']), d)
            n = choi(lambda h: c['p']*h+(1-c['p'])*np.trace(h)*c['reset'], d)
            self.assertGreater(np.linalg.eigvalsh(m).min(), -1e-12)
            np.testing.assert_allclose(m, n, atol=2e-13)

    def test_generic_complex_pair_agrees_on_correlated_inputs(self):
        rng = np.random.default_rng(222)
        for da, db in ((2, 2), (3, 2)):
            states, effects, basis = interior_frames(da)
            c = finite_contract(states, effects, basis)
            for _ in range(4):
                v = rng.normal(size=da*db)+1j*rng.normal(size=da*db)
                rho = np.outer(v, v.conj())/np.vdot(v, v)
                np.testing.assert_allclose(extended_mp(rho, effects, c['outputs'], da, db),
                    extended_noisy_identity(rho, c['reset'], c['p'], db), atol=2e-13)

    def test_axes_full_complex_channel_equality(self):
        effects, states, p = axis_contract([X, Y, Z])
        self.assertAlmostEqual(p, 1/3)
        m = choi(lambda h: measure_prepare(h, effects, states), 2)
        n = choi(lambda h: p*h+(1-p)*np.trace(h)*I/2, 2)
        np.testing.assert_allclose(m, n, atol=2e-14)
        self.assertGreater(np.linalg.eigvalsh(m).min(), -1e-14)

    def test_real_failure_is_only_in_context_for_declared_pair(self):
        r = real_axis_witness()
        self.assertLess(r['isolated_max_residual'], 1e-14)
        self.assertLess(r['hidden_gap_formula_residual'], 1e-14)
        self.assertAlmostEqual(r['noisy_yes'], 3/4)
        self.assertAlmostEqual(r['measured_yes'], 1/2)
        self.assertAlmostEqual(r['trace_distance'], 1/4)
        self.assertGreaterEqual(min(r['output_min_eigenvalues']), -1e-14)

    def test_product_statistics_miss_exact_process_difference(self):
        effects, states, p = axis_contract([X, Z])
        rho = bell_state()
        gap = extended_noisy_identity(rho, I/2, p, 2)-extended_mp(rho, effects, states, 2, 2)
        for a in (I, X, Z):
            for b in (I, X, Z):
                self.assertAlmostEqual(abs(np.trace(np.kron(a, b) @ gap)), 0)
        self.assertAlmostEqual(np.trace(np.kron(J, J) @ gap).real, p)

    def test_identity_preserving_real_kraus_lemma_finite_certificates(self):
        for d in (2, 3, 4, 5):
            rank, residual = pure_preservation_constraint_rank(d)
            self.assertEqual(rank, d*d-1)
            self.assertLess(residual, 1e-13)


def report():
    certificates = []
    for d in (2, 3):
        states, effects, basis = interior_frames(d)
        c = finite_contract(states, effects, basis)
        certificates.append(dict(complex_hilbert_dimension=d, operational_dimension=d*d,
            identity_mixture_probability=c['p'], minimum_actual_probability=float(c['w'].min()),
            algebraic_identity_residual=c['algebraic_residual'],
            max_local_map_residual=float(max(np.linalg.norm(measure_prepare(h, effects, c['outputs'])-
                (c['p']*h+(1-c['p'])*np.trace(h)*c['reset'])) for h in basis))))
    return {'round': 222, 'date': '2026-09-20',
        'scope': 'Finite-dimensional causal convex operational framework with identity, discard/prepare, conditioning and random control; no tomography presupposed.',
        'analytic_result': 'L iff context stability of all locally equal channels iff context stability of one constructed noisy-identity / measure-prepare pair per system, for every partner.',
        'quantifiers': 'Complete local input-output interface, all allowed correlated inputs and future partners. Fixed repeatable channels, not unmodelled device memory.',
        'inherited_selection': 'Round189 F+U+C+L selects complex matrix state cones. Ordinary real QM remains excluded from that branch.',
        'generic_construction_certificates': certificates,
        'concrete_real_pair': real_axis_witness(),
        'complex_axis_pair': {'identity_weight': '1/3', 'full_complex_channel_equality': True},
        'reversible_only_insufficient': {'analytic': 'A real CPTP map identical on all local states to an orthogonal channel has that same complete extension; nevertheless real QM fails L.',
            'finite_nullities': {str(d): d*d-pure_preservation_constraint_rank(d)[0] for d in (2, 3, 4, 5)}},
        'cognitive_status': 'Candidate modular-calibration principle, not derived from SoCA or adopted as a new axiom. The theorem supplies an operational test, not an independent necessity argument.',
        'sources': [{'url': 'https://arxiv.org/html/1011.6451', 'location': 'III.1.4, equation (5)', 'use': 'Known implication from local tomography to local process characterization'},
                    {'url': 'https://arxiv.org/html/1907.07043', 'location': 'Proposition 1', 'use': 'Known distinction between local action and complete transformation'}],
        'claims_original_discovery_of_process_tomography': False,
        'next': 'Keep L / modular calibration explicit; audit U quantifiers against capability preservation and remote ensemble choice without assuming a complete source or same operational dimension.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CalibrationTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        Path(__file__).with_name('local_calibration_contract_results.json').write_text(
            json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
