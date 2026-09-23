"""Round 391: infer common sharp input information from actual record channels.

Applies the existing correctable-algebra theorem to the round-250 interaction.
No spatial geometry or preferred input pointer observable is given to the solver.
Default read-only; --write-results creates the result file exclusively.
"""
import argparse
import json
from pathlib import Path
import platform
import unittest

import numpy as np


HERE = Path(__file__).resolve().parent
TARGET = HERE / 'shared_sharp_record_audit_results.json'
TOL = 2e-10
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)
I4 = np.eye(4, dtype=complex)
S = np.kron(Z, Z)  # Expected answer, used only in independent verification.


def interaction_isometry():
    """Two parent qubits, two blank |+> readers, all four parent-reader CZs."""
    v = np.zeros((16, 4), dtype=complex)
    for parent in range(4):
        a, b = divmod(parent, 2)
        for c in range(2):
            for d in range(2):
                v[4*parent+2*c+d, parent] = (-1.)**((a+b)*(c+d))/2
    return v


def marginal_kraus(operators, dimensions, keep):
    keep = tuple(keep)
    discard = tuple(i for i in range(len(dimensions)) if i not in keep)
    dk = int(np.prod([dimensions[i] for i in keep]))
    de = int(np.prod([dimensions[i] for i in discard]))
    din = operators[0].shape[1]
    answer = []
    for v in operators:
        tensor = v.reshape(tuple(dimensions)+(din,))
        tensor = tensor.transpose(keep+discard+(len(dimensions),)).reshape(dk, de, din)
        answer.extend(tensor[:, j, :] for j in range(de))
    return answer


def apply(operators, a):
    d = operators[0].shape[0]
    return sum((k @ a @ k.conj().T for k in operators), np.zeros((d, d), dtype=complex))


def adjoint(operators, a):
    d = operators[0].shape[1]
    return sum((k.conj().T @ a @ k for k in operators), np.zeros((d, d), dtype=complex))


def products(operators):
    return [a.conj().T @ b for a in operators for b in operators]


def commutant(generators, d, tol=1e-9):
    """Solve [A,G]=0 by complex linear equations; columns use column vec."""
    identity = np.eye(d, dtype=complex)
    system = np.vstack([np.kron(g.T, identity)-np.kron(identity, g) for g in generators])
    _, singular, vh = np.linalg.svd(system, full_matrices=False)
    rank = int(np.sum(singular > tol))
    vectors = vh.conj().T[:, rank:]
    return [vectors[:, j].reshape(d, d, order='F') for j in range(vectors.shape[1])], singular


def recoverable(operators):
    return commutant(products(operators), operators[0].shape[1])[0]


def common(channels):
    return commutant([a for k in channels for a in products(k)], channels[0][0].shape[1])[0]


def projector(basis):
    v = np.column_stack([a.reshape(-1, order='F') for a in basis])
    q, _ = np.linalg.qr(v)
    return q @ q.conj().T


def span_error(a, basis):
    v = a.reshape(-1, order='F')
    return float(np.linalg.norm(v-projector(basis)@v))


def max_commutator(left, right=None):
    right = left if right is None else right
    return max(float(np.linalg.norm(a@b-b@a)) for a in left for b in right)


def noisy_readers(v, p):
    flip = np.kron(I4, np.kron(Z, I2))
    return [np.sqrt(1-p)*v, np.sqrt(p)*flip@v]


def tetrahedral_channel():
    directions = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])/np.sqrt(3)
    effects, operators = [], []
    for j, n in enumerate(directions):
        m = (I2+n[0]*X+n[1]*Y+n[2]*Z)/4
        vals, vecs = np.linalg.eigh(m)
        ket = np.zeros(16, dtype=complex)
        ket[4*j+j] = 1
        operators.append(np.outer(ket, vecs[:, -1].conj())*np.sqrt(vals[-1]))
        effects.append(m)
    return operators, effects


def input_change():
    a = np.kron(X, Z)
    return np.cos(.37)*I4-1j*np.sin(.37)*a


def trace_norm(a):
    return float(np.abs(np.linalg.eigvalsh((a+a.conj().T)/2)).sum())


def report():
    v = interaction_isometry()
    channels = [marginal_kraus([v], (4, 2, 2), (j,)) for j in range(3)]
    algebras = [recoverable(k) for k in channels]
    c = common(channels[1:])
    w = input_change()
    changed = [marginal_kraus([v@w], (4, 2, 2), (j,)) for j in (1, 2)]
    c_changed = common(changed)
    noisy = []
    for p in (0., 1e-6, .03, .2, .5, 1.):
        full = noisy_readers(v, p)
        readers = [marginal_kraus(full, (4, 2, 2), (j,)) for j in (1, 2)]
        noisy.append(dict(p=p, first_reader_exact_dimension=len(recoverable(readers[0])),
                          common_exact_dimension=len(common(readers)),
                          diamond_distance_first_reader_from_ideal=2*p,
                          promised_parity_minimax_error=min(p, 1-p)))
    tk, effects = tetrahedral_channel()
    tc = [marginal_kraus(tk, (4, 4), (j,)) for j in (0, 1)]
    return dict(
        round=391, scientific_base_through_round=390,
        hypothesis='Actual record channels determine a maximal common sharp input algebra without a preselected pointer.',
        scope=dict(ontology='All output fragments and memory carriers belong to the cognitive whole.',
                   established_theorem='Beny-Kempf-Kribs correctable algebra, applied with full input code P=I.',
                   new_project_interface='Compute the maximal common sharp algebra from existing interaction channels.',
                   additional_inputs=['Actual channel on all admitted unknown inputs',
                                      'Fixed complete input domain', 'Specified disjoint accessible output fragments'],
                   assumptions_removed=['A preselected sharp pointer algebra within the specified channel'],
                   not_derived=['Choice of natural channel', 'Fragment accessibility', 'Actual location or motion',
                                'Spatial topology or dimension', 'Future record stability', 'Finite-data certification of exact algebra']),
        existing_cz_interaction=dict(input_dimension=4, output_dimensions=[4, 2, 2], blank_reader_qubits=2,
                                     parent_reader_CZ_gates=4,
                                     recoverable_complex_dimensions=[len(a) for a in algebras],
                                     common_sharp_complex_dimension=len(c), full_joint_complex_dimension=len(recoverable([v])),
                                     inferred_common_generator='Z tensor Z, derived after solving the channel constraints',
                                     maximal_common_commutator=max_commutator(c),
                                     parent_to_reader_commutator=max_commutator(algebras[0], algebras[1]),
                                     expected_span_projector_error=float(np.linalg.norm(projector(c)-projector([I4, S])))),
        input_basis_covariance=dict(expected_generator='W^dagger (Z tensor Z) W',
                                    inferred_dimension=len(c_changed),
                                    span_error=span_error(w.conj().T@S@w, c_changed)),
        noisy_reader_cases=noisy,
        unsharp_consensus=dict(outcomes=4, reader_dimensions=[4, 4],
                               same_outcome_with_probability=1., input_effect_span_rank=4,
                               effect_max_commutator_frobenius=max_commutator(effects),
                               reader_recoverable_dimensions=[len(recoverable(k)) for k in tc],
                               common_sharp_dimension=len(common(tc))),
        numerical_policy=dict(svd_absolute_tolerance=1e-9, smallest_nonzero_noise_test=1e-6,
                              arbitrary_small_noise_proved_analytically=True,
                              numerical_rank_not_exact_experimental_certificate=True),
        sources=['https://arxiv.org/pdf/0705.1574', 'https://arxiv.org/pdf/0802.0685'])


class Checks(unittest.TestCase):
    def test_01_existing_cz_channel_and_marginals_are_trace_preserving(self):
        v = interaction_isometry()
        np.testing.assert_allclose(v.conj().T@v, I4, atol=TOL)
        for j in range(3):
            k = marginal_kraus([v], (4, 2, 2), (j,))
            np.testing.assert_allclose(adjoint(k, np.eye(k[0].shape[0])), I4, atol=TOL)

    def test_02_infer_maximal_algebras_and_retain_noncommuting_internal_information(self):
        v = interaction_isometry()
        k = [marginal_kraus([v], (4, 2, 2), (j,)) for j in range(3)]
        a = [recoverable(channel) for channel in k]
        self.assertEqual([len(x) for x in a], [8, 2, 2])
        self.assertEqual(len(recoverable([v])), 16)
        c = common(k[1:])
        np.testing.assert_allclose(projector(c), projector([I4, S]), atol=TOL)
        self.assertLess(max_commutator(c), TOL)
        self.assertLess(max_commutator(a[0], a[1]), TOL)
        for observable in (np.kron(Z, I2), np.kron(X, X)):
            self.assertLess(span_error(observable, a[0]), TOL)
        self.assertGreater(np.linalg.norm(np.kron(Z, I2)@np.kron(X, X)-np.kron(X, X)@np.kron(Z, I2)), 1.)

    def test_03_sharp_readout_intertwines_on_the_full_unknown_input_domain(self):
        v = interaction_isometry()
        for sign in (-1, 1):
            p = (I4+sign*S)/2
            e = (I2+sign*X)/2
            for j in (1, 2):
                effect = np.kron(I4, np.kron(e, I2) if j == 1 else np.kron(I2, e))
                np.testing.assert_allclose(effect@v, v@p, atol=TOL)
                k = marginal_kraus([v], (4, 2, 2), (j,))
                np.testing.assert_allclose(adjoint(k, e), p, atol=TOL)

    def test_04_kraus_representation_does_not_change_the_answer(self):
        k = marginal_kraus([interaction_isometry()], (4, 2, 2), (1,))
        n = len(k)
        u = np.exp(2j*np.pi*np.outer(np.arange(n), np.arange(n))/n)/np.sqrt(n)
        changed = [sum(u[i,j]*k[j] for j in range(n)) for i in range(n)]
        np.testing.assert_allclose(projector(recoverable(k)), projector(recoverable(changed)), atol=TOL)

    def test_05_input_coordinate_change_transforms_inferred_algebra(self):
        w = input_change()
        v = interaction_isometry()@w
        k = [marginal_kraus([v], (4, 2, 2), (j,)) for j in (1, 2)]
        np.testing.assert_allclose(w.conj().T@w, I4, atol=TOL)
        np.testing.assert_allclose(projector(common(k)), projector([I4, w.conj().T@S@w]), atol=TOL)
        self.assertGreater(np.linalg.norm(w.conj().T@S@w-S), .1)

    def test_06_arbitrarily_weak_positive_noise_destroys_exact_sharp_recovery(self):
        v = interaction_isometry()
        for p, expected in ((0., 2), (1e-6, 1), (.03, 1), (.2, 1), (.5, 1), (1., 2)):
            full = noisy_readers(v, p)
            k = [marginal_kraus(full, (4, 2, 2), (j,)) for j in (1, 2)]
            self.assertEqual(len(recoverable(k[0])), expected)
            self.assertEqual(len(common(k)), expected)
            self.assertEqual(len(recoverable(k[1])), 2)
            self.assertEqual(len(recoverable(full)), 16)

    def test_07_noise_distance_and_promised_parity_discrimination_are_distinct_tasks(self):
        v = interaction_isometry()
        ideal = marginal_kraus([v], (4, 2, 2), (1,))
        plus, minus = (I4+S)/4, (I4-S)/4
        for p in (1e-6, .03, .2, .5, 1.):
            k = marginal_kraus(noisy_readers(v, p), (4, 2, 2), (1,))
            outp, outm = apply(k, plus), apply(k, minus)
            self.assertAlmostEqual(trace_norm(outp-apply(ideal, plus)), 2*p)
            helstrom_error = .5-trace_norm(outp-outm)/4
            self.assertAlmostEqual(helstrom_error, min(p, 1-p))
            eplus = (I2+X)/2
            self.assertAlmostEqual(float(np.trace((I2-eplus)@outp).real), p)

    def test_08_unsharp_consensus_is_not_common_sharp_input_information(self):
        operators, effects = tetrahedral_channel()
        np.testing.assert_allclose(adjoint(operators, np.eye(16)), I2, atol=TOL)
        k = [marginal_kraus(operators, (4, 4), (j,)) for j in (0, 1)]
        self.assertEqual(len(common(k)), 1)
        self.assertEqual([len(recoverable(c)) for c in k], [1, 1])
        self.assertEqual(np.linalg.matrix_rank(np.column_stack([e.reshape(-1) for e in effects])), 4)
        self.assertGreater(max_commutator(effects), .1)
        for j, effect in enumerate(effects):
            output_effect = np.zeros((4, 4), dtype=complex)
            output_effect[j,j] = 1
            for channel in k:
                np.testing.assert_allclose(adjoint(channel, output_effect), effect, atol=TOL)
        rho = (I2+.2*X+.3*Y+.4*Z)/2
        output = apply(operators, rho)
        equal = np.diag([int(i//4 == i%4) for i in range(16)])
        self.assertAlmostEqual(float(np.trace(equal@output).real), 1.)

    def test_09_access_and_input_domain_are_essential(self):
        v = interaction_isometry()
        reader = marginal_kraus([v], (4, 2, 2), (1,))
        # A code promised to be even has a constant reader: no logical qubit there.
        encoding = np.eye(4, dtype=complex)[:, [0, 3]]
        restricted = [k@encoding for k in reader]
        self.assertEqual(len(recoverable(restricted)), 1)
        # Both disjoint input factors survive the identity channel, but no
        # nontrivial input sharp information is separately available in both.
        channels = [marginal_kraus([I4], (2, 2), (j,)) for j in (0, 1)]
        self.assertEqual([len(recoverable(k)) for k in channels], [4, 4])
        self.assertEqual(len(common(channels)), 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    run = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not run.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result['checks'] = dict(run=run.testsRun, failures=len(run.failures), errors=len(run.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
