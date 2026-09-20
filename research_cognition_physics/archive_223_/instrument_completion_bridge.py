"""Round 228: finite instruments compiled into a unitary and ancilla POVM.

Actual exact permission is conditional on round 227 C_path or closed groups.
The base contract gives arbitrary approximation, with a proved 2*delta bound.
NumPy checks entire Choi matrices, including rectangular and zero branches.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest

import numpy as np

from reversible_dynamics_bridge import dagger, sample_generator, sample_state, unitary
from tensor_process_bridge import apply_kraus, choi, extend_left, maximally_entangled
from steering_permission_bridge import positive_power


def random_kraus(a, b, count):
    rng = np.random.default_rng(22800+100*a+10*b+count)
    raw = [rng.normal(size=(b, a))+1j*rng.normal(size=(b, a)) for _ in range(count)]
    normalizer = positive_power(sum(dagger(k)@k for k in raw), -.5)
    return [k@normalizer for k in raw]


def auxiliary_dimension(base_dim, count):
    """Use powers of an existing nontrivial type; no new integer type assumed."""
    if base_dim < 2:
        raise ValueError('Use an existing nontrivial object as the base')
    e = base_dim
    while e < count:
        e *= base_dim
    return e


def compile_instrument(kraus, labels, outcomes, e):
    if not kraus or len(labels) != len(kraus):
        raise ValueError('Each Kraus operator needs one label')
    b, a = kraus[0].shape
    if any(k.shape != (b, a) for k in kraus):
        raise ValueError('All branches must have the same input/output type')
    if e < len(kraus) or any(j < 0 or j >= outcomes for j in labels):
        raise ValueError('Insufficient environment or invalid outcome label')
    if np.linalg.norm(sum(dagger(k)@k for k in kraus)-np.eye(a)) > 1e-10:
        raise ValueError('The sum of the instrument must be trace preserving')
    n = a*b*e
    # Columns describe |i>_A |0>_B |0>_E -> |0>_A sum_l K_l|i>_B |l>_E.
    w = np.zeros((a, b, e, a), dtype=complex)
    for ell, k in enumerate(kraus):
        w[0, :, ell, :] = k
    w = w.reshape(n, a)
    left, _, _ = np.linalg.svd(w, full_matrices=True)
    prepared_columns = np.arange(a)*b*e
    other_columns = np.setdiff1d(np.arange(n), prepared_columns)
    u = np.zeros((n, n), dtype=complex)
    u[:, prepared_columns] = w
    u[:, other_columns] = left[:, a:]
    # Complete the ancilla POVM on the unused subspace, assigning it outcome 0.
    environment_labels = np.zeros(e, dtype=int)
    environment_labels[:len(labels)] = labels
    return u, environment_labels


def implemented_branches(u, a, b, e, environment_labels, outcomes):
    prepared_columns = np.arange(a)*b*e
    tensor = u[:, prepared_columns].reshape(a, b, e, a)
    return [[tensor[old_a, :, ell, :] for old_a in range(a)
             for ell in range(e) if environment_labels[ell] == outcome]
            for outcome in range(outcomes)]


def branch_action(branch, x, b):
    return sum((k@x@dagger(k) for k in branch), np.zeros((b, b), dtype=complex))


def branch_chois(branches, a, b):
    return [choi(lambda x, branch=branch: branch_action(branch, x, b), a)
            for branch in branches]


def complete_event(kraus):
    b, a = kraus[0].shape
    gap = np.eye(a)-sum(dagger(k)@k for k in kraus)
    vals, vecs = np.linalg.eigh(gap)
    if vals.min() < -1e-10:
        raise ValueError('Event is not trace nonincreasing')
    root = (vecs*np.sqrt(np.maximum(vals, 0)))@dagger(vecs)
    output = np.eye(b)[:, 0]
    failure = [np.outer(output, root[j, :]) for j in range(a)]
    return kraus+failure, [0]*len(kraus)+[1]*a


@lru_cache(None)
def certificate(a, b, count, outcomes=4):
    kraus = random_kraus(a, b, count)
    labels = [j % (outcomes-1) for j in range(count)]
    e = auxiliary_dimension(a, count)
    u, environment_labels = compile_instrument(kraus, labels, outcomes, e)
    expected = [[k for k, j in zip(kraus, labels) if j == outcome] for outcome in range(outcomes)]
    actual = implemented_branches(u, a, b, e, environment_labels, outcomes)
    targets, got = branch_chois(expected, a, b), branch_chois(actual, a, b)
    return {'input_dim': a, 'output_dim': b, 'kraus_count': count,
            'environment_dim': e, 'unitary_dim': len(u), 'outcomes': outcomes,
            'unitarity_error': float(np.linalg.norm(dagger(u)@u-np.eye(len(u)))),
            'max_branch_choi_error': float(max(np.linalg.norm(x-y) for x, y in zip(targets, got))),
            'tp_error': float(np.linalg.norm(sum(dagger(k)@k for branch in actual for k in branch)-np.eye(a))),
            'zero_outcome_choi_norm': float(np.linalg.norm(got[-1])),
            'unused_environment_dimensions': e-count}


def approximate_certificate(epsilon):
    a, b, e, outcomes = 2, 3, 4, 3
    kraus = random_kraus(a, b, 3)
    u, labels = compile_instrument(kraus, [0, 1, 1], outcomes, e)
    k = sample_generator(len(u), seed=2283)
    k /= np.linalg.norm(k, 2)
    perturbed = unitary(k, epsilon)@u
    delta = np.linalg.norm(perturbed-u, 2)
    target = branch_chois(implemented_branches(u, a, b, e, labels, outcomes), a, b)
    approx = branch_chois(implemented_branches(perturbed, a, b, e, labels, outcomes), a, b)
    error = sum(np.abs(np.linalg.eigvalsh((x-y)/a)).sum() for x, y in zip(target, approx))
    return {'epsilon': epsilon, 'unitary_operator_norm_error': float(delta),
            'flagged_trace_norm_error_on_maximally_entangled_input': float(error),
            'analytic_diamond_norm_upper_bound': float(2*delta),
            'scope': 'Perturbed unitaries test the stability bound; no claim that these particular perturbations belong to an unknown allowed gate set.'}


class InstrumentCompletionTests(unittest.TestCase):
    def test_rectangular_and_padded_instruments_match_every_branch_choi(self):
        for a, b, count in ((2, 3, 5), (3, 2, 4), (3, 4, 6), (2, 2, 3)):
            c = certificate(a, b, count)
            self.assertLess(c['unitarity_error'], 2e-12)
            self.assertLess(c['max_branch_choi_error'], 2e-12)
            self.assertLess(c['tp_error'], 2e-12)
            self.assertEqual(c['zero_outcome_choi_norm'], 0)

    def test_rank_one_isometry_and_one_nonzero_outcome(self):
        c = certificate(2, 3, 1)
        self.assertLess(c['unitarity_error'], 2e-12)
        self.assertLess(c['max_branch_choi_error'], 2e-12)

    def test_direct_joint_state_calculation_preserves_an_entangled_reference(self):
        a, b, e = 2, 3, 4
        kraus = random_kraus(a, b, 3)
        labels_in = [0, 1, 1]
        u, labels = compile_instrument(kraus, labels_in, 3, e)
        reference = sample_state(a*a, seed=2285)
        reference = .6*maximally_entangled(a)+.4*reference
        embedding = np.zeros((a*b*e, a))
        embedding[np.arange(a)*b*e, np.arange(a)] = 1
        preparation = np.kron(embedding, np.eye(a))
        joint_u = np.kron(u, np.eye(a))
        final = joint_u@preparation@reference@dagger(preparation)@dagger(joint_u)
        tensor = final.reshape(a, b, e, a, a, b, e, a)
        for outcome in range(3):
            actual = np.zeros((b, a, b, a), dtype=complex)
            for old_a in range(a):
                for ell in np.flatnonzero(labels == outcome):
                    actual += tensor[old_a, :, ell, :, old_a, :, ell, :]
            branch = [k for k, j in zip(kraus, labels_in) if j == outcome]
            expected = extend_left(lambda x: branch_action(branch, x, b), reference, a, a)
            np.testing.assert_allclose(actual.reshape(b*a, b*a), expected, atol=2e-13)

    def test_amplitude_damping_click_probabilities_and_states(self):
        gamma = .37
        kraus = [np.diag([1, np.sqrt(1-gamma)]), np.array([[0, np.sqrt(gamma)], [0, 0]])]
        u, labels = compile_instrument(kraus, [0, 1], 2, 2)
        branches = implemented_branches(u, 2, 2, 2, labels, 2)
        excited = np.diag([0., 1.])
        np.testing.assert_allclose(branch_action(branches[0], excited, 2), np.diag([0, 1-gamma]), atol=1e-14)
        np.testing.assert_allclose(branch_action(branches[1], excited, 2), np.diag([gamma, 0]), atol=1e-14)

    def test_swapping_to_the_designated_auxiliary_transfers_complex_povms(self):
        d = 2
        swap = np.eye(d*d).reshape(d, d, d, d).transpose(1, 0, 2, 3).reshape(d*d, d*d)
        y_plus = np.array([[1, -1j], [1j, 1]])/2
        blank = np.diag([1., 0.])
        state = y_plus
        final = swap@np.kron(state, blank)@dagger(swap)
        probability = np.trace(np.kron(np.eye(d), y_plus)@final).real
        self.assertAlmostEqual(probability, 1)
        self.assertAlmostEqual(np.trace(y_plus.T@state).real, 0)

    def test_trace_nonincreasing_event_can_be_completed_with_a_failure_record(self):
        event = [.6*k for k in random_kraus(3, 2, 2)]
        kraus, labels = complete_event(event)
        u, environment_labels = compile_instrument(kraus, labels, 2, 9)
        actual = implemented_branches(u, 3, 2, 9, environment_labels, 2)
        expected = choi(lambda x: apply_kraus(event, x), 3)
        got = choi(lambda x: branch_action(actual[0], x, 2), 3)
        np.testing.assert_allclose(got, expected, atol=2e-13)

    def test_invalid_trace_and_environment_contracts_are_rejected(self):
        with self.assertRaises(ValueError):
            compile_instrument([2*np.eye(2)], [0], 1, 2)
        with self.assertRaises(ValueError):
            compile_instrument(random_kraus(2, 2, 3), [0, 1, 1], 2, 2)
        with self.assertRaises(ValueError):
            complete_event([2*np.eye(2)])

    def test_misreading_the_record_changes_the_instrument_even_if_channel_agrees(self):
        kraus = random_kraus(2, 3, 3)
        u, labels = compile_instrument(kraus, [0, 1, 1], 2, 4)
        correct = branch_chois(implemented_branches(u, 2, 3, 4, labels, 2), 2, 3)
        wrong = branch_chois(implemented_branches(u, 2, 3, 4, 1-labels, 2), 2, 3)
        np.testing.assert_allclose(sum(correct), sum(wrong), atol=1e-14)
        self.assertGreater(np.linalg.norm(correct[0]-wrong[0]), .1)

    def test_approximation_bound_covers_flagged_entangled_outputs(self):
        for epsilon in (.1, .01, .001):
            c = approximate_certificate(epsilon)
            self.assertLessEqual(c['flagged_trace_norm_error_on_maximally_entangled_input'], c['analytic_diamond_norm_upper_bound']+1e-12)


def report():
    return {'round': 228, 'date': '2026-09-20',
            'exact_result': 'F+U+C+P with C_path (or closed reversible groups) yields every finite CP instrument between existing object types, hence all POVMs and CPTP channels. Actual events were already CP by round225.',
            'baseline_result': 'Without that extra regularity, every such instrument is arbitrarily approximable in diamond norm with its classical outcome retained.',
            'construction_inputs': ['round227 unitary access with exact/approximate scope', 'round226 full POVMs on designated steering auxiliaries', 'F finite composition, preparation, discard and conditioning'],
            'instrument_certificates': [certificate(a, b, count) for a, b, count in ((2, 3, 5), (3, 2, 4), (3, 4, 6), (2, 2, 3), (2, 3, 1))],
            'approximation_certificates': [approximate_certificate(epsilon) for epsilon in (.1, .01, .001)],
            'runtime': {'python': platform.python_version(), 'numpy': np.__version__},
            'not_derived': ['existence of every integer-dimensional elementary object', 'a finite control cost bound', 'a chosen Hamiltonian, time scale, spacetime or gravity'],
            'next': 'Audit whether the exact regularity condition can be weakened within the full composite theory; distinguish operational completeness from selection of local dynamics.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(InstrumentCompletionTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    data = report()
    data['checks'] = {'run': checks.testsRun, 'failures': 0, 'errors': 0}
    if args.write_results:
        target = Path(__file__).with_name('instrument_completion_bridge_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8')) != data:
            raise RuntimeError('Existing result differs; review before replacing it.')
        target.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
