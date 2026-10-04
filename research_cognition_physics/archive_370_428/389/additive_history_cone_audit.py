"""Round 389: additive PSD summaries of independent history configurations.

This audits an order representation, not spacetime or a quantum reconstruction.
Run read-only by default; --write-results creates, never overwrites, its JSON.
"""
import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np


HERE = Path(__file__).resolve().parent
TARGET = HERE / 'additive_history_cone_audit_results.json'
TOL = 2e-10
I2 = np.eye(2, dtype=complex)
P0 = np.diag([1., 0.]).astype(complex)
P1 = I2 - P0
PLUS = np.ones((2, 2), dtype=complex) / 2
GENERATORS = (P0, P1, PLUS)


def summary(generators, counts):
    return sum((a * int(n) for a, n in zip(generators, counts)),
               np.zeros_like(generators[0]))


def psd(a):
    return bool(np.linalg.eigvalsh(a).min() >= -TOL)


def component_leq(n, m):
    return all(x <= y for x, y in zip(n, m))


def private_vector(generators, index):
    """Numerical diagnostic; the note supplies the exact support proof."""
    a = generators[index]
    b = sum((x for j, x in enumerate(generators) if j != index),
            np.zeros_like(a))
    vals, vectors = np.linalg.eigh(b)
    kernel = vectors[:, vals < 1e-11]
    if kernel.shape[1] == 0:
        return None
    restricted = kernel.conj().T @ a @ kernel
    av, au = np.linalg.eigh(restricted)
    if av[-1] <= 1e-11:
        return None
    return kernel @ au[:, -1]


def dominance_witness(generators, index):
    a = generators[index]
    b = sum((x for j, x in enumerate(generators) if j != index),
            np.zeros_like(a))
    values, vectors = np.linalg.eigh(b)
    support = vectors[:, values > 1e-11]
    assert support.shape[1] > 0
    projector = support @ support.conj().T
    assert np.linalg.norm(a - projector @ a @ projector) < TOL
    inverse_root = (support / np.sqrt(values[values > 1e-11])) @ support.conj().T
    threshold = float(np.linalg.eigvalsh(inverse_root @ a @ inverse_root)[-1])
    count = math.floor(threshold + 1e-12) + 1
    old = np.zeros(len(generators), dtype=int)
    new = np.full(len(generators), count, dtype=int)
    old[index], new[index] = 1, 0
    difference = summary(generators, new - old)
    return old, new, difference, count


def exact_qubit_difference_psd(delta):
    a, b, c = map(int, delta)
    off = Fraction(c, 2)
    left, right = Fraction(a) + off, Fraction(b) + off
    return left >= 0 and right >= 0 and left * right - off * off >= 0


def finite_table():
    states = list(itertools.product(range(3), repeat=3))
    totals = dict(pairs=0, true_comparisons=0, psd_comparisons=0,
                  false_comparisons=0, lost_comparisons=0)
    for n, m in itertools.product(states, repeat=2):
        actual = component_leq(n, m)
        coded = exact_qubit_difference_psd(np.subtract(m, n))
        numerical = psd(summary(GENERATORS, np.subtract(m, n)))
        assert coded == numerical
        totals['pairs'] += 1
        totals['true_comparisons'] += int(actual)
        totals['psd_comparisons'] += int(coded)
        totals['false_comparisons'] += int(coded and not actual)
        totals['lost_comparisons'] += int(actual and not coded)
    return totals


def nonorthogonal_family(d):
    w = np.eye(d, dtype=complex)
    for j in range(1, d):
        w[j-1, j] = .3 + .2j
    return tuple(np.outer(w[:, i], w[:, i].conj()) for i in range(d))


def ledger_operators():
    # Two qubits per branch; use counts 0, 1, 2 without wrapping.
    shift = np.roll(np.eye(4, dtype=complex), 1, axis=0)
    ops = []
    for i in range(3):
        factors = [shift if j == i else np.eye(4) for j in range(3)]
        ops.append(np.kron(np.kron(factors[0], factors[1]), factors[2]))
    return ops


class Checks(unittest.TestCase):
    def test_01_exact_history_table_and_strict_witness(self):
        table = finite_table()
        self.assertEqual(table['pairs'], 729)
        self.assertEqual(table['true_comparisons'], 216)
        self.assertEqual(table['lost_comparisons'], 0)
        self.assertGreater(table['false_comparisons'], 0)
        n, m = (0, 0, 1), (2, 2, 0)
        self.assertFalse(component_leq(n, m))
        self.assertFalse(component_leq(m, n))
        np.testing.assert_allclose(np.linalg.eigvalsh(summary(GENERATORS, np.subtract(m, n))), [1, 2])

    def test_02_injective_summary_still_has_wrong_order(self):
        columns = np.array([[a[0, 0].real, a[1, 1].real,
                             a[0, 1].real, a[0, 1].imag] for a in GENERATORS]).T
        self.assertEqual(np.linalg.matrix_rank(columns), 3)
        for n in itertools.product(range(3), repeat=3):
            x = summary(GENERATORS, n)
            c = 2 * x[0, 1].real
            decoded = (x[0, 0].real-c/2, x[1, 1].real-c/2, c)
            np.testing.assert_allclose(decoded, n)

    def test_03_private_vectors_are_absent_and_give_integer_witnesses(self):
        for index in range(3):
            self.assertIsNone(private_vector(GENERATORS, index))
            n, m, difference, _ = dominance_witness(GENERATORS, index)
            self.assertTrue(psd(difference))
            self.assertFalse(component_leq(n, m))

    def test_04_dimension_bound_is_attained_without_commutativity(self):
        for d in (2, 3, 4):
            generators = nonorthogonal_family(d)
            private = [private_vector(generators, i) for i in range(d)]
            self.assertEqual(np.linalg.matrix_rank(np.column_stack(private)), d)
            self.assertGreater(np.linalg.norm(generators[0] @ generators[1]
                                             - generators[1] @ generators[0]), .1)
            for i, v in enumerate(private):
                for j, a in enumerate(generators):
                    value = float(np.vdot(v, a @ v).real)
                    self.assertGreater(value, 0) if i == j else self.assertLess(abs(value), TOL)
            for delta in itertools.product(range(-2, 3), repeat=d):
                self.assertEqual(psd(summary(generators, delta)), min(delta) >= 0)

    def test_05_private_directions_allow_shared_support_and_higher_rank(self):
        generators = []
        for i in range(3):
            a = np.zeros((4, 4), dtype=complex)
            a[i, i], a[3, 3] = 1, 1+i/4
            a[i, 3], a[3, i] = .2j, -.2j
            generators.append(a)
        self.assertTrue(all(np.linalg.matrix_rank(a) == 2 for a in generators))
        self.assertTrue(all(private_vector(generators, i) is not None for i in range(3)))
        for delta in itertools.product(range(-2, 3), repeat=3):
            self.assertEqual(psd(summary(generators, delta)), min(delta) >= 0)

    def test_06_dominance_uses_support_not_an_illegal_full_inverse(self):
        padded = tuple(np.pad(a, ((0, 1), (0, 1))) for a in GENERATORS)
        _, _, difference, _ = dominance_witness(padded, 2)
        np.testing.assert_allclose(np.linalg.eigvalsh(difference), [0, 1, 2])
        self.assertIsNone(private_vector(padded, 2))

    def test_07_quantum_ledger_is_readable_and_preserves_unknown_reference(self):
        operators = ledger_operators()
        for a in operators:
            np.testing.assert_allclose(a.conj().T @ a, np.eye(64))
        for a, b in itertools.combinations(operators, 2):
            np.testing.assert_allclose(a @ b, b @ a)
        rng = np.random.default_rng(389)
        z = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
        rho = z @ z.conj().T
        rho /= np.trace(rho)
        initial = np.zeros((64, 64), dtype=complex)
        initial[0, 0] = 1
        for n in ((0, 0, 1), (2, 2, 0), (2, 1, 2)):
            u = np.eye(64, dtype=complex)
            for a, count in zip(operators, n):
                u = np.linalg.matrix_power(a, count) @ u
            full_u = np.kron(u, np.eye(4))
            full_rho = full_u @ np.kron(initial, rho) @ full_u.conj().T
            index = n[0]*16 + n[1]*4 + n[2]
            expected_log = np.zeros((64, 64), dtype=complex)
            expected_log[index, index] = 1
            np.testing.assert_allclose(full_rho, np.kron(expected_log, rho), atol=TOL)
            reduced = np.trace(full_rho.reshape(64, 4, 64, 4), axis1=0, axis2=2)
            np.testing.assert_allclose(reduced, rho, atol=TOL)
            pointer = np.kron(expected_log, np.eye(4))
            np.testing.assert_allclose(pointer @ full_rho @ pointer, full_rho, atol=TOL)

    def test_08_false_comparison_survives_bounded_matrix_errors(self):
        rng = np.random.default_rng(390)
        delta = .1
        ideal = 2*P0+2*P1-PLUS
        for _ in range(24):
            errors = []
            for _ in range(3):
                z = rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2))
                e = z+z.conj().T
                errors.append(delta*e/np.linalg.norm(e, 2))
            observed = ideal+2*errors[0]+2*errors[1]-errors[2]
            self.assertGreaterEqual(np.linalg.eigvalsh(observed)[0], 1-5*delta-TOL)
        physical = ((1-delta)*P0, (1-delta)*P1, (1+delta)*PLUS)
        self.assertTrue(all(psd(a) for a in physical))
        np.testing.assert_allclose(np.linalg.eigvalsh(summary(physical, (2, 2, -1))), [.7, 1.8])
        # Unconstrained Hermitian calibration errors can saturate the bound.
        perturbed = (P0-delta*PLUS, P1-delta*PLUS, PLUS+delta*PLUS)
        # Saturation concerns Hermitian calibration errors; these first two
        # matrices need not themselves be physical positive generators.
        np.testing.assert_allclose(np.linalg.eigvalsh(summary(perturbed, (2, 2, -1))), [.5, 2])

    def test_09_small_positive_leakage_hides_failure_until_longer_histories(self):
        epsilon = Fraction(1, 16)
        def exact_noisy(delta):
            a, b = delta
            return a+(a+b)*epsilon >= 0 and b+(a+b)*epsilon >= 0
        for delta in itertools.product(range(-16, 17), repeat=2):
            self.assertEqual(exact_noisy(delta), min(delta) >= 0)
        self.assertTrue(exact_noisy((17, -1)))
        self.assertFalse(exact_noisy((16, -1)))
        noisy = (P0+float(epsilon)*I2, P1+float(epsilon)*I2)
        self.assertIsNone(private_vector(noisy, 0))
        self.assertIsNone(private_vector(noisy, 1))
        np.testing.assert_allclose(np.linalg.eigvalsh(summary(noisy, (18, -1))), [1/16, 305/16])

    def test_10_congruence_recalibration_cannot_remove_false_order(self):
        w = np.array([[1, .2+.3j], [.1j, 1.1]], dtype=complex)
        changed = tuple(w @ a @ w.conj().T for a in GENERATORS)
        for delta in itertools.product(range(-2, 3), repeat=3):
            self.assertEqual(psd(summary(changed, delta)),
                             exact_qubit_difference_psd(delta))
        self.assertGreater(np.linalg.eigvalsh(summary(changed, (2, 2, -1)))[0], 0)


def results(checks):
    table = finite_table()
    witnesses = []
    for i in range(3):
        old, new, difference, count = dominance_witness(GENERATORS, i)
        witnesses.append(dict(index=i, earlier_candidate=old.tolist(),
                              later_candidate=new.tolist(), positive_multiplier=count,
                              eigenvalues=np.linalg.eigvalsh(difference).tolist()))
    return dict(
        round=389, scientific_base_through_round=388,
        hypothesis='A fixed additive PSD summary reflects all independent finite history extensions.',
        scope=dict(object='Global record configurations / consistent cuts of independent chains, not individual spacetime events.',
                   additional_inputs=['Independent repeatable recorded updates', 'Fixed PSD increment per update type',
                                      'Additive unnormalized summary', 'Order tested by PSD differences'],
                   conclusion='Private kernel directions are necessary and sufficient; k <= matrix size D.',
                   excludes='The stated cumulative-history-to-cone identification.',
                   does_not_exclude=['A local event displacement cone', 'Physical 3+1 spacetime',
                                     'Nonlinear or history-dependent encodings', 'The full cognitive programme']),
        finite_table=table,
        strict_witness=dict(n=[0, 0, 1], m=[2, 2, 0], difference=[[1.5, -.5], [-.5, 1.5]],
                            eigenvalues=[1., 2.], summary_real_rank=3, summaries_injective=True,
                            per_generator_error=.1, certified_margin=.5),
        constructive_witnesses=witnesses,
        positive_cases=dict(noncommuting_saturated_dimensions=[2, 3, 4],
                            higher_rank_case=dict(matrix_size=4, independent_histories=3, generator_ranks=[2, 2, 2])),
        ledger=dict(branches=3, maximum_local_count=2, register_dimensions=[4, 4, 4],
                    record_qubits=6, configurations_used=27,
                    log_dimension=64, test_payload_reference_dimension=4,
                    complete_test_dimension=256, updates_across_two_witness_preparations=5,
                    pointer_readout='Existing log basis; output copying would require additional blank memory.',
                    global_schedule_autonomously_derived=False),
        leakage=dict(epsilon=1/16, all_integer_differences_through=16,
                     first_false_ratio=17, strict_witness_ratio=18,
                     strict_eigenvalues=[1/16, 305/16]),
        sources=['https://arxiv.org/pdf/1407.4095'],
        checks=checks, runtime=dict(python=platform.python_version(), numpy=np.__version__))


def report():
    """Scientific payload consumed by the existing read-only audit chain."""
    return {k: v for k, v in results(None).items() if k not in ('checks', 'runtime')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    run = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks))
    if not run.wasSuccessful():
        raise SystemExit(1)
    report = results(dict(run=run.testsRun, failures=len(run.failures), errors=len(run.errors)))
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
