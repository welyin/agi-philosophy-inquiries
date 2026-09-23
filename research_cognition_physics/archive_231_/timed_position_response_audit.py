"""Round 396: timed BB84 responses do not by themselves certify presence.

This instantiates a known position-verification attack, not a new teleportation
theorem. Geometry, signal speed, clocks, operations, and Bell preparation are
explicit inputs. Output is a destructive measurement with two classical copies;
the conditional state of any internal reference must agree with the honest map.
"""
import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET = Path(__file__).with_name('timed_position_response_audit_results.json')
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
H = np.array([[1, 1], [1, -1]], complex)/np.sqrt(2)
PHI = np.eye(2, dtype=complex)/np.sqrt(2)


def pauli(u, v):
    return np.linalg.matrix_power(X, u) @ np.linalg.matrix_power(Z, v)


def resource_vectors(noise):
    if not 0 <= noise <= 1:
        raise ValueError('White-noise fraction must lie in [0,1].')
    return [(np.sqrt(1-noise)*PHI)] + [
        np.sqrt(noise/4)*np.eye(4, dtype=complex)[j].reshape(2, 2)
        for j in range(4)]


def raw_rows(basis, noise=0.):
    """Direct contraction of the Bell bra, resource ket, and right measurement."""
    rows = []
    for u, v, y in itertools.product(range(2), repeat=3):
        bell = PHI @ pauli(u, v).T  # (I tensor P)|Phi>, indices Q,a.
        bra = basis[:, y].conj()
        for resource in resource_vectors(noise):
            transfer = np.einsum('qa,ab->bq', bell.conj(), resource)
            rows.append((u, v, y, bra @ transfer))
    return rows


def output_flag(x):
    flag = np.zeros((4, 4), complex)
    flag[3*x, 3*x] = 1
    return flag


def honest_rows(theta):
    basis = I if theta == 0 else H
    return [(x, basis[:, x].conj()) for x in range(2)]


def attack_rows(theta, noise=0.):
    basis = I if theta == 0 else H
    return [(y ^ (u if theta == 0 else v), row)
            for u, v, y, row in raw_rows(basis, noise)]


def channel(rows, operator):
    return sum((row @ operator @ row.conj()) * output_flag(x) for x, row in rows)


def referenced_output(rows, input_vector, reference_dimension):
    state = input_vector.reshape(2, reference_dimension)
    return sum(np.kron(output_flag(x), np.outer(row @ state, (row @ state).conj()))
               for x, row in rows)


def matrix_units():
    for i, j in itertools.product(range(2), repeat=2):
        yield np.outer(I[:, i], I[:, j])


def exact_audit():
    maximum = branch_error = tp_error = 0.
    for theta in range(2):
        rows = attack_rows(theta)
        tp_error = max(tp_error, float(np.linalg.norm(
            sum(np.outer(row.conj(), row) for _, row in rows)-I)))
        for operator in matrix_units():
            maximum = max(maximum, float(np.linalg.norm(
                channel(rows, operator)-channel(honest_rows(theta), operator))))
        for u, v, y, row in raw_rows(I if theta == 0 else H):
            # Zero-weight white-noise rows are not individual Bell outcomes.
            if np.linalg.norm(row) < 1e-15:
                continue
            x = y ^ (u if theta == 0 else v)
            ideal = honest_rows(theta)[x][1]
            branch_error = max(branch_error, float(np.linalg.norm(
                np.outer(row.conj(), row)-np.outer(ideal.conj(), ideal)/4)))
    return dict(channel_basis_error=maximum, branch_effect_error=branch_error,
                trace_preservation_error=tp_error, classical_outputs='identical x at both verifiers',
                input_qubit_is_measured_not_preserved=True)


def reference_audit():
    rng = np.random.default_rng(396)
    error = 0.
    for dimension in (2, 3, 5):
        state = rng.normal(size=2*dimension)+1j*rng.normal(size=2*dimension)
        state /= np.linalg.norm(state)
        for theta in range(2):
            actual = referenced_output(attack_rows(theta), state, dimension)
            ideal = referenced_output(honest_rows(theta), state, dimension)
            error = max(error, float(np.linalg.norm(actual-ideal)))
    return dict(reference_dimensions=[2, 3, 5], full_classical_reference_error=error,
                arbitrary_reference_equality_proved_by_row_identity=True)


def no_message_audit():
    rng = np.random.default_rng(1396)
    state = rng.normal(size=6)+1j*rng.normal(size=6)
    state = (state/np.linalg.norm(state)).reshape(2, 3)
    rho_r = state.T @ state.conj()
    maximum = 0.
    for theta in range(2):
        totals = [np.zeros((3, 3), complex) for _ in range(2)]
        for _, _, y, row in raw_rows(I if theta == 0 else H):
            vector = row @ state
            totals[y] += np.outer(vector, vector.conj())
        maximum = max(maximum, *(float(np.linalg.norm(value-rho_r/2)) for value in totals))
    return dict(right_outcome_and_reference_product_error=maximum,
                right_uninformed_bb84_success=.5, superluminal_message_required=False)


def noise_audit():
    cases = []
    for noise in (0., .02, .1, .5, 1.):
        maximum = extremal = 0.
        for theta in range(2):
            for operator in matrix_units():
                ideal = channel(honest_rows(theta), operator)
                random = np.trace(operator)*(output_flag(0)+output_flag(1))/2
                actual = channel(attack_rows(theta, noise), operator)
                maximum = max(maximum, float(np.linalg.norm(actual-((1-noise)*ideal+noise*random))))
            eigenstate = (I if theta == 0 else H)[:, 0]
            operator = np.outer(eigenstate, eigenstate.conj())
            difference = channel(attack_rows(theta, noise), operator)-channel(honest_rows(theta), operator)
            extremal = max(extremal, float(np.linalg.norm(difference, ord='nuc')/2))
        cases.append(dict(resource_white_noise=noise, channel_mixture_error=maximum,
                          proved_half_diamond_error=noise/2,
                          extremal_trace_distance=extremal, bb84_wrong_bit_probability=noise/2))
    return cases


def schedule(left, right, a, b, speed, delay_a, delay_b, correction, honest_delay):
    if not (0 < a < left and 0 < b < right and speed > 0):
        raise ValueError('Interceptors must lie strictly between target and verifiers.')
    if min(delay_a, delay_b, correction, honest_delay) < 0:
        raise ValueError('Processing times cannot be negative.')
    target_time = max(left, right)/speed
    local_a = target_time-a/speed+delay_a
    local_b = target_time-b/speed+delay_b
    # Messages are sent independently; neither waits to receive the other.
    received_a = local_b+(a+b)/speed
    received_b = local_a+(a+b)/speed
    reply_left = max(local_a, received_a)+correction+(left-a)/speed
    reply_right = max(local_b, received_b)+correction+(right-b)/speed
    deadlines = [target_time+left/speed+honest_delay, target_time+right/speed+honest_delay]
    excess = [reply_left-deadlines[0], reply_right-deadlines[1]]
    predicted = [max(delay_b, delay_a-2*a/speed)+correction-honest_delay,
                 max(delay_a, delay_b-2*b/speed)+correction-honest_delay]
    return dict(return_times=[reply_left, reply_right], honest_deadlines=deadlines,
                lateness=excess, formula_error=float(np.max(np.abs(np.array(excess)-predicted))),
                both_on_time=bool(max(excess) <= 1e-12),
                interceptor_positions=[-a, b], signal_distance=a+b)


def timing_audit():
    cases = [schedule(3., 5., a, b, c, da, db, dc, dh)
             for a, b, c in itertools.product((.2, 1., 2.8), (.3, 2., 4.5), (.5, 1., 2.))
             for da, db, dc, dh in ((0., 0., 0., 0.), (.002, .001, .001, .003),
                                    (.002, .001, .001, 0.), (3., .001, .1, 4.))]
    return dict(checked_schedules=len(cases), formula_error=max(x['formula_error'] for x in cases),
                exact_deadline_case=schedule(3., 5., 1., 2., 1., 0., 0., 0., 0.),
                positive_processing_case=schedule(3., 5., 1., 2., 1., .002, .001, .001, .003),
                zero_slack_rejects_positive_delay=schedule(3., 5., 1., 2., 1., .002, .001, .001, 0.))


def tilted_basis_audit():
    # This single Bell / outcome-relabel strategy is not universal for bases.
    observable = (X+Z)/np.sqrt(2)
    values, vectors = np.linalg.eigh(observable)
    basis = vectors[:, ::-1]  # outcome 0 has eigenvalue +1.
    success = 0.
    for source_bit in range(2):
        state = basis[:, source_bit]
        for u, v, y, row in raw_rows(basis):
            transformed = pauli(u, v).conj().T @ observable @ pauli(u, v)
            overlap = float(np.trace(observable @ transformed).real/2)
            guessed = y ^ int(overlap < -1e-12)
            if guessed == source_bit:
                success += abs(row @ state)**2/2
    return dict(tilted_axis=['X/sqrt(2)', 'Z/sqrt(2)'],
                mean_eigenstate_success=float(success),
                exact_one_bell_strategy_generalized_to_all_measurements=False)


@lru_cache(None)
def report():
    return dict(round=396,
                scope='A known BB84 position-verification attack with explicit channel, internal-reference, noise, and timing checks. Given effective line geometry, propagation speed, clocks and Bell resources; no dimension or universal cognition countermodel is derived.',
                exact_instrument=exact_audit(), reference=reference_audit(),
                before_classical_messages=no_message_audit(), noisy_resources=noise_audit(),
                timing=timing_audit(), non_pauli_negative_control=tilted_basis_audit(),
                resource_ledger=dict(bell_pairs_consumed_per_trial=1,
                    cross_interceptor_bits_each_direction=2, final_reply_bits_total=2,
                    challenge_qubits=1, basis_challenge_bits=1,
                    bell_resource_returned=False, replenishment_for_repetition_required=True),
                trusted_unique_local_marker_disproved=False,
                quantum_identity_authentication_disproved=False,
                response_and_timing_suffice_for_presence=False,
                full_position_generation_completed=False)


class Audit(unittest.TestCase):
    def test_complete_timed_challenge_measurement_instrument(self):
        for key in ('channel_basis_error', 'branch_effect_error', 'trace_preservation_error'):
            self.assertLess(report()['exact_instrument'][key], 3e-12)

    def test_unknown_complex_input_with_internal_reference(self):
        self.assertLess(report()['reference']['full_classical_reference_error'], 3e-12)

    def test_no_outcome_information_before_correction_messages(self):
        self.assertLess(report()['before_classical_messages']['right_outcome_and_reference_product_error'], 3e-12)

    def test_noisy_resource_full_channel_and_worst_case_witness(self):
        for case in report()['noisy_resources']:
            self.assertLess(case['channel_mixture_error'], 3e-12)
            self.assertAlmostEqual(case['extremal_trace_distance'], case['proved_half_diamond_error'])

    def test_timelike_message_schedule_and_nonzero_processing_limits(self):
        data = report()['timing']
        self.assertLess(data['formula_error'], 3e-12)
        self.assertTrue(data['exact_deadline_case']['both_on_time'])
        self.assertTrue(data['positive_processing_case']['both_on_time'])
        self.assertFalse(data['zero_slack_rejects_positive_delay']['both_on_time'])

    def test_this_single_pair_attack_is_not_a_universal_basis_claim(self):
        self.assertAlmostEqual(report()['non_pauli_negative_control']['mean_eigenstate_success'], .75)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    checked = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result = dict(report())
    result['checks'] = dict(run=checked.testsRun, failures=len(checked.failures), errors=len(checked.errors))
    result['runtime'] = dict(python=platform.python_version(), numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
