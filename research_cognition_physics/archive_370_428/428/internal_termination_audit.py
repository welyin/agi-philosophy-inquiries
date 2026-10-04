"""Round 428: fixed quantum processor, internal retry records and resource scope.

Finite tests reproduce Vidal--Masanes--Cirac's known U(1) processor. The
countable-termination/no-programming implication is proved in the note, not
inferred from these matrices. Only NumPy and the standard library are used.
"""
import argparse
import io
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

TARGET = Path(__file__).with_name('internal_termination_audit_results.json')
OBS = {}
I = np.eye(2, dtype=complex)
ZERO = np.diag([1., 0.]).astype(complex)


def phase(theta):
    return np.diag([1., np.exp(1j*theta)])


def bank(theta, count):
    out = np.ones(1, dtype=complex)
    for k in range(count):
        out = np.kron(out, np.array([1., np.exp(1j*(2**k)*theta)])/math.sqrt(2))
    return out


def processor(count):
    """Kraus operators D x (D tensor P); the gates never depend on theta.

    Keep all readout bits coherently. Flip bit k if the data bit is one and
    every earlier readout bit is one (failure). Final tracing implements the
    measurement-and-feedback instrument, including unused program bits.
    """
    size = 2**count
    perm = np.arange(2*size)
    for k in range(count):
        bit = 1 << (count-1-k)
        prefix_mask = ((1 << k)-1) << (count-k)
        active = (perm >= size) & ((perm % size & prefix_mask) == prefix_mask)
        perm[active] ^= bit
    assert len(np.unique(perm)) == 2*size
    inverse = np.argsort(perm)
    operators = []
    for record in range(size):
        operator = np.zeros((2, 2*size), dtype=complex)
        operator[np.arange(2), inverse[np.arange(2)*size+record]] = 1
        operators.append(operator)
    return operators


def effective(theta, count):
    prepare = np.kron(I, bank(theta, count)[:, None])
    return [operator @ prepare for operator in processor(count)]


def channel(operators, rho, reference=1):
    out = np.zeros_like(rho, dtype=complex)
    for operator in operators:
        lifted = np.kron(operator, np.eye(reference))
        out += lifted @ rho @ lifted.conj().T
    return out


def choi(operators):
    bell = np.array([1., 0., 0., 1.])/math.sqrt(2)
    return channel(operators, np.outer(bell, bell), reference=2)


def trace_distance(a, b):
    delta = a-b
    return float(np.sum(np.abs(np.linalg.eigvalsh((delta+delta.conj().T)/2)))/2)


def first_stop(record, count):
    for k in range(count):
        if not (record & (1 << (count-1-k))):
            return k+1
    return None


class Audit(unittest.TestCase):
    def test_01_single_gate_and_unknown_reference(self):
        rng = np.random.default_rng(42801)
        worst = 0.
        for theta in (-2.1, 0., .37, math.pi, 2.7):
            operators = effective(theta, 1)
            targets = [phase(theta)/math.sqrt(2),
                       np.exp(1j*theta)*phase(-theta)/math.sqrt(2)]
            for actual, expected in zip(operators, targets):
                worst = max(worst, float(np.linalg.norm(actual-expected)))
                self.assertTrue(np.allclose(actual, expected, atol=1e-13))
            for _ in range(4):
                psi = rng.normal(size=6)+1j*rng.normal(size=6)
                psi /= np.linalg.norm(psi)
                rho = np.outer(psi, psi.conj())
                for operator, target in zip(operators, targets):
                    output = channel([operator], rho, reference=3)
                    self.assertAlmostEqual(float(np.trace(output).real), .5)
                    self.assertLess(trace_distance(output, channel([target], rho, 3)), 1e-13)
        OBS['single_gate_matrix_error'] = worst

    def test_02_fixed_feedback_circuit_and_internal_records(self):
        worst = 0.
        cases = 0
        for count in range(1, 7):
            for theta in (-.91, .17, .73, 2.1):
                operators = effective(theta, count)
                self.assertTrue(np.allclose(sum(k.conj().T @ k for k in operators), I))
                target_choi = choi([phase(theta)])
                for stop in range(1, count+1):
                    group = [k for m, k in enumerate(operators) if first_stop(m, count) == stop]
                    error = float(np.linalg.norm(choi(group)-2.**(-stop)*target_choi))
                    worst = max(worst, error)
                    self.assertLess(error, 1e-12)
                failure_target = 2.**(-count)*choi([phase(-(2**count-1)*theta)])
                self.assertLess(np.linalg.norm(choi([operators[-1]])-failure_target), 1e-12)
                cases += 1
        OBS['feedback_circuit'] = dict(cases=cases, maximum_first_stop_instrument_error=worst,
                                       target_independent_gates=True, full_reference_choi_checked=True)

    def test_03_trace_defect_completion_on_all_program_inputs(self):
        # This audit uses the full processor, not just the valid phase programs.
        # Completion erases a failed output; it differs from returning that output.
        count = 3
        operators = processor(count)
        success, fail = operators[:-1], operators[-1]
        effect = sum(k.conj().T @ k for k in success)
        defect = np.eye(effect.shape[0])-effect
        self.assertTrue(np.allclose(defect, fail.conj().T @ fail))
        completion = []
        for row in fail:
            operator = np.zeros_like(fail)
            operator[0] = row
            completion.append(operator)
        all_ops = success+completion
        completeness_error = float(np.linalg.norm(sum(k.conj().T @ k for k in all_ops)-np.eye(effect.shape[0])))
        self.assertLess(completeness_error, 1e-13)
        rng = np.random.default_rng(42803)
        psi = rng.normal(size=effect.shape[0])+1j*rng.normal(size=effect.shape[0])
        psi /= np.linalg.norm(psi)
        x = np.outer(psi, psi.conj())
        apply = lambda ops: sum(k @ x @ k.conj().T for k in ops)
        actual = apply(all_ops)
        expected = apply(success)+np.trace(defect @ x)*ZERO
        self.assertTrue(np.allclose(actual, expected))
        self.assertAlmostEqual(float(np.trace(actual).real), 1.)
        self.assertGreaterEqual(float(np.linalg.eigvalsh(actual).min()), -1e-13)
        OBS['trace_defect_completion'] = dict(input_dimension=effect.shape[0],
            defect_rank=int(np.linalg.matrix_rank(defect)), completeness_error=completeness_error,
            tested_on_entangled_data_program_input=True, completion_returns_fixed_failure_state=True)

    def test_04_unconditional_error_including_failure(self):
        rng = np.random.default_rng(42804)
        plus = np.ones(2)/math.sqrt(2)
        plus_state = np.outer(plus, plus)
        worst = 0.
        rows = []
        for count in (1, 2, 4, 6):
            for theta in (.13, .73, 1.91):
                operators = effective(theta, count)
                q = 2.**(-count)
                ideal = [phase(theta)]
                expected = q*abs(math.sin(2**(count-1)*theta))
                observed = trace_distance(channel(operators, plus_state), channel(ideal, plus_state))
                worst = max(worst, abs(observed-expected))
                self.assertAlmostEqual(observed, expected, places=12)
                self.assertAlmostEqual(trace_distance(choi(operators), choi(ideal)), expected, places=12)
                for _ in range(5):
                    matrix = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
                    rho = matrix @ matrix.conj().T
                    rho /= np.trace(rho)
                    actual = channel(operators, rho, reference=4)
                    target = channel(ideal, rho, reference=4)
                    mixture = (1-q)*target+q*channel([phase(-(2**count-1)*theta)], rho, 4)
                    self.assertLess(np.linalg.norm(actual-mixture), 1e-12)
                    self.assertLessEqual(trace_distance(actual, target), expected+1e-12)
                if theta == .73:
                    rows.append(dict(program_qubits=count, failure_probability=q,
                                     full_half_diamond_error=expected, saturating_plus_input_error=observed))
        OBS['cutoff_channel'] = dict(maximum_saturating_witness_error=worst, examples=rows,
                                    arbitrary_reference_bound_proved_in_note=True,
                                    random_inputs_not_used_as_diamond_proof=True)

    def test_05_program_overlap_and_source_information(self):
        worst = 0.
        rows = []
        theta, phi = .13, .83
        delta = theta-phi
        for count in (1, 2, 4, 8, 10):
            overlap = np.vdot(bank(phi, count), bank(theta, count))
            formula = np.expm1(1j*(2**count)*delta)/(2**count*np.expm1(1j*delta))
            worst = max(worst, float(abs(overlap-formula)))
            self.assertLess(abs(overlap-formula), 1e-12)
            upper = 1/(2**count*abs(math.sin(delta/2)))
            self.assertLessEqual(abs(overlap), upper+1e-12)
            rows.append(dict(program_qubits=count, absolute_overlap=float(abs(overlap)),
                             analytic_upper_bound=upper))
        # A concrete trace-distance obstruction to making two complete program
        # prefixes from one nonorthogonal source token by a fixed channel.
        source_overlap = abs(np.vdot(bank(0., 1), bank(math.pi/2, 1)))
        output_overlap = abs(np.vdot(bank(0., 2), bank(math.pi/2, 2)))
        source_distance = math.sqrt(max(0., 1-source_overlap**2))
        output_distance = math.sqrt(max(0., 1-output_overlap**2))
        self.assertGreater(output_distance, source_distance+.2)
        OBS['program_sources'] = dict(overlap_identity_error=worst, prefix_overlaps=rows,
            single_token_trace_distance=source_distance, two_token_prefix_trace_distance=output_distance,
            exact_generation_would_violate_trace_distance_contractivity=True,
            pairwise_orthogonal_limit_is_analytic_not_a_finite_sample_inference=True)

    def test_06_consumed_length_and_prepared_capacity(self):
        rows = []
        for count in (1, 2, 4, 8):
            operators = effective(.73, count)
            consumed = 0.
            success = 0.
            for record, operator in enumerate(operators):
                probability = float(np.trace(operator.conj().T @ operator).real/2)
                stop = first_stop(record, count)
                consumed += probability*(count if stop is None else stop)
                success += probability*(stop is not None)
            self.assertAlmostEqual(success, 1-2.**(-count), places=12)
            self.assertAlmostEqual(consumed, 2*(1-2.**(-count)), places=12)
            rows.append(dict(prepared_program_qubits=count, program_dimension=2**count,
                             mean_consumed_program_qubits=consumed, heralded_success_probability=success))
        OBS['resource_accounting'] = dict(examples=rows, asymptotic_mean_consumed_qubits=2,
            source_preparation_cost_derived=False, mean_consumption_not_a_capacity_bound=True,
            unused_program_registers_and_retry_records_remain_inside_whole=True)


def run():
    output = io.StringIO()
    tests = unittest.TextTestRunner(stream=output).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():
        raise RuntimeError(output.getvalue())
    return dict(round=428, baseline_round=427, status='internal_termination_programming_scope_audit',
        tests_run=tests.testsRun, failures=len(tests.failures), errors=len(tests.errors),
        python=platform.python_version(), numpy=np.__version__, observations=OBS,
        scope=dict(internal_generation_is_user_authorized_candidate=True,
            countable_first_termination_reduction_proved=True,
            separable_normal_program_contract_explicit=True,
            all_target_information_accounted_at_initial_interface=True,
            growing_memory_and_random_retries_allowed=True,
            common_deadline_or_fixed_hamiltonian_required=False,
            exact_universal_almost_sure_processor_under_contract=False,
            candidate_principle_itself_refuted=False,
            known_phase_processor_claimed_as_new_discovery=False,
            full_cognitive_implementation_completed=False,
            actual_position_interface_still_required=True,
            three_dimensional_space_unconditionally_derived=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.check:
        assert json.loads(TARGET.read_text(encoding='utf-8')) == result
    else:
        with TARGET.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
